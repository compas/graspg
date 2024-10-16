!***********************************************************************
!                                                                      *
      SUBROUTINE SETHAM_POT (EOL, JBLOCK)
!----------------------------------------------------------------------*
! CSFG version of SETHAM_POT                                           *
!      Written by Chongyang Chen, Fudan university, Dec 2023           *
!                                                                      *
! This subroutine calls SETHAM_CYC, performs the following calculations*
!   0) IF (LREADMCP) the first call of setham_pot, reading mcpfiles    *
!   1) IF (LABTVFRST):                                                 *
!        Construct NXAOPT, NYAOPT, NDAOPT arrays                       *
!   2) IF (LCHM) Construct the H-matrix, diagonalize them in Matrixmpi *
!   3) IF (LCPOT) Calculate XAopt and YAopt, and DAopt                 *
!     These data are used to calculate the potentials.                 *
!                                                                      *
!     Arrays: VCOEFP, KVCOFP, NVKPG, LABPV, and ICIRPV are not used to *
!     save memory                                                      *      
!     Arrays: TCOEFP, NCRPST, LABPT, and ICIRPT are not used anymore   *      
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE, LONG, BYTE
      USE parameter_def,   ONLY: KEYORB
      USE csfg_memory_man
      use symexpand_mod
      use symmatrix_mod, CUTOFF0=>CUTOFF
      use, intrinsic :: iso_fortran_env
      !use csfg_scf_C   !,  ONLY: NPASTCSF, NPASTCSFG
      use csfg_scf_C    ,  ONLY: NXCOFOPT, NYCOFOPT, XAOPT, YAOPT, &
                                 NXAOPT,   NXAIORB,  JXIPOS,       &
                                 NDCOFOPT, NDAOPT,   DAOPT,        &
                                 XAWRK,    YAWRK
      use csfg_tv_C
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      USE hmat_C
      !USE orb_C, ONLY: NCFTr
      USE iounit_C
      USE iccu_C
      USE mpi_C
      
!CYC      USE decide_C,  ONLY: LFORDR
      USE mcp_C,           ONLY: LFORDR
      use csfg_decide_C,   ONLY: LSYMD
      USE debug_C
      USE buffer_C
      USE hblock_C

      USE bcore_C, ONLY: icore
      USE bilst_C
      USE coeils_C
      USE def_C
      USE keilst_C
      USE prnt_C
      USE stat_C
      USE stor_C
      USE wave_C
      USE cteilsrk_C
      USE blim_c
      USE orb_C, ONLY: NW

      USE symmatrix_mod, ONLY: NCFTr
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I
      USE convrt_I
      USE ichop_I
      IMPLICIT NONE
      EXTERNAL CORD
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL             :: EOL
      INTEGER             :: JBLOCK
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW) :: tshell
      REAL(DOUBLE) :: tgrl1, tgrl2, tegral
      INTEGER, PARAMETER :: KEY = KEYORB
!
!     Matrix elements smaller than CUTOFF are not accumulated
!     CUTOFF is set within symmatrix_mod.f90
!
      INTEGER :: I, J, K, NELC, IC, IR, IA, IB,                 &
                 itype, iia, ICC, IOUTF, IRR, NCFGEN, NCFSMALL, &
                 NCFTOT, IC0, IR0, IOUTT, IORB, NXCOF
      INTEGER :: ICUT, NB, NCOLG, NROWG, ICOLG, IROWG, NCOLS, NROWS
      REAL(DOUBLE) :: endtime, starttime

      INTEGER :: NPAIRtmp(2), NELMNTCC
      INTEGER :: IENDW
!-----------------------------------------------------------------------
      !write(*,*)'LSYMD,  LFORDR =', LSYMD,  LFORDR 
      NB = JBLOCK
      IF (NB == 1) THEN
        IF (LREADMCP) THEN
          REWIND(IGMCP31)
          REWIND(IGMCP32)
          READ(IGMCP31)
          READ(IGMCP32)
        ENDIF 
        NUMHT  = 0
        NUMHP  = 0
        NUMHVK = 0
! One-body integral initiallization
        FRSTCO = .TRUE.
        NCOEI = 0
! Initialization for NDCOFOPT, NYCOFOPT, NXCOFOPT
        IF (LABTVFRST) THEN
          NDCOFOPT = 0
          NYCOFOPT = 0
          NXCOFOPT = 0
! Initialization for DAOPT, YAOPT, XAOPT
        ELSEIF (LCPOT) THEN
          DAOPT = 0.D0
          !YAOPT = 0.D0
          !XAOPT = 0.D0
          ! XAWRK is allocated in VECTOR_XYA 
          XAWRK = 0.D0
          ! YAWRK is allocated in VECTOR_XYA 
          YAWRK = 0.D0
        ENDIF
      ENDIF
! Initialization the beginning CSF / CSFG of JBLOCK
      IF (NB == 1) THEN
        NCSFPAST = 0
        NCFGPAST = 0
      ELSE
        NCSFPAST = SUM(TotCSFs_perblock(3,1:NB-1))
        NCFGPAST = SUM(TotCSFs_perblock(2,1:NB-1))
      ENDIF

! Initialization for the CSFG (COLUMNS) for MYID process: ICGPDET and
! Initialization for the detailed interacting CSFG pairs ICGPEND
! Both are in CSFG list.
      NPAIRPAST(1) = SUM(NPAIRS(1, 0:NB-1)) 

! Initialization for IRGPDET in CSFG list.
      NPAIRPAST(2) = SUM(NPAIRS(2, 0:NB-1)) 

! Initialization for NODECOLS and IENDCTV (in normal CSF list),
! NCOLBLK(:) record the values of NODENCOLS.
! NCOLPAST differs from NCSFPAST
      NCOLPAST = SUM(NCOLBLK(0:NB-1))

!! Initialization for IROWTV, NZPAST: Non-Zero matrixelements
!! Initilized but never used again?
!      NZPAST = SUM(NZBLK(0:NB-1))

! The number of CSFG for the JBLOCK-th block.
      NCFSMALL = TotCSFs_perblock(2,NB)
      NCFTr = TotCSFs_perblock(3,NB)

! Allocate working arrays...
! NCOLG: CSFG-Columns calculated by MYID process
      NCOLG = NPAIRS(1, NB)

! NCOLS: CSF-Columns calculated by MYID process
      NCOLS = NCOLBLK(NB)

! IROWSYM: Row index of the non-zero elements for the MAXSPAN COLUMN-CSFs  
      ALLOCATE(IROWSYM(NCFTr,MAXSPAN))
! Matrixelement values of the ICOLG-CSFG interacting all the 
! IROWG (1 - ICOLG) CSFGs 
      ALLOCATE(EMTSYM(NCFTr,MAXSPAN))

! Sub-block of matrix generated by each ICOLG - IROWG CSFGs pair.
      ALLOCATE(EMTBLOCK(MAXSPAN,MAXSPAN))

! Number of the non-zero elements belong to the ICOLG CSFG
      ALLOCATE(NELCSYM(MAXSPAN))

! NODECOLS used in spicmvsym.f, allocated in SetArrs_csfg.f90, 
! deallocated in matrix_csfg.f90
!      ALLOCATE(NODECOLS(NCOLS))

! Open scratch file to save the CSFG Hamiltonian matrixelements for
! check.
      !IOUTF = 90 + jblock
      !IOUTF = 91
      !WRITE(IOUTF, *)'JBLOCK = ', JBLOCK 
      !IOUTT = 9900 + jblock

! Scratch file stores the IROWs for nonzero elements in sparse
! structure.
      !OPEN (IOUTT, STATUS = 'SCRATCH', FORM = 'UNFORMATTED')

! Now loop over CSFG-Columns of the Hamiltonian matrix - distributed by 
! the same ways of rangular_csfg_mpi.

! For check
      NPAIRtmp(1:2) = 0
!
      ICUT = ICCUT(NB)
      NROWGPAST = NPAIRPAST(2)
      NODENCOLS = 0

! Record how many times returning non-zero results by ONESCALAR and RKCO
      NONESCALAR = 0
      NRKCO = 0

! Total number of non-zero matrixelements of JBLOCK
      NELMNTCC = 0
!
! Allocate storage for the arrays in BUFFER
!
      CALL ALCBUF (1)
      
      IF (.NOT. LREADMCP) THEN
! Maximum size of buffer is obtained as being MAXNV, then no requrement
! to call ALCBUF again.
        CALL RALLOC(LABEL, 5, MAXNV, 'LABEL', 'SETHAM_POT')
        CALL RALLOC(COEFF, MAXNV, 'COEFF', 'SETHAM_POT')        
      ENDIF 

! Sparse structure of H-matrix, initializaiton
      IENDW = 0
      ICTOT = 0 ! Column counter in the fully expanded CSF
      DO 10 ICOLG = 1, NCOLG
        ! For check
        NPAIRtmp(1) = NPAIRtmp(1) + 1
        ! ICC is the column index in CSFG list, ICC is the BLOCK's
        ! index, not in the whole CSFG list file (rcsfg.inp)
        ICC = ICGPDET(NPAIRPAST(1) + ICOLG)
        ! column number in the CSFG list, when only rcsfg.inp is used.
        IC = ICC
        ! column number in the expanded CSFs file: rcsf.inp (block),
        ! not used again
        !IC  = MAP(ICC)
        
        ! ICTOT is used to map the .c list and H-matrix
        IF (ICC.GT.1) ICTOT = MAP1(ICC-1, JBLOCK)
        ! Full matrixelements of the diagonal BLOCK generated within
        ! the IC-th CSFG with IC <= ICCUT should be all kept in
        ! the zero-first calculation.
        IF (ICC.LE.ICUT) THEN
          LICCUT = .FALSE.
        ELSE
          LICCUT = .TRUE.
        ENDIF

        ! Total number of non-zero matrixelements of ICC-th CSFG
        NELC = 0 
        ! Number of non-zeros of every column spanned by ICC CSFG
        NELCSYM(:) = 0 

! Loop over IROWG rows for the current ICOLG column
        ! Row counter in the fully expanded CSF and H-matrix 
        IRTOT = 0

        IF (ICOLG == 1) THEN
          NROWG = ICGPEND(NPAIRPAST(1)+ICOLG)
        ELSE
          NROWG = ICGPEND(NPAIRPAST(1)+ICOLG) -  &
                  ICGPEND(NPAIRPAST(1)+ICOLG - 1)
        ENDIF
        DO 85 IROWG = 1, NROWG
          ! For check
          NPAIRtmp(2) = NPAIRtmp(2) + 1

          ! row index in CSFG list (block) 
          IRR = IRGPDET(NROWGPAST + IROWG)

          ! Use only CSFG list
          IR = IRR
          !row number in the expanded CSFs file: rcsf.inp, not used
          !anymore
          !IR = MAP(IRR)

          IF (IRR.GT.1) IRTOT = MAP1(IRR-1, JBLOCK)

! Zero-First calculation
          IF (LFORDR .AND. (IRR .GT. ICUT)) THEN
             IF (IRR.NE.ICC) CYCLE
          END IF

! Obtain the the generalized weights for the sub-block spanned by the 
! ICC-IRR CSFGs pair.They are calculated repeatedly by DSUBRS in
! GRASP2018.
          IF (LCPOT) CALL CALDDRS(EOL, JBLOCK, ICC, IRR)
!
! Now, calculate the sub-block of H-matrix generated by the ICC-IRR CSFG
! pair. Only CSFG list is used, IC = ICC and IR = IRR. 
!          
          CALL SETHAM_CYC(ICC, IRR, IC, IR)

! Total number of V-K-Coefficients generated by one ICG-IRG pair, the
! demension is SUM(NPAIRS2,1:NBLOCK)), these VK terms use the same
! drs, i.e. the generalized weights.
! NTPV is set by savep_vk
        !IF (LABTVFRST) &
        !  NVKPG(NPAIRPAST(2) + NPAIRtmp(2)) = NTPV

   85   CONTINUE ! End IR/IRR loop

        NROWGPAST = NPAIRPAST(2) + ICGPEND(NPAIRPAST(1)+ICOLG)
        IF (LCHM .AND. NELCSYM(NTYPE(2, ICC+NCFGPAST)).gt.0) then
          DO I = 1, NTYPE(2, ICC+NCFGPAST)
            !WRITE (IOUTT)                                              &
            !       NELCSYM(I), (EMTSYM(J,I),  J = 1, NELCSYM(I)),      &
            !                   (IROWSYM(J,I), J = 1, NELCSYM(I))

            EMT(IENDW+1:IENDW+NELCSYM(I)) =  EMTSYM(1:NELCSYM(I), I)
            IROW(IENDW+1:IENDW+NELCSYM(I))= IROWSYM(1:NELCSYM(I), I)
            IENDW = IENDW + NELCSYM(I)

            NODENCOLS = NODENCOLS + 1
            NODECOLS(NODENCOLS) = ICTOT + I
            IENDC(ICTOT + I) = IENDW
            NELMNTCC = NELMNTCC + NELCSYM(I)
            NELC = NELC + NELCSYM(I)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output H-Matrix for check
            !write(IOUTF,'(2i12,5000(1pe22.14))') &
            !      ICTOT+I, 0        ,  (EMTSYM(J,I),  J = 1,NELCSYM(I))
            !write(IOUTF,'(2i12,5000i22)') &
            !      ICTOT+I, NELCSYM(I), (IROWSYM(J,I), J = 1, NELCSYM(I))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ENDDO 
        ENDIF 
   10 CONTINUE ! End IC/ICC loop

! All the column CSFs indexes of MYID process
      IF (LREADMCP) ICOLIDX(NCOLPAST+1:NCOLPAST+NCOLS) =   &
                     NODECOLS(1:NODENCOLS)
! Check
      IF (LCHM .AND. NODENCOLS .NE. NCOLS) THEN
        WRITE(*,*)'MYID, NB, NCOLS, NODENCOLS =', &
                   MYID, NB, NCOLS, NODENCOLS
        STOP "ERROR! NODENCOLS .NE. NCOLS ..."
      ENDIF
! if symmatrix_mod :: CUTOFF=1.0D-10, the number of non-zero elements 
! counted by rangular_csfg_mpi will be larger than NELMNTCC.
      IF (LCHM .AND. (NELMNTCC.GT.NELMNT .OR. IENDW.NE.NELMNTCC)) THEN
        WRITE(*,*)'MYID, NB, NELMNTCC, IENDW, NELMNT =', &
                   MYID, NB, NELMNTCC, IENDW, NELMNT
        STOP "ERROR! NELMNTCC .GT. NELMNT ..." 
      ENDIF
!
! Block ending lines of csfg_mcpXXX.31 and csfg_mcpXXX.32: "0 0 0"
!
      IF (LREADMCP) THEN  
        READ(IGMCP31)I, J, K
        IF (I /= 0 .OR. J /= 0 .OR.  K /= 0) &
          STOP "Sth error in reading csfg_mcpXXX.31 ..."
        READ(IGMCP32)I, J, K
        IF (I /= 0 .OR. J /= 0 .OR.  K /= 0) &
          STOP "Sth error in reading csfg_mcpXXX.32 ..."
      ENDIF

! Record potential data for every Jp block
      IF (LABTVFRST) THEN
        NBPT(NB)   = NTPT  - SUM(NBPT(0:NB-1))
        NBCRPT(NB) = NCRPT - SUM(NBCRPT(0:NB-1))

        NBPV(NB)   = NTPV  - SUM(NBPV(0:NB-1))
        NBCRPV(NB) = NCRPV - SUM(NBCRPV(0:NB-1))
      ENDIF
!
! Check
      IF (NPAIRtmp(1) /= NPAIRS(1,NB) .OR. NPAIRtmp(2) /= NPAIRS(2,NB)) THEN
        WRITE(*,*)"NPAIRtmp(1), NPAIRS(1,NB) =",NPAIRtmp(1),NPAIRS(1,NB)
        WRITE(*,*)"NPAIRtmp(2), NPAIRS(2,NB) =",NPAIRtmp(1),NPAIRS(2,NB)
        STOP "NPAIRtmp(1) /= NPAIRS(1,NB) &
              .OR. NPAIRtmp(2) /= NPAIRS(2,NB) "
      ENDIF

! Release LABEL,COEFF?
      CALL ALCBUF (3)

!     ...Locals
!  Deallocate arrays related to the CSFG part
      DEALLOCATE(EMTBLOCK)
!      DEALLOCATE(MAP)

! Deallocated in matrix_csfg.f90
!      DEALLOCATE(NODECOLS)
      DEALLOCATE(NELCSYM)
      DEALLOCATE(EMTSYM)
      DEALLOCATE(IROWSYM)

! Close the scratch file
!      CLOSE(UNIT=IOUTF)
      !WRITE(IOUTF,*)
!
! Merge NXAIORB to obtain NXAOPT, and reset JXIPOS values
!
      IF (LABTVFRST .AND. JBLOCK == NBLOCK) THEN
        DO J = 1, NW
          DO IORB = 1, NW
            JXIPOS(IORB,J) = JXIPOS(IORB-1,J) + JXIPOS(IORB,J)
          ENDDO
          NXCOFOPT(J) = JXIPOS(NW,J)
        ENDDO

        NXCOF = MAXVAL(NXCOFOPT(1:NW))
        IF (NXCOF > SIZE(NXAOPT,DIM=1)) THEN
          CALL RALLOC (NXAOPT, NXCOF, NW, 'NXAOPT', 'SETCOF_V_CSFG')
        ENDIF

        DO J = 1, NW
          DO IORB = 1, NW
             I = JXIPOS(IORB,J) - JXIPOS(IORB-1,J)
             NXAOPT(JXIPOS(IORB-1,J)+1 : JXIPOS(IORB,J), J) = &
                NXAIORB(1:I, IORB, J)
          ENDDO
        ENDDO
      ENDIF

!***********************************************************************

      RETURN
      END SUBROUTINE SETHAM_POT
