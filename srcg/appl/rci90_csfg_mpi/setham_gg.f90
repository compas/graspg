!***********************************************************************
!                                                                      *
      SUBROUTINE SETHAM (myiddumm, nprocsdumm, jblock, ELSTO,ICSTRT,   &
                                               nelmntt, atwinv,slf_en)
!                                                                      *
!   Sets up the Hamiltonian matrix and determines the average energy.  *
!
!   Serial I/O moved out; able to run on single/multiple processors
!   For this purpose a new common /setham_to_genmat2/ is created
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, CONVRT, DALLOC, ICHOP, RKCO, TNSRJJ.  *
!               [RCI92]: BRINT, IABINT, KEINT, RKINTC, VINT, VPINT.    *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 30 Oct 1992   *
!   Block version by Xinghong He          Last revision: 15 Jun 1998   *
!   Modify  by G Gaigalas                         May 2021             *
!                                                                      *
!                                                                      *
!   Last modification for CSFG by Chongyang Chen  Dec 2023             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE, LONG, BYTE
      USE parameter_def,   ONLY: NNNP, KEYORB
      USE memory_man
!cjb iso_fortran_env :: output_unit, error_unit
      use, intrinsic :: iso_fortran_env
      use symexpand_mod
      use symmatrix_mod, CUTOFF0=>CUTOFF
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      USE bcore_C, ONLY: icore
      USE bilst_C
      USE buffer_C
      USE coeils_C
      USE decide_C
      USE debug_C
      USE def_C
      USE eigv_C
!GG      USE iccu_C
      USE foparm_C
      USE grid_C
      USE hmat_C
      USE keilst_C
      USE ncdist_C
      USE orb_C, ONLY: NCF, NW, IQA
      USE prnt_C
      USE stat_C
      USE stor_C
      USE tatb_C
      USE vinlst_C
      USE vpilst_C
      USE wave_C
      USE cteilsrk_C
      USE blim_c
!CYC      USE where_C
      USE eigvec1_C
!      USE mpi_C,   ONLY: ierr, MPI_COMM_WORLD, MPI_INTEGER
      USE mpi_C
!      USE MPI

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I
      USE iabint_I
      USE brint1_I
      USE brint2_I
      USE brint3_I
      USE brint4_I
      USE brint5_I
      USE brint6_I
      USE convrt_I
      USE ichop_I
      USE keint_I
      USE rkintc_I
      USE vint_I
      USE vpint_I
      IMPLICIT NONE
      EXTERNAL BREID,CORD
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER(LONG)       :: NELMNTT
      INTEGER             :: JBLOCK, ICSTRT
      INTEGER, INTENT(IN) :: MYIDDUMM, NPROCSDUMM
      REAL(DOUBLE)        :: ELSTO, ATWINV, endtime, starttime
      REAL(DOUBLE), DIMENSION(*) :: slf_en
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
      INTEGER :: ipi, ipj, inc1, inc2, kt, ipt, incor, ncoec, nctec, &
                 i, j, nmcbp, ncore, ic, nelc, irstart, ir, ia, ib,  &
                 itype, nctei, iia, ICC, ioutf, IRR, NCFGEN, NCFSMALL, &
                 NCFTOT, IC0, IR0
      REAL(DOUBLE) :: elemnt, precoeff, tcoeff, vcoeff, contr

      INTEGER(BYTE) :: IQAICC(NW), IQAIRR(NW), IQADiff(NW)
!-----------------------------------------------------------------------
      nelmnt = nelmntt

! Call findtype to initialize the values in symmatrix_C.mod
      NCFSMALL=TotCSFs_perblock(2,jblock)
!Check: 
      IF (NCF.NE.NCFSMALL) THEN
        IF (myid .eq. 0) THEN
          WRITE(*,*)"Error, Unexpected NCF.NE.NCFSMALL within setham_gg"
        ENDIF
        STOP "Error, Unexpected NCF.NE.NCFSMALL within setham_gg ..."
      ENDIF

      ALLOCATE(MAP(NCFSMALL))
!GG      CALL alloc (MAP, NCFSMALL, 'MAP', 'SETHAM')
      DO I = 1,NCFSMALL
         MAP(I)=map1(i,jblock)
      END DO

!  Call findtype in order to define the type of CSF and largest value
!  for the orbitals of the CSFG. Note that we define the type
!  for all CSFs in the added list although we only will need the
!  information for the CSFs defined by the mapping. [Old version, TWO
!  lists were used throughout the codes.]

!  New version, the smallest list is used throughout the codes to save
!  memory.
!
      ALLOCATE(NTYPE(6,NCF))
!GG      CALL alloc (NTYPE,NCF, 'NTYPE','SETHAM')
      if (MYIDDUMM.eq.0) then
         CALL FINDTYPE(NCFSMALL,NCFGEN,NCFTOT)
         if (NCFGEN.ne.TotCSFs_perblock(1,jblock).or.                &
             NCFTOT.ne.TotCSFs_perblock(3,jblock).or.                &
             NCFTr.NE.NCFTOT) then
           write(*,*)'NCFGEN=',NCFGEN,                               &
             ' TotCSFs_perblock(1,jblock)=',TotCSFs_perblock(1,jblock)
           write(*,*)'NCFTOT=',NCFTOT,                               &
             ' TotCSFs_perblock(3,jblock)=',TotCSFs_perblock(3,jblock)
           write(*,*)'NCFTr=',NCFTr, ' NCFTOT=',NCFTOT
           stop 'Unexpected error in setham_gg ===1=== ...'
         endif
      endif
      !write(*,*)'Now in setham_gg ===2=== '
! Broadcast ...
      CALL MPI_Bcast(ncsfDF1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(MAXSPAN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NORBGEN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NCFTOT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NTYPE(:,:),6*NCF,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
! Allocate working arrays...
      ALLOCATE(IROWSYM(NCFTOT,MAXSPAN))
!GG      CALL DALLOC (IROWSYM,NCFTOT,MAXSPAN, 'IROWSYM', 'SETHAM')
      ALLOCATE(EMTSYM(NCFTOT,MAXSPAN))
!GG      CALL DALLOC (EMTSYM,NCFTOT,MAXSPAN, 'EMTSYM', 'SETHAM')
      ALLOCATE(EMTBLOCK(MAXSPAN,MAXSPAN))
!GG      CALL DALLOC (EMTBLOCK,MAXSPAN,MAXSPAN, 'EMTBLOCK', 'SETHAM')
      ALLOCATE(NELCSYM(MAXSPAN))
!GG      CALL DALLOC (NELCSYM,MAXSPAN, 'NELCSYM', 'SETHAM')

! Node-ID for the IC-th column data
!CC      ALLOCATE(ICNODES(NCFTOT))
! Number of the non-zero matrixelement of the IC-th column
!CC      ALLOCATE(NZMCOL(NCFTOT))
!CC      ICNODES = 0
!CC      NZMCOL = 0
      NODENCOLS = 0
! NODECOLS used in spicmvsym.f, deallocated in maneig.f
      ALLOCATE(NODECOLS(NCFTOT))
!GG      CALL DALLOC (NODECOLS,NCFTOT, 'NODECOLS', 'SETHAM')

      ! Positions of the FICTIOUS CSFs
      IRFICT = NCF + 1
      ICFICT = NCF + 6

! Continue the original setham_gg ...
      nelmnt = nelmntt
      EAV = 0.D0
!
!*************************
!
!     CSF(i) = CSFG that spans N(i) CSFs
!     <CSF   |H|CSF(i)> generates an 1 x N(i) array of elements
!     <CSF(i)|H|CSF(i)> generates an upper triangular N(i) x N(i) array
!     <CSF(i)|H|CSF(j)> generates an N(i) x N(j) array of elements
!
!     Allocate an array of some largest dimension Nmax x Nmax as
!     computed by
!     FINDTYPE. Store the elements in this array
!
!     Aoocate a double precision and an integer array with the largest
!     dimension NCSFtot x Nmax for storing the matrix elements and the
!     row positions.
!     NCSFtot is computed by FINDTYPE. In reality we can take NCSFtot
!     smaller as the matrix will be sparse.
!     10 000 000 x 400 = 32 Gb is very much nore than we need.
!
!     As the computation of the interactions between the CSFGs
!     proceeds fill the arrays above with non-zero elements. We can then
!     write to rci.res in the same format as before
!
!
      ATWINV = 1.D0/EMN
      IPRERUN = 0   ! Pre-run has been disabled
!
!   Allocate storage to arrays in COMMON/BUFFER/; these are
!   used for the Coulomb and transverse two-electron integrals
!
      CALL ALCBUF (1)

!     ...Locals
!      CALL alloc (pntemt, NCF, 8) !!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!
!      CALL alloc (EMT, NCFTr, 'EMT', 'SETHAM')  ! Not used.
!      CALL alloc (pnirow, NCF, 4) !!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!
!      CALL alloc (IROW, NCFTr, 'IROW', 'SETHAM') ! Not used.
!      EMT(1:NCFTr) = 0.D0
!      IROW(1:NCFTr) = 0
!
      INC1 = 1
      INC2 = 1
!
!   Initialisations for contributions from the Dirac-Coulomb
!   operator
!
      KT  = 0
      IPT = 1
      INCOR = 1

      NCOEC = 0
!
      NCTEC   = 0

      IF (LTRANS) THEN

!        ...Initialisations for transverse interaction correction
        DO 2 I = 1, NW
          ICORE(I) = 0
          DO J = 1, NCF ! Here NCF is number of <name>.g. It is OK.
            IF (ICHOP (I,J) .LE. 0) GOTO 2
          ENDDO
          ICORE(I) = 1 ! core-orbital labels
    2   CONTINUE

        NMCBP = 0  ! Count for mcbp calculations, just for Log.
        NCORE = 0
      ENDIF

      IF (LFORDR .AND. MYIDDUMM.EQ.0) WRITE(*,*)'ICCUT=',ICCUT
!
! Open scratch file to save the CSFG Hamiltonian matrixelements.
!
!CC      itmpf=999
      ioutf = 90 + jblock

      !open(itmpf,status='scratch',form='unformatted')
! Loop over rows of the Hamiltonian matrix - distributed
      !starttime = MPI_Wtime()
      ICTOT = 0 ! Column counter in the fully expanded CSF
      !DO 10 IC = icstrt, NCF, NPROCS
      DO 10 ICC = icstrt, NCFSMALL, NPROCSDUMM
        !EMTSYM = 0.D0
        IQAICC(1:NW) = IQA(1:NW, ICC) 
        ! Smallest list used. There is no CSF added.           
        ! IC = MAP(ICC)  ! column number in the expanded CSFs file: name.c
        IC = ICC 
        IC0 = IC ! For check

        ! ICTOT is used to map the .c list and H-matrix
        IF (ICC.GT.1) ICTOT = MAP(ICC-1)

        ! Full matrixelements of the diagonal BLOCK generated within
        ! the IC-th CSFG with IC <= ICCUT should be all kept in
        ! the zero-first calculation.
        IF (ICC.LE.ICCUT) THEN
          LICCUT = .FALSE.
        ELSE
          LICCUT = .TRUE.
        ENDIF

        ! Total number of non-zeros of ICC CSFG
        NELC = 0 
        ! Number of non-zeros of every column generated by ICC CSFG
        NELCSYM = 0 

! Loop over rows for the current column
        ! Row counter in the fully expanded CSF and H-matrix 
        IRTOT = 0
        !DO 85 IR = 1, IC
         DO 85 IRR = 1, ICC
          ! Smallest list used. There is no CSF added.
          !IR = MAP(IRR)  !row number in the expanded CSFs file: name.c
          IR = IRR
          IR0 = IR  ! For check
          !IF (IC.EQ.2756) WRITE(*,*)"IC, IR = ", IC, IR
          IF (IRR.GT.1) IRTOT = MAP(IRR-1)
! zero-first calculation
          IF (LFORDR .AND. (IRR .GT. ICCUT)) THEN
             IF (IRR.NE.ICC) CYCLE
          END IF
!CYC: 2022/12/01
          IQAIRR(1:NW) = IQA(1:NW, IRR)
          IQADiff = IQAICC - IQAIRR
          IF (COUNT(IQADiff /= 0) > 4) CYCLE
!CYC: 2022/12/01
!
          EMTBLOCK = 0.D0     ! accumulates various contributions to H
          LTRANSPOSE = .false.
          LTRANSFER = .false.
          IF (NTYPE(1,IC).EQ.1 .AND. NTYPE(1,IR).EQ.1) THEN
            ! IC = IR (Diagonal matrixelement) is also included here.
            CALL MATRIXBLOCK11(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,   &
                               NCORE,ELSTO)

          ELSEIF ((NTYPE(1,IC).EQ.2).AND.(NTYPE(1,IR).EQ.1)) THEN
            CALL MATRIXBLOCK12(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,   &
                               NCORE,ELSTO)
          ! Zero - First calculation, in some cases IC-CSF is TYPE 1,
          ! but IR-CSF is other TYPE (2, 3, 4, 5)
          ELSEIF ((NTYPE(1,IC).EQ.1).AND.(NTYPE(1,IR).EQ.2)) THEN
            LTRANSPOSE=.TRUE.
            CALL MATRIXBLOCK12(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,   &
                               NCORE,ELSTO)

          ELSEIF ((NTYPE(1,IC).EQ.2).AND.(NTYPE(1,IR).EQ.2)) THEN
             IF (IC.EQ.IR) THEN   ! Diagonal block of type 2
              CALL MATRIXBLOCK2(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,  &
                      NCORE,ELSTO)
             ELSE IF (IC.NE.IR) THEN
              CALL MATRIXBLOCK22(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)
            END IF

          ELSEIF ((NTYPE(1,IC).EQ.3).AND.(NTYPE(1,IR).EQ.1)) THEN
              CALL MATRIXBLOCK13(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)
          ! Zero - First calculation, in some cases IC-CSF is TYPE 1,
          ! but IR-CSF is other TYPE (2, 3, 4, 5)
          ELSEIF ((NTYPE(1,IC).EQ.1).AND.(NTYPE(1,IR).EQ.3)) THEN
              LTRANSPOSE=.TRUE.
              CALL MATRIXBLOCK13(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)

          ELSEIF (((NTYPE(1,IC).EQ.3).AND.(NTYPE(1,IR).EQ.2)) .OR.   &
                  ((NTYPE(1,IC).EQ.2).AND.(NTYPE(1,IR).EQ.3))) THEN
              IF ((NTYPE(1,IC).EQ.2).AND.(NTYPE(1,IR).EQ.3))         &
                 LTRANSPOSE=.TRUE.
              CALL MATRIXBLOCK23(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)

          ELSEIF ((NTYPE(1,IC).EQ.3).AND.(NTYPE(1,IR).EQ.3)) THEN
             IF (IC.EQ.IR) THEN   ! Diagonal block of type 3
              CALL MATRIXBLOCK3(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,  &
                      NCORE,ELSTO)
             ELSE
              CALL MATRIXBLOCK33(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)
             END IF

          ELSEIF ((NTYPE(1,IC).EQ.4).AND.(NTYPE(1,IR).EQ.1)) THEN
              CALL MATRIXBLOCK14(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)
         ! Zero - First calculation, in some cases IC-CSF is TYPE 1,
          ! but IR-CSF is other TYPE (2, 3, 4, 5)
          ELSEIF ((NTYPE(1,IC).EQ.1).AND.(NTYPE(1,IR).EQ.4)) THEN
              LTRANSPOSE=.TRUE.
              CALL MATRIXBLOCK14(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)

          ELSEIF (((NTYPE(1,IC).EQ.4).AND.(NTYPE(1,IR).EQ.2)) .OR.   &
                  ((NTYPE(1,IC).EQ.2).AND.(NTYPE(1,IR).EQ.4))) THEN
              IF ((NTYPE(1,IC).EQ.2).AND.(NTYPE(1,IR).EQ.4))         &
                 LTRANSPOSE=.TRUE.
              CALL MATRIXBLOCK24(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)

          ELSEIF (((NTYPE(1,IC).EQ.4).AND.(NTYPE(1,IR).EQ.3)) .OR.   &
                  ((NTYPE(1,IC).EQ.3).AND.(NTYPE(1,IR).EQ.4))) THEN
              IF ((NTYPE(1,IC).EQ.3).AND.(NTYPE(1,IR).EQ.4))         &
                 LTRANSPOSE=.TRUE.
             CALL MATRIXBLOCK34(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,  &
                      NCORE,ELSTO)

          ELSEIF ((NTYPE(1,IC).EQ.4).AND.(NTYPE(1,IR).EQ.4)) THEN
            IF (IC.EQ.IR) THEN
              CALL MATRIXBLOCK4(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,  &
                      NCORE,ELSTO)
            ELSE
              CALL MATRIXBLOCK44(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)
            ENDIF

          ELSEIF ((NTYPE(1,IC).EQ.5).AND.(NTYPE(1,IR).EQ.1)) THEN
              CALL MATRIXBLOCK15(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)
          ! Zero - First calculation, in some cases IC-CSF is TYPE 1,
          ! but IR-CSF is other TYPE (2, 3, 4, 5)
          ELSEIF ((NTYPE(1,IC).EQ.1).AND.(NTYPE(1,IR).EQ.5)) THEN
              LTRANSPOSE=.TRUE.
              CALL MATRIXBLOCK15(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)

          ELSEIF (((NTYPE(1,IC).EQ.5).AND.(NTYPE(1,IR).EQ.2)) .OR.   &
                  ((NTYPE(1,IC).EQ.2).AND.(NTYPE(1,IR).EQ.5))) THEN
              IF ((NTYPE(1,IC).EQ.2).AND.(NTYPE(1,IR).EQ.5))         &
                 LTRANSPOSE=.TRUE.
              CALL MATRIXBLOCK25(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)

          ELSEIF (((NTYPE(1,IC).EQ.5).AND.(NTYPE(1,IR).EQ.3)) .OR.   &
                  ((NTYPE(1,IC).EQ.3).AND.(NTYPE(1,IR).EQ.5))) THEN
              IF ((NTYPE(1,IC).EQ.3).AND.(NTYPE(1,IR).EQ.5))         &
                 LTRANSPOSE=.TRUE.
              CALL MATRIXBLOCK35(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)

          ELSEIF (((NTYPE(1,IC).EQ.5).AND.(NTYPE(1,IR).EQ.4)) .OR.   &
                  ((NTYPE(1,IC).EQ.4).AND.(NTYPE(1,IR).EQ.5))) THEN
              IF ((NTYPE(1,IC).EQ.4).AND.(NTYPE(1,IR).EQ.5))         &
                 LTRANSPOSE=.TRUE.
              CALL MATRIXBLOCK45(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)

          ELSEIF ((NTYPE(1,IC).EQ.5).AND.(NTYPE(1,IR).EQ.5)) THEN
            IF (IC.EQ.IR) THEN
              CALL MATRIXBLOCK5(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,  &
                      NCORE,ELSTO)
            ELSE
              CALL MATRIXBLOCK55(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP, &
                      NCORE,ELSTO)
            ENDIF
          ENDIF
  85   CONTINUE ! End IR/IRR loop
      IF (IC.NE.IC0 .OR. IR.NE.IR0) THEN
       WRITE(*,*)'IC,IC0,IR,IR0=',IC,IC0,IR,IR0
       STOP "Unexpected IC.NE.IC0 .OR. IR.NE.IR0"
      ENDIF
      IF (NELCSYM(NTYPE(2, IC)).gt.0) then
        DO I = 1, NTYPE(2, IC)
         ! SLF_EN, calculated within matrix.f
         !EMTSYM(NELCSYM(I),I)=EMTSYM(NELCSYM(I),I)+SLF_EN(ICTOT+I)
         ! Smallest list version, NQEDMAX should be not larger than nmaxgen
         EMTSYM(NELCSYM(I),I)=EMTSYM(NELCSYM(I),I)+SLF_EN(IC)
         WRITE (imcdf)                                               &
                 NELCSYM(I), ELSTO, (EMTSYM(IR,I),IR=1,NELCSYM(I)),  &
                                    (IROWSYM(IR,I), IR = 1, NELCSYM(I))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output H-Matrix for check
!         write(ioutf,'(2i12,5000(1pe12.4))') &
!      ICTOT+I, 0         , (EMTSYM(IR,I),  IR = 1,NELCSYM(I))
!         write(ioutf,'(2i12,5000i12)') &
!      ICTOT+I, NELCSYM(I), (IROWSYM(IR,I), IR = 1, NELCSYM(I))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         NODENCOLS = NODENCOLS + 1
         NODECOLS(NODENCOLS) = ICTOT+I

! This EAV (and the above EMTSYM) does not have ELSTO.
         EAV = EAV + EMTSYM(NELCSYM(I),I)
! Update the counter of total non-zero matrixelement.
         NELMNT = NELMNT + NELCSYM(I)
         NELC = NELCSYM(I) + NELC
        ENDDO
      ENDIF
      IF (MOD (ICC, 100) .EQ. 0 ) THEN
      !IF (MOD (ICC, 1000) .EQ. NPROCS ) THEN
        PRINT *, 'NCSF(G)s:', NCFSMALL, ' ICSF(G):', ICC,        &
             ' ICSF:', ICTOT + NTYPE(2,IC), ': ',                &
              NELC, 'nonzero elements;', ' block = ', jblock,    &
             ' MYID =',MYID
      ENDIF

   10 CONTINUE ! End IC/ICC loop
!
!   Deallocate storage for the arrays in /BUFFER/
!
      CALL ALCBUF (3)

!     ...Locals
      !CALL DALLOC (EMT, 'EMT', 'SETHAM')
      !CALL DALLOC (IROW, 'IROW', 'SETHAM')
!  Deallocate arrays related to the CSFG part

! Deallocate in maneig.f
!GG      CALL DALLOC(IROWSYM, 'IROWSYM', 'SETHAM')
!GG      CALL DALLOC(EMTSYM, 'EMTSYM', 'SETHAM')
!GG      CALL DALLOC(EMTBLOCK, 'EMTBLOCK', 'SETHAM')
!GG      CALL DALLOC(MAP, 'MAP', 'SETHAM')
!GG      CALL DALLOC(NELCSYM, 'NELCSYM', 'SETHAM')
      DEALLOCATE(NTYPE)
      DEALLOCATE(IROWSYM)
      DEALLOCATE(EMTSYM)
      DEALLOCATE(EMTBLOCK)
      DEALLOCATE(MAP)
      DEALLOCATE(NELCSYM)

!  Fill the common block /setham_to_genmat2/ for use in genmat2

      CUTOFFtmp = CUTOFF0
      NCOEItmp = NCOEI
      NCOECtmp = NCOEC
      NCTEItmp = NCTEI
      NCTECtmp = NCTEC
      NTPItmp = NTPI
      NMCBPtmp = NMCBP
      NCOREtmp = NCORE
      NVPItmp = NVPI
      NKEItmp = NKEI
      NVINTItmp = NVINTI
      NELMNTtmp = NELMNT
      !NCFtmp = NCF
      NCFtmp = NCFTr

! Ensure the H matrix obtained in all nodes
      !CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
      !flush output_unit
      !endtime = MPI_Wtime()
      !if (myiddumm.eq.0)                                             &
      !  print *, 'Time for calculation of H matrix (m) = ',          &
      !         (endtime-starttime)/60.0
      !flush output_unit
!      write(*,*)'MYID, CSFG NELMNT=', MYIDDUMM, NELMNT
      flush output_unit
!***********************************************************************
      !WRITE(*,*)"NTPI=",NTPITMP
      !WRITE(*,*)"NDTPA=",NDTPA
!      IF (MYIDDUMM.EQ.0) WRITE(*,*)"Now, exiting setham_gg ..."

      RETURN
      END SUBROUTINE SETHAM
