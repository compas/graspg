!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIX (ncore, j2max)
!                                                                      *
!   This SUBROUTINE calls routines to  form the  Hamiltonian  matrix   *
!   and to diagonalise it.   The total angular momenta and parity of   *
!   the ASF are found;  the eigenvectors are chosen so that the sign   *
!   of the largest element is positive.                                *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, HMOUT, IQ              *
!                        WGHTD5.                                       *
!               [RCI92]: MANEIG, QED, SETHAM.                          *
!                                                                      *
!                                         Last revision: 28 Dec 1992   *
!   Modified by Xinghong He               Last revision: 31 Dec 1997   *
!   Block version Xinghong He             Last revision: 12 Jun 1998   *
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
      USE vast_kind_param, ONLY: DOUBLE, LONG
      USE parameter_def,   ONLY: NNNW
      USE memory_man
      USE decide_C
      USE DEBUG_C
      USE def_C
      USE eigv_C
      USE hblock_C
      USE hmat_C
      USE iounit_C
      USE jlabl_C, LABJ=>JLBR, LABP=>JLBP
      USE orb_C, ONLY: NCF, nw, iqa
      USE prnt_C
      USE stat_C
      USE wave_C
!CYC      USE where_C
      USE blim_C
      USE eigvec1_C
      USE iccu_C,   ONLY: iccutblk
      USE foparm_C
      USE cteilsrk_C
      USE coeils_C
      USE bilst_C
      USE vpilst_C
      USE keilst_C
      USE vinlst_C
      USE fposition_C
      USE symexpand_mod
      use symmatrix_mod, ONLY: NCFTr, IMCDF, NCFBLKTr
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE auxblk_I
!CYC      USE lodcslmpi_I
      USE lodcslmpi_csfg_I
      USE qed_slfen_I
      USE genmat_I
      USE genmat2_I
      USE hmout_I
      USE iq_I
      USE maneig_I
      USE engout_I
      USE wghtd5_I
      USE qed_I
      USE mpi_C
      USE itjpo_I
      USE ispar_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: ncore
      INTEGER,  INTENT(IN):: j2max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW) :: slfint
      CHARACTER(LEN=8) :: CNUM
      REAL(DOUBLE) :: atwinv, elsto, eau, ecm, eev, elemnt
      REAL(DOUBLE), DIMENSION(:), pointer :: slf_en, ucf, etot
      INTEGER(LONG) :: nelmnt_a
      INTEGER :: iiatjpo, iiaspar
      INTEGER :: i, j, irestart, ncminpas, jblock, ic, ip, mode
      INTEGER :: iflagda
      INTEGER :: ipar4, jblock4, jpo4, NCF4, NVEC4
      REAL(DOUBLE)  :: endtime, starttime
!-----------------------------------------------
!     ...Common to all blocks - place here to save CPU time
      CALL auxblk (j2max, atwinv)

!***************************************************************
!      Loop over blocks
!***************************************************************
      ncminpas = 0
      DO 100 jblock = 1, nblock
         NCF    =   NCFblk(jblock)
         NCFTr  =   NCFBLKTr(JBLOCK)
         nvec   =   nevblk(jblock)
         nvecmx = ncmaxblk(jblock)
         iccut  = iccutblk(jblock)
         !.. Determine position of the previous block in the .res file
!
! CYC: By now, the RESTART mode has not yet been implemented. 
!

         nposition = 7 + nw + nw  ! File position of the previous block
                                  ! in the .res file
! CYC: Only H matrixelements for the present BLOCK are saved. 
!         DO i = 1, jblock - 1
!            j = (NCFblk(i) - myid - 1 + nprocs) / nprocs
!            IF (NCFblk(i) .LT. nprocs) j = NCFblk(i) / (myid+1)
!            nposition = nposition + j + 1
!         ENDDO

         !.. SETHAM does not handle this extrem case
         IF (nprocs .GT. NCF)                                          &
                         CALL stopmpi ('matrix: too many nodes', myid)

!        ...Obtain ivec() from iccmin()
         IF (nvec .GT. 0) THEN
            CALL alloc (ivec, nvec, 'IVEC', 'MATRIX')
            DO i = 1, nvec
               ivec(i) = iccmin(i+ncminpas)
            ENDDO
            ncminpas = ncminpas + nvec
         ENDIF

!        ...These 3 were allocated in lodcsh2 and deallocated at the end
!        ... of this routine and in the setham. In this block version,
!        ... both allocation and deallocation are placed here. See the
!        ... following goto 80 for reason.

! For CSFG-version, TEN FICTIOUS CSFs are needed. 
! They are constructed as needed just before ONESCALAR and RKCO_GG calls.
         CALL ALLOC (IQA, NNNW, NCF+10, 'IQA', 'MATRIX')
         CALL ALLOC (JQSA, NNNW, 3, NCF+10, 'JQSA', 'MATRIX')
         CALL ALLOC (JCUPA, NNNW, NCF+10, 'JCUPA', 'MATRIX')
         CALL ALLOC (SLF_EN,NCF,'SLF_EN', 'MATRIX')
         CALL ALLOC (UCF, nw,'UCF', 'MATRIX')
         iflagda = 0
          do ic=1,NCF
               SLF_EN(IC) = 0.0
          enddo

!      ...Load CSF list of the current block
        CALL lodcslmpi_csfg (21, ncore, jblock)
         iiatjpo = ITJPO (1)
         iiaspar = ISPAR (1)
         IF (LSE) THEN
            IF (myid .EQ. 0) THEN
              PRINT *, 'Entering QED ...'

              WRITE (24,*)
              WRITE (24,*) ' Self Energy Corrections: '
              WRITE (24,*)
            endif
            CALL QED_SLFEN (SLFINT)
            DO IC = 1, NCF
               ELEMNT = 0.0D00
               !DO I = 1,NW
               ! Include only the contributions from labelling orbitals.
               DO I = 1, nonsym
                  ELEMNT = ELEMNT+IQ (I,IC)*SLFINT(I)
               ENDDO
               SLF_EN(IC) = ELEMNT
            ENDDO

            IF (myid .EQ. 0) THEN
               WRITE (24,*)
               WRITE (24,*) 'Self-energy corrections estimated'       &
                       //' --- these will influence the data'
               WRITE (24,*) ' in the RCI92 MIXing coefficients File.'
               endif
            ENDIF

! zou
         IF (nvec <= 0) THEN
            NVEC4 = nvec
            !NCF4 = NCF
            NCF4 = NCFTr
            jblock4 = jblock
            IF (myid .EQ. 0) WRITE (25) jblock4, NCF4, NVEC4, 999, 999
!           ...Generate H matrix of the current block and write to disk
!           ...eav, nelmnt are also obtained from genmat
            CALL genmat (atwinv, jblock, myid, nprocs, elsto, irestart,&
                         slf_en)
               ! This call is optional, meaning the matrix of this block
               ! does not need to be computed. But don't comment it out
               ! since other programs may assume thet existence of them.
            CALL genmat2 (irestart, nelmnt_a, elsto)
            GOTO 80 ! need to clear memory
         ENDIF
!        ------------------------
         CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
         starttime = MPI_Wtime()
         !if (myid == 0) print *, "Entering genmat ...."
         CALL genmat (atwinv, jblock, myid, nprocs, elsto, irestart,   &
                      slf_en)
         !if (myid == 0) print *, "Exiting  genmat ...."
         CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
         endtime = MPI_Wtime()
         if (myid.eq.0)                                                &
           print *, 'Time for calculation of H matrix (m) = ',         &
                  (endtime-starttime)/60.0
         !if (myid.eq.0) write(*,*)"Now calling genmat2 ..."
         CALL genmat2 (irestart, nelmnt_a, elsto)
         !if (myid.eq.0) write(*,*)"Now existing from genmat2 ..."
!
!   Allocate and deallocate memory for the mixing coefficients
!   from the prerun
!
      !IF (IPRERUN.EQ.1) CALL ALLOC (EVEC1,NCF*NVEC, 'EVEC1', 'MATRIX')
      IF (IPRERUN.EQ.1) CALL ALLOC (EVEC1,NCFTr*NVEC, &
                                    'EVEC1', 'MATRIX')
      IF (IPRERUN.EQ.2) CALL DALLOC (EVEC1, 'EVEC1', 'MATRIX')

! Since maneig needs both nelmnt and nelmnt_a, # of non-zero
! matrix elements computed by current node and the total #
! from all nodes, a new var nelmnt_a is created rathan than
! using the one in the common block.

!CYC: Deallocate arrays IQA, JQSA, JCUPA, SLF_EN before maneig.
        if (iflagda.eq.0) then
!GG           CALL dalloc (PNTRIQ)
!GG           CALL dalloc (PNTJQS)
!GG           CALL dalloc (PNJCUP)
!GG           CALL dalloc (PNTSLF)
           CALL dalloc (IQA, 'IQA', 'MATRIX')
           CALL dalloc (JQSA, 'JQSA', 'MATRIX')
           CALL dalloc (JCUPA, 'JCUPA', 'MATRIX')
           CALL dalloc (SLF_EN, 'SLF_EN', 'MATRIX')
           !CALL dalloc (PNTUCF)
           iflagda = 1
        endif
        if (myid .eq. 0) write (*,*)
        if (myid .eq. 0) write (*,*) "Calling maneig ..."
        starttime = MPI_Wtime()

        CALL MANEIG (iiatjpo, iiaspar, nelmnt_a)
!        if (myid .eq. 0) write (*,*) 'Exit from maneig'
        CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

        endtime = MPI_Wtime()
        if (myid.eq.0)                                                &
          print *, &
         'Time for maneig (reading and diagonalizing matrix) (m) = ', &
         (endtime-starttime)/60.0
!
!   Write out eigenvalues (ENGOUT), dominant components of the
!   eigenvectors (WGHTD5) to stream 24; write out ASF symmetries,
!   eigenvalues  and eigenvectors to RCI92 mixing coefficients file.
!   EAV and ELSTO are added back to energy here
!
        IF (myid .EQ. 0) THEN
! ELSTO has never been in Hamiltonian matrix, yet it was
! added to EAV which was later substracted from H. Thus at
! this point, EAV is correct (it has ELSTO added), EVAL()
! need ELSTO and the correct EAV.
!           IF (MYID.EQ.0) WRITE(*,*)"Matrix: ELSTO = ", ELSTO
           IF (NCF > 1) then
              DO i = 1, NVEC
                 EVAL(i) = EVAL(i) + ELSTO
              ENDDO
           END IF

           CALL ENGOUT (EAV,EVAL,IiATJPO,iIASPAR,IVEC,NVEC,3)
           CALL WGHTD5 (iiatjpo, iiaspar)

!          ...Write ASF symmetries, eigenvalues, and eigenvectors to RCI92
!          ...MIXing coefficients File; close the file; print a report
           jblock4 = jblock
           !NCF4 = NCF
           NCF4 = NCFTr
           nvec4 = nvec
           jpo4 = iiatjpo
           ipar4 = iiaspar
           WRITE (25) jblock4, NCF4, nvec4, jpo4, ipar4
           WRITE (25) (ivec(i), i = 1, nvec)
           WRITE (25) EAV,(EVAL(I),I = 1,NVEC)
           WRITE (25) ((EVEC(NCFTr*(J-1)+I),I=1,NCFTr),J=1,NVEC)

           PRINT *, 'RCI90 MIXing coefficients File generated.'
           PRINT *
        ENDIF
!
!   Save the mixing coefficients from the prerun
!
      IF (IPRERUN .EQ. 1) THEN
         DO J = 1, NVEC
            DO I = 1, NCFTr
               EVEC1(I+(J-1)*NCFTr) = EVEC(I+(J-1)*NCFTr)
            ENDDO
         ENDDO
      ENDIF
!
!   Estimate diagonal self-energy corrections; the stored
!   eigenvalues and eigenvectors are not modified by these
!   estimates
!
!      IF (.not.LSE) THEN
!         IF (myid .EQ. 0)            PRINT *, 'Entering QED ...'
!            CALL ALLOC (ETOT,NVEC,'ETOT', 'MATRIX')
!         IF (myid .EQ. 0) THEN
!            WRITE (24,*)
!            WRITE (24,*) ' Self Energy Corrections: '
!            WRITE (24,*)
!            WRITE (24,301)
!            WRITE (24,*)
!         END IF
!  301 FORMAT (' Level  J Parity',7X,'Hartrees',14X,'Kaysers',          &
!               16X,'eV' )
!  302 FORMAT (1I3,2X,2A4,1P,3D22.14)
!         DO J = 1, nvec
!            CALL QED (j,SLFINT,UCF)
!            ELEMNT = 0.0D00
!            IC = IVEC(J)
!            DO I = 1,NW
!               ELEMNT = ELEMNT+UCF(I)*SLFINT(I)
!              ELEMNT = ELEMNT+IQ (I,IC)*SLFINT(I)
!            ENDDO
!            ETOT(J) = EVAL(J)+ELEMNT
!
!            EAU = ELEMNT
!            ECM = EAU*AUCM
!            EEV = EAU*AUEV
!            IP = (IIASPAR+3)/2
!            IF (myid .EQ. 0)                                           &
!               WRITE (24,302) j,LABJ(IiATJPO),LABP(IP),EAU,ECM,EEV
!
!         ENDDO
!         IF (myid .EQ. 0) THEN
!            WRITE (24,*)
!            WRITE (24,*) 'Self-energy corrections estimated'           &
!                    //' --- these do not influence the data'
!            WRITE (24,*) ' in the RCI92 MIXing coefficients File.'
!            CALL ENGOUT (EAV,ETOT,IiATJPO,iIASPAR,IVEC,NVEC,MODE)
! zou       CALL ENGOUT (EAV+elsto,ETOT,IiATJPO,iIASPAR,IVEC,NVEC,MODE)
!         ENDIF
!         CALL dalloc (ETOT, 'ETOT', 'MATRIX')
!      ENDIF
!        ...Locals
         CALL dalloc (ivec, 'IVEC', 'MATRIX')
!        ...Allocated in maneig
         CALL dalloc (eval, 'EVAL', 'MATRIX')
         CALL dalloc (evec, 'EVEC', 'MATRIX')

   80    CONTINUE

!        ...Locals
        if (iflagda.eq.0) then
            CALL dalloc (IQA, 'IQA', 'MATRIX')
            CALL dalloc (JQSA, 'JQSA', 'MATRIX')
            CALL dalloc (JCUPA, 'JCUPA', 'MATRIX')
            CALL dalloc (SLF_EN, 'SLF_EN', 'MATRIX')
            CALL dalloc (UCF, 'UCF', 'MATRIX')
        else
            CALL dalloc (UCF, 'UCF', 'MATRIX')
        endif

  100 CONTINUE
!
!   Close the restart files; nothing will be added to them now
!
      CLOSE (imcdf)
!     ...Clean up
      CALL dalloc (NCFblk, 'NCFBLK', 'MATRIX')
      CALL dalloc (nevblk, 'NEVBLK', 'MATRIX')
      CALL dalloc (ncmaxblk, 'NCMAXBLK', 'MATRIX')
      CALL dalloc (iccutblk, 'ICUTTBLK', 'MATRIX')
      CALL dalloc (iccmin, 'ICCMIN', 'MATRIX') ! allocated in items as pnccmin

      CALL dalloc (VALTEIRK, 'VALTEIRK', 'MATRIX') ! allocated in genintrk
      CALL dalloc (INDTEIRK, 'INDTEIRK', 'MATRIX') ! allocated in genintrk
!
!   Deallocate storage for the integral lists from the
!   Dirac-Coulomb operator; the storage was allocated
!   in IABINT and RKINTC
!
      IF (NCOEI .GT. 0) THEN
         CALL DALLOC (INDOEI, 'INDOEI', 'MATRIX')
         CALL DALLOC (VALOEI, 'VALOEI', 'MATRIX')
      ENDIF
!
!   Deallocate storage for the integral lists from the
!   transverse photon interaction operator; this storage
!   was allocated in BRINT1, brint2,...brint6
!
      IF (LTRANS) THEN
         IF (NTPI(1) .GT. 0) THEN
            CALL DALLOC (INDTP1, 'INDTP1', 'MATRIX')
            CALL DALLOC (VALTP1, 'VALTP1', 'MATRIX')
         ENDIF
         IF (NTPI(2) .GT. 0) THEN
            CALL DALLOC (INDTP2, 'INDTP2', 'MATRIX')
            CALL DALLOC (VALTP2, 'VALTP2', 'MATRIX')
         ENDIF
         IF (NTPI(3) .GT. 0) THEN
            CALL DALLOC (INDTP3, 'INDTP3', 'MATRIX')
            CALL DALLOC (VALTP3, 'VALTP3', 'MATRIX')
         ENDIF
         IF (NTPI(4) .GT. 0) THEN
            CALL DALLOC (INDTP4, 'INDTP4', 'MATRIX')
            CALL DALLOC (VALTP4, 'VALTP4', 'MATRIX')
         ENDIF
         IF (NTPI(5) .GT. 0) THEN
            CALL DALLOC (INDTP5, 'INDTP5', 'MATRIX')
            CALL DALLOC (VALTP5, 'VALTP5', 'MATRIX')
         ENDIF
         IF (NTPI(6) .GT. 0) THEN
            CALL DALLOC (INDTP6, 'INDTP6', 'MATRIX')
            CALL DALLOC (VALTP6, 'VALTP6', 'MATRIX')
         ENDIF
      ENDIF
!
!   Deallocate storage for the nuclear motional energy integral
!   lists; this was allocated in KEINT and VINT
!
      IF (LNMS) THEN
         IF (NKEI .GT. 0) THEN
            CALL DALLOC (INDKEI, 'INDKEI', 'MATRIX')
            CALL DALLOC (VALKEI, 'VALKEI', 'MATRIX')
         ENDIF
      ENDIF
      IF (LSMS) THEN
         IF (NVINTI .GT. 0) THEN
            CALL DALLOC (INDTEI, 'INDTEI', 'MATRIX')
            CALL DALLOC (VALTEI, 'VALTIE', 'MATRIX')
         ENDIF
      ENDIF
!
!   Deallocate storage for the vacuum-polarisation integral list;
!   this was allocated in VPINT
!
      IF (LVP) THEN
         IF (NVPI .GT. 0) THEN
            CALL DALLOC (INDVPI, 'INDVPI', 'MATRIX')
            CALL DALLOC (VALVPI, 'VALVPI', 'MATRIX')
         ENDIF
      ENDIF

      CALL dalloc (PF, 'PF', 'MATRIX') ! lodrwf or lodres
      CALL dalloc (QF, 'QF', 'MATRIX') ! lodrwf or lodres

      RETURN

      END SUBROUTINE MATRIX
