!***********************************************************************
!                                                                      *
!cjb  myid, nprocs = NOT args
      SUBROUTINE mcpmpi (nb, RESTRT, fhead)
!                                                                      *
!   This routine controls the computation  and storage of the values   *
!   and all indices of the angular coefficients                        *
!                                                                      *
!                                       k                              *
!                   T  (ab)            V  (abcd)                       *
!                    rs                 rs                             *
!                                                                      *
!   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
!   c and d are orbital sequence numbers.  r and s are configuration   *
!   state function indices.                                            *
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, RKCO_GG,       *
!                        TNSRJJ.                                       *
!               [GENMCP]: FNDBEG, SETSDA, SORT.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 28 Sep 1993   *
!   Modified by C. Froese Fischer for block computation.               *
!   Modified by G. Gaigalas and J. Bieron for new spin-angular         *
!   integration                                        01 April 2012   *
!                                                                      *
!***********************************************************************
!                                                                      *
!   Modified for CSFG list                                             *
!   Chongyang Chen,  Fudan University, Shanghai,           Oct  2023   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,   ONLY:  NNNW, KEYORB
      USE csfg_memory_man
!CYC      USE mpi_C,           ONLY:  myid, nprocs
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE outsdampi_I
      USE fndbeg_I
      USE alcbuf_I
      USE rkco_GG_I
      USE setsda_csfg_I
      USE sort_I
!-----------------------------------------------
!   C O M M O N  B L O C K S
!-----------------------------------------------
      USE BUFFER_C,  ONLY: NVCOEF, LABEL, COEFF
      USE DEBUG_C,   ONLY: LDBPA
      USE DEFAULT_C, ONLY: NDEF
      USE iccu_C,    ONLY: ICCUT
      USE MCP_C,     ONLY: KMAX, DIAG, LFORDR
      USE ORB_C,     ONLY: NCF, IQA, NW
      USE STAT_C,    ONLY: JQSA, JCUPA
!CYC
      Use csfg_decide_C, ONLY : LSYMD
      Use symexpand_mod, ONLY : TotCSFs_perblock, ncsfDF1
      Use symmatrix_mod, KMAXTmp=>KMAX
      Use symmatrix_restart_C, ONLY : ncsfDF1_core 
      Use csfg_tv_C

      IMPLICIT NONE
      EXTERNAL CORD
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NB
!     INTEGER  :: MYID
!     INTEGER  :: NPROCS
      LOGICAL , INTENT(IN) :: RESTRT
      CHARACTER(len=*), INTENT(IN) :: FHEAD
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: KEYSQ = KEYORB*KEYORB
!GG      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-10
!GG      REAL(DOUBLE), PARAMETER :: CUTOFF = 0.0D0
!cjb  REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-08
      REAL(DOUBLE), PARAMETER :: CUTOFF0 = 1.0D-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      !INTEGER, DIMENSION(:), pointer :: LLISTV
      !INTEGER, DIMENSION(:), allocatable :: MAPCS
      INTEGER :: JASTRT,JBSTRT,NPOS,K,JA,JB,IA,IB,LAB,I,IC, &
         ID, NSWAP, ISWAP, LAC, LBD, NTGI
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      REAL(DOUBLE) :: VCOEFF
      LOGICAL :: F0INT, LINCR

      INTEGER :: NCFSMALL, NCFGEN, NCFTOT, IOUTF, IOUTT
      INTEGER(BYTE) :: IQAICC(NW), IQAIRR(NW), IQADiff(NW)
      INTEGER :: NELMNTCC, IUNITF, NUMOS, NUMRKCO, NUMV(0:KMAX)
      INTEGER :: NDIM, IPOS, NUMP, NCOLG
      
!=======================================================================
!
      IOUTF = 90 + NB
      IOUTT = 9900 + NB
!
! Scratch file stores the IROWs for nonzero elements in sparse
! structure.
      OPEN (IOUTT, STATUS = 'SCRATCH', FORM = 'UNFORMATTED')
!
! Set NTYPE(6,JA/JB)--BEGIN
!

! Call findtype to initialize the values in symmatrix_C.mod
      NCFSMALL=TotCSFs_perblock(2,NB)
! Check: 
      IF (NCF.NE.NCFSMALL) THEN
        IF (myid .eq. 0) THEN
          WRITE(*,*)"Error, Unexpected NCF.NE.NCFSMALL within setham_gg"
        ENDIF
        STOP "Error, Unexpected NCF.NE.NCFSMALL within setham_gg ..."
      ENDIF

! Call findtype in order to define the type of CSF and largest value
! for the orbitals of the CSFG. 
!
      ALLOCATE(NTYPE(6,NCF))
      IF (MYID.EQ.0) THEN
         CALL FINDTYPE(NCFSMALL,NCFGEN,NCFTOT)
         if (NCFGEN.ne.TotCSFs_perblock(1,NB).or.                &
             NCFTOT.ne.TotCSFs_perblock(3,NB).or.                &
             NCFTr.NE.NCFTOT) then
           write(*,*)'NCFGEN=',NCFGEN,                           &
             ' TotCSFs_perblock(1,NB)=',TotCSFs_perblock(1,NB)
           write(*,*)'NCFTOT=',NCFTOT,                           &
             ' TotCSFs_perblock(3,NB)=',TotCSFs_perblock(3,NB)
           write(*,*)'NCFTr=',NCFTr, ' NCFTOT=',NCFTOT
           stop 'Unexpected error in mcpmpi_gg ===1=== ...'
         endif
      ENDIF
! Broadcast ...
      CALL MPI_Bcast(ncsfDF1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(MAXSPAN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NORBGEN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NCFTOT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(NTYPE(:,:),6*NCF,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
!
! Set NTYPE(6,JA/JB)-- END
!
!=======================================================================
!
! Allocate storage to the array that stores the lengths
! of the lists of V coefficients of each multipolarity
!
      CALL ALLOC(LLISTV, 0, kmax, 'LLISTV', 'MCPMPI')
!
! If this is a restart, determine the pair of CSFs with which
! the computation should begin.
! If not restart, initialize the counters for the total number of
! T coefficients and V coefficients of each multipolarity.
!
      IF (RESTRT) THEN
         PRINT *, 'Not ready for RESTART'
         STOP
         CALL FNDBEG (JASTRT, JBSTRT, npos, LLISTT, LLISTV)
      ELSE
         JASTRT = 1 + myid
         NPOS = 0
         LLISTT = 0
         DO K = 0, KMAX
            LLISTV(K) = 0
         ENDDO
      ENDIF
!
! Allocate storage for the arrays in BUFFER
!
      CALL ALCBUF (1)
!
! Write header to  .dbg  file if appropriate
!
      IF (LDBPA(2) .OR. LDBPA(3)) WRITE (99,300)
!
! JA and JB respectively refer to the initial and final states
! in the list of NCF configurations
!
! NCF:  NCSFsmall, Number of CSFGs in rcsfg.inp ....
! JA-th CSFG, Column 
!
! Initialize arrays to store all the spin-angular coefficiests
! involving the Column JA-th CSFG
!
! IROWSYM(1:NCFTr,NTYPE(1,JA)) to store the Sparse structure of
! the H-matrix, i.e., the data to write into mcpXXX.30, at the end of
! the DO 5 LOOP:  "JA, JB, npos" 
      ALLOCATE(IROWSYM(NCFTOT,MAXSPAN))

! For each sub-block generated by JB -- JA CSFG pairs, record the
! positions of nonzero elements.
      ALLOCATE(LNonZero(MAXSPAN,MAXSPAN))
!
! For each sub-block generated by JB -- JA CSFG pairs, record the
! nonzero numbers of the columns expanded by the JA-th CSFG.
      ALLOCATE(NELCSYM(MAXSPAN)) ! arrays of npos ...

! Labelling CSFs number, are used by rmcdhf_mpi_csfg, for
! estimations of eigen-pairs. 
      ncsfDF1_core = 0
!
! The number of Columns calculated by this MPI process
      NODENCOLS = 0
! The detailed Column index of the NODENCOLS, they will be used in
! rmcdhf_mpi_csfg.
      ALLOCATE(NODECOLS(NCFTOT))

! Positions of the FICTIOUS CSFs
      IRFICT = NCF + 1
      ICFICT = NCF + 6

! Begin the loops:
! Number of calling ONESCALAR, in which the returned coefficients
! should be recorded.
      NONESCALAR = 0

! Number of calling RKCO_GG, in which the returned coefficients
! should be recorded. 
      NRKCO = 0
!
! Column counter in the fully expanded CSF list
      ICTOT = 0
!
! The total number of non-zero matrixelement calculated by MYID node.
      NELMNTCC = 0

! Number of non-zero JA - JB pairs.
      IPOS = 0
      NDIM = NCFTr
! Deatailed ICG - IRG indexes for the interacting CSFG pairs
      CALL ALLOC(ICGPDET, NDIM, 'ICGPDET', 'MCPMPI_CSFG')
      CALL ALLOC(IRGPDET, NDIM, 'IRGPDET', 'MCPMPI_CSFG')

! JA-th CSFG, Column index
      NCOLG = 0      
      DO 5 JA = JASTRT, NCFSMALL, nprocs
         NCOLG = NCOLG + 1
         IF (JA.GT.1) ICTOT = MAP(JA-1)
         IQAICC(1:NW) = IQA(1:NW, JA)

         ! Full matrixelements of the diagonal BLOCK generated within
         ! the JA-th CSFG with JA <= ICCUT should be all kept in
         ! the zero-first calculation.
         IF (JA.LE.ICCUT(nb)) THEN
           LICCUT = .FALSE.
         ELSE
           LICCUT = .TRUE.
         ENDIF

         ! Number of non-zeros of every column generated by JA CSFG
         NELCSYM = 0
 
         ! Row counter in the fully expanded CSF and H-matrix 
         IRTOT = 0

         ! JB-th CSFG, Row index
         JBSTRT = 1

         DO 4 JB = JBSTRT, JA
           IF (JB.GT.1) IRTOT = MAP(JB-1)

           ! Zero-first approximation
           IF (DIAG .OR. (LFORDR .AND. (JB .GT. ICCUT(nb))) ) THEN
             IF (JB /= JA) CYCLE
           ENDIF

           !H contains only ONE- and TWO-body interactions
           IQAIRR(1:NW) = IQA(1:NW, JB)
           IQADiff = IQAICC - IQAIRR
           IF (COUNT(IQADiff /= 0) > 4) CYCLE

           ! Initialization
           LTRANSPOSE = .FALSE.
           LTRANSFER = .false.
           LNonZero = .false.
           IF (NTYPE(1,JA) == 1 .AND. NTYPE(1,JB) == 1) THEN
             IF (JA == JB) THEN
               IUNITF = 2001
               CALL SPINANGULAR1(JA, JB, IUNITF, NPOS)
             ELSE
               IUNITF = 2011
               CALL SPINANGULAR11(JA, JB, IUNITF, NPOS)
             ENDIF

           ELSEIF (NTYPE(1,JA) == 2 .AND. NTYPE(1,JB) == 1) THEN
             IUNITF = 2012
             CALL SPINANGULAR12(JA, JB, IUNITF, NPOS)
           ELSEIF (NTYPE(1,JA) == 1 .AND. NTYPE(1,JB) == 2) THEN
             IUNITF = 2012
             LTRANSPOSE = .TRUE.
             CALL SPINANGULAR12(JA, JB, IUNITF, NPOS)

           ELSEIF (NTYPE(1,JA) == 3 .AND. NTYPE(1,JB) == 1) THEN
             IUNITF = 2013
             CALL SPINANGULAR13(JA, JB, IUNITF, NPOS)
           ELSEIF (NTYPE(1,JA) == 1 .AND. NTYPE(1,JB) == 3) THEN
             IUNITF = 2013
             LTRANSPOSE = .TRUE.
             CALL SPINANGULAR13(JA, JB, IUNITF, NPOS)

           ELSEIF (NTYPE(1,JA) == 4 .AND. NTYPE(1,JB) == 1) THEN
             IUNITF = 2014
             CALL SPINANGULAR14(JA, JB, IUNITF, NPOS)
           ELSEIF (NTYPE(1,JA) == 1 .AND. NTYPE(1,JB) == 4) THEN
             IUNITF = 2014
             LTRANSPOSE = .TRUE.
             CALL SPINANGULAR14(JA, JB, IUNITF, NPOS)

           ELSEIF (NTYPE(1,JA) == 5 .AND. NTYPE(1,JB) == 1) THEN
             IUNITF = 2015
             CALL SPINANGULAR15(JA, JB, IUNITF, NPOS)
           ELSEIF (NTYPE(1,JA) == 1 .AND. NTYPE(1,JB) == 5) THEN 
             IUNITF = 2015
             LTRANSPOSE = .TRUE.
             CALL SPINANGULAR15(JA, JB, IUNITF, NPOS)

           ELSEIF (NTYPE(1,JA) == 2 .AND. NTYPE(1,JB) == 2) THEN
             IF (JA == JB) THEN
               IUNITF = 2002
               CALL SPINANGULAR2(JA, JB, IUNITF, NPOS)
             ELSE
               IUNITF = 2022
               CALL SPINANGULAR22(JA, JB, IUNITF, NPOS)
             ENDIF 

           ELSEIF (NTYPE(1,JA) == 3 .AND. NTYPE(1,JB) == 2) THEN
             IUNITF = 2023
             CALL SPINANGULAR23(JA, JB, IUNITF, NPOS)
           ELSEIF (NTYPE(1,JA) == 2 .AND. NTYPE(1,JB) == 3) THEN
             IUNITF = 2023
             LTRANSPOSE = .TRUE.
             CALL SPINANGULAR23(JA, JB, IUNITF, NPOS)

           ELSEIF (NTYPE(1,JA) == 4 .AND. NTYPE(1,JB) == 2) THEN
             IUNITF = 2024
             CALL SPINANGULAR24(JA, JB, IUNITF, NPOS)
           ELSEIF (NTYPE(1,JA) == 2 .AND. NTYPE(1,JB) == 4) THEN
             IUNITF = 2024
             LTRANSPOSE = .TRUE.
             CALL SPINANGULAR24(JA, JB, IUNITF, NPOS)

           ELSEIF (NTYPE(1,JA) == 5 .AND. NTYPE(1,JB) == 2) THEN
             IUNITF = 2025
             CALL SPINANGULAR25(JA, JB, IUNITF, NPOS)
           ELSEIF (NTYPE(1,JA) == 2 .AND. NTYPE(1,JB) == 5) THEN
             IUNITF = 2025
             LTRANSPOSE = .TRUE.
             CALL SPINANGULAR25(JA, JB, IUNITF, NPOS)

           ELSEIF (NTYPE(1,JA) == 3 .AND. NTYPE(1,JB) == 3) THEN
             IF (JA == JB) THEN
               IUNITF = 2003
               CALL SPINANGULAR3(JA, JB, IUNITF, NPOS)
             ELSE
               IUNITF = 2033
               CALL SPINANGULAR33(JA, JB, IUNITF, NPOS)
             ENDIF

           ELSEIF (NTYPE(1,JA) == 4 .AND. NTYPE(1,JB) == 3) THEN
             IUNITF = 2034
             CALL SPINANGULAR34(JA, JB, IUNITF, NPOS)
           ELSEIF (NTYPE(1,JA) == 3 .AND. NTYPE(1,JB) == 4) THEN 
             IUNITF = 2034
             LTRANSPOSE = .TRUE.
             CALL SPINANGULAR34(JA, JB, IUNITF, NPOS)

           ELSEIF (NTYPE(1,JA) == 5 .AND. NTYPE(1,JB) == 3) THEN
             IUNITF = 2035
             CALL SPINANGULAR35(JA, JB, IUNITF, NPOS)
           ELSEIF (NTYPE(1,JA) == 3 .AND. NTYPE(1,JB) == 5) THEN 
             IUNITF = 2035
             LTRANSPOSE = .TRUE.
             CALL SPINANGULAR35(JA, JB, IUNITF, NPOS)

           ELSEIF (NTYPE(1,JA) == 4 .AND. NTYPE(1,JB) == 4) THEN
             IF (JA == JB) THEN
               IUNITF = 2004
               CALL SPINANGULAR4(JA, JB, IUNITF, NPOS)
             ELSE
               IUNITF = 2044
               CALL SPINANGULAR44(JA, JB, IUNITF, NPOS)
             ENDIF

           ELSEIF (NTYPE(1,JA) == 5 .AND. NTYPE(1,JB) == 4) THEN
             IUNITF = 2045
             CALL SPINANGULAR45(JA, JB, IUNITF, NPOS)
           ELSEIF (NTYPE(1,JA) == 4 .AND. NTYPE(1,JB) == 5) THEN 
             IUNITF = 2045
             LTRANSPOSE = .TRUE.
             CALL SPINANGULAR45(JA, JB, IUNITF, NPOS)

           ELSEIF (NTYPE(1,JA) == 5 .AND. NTYPE(1,JB) == 5) THEN
             IF (JA == JB) THEN
               IUNITF = 2005
               CALL SPINANGULAR5(JA, JB, IUNITF, NPOS)
             ELSE
               IUNITF = 2055
               CALL SPINANGULAR55(JA, JB, IUNITF, NPOS)
             ENDIF
           ELSE
             STOP "Sth Error within NTYPE(1,:) in mcpmpi_csfg ..." 
           ENDIF
!
! Record the positions of non-zero CSFG-pairs, JA (IC) -- JB (IR). 
           IF (LTRANSFER) THEN
             IPOS = IPOS + 1
             IF (IPOS > NDIM) THEN
               NDIM = NDIM * 2
               CALL RALLOC(IRGPDET, NDIM, 'IRGPDET', 'MCPMPI_CSFG')
             ENDIF
             IRGPDET(IPOS) = JB
           ENDIF
    4    CONTINUE
!
! Record the data for JA-th CSFG
!
! Section 1:
! Record the end position of the interacting CSFG-pairs list for JA-th
! CSFG.
!
         ICGPDET(NCOLG) = JA
         ICGPEND(NCOLG) = IPOS 

! Section 2:
! The columns included by the MYID MPI process.
         DO I = 1, NTYPE(2, JA)
           NODENCOLS = NODENCOLS + 1
           NODECOLS(NODENCOLS) = ICTOT + I
! The labelling CSFs (before any correlation CSFG) included by MYID
           IF (NODECOLS(NODENCOLS) .LE. ncsfDF1) &
             ncsfDF1_core = ncsfDF1_core + 1
!
! Diagonal IROWSYM(NELCSYM(I), I) should be equal to Index of Column
!
           IF (NODECOLS(NODENCOLS) .NE. IROWSYM(NELCSYM(I), I)) THEN
             WRITE(*,*)NODECOLS(NODENCOLS),IROWSYM(NELCSYM(I), I)
             STOP "NODECOLS(NODENCOLS) .NE. IROWSYM(NELCSYM(I), I) ..."
           ENDIF
!
! Output for check, by comparing to those calculated by rci_mpi_csfg
!
           !WRITE(9900, *)ICTOT+I, NELCSYM(I), &
           !              '   JA = ', JA, '  JA*3+5=', JA*3+5
           !WRITE(IOUTF, '(2i8,5000I8)')ICTOT+I, &
           !      NELCSYM(I), (IROWSYM(JB,I), JB = 1,NELCSYM(I))
!
! Section 3:
! The total number of non-zero matrixelement calculated by MYID node.
           NELMNTCC = NELMNTCC + NELCSYM(I)
!
! Section 4:
! Output results into the SCRATCH file, read in setsda_csfg.f90
           WRITE(IOUTT)ICTOT+I, &
                 NELCSYM(I), (IROWSYM(JB,I), JB = 1,NELCSYM(I))
         ENDDO

         if (mod(JA-1,100) == 0) PRINT '(A4,I7,A9,I10,A9,I3,A8,I3)',  &
            'Row ',JA,' nnonz = ',npos,' block = ',nb,' myid = ',myid

    5 CONTINUE
!
! Total number of interacting CSFG-pairs (non-zero block matrixelement).
      NPAIRS(1,NB) = NCOLG
      NPAIRS(2,NB) = IPOS
!
! Spin-angular coefficients for block NB are calculated, write the
! sparse structure into mcpXXX.30
!
!   Write out a report for this run
!
      CALL MPI_Allreduce (NPAIRS(2,NB), NUMP, 1, MPI_INTEGER,  &
         MPI_SUM, MPI_COMM_WORLD, ierr)
      IF (MYID == 0) WRITE (*,*) NUMP, ' Interacting CSFG-Pairs'

      CALL MPI_Allreduce (NONESCALAR, NUMOS, 1, MPI_INTEGER,  &
         MPI_SUM, MPI_COMM_WORLD, ierr)
      IF (MYID == 0) WRITE (*,*) NUMOS, ' Effective calls ONESCALAR'
      !WRITE (*,*) LLISTT, ' T coefficients generated;'

      CALL MPI_Allreduce (NRKCO, NUMRKCO, 1, MPI_INTEGER,  &
         MPI_SUM, MPI_COMM_WORLD, ierr)
      IF (MYID == 0) WRITE (*,*) NUMRKCO, ' Effective calls RKCO_GG'

      CALL MPI_Allreduce (LLISTV(0:KMAX), NUMV(0:KMAX), KMAX+1, &
                          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      IF (MYID == 0) THEN
        DO K = 0, KMAX
          WRITE (*,*) NUMV(K),' V(k=',K,') coefficients generated;'
        ENDDO
      ENDIF
      CALL DALLOC(LLISTV, 'LLISTV', 'MCPMPI')
!
!   Set up sparse structure definition arrays in file 30
!
      REWIND(IOUTT)
      CALL SETSDA_CSFG (outsdampi, IOUTT, NELMNTCC, LDBPA(4),          &
                                          NB, MYID, NPROCS, FHEAD)
! Close the scratch file
      CLOSE(UNIT=IOUTF)

!
! Deallocate the arrays
!
      DEALLOCATE(NTYPE)
      DEALLOCATE(IROWSYM) 
      DEALLOCATE(LNonZero)
      DEALLOCATE(NELCSYM)
      DEALLOCATE(NODECOLS)
!
! Deallocate storage that is no longer required
!
      CALL DALLOC (iqa, 'IQA', 'MCPMPI')
      CALL DALLOC (jqsa, 'JQSA', 'MCPMPI')
      CALL DALLOC (jcupa, 'JCUPA', 'MCPMPI')
      CALL ALCBUF (3)

      RETURN
  300 FORMAT (/'From MCPMPI:')
      END SUBROUTINE mcpmpi
