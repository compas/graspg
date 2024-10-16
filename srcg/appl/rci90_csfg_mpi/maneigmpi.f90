!***********************************************************************
!                                                                      *
      SUBROUTINE MANEIG(IATJPO, IASPAR, NELMNT_a)
!                                                                      *
!   This module  manages the  operation of the  eigensolvers and the   *
!   storage of the eigenpairs.  There are two principal branches:      *
!                                                                      *
!      (1) Matrix of order 1: the trivial case                         *
!      (2) Matrix of order greater than 1: there are two branches      *
!             (i) Matrices of order less than or equal to IOLPCK:      *
!                 eigenpairs are found using LAPACK SUBROUTINEs        *
!            (ii) Matrices of order greater than IOLPCK: eigenpairs    *
!                 are found using DVDSON; this involves up to three    *
!                 steps:                                               *
!                    (a) The matrix is analysed to determine its       *
!                        block structure (only irreducibe matrices     *
!                        are correctly treated by DVDSON)              *
!                    (b) Eigenpairs are extracted for each block       *
!                    (c) The appropriate eigenpairs are selected and   *
!                        stored                                        *
!                 Different methods of storage and different           *
!                 versions of the matrix-vector multiply are used      *
!                 depending upon the order and density of the matrix   *
!
!cjb comments [re-introduced] from GRASP-2013
!   iatjpo - Output. (2j+1) of the dominant coefficients
!   iaspar - Output, Parity (1 or -1) of the dominant coefficients
!   nelmnt - Number of non-zero matrix elements on the current process *
!            stored in the common block                                *
!   nelmnt_a - The maximum number of non-zero matrix elements for all  *
!             nodes                                                    *
!                                                                      *
!     value  meaning  of  method used IV                               *
!       1    LAPACK                  =>  NCF < IOLPCK                  *
!       2    Davidson, Dense-Memory  =>  dense requires less memory    *
!       3    Davidson, Sparse-Memory =>  sparse requires less memory   *
!       4    Davidson, Sparse-Disk   =>  memory requirement too large  *
!                                        for at least one  process     *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, ISPAR, ITJPO, posfile,         *
!                        RALLOC.                                       *
!               [RCI92]: DNICMV, SPICMVsym, SPODMV.                    *
!               [DVDSON]: DVDSON.                                      *
!               [AUXBLAS]: DINIT/SINIT.                                *
!               [BLAS]: DCOPY/SCOPY, DSWAP/SSWAP.                      *
!               [LAPACK]: DSPEVX/SSPEVX.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 27 Sep 1993   *
!   Modified Misha Saparov                                  Feb 1997   *
!            Charlotte F. Fischer                           May 1997   *
!                 Except for the disk version, all matrices have       *
!                 diagonals  shifted by EAV                            *
!  All arrays allocated here are de-allocated except pneval and pnevec *
!  which will be de-allocated in matrix.                               *
!  Block Version  By Xinghong He          Last revision: 18 Jun 1998   *
!   Modify  by G Gaigalas                         May 2021             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:38   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!                                                                      *
!   Last modification for CSFG by Chongyang Chen  Aug 2024             *
!                                                                      *
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE, LONG
      USE memory_man
      USE eigv_C
      USE fposition_C
      USE hmat_C
      USE prnt_C
      USE WCHBLK_C
      USE iounit_C
      USE mpi_C
      USE fposition_C
      use symexpand_mod
      use symmatrix_mod
      use csfg_decide_C, ONLY: MaxMemPerProcs 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dnicmv_I
      USE spicmv2_I
      USE spodmv_I
      USE posfile_I
      USE dinit_I
      USE dspevx_I
      USE iniestsd_I
      USE gdvd_I
      use iniestmpi_I
      USE iniestdm_I
      USE itjpo_I
      USE ispar_I
      USE spicmvmpi_I
      IMPLICIT NONE
!-----------------------------------------------
!  E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      REAL(DOUBLE) :: DLAMCH
      EXTERNAL     :: DLAMCH
!      EXTERNAL     :: DNICMV,SPICMVsym,SPODMV
      EXTERNAL     :: SPICMVsym
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(OUT) :: IATJPO
      INTEGER, INTENT(OUT) :: IASPAR
      INTEGER(LONG)        :: NELMNT_a
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER      :: IOLPCK0 = 1000
! GG      REAL(DOUBLE), PARAMETER :: ABSTOL = 1.0D-10
!cjb  NINCOR
!cjb  INTEGER, PARAMETER      :: NINCOR = 1         ! To enforce DISK
!CYC      INTEGER(LONG), PARAMETER      :: NINCOR = 268435456 ! = 2 GB or  memory
      INTEGER(LONG)  :: NINCOR != 268435456 ! = 2 GB or  memory
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!cff   ... variable for deciding between sparse or dense
!     Memory need:
!     sparse : nstore_s = nelmnt + (isize*(nelmnt + NCF+1))/8
!     dense  : nstore_d = (NCF*NCF+1))/(2*nprocs)    (on average)
!     Use Davidson (disk)  if nstore_s > NINCOR
!     Use Davidson (memory, sparse) if nstore_s < nstore_d
!     Use Davidson (memory, dense) if nstore_d < nstore_s
!-----------------------------------------------
      INTEGER :: LWORK, LIWORK
      INTEGER(LONG) :: NSTORE_s, NSTORE_d
      INTEGER :: IOFSET, NCFT
      INTEGER :: NROWS, I, NDENSE, NCFDUM, ICCUTDUM, MYIDDUM, &
         NPROCSDUM, NELC, IR, NVECMN, NVEX, M, INFO, LOC, NBRKEV,   &
         IMV, NDENSE_L, LIM, MAXITR, MBLOCK, NEND, &
         ILOW, IHIGH, NIV, NLOOPS, NMV, J, IA,   &
         idummy, IOLPCK, ICOL
      INTEGER :: iitmp, irestart_gdvd, istep, ncount0, ncount1, NCF
      INTEGER :: ncount2, ncount_max, ncount_rate, nminutes, nmv0
      INTEGER :: nseconds, nstep, AllocateStatus, DeAllocateStatus
      REAL(DOUBLE) :: DSIZE, HSIZE, SYSMemT, WSIZE, ELSTO,  DUMMY, &
          DIATMP, CRITE, CRITC, CRITR, ORTHO, DMUNGO, AMAX, WA, ABSTOL
! GG          DIATMP, CRITE, CRITC, CRITR, ORTHO, DMUNGO, AMAX, WA
      LOGICAL :: HIEND, LDISC, SPARSE
      CHARACTER(LEN=8) :: CNUM
      REAL(DOUBLE), DIMENSION(:), pointer :: w, z, work, diag
      INTEGER, DIMENSION(:), pointer :: iwork, ifail, jwork

      REAL(DOUBLE)  :: endtime, starttime
!CYC-2023/12/04
      INTEGER(LONG)        :: NELMNTMax
      INTEGER(LONG), DIMENSION(:), POINTER   :: NELMNTArr(:)
!-----------------------------------------------------------------------
      NINCOR = 2147483648d0
      IOLPCK = min(ncsfDF1,IOLPCK0)
      If (IOLPCK .LT. NVECMX) then
        write(*,*)"ncsfDF1=", ncsfDF1
        write(*,*)"NVECMX=", NVECMX 
        stop "Error, IOLPCK .LT. nvecmx, search too many levels ..."
      endif
      NCF = NCFTr
      ABSTOL = 2*DLAMCH('S')
!      IF (MYID == 0) WRITE (6, *) 'Calling maneig...'
      NCFt = NCF
      NMV = 0
!  (nrows+1) is the number of records of the present block's .res file

      NROWS = (NCF - MYID - 1 + NPROCS)/NPROCS
      IF (NCF < NPROCS) NROWS = NCF/(MYID + 1)
      !CALL posfile (1, imcdf, nrows+1)
      CALL POSFILE (0, IMCDF, NPOSITION)

      IF (NCF == 1) THEN
!-----------------------------------------------------------------------
!
! (1) - Trivial    NCF = 1
!
!-------------------------------------------------------
         IF (myid .EQ. 0) WRITE (24, *) 'Trivial eigenvalue problem.'

!   Matrix of order 1: the trivial case; we assume that the value
!   of EAV is available

         CALL ALLOC (EVAL, 1,'EVAL', 'MANEIG' )
         CALL ALLOC (EVEC, 1, 'EVECO', 'MANEIG')
         EVAL(1) = 0.D0
         EVEC(1) = 1.D0

         ! Still read through the .res file
!GG
!GG Gediminas NIST 2005.11.03
!GG         READ (imcdf)
         DO I = 1, NROWS + 1
            READ (IMCDF)
         END DO

      ELSE                                       !if-2
!-----------------------------------------------------------------------
!
! (2) - Non trivial
!
!-------------------------------------------------------
!         if(myid==0) write(*,*) "NCSF=",NCF, " iolpck=", iolpck
         IF (NCF <= IOLPCK) THEN
!-----------------------------------------------------------------------
!
! (2.1) - LAPACK    Dense, Memory,
!
!-------------------------------------------------------
            IF (MYID == 0) THEN
               WRITE (6, *) &
                  'LAPACK routine DSPEVX selected for eigenvalue problem.'
               WRITE (24, *) &
                  'LAPACK routine DSPEVX selected for eigenvalue problem.'
            ENDIF

!   Allocate storage for the dense representation of the matrix
!   and initialize emt

            NDENSE = (NCF*(NCF + 1))/2
            CALL ALLOC (EMT, NDENSE, 'EMT', 'MANEIG')
            CALL DINIT (NDENSE, 0.0D00, EMT, 1)

!   Read the matrix into position from the disc file; it's already
!   been properly positioned.

            CALL ALLOC (WORK, NCF,'WORK', 'MANEIG' )
            CALL ALLOC (IROW, NCF, 'IROW', 'MANEIG')
            READ (IMCDF) NCFDUM, ICCUTDUM, MYIDDUM, NPROCSDUM
            IF (NCF/=NCFDUM .OR. MYID/=MYIDDUM .OR. NPROCSDUM/=NPROCS) STOP &
               'maneig:1'
!CYC: Modifications for CSFG
!            DO I = MYID + 1, NCF, NPROCS
            DO I = 1, NODENCOLS
               ICOL = NODECOLS(I)
!               IOFSET = (I*(I - 1))/2
               IOFSET = (ICOL*(ICOL-1))/2
               READ (IMCDF) NELC, ELSTO, (WORK(IR),IR=1,NELC), (IROW(IR),IR=1,&
                  NELC)
               ! In the row-mode of the lower triangle,
               ! diagonal is the last one
               ! Check
                     IF (IROW(NELC).NE.NODECOLS(I)) THEN
                      WRITE(*,*)'Wrong column index....',              &
                'I,NODECOLS(I), IROW(NELC)=', I,NODECOLS(I),IROW(NELC)
      STOP 'Unexpected IROW(NELC).NE.NODECOLS(I) in maneig.f Line 531 '
                     ENDIF

               DO IR = 1, NELC - 1
                  EMT(IOFSET+IROW(IR)) = WORK(IR)
               END DO
               EMT(IOFSET+IROW(NELC)) = WORK(NELC) - EAV

            END DO

            CALL DALLOC (WORK, 'WORK', 'MANEIG')
            CALL DALLOC (IROW, 'IROW', 'MANEIG')

! Let each node have a complete copy of EMT

            CALL gdsummpi (EMT, NDENSE)

!   Find the eigenpairs
!
!    ivec() - serial numbers of eigenstates of the current block
!    iccmin() - serial numbers of eigenstates of all blocks.
!    nvecmn - minimum serial number of the eigenstates of the block
!    NVECMX - maximum .............
!    nvex - clear from def: NVECMX-NVECMN+1

            NVECMN = NCF
            DO I = 1, NVEC
               NVECMN = MIN(NVECMN,IVEC(I))
            END DO
            NVEX = NVECMX - NVECMN + 1
            CALL ALLOC (W, NVEX, 'W', 'MANEIG')
            CALL ALLOC (Z, NCF*NVEX,'Z', 'MANEIG' )
            CALL ALLOC (WORK, NCF*8,'WORK', 'MANEIG' )
            CALL ALLOC (IWORK, NCF*5,'IWORK', 'MANEIG' )
            CALL ALLOC (IFAIL, NVEX, 'IFAIL', 'MANEIG')
            CALL DSPEVX ('V', 'I', 'U', NCF, EMT, DUMMY, DUMMY, NVECMN, NVECMX&
               , ABSTOL, M, W, Z, NCF, WORK, IWORK, IFAIL, INFO)
            IF (INFO /= 0)                                             &
                CALL STOPMPI('maneig: Failure in DSPEVX [LAPACK]',myid)
            CALL DALLOC (WORK, 'WORK', 'MANEIG')
            CALL DALLOC (IWORK, 'IWORK', 'MANEIG')
            CALL DALLOC (IFAIL, 'IFAIL', 'MANEIG')
            CALL DALLOC (EMT, 'EMT', 'MANEIG')

!   Store the eigenpairs in their proper positions EVAL() and EVEC()

            CALL ALLOC (EVAL, NVEC,'EVAL', 'MANEIG')
            CALL ALLOC (EVEC, NCF*NVEC, 'EVEC', 'MANEIG')

            DO I = 1, NVEC
               LOC = IVEC(I)
               EVAL(I) = W(LOC - NVECMN + 1)
               IOFSET = NCF*(I - 1)
               LOC = NCF*(LOC - NVECMN)
               CALL DCOPY (NCF, Z(LOC + 1), 1, EVEC(IOFSET+1), 1)
            END DO
            CALL DALLOC (W, 'W', 'MANEIG')
            CALL DALLOC (Z, 'Z', 'MANEIG')

         ELSE
!-----------------------------------------------------------------------

! (2.2) - DVDSON --- preparation work
!
!-------------------------------------------------------
            IF (myid == 0) WRITE (24,*)'DVDSON routine selected' &
                                      //' for eigenvalue problem;'
!--------------------------------------------------------------
!            IF (myid == 0) print *,'zou: nelmnt',NELMNT
!            print *, 'myid, nprocs,nincor', myid, nprocs, nincor

!           If on any process more memory is needed than NINCOR, use disc.
!           If not, should the matrix be sparse or dense (distributed)?
!           sparse requires = nelmnt + (isize*(nelmnt+NCF+1))/8
!           dense  storage = (NCF*(NCF+1)/2/nproc - on average
!           (because of integer arithmetic nproc-1 is added to numerator

!           This is a comparion of matrix elements only.
!
            NSTORE_s = NELMNT_a + ((NELMNT_a + NCF+1)*4 )/8
            NSTORE_d = ((NCF+nprocs-1)/nprocs)
            NSTORE_d = NSTORE_d*((NCF+1)/2)
!            IF (myid .EQ. 0) write(*,*) "nstore_s=", nstore_s,   &
!                             "nstore_d=", nstore_d, 'nincor=', nincor
!===================================================================
! CYC:
!  Determine imethod: Sparse memory / Sparse disk
!            ! Total physical memory
!            SYSMemT = SYSTEMMem()  ! GB

            ! Non-zero H-matrix size, and the row-index array size
            HSIZE = (NELMNT*(8 + 4)) / 1.0D9  ! GB

            ! Working array size (used in dvdson module)
            !LIM = 3*NVECMX + 30
            LIM = MIN(NCF,2*NVECMX + 60)
            LWORK = 2*NCFt*LIM+LIM*LIM*2+11*LIM+NVECMX
            IF (LWORK < 0 ) THEN
              WRITE(*,*)"LWORK Overflow, the calculation is too" &
                         // " large-scale: NCF is too large,"   &
                         // " or too many levels to be resolved."
              STOP "Error, Overflow occurs 1, &
                    out of the ability of the present version!"
            ENDIF

            ! Additional size for other work-arrays, CHECK needed
            LWORK = LWORK + min(10, NVECMX) * NCFt
            WSIZE = 1.0d-09 * LWORK * 8         ! GB

            ! Dense matrix size for each process
            DSIZE = 8 * NCFt * (NCFt + 1) / 2 / 1.0D9 / NPROCS ! GB

            ! Determine SPARSE / DENSE, Memory / DISC
!CYC-2023/12/05: 
! Set always SPARSE = .TRUE., determine LDISC using MaxMemPerProcs
            IF (MYID .EQ. 0) THEN
              SPARSE = .TRUE.   
              IF ((HSIZE + WSIZE) < MaxMemPerProcs) THEN
                LDISC = .FALSE.
              ELSE
                LDISC = .TRUE.
              ENDIF
              !IF ((WSIZE + DSIZE)*NPROCS .GT. 0.8*SYSMemT) THEN
              !  SPARSE = .TRUE.
              !  IF ((SYSMemT-(HSIZE + WSIZE)*NPROCS).GT.50.0) THEN
              !    ! System memory is enough for SPARSE / MEMORY
              !    ! At least 50 GB memory left for other arrays /
              !    ! system.
              !    LDISC = .FALSE.
              !  ELSE
              !    LDISC = .TRUE.
              !  ENDIF
              !ELSE
              !   ! System memory is enough for DENSE matrix loaded into
              !   ! MEMORY.
              !   SPARSE = .FALSE.
              !   LDISC  = .FALSE.
              !ENDIF
              WRITE(*,*)"Estimation of Memory Needed:"
              WRITE(*,'(1x,a15, i10      )')"NCSFs         =", NCFt
              WRITE(*,'(1x,a15, f10.2, a4, 1x, a)')"HSIZE (Total) =", &
                       HSIZE*NPROCS, ' GB,', &
              'which is nearly constant nomatter how many cores used.'
              WRITE(*,'(1x,a15, f10.2, a4, 1x, a)')"WSIZE (Total) =", &
                       WSIZE*NPROCS, ' GB,', &
              'which increases linearly with the cores used.'
              WRITE(*,'(1x,a15, f10.2, a4)')"HWSIZE(Total) =", &
                       (HSIZE + WSIZE)*NPROCS, ' GB '
              !WRITE(*,'(a8, f10.2, a4)')"DWSIZE=", DSIZE + WSIZE, ' GB'
              !WRITE(*,'(a8, f10.2, a4)')"SSIZE =", SYSMemT/nprocs, ' GB'
              !WRITE(*,'(a8, f10.2, a4)')"MemLF1=",                    &
              !         SYSMemT -  (HSIZE + WSIZE) * NPROCS, ' GB'
              ! Set SPARSE always .TRUE. In case of DENSE/Memory model
              ! can be employed, the calculation would be small. Though,
              ! in such cases, DENSE/Memory model could save time, it
              ! doesn't matter much.
              ! Temporary use
              !WRITE(*,*)"SPARSE=", SPARSE
              !WRITE(*,*)"LDISC =", LDISC
            ENDIF

            CALL MPI_Bcast (SPARSE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            CALL MPI_Bcast (LDISC,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            call MPI_Barrier(MPI_COMM_WORLD, ierr)

            ALLOCATE(NELMNTArr(1:NPROCS), STAT = AllocateStatus)
            CALL MPI_GATHER(NELMNT, 1, MPI_LONG, NELMNTArr, 1,  &
                            MPI_LONG, 0, MPI_COMM_WORLD, ierr)
            NELMNTTMP = SUM(NELMNTArr(1:NPROCS))
            IF (MYID == 0) WRITE(*,*) &
               "Total Number of Non-Zero matrix elements:", NELMNTTMP
            IF (MYID == 0) WRITE(*,*) &
               "Sparsity of H-Matrix:", (2.0d0*NELMNTTMP-NCFt)/NCFt/NCFt
!===================================================================
!CYC-2023/12/04
! Check LDISC again, to avoid overflow 
            IF (.NOT.LDISC) THEN
              IF (MYID == 0) THEN
                NELMNTMax = MAXVAL(NELMNTArr)
                IF (NELMNTMax >= NINCOR) THEN
                  WRITE(*, *)"NELMNTMax, NINCOR =", NELMNTMax, NINCOR
                  WRITE(*, *)"Load Sparse H-matrix from Disc ..."
                  LDISC = .TRUE. 
                ELSE
                  WRITE(*, *)"Load Sparse H-matrix from Memory ..."
                ENDIF
              ENDIF
              CALL MPI_Bcast (LDISC,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            ELSE
              WRITE(*, *)"Load Sparse H-matrix from Disc ..."
            ENDIF 
            IF (MYID == 0) WRITE(*,*)
!CYC-2023/12/04 End
!===================================================================
!            LDISC = .FALSE.
!            IF ( NSTORE_s > NINCOR) THEN
!              .. we need to use disk
!               LDISC =.TRUE.
!               SPARSE=.TRUE.
!            ELSE

!               IF (NSTORE_s .LT. NSTORE_d ) THEN
!                  SPARSE = .TRUE.
!               ELSE
!                  SPARSE = .FALSE.
!               ENDIF
!            ENDIF

            CALL ALLOC (DIAG,NCF,'DIAG', 'MANEIG')
            CALL dinit (NCF, 0.d0, diag, 1)
            IF (LDISC) THEN
!-----------------------------------------------------------------------
!
! (2.2.1) - DVDSON --- Disk, load diagonal
!
!-------------------------------------------------------
               !IF (myid == 0)  THEN
               !   WRITE (*, *) ' matrix stored on disc;'
               !   WRITE (24, *) ' matrix stored on disc;'
               !END IF

!   Disk storage; necessarily sparse; one column of the matrix in
!   memory

               IMV = 1

!   Load diagonal - Each node will have the same, complete copy
!   after this if block
               print *, myid, 'Reading diagonal, unit imcdf =', imcdf
               READ (IMCDF) NCFDUM, ICCUTDUM, MYIDDUM, NPROCSDUM
               IF (NCF/=NCFDUM .OR. MYID/=MYIDDUM .OR. NPROCSDUM/=NPROCS) STOP &
                  'maneig:2'
!===================================================================
! CYC: CSFG-version implementation
               ! DO I = myid + 1, NCF, nprocs
               DO I = 1, NODENCOLS
                  READ (IMCDF) NELC,ELSTO,(dummy,IR=2,NELC), diatmp    &
                                       , (idummy, ir=1,nelc)
                  IF (IDUMMY.NE.NODECOLS(I)) THEN
                    WRITE(*,*)'Wrong column index....',                &
      'I,NODECOLS(I), IROW(NELC+IOFSET)=', I,NODECOLS(I),IDUMMY
      STOP 'Unexpected IDUMMY.NE.NODECOLS(I) in maneig.f ...'
                  ENDIF
                  !DIAG(I) = diatmp - EAV
                  DIAG(NODECOLS(I)) = diatmp - EAV
               ENDDO
               !print *, ' finished'
! CYC: CSFG-version implementation
!===================================================================
            ELSE
!-----------------------------------------------------------------------
!
! (2.2.2) - DVDSON --- Memory, load all
!
!-------------------------------------------------------
!
!   Core storage; load matrix into memory
!
!               LDISC = .FALSE.
               IF (SPARSE) THEN
!-----------------------------------------------------------------------

! (2.2.2.1) - DVDSON --- Memory, load all, sparse
!
!-------------------------------------------------------
                  !IF (MYID == 0) THEN
                  !   WRITE (24, *) ' sparse matrix stored in memory'
                  !   WRITE ( *, *) ' sparse matrix stored in memory'
                  !END IF

                  IMV = 2
                  CALL ALLOC (EMT, NELMNT, 'EMT', 'MANEIG')
                  CALL ALLOC (IROW, NELMNT, 'IROW', 'MANEIG')
!                  CALL ALLOC (IENDC, NCF+1, 'IENDC', 'MANEIG')
                  ALLOCATE(IENDC(0:NCF+1), STAT = AllocateStatus)
                  IF (AllocateStatus /= 0)                             &
                    STOP "*** Not enough memory for IENDC***"
                  IOFSET = 0
                  IENDC(0) = 0
                  READ (IMCDF) NCFDUM, ICCUTDUM, MYIDDUM, NPROCSDUM
                  IF (NCF/=NCFDUM .OR. MYID/=MYIDDUM .OR. NPROCSDUM/=NPROCS) &
                     STOP 'maneig:3'
!                  DO I = MYID + 1, NCF, NPROCS
                  DO I = 1, NODENCOLS
                     READ (IMCDF) NELC, ELSTO,                         &
                                  (EMT(IR+IOFSET),IR=1,NELC),          &
                                  (IROW(IR+IOFSET),IR=1,NELC)
                     IF (IROW(NELC+IOFSET).NE.NODECOLS(I)) THEN
                        WRITE(*,*)'Wrong column index....',          &
                 'I,NODECOLS(I), IROW(NELC+IOFSET)=', I,NODECOLS(I), &
                  IROW(NELC+IOFSET)
      STOP 'Unexpected IROW(NELC+IOFSET).NE.NODECOLS(I) in maneig.f ...'
                     ENDIF
                     EMT(NELC+IOFSET) = EMT(NELC+IOFSET) - EAV
!                     DIAG(I) = EMT(NELC+IOFSET)
                     DIAG(NODECOLS(I)) = EMT(NELC+IOFSET)
                     IOFSET = IOFSET + NELC
!                     IENDC(I) = IOFSET
                     IENDC(NODECOLS(I)) = IOFSET
                  END DO
               ELSE
!-----------------------------------------------------------------------
!
! (2.2.2.2) - DVDSON --- Memory, load all, dense
!
!-------------------------------------------------------
                  IF (myid .EQ. 0) THEN
                     WRITE (24,*) ' full matrix stored in in memory'
                     WRITE ( *,*) ' full matrix stored in in memory'
                  END IF

                  IMV = 3
!===================================================================
! CYC: CSFG-version implementation
! Find NDENSE_L, the number of elements on the node (dense form)

                  NDENSE_L = 0
                  !DO i = myid + 1, NCF, nprocs
                  DO I = 1, NODENCOLS
                     !NDENSE_L = NDENSE_L + i
                     NDENSE_L = NDENSE_L + NODECOLS(I)
                  ENDDO

                  CALL ALLOC (EMT, NDENSE_L, 'EMT', 'MANEIG')
                  CALL DINIT (NDENSE_L, 0.0D00, EMT, 1)
                  CALL ALLOC (WORK, NCF, 'WORK', 'MANEIG')
                  CALL ALLOC (IROW, NCF, 'IROW', 'MANEIG')

                  READ (IMCDF) NCFDUM, ICCUTDUM, MYIDDUM, NPROCSDUM
                  IF (NCF/=NCFDUM .OR. MYID/=MYIDDUM .OR. NPROCSDUM/=NPROCS) &
                     STOP 'maneig:4'

                  IOFSET = 0
!                  DO I = MYID + 1, NCF, NPROCS
                  DO I = 1, NODENCOLS
                     READ (IMCDF) NELC, ELSTO,                         &
                                  (WORK(IR),IR=1,NELC),                &
                                  (IROW(IR),IR=1,NELC)
                     ! Check
                     IF (IROW(NELC).NE.NODECOLS(I)) THEN
                      WRITE(*,*)'Wrong column index....',             &
      'I,NODECOLS(I), IROW(NELC)=', I,NODECOLS(I),IROW(NELC)
      STOP 'Unexpected IROW(NELC).NE.NODECOLS(I) in maneig.f Line 531 '
                     ENDIF
                     WORK(NELC) = WORK(NELC) - EAV
!                     DIAG(I) = WORK(NELC)
                     DIAG(NODECOLS(I)) = WORK(NELC)
                     DO IR = 1, NELC
                        EMT(IOFSET+IROW(IR)) = WORK(IR)
                     END DO
!                     IOFSET = IOFSET + I
                     IOFSET = IOFSET + NODECOLS(I)
                  END DO
                  CALL DALLOC (WORK, 'WORK', 'MANEIG')
                  CALL DALLOC (IROW, 'IROW', 'MANEIG')

               ENDIF
! CYC: CSFG-version implementation
!===================================================================
!               ...Memory mode - sparse or dense
!-----------------------------------------------------------------------
!  (2.2.2.3e)          *** E n d   o f   D V D S O N   m e m o r y
!-----------------------------------------------------------------------
            ENDIF
!            ...Disk or Memory
!-----------------------------------------------------------------------
!  (2.2.3e)      *** E n d   o f   D V D S O N
!-----------------------------------------------------------------------

! Make diagonals global, no matter it is disk or memory mode

            CALL gdsummpi (DIAG, NCF)
!            if (myid.eq.0) print *, 'Got diag, imv=',imv
!
!   Allocate storage for workspace; see the header of DVDSON for
!   the expression below; the value of LIM can be reduced to NVECMX
!   plus a smaller number if storage is severely constrained
!
            LIM = MIN(NCF,2*NVECMX + 60)
            LWORK = 2*NCFT*LIM + LIM*LIM*2 + 11*LIM + NVECMX
            CALL ALLOC (WORK, LWORK, 'WORK', 'MANEIG')
            work(1:lwork) = 0.0d0
            LIWORK = 6*LIM + NVECMX
            CALL ALLOC (IWORK, LIWORK, 'IWORK', 'MANEIG')

            CRITE = 1.0D-17
            CRITC = 1.0D-09
            CRITR = 1.0D-09
            ORTHO = MAX(1D-09,CRITR)

!CYC:
            MAXITR = MAX (NVECMX*1000,200)
            !maxitr = MIN (nvect*100, NCF)  ! FROM RSCFVU !!!
            CALL ALLOC (JWORK, LIM,'JWORK', 'MANEIG' )

            CALL ALLOC (EVAL, NVECMX, 'EVAL', 'MANEIG')
            CALL ALLOC (EVEC, NCFT*NVECMX, 'EVEC', 'MANEIG')

            DMUNGO = 10.D99
            CALL DINIT (NVECMX, DMUNGO, EVAL, 1)

!   Compute the eigenpairs in each block
            NVEX = NVECMX
            IF (LDISC) THEN
               MBLOCK = NVEX
            ELSE
               MBLOCK = 1
            ENDIF
            NEND = NCF*NVEX
            IF (NEND < 0 ) THEN
              WRITE(*,*)"NEND Overflow, the calculation is too" &
                         // " large-scale: NCF is too large,"   &
                         // " or too many levels to be resolved."
              STOP "Error, Overflow occurs 2, &
                    out of the ability of the present version!"
            ENDIF
            ILOW = 1
            IHIGH = NVEX
            NIV = NVEX
            NMV = 0
!************************************************************************
!
!   Call Davidson eigensolver
!
            CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

            SELECT CASE (IMV)
            CASE (1)
!******************** sparse and matrix on disk **********************
              IF (myid .EQ. 0) print *, ' Sparse - Disk, iniestsd'
                CALL POSFILE (0, IMCDF, NPOSITION)! was within iniestsd before
              IF (myid .eq. 0) print *, 'Calling iniestsd ',           &
                'IOLPCK=', IOLPCK, ' NCF=', NCF
              CALL INIESTSD (IOLPCK, NCF, MYID, NPROCS, NIV, WORK, IMCDF, EAV)
!             print *, 'Returned from iniestsd '
              if (NCF.gt. IOLPCK) then
                CALL alloc (emt, NCF, 'EMT', 'MANEIGMPI')
                CALL alloc (irow, NCF, 'IROW', 'MANEIGMPI')
                if (myid == 0) print *, &
                  'Calling GDVD to diagonalize the matrix ...'
                starttime = MPI_Wtime()
                CALL GDVD(SPODMV,NCF,LIM,DIAG,ILOW,IHIGH,&
                  JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR,   &
                  WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,               &
                  NMV,IERR)
                endtime = MPI_Wtime()
                if (myid.eq.0)                                        &
                   print *, 'Time used by gdvd (m) = ',               &
                          (endtime-starttime)/60.0,                   &
                          '  NLOOPS = ', NLOOPS
                CALL dalloc (emt, 'EMT', 'MANEIG')
                CALL dalloc (irow, 'IROW', 'MANEIG')
              end if
            CASE (2)
!******************** sparse and matrix in memory ********************
              IF (myid .EQ. 0) print *, 'Sparse - Memory, iniestmpi'
              CALL INIESTMPI (IOLPCK, NCF,NIV,WORK,EMT,IENDC,IROW)
!              IF (myid .EQ. 0) print *, 'Exit from iniestmpi ......'
              if (NCF.gt. IOLPCK) then
                if (myid == 0) print *, &
                  'Calling GDVD to diagonalize the matrix ...'
                starttime = MPI_Wtime()
                CALL GDVD (SPICMVsym,NCF,LIM,DIAG,ILOW,IHIGH,       &
                 JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR,  &
                 WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,              &
                 NMV,IERR)
                endtime = MPI_Wtime()
              endif

              if (myid.eq.0)                                        &
                 print *, 'Time used by gdvd (m) = ',               &
                          (endtime-starttime)/60.0,                 &
                          '  NLOOPS = ', NLOOPS

              CALL DALLOC (EMT, 'EMT', 'MANEIG')
              CALL DALLOC (IROW, 'IROW', 'MANEIG')
              CALL DALLOC (IENDC, 'IENDC', 'MANEIG')
            CASE (3)
!*************************** dense and in memory **********************
              IF (myid .EQ. 0) print *, ' Dense - Memory, iniestdm'
              CALL INIESTDM (IOLPCK,NCF,NIV,WORK,EMT)
              if (NCF.gt. IOLPCK) then
                if (myid == 0) print *, &
                  'Calling GDVD to diagonalize the matrix ...'
                starttime = MPI_Wtime()
                CALL GDVD (DNICMV,NCF,LIM,DIAG,ILOW,IHIGH,  &
                    JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR, &
                    WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,             &
                    NMV,IERR)
                endtime = MPI_Wtime()
              end if

              if (myid.eq.0)                                        &
                 print *, 'Time used by gdvd (m) = ',               &
                          (endtime-starttime)/60.0,                 &
                          '  NLOOPS = ', NLOOPS

              CALL DALLOC (EMT, 'EMT', 'MANEIG')
            END SELECT
!************************************************************************
            CALL DALLOC (DIAG, 'DIAG', 'MANEIG')
            CALL DALLOC (IWORK, 'IWORK', 'MANEIG')
            CALL DALLOC (JWORK, 'JWORK', 'MANEIG')
            IF (myid .EQ. 0) THEN
               WRITE (24, *) ' ', NLOOPS, ' iterations;'
               WRITE (24, *) ' ', NMV, ' matrix-vector multiplies.'
            ENDIF
            IF (IERR /= 0) THEN
               WRITE (ISTDE, *) 'MANEIG: Returned from DVDSON with'
               WRITE (ISTDE, *) ' IERR = ', IERR, '.'
               CALL stopmpi ('maneig: DVDSON wrong', myid)
            ENDIF

!   Put the eigenpairs in order, overwriting as necessary
            CALL DCOPY (NVEX, WORK(NEND+1), 1, EVAL, 1)
            CALL DCOPY (NCFt*NVEX, WORK(1), 1, EVEC, 1)
            CALL DALLOC (WORK, 'WORK', 'MANEIG')

!   Rearrange and reallocate storage for the eigenpairs
!   as necessary

            IF (NVEC < NVECMX) THEN
               CALL ALLOC (IWORK, NVECMX, 'IWORK', 'MANEIG')
               DO I = 1, NVECMX
                  IWORK(I) = I
               END DO
               DO I = 1, NVEC
                  IOFSET = IVEC(I)
                  LOC = IWORK(I)
                  IF (IOFSET == LOC) CYCLE
                  CALL DSWAP (1, EVAL(IOFSET), 1, EVAL(I), 1)
                  IWORK(I) = IWORK(IOFSET)
                  IWORK(IOFSET) = LOC
                  IOFSET = NCF*(IOFSET - 1)
                  LOC = NCF*(I - 1)
                  CALL DSWAP (NCF, EVEC(IOFSET+1), 1, EVEC(LOC+1), 1)
               END DO
               CALL DALLOC (IWORK, 'IWORK', 'MANEIG')
               CALL RALLOC (EVAL,  NVEC, 'EVAL', 'MANEIG')
               CALL RALLOC (EVEC, NCF*NVEC, 'EVEC', 'MANEIG' )
            ENDIF

         ENDIF
! (2.3e)              *** E N D   O F    N O N - T R I V I A L   C A S E

      ENDIF
! (3e)               *** E N D   O F    A L L
      DEALLOCATE(NODECOLS)

!--------------------------------------------------------------------
! Only the following quantities are needed after this routine is
! finished:
!   eval(), evec(), iatjpo, iaspar
!--------------------------------------------------------------------
!
!   Clean up eigenvectors; determine their J/P values
!
      DO J = 1, NVEC

!         Find the dominant component of each eigenvector

         IOFSET = (J - 1)*NCF
         AMAX = 0.D0
         DO I = 1, NCF
            WA = ABS(EVEC(I+IOFSET))
            IF (WA <= AMAX) CYCLE
            AMAX = WA
            IA = I
         END DO

!          Find the angular momentum and parity of the dominant component
! Removed into matrix.f
!         IATJPO = ITJPO(IA)
!         IASPAR = ISPAR(IA)

!          Change sign of eigenvactor if dominant component is negative

         IF (EVEC(IA+IOFSET) >= 0.D0) CYCLE
         EVEC(1+IOFSET:NCF+IOFSET) = -EVEC(1+IOFSET:NCF+IOFSET)
      END DO

      RETURN
      END SUBROUTINE MANEIG
