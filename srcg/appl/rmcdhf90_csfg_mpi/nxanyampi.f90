!***********************************************************************
      SUBROUTINE NXANYAmpi
!                                                                      *
!   Called by SCFMPI subroutine, make every node has the same arrays   *
!   NXCOFopt, NYCOFopt, NXAopt and NYAopt, then XPOT and YPOT
!   subroutines could be parrallelized. 

!   Written by Chongyang Chen               Last update:    Jul 2022   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,    ONLY: KEYORB
      USE csfg_memory_man
      USE def_C
      USE mpi_C
      USE orb_C, ONLY: NP,NH,NW
      USE csfg_scf_C
      use, intrinsic :: iso_fortran_env
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB*KEYORB*KEYORB
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I_MPI, IORB, IXY
      INTEGER :: I, J, K, NUM, NUMT, NTOT
      INTEGER :: ITHISB, ITHISC, ICA, ICB, ICC
      INTEGER :: NXYDIM0, NXYDMAX, NXYD0
      INTEGER, DIMENSION(:), POINTER :: NXYDBUF
      INTEGER, DIMENSION(:), POINTER :: NXYABUF 
      INTEGER, DIMENSION(:), POINTER :: NXYATA, NXYATB, NXYATC 
      LOGICAL :: FLAGIORB(NW)
!-----------------------------------------------------------------------
      I_MPI = MPI_INTEGER
!
! Consolidate NXCOFOPT, NXAOPT, NYCOFOPT, and NYAOPT from all nodes
!
!      IF (NPROCS == 1 ) RETURN
      IF (NPROCS == 1 ) GOTO 100

      ! IXY = 1: Gather NXCOFOPT and NXAOPT
      ! IXY = 2: Gather NYCOFOPT and NYAOPT
      DO IXY = 1, 2 
        ! The maximum value of NXCOFOPT / NYCOFOPT at each node 
        CALL ALLOC(NXYDBUF, nprocs, 'NXYDBUF','NXANYAmpi')
        IF (IXY == 1 ) THEN 
          NXYDIM0 = MAXVAL(NXCOFOPT)
        ELSE
          NXYDIM0 = MAXVAL(NYCOFOPT)
        ENDIF 
        ! Gather the maximum value of NXCOFOPT / NYCOFOPT 
        CALL MPI_GATHER(NXYDIM0, 1, I_MPI, NXYDBUF,            &
             1, I_MPI, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_Bcast (NXYDBUF, nprocs, I_MPI, 0,             &
                              MPI_COMM_WORLD, ierr)
        ! The maximum value of NXCOFOPT / NYCOFOPT over all nodes
        NXYDMAX = MAXVAL(NXYDBUF)
        IF (NXYDMAX.LE.0) RETURN
        IF (MYID == 0 .and. IXY == 1) &
          WRITE(*,'(1X,A,I12)')'MAX NXCOF elements before merge:', NXYDMAX
        IF (MYID == 0 .and. IXY == 2) &
          WRITE(*,'(1X,A,I12)')'MAX NYCOF elements before merge:', NXYDMAX
        ! Working arrays
        NTOT = NXYDMAX*NPROCS
        IF (NTOT .LT. 0) THEN
          WRITE(*,*)"Eorror, NTOT of nxanya.f90 exceeds 2^31 ..."
          STOP "Eorror, NTOT of nxanya.f90 exceeds 2^31 ..."
        ENDIF
        CALL ALLOC (NXYABUF, NTOT, 'NXYABUF', 'NXANYAmpi')
        IF (MYID == 0) THEN
          CALL ALLOC (NXYATA, NTOT, 'NXYATA', 'NXANYAmpi') 
          CALL ALLOC (NXYATB, NTOT, 'NXYATB', 'NXANYAmpi')
          CALL ALLOC (NXYATC, NXYDMAX, 'NXYATC', 'NXANYAmpi')
        ENDIF

        IF (IXY == 1) THEN 
          NXYD0 = SIZE(NXAOPT,DIM=1)
        ELSE
          NXYD0 = SIZE(NYAOPT,DIM=1)
        ENDIF 

        DO J = 1, NW
          IF (LSKIPPOT(J)) CYCLE
          IF (IXY == 1) THEN
            NXYDIM0 = NXCOFOPT(J)
          ELSE
            NXYDIM0 = NYCOFOPT(J)
          ENDIF
          CALL MPI_GATHER(NXYDIM0, 1, I_MPI, NXYDBUF,             &
             1, I_MPI, 0, MPI_COMM_WORLD, ierr)
          CALL MPI_Bcast (NXYDBUF, NPROCS, I_MPI, 0,             &
                              MPI_COMM_WORLD, ierr)
          NXYDIM0 = MAXVAL(NXYDBUF)
          ! Add the needed data to make sure NXAOPT / NYAOPT of every 
          ! node has at least NXYDIM0 elements.
          IF (IXY == 1) THEN
            IF (NXYDIM0 > SIZE(NXAOPT,DIM=1))                      &
              CALL RALLOC (NXAOPT, NXYDIM0, NW, 'NXAOPT', 'NXANYA')
          ELSE
            IF (NXYDIM0 > SIZE(NYAOPT,DIM=1))                      &
              CALL RALLOC (NYAOPT, NXYDIM0, NW, 'NYAOPT', 'NXANYA')
          ENDIF

          ! Gather the LABELs from all nodes
          IF (IXY == 1) THEN
            CALL MPI_GATHER(NXAOPT(1:NXYDIM0,J), NXYDIM0, I_MPI, NXYABUF,  &
                          NXYDIM0, I_MPI, 0, MPI_COMM_WORLD, ierr)
          ELSE
            CALL MPI_GATHER(NYAOPT(1:NXYDIM0,J), NXYDIM0, I_MPI, NXYABUF,  &
                          NXYDIM0, I_MPI, 0, MPI_COMM_WORLD, ierr)
          ENDIF

          ! Perform the mergement at node 0 
          IF (MYID.EQ.0) THEN
            NUMT = NXYDBUF(1)
            ! Initialize arrays NXYATA and NXYATB
            NXYATA(1:NUMT) = NXYABUF(1:NUMT)
            NXYATB(1:NUMT) = NXYABUF(1:NUMT)

            DO I = 1, NPROCS-1
              ! Here, I = MYID
              NUM = NXYDBUF(I+1)
              ! data of node i, (counting from node 0) 
              ! Initialize array NXYATC
              NXYATC(1:NUM) = NXYABUF(I*NXYDIM0+1:I*NXYDIM0+NUM)
              ICA = 0
              ICB = 1
              ICC = 1
              DO WHILE (ICB<=NUMT .AND. ICC<=NUM) 
                ITHISB = NXYATB(ICB)
                ITHISC = NXYATC(ICC)
                IF (ITHISB == ITHISC) THEN
                  ICA = ICA + 1
                  NXYATA(ICA) = ITHISC
                  ICB = ICB + 1
                  ICC = ICC + 1
                ELSEIF (ITHISB < ITHISC) THEN
                  ICA = ICA + 1
                  NXYATA(ICA) = ITHISB
                  ICB = ICB + 1
                ELSE
                  ICA = ICA + 1
                  NXYATA(ICA) = ITHISC
                  ICC = ICC + 1
                ENDIF
              ENDDO

              IF (ICB <= NUMT) THEN
                DO K = 1, NUMT - ICB + 1
                   ICA = ICA + 1
                   NXYATA(ICA) = NXYATB(ICB+K-1)
                ENDDO 
              ENDIF
              IF (ICC <= NUM) THEN
                DO K = 1, NUM - ICC + 1
                  ICA = ICA + 1
                  NXYATA(ICA) = NXYATC(ICC+K-1)
                ENDDO
              ENDIF
              ! Merged data of node i
              NUMT = ICA
              NXYATB(1:NUMT) = NXYATA(1:NUMT) 
            ENDDO

            IF (IXY == 1) THEN
              NXCOFOPT(J) = NUMT
              IF (NUMT > NXYD0) THEN
                NXYD0 = 1.25*NUMT
                CALL RALLOC(NXAOPT, NXYD0, NW, 'NXAOPT','NXANYAmpi')
              ENDIF
              NXAOPT(1:NUMT,J) = NXYATA(1:NUMT)
            ELSE
              NYCOFOPT(J) = NUMT
              IF (NUMT > NXYD0) THEN
                NXYD0 = 1.25*NUMT
                CALL RALLOC(NYAOPT, NXYD0, NW, 'NYAOPT','NXANYAmpi')
              ENDIF
              NYAOPT(1:NUMT,J) = NXYATA(1:NUMT)
            ENDIF
          ENDIF
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
        ENDDO

        IF (IXY == 1) THEN
          CALL MPI_Bcast (NXCOFOPT, NW, MPI_INTEGER, 0,              &
                                    MPI_COMM_WORLD, ierr)
          NXYDMAX = MAXVAL(NXCOFOPT)
          CALL RALLOC(NXAOPT, NXYDMAX, NW, 'NXAOPT', 'NXANYAmpi')
          CALL MPI_Bcast (NXAOPT, NXYDMAX*NW, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          CALL RALLOC (XAOPT, NXYDMAX, NW, 'XAOPT','NXANYAmpi')

          ! Update JXIPOS(0:NNNW, NNNW)
          DO J = 1, NW
             IF (LSKIPPOT(J)) CYCLE
             FLAGIORB = .FALSE.
             DO K = 1, NXCOFOPT(J)
               IORB = NXAOPT(K,J)/KEY
               JXIPOS(IORB,J) = K
               FLAGIORB(IORB) = .TRUE.
             ENDDO
             DO K = 1, NW
               IF (.NOT.FLAGIORB(K)) JXIPOS(K,J) = JXIPOS(K-1,J) 
             ENDDO
          ENDDO
        ELSE
          CALL MPI_Bcast (NYCOFOPT, NW, MPI_INTEGER, 0,              &
                                    MPI_COMM_WORLD, ierr)
          NXYDMAX = MAXVAL(NYCOFOPT)
          CALL RALLOC(NYAOPT, NXYDMAX, NW, 'NYAOPT', 'NXANYAmpi')
          CALL MPI_Bcast (NYAOPT, NXYDMAX*NW, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          CALL RALLOC (YAOPT, NXYDMAX, NW, 'YAOPT','NXANYAmpi')
        ENDIF

        CALL DALLOC (NXYDBUF, 'NXYDBUF','NXANYAmpi')
        CALL DALLOC (NXYABUF, 'NXYABUF', 'NXANYAmpi')
        IF (MYID == 0) THEN
          CALL DALLOC (NXYATA,  'NXYATA', 'NXANYAmpi')
          CALL DALLOC (NXYATB,  'NXYATB', 'NXANYAmpi')
          CALL DALLOC (NXYATC,  'NXYATC', 'NXANYAmpi')
        ENDIF
      ENDDO

100   CONTINUE
!      IF (MYID.EQ.0) THEN
!        DO I = 1, NW
!          IF (LSKIPPOT(I)) CYCLE
!          WRITE(*,'(I2,A2,2I10)')NP(I),NH(I),NXCOFOPT(I), NYCOFOPT(I)
!        ENDDO
!        WRITE(*,*)
!      ENDIF

! CYC 2024/03
! Output NXCOFOPT, NXAOPT and JXIPOS for check
!      IF (MYID == 0) THEN
!        DO I = 1, NW
!          IF (LSKIPPOT(I)) CYCLE
!          WRITE(1948,'(I2,A2,I10)')NP(I),NH(I),NXCOFOPT(I)
!          WRITE(1948,'(10I12)')JXIPOS(1:NW,I)
!          WRITE(1948,*)
!          !WRITE(1948,'(10I12)')NXAOPT(1:NXCOFOPT(I),I)
!        ENDDO 
!        CLOSE(1948)
!      ENDIF
!      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      flush output_unit

!      IF (MYID.EQ.0) WRITE(*,*)"Finished NXANYAmpi ..."
      RETURN
      END SUBROUTINE NXANYAmpi
