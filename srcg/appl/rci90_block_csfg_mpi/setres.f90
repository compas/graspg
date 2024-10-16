!***********************************************************************
!                                                                      *
      SUBROUTINE SETRES(ISOFILE, RWFFILE, IDBLK)
!                                                                      *
!   Open, check, load data from the  .res  file.                       *
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, OPENFL.                                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
!   Modified by Xinghong                  Last revision: 23 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE memory_man
      USE default_C
!CYC      USE where_C
      USE hblock_C
      USE iccu_C
      USE iounit_C
      USE mpi_C

      USE def_C, ONLY : NELEC
      USE orb_C
      use symmatrix_mod,       ONLY: IMCDF, NCFTOTTr
      use symmatrix_restart_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE openfl_I
      USE lodres_I
      USE getcid_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: ISOFILE*(*)
      CHARACTER  :: RWFFILE*(*)
      CHARACTER(LEN=8), DIMENSION(*)  :: IDBLK
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER*11, PARAMETER :: FORM = 'UNFORMATTED'
      CHARACTER*7, PARAMETER :: STATUS = 'UNKNOWN'
      CHARACTER*6, PARAMETER :: RESTITLE = 'R92RES'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IOS, I, J, NBLOCK_OLD
      LOGICAL :: FOUND
      CHARACTER :: R92RES*6, DEFNAM*11, IDSTRING*3, G92MIX*6

! CYC: Restart Mode
      CHARACTER*13 :: CINFNCOL
      REAL(DOUBLE) :: EAV_OLD
      INTEGER  :: NELEC_OLD, NB_OLD, NCFTOT_OLD, NW_OLD
      INTEGER  :: NVECSIZ_OLD, NVECTOT_OLD
!-----------------------------------------------
!
! Compose the "rciXXX.res" file names (for each node)
!
      WRITE (idstring, '(I3.3)') myid
      DEFNAM = 'rci' // idstring // '.res'

! CYC: Restart Mode
! Compose the "ncolsXXX.txt" file names (for each node), storing/reload
! the total number of NCols (NODENCOLS) and the detailed COL-number
! (NODECOLS[])
      CINFNCOL = 'ncols' // idstring // '.txt'

!
! Ask if this is a restart
!
      IF (myid .EQ. 0) THEN
         IF (NDEF /= 0) THEN
            WRITE (ISTDE, *) 'Restarting RCI90 ?'
            RESTRT = GETYN()
         ELSE
            RESTRT = .FALSE.
         ENDIF
!         IF (RESTRT) THEN
!           WRITE(734,'(a)') 'y            ! Restarting RCI90 ?'
!         ELSE
!           WRITE(734,'(a)') 'n            ! Restarting RCI90 ?'
!         END IF
      END IF
      CALL MPI_Bcast (RESTRT, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

!
! Do some settings and checks
!
      IF (RESTRT) THEN
!         ...Restart, make sure file exist
         INQUIRE(FILE=DEFNAM, EXIST=FOUND)
         IF (.NOT. FOUND) THEN
            CALL stopmpi ('setres: .res not exist', myid)
         ENDIF

         INQUIRE (FILE = CINFNCOL, EXIST = FOUND)
         IF (.NOT. FOUND) THEN
            CALL stopmpi ('setres: ncolsXXX.txt does not exist', myid)
         ENDIF
      ENDIF
!
! Open the .res file
!
      CALL OPENFL (IMCDF, DEFNAM, FORM, STATUS, IERR)
      IF (IERR .NE. 0) THEN
         CALL stopmpi ('setres: Error openning .res file', myid)
      ENDIF

      CALL OPENFL (IFNCOL, CINFNCOL, 'FORMATTED', STATUS, IERR)
      IF (IERR .NE. 0) THEN
         CALL stopmpi ('setres: Error openning ncolsXXX.txt file', myid)
      ENDIF

!-----------------------------------------------------------------------
! CYC:
! Try to read old.cm
!
      LOLDCM = .FALSE.
      IF (MYID == 0) THEN
        INQUIRE(FILE=TRIM(FOLDCM), EXIST=FOUND)
        IF (FOUND) THEN
          OPEN(1457,FILE=TRIM(FOLDCM),FORM='UNFORMATTED',STATUS ='OLD')
          READ (1457, IOSTAT=IOS) G92MIX
          IF (IOS/=0 .OR. G92MIX/='G92MIX') THEN
            WRITE (ISTDE, *) 'Not a GRASP92 MIXing Coefficients File;'
            Close(1457)
            STOP 'ERROR! Not a GRASP92 MIXing Coefficients File ...'
          ENDIF
          READ (1457) NELEC_OLD, NCFTOT_OLD, NW_OLD, NVECTOT_OLD, &
               NVECSIZ_OLD, NBLOCK_OLD
          IF (NBLOCK_OLD .GT. 1) Then
                Write(*,*)"NBlock_OLD > 1, Error!!!"
                STOP
          ENDIF
          READ (1457) NB_OLD, NCFBLK_OLD, NEVBLK_OLD, IATJPO_OLD, &
                      IASPAR_OLD
          IF (NELEC_OLD/=NELEC .OR. NW_OLD/=NW .OR. NCFBLK_OLD/=NCFTOTTr) &
            CALL STOPMPI &
                 ('lodres: NELEC/NCF/NW does not match for old.cm',myid)
          ! Allocate the arrays neede for this block
          ALLOCATE (EVAL_OLD(NEVBLK_OLD))
          ALLOCATE (EVEC_OLD(NCFBLK_OLD*NEVBLK_OLD))
          ALLOCATE (IVEC_OLD(NEVBLK_OLD))
          EVEC_OLD(:NEVBLK_OLD*NCFBLK_OLD) = 0.D0

          READ (1457) (IVEC_OLD(I),I=1,NEVBLK_OLD)
          READ (1457) EAV_OLD, (EVAL_OLD(I),I=1,NEVBLK_OLD)
          READ (1457)((EVEC_OLD(I+(J-1)*NCFBLK_OLD),I=1,NCFBLK_OLD), &
                       J=1,NEVBLK_OLD)
          Close(1457)
        ENDIF
      ENDIF

      CALL MPI_Bcast (FOUND, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      IF (FOUND) THEN
         LOLDCM = .TRUE.
         CALL MPI_Bcast (NEVBLK_OLD,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_Bcast (NCFBLK_OLD,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_Bcast (IATJPO_OLD,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_Bcast (IASPAR_OLD,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         IF (myid .NE. 0) THEN
           ALLOCATE (EVAL_OLD(NEVBLK_OLD))
           ALLOCATE (EVEC_OLD(NCFBLK_OLD*NEVBLK_OLD))
           !ALLOCATE (IVEC_OLD(NEVBLK_OLD))
         ENDIF

         CALL MPI_Bcast (EVAL_OLD, NEVBLK_OLD, MPI_DOUBLE_PRECISION, 0,&
                         MPI_COMM_WORLD, ierr)
         CALL MPI_Bcast (EVEC_OLD, NEVBLK_OLD*NCFBLK_OLD, &
                         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
      ENDIF
!
! End reload old.cm
!-----------------------------------------------------------------------
!
! If restart, load the contents. Otherwise generate them via getcid
!
! But first of all, iccutblk() is needed in both cases
!
      CALL ALLOC (ICCUTBLK, NBLOCK, 'ICCUTBLK', 'SETRES')

      IF (RESTRT) THEN
!        ...Check the signature of the file
         READ (IMCDF, IOSTAT=IOS) R92RES
         IF (IOS/=0 .OR. R92RES/=RESTITLE) THEN
            CLOSE(IMCDF)
            CALL stopmpi ('setres: Not RCI92 .res file', myid)
         ENDIF

!         ...Read and check restart information
         CALL LODRES

      ELSE

!         ...Write the file header
!         ...Generate the first part of the .res file
         WRITE (IMCDF) RESTITLE
         CALL GETCID (ISOFILE, RWFFILE, IDBLK)

      ENDIF

      RETURN
      END SUBROUTINE SETRES
