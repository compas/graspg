!***********************************************************************
!                                                                      *
      SUBROUTINE SETLABORB (FNAME)
!                                                                      *
!   Set nonsym, i.e., the number of labeling orbitals.                *
!                                                                      *
!   Written by Chongyang Chen                               Nov 2023   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE IOUNIT_C,      ONLY: ISTDE
      USE SYMEXPAND_MOD, ONLY: nonsyminf, strlaborb
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(*)  :: FNAME
!-----------------------------------------------

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IERR
!-----------------------------------------------

      OPEN (2023,FILE=TRIM(FNAME), STATUS='OLD', &
               FORM='FORMATTED',IOSTAT=IERR)
      IF (IERR == 0) THEN
        READ(2023, *) nonsyminf
        READ(2023, '(A)')strlaborb
        CLOSE(2023)
!PERJ        WRITE(*,*)'nonsyminf =', nonsyminf
      ELSE
        WRITE (ISTDE, *) TRIM(FNAME) // &
                         ' not exists, Stop Error!'
        STOP
      ENDIF
!
      RETURN
      END SUBROUTINE SETLABORB
