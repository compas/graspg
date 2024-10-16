      SUBROUTINE TESTFILES(FDNAME)
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*),INTENT(IN) :: FDNAME
!
      INTEGER :: ISTATE

      ISTATE = ACCESS (TRIM(FDNAME), ' ')
      IF (ISTATE .NE. 0) THEN
        WRITE(*, '(A31, A)')"Directory/File does not exist: ", TRIM(FDNAME)
        STOP "ERROR! Directory/File does not exist ..."
      ENDIF

      RETURN
      END
