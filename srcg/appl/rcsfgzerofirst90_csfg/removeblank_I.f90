      MODULE removeblank_I
      INTERFACE
      SUBROUTINE removeblank (linein, m, lineout, n, ICSF)
      CHARACTER(LEN=256), INTENT(IN) :: LINEIN 
      INTEGER,  INTENT(IN) :: M, ICSF    
      CHARACTER(LEN=256), INTENT(OUT) :: LINEOUT
      INTEGER,  INTENT(OUT) :: N
      END SUBROUTINE
      END INTERFACE
      END MODULE
