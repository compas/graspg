      MODULE rmeajj9_I
      INTERFACE
!
      SUBROUTINE RMEAJJ9(IT,LQ,J,ITS,LQS,J1S,COEF)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER,      INTENT(IN)  :: IT, LQ, J, ITS, LQS, J1S
      REAL(DOUBLE), INTENT(OUT) :: COEF
      END SUBROUTINE
      END INTERFACE
      END MODULE
