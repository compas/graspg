      MODULE locateall_I
      INTERFACE
      SUBROUTINE locateall (ILABWORK,NUMB,NUME,IARR,LFOUND,LOC)
      USE vast_kind_param, ONLY: LONG
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ILABWORK, NUMB, NUME
      INTEGER :: IARR(*)
      LOGICAL, INTENT(OUT) :: LFOUND
      INTEGER, INTENT(OUT) :: LOC
      END SUBROUTINE
      END INTERFACE
      END MODULE
