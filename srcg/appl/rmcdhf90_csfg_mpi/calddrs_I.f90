      MODULE calddrs_I
      INTERFACE
      SUBROUTINE calddrs(EOL, JBLOCK, ICG, IRG)
      USE vast_kind_param, ONLY: DOUBLE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ICG, IRG, JBLOCK
      LOGICAL, INTENT(IN) :: EOL

      END SUBROUTINE
      END INTERFACE
      END MODULE
