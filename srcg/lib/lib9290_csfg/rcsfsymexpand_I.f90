      MODULE rcsfsymexpand_I
      INTERFACE
      SUBROUTINE rcsfsymexpand (NAME, FILE_RL, FILE_OUT, IFLAG)
      CHARACTER(*), INTENT(IN) :: NAME
      CHARACTER(*), INTENT(IN) :: FILE_RL
      CHARACTER(*), INTENT(IN) :: FILE_OUT
      INTEGER,      INTENT(IN) :: IFLAG
      END SUBROUTINE
      END INTERFACE
      END MODULE
