      MODULE lasa1_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
      SUBROUTINE lasa1 (FIL, RAD, POP, SKAL, SLUT)
      integer, INTENT(IN) :: FIL
      character (LEN = 200), INTENT(OUT) :: RAD
      integer, DIMENSION(25,0:10,0:1) :: POP
      integer :: SKAL
      logical, INTENT(INOUT) :: SLUT
      END SUBROUTINE
      END INTERFACE
      END MODULE
