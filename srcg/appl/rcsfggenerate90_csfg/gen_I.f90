      MODULE gen_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
      SUBROUTINE gen (ANSATS, POSN, POSL, SKAL, CF, FIRST, MINJ, MAXJ, PAR)
      integer, DIMENSION(25,0:10,0:1), INTENT(IN) :: ANSATS
      integer, DIMENSION(220), INTENT(IN) :: POSN
      integer, DIMENSION(220), INTENT(IN) :: POSL
      integer, INTENT(IN) :: SKAL
      integer, INTENT(INOUT) :: CF
      logical, INTENT(IN) :: FIRST
      integer, INTENT(IN) :: MINJ
      integer, INTENT(IN) :: MAXJ
      integer :: PAR
      END SUBROUTINE
      END INTERFACE
      END MODULE
