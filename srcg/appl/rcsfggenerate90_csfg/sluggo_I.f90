      MODULE sluggo_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
      SUBROUTINE sluggo (I, J, VARMAX, VARUPP, VARNED, ANSATS, ORG, LOCK, LOW&
         , START, STOPP)
      integer, INTENT(IN) :: I
      integer, INTENT(IN) :: J
      integer, INTENT(IN) :: VARMAX
      integer, DIMENSION(25,0:10), INTENT(INOUT) :: VARUPP
      integer, DIMENSION(25,0:10), INTENT(INOUT) :: VARNED
      integer, DIMENSION(25,0:10,0:1), INTENT(IN) :: ANSATS
      integer, DIMENSION(25,0:10), INTENT(IN) :: ORG
      logical, INTENT(IN) :: LOCK
      integer, DIMENSION(25,0:10), INTENT(IN) :: LOW
      integer, INTENT(OUT) :: START
      integer, INTENT(OUT) :: STOPP
      END SUBROUTINE
      END INTERFACE
      END MODULE