      MODULE nine_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90 4.3E  10:07:58  11/16/01
      SUBROUTINE nine (J1, J2, J3, L1, L2, L3, K1, K2, K3, I, IN, AA)
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER, INTENT(IN) :: J1
      INTEGER, INTENT(IN) :: J2
      INTEGER, INTENT(IN) :: J3
      INTEGER, INTENT(IN) :: L1
      INTEGER, INTENT(IN) :: L2
      INTEGER, INTENT(IN) :: L3
      INTEGER, INTENT(IN) :: K1
      INTEGER, INTENT(IN) :: K2
      INTEGER, INTENT(IN) :: K3
      INTEGER, INTENT(IN) :: I
      INTEGER, INTENT(OUT) :: IN
      REAL(DOUBLE), INTENT(INOUT) :: AA
!VAST.../CONSTS/ ZERO(IN)
!VAST...Calls: ITTK, NINE0, SIXJ
      END SUBROUTINE
      END INTERFACE
      END MODULE
