!=======================================================================
!                                                                      *
! Pack LABEL from IA, IB, IC, ID                                       *
!                                                                      *
!=======================================================================
      FUNCTION LABPACK(IA0, IB0, IC0, ID0)

      USE parameter_def,   ONLY: KEYORB

      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,  INTENT(IN)  :: IA0, IB0, IC0, ID0
      INTEGER               :: LABPACK
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: KEYSQ = KEYORB*KEYORB
      INTEGER :: IA, IB, IC, ID
      INTEGER :: LAB, ISWAP, LAC, LBD
!-----------------------------------------------
      IA = IA0; IB = IB0; IC = IC0; ID = ID0
      IF (IA > IC) THEN
        ISWAP = IC
        IC = IA
        IA = ISWAP
      ENDIF
      IF (IB > ID) THEN
        ISWAP = ID
        ID = IB
        IB = ISWAP
      ENDIF
      IF (IA == IB .AND. IC > ID) THEN
        ISWAP = IC
        IC = ID
        ID = ISWAP
      ENDIF
      LAC = IA * KEY + IC
      LBD = IB * KEY + ID
      IF (LAC .LT. LBD) THEN
        LABPACK = LAC * KEYSQ + LBD
      ELSE
        LABPACK = LBD * KEYSQ + LAC
      ENDIF
      
      RETURN
      END FUNCTION LABPACK
