!***********************************************************************
!                                                                      *
      SUBROUTINE twobody_DC(TEGRAL)
!                                                                      *
!   This subroutine calculates the two-body DC, MS and VP contributions* 
!                                                                      *
!   Modified by CHONG-YANG CHEN                              JUNE 2020 *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE symmatrix_mod, ONLY : LABV
!-----------------------------------------------
!   C O M M O N  B L O C K S
!-----------------------------------------------
      USE decide_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) :: TEGRAL
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: IA, IB, IC, ID, K, ISW
!-----------------------------------------------
!   Accumulate the contributions from the two-electron
!   Coulomb operator.

      TEGRAL = 0.0d0
      IA = LABV(1)
      IB = LABV(2)
      IC = LABV(3)
      ID = LABV(4)
      K  = LABV(5)

      !IF (IA > IC) THEN
      !  ISW = IA
      !  IA  = IC
      !  IC  = ISW
      !ENDIF
      !IF (IB > ID) THEN
      !  ISW = IB
      !  IB  = ID
      !  ID  = ISW
      !ENDIF

      CALL RKINTC (IA, IB, IC, ID, K, TEGRAL)
      RETURN
      END SUBROUTINE twobody_DC
