!***********************************************************************
!                                                                      *
      SUBROUTINE twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
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
      REAL(DOUBLE) :: ATWINV, VCOEFF, EMTTMP
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: TGRL1, TGRL2, TEGRAL
!-----------------------------------------------
!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation; the latter
!   is computed first because the orbital indices may be
!   permuted by RKINTC

      EMTTMP = 0.0d0
      IF (LSMS) THEN
        IF (LABV(5) .EQ. 1) THEN
          CALL VINT (LABV(1), LABV(3), TGRL1)
          CALL VINT (LABV(2), LABV(4), TGRL2)
          EMTTMP = EMTTMP - TGRL1*TGRL2*ATWINV*VCOEFF
        ENDIF
      ENDIF
      CALL RKINTC (LABV(1), LABV(2), LABV(3), LABV(4), LABV(5), TEGRAL)
     
      EMTTMP = EMTTMP + TEGRAL*VCOEFF

      RETURN
      END SUBROUTINE twobody_DC_SMS
