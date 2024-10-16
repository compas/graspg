!***********************************************************************
!                                                                      *
      SUBROUTINE onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
!                                                                      *
!   This subroutine calculates the one electron matrix elements.       *         
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
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
      INTEGER      :: IA, IB
      REAL(DOUBLE) :: TCOEFF, EMTTMP, ATWINV
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: TEGRAL
!-----------------------------------------------------------------------
!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
!
      EMTTMP=0.0D0
      CALL IABINT(IA,IB,TEGRAL)
      EMTTMP = EMTTMP + TEGRAL*TCOEFF

      IF (LNMS) THEN
         CALL KEINT(IA,IB,TEGRAL)
         EMTTMP = EMTTMP + TEGRAL*ATWINV*TCOEFF
      ENDIF

      IF (LVP) THEN
         CALL VPINT(IA,IB,TEGRAL)
         EMTTMP = EMTTMP + TEGRAL*TCOEFF
      ENDIF

      RETURN
      END SUBROUTINE onebody_DC_MS_VP
