
!======================================================================
! Written by Chongyang Chen, Fudan University, Shanghai, Feb 2024
!
! Calculate d_{rs}, r, s are the CSFs generated by the ICG-IRG CSFG pair
! 
!======================================================================
      SUBROUTINE calddrs(EOL, JBLOCK, ICG, IRG)
      USE vast_kind_param, ONLY: DOUBLE

      use mcp_C,          ONLY: LFORDR
      use symmatrix_mod,  ONLY: LICCUT
      use csfg_decide_C,  ONLY: LSYMD
      use symexpand_mod,  ONLY: MAP1
      use csfg_scf_C,     ONLY: DDRS

      USE dsubrs_I

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ICG, IRG, JBLOCK
      LOGICAL, INTENT(IN) :: EOL
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: I, J, M, N, IR0, IC0, ICEND, IREND
!-----------------------------------------------
      IC0 = 0
      IR0 = 0
      IF (ICG > 1) IC0 = MAP1(ICG-1, JBLOCK)
      IF (IRG > 1) IR0 = MAP1(IRG-1, JBLOCK)

      ICEND = MAP1(ICG, JBLOCK)
      IREND = MAP1(IRG, JBLOCK)
      
      DDRS(1:IREND-IR0, 1:ICEND-IC0) = 0.0D0

      DO J = IC0 + 1, ICEND
        N = J - IC0
        DO I = IR0 + 1, IREND
          M = I - IR0

          ! IF IRG == ICG, Up triangle needed. 
          IF (IRG == ICG .AND. I > J) EXIT
          
          IF ((.NOT.LFORDR) .OR. (.NOT.LICCUT) .OR. IRG.NE.ICG .OR. LSYMD) THEN
            DDRS(M,N) = DSUBRS(EOL,I,J,JBLOCK)

          ELSEIF (M == N) THEN
          ! Here: LFORDR = .TRUE. (ZF); LICCUT = .TRUE. (ICG > ICCUT);
          ! IRG = ICG; LSYMD = .FALSE.; and M = N : 
          ! Zero-first approxiamtion same as GRASP2018.
            DDRS(M,N) = DSUBRS(EOL,I,J,JBLOCK)
          ENDIF
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE CALDDRS