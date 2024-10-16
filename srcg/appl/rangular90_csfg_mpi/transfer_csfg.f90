!**********************************************************************
!                                                                     *      
      SUBROUTINE TRANSFER_CSFG(IC, IR)
!                                                                     *
!  This subroutine transfers the non-zero data in LNonZero            *
!  with the corresponding ACCUMULATED row position to the arrays      *
!  IROWSYM, recording the rows of the non-zero element                *
!                                                                     *
!  Written by Chongyang Chen , Fudan university,       Oct 2023       *
!                                                                     *
!********************************************************************** 
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------        
!CYC      USE decide_C,        ONLY: LFORDR
      USE MCP_C,           ONLY: LFORDR
      use symmatrix_mod,   ONLY: CUTOFF, IRTOT, LICCUT, NTYPE,         &
                                 NELCSYM, IROWSYM, LNonZero
      use csfg_decide_C,   ONLY: LSYMD
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER IC, IR
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER I, J, M, N, IR0
!-----------------------------------------------
      M=NTYPE(2,IR)
      N=NTYPE(2,IC)
      !IF (IC == 12) WRITE(*,*)'IC, IR=',IC,IR, ' IRTOT=',IRTOT
! Needing more logical tests !!!!!
      IF ((.NOT.LFORDR) .OR. (.NOT.LICCUT) .OR. IC.NE.IR .OR. LSYMD)THEN
        ! Include all of the non-zero matrixelements
        ! LSYMD: Diagonal-Block interaction
        DO J = 1, N
           IR0 = IRTOT
           DO I = 1, M
              IR0 = IR0 + 1
              IF (LNonZero(I,J)) THEN
                 NELCSYM(J) = NELCSYM(J) + 1
                 IROWSYM(NELCSYM(J),J)=IR0
              ENDIF
           ENDDO
        ENDDO
        RETURN
      ELSE
      ! Discard the CI within the Diagonal-Block
        DO I = 1, M
          J = I
          NELCSYM(J) = NELCSYM(J) + 1
          IROWSYM(NELCSYM(J),J)= IRTOT + I
        ENDDO
        RETURN
      ENDIF

      RETURN
      END SUBROUTINE TRANSFER_CSFG
