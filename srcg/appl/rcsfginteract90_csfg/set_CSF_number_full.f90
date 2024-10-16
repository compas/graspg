!***********************************************************************
!                                                                      *
      SUBROUTINE SET_CSF_number_full
!                                                                      *
!                                                                      *
!                                                                      *
!   Written by  G. Gaigalas                   NIST, December 2015      *
!                                                                      * 
!   CYC: 2023/07/17:                                                   * 
!   Copied from set_csf_number for full-CSF                            *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!    M O D U L E S
!-----------------------------------------------
      USE rang_Int_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER   :: IOS, I_num, I_Num_BLK
      CHARACTER :: RECORD*15, S_closed*181, S_orbitals*1070
!-----------------------------------------------
!
      NUM_in_BLK_full(1:20) = 0
      I_num = 0
      I_Num_BLK = 1
      DO
         READ (20, '(A)', IOSTAT=IOS) S_orbitals
         IF (S_orbitals(1:2) == ' *') THEN
            NUM_in_BLK_full(I_Num_BLK) = I_num
            I_Num_BLK = I_Num_BLK + 1
            I_num = 0
            READ (20, '(A)') S_orbitals
         ENDIF
         IF (IOS == 0) THEN
            I_num = I_num + 1
!     Read the J_sub and v quantum numbers
            READ (20, '(A)') S_orbitals
!     Read the X, J, and (sign of) P quantum numbers
            READ (20, '(A)') S_orbitals
            CYCLE
         END IF
         EXIT
      END DO
      NUM_in_BLK_full(I_Num_BLK) = I_num
!      WRITE(*,*)"NUM_in_BLK_full =",NUM_in_BLK_full(1:I_Num_BLK)
      I_Num_BLK = I_Num_BLK + 1
!---------------------------------------------------
      I_Num_BLK = 1
      REWIND(20)
      READ (20, '(1A15)') RECORD
!     Closed orbitals
      READ (20, '(A)') S_closed
      READ (20, '(1A15)') RECORD
!     Peel orbitals
      READ (20, '(A)') S_orbitals
      READ (20, '(1A7)') RECORD
      RETURN
      END SUBROUTINE SET_CSF_number_full
