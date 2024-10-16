!***********************************************************************
!                                                                      *
      SUBROUTINE SET_CSF_numberF
!                                                                      *
!                                                                      *
!                                                                      *
!   Written by  G. Gaigalas                   NIST, December 2015      *
!                                                                      *
!***********************************************************************
!   Modified from SET_CSF_number.f90 for reading the full CSFs set,
!   initialize the block numbers for the working arrys. 
!   Chongyang Chen, Fudan University, Shanghai, China
!***********************************************************************
!-----------------------------------------------
!    M O D U L E S
!-----------------------------------------------
      USE symzf_mod
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
      CHARACTER :: RECORD*1, S_closed*1, S_orbitals*2
!-----------------------------------------------
!
      NUM_in_BLKF(1:20) = 0
      I_num = 0
      I_Num_BLK = 1
      DO
         READ (20, '(A2)', IOSTAT=IOS) S_orbitals
         IF (IOS == 0) THEN
            IF (S_orbitals(1:2) == ' *') THEN
               NUM_in_BLKF(I_Num_BLK) = I_num
               I_Num_BLK = I_Num_BLK + 1
               I_num = 0
               READ (20, '(A2)') S_orbitals
            ENDIF
            I_num = I_num + 1
!     Read the J_sub and v quantum numbers
            READ (20, '(A2)') S_orbitals
!     Read the X, J, and (sign of) P quantum numbers
            READ (20, '(A2)') S_orbitals
            CYCLE 
         END IF
         EXIT
      END DO
      NUM_in_BLKF(I_Num_BLK) = I_num
      NBlockF = I_Num_BLK
!cyc      I_Num_BLK = I_Num_BLK + 1
!---------------------------------------------------
      I_Num_BLK = 1
      REWIND(20)
! Read the first FIVE lines
      READ (20, '(1A1)') RECORD
!     Closed orbitals
      READ (20, '(A1)') S_closed
      READ (20, '(1A1)') RECORD
!     Peel orbitals
      READ (20, '(A2)') S_orbitals
      READ (20, '(1A1)') RECORD
      RETURN  
      END SUBROUTINE SET_CSF_numberF
