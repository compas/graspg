!***********************************************************************
!                                                                      *
      SUBROUTINE LODCSL_Full(NEXT_BLOCKF) 
!                                                                      *
!   Loads the data from the  .csl  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!                                                                      *
!                                                                      *
!   Written by  G. Gaigalas                        Vilnius, May 2016   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE symzf_mod

!      USE BLK_C,            only: NBLOCK,NCFBLK
!      USE rang_Int_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL,  INTENT(OUT) :: NEXT_BLOCKF 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER            :: IOS, NCF, NCFD
      CHARACTER(LEN=256) :: RECORD 
!-----------------------------------------------
!
!   Initial allocation for arrays with a dimension dependent
!   on the number of CSFs; the initial allocation must be
!   greater than 1
!
      NEXT_BLOCKF = .TRUE.
      NCF = 0 
      DO
         NCF = NCF + 1 
         READ (20, '(A)', IOSTAT=IOS) RECORD 
         IF (IOS == 0) THEN 
            IF (RECORD(1:2) == ' *') THEN 
               RETURN
            ENDIF 
            !C_shell(NCF) = RECORD
            FCONF1(NCF) = RECORD 
!
!   Read the J_sub and v quantum numbers
!
            READ (20, '(A)', IOSTAT=IOS) RECORD 
            !C_quant(NCF) = RECORD
            FCONF2(NCF) = RECORD
!
!   Read the X, J, and (sign of) P quantum numbers
!
            READ (20, '(A)') RECORD
            !C_coupl(NCF) = RECORD
            FCONF3(NCF) = RECORD
         ELSE
            EXIT
         ENDIF
      END DO
!   The last block
!      NBLOCK = NBLOCK + 1 
!      NCFBLK(NBLOCK) = NCF - 1
!      NotFound = NCFBLK(NBLOCK)
!      Found(1:NotFound) = 0
!      NEXT_BLOCKZ = .FALSE.
      RETURN  
!
      END SUBROUTINE LODCSL_Full
