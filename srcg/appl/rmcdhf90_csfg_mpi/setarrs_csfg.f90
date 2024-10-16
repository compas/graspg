!*******************************************************************
!                                                                  *
      SUBROUTINE SetArrs_csfg(JBLOCK)
!                                                                  *
!   Written by  Chongyang Chen, Fudan university,  Jan 2024        *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE csfg_memory_man
      USE hmat_C,     ONLY: EMT, IROW, IENDC, NELMNT
      USE orb_C,      ONLY: NCF 
      USE mpi_C
      USE hblock_C 
      USE csfg_tv_C
      USE symexpand_mod
      USE symmatrix_mod, CUTOFF0=>CUTOFF
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: JBLOCK
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER :: I, NCOLS
!-----------------------------------------------
      NCF = NCFBLKTr(JBLOCK)

      NELMNT = INT8(NZBLK(JBLOCK))
      CALL ALLOC (IROW, NELMNT, 'IROW', 'MATRIX_CSFG')
      IROW(:) = 0
      CALL ALLOC (EMT, NELMNT, 'EMT', 'MATRIX_CSFG')
      EMT(:) = 0.0D0
      NCOLS = NCOLBLK(JBLOCK)
      CALL ALLOC (IENDC, 0, NCF, 'IENDC', 'MATRIX_CSFG')
      IENDC(:) = 0

! NODENCOLS will be set again in setham_pot.f90, and ensure it is equal
! to NCOLS, i.e., NCOLBLK(JBLOCK).
      NODENCOLS = NCOLS

! NODECOLS used in spicmvsym.f, allocated in SetArrs_csfg.f90, 
! deallocated in matrix_csfg.f90, set values in setham_pot.f90, checked
! in checkh_csfg.f90
      ALLOCATE(NODECOLS(NCOLS)) 

      RETURN
      END SUBROUTINE SetArrs_csfg
