!***********************************************************************
!                                                                      *
      SUBROUTINE SETDBG
!                                                                      *
!   This subroutine sets the arrays that control debug printout from   *
!   the radial and angular modules of the GRASP92 suite.               *
!                                                                      *
!                                                                      *
!   Written by Farid A Parpia               Last update: 21 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE debug_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!
      LDBPA = .FALSE.
      LDBPG = .FALSE.

      LDBPR = .FALSE.

      RETURN
      END SUBROUTINE SETDBG
