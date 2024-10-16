!***********************************************************************
!                                                                      *
      SUBROUTINE GENINTRKWRAP_CSFG(MYID, NPROCS, J2MAX)
!
!   Written by     Xinghong He            Last revision: 12 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:42:46   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE CTEILSRK_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE genintrk_csfg_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: MYID
      INTEGER  :: NPROCS
      INTEGER  :: J2MAX
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL  :: LFRST = .TRUE.
      INTEGER  :: N = 0
!-----------------------------------------------
      IF (.NOT. LFRST) THEN 
        INDTEIRK = 0
        VALTEIRK = 0.0d0 
      ENDIF  
      CALL GENINTRK_CSFG (LFRST, MYID, NPROCS, N, J2MAX)

! Gather integrals (and their indeces) from- and send to- all nodes
      CALL GISUMMPI (INDTEIRK, N)
      CALL GDSUMMPI (VALTEIRK, N)

      RETURN
      END SUBROUTINE GENINTRKWRAP_CSFG
