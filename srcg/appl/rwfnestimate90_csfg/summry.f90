

!***********************************************************************
!                                                                      *
      SUBROUTINE SUMMRY(NUNIT)
!                                                                      *
!   Prints  a summary of the complete  list of subshell radial wave-   *
!   functions on NUNIT.                                                *
!                                                                      *
!   Call(s) to: [LIB92]:                                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!   Modified by Yanting Li, output <r> for each orbital, 19 Aug 2022 
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE GRID_C
      USE ORB_C
      USE WAVE_C
      USE WHFROM_C, ONLY: SOURCE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!-----add use of rint_I, Yanting
      USE rint_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NUNIT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, LENTH
!-----------------------------------------------
!
!
      WRITE (NUNIT, 300)
!
      DO I = 1, NW
         LENTH = LEN_TRIM(SOURCE(I))
!-------Replace P(2) and Q(2) output with <r> (RINT(I,I,1))
         WRITE (NUNIT, 301) NP(I), NH(I), E(I), PZ(I), GAMA(I), &
            RINT(I,I,1), MF(I), SOURCE(I)(1:LENTH)
!         WRITE (NUNIT,302) SOURCE(I)(1:LENTH)
      END DO
!
      RETURN
!
  300 FORMAT('Shell',6X,'e',11X,'p0',8X,'gamma',8X,'<r>',6X,'MTP',&
         '  SRC'/)
  301 FORMAT(1X,I2,A2,4D12.4,I5,2X,A3)
      RETURN
!  302 FORMAT ('      Source: ',A)
!
      END SUBROUTINE SUMMRY
