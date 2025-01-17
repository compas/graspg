!***********************************************************************
!                                                                      *
      INTEGER FUNCTION csfg_JCUP (LOC, ICSF)
!                                                                      *
!   JCUP is the 2J+1 value of the LOCth nontrivial intermediate ang-   *
!   ular momentum in CSF  ICSF.                                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  08:16:51   2/21/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def, ONLY: NNNW
      USE IOUNIT_C,      ONLY: ISTDE
      use orb_C,         ONLY: NW, NCF
      use stat_C,        ONLY: JCUPA
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: LOC
      INTEGER             :: ICSF
!-----------------------------------------------
!CYC:
! For symbolic-CSF list, some additional FICRIOUS ones are needed.
      IF (LOC>=1 .AND. LOC<=NNNW-1) THEN
         !IF (ICSF>=1 .AND. ICSF<=NCF) THEN
         IF (ICSF>=1 .AND. ICSF<=NCF+10) THEN
            csfg_jcup = jcupa(loc,icsf)
         ELSE
            WRITE (ISTDE, *) 'csfg_JCUP: Argument ICSF is out of range.'
            WRITE (ISTDE, *) 'ICSF, NCF+10 =', ICSF, NCF+10
            STOP
         ENDIF
      ELSE
         WRITE (ISTDE, *) 'csfg_JCUP: Argument LOC is out of range.'
         STOP
      ENDIF
!
      RETURN
      END FUNCTION csfg_JCUP
