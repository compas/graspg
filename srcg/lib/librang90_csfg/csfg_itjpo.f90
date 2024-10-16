!***********************************************************************
!                                                                      *
      INTEGER FUNCTION csfg_ITJPO (ICSF)
!                                                                      *
!   ITJPO is the value of 2J+1 for CSF number ICSF.                    *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!     Last modification by Chongyang Chen   Nov 2023                   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:45   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def, ONLY: NNNW
      USE STAT_C,        ONLY: JCUPA
      USE IOUNIT_C,      ONLY: ISTDE
      USE orb_C,         ONLY: NCF
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER       :: ICSF
!-----------------------------------------------
!CYC:
! For symbolic-CSF list, some additional FICRIOUS ones are needed.
      !IF (ICSF>=1 .AND. ICSF<=NCF) THEN
      IF (ICSF>=1 .AND. ICSF<=NCF+10) THEN
        csfg_itjpo = jcupa(NNNW,icsf)
        IF (csfg_ITJPO > 127) csfg_ITJPO = 256 - csfg_ITJPO
        csfg_ITJPO = IABS (csfg_ITJPO)
      ELSE
         WRITE (ISTDE, *) 'ITJPO: Argument ICSF is out of range.'
         WRITE (ISTDE, *) 'ICSF =',ICSF,' NCF+10 =',NCF+10
         STOP
      ENDIF
!
      RETURN
      END FUNCTION csfg_ITJPO
