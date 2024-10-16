!***********************************************************************
!                                                                      *
      INTEGER FUNCTION csfg_ISPAR (ICSF)
!                                                                      *
!   ISPAR is the value of P for CSF number ICSF.                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:41   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def, ONLY: NNNW
      USE STAT_C,        ONLY: JCUPA
      USE IOUNIT_C,      ONLY: ISTDE
      USE ORB_C,         ONLY: NCF
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER       :: ICSF
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!CYC:
! For symbolic-CSF list, some additional FICRIOUS ones are needed.
      !IF (ICSF>=1 .AND. ICSF<=NCF) THEN
      IF (ICSF>=1 .AND. ICSF<=NCF+10) THEN
         csfg_ispar = jcupa(NNNW,icsf)
         IF (csfg_ISPAR > 127) csfg_ISPAR = csfg_ISPAR - 256
         csfg_ISPAR = SIGN(1,csfg_ISPAR)
      ELSE
         WRITE (ISTDE, *) 'ISPAR: Argument ICSF is out of range.'
         WRITE (ISTDE, *) 'ICSF =',ICSF,' NCF+10 =',NCF+10
         STOP
      ENDIF
!
      RETURN
      END FUNCTION csfg_ISPAR
