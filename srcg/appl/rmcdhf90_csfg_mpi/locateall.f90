      SUBROUTINE LOCATEALL(ILABWORK,NUMB,NUME,IARR,LFOUND,LOC)
!======================================================================
! Written by Chongyang Chen, Fudan University, Shanghai, Dec 2021
!
! IARR is one sorted arrays
!
!======================================================================

      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ILABWORK, NUMB, NUME
      INTEGER :: IARR(*)

      LOGICAL, INTENT(OUT) :: LFOUND
      INTEGER, INTENT(OUT) :: LOC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: JU, JL, JM

      LFOUND = .FALSE.
      IF ( NUMB > NUME ) THEN
        WRITE(*,*)"Unexpected NUMB > NUME in LOCATEALL ..."
        STOP "Unexpected NUMB > NUME in LOCATEALL ..."
        RETURN
      ELSEIF ( NUMB == NUME ) THEN
        IF (ILABWORK == IARR(NUMB)) THEN
          LFOUND = .TRUE.
          LOC = NUMB
        ELSE
          WRITE(*,*)"Unexpected ILABWORK != IARR(NUMB) in LOCATEALL ..."
          STOP "Unexpected ILABWORK != IARR(NUMB) in LOCATEALL ..."
          RETURN
        ENDIF

      ELSE
        IF (ILABWORK > IARR(NUME)) THEN
          WRITE(*,*)"Unexpected ILABWORK > IARR(NUME) in LOCATEALL ..."
          STOP "Unexpected ILABWORK > IARR(NUME) in LOCATEALL ..."
          RETURN
        ELSEIF (ILABWORK < IARR(NUMB)) THEN
          WRITE(*,*)"Unexpected ILABWORK < IARR(NUMB) in LOCATEALL ..."
          STOP "Unexpected ILABWORK < IARR(NUMB) in LOCATEALL ..."
          RETURN
        ELSE
          JL = NUMB
          JU = NUME

  1       IF (JU-JL .GT. 1) THEN
            JM = (JU + JL)/2
            !JM = JL + (ILABWORK - IARR(JL)) / (IARR(JU)-IARR(JL)) * (JU-JL)
            IF (IARR(JM) .GT. ILABWORK) THEN
              JU = JM
            ELSE
              JL = JM
            ENDIF
            GOTO 1

          ELSE
            IF (ILABWORK .EQ. IARR(JU)) THEN
              LFOUND = .TRUE.
              LOC = JU
            ELSEIF (ILABWORK .EQ. IARR(JL)) THEN
              LFOUND = .TRUE.
              LOC = JL
            ELSE
              !WRITE(*,*)"Unexpected ILABWORK .NOT. FOUND in LOCATEALL "
              !WRITE(*,*)'JU, JL, IARR(JU), IARR(JL), ILABWORK =', &
              !          JU, JL, IARR(JU), IARR(JL), ILABWORK 
              !STOP "Unexpected ILABWORK .NOT. FOUND in LOCATEALL "
              LFOUND = .FALSE.
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE LOCATEALL
