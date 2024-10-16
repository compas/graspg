
      SUBROUTINE LOCATE(ILABWORK,NUM,IARR,LFOUND,LOC)
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
      INTEGER, INTENT(IN) :: ILABWORK, NUM
      INTEGER :: IARR(*)
      LOGICAL, INTENT(OUT) :: LFOUND
      INTEGER, INTENT(OUT) :: LOC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: JU, JL, JM

      LFOUND = .FALSE.
      IF ( NUM == 0 ) THEN
        LOC = 0

      ELSEIF ( NUM == 1 ) THEN
        IF (ILABWORK == IARR(1)) THEN
          LFOUND = .TRUE.
          LOC = 1
        ELSEIF (ILABWORK > IARR(1)) THEN
          LOC = 1
        ELSE
          LOC = 0
        ENDIF

      ELSEIF ( NUM > 1 ) THEN
        IF (ILABWORK > IARR(NUM)) THEN
          LOC = NUM
        ELSEIF (ILABWORK < IARR(1)) THEN
          LOC = 0
        ELSE
          JL = 1
          JU = NUM

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
              LOC = JL
            ENDIF
          ENDIF
        ENDIF
      ELSE
        WRITE(*,*)"Unexpected Num < 0 in LOCATE ... "
        STOP "Unexpected Num < 0 in LOCATE ... "
      ENDIF

      RETURN
      END SUBROUTINE LOCATE
