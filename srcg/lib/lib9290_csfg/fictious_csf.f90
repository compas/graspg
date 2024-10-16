!**********************************************************************
!                                                                     *      
      SUBROUTINE FICTIOUS_CSF(INTYPE, ICSFN, ICSFO, IORBN1, &
                              IORBO1,IORBN2,IORBO2)
!                                                                     *
!  This subroutine builds the TWO fictious CSFGs used                 *
!  within ONESCALAR and RKCO_GG.                                      *
!                                                                     *
!  Written by Chongyang Chen, Jan 2022                                *
!********************************************************************** 
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------        
      USE orb_C,           ONLY: IQA
      USE stat_C,          ONLY: JQSA, JCUPA
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER INTYPE, ICSFN, ICSFO, IORBN1, IORBN2,IORBO1, IORBO2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER I, J, M, N, ICSF
      INTEGER THZERO(3)
      DATA THZERO/0,0,1/ 
!-----------------------------------------------
!
      IQA(:,ICSFN) = IQA(:,ICSFO)
      JQSA(:,:,ICSFN) = JQSA(:,:,ICSFO)
      JCUPA(:,ICSFN) = JCUPA(:,ICSFO)

      ! Copy ICSFO to ICSFN
      IF (IORBN1.EQ.IORBO1 .AND. IORBN2.EQ.IORBO2) RETURN

      SELECT CASE (INTYPE)

        CASE(2)
          IQA(IORBN1, ICSFN) = 1 
          JQSA(IORBN1,1:3,ICSFN) = JQSA(IORBO1,1:3,ICSFO)
          JCUPA(IORBN1,ICSFN) = JCUPA(IORBO1,ICSFO)

          IQA(IORBO1, ICSFN) = 0
          JQSA(IORBO1,1:3,ICSFN) = THZERO
          JCUPA(IORBO1,ICSFN) = 0

        !CASE(3)
        !  RETURN

        CASE(4)
          IQA(IORBN1, ICSFN) = 1 
          JQSA(IORBN1,1:3,ICSFN) = JQSA(IORBO1,1:3,ICSFO)
          JCUPA(IORBN1,ICSFN) = JCUPA(IORBO1,ICSFO)

          IF (IORBN2.EQ.0) THEN
            ! One orbital is shifted
            IQA(IORBO1, ICSFN) = 0
            JQSA(IORBO1,1:3,ICSFN) = THZERO
            JCUPA(IORBO1,ICSFN) = 0
          ELSE
          ! Two electron are shifted adjacently 
            IF (IORBN1+1 .NE. IORBN2) THEN
              WRITE(*,*)'Unexpected IORBN1+1 .NE. IORBN2 ...'
              STOP 'Unexpected IORBN1+1 .NE. IORBN2 ...'
            ENDIF
            IQA(IORBN2, ICSFN) = 1
            JQSA(IORBN2,1:3,ICSFN) = JQSA(IORBO2,1:3,ICSFO)
            JCUPA(IORBN2,ICSFN) = JCUPA(IORBO2,ICSFO)
 
            IF (IORBO1.NE.IORBN1 .AND. IORBO1.NE.IORBN2) THEN
              IQA(IORBO1, ICSFN) = 0
              JQSA(IORBO1,1:3,ICSFN) = THZERO
              JCUPA(IORBO1,ICSFN) = 0
            ENDIF
            IF (IORBO2.NE.IORBN1 .AND. IORBO2.NE.IORBN2) THEN
              IQA(IORBO2, ICSFN) = 0
              JQSA(IORBO2,1:3,ICSFN) = THZERO
              JCUPA(IORBO2,ICSFN) = 0
            ENDIF
          ENDIF

        CASE(5) 
          IQA(IORBN1, ICSFN) = 2
          JQSA(IORBN1,1:3,ICSFN) = JQSA(IORBO1,1:3,ICSFO)
          JCUPA(IORBN1,ICSFN) = JCUPA(IORBO1,ICSFO)

          IQA(IORBO1, ICSFN) = 0
          JQSA(IORBO1,1:3,ICSFN) = THZERO
          JCUPA(IORBO1,ICSFN) = 0
          
        CASE DEFAULT
          WRITE(*,*)"Error! INTYPE WRONG ..."
          STOP "Error! INTYPE WRONG ..."
      END SELECT
  
      RETURN
      END SUBROUTINE FICTIOUS_CSF

