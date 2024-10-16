!=======================================================================
! Print the constructed FICTIOUS CSFs for check.                       *
! Written by Chongyang Chen,      Fudan University,          Jan 2022  *
!=======================================================================
      SUBROUTINE PRINTFICTCSF(ICSYM, IRSYM, ICSF)

      USE orb_C,         ONLY: IQA, NW
      USE stat_C,        ONLY: JQSA, JCUPA
      USE symmatrix_mod, ONLY: MAP

      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER    :: ICSYM, IRSYM, ICSF
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER   :: I

      WRITE(*,*)'ICSYM, IRSYM, MAP(ICSYM), MAP(IRSYM),ICSF=', &
                 ICSYM, IRSYM, MAP(ICSYM), MAP(IRSYM),ICSF
      WRITE(*,'(250(I4))')(I, I=1,NW+1)
      WRITE(*,'(250(I4))')(IQA(I,ICSF),I=1,NW)
      WRITE(*,'(250(I4))')(JQSA(I,1,ICSF),I=1,NW+1)
      WRITE(*,'(250(I4))')(JQSA(I,2,ICSF),I=1,NW+1)
      WRITE(*,'(250(I4))')(JQSA(I,3,ICSF),I=1,NW+1)
      WRITE(*,'(250(I4))')(JCUPA(I,ICSF),I=1,NW+1)

      END SUBROUTINE PRINTFICTCSF

