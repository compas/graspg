!***********************************************************************
!                                                                      *
      SUBROUTINE REMOVEBLANK(LINEIN, M, LINEOUT, N, ICSF) 
!                                                                      *
!   Remove the blanks within string Linein                             *
!                                                                      *
!                                                                      *
!   Written by  Chongyang Chen      Fudan University,       Jan 2021   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!      USE symzf_mod

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=256), INTENT(IN) :: LINEIN
      INTEGER,  INTENT(IN) :: M, ICSF 
      CHARACTER(LEN=256), INTENT(OUT) :: LINEOUT
      INTEGER,  INTENT(OUT) :: N 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER            :: I
!-----------------------------------------------
!
      n = 0
      lineout = ' '
      do i = 1, m
       if (linein(i:i).ne.' ') then
         n = n + 1
         lineout(n:n) = linein(i:i)
       endif
      enddo
!      if (len_trim(lineout).eq.0) then 
!       write(*,*)"ICSF=", ICSF, " m=", m
!       write(*,'(a)') trim(linein)
!      endif
      RETURN  
!
      END SUBROUTINE REMOVEBLANK
