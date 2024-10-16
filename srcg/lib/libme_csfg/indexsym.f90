!***********************************************************************
!                                                                      *
      INTEGER function indexsym(char2)
!                                                                      *
!cyc                                                                   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(LEN=2) :: char2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(LEN=2) :: sym(21)
      INTEGER :: i
!-----------------------------------------------
!
! Intitialize symmetry and number

      sym(1)  = 's '
      sym(2)  = 'p-'
      sym(3)  = 'p '
      sym(4)  = 'd-'
      sym(5)  = 'd '
      sym(6)  = 'f-'
      sym(7)  = 'f '
      sym(8)  = 'g-'
      sym(9)  = 'g '
      sym(10) = 'h-'
      sym(11) = 'h '
      sym(12) = 'i-'
      sym(13) = 'i '
      sym(14) = 'k-'
      sym(15) = 'k '
      sym(16) = 'l-'
      sym(17) = 'l '
      sym(18) = 'm-'
      sym(19) = 'm '
      sym(20) = 'n-'
      sym(21) = 'n '
      do i=1,21
         if (char2.eq.sym(i)) then
            indexsym=i
            return
         endif
      enddo

      stop 'Error, Not found the symetery in the list...'
      END function indexsym
