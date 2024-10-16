      subroutine removeblank(line1,m,line2,n)
!...Translated by  Gediminas Gaigalas  04/21/2021

      implicit none
      character*300 line1
      character*300 line2
      integer m,n,i

      n = 0
      line2 = ' '
      do i = 1, m
       if (line1(i:i).ne.' ') then
         n = n + 1
         line2(n:n) = line1(i:i)
       endif
      enddo

      return
      end

