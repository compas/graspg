!     last edited Januar 2, 1997
! CYC modifications for excitations up to n = 25 (2020).
      subroutine blanda(org, varmax, lock, minj, maxj, skal, nmax, low, posn, &
         posl, lim, dubbel, first)
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
!...Switches:
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use slug_I
      use gen_I
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: varmax
      integer  :: minj
      integer  :: maxj
      integer  :: skal
      integer , intent(in) :: nmax
      logical  :: first
      integer  :: org(25,0:10)
      integer  :: low(25,0:10)
      integer  :: posn(220) ! Subshell numbers up to n=25
      integer  :: posl(220) ! Subshell numbers
      integer , intent(in) :: lim(25)
      logical  :: lock(25,0:10)
      logical  :: dubbel(25,0:10)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: fil_1 = 7
      integer, parameter :: fil_2 = 8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(25,0:10) :: antel, start
      integer :: cf
      integer , dimension(25,0:10,0:1) :: ansats
      integer , dimension(25,0:10) :: varupp, varned
      integer :: an10, an20, an21, an30, an31, an32, an40, an41, an42, an43, k&
         , an50, an51, an52, an53, an54, an60, an61, an62, an63, an64, an65, &
         an70, an71, an72, an73, an74, an75, an76
      integer , dimension(25,0:10) :: stopp
      integer :: an80, an81, an82, an83, an84, an85, an86, an87, an90, an91, &
         an92, an93, an94, an95, an96, an97, an98, ana0, ana1, ana2, ana3, ana4&
         , ana5, ana6, ana7, ana8, ana9, plus21, plus31, plus32, plus41, plus42&
         , plus43, plus51, plus52, plus53, plus54, plus61, plus62, plus63, &
         plus64, plus65, plus71, plus72, plus73, plus74, plus75, plus76, plus81&
         , plus82, plus83, plus84, plus85, plus86, plus87, plus91, plus92, &
         plus93, plus94, plus95, plus96, plus97, plus98, plusa1, plusa2, plusa3&
         , plusa4, plusa5, plusa6, plusa7, plusa8, plusa9, par0, par, ress, &
         resl, i, j, antal
      integer , dimension(25,0:10) :: steg
      integer :: dum, ras1, ras3, elar, rasett, rastre
      integer , dimension(25,0:10) :: ras
      integer :: plusba, plusca, plusda, plusea, plusfa, plusb1, plusb2, plusb3&
         , plusb4, plusb5, plusb6, plusb7, plusb8, plusb9, plusc1, plusc2, &
         plusc3, plusc4, plusc5, plusc6, plusc7, plusc8, plusc9, plusd1, plusd2&
         , plusd3, plusd4, plusd5, plusd6, plusd7, plusd8, plusd9, pluse1, &
         pluse2, pluse3, pluse4, pluse5, pluse6, pluse7, pluse8, pluse9, plusf1&
         , plusf2, plusf3, plusf4, plusf5, plusf6, plusf7, plusf8, plusf9, anba&
         , anca, anda, anea, anfa, anb0, anb1, anb2, anb3, anb4, anb5, anb6, &
         anb7, anb8, anb9, anc0, anc1, anc2, anc3, anc4, anc5, anc6, anc7, anc8&
         , anc9, and0, and1, and2, and3, and4, and5, and6, and7, and8, and9, &
         ane0, ane1, ane2, ane3, ane4, ane5, ane6, ane7, ane8, ane9, anf0, anf1&
         , anf2, anf3, anf4, anf5, anf6, anf7, anf8, anf9
!CYC-2020
      integer :: ng16, nh17, ni18, nj19, nk20, nl21, nm22, nn23, no24, np25
! n = 16
      integer :: ang0, ang1, ang2, ang3, ang4, ang5, ang6, ang7, ang8, ang9, anga
      integer :: plusg1, plusg2, plusg3, plusg4, plusg5, plusg6, plusg7, plusg8, plusg9, plusga

! n = 17
      integer :: anh0, anh1, anh2, anh3, anh4, anh5, anh6, anh7, anh8, anh9, anha
      integer :: plush1, plush2, plush3, plush4, plush5, plush6, plush7, plush8, plush9, plusha

! n = 18
      integer :: ani0, ani1, ani2, ani3, ani4, ani5, ani6, ani7, ani8, ani9, ania
      integer :: plusi1, plusi2, plusi3, plusi4, plusi5, plusi6, plusi7, plusi8, plusi9, plusia

! n = 19
      integer :: anj0, anj1, anj2, anj3, anj4, anj5, anj6, anj7, anj8, anj9, anja
      integer :: plusj1, plusj2, plusj3, plusj4, plusj5, plusj6, plusj7, plusj8, plusj9, plusja

! n = 20
      integer :: ank0, ank1, ank2, ank3, ank4, ank5, ank6, ank7, ank8, ank9, anka
      integer :: plusk1, plusk2, plusk3, plusk4, plusk5, plusk6, plusk7, plusk8, plusk9, pluska

! n = 21
      integer :: anl0, anl1, anl2, anl3, anl4, anl5, anl6, anl7, anl8, anl9, anla
      integer :: plusl1, plusl2, plusl3, plusl4, plusl5, plusl6, plusl7, plusl8, plusl9, plusla

! n = 22
      integer :: anm0, anm1, anm2, anm3, anm4, anm5, anm6, anm7, anm8, anm9, anma
      integer :: plusm1, plusm2, plusm3, plusm4, plusm5, plusm6, plusm7, plusm8, plusm9, plusma

! n = 23
      integer :: ann0, ann1, ann2, ann3, ann4, ann5, ann6, ann7, ann8, ann9, anna
      integer :: plusn1, plusn2, plusn3, plusn4, plusn5, plusn6, plusn7, plusn8, plusn9, plusna

! n = 24
      integer :: ano0, ano1, ano2, ano3, ano4, ano5, ano6, ano7, ano8, ano9, anoa
      integer :: pluso1, pluso2, pluso3, pluso4, pluso5, pluso6, pluso7, pluso8, pluso9, plusoa

! n = 25
      integer :: anp0, anp1, anp2, anp3, anp4, anp5, anp6, anp7, anp8, anp9, anpa
      integer :: plusp1, plusp2, plusp3, plusp4, plusp5, plusp6, plusp7, plusp8, plusp9, pluspa
!CYC-2020 end

      logical :: finns, napp
!-----------------------------------------------

      cf = 0
      antal = 0
      par0 = 0
      finns = .FALSE.
      do i = 1, nmax
         do j = 0, min(10,i - 1)
            if (dubbel(i,j)) then
               steg(i,j) = -2
            else
               steg(i,j) = -1
            endif
            antal = antal + org(i,j)
            par0 = mod(par0 + j*org(i,j),2)
         end do
      end do
      if (nmax < 25) then
         do i = nmax + 1, 25
            steg(i,:min(10,i-1)) = -1
         end do
      endif
!     1s
      call slug (1, 0, varmax, varupp, varned, ansats, org, lock(1,0), dubbel, &
         low, start(1,0), stopp(1,0))
      do an10 = start(1,0), stopp(1,0), steg(1,0)
         antel(1,0) = an10
         if (antel(1,0)>antal .or. antel(1,0)<lim(1)) cycle
         ansats(1,0,0) = an10
!     2s
         call slug (2, 0, varmax, varupp, varned, ansats, org, lock(2,0), &
            dubbel, low, start(2,0), stopp(2,0))
         do an20 = start(2,0), stopp(2,0), steg(2,0)
            antel(2,0) = an20 + antel(1,0)
            if (antel(2,0) > antal) cycle
            ansats(2,0,0) = an20
!     2p
            call slug (2, 1, varmax, varupp, varned, ansats, org, lock(2,1), &
               dubbel, low, start(2,1), stopp(2,1))
            do an21 = start(2,1), stopp(2,1), steg(2,1)
               antel(2,1) = an21 + antel(2,0)
               if (antel(2,1)>antal .or. antel(2,1)<lim(2)) cycle
               do plus21 = min(an21,4), max(an21 - 2,0), -1
                  ansats(2,1,1) = plus21
                  ansats(2,1,0) = an21 - plus21
!     3s
                  call slug (3, 0, varmax, varupp, varned, ansats, org, lock(3,&
                     0), dubbel, low, start(3,0), stopp(3,0))
                  do an30 = start(3,0), stopp(3,0), steg(3,0)
                     antel(3,0) = an30 + antel(2,1)
                     if (antel(3,0) > antal) cycle
                     ansats(3,0,0) = an30
!     3p
                     call slug (3, 1, varmax, varupp, varned, ansats, org, lock&
                        (3,1), dubbel, low, start(3,1), stopp(3,1))
                     do an31 = start(3,1), stopp(3,1), steg(3,1)
                        antel(3,1) = an31 + antel(3,0)
                        if (antel(3,1) > antal) cycle
                        do plus31 = min(an31,4), max(an31 - 2,0), -1
                           ansats(3,1,1) = plus31
                           ansats(3,1,0) = an31 - plus31
!     3d
                           call slug (3, 2, varmax, varupp, varned, ansats, org&
                              , lock(3,2), dubbel, low, start(3,2), stopp(3,2))
                           do an32 = start(3,2), stopp(3,2), steg(3,2)
                              antel(3,2) = an32 + antel(3,1)
                              if (antel(3,2)>antal .or. antel(3,2)<lim(3)) &
                                 cycle
                              do plus32 = min(an32,6), max(an32 - 4,0), -1
                                 ansats(3,2,1) = plus32
                                 ansats(3,2,0) = an32 - plus32
!     4s
                                 call slug (4, 0, varmax, varupp, varned, &
                                    ansats, org, lock(4,0), dubbel, low, start(&
                                    4,0), stopp(4,0))
                                 do an40 = start(4,0), stopp(4,0), steg(4,0)
                                    antel(4,0) = an40 + antel(3,2)
                                    if (antel(4,0) > antal) cycle
                                    ansats(4,0,0) = an40
!     4p
                                    call slug (4, 1, varmax, varupp, varned, &
                                       ansats, org, lock(4,1), dubbel, low, &
                                       start(4,1), stopp(4,1))
                                    do an41 = start(4,1), stopp(4,1), steg(4,1)
                                    antel(4,1) = an41 + antel(4,0)
                                    if (antel(4,1) > antal) cycle
                                    do plus41 = min(an41,4), max(an41 - 2,0), &
                                       -1
                                    ansats(4,1,1) = plus41
                                    ansats(4,1,0) = an41 - plus41
!     4d
                                    call slug (4, 2, varmax, varupp, varned, &
                                       ansats, org, lock(4,2), dubbel, low, &
                                       start(4,2), stopp(4,2))
                                    do an42 = start(4,2), stopp(4,2), steg(4,2)
                                    antel(4,2) = an42 + antel(4,1)
                                    if (antel(4,2) > antal) cycle
                                    do plus42 = min(an42,6), max(an42 - 4,0), &
                                       -1
                                    ansats(4,2,1) = plus42
                                    ansats(4,2,0) = an42 - plus42
!     4f
                                    call slug (4, 3, varmax, varupp, varned, &
                                       ansats, org, lock(4,3), dubbel, low, &
                                       start(4,3), stopp(4,3))
                                    do an43 = start(4,3), stopp(4,3), steg(4,3)
                                    antel(4,3) = an43 + antel(4,2)
                                    if (antel(4,3)>antal .or. antel(4,3)<lim(4)&
                                       ) cycle
                                    do plus43 = min(an43,8), max(an43 - 6,0), &
                                       -1
                                    ansats(4,3,1) = plus43
                                    ansats(4,3,0) = an43 - plus43
!     5s
                                    call slug (5, 0, varmax, varupp, varned, &
                                       ansats, org, lock(5,0), dubbel, low, &
                                       start(5,0), stopp(5,0))
                                    do an50 = start(5,0), stopp(5,0), steg(5,0)
                                    antel(5,0) = an50 + antel(4,3)
                                    if (antel(5,0) > antal) cycle
                                    ansats(5,0,0) = an50
!     5p
                                    call slug (5, 1, varmax, varupp, varned, &
                                       ansats, org, lock(5,1), dubbel, low, &
                                       start(5,1), stopp(5,1))
                                    do an51 = start(5,1), stopp(5,1), steg(5,1)
                                    antel(5,1) = an51 + antel(5,0)
                                    if (antel(5,1) > antal) cycle
                                    do plus51 = min(an51,4), max(an51 - 2,0), &
                                       -1
                                    ansats(5,1,1) = plus51
                                    ansats(5,1,0) = an51 - plus51
!     5d
                                    call slug (5, 2, varmax, varupp, varned, &
                                       ansats, org, lock(5,2), dubbel, low, &
                                       start(5,2), stopp(5,2))
                                    do an52 = start(5,2), stopp(5,2), steg(5,2)
                                    antel(5,2) = an52 + antel(5,1)
                                    if (antel(5,2) > antal) cycle
                                    do plus52 = min(an52,6), max(an52 - 4,0), &
                                       -1
                                    ansats(5,2,1) = plus52
                                    ansats(5,2,0) = an52 - plus52

!     5f
                                    call slug (5, 3, varmax, varupp, varned, &
                                       ansats, org, lock(5,3), dubbel, low, &
                                       start(5,3), stopp(5,3))
                                    do an53 = start(5,3), stopp(5,3), steg(5,3)
                                    antel(5,3) = an53 + antel(5,2)
                                    if (antel(5,3) > antal) cycle
                                    do plus53 = min(an53,8), max(an53 - 6,0), &
                                       -1
                                    ansats(5,3,1) = plus53
                                    ansats(5,3,0) = an53 - plus53
!     5g
                                    call slug (5, 4, varmax, varupp, varned, &
                                       ansats, org, lock(5,4), dubbel, low, &
                                       start(5,4), stopp(5,4))
                                    do an54 = start(5,4), stopp(5,4), steg(5,4)
                                    antel(5,4) = an54 + antel(5,3)
                                    if (antel(5,4)>antal .or. antel(5,4)<lim(5)&
                                       ) cycle
                                    do plus54 = min(an54,10), max(an54 - 8,0), &
                                       -1
                                    ansats(5,4,1) = plus54
                                    ansats(5,4,0) = an54 - plus54
!     6s
                                    call slug (6, 0, varmax, varupp, varned, &
                                       ansats, org, lock(6,0), dubbel, low, &
                                       start(6,0), stopp(6,0))
                                    do an60 = start(6,0), stopp(6,0), steg(6,0)
                                    antel(6,0) = an60 + antel(5,4)
                                    if (antel(6,0)>antal .or. ansats(5,4,1)>2) &
                                       cycle
                                    ansats(6,0,0) = an60
!     6p
                                    call slug (6, 1, varmax, varupp, varned, &
                                       ansats, org, lock(6,1), dubbel, low, &
                                       start(6,1), stopp(6,1))
                                    do an61 = start(6,1), stopp(6,1), steg(6,1)
                                    antel(6,1) = an61 + antel(6,0)
                                    if (antel(6,1) > antal) cycle
                                    do plus61 = min(an61,4), max(an61 - 2,0), &
                                       -1
                                    ansats(6,1,1) = plus61
                                    ansats(6,1,0) = an61 - plus61
!     6d
                                    call slug (6, 2, varmax, varupp, varned, &
                                       ansats, org, lock(6,2), dubbel, low, &
                                       start(6,2), stopp(6,2))
                                    do an62 = start(6,2), stopp(6,2), steg(6,2)
                                    antel(6,2) = an62 + antel(6,1)
                                    if (antel(6,2) > antal) cycle
                                    do plus62 = min(an62,6), max(an62 - 4,0), &
                                       -1
                                    ansats(6,2,1) = plus62
                                    ansats(6,2,0) = an62 - plus62
!     6f
                                    call slug (6, 3, varmax, varupp, varned, &
                                       ansats, org, lock(6,3), dubbel, low, &
                                       start(6,3), stopp(6,3))
                                    do an63 = start(6,3), stopp(6,3), steg(6,3)
                                    antel(6,3) = an63 + antel(6,2)
                                    if (antel(6,3) > antal) cycle
                                    do plus63 = min(an63,8), max(an63 - 6,0), &
                                       -1
                                    ansats(6,3,1) = plus63
                                    ansats(6,3,0) = an63 - plus63
!     6g
                                    call slug (6, 4, varmax, varupp, varned, &
                                       ansats, org, lock(6,4), dubbel, low, &
                                       start(6,4), stopp(6,4))
                                    do an64 = start(6,4), stopp(6,4), steg(6,4)
                                    antel(6,4) = an64 + antel(6,3)
                                    if (antel(6,4) > antal) cycle
                                    do plus64 = min(an64,10), max(an64 - 8,0), &
                                       -1
                                    ansats(6,4,1) = plus64
                                    ansats(6,4,0) = an64 - plus64
!     6h
                                    call slug (6, 5, varmax, varupp, varned, &
                                       ansats, org, lock(6,5), dubbel, low, &
                                       start(6,5), stopp(6,5))
                                    do an65 = start(6,5), stopp(6,5), steg(6,5)
                                    antel(6,5) = an65 + antel(6,4)
                                    if (.not.(antel(6,5)<=antal .and. ansats(6,&
                                       4,1)<=2 .and. antel(6,5)>=lim(6))) &
                                       cycle
                                    do plus65 = min(an65,12), max(an65 - 10,0)&
                                       , -1
                                    ansats(6,5,1) = plus65
                                    ansats(6,5,0) = an65 - plus65
!     7s
                                    call slug (7, 0, varmax, varupp, varned, &
                                       ansats, org, lock(7,0), dubbel, low, &
                                       start(7,0), stopp(7,0))
                                    do an70 = start(7,0), stopp(7,0), steg(7,0)
                                    antel(7,0) = an70 + antel(6,5)
                                    if (.not.(antel(7,0)<=antal .and. ansats(6,&
                                       5,1)<=2 .and. ansats(6,5,0)<=2)) cycle
                                    ansats(7,0,0) = an70
!     7p
                                    call slug (7, 1, varmax, varupp, varned, &
                                       ansats, org, lock(7,1), dubbel, low, &
                                       start(7,1), stopp(7,1))
                                    do an71 = start(7,1), stopp(7,1), steg(7,1)
                                    antel(7,1) = an71 + antel(7,0)
                                    if (antel(7,1) > antal) cycle
                                    do plus71 = min(an71,4), max(an71 - 2,0), &
                                       -1
                                    ansats(7,1,1) = plus71
                                    ansats(7,1,0) = an71 - plus71
!     7d
                                    call slug (7, 2, varmax, varupp, varned, &
                                       ansats, org, lock(7,2), dubbel, low, &
                                       start(7,2), stopp(7,2))
                                    do an72 = start(7,2), stopp(7,2), steg(7,2)
                                    antel(7,2) = an72 + antel(7,1)
                                    if (antel(7,2) > antal) cycle
                                    do plus72 = min(an72,6), max(an72 - 4,0), &
                                       -1
                                    ansats(7,2,1) = plus72
                                    ansats(7,2,0) = an72 - plus72
!     7f
                                    call slug (7, 3, varmax, varupp, varned, &
                                       ansats, org, lock(7,3), dubbel, low, &
                                       start(7,3), stopp(7,3))
                                    do an73 = start(7,3), stopp(7,3), steg(7,3)
                                    antel(7,3) = an73 + antel(7,2)
                                    if (antel(7,3) > antal) cycle
                                    do plus73 = min(an73,8), max(an73 - 6,0), &
                                       -1
                                    ansats(7,3,1) = plus73
                                    ansats(7,3,0) = an73 - plus73
!     7g
                                    call slug (7, 4, varmax, varupp, varned, &
                                       ansats, org, lock(7,4), dubbel, low, &
                                       start(7,4), stopp(7,4))
                                    do an74 = start(7,4), stopp(7,4), steg(7,4)
                                    antel(7,4) = an74 + antel(7,3)
                                    if (antel(7,4) > antal) cycle
                                    do plus74 = min(an74,10), max(an74 - 8,0), &
                                       -1
                                    ansats(7,4,1) = plus74
                                    ansats(7,4,0) = an74 - plus74
!     7h
                                    call slug (7, 5, varmax, varupp, varned, &
                                       ansats, org, lock(7,5), dubbel, low, &
                                       start(7,5), stopp(7,5))
                                    do an75 = start(7,5), stopp(7,5), steg(7,5)
                                    antel(7,5) = an75 + antel(7,4)
                                    if (antel(7,5)>antal .or. ansats(7,4,1)>2) &
                                       cycle
                                    do plus75 = min(an75,12), max(an75 - 10,0)&
                                       , -1
                                    ansats(7,5,1) = plus75
                                    ansats(7,5,0) = an75 - plus75
!     7i
                                    call slug (7, 6, varmax, varupp, varned, &
                                       ansats, org, lock(7,6), dubbel, low, &
                                       start(7,6), stopp(7,6))
                                    do an76 = start(7,6), stopp(7,6), steg(7,6)
                                    antel(7,6) = an76 + antel(7,5)
                                    if (.not.(antel(7,6)<=antal .and. ansats(7,&
                                       5,1)<=2 .and. ansats(7,5,0)<=2 .and. &
                                       antel(7,6)>=lim(7))) cycle
                                    do plus76 = min(an76,14), max(an76 - 12,0)&
                                       , -1
                                    ansats(7,6,1) = plus76
                                    ansats(7,6,0) = an76 - plus76
!     8s
                                    call slug (8, 0, varmax, varupp, varned, &
                                       ansats, org, lock(8,0), dubbel, low, &
                                       start(8,0), stopp(8,0))
                                    do an80 = start(8,0), stopp(8,0), steg(8,0)
                                    antel(8,0) = an80 + antel(7,6)
                                    if (.not.(antel(8,0)<=antal .and. ansats(7,&
                                       6,1)<=2 .and. ansats(7,6,0)<=2)) cycle
                                    ansats(8,0,0) = an80
!     8p
                                    call slug (8, 1, varmax, varupp, varned, &
                                       ansats, org, lock(8,1), dubbel, low, &
                                       start(8,1), stopp(8,1))
                                    do an81 = start(8,1), stopp(8,1), steg(8,1)
                                    antel(8,1) = an81 + antel(8,0)
                                    if (antel(8,1) > antal) cycle
                                    do plus81 = min(an81,4), max(an81 - 2,0), &
                                       -1
                                    ansats(8,1,1) = plus81
                                    ansats(8,1,0) = an81 - plus81
!     8d
                                    call slug (8, 2, varmax, varupp, varned, &
                                       ansats, org, lock(8,2), dubbel, low, &
                                       start(8,2), stopp(8,2))
                                    do an82 = start(8,2), stopp(8,2), steg(8,2)
                                    antel(8,2) = an82 + antel(8,1)
                                    if (antel(8,2) > antal) cycle
                                    do plus82 = min(an82,6), max(an82 - 4,0), &
                                       -1
                                    ansats(8,2,1) = plus82
                                    ansats(8,2,0) = an82 - plus82
!     8f
                                    call slug (8, 3, varmax, varupp, varned, &
                                       ansats, org, lock(8,3), dubbel, low, &
                                       start(8,3), stopp(8,3))
                                    do an83 = start(8,3), stopp(8,3), steg(8,3)
                                    antel(8,3) = an83 + antel(8,2)
                                    if (antel(8,3) > antal) cycle
                                    do plus83 = min(an83,8), max(an83 - 6,0), &
                                       -1
                                    ansats(8,3,1) = plus83
                                    ansats(8,3,0) = an83 - plus83
!     8g
                                    call slug (8, 4, varmax, varupp, varned, &
                                       ansats, org, lock(8,4), dubbel, low, &
                                       start(8,4), stopp(8,4))

                                    do an84 = start(8,4), stopp(8,4), steg(8,4)
                                    antel(8,4) = an84 + antel(8,3)
                                    if (antel(8,4) > antal) cycle
                                    do plus84 = min(an84,10), max(an84 - 8,0), &
                                       -1
                                    ansats(8,4,1) = plus84
                                    ansats(8,4,0) = an84 - plus84
!     8h
                                    call slug (8, 5, varmax, varupp, varned, &
                                       ansats, org, lock(8,5), dubbel, low, &
                                       start(8,5), stopp(8,5))
                                    do an85 = start(8,5), stopp(8,5), steg(8,5)
                                    antel(8,5) = an85 + antel(8,4)
                                    if (antel(8,5)>antal .or. ansats(8,4,1)>2) &
                                       cycle
                                    do plus85 = min(an85,12), max(an85 - 10,0)&
                                       , -1
                                    ansats(8,5,1) = plus85
                                    ansats(8,5,0) = an85 - plus85
!     8i
                                    call slug (8, 6, varmax, varupp, varned, &
                                       ansats, org, lock(8,6), dubbel, low, &
                                       start(8,6), stopp(8,6))
                                    do an86 = start(8,6), stopp(8,6), steg(8,6)
                                    antel(8,6) = an86 + antel(8,5)
                                    if (.not.(antel(8,6)<=antal .and. ansats(8,&
                                       5,1)<=2 .and. ansats(8,5,0)<=2)) cycle
                                    do plus86 = min(an86,14), max(an86 - 12,0)&
                                       , -1
                                    ansats(8,6,1) = plus86
                                    ansats(8,6,0) = an86 - plus86
!     8k
                                    call slug (8, 7, varmax, varupp, varned, &
                                       ansats, org, lock(8,7), dubbel, low, &
                                       start(8,7), stopp(8,7))
                                    do an87 = start(8,7), stopp(8,7), steg(8,7)
                                    antel(8,7) = an87 + antel(8,6)
                                    if (.not.(antel(8,7)<=antal .and. ansats(8,&
                                       6,1)<=2 .and. ansats(8,6,0)<=2 .and. &
                                       antel(8,7)>=lim(8))) cycle
                                    do plus87 = min(an87,16), max(an87 - 14,0)&
                                       , -1
                                    ansats(8,7,1) = plus87
                                    ansats(8,7,0) = an87 - plus87
!     9s
                                    call slug (9, 0, varmax, varupp, varned, &
                                       ansats, org, lock(9,0), dubbel, low, &
                                       start(9,0), stopp(9,0))
                                    do an90 = start(9,0), stopp(9,0), steg(9,0)
                                    antel(9,0) = an90 + antel(8,7)
                                    if (.not.(antel(9,0)<=antal .and. ansats(8,&
                                       7,1)<=2 .and. ansats(8,7,0)<=2)) cycle
                                    ansats(9,0,0) = an90
!     9p
                                    call slug (9, 1, varmax, varupp, varned, &
                                       ansats, org, lock(9,1), dubbel, low, &
                                       start(9,1), stopp(9,1))
                                    do an91 = start(9,1), stopp(9,1), steg(9,1)
                                    antel(9,1) = an91 + antel(9,0)
                                    if (antel(9,1) > antal) cycle
                                    do plus91 = min(an91,4), max(an91 - 2,0), &
                                       -1
                                    ansats(9,1,1) = plus91
                                    ansats(9,1,0) = an91 - plus91
!     9d
                                    call slug (9, 2, varmax, varupp, varned, &
                                       ansats, org, lock(9,2), dubbel, low, &
                                       start(9,2), stopp(9,2))
                                    do an92 = start(9,2), stopp(9,2), steg(9,2)
                                    antel(9,2) = an92 + antel(9,1)
                                    if (antel(9,2) > antal) cycle
                                    do plus92 = min(an92,6), max(an92 - 4,0), &
                                       -1
                                    ansats(9,2,1) = plus92
                                    ansats(9,2,0) = an92 - plus92
!     9f
                                    call slug (9, 3, varmax, varupp, varned, &
                                       ansats, org, lock(9,3), dubbel, low, &
                                       start(9,3), stopp(9,3))
                                    do an93 = start(9,3), stopp(9,3), steg(9,3)
                                    antel(9,3) = an93 + antel(9,2)
                                    if (antel(9,3) > antal) cycle
                                    do plus93 = min(an93,8), max(an93 - 6,0), &
                                       -1
                                    ansats(9,3,1) = plus93
                                    ansats(9,3,0) = an93 - plus93
!     9g
                                    call slug (9, 4, varmax, varupp, varned, &
                                       ansats, org, lock(9,4), dubbel, low, &
                                       start(9,4), stopp(9,4))
                                    do an94 = start(9,4), stopp(9,4), steg(9,4)
                                    antel(9,4) = an94 + antel(9,3)
                                    if (antel(9,4) > antal) cycle
                                    do plus94 = min(an94,10), max(an94 - 8,0), &
                                       -1
                                    ansats(9,4,1) = plus94
                                    ansats(9,4,0) = an94 - plus94
!     9h
                                    call slug (9, 5, varmax, varupp, varned, &
                                       ansats, org, lock(9,5), dubbel, low, &
                                       start(9,5), stopp(9,5))
                                    do an95 = start(9,5), stopp(9,5), steg(9,5)
                                    antel(9,5) = an95 + antel(9,4)
                                    if (antel(9,5)>antal .or. ansats(9,4,1)>2) &
                                       cycle
                                    do plus95 = min(an95,12), max(an95 - 10,0)&
                                       , -1
                                    ansats(9,5,1) = plus95
                                    ansats(9,5,0) = an95 - plus95
!     9i
                                    call slug (9, 6, varmax, varupp, varned, &
                                       ansats, org, lock(9,6), dubbel, low, &
                                       start(9,6), stopp(9,6))
                                    do an96 = start(9,6), stopp(9,6), steg(9,6)
                                    antel(9,6) = an96 + antel(9,5)
                                    if (.not.(antel(9,6)<=antal .and. ansats(9,&
                                       5,1)<=2 .and. ansats(9,5,0)<=2)) cycle
                                    do plus96 = min(an96,14), max(an96 - 12,0)&
                                       , -1
                                    ansats(9,6,1) = plus96
                                    ansats(9,6,0) = an96 - plus96
!     9k
                                    call slug (9, 7, varmax, varupp, varned, &
                                       ansats, org, lock(9,7), dubbel, low, &
                                       start(9,7), stopp(9,7))
                                    do an97 = start(9,7), stopp(9,7), steg(9,7)
                                    antel(9,7) = an97 + antel(9,6)
                                    if (.not.(antel(9,7)<=antal .and. ansats(9,&
                                       6,1)<=2 .and. ansats(9,6,0)<=2)) cycle
                                    do plus97 = min(an97,16), max(an97 - 14,0)&
                                       , -1
                                    ansats(9,7,1) = plus97
                                    ansats(9,7,0) = an97 - plus97
!     9l
                                    call slug (9, 8, varmax, varupp, varned, &
                                       ansats, org, lock(9,8), dubbel, low, &
                                       start(9,8), stopp(9,8))
                                    do an98 = start(9,8), stopp(9,8), steg(9,8)
                                    antel(9,8) = an98 + antel(9,7)
                                    if (.not.(antel(9,8)<=antal .and. ansats(9,&
                                       7,1)<=2 .and. ansats(9,7,0)<=2 .and. &
                                       antel(9,8)>=lim(9))) cycle
                                    do plus98 = min(an98,18), max(an98 - 16,0)&
                                       , -1
                                    ansats(9,8,1) = plus98
                                    ansats(9,8,0) = an98 - plus98
!     10s
                                    call slug (10, 0, varmax, varupp, varned, &
                                       ansats, org, lock(10,0), dubbel, low, &
                                       start(10,0), stopp(10,0))
                                    do ana0 = start(10,0), stopp(10,0), steg(10&
                                       ,0)
                                    antel(10,0) = ana0 + antel(9,8)
                                    if (.not.(antel(10,0)<=antal .and. ansats(9&
                                       ,8,1)<=2 .and. ansats(9,8,0)<=2)) cycle
                                    ansats(10,0,0) = ana0
!     10p
                                    call slug (10, 1, varmax, varupp, varned, &
                                       ansats, org, lock(10,1), dubbel, low, &
                                       start(10,1), stopp(10,1))
                                    do ana1 = start(10,1), stopp(10,1), steg(10&
                                       ,1)
                                    antel(10,1) = ana1 + antel(10,0)
                                    if (antel(10,1) > antal) cycle
                                    do plusa1 = min(ana1,4), max(ana1 - 2,0), &
                                       -1
                                    ansats(10,1,1) = plusa1
                                    ansats(10,1,0) = ana1 - plusa1
!     10d
                                    call slug (10, 2, varmax, varupp, varned, &
                                       ansats, org, lock(10,2), dubbel, low, &
                                       start(10,2), stopp(10,2))
                                    do ana2 = start(10,2), stopp(10,2), steg(10&
                                       ,2)
                                    antel(10,2) = ana2 + antel(10,1)
                                    if (antel(10,2) > antal) cycle
                                    do plusa2 = min(ana2,6), max(ana2 - 4,0), &
                                       -1
                                    ansats(10,2,1) = plusa2
                                    ansats(10,2,0) = ana2 - plusa2
!     10f
                                    call slug (10, 3, varmax, varupp, varned, &
                                       ansats, org, lock(10,3), dubbel, low, &
                                       start(10,3), stopp(10,3))
                                    do ana3 = start(10,3), stopp(10,3), steg(10&
                                       ,3)
                                    antel(10,3) = ana3 + antel(10,2)
                                    if (antel(10,3) > antal) cycle
                                    do plusa3 = min(ana3,8), max(ana3 - 6,0), &
                                       -1
                                    ansats(10,3,1) = plusa3
                                    ansats(10,3,0) = ana3 - plusa3
!     10g
                                    call slug (10, 4, varmax, varupp, varned, &
                                       ansats, org, lock(10,4), dubbel, low, &
                                       start(10,4), stopp(10,4))
                                    do ana4 = start(10,4), stopp(10,4), steg(10&
                                       ,4)
                                    antel(10,4) = ana4 + antel(10,3)
                                    if (antel(10,4) > antal) cycle
                                    do plusa4 = min(ana4,10), max(ana4 - 8,0), &
                                       -1
                                    ansats(10,4,1) = plusa4
                                    ansats(10,4,0) = ana4 - plusa4
!     10h
                                    call slug (10, 5, varmax, varupp, varned, &
                                       ansats, org, lock(10,5), dubbel, low, &
                                       start(10,5), stopp(10,5))
                                    do ana5 = start(10,5), stopp(10,5), steg(10&
                                       ,5)
                                    antel(10,5) = ana5 + antel(10,4)
                                    if (antel(10,5)>antal .or. ansats(10,4,1)>2&
                                       ) cycle
                                    do plusa5 = min(ana5,12), max(ana5 - 10,0)&
                                       , -1
                                    ansats(10,5,1) = plusa5
                                    ansats(10,5,0) = ana5 - plusa5
!     10i
                                    call slug (10, 6, varmax, varupp, varned, &
                                       ansats, org, lock(10,6), dubbel, low, &
                                       start(10,6), stopp(10,6))
                                    do ana6 = start(10,6), stopp(10,6), steg(10&
                                       ,6)
                                    antel(10,6) = ana6 + antel(10,5)
                                    if (.not.(antel(10,6)<=antal .and. ansats(&
                                       10,5,1)<=2 .and. ansats(10,5,0)<=2)) &
                                       cycle
                                    do plusa6 = min(ana6,14), max(ana6 - 12,0)&
                                       , -1
                                    ansats(10,6,1) = plusa6
                                    ansats(10,6,0) = ana6 - plusa6
!     10k
                                    call slug (10, 7, varmax, varupp, varned, &
                                       ansats, org, lock(10,7), dubbel, low, &
                                       start(10,7), stopp(10,7))
                                    do ana7 = start(10,7), stopp(10,7), steg(10&
                                       ,7)
                                    antel(10,7) = ana7 + antel(10,6)
                                    if (.not.(antel(10,7)<=antal .and. ansats(&
                                       10,6,1)<=2 .and. ansats(10,6,0)<=2)) &
                                       cycle
                                    do plusa7 = min(ana7,16), max(ana7 - 14,0)&
                                       , -1
                                    ansats(10,7,1) = plusa7
                                    ansats(10,7,0) = ana7 - plusa7
!     10l
                                    call slug (10, 8, varmax, varupp, varned, &
                                       ansats, org, lock(10,8), dubbel, low, &
                                       start(10,8), stopp(10,8))
                                    do ana8 = start(10,8), stopp(10,8), steg(10&
                                       ,8)
                                    antel(10,8) = ana8 + antel(10,7)
                                    if (.not.(antel(10,8)<=antal .and. ansats(&
                                       10,7,1)<=2 .and. ansats(10,7,0)<=2)) &
                                       cycle
                                    do plusa8 = min(ana8,18), max(ana8 - 16,0)&
                                       , -1
                                    ansats(10,8,1) = plusa8
                                    ansats(10,8,0) = ana8 - plusa8
!     10m
                                    call slug (10, 9, varmax, varupp, varned, &
                                       ansats, org, lock(10,9), dubbel, low, &
                                       start(10,9), stopp(10,9))
                                    do ana9 = start(10,9), stopp(10,9), steg(10&
                                       ,9)
                                    antel(10,9) = ana9 + antel(10,8)
                                    if (.not.(antel(10,9)<=antal .and. ansats(&
                                       10,8,1)<=2 .and. ansats(10,8,0)<=2&
                                        .and. antel(10,9)>=lim(10))) cycle
                                    do plusa9 = min(ana9,20), max(ana9 - 18,0)&
                                       , -1
                                    ansats(10,9,1) = plusa9
                                    ansats(10,9,0) = ana9 - plusa9
!     11s
                                    call slug (11, 0, varmax, varupp, varned, &
                                       ansats, org, lock(11,0), dubbel, low, &
                                       start(11,0), stopp(11,0))
                                    do anb0 = start(11,0), stopp(11,0), steg(11&
                                       ,0)
                                    antel(11,0) = anb0 + antel(10,9)
                                    if (.not.(antel(11,0)<=antal .and. ansats(&
                                       10,9,1)<=2 .and. ansats(10,9,0)<=2)) &
                                       cycle
                                    ansats(11,0,0) = anb0
!     11p
                                    call slug (11, 1, varmax, varupp, varned, &
                                       ansats, org, lock(11,1), dubbel, low, &
                                       start(11,1), stopp(11,1))
                                    do anb1 = start(11,1), stopp(11,1), steg(11&
                                       ,1)
                                    antel(11,1) = anb1 + antel(11,0)
                                    if (antel(11,1) > antal) cycle
                                    do plusb1 = min(anb1,4), max(anb1 - 2,0), &
                                       -1
                                    ansats(11,1,1) = plusb1
                                    ansats(11,1,0) = anb1 - plusb1
!     11d
                                    call slug (11, 2, varmax, varupp, varned, &
                                       ansats, org, lock(11,2), dubbel, low, &
                                       start(11,2), stopp(11,2))
                                    do anb2 = start(11,2), stopp(11,2), steg(11&
                                       ,2)
                                    antel(11,2) = anb2 + antel(11,1)
                                    if (antel(11,2) > antal) cycle
                                    do plusb2 = min(anb2,6), max(anb2 - 4,0), &
                                       -1
                                    ansats(11,2,1) = plusb2
                                    ansats(11,2,0) = anb2 - plusb2
!     11f
                                    call slug (11, 3, varmax, varupp, varned, &
                                       ansats, org, lock(11,3), dubbel, low, &
                                       start(11,3), stopp(11,3))
                                    do anb3 = start(11,3), stopp(11,3), steg(11&
                                       ,3)
                                    antel(11,3) = anb3 + antel(11,2)
                                    if (antel(11,3) > antal) cycle
                                    do plusb3 = min(anb3,8), max(anb3 - 6,0), &
                                       -1
                                    ansats(11,3,1) = plusb3
                                    ansats(11,3,0) = anb3 - plusb3
!     11g
                                    call slug (11, 4, varmax, varupp, varned, &
                                       ansats, org, lock(11,4), dubbel, low, &
                                       start(11,4), stopp(11,4))
                                    do anb4 = start(11,4), stopp(11,4), steg(11&
                                       ,4)
                                    antel(11,4) = anb4 + antel(11,3)
                                    if (antel(11,4) > antal) cycle
                                    do plusb4 = min(anb4,10), max(anb4 - 8,0), &
                                       -1
                                    ansats(11,4,1) = plusb4
                                    ansats(11,4,0) = anb4 - plusb4
!     11h
                                    call slug (11, 5, varmax, varupp, varned, &
                                       ansats, org, lock(11,5), dubbel, low, &
                                       start(11,5), stopp(11,5))
                                    do anb5 = start(11,5), stopp(11,5), steg(11&
                                       ,5)
                                    antel(11,5) = anb5 + antel(11,4)
                                    if (antel(11,5)>antal .or. ansats(11,4,1)>2&
                                       ) cycle
                                    do plusb5 = min(anb5,12), max(anb5 - 10,0)&
                                       , -1
                                    ansats(11,5,1) = plusb5
                                    ansats(11,5,0) = anb5 - plusb5
!     11i
                                    call slug (11, 6, varmax, varupp, varned, &
                                       ansats, org, lock(11,6), dubbel, low, &
                                       start(11,6), stopp(11,6))
                                    do anb6 = start(11,6), stopp(11,6), steg(11&
                                       ,6)
                                    antel(11,6) = anb6 + antel(11,5)
                                    if (.not.(antel(11,6)<=antal .and. ansats(&
                                       11,5,1)<=2 .and. ansats(11,5,0)<=2)) &
                                       cycle
                                    do plusb6 = min(anb6,14), max(anb6 - 12,0)&
                                       , -1
                                    ansats(11,6,1) = plusb6
                                    ansats(11,6,0) = anb6 - plusb6
!     11k
                                    call slug (11, 7, varmax, varupp, varned, &
                                       ansats, org, lock(11,7), dubbel, low, &
                                       start(11,7), stopp(11,7))
                                    do anb7 = start(11,7), stopp(11,7), steg(11&
                                       ,7)
                                    antel(11,7) = anb7 + antel(11,6)
                                    if (.not.(antel(11,7)<=antal .and. ansats(&
                                       11,6,1)<=2 .and. ansats(11,6,0)<=2)) &
                                       cycle
                                    do plusb7 = min(anb7,16), max(anb7 - 14,0)&
                                       , -1
                                    ansats(11,7,1) = plusb7
                                    ansats(11,7,0) = anb7 - plusb7
!     11l
                                    call slug (11, 8, varmax, varupp, varned, &
                                       ansats, org, lock(11,8), dubbel, low, &
                                       start(11,8), stopp(11,8))
                                    do anb8 = start(11,8), stopp(11,8), steg(11&
                                       ,8)
                                    antel(11,8) = anb8 + antel(11,7)
                                    if (.not.(antel(11,8)<=antal .and. ansats(&
                                       11,7,1)<=2 .and. ansats(11,7,0)<=2)) &
                                       cycle
                                    do plusb8 = min(anb8,18), max(anb8 - 16,0)&
                                       , -1
                                    ansats(11,8,1) = plusb8
                                    ansats(11,8,0) = anb8 - plusb8
!     11m
                                    call slug (11, 9, varmax, varupp, varned, &
                                       ansats, org, lock(11,9), dubbel, low, &
                                       start(11,9), stopp(11,9))
                                    do anb9 = start(11,9), stopp(11,9), steg(11&
                                       ,9)
                                    antel(11,9) = anb9 + antel(11,8)
                                    if (.not.(antel(11,9)<=antal .and. ansats(&
                                       11,8,1)<=2 .and. ansats(11,8,0)<=2)) &
                                       cycle
                                    do plusb9 = min(anb9,20), max(anb9 - 18,0)&
                                       , -1
                                    ansats(11,9,1) = plusb9
                                    ansats(11,9,0) = anb9 - plusb9
!     11n
                                    call slug (11, 10, varmax, varupp, varned, &
                                       ansats, org, lock(11,10), dubbel, low, &
                                       start(11,10), stopp(11,10))
                                    do anba = start(11,10), stopp(11,10), steg(&
                                       11,10)
                                    antel(11,10) = anba + antel(11,9)
                                    if (.not.(antel(11,10)<=antal .and. ansats(&
                                       11,9,1)<=2 .and. ansats(11,9,0)<=2&
                                        .and. antel(11,10)>=lim(11))) cycle
                                    do plusba = min(anba,22), max(anba - 20,0)&
                                       , -1
                                    ansats(11,10,1) = plusba
                                    ansats(11,10,0) = anba - plusba
!     12s
                                    call slug (12, 0, varmax, varupp, varned, &
                                       ansats, org, lock(12,0), dubbel, low, &
                                       start(12,0), stopp(12,0))
                                    do anc0 = start(12,0), stopp(12,0), steg(12&
                                       ,0)
                                    antel(12,0) = anc0 + antel(11,10)
                                    if (.not.(antel(12,0)<=antal .and. ansats(&
                                       11,10,1)<=2 .and. ansats(11,10,0)<=2)) &
                                       cycle
                                    ansats(12,0,0) = anc0
!     12p
                                    call slug (12, 1, varmax, varupp, varned, &
                                       ansats, org, lock(12,1), dubbel, low, &
                                       start(12,1), stopp(12,1))
                                    do anc1 = start(12,1), stopp(12,1), steg(12&
                                       ,1)
                                    antel(12,1) = anc1 + antel(12,0)
                                    if (antel(12,1) > antal) cycle
                                    do plusc1 = min(anc1,4), max(anc1 - 2,0), &
                                       -1
                                    ansats(12,1,1) = plusc1
                                    ansats(12,1,0) = anc1 - plusc1
!     12d
                                    call slug (12, 2, varmax, varupp, varned, &
                                       ansats, org, lock(12,2), dubbel, low, &
                                       start(12,2), stopp(12,2))
                                    do anc2 = start(12,2), stopp(12,2), steg(12&
                                       ,2)
                                    antel(12,2) = anc2 + antel(12,1)
                                    if (antel(12,2) > antal) cycle
                                    do plusc2 = min(anc2,6), max(anc2 - 4,0), &
                                       -1
                                    ansats(12,2,1) = plusc2
                                    ansats(12,2,0) = anc2 - plusc2
!     12f
                                    call slug (12, 3, varmax, varupp, varned, &
                                       ansats, org, lock(12,3), dubbel, low, &
                                       start(12,3), stopp(12,3))
                                    do anc3 = start(12,3), stopp(12,3), steg(12&
                                       ,3)
                                    antel(12,3) = anc3 + antel(12,2)
                                    if (antel(12,3) > antal) cycle
                                    do plusc3 = min(anc3,8), max(anc3 - 6,0), &
                                       -1
                                    ansats(12,3,1) = plusc3
                                    ansats(12,3,0) = anc3 - plusc3
!     12g
                                    call slug (12, 4, varmax, varupp, varned, &
                                       ansats, org, lock(12,4), dubbel, low, &
                                       start(12,4), stopp(12,4))
                                    do anc4 = start(12,4), stopp(12,4), steg(12&
                                       ,4)
                                    antel(12,4) = anc4 + antel(12,3)
                                    if (antel(12,4) > antal) cycle
                                    do plusc4 = min(anc4,10), max(anc4 - 8,0), &
                                       -1
                                    ansats(12,4,1) = plusc4
                                    ansats(12,4,0) = anc4 - plusc4
!     12h
                                    call slug (12, 5, varmax, varupp, varned, &
                                       ansats, org, lock(12,5), dubbel, low, &
                                       start(12,5), stopp(12,5))
                                    do anc5 = start(12,5), stopp(12,5), steg(12&
                                       ,5)
                                    antel(12,5) = anc5 + antel(12,4)
                                    if (antel(12,5)>antal .or. ansats(12,4,1)>2&
                                       ) cycle
                                    do plusc5 = min(anc5,12), max(anc5 - 10,0)&
                                       , -1
                                    ansats(12,5,1) = plusc5
                                    ansats(12,5,0) = anc5 - plusc5
!     12i
                                    call slug (12, 6, varmax, varupp, varned, &
                                       ansats, org, lock(12,6), dubbel, low, &
                                       start(12,6), stopp(12,6))
                                    do anc6 = start(12,6), stopp(12,6), steg(12&
                                       ,6)
                                    antel(12,6) = anc6 + antel(12,5)
                                    if (.not.(antel(12,6)<=antal .and. ansats(&
                                       12,5,1)<=2 .and. ansats(12,5,0)<=2)) &
                                       cycle
                                    do plusc6 = min(anc6,14), max(anc6 - 12,0)&
                                       , -1
                                    ansats(12,6,1) = plusc6
                                    ansats(12,6,0) = anc6 - plusc6
!     12k
                                    call slug (12, 7, varmax, varupp, varned, &
                                       ansats, org, lock(12,7), dubbel, low, &
                                       start(12,7), stopp(12,7))
                                    do anc7 = start(12,7), stopp(12,7), steg(12&
                                       ,7)
                                    antel(12,7) = anc7 + antel(12,6)
                                    if (.not.(antel(12,7)<=antal .and. ansats(&
                                       12,6,1)<=2 .and. ansats(12,6,0)<=2)) &
                                       cycle
                                    do plusc7 = min(anc7,16), max(anc7 - 14,0)&
                                       , -1
                                    ansats(12,7,1) = plusc7
                                    ansats(12,7,0) = anc7 - plusc7
!     12l
                                    call slug (12, 8, varmax, varupp, varned, &
                                       ansats, org, lock(12,8), dubbel, low, &
                                       start(12,8), stopp(12,8))
                                    do anc8 = start(12,8), stopp(12,8), steg(12&
                                       ,8)
                                    antel(12,8) = anc8 + antel(12,7)
                                    if (.not.(antel(12,8)<=antal .and. ansats(&
                                       12,7,1)<=2 .and. ansats(12,7,0)<=2)) &
                                       cycle
                                    do plusc8 = min(anc8,18), max(anc8 - 16,0)&
                                       , -1
                                    ansats(12,8,1) = plusc8
                                    ansats(12,8,0) = anc8 - plusc8
!     12m
                                    call slug (12, 9, varmax, varupp, varned, &
                                       ansats, org, lock(12,9), dubbel, low, &
                                       start(12,9), stopp(12,9))
                                    do anc9 = start(12,9), stopp(12,9), steg(12&
                                       ,9)
                                    antel(12,9) = anc9 + antel(12,8)
                                    if (.not.(antel(12,9)<=antal .and. ansats(&
                                       12,8,1)<=2 .and. ansats(12,8,0)<=2)) &
                                       cycle
                                    do plusc9 = min(anc9,20), max(anc9 - 18,0)&
                                       , -1
                                    ansats(12,9,1) = plusc9
                                    ansats(12,9,0) = anc9 - plusc9
!     12n
                                    call slug (12, 10, varmax, varupp, varned, &
                                       ansats, org, lock(12,10), dubbel, low, &
                                       start(12,10), stopp(12,10))
                                    do anca = start(12,10), stopp(12,10), steg(&
                                       12,10)
                                    antel(12,10) = anca + antel(12,9)
                                    if (.not.(antel(12,10)<=antal .and. ansats(&
                                       12,9,1)<=2 .and. ansats(12,9,0)<=2&
                                        .and. antel(12,10)>=lim(12))) cycle
                                    do plusca = min(anca,22), max(anca - 20,0)&
                                       , -1
                                    ansats(12,10,1) = plusca
                                    ansats(12,10,0) = anca - plusca
!     13s
                                    call slug (13, 0, varmax, varupp, varned, &
                                       ansats, org, lock(13,0), dubbel, low, &
                                       start(13,0), stopp(13,0))
                                    do and0 = start(13,0), stopp(13,0), steg(13&
                                       ,0)
                                    antel(13,0) = and0 + antel(12,10)
                                    if (.not.(antel(13,0)<=antal .and. ansats(&
                                       12,10,1)<=2 .and. ansats(12,10,0)<=2)) &
                                       cycle
                                    ansats(13,0,0) = and0
!     13p
                                    call slug (13, 1, varmax, varupp, varned, &
                                       ansats, org, lock(13,1), dubbel, low, &
                                       start(13,1), stopp(13,1))
                                    do and1 = start(13,1), stopp(13,1), steg(13&
                                       ,1)
                                    antel(13,1) = and1 + antel(13,0)
                                    if (antel(13,1) > antal) cycle
                                    do plusd1 = min(and1,4), max(and1 - 2,0), &
                                       -1
                                    ansats(13,1,1) = plusd1
                                    ansats(13,1,0) = and1 - plusd1
!     13d
                                    call slug (13, 2, varmax, varupp, varned, &
                                       ansats, org, lock(13,2), dubbel, low, &
                                       start(13,2), stopp(13,2))
                                    do and2 = start(13,2), stopp(13,2), steg(13&
                                       ,2)
                                    antel(13,2) = and2 + antel(13,1)
                                    if (antel(13,2) > antal) cycle
                                    do plusd2 = min(and2,6), max(and2 - 4,0), &
                                       -1
                                    ansats(13,2,1) = plusd2
                                    ansats(13,2,0) = and2 - plusd2
!     13f
                                    call slug (13, 3, varmax, varupp, varned, &
                                       ansats, org, lock(13,3), dubbel, low, &
                                       start(13,3), stopp(13,3))
                                    do and3 = start(13,3), stopp(13,3), steg(13&
                                       ,3)
                                    antel(13,3) = and3 + antel(13,2)
                                    if (antel(13,3) > antal) cycle
                                    do plusd3 = min(and3,8), max(and3 - 6,0), &
                                       -1
                                    ansats(13,3,1) = plusd3
                                    ansats(13,3,0) = and3 - plusd3
!     13g
                                    call slug (13, 4, varmax, varupp, varned, &
                                       ansats, org, lock(13,4), dubbel, low, &
                                       start(13,4), stopp(13,4))
                                    do and4 = start(13,4), stopp(13,4), steg(13&
                                       ,4)
                                    antel(13,4) = and4 + antel(13,3)
                                    if (antel(13,4) > antal) cycle
                                    do plusd4 = min(and4,10), max(and4 - 8,0), &
                                       -1
                                    ansats(13,4,1) = plusd4
                                    ansats(13,4,0) = and4 - plusd4
!     13h
                                    call slug (13, 5, varmax, varupp, varned, &
                                       ansats, org, lock(13,5), dubbel, low, &
                                       start(13,5), stopp(13,5))
                                    do and5 = start(13,5), stopp(13,5), steg(13&
                                       ,5)
                                    antel(13,5) = and5 + antel(13,4)
                                    if (antel(13,5)>antal .or. ansats(13,4,1)>2&
                                       ) cycle
                                    do plusd5 = min(and5,12), max(and5 - 10,0)&
                                       , -1
                                    ansats(13,5,1) = plusd5
                                    ansats(13,5,0) = and5 - plusd5
!     13i
                                    call slug (13, 6, varmax, varupp, varned, &
                                       ansats, org, lock(13,6), dubbel, low, &
                                       start(13,6), stopp(13,6))
                                    do and6 = start(13,6), stopp(13,6), steg(13&
                                       ,6)
                                    antel(13,6) = and6 + antel(13,5)
                                    if (.not.(antel(13,6)<=antal .and. ansats(&
                                       13,5,1)<=2 .and. ansats(13,5,0)<=2)) &
                                       cycle
                                    do plusd6 = min(and6,14), max(and6 - 12,0)&
                                       , -1
                                    ansats(13,6,1) = plusd6
                                    ansats(13,6,0) = and6 - plusd6
!     13k
                                    call slug (13, 7, varmax, varupp, varned, &
                                       ansats, org, lock(13,7), dubbel, low, &
                                       start(13,7), stopp(13,7))
                                    do and7 = start(13,7), stopp(13,7), steg(13&
                                       ,7)
                                    antel(13,7) = and7 + antel(13,6)
                                    if (.not.(antel(13,7)<=antal .and. ansats(&
                                       13,6,1)<=2 .and. ansats(13,6,0)<=2)) &
                                       cycle
                                    do plusd7 = min(and7,16), max(and7 - 14,0)&
                                       , -1
                                    ansats(13,7,1) = plusd7
                                    ansats(13,7,0) = and7 - plusd7
!     13l
                                    call slug (13, 8, varmax, varupp, varned, &
                                       ansats, org, lock(13,8), dubbel, low, &
                                       start(13,8), stopp(13,8))
                                    do and8 = start(13,8), stopp(13,8), steg(13&
                                       ,8)
                                    antel(13,8) = and8 + antel(13,7)
                                    if (.not.(antel(13,8)<=antal .and. ansats(&
                                       13,7,1)<=2 .and. ansats(13,7,0)<=2)) &
                                       cycle
                                    do plusd8 = min(and8,18), max(and8 - 16,0)&
                                       , -1
                                    ansats(13,8,1) = plusd8
                                    ansats(13,8,0) = and8 - plusd8
!     13m
                                    call slug (13, 9, varmax, varupp, varned, &
                                       ansats, org, lock(13,9), dubbel, low, &
                                       start(13,9), stopp(13,9))
                                    do and9 = start(13,9), stopp(13,9), steg(13&
                                       ,9)
                                    antel(13,9) = and9 + antel(13,8)
                                    if (.not.(antel(13,9)<=antal .and. ansats(&
                                       13,8,1)<=2 .and. ansats(13,8,0)<=2)) &
                                       cycle
                                    do plusd9 = min(and9,20), max(and9 - 18,0)&
                                       , -1
                                    ansats(13,9,1) = plusd9
                                    ansats(13,9,0) = and9 - plusd9
!     13n
                                    call slug (13, 10, varmax, varupp, varned, &
                                       ansats, org, lock(13,10), dubbel, low, &
                                       start(13,10), stopp(13,10))
                                    do anda = start(13,10), stopp(13,10), steg(&
                                       13,10)
                                    antel(13,10) = anda + antel(13,9)
                                    if (.not.(antel(13,10)<=antal .and. ansats(&
                                       13,9,1)<=2 .and. ansats(13,9,0)<=2&
                                        .and. antel(13,10)>=lim(13))) cycle
                                    do plusda = min(anda,22), max(anda - 20,0)&
                                       , -1
                                    ansats(13,10,1) = plusda
                                    ansats(13,10,0) = anda - plusda
!     14s
                                    call slug (14, 0, varmax, varupp, varned, &
                                       ansats, org, lock(14,0), dubbel, low, &
                                       start(14,0), stopp(14,0))
                                    do ane0 = start(14,0), stopp(14,0), steg(14&
                                       ,0)
                                    antel(14,0) = ane0 + antel(13,10)
                                    if (.not.(antel(14,0)<=antal .and. ansats(&
                                       13,10,1)<=2 .and. ansats(13,10,0)<=2)) &
                                       cycle
                                    ansats(14,0,0) = ane0
!     14p
                                    call slug (14, 1, varmax, varupp, varned, &
                                       ansats, org, lock(14,1), dubbel, low, &
                                       start(14,1), stopp(14,1))
                                    do ane1 = start(14,1), stopp(14,1), steg(14&
                                       ,1)
                                    antel(14,1) = ane1 + antel(14,0)
                                    if (antel(14,1) > antal) cycle
                                    do pluse1 = min(ane1,4), max(ane1 - 2,0), &
                                       -1
                                    ansats(14,1,1) = pluse1
                                    ansats(14,1,0) = ane1 - pluse1
!     14d
                                    call slug (14, 2, varmax, varupp, varned, &
                                       ansats, org, lock(14,2), dubbel, low, &
                                       start(14,2), stopp(14,2))
                                    do ane2 = start(14,2), stopp(14,2), steg(14&
                                       ,2)
                                    antel(14,2) = ane2 + antel(14,1)
                                    if (antel(14,2) > antal) cycle
                                    do pluse2 = min(ane2,6), max(ane2 - 4,0), &
                                       -1
                                    ansats(14,2,1) = pluse2
                                    ansats(14,2,0) = ane2 - pluse2
!     14f
                                    call slug (14, 3, varmax, varupp, varned, &
                                       ansats, org, lock(14,3), dubbel, low, &
                                       start(14,3), stopp(14,3))
                                    do ane3 = start(14,3), stopp(14,3), steg(14&
                                       ,3)
                                    antel(14,3) = ane3 + antel(14,2)
                                    if (antel(14,3) > antal) cycle
                                    do pluse3 = min(ane3,8), max(ane3 - 6,0), &
                                       -1
                                    ansats(14,3,1) = pluse3
                                    ansats(14,3,0) = ane3 - pluse3
!     14g
                                    call slug (14, 4, varmax, varupp, varned, &
                                       ansats, org, lock(14,4), dubbel, low, &
                                       start(14,4), stopp(14,4))
                                    do ane4 = start(14,4), stopp(14,4), steg(14&
                                       ,4)
                                    antel(14,4) = ane4 + antel(14,3)
                                    if (antel(14,4) > antal) cycle
                                    do pluse4 = min(ane4,10), max(ane4 - 8,0), &
                                       -1
                                    ansats(14,4,1) = pluse4
                                    ansats(14,4,0) = ane4 - pluse4
!     14h
                                    call slug (14, 5, varmax, varupp, varned, &
                                       ansats, org, lock(14,5), dubbel, low, &
                                       start(14,5), stopp(14,5))
                                    do ane5 = start(14,5), stopp(14,5), steg(14&
                                       ,5)
                                    antel(14,5) = ane5 + antel(14,4)
                                    if (antel(14,5)>antal .or. ansats(14,4,1)>2&
                                       ) cycle
                                    do pluse5 = min(ane5,12), max(ane5 - 10,0)&
                                       , -1
                                    ansats(14,5,1) = pluse5
                                    ansats(14,5,0) = ane5 - pluse5
!     14i
                                    call slug (14, 6, varmax, varupp, varned, &
                                       ansats, org, lock(14,6), dubbel, low, &
                                       start(14,6), stopp(14,6))
                                    do ane6 = start(14,6), stopp(14,6), steg(14&
                                       ,6)
                                    antel(14,6) = ane6 + antel(14,5)
                                    if (.not.(antel(14,6)<=antal .and. ansats(&
                                       14,5,1)<=2 .and. ansats(14,5,0)<=2)) &
                                       cycle
                                    do pluse6 = min(ane6,14), max(ane6 - 12,0)&
                                       , -1
                                    ansats(14,6,1) = pluse6
                                    ansats(14,6,0) = ane6 - pluse6
!     14k
                                    call slug (14, 7, varmax, varupp, varned, &
                                       ansats, org, lock(14,7), dubbel, low, &
                                       start(14,7), stopp(14,7))
                                    do ane7 = start(14,7), stopp(14,7), steg(14&
                                       ,7)
                                    antel(14,7) = ane7 + antel(14,6)
                                    if (.not.(antel(14,7)<=antal .and. ansats(&
                                       14,6,1)<=2 .and. ansats(14,6,0)<=2)) &
                                       cycle
                                    do pluse7 = min(ane7,16), max(ane7 - 14,0)&
                                       , -1
                                    ansats(14,7,1) = pluse7
                                    ansats(14,7,0) = ane7 - pluse7
!     14l
                                    call slug (14, 8, varmax, varupp, varned, &
                                       ansats, org, lock(14,8), dubbel, low, &
                                       start(14,8), stopp(14,8))
                                    do ane8 = start(14,8), stopp(14,8), steg(14&
                                       ,8)
                                    antel(14,8) = ane8 + antel(14,7)
                                    if (.not.(antel(14,8)<=antal .and. ansats(&
                                       14,7,1)<=2 .and. ansats(14,7,0)<=2)) &
                                       cycle
                                    do pluse8 = min(ane8,18), max(ane8 - 16,0)&
                                       , -1
                                    ansats(14,8,1) = pluse8
                                    ansats(14,8,0) = ane8 - pluse8
!     14m
                                    call slug (14, 9, varmax, varupp, varned, &
                                       ansats, org, lock(14,9), dubbel, low, &
                                       start(14,9), stopp(14,9))
                                    do ane9 = start(14,9), stopp(14,9), steg(14&
                                       ,9)
                                    antel(14,9) = ane9 + antel(14,8)
                                    if (.not.(antel(14,9)<=antal .and. ansats(&
                                       14,8,1)<=2 .and. ansats(14,8,0)<=2)) &
                                       cycle
                                    do pluse9 = min(ane9,20), max(ane9 - 18,0)&
                                       , -1
                                    ansats(14,9,1) = pluse9
                                    ansats(14,9,0) = ane9 - pluse9
!     14n
                                    call slug (14, 10, varmax, varupp, varned, &
                                       ansats, org, lock(14,10), dubbel, low, &
                                       start(14,10), stopp(14,10))
                                    do anea = start(14,10), stopp(14,10), steg(&
                                       14,10)
                                    antel(14,10) = anea + antel(14,9)
                                    if (.not.(antel(14,10)<=antal .and. ansats(&
                                       14,9,1)<=2 .and. ansats(14,9,0)<=2&
                                        .and. antel(14,10)>=lim(14))) cycle
                                    do plusea = min(anea,22), max(anea - 20,0)&
                                       , -1
                                    ansats(14,10,1) = plusea
                                    ansats(14,10,0) = anea - plusea
!     15s
                                    call slug (15, 0, varmax, varupp, varned, &
                                       ansats, org, lock(15,0), dubbel, low, &
                                       start(15,0), stopp(15,0))
                                    do anf0 = start(15,0), stopp(15,0), steg(15&
                                       ,0)
                                        
                                    antel(15,0) = anf0 + antel(14,10)
                                    if (.not.(antel(15,0)<=antal .and. ansats(&
                                       14,10,1)<=2 .and. ansats(14,10,0)<=2)) &
                                       cycle
                                    ansats(15,0,0) = anf0
!     15p
                                    call slug (15, 1, varmax, varupp, varned, &
                                       ansats, org, lock(15,1), dubbel, low, &
                                       start(15,1), stopp(15,1))
                                    do anf1 = start(15,1), stopp(15,1), steg(15&
                                       ,1)
                                    antel(15,1) = anf1 + antel(15,0)
                                    if (antel(15,1) > antal) cycle
                                    do plusf1 = min(anf1,4), max(anf1 - 2,0), &
                                       -1
                                    ansats(15,1,1) = plusf1
                                    ansats(15,1,0) = anf1 - plusf1
!     15d
                                    call slug (15, 2, varmax, varupp, varned, &
                                       ansats, org, lock(15,2), dubbel, low, &
                                       start(15,2), stopp(15,2))
                                    do anf2 = start(15,2), stopp(15,2), steg(15&
                                       ,2)
                                    antel(15,2) = anf2 + antel(15,1)
                                    if (antel(15,2) > antal) cycle
                                    do plusf2 = min(anf2,6), max(anf2 - 4,0), &
                                       -1
                                    ansats(15,2,1) = plusf2
                                    ansats(15,2,0) = anf2 - plusf2
!     15f
                                    call slug (15, 3, varmax, varupp, varned, &
                                       ansats, org, lock(15,3), dubbel, low, &
                                       start(15,3), stopp(15,3))
                                    do anf3 = start(15,3), stopp(15,3), steg(15&
                                       ,3)
                                    antel(15,3) = anf3 + antel(15,2)
                                    if (antel(15,3) > antal) cycle
                                    do plusf3 = min(anf3,8), max(anf3 - 6,0), &
                                       -1
                                    ansats(15,3,1) = plusf3
                                    ansats(15,3,0) = anf3 - plusf3
!     15g
                                    call slug (15, 4, varmax, varupp, varned, &
                                       ansats, org, lock(15,4), dubbel, low, &
                                       start(15,4), stopp(15,4))
                                    do anf4 = start(15,4), stopp(15,4), steg(15&
                                       ,4)
                                    antel(15,4) = anf4 + antel(15,3)
                                    if (antel(15,4) > antal) cycle
                                    do plusf4 = min(anf4,10), max(anf4 - 8,0), &
                                       -1
                                    ansats(15,4,1) = plusf4
                                    ansats(15,4,0) = anf4 - plusf4
!     15h
                                    call slug (15, 5, varmax, varupp, varned, &
                                       ansats, org, lock(15,5), dubbel, low, &
                                       start(15,5), stopp(15,5))
                                    do anf5 = start(15,5), stopp(15,5), steg(15&
                                       ,5)
                                    antel(15,5) = anf5 + antel(15,4)
                                    if (antel(15,5)>antal .or. ansats(15,4,1)>2&
                                       ) cycle
                                    do plusf5 = min(anf5,12), max(anf5 - 10,0)&
                                       , -1
                                    ansats(15,5,1) = plusf5
                                    ansats(15,5,0) = anf5 - plusf5
      !if (plusf5.ne.0) write(*,*)'Now 15h plusf5=', plusf5
!     15i
                                    call slug (15, 6, varmax, varupp, varned, &
                                       ansats, org, lock(15,6), dubbel, low, &
                                       start(15,6), stopp(15,6))
                                    do anf6 = start(15,6), stopp(15,6), steg(15&
                                       ,6)
                                    antel(15,6) = anf6 + antel(15,5)
                                    if (.not.(antel(15,6)<=antal .and. ansats(&
                                       15,5,1)<=2 .and. ansats(15,5,0)<=2)) &
                                       cycle
                                    do plusf6 = min(anf6,14), max(anf6 - 12,0)&
                                       , -1
                                    ansats(15,6,1) = plusf6
                                    ansats(15,6,0) = anf6 - plusf6
!     15k
                                    call slug (15, 7, varmax, varupp, varned, &
                                       ansats, org, lock(15,7), dubbel, low, &
                                       start(15,7), stopp(15,7))
                                    do anf7 = start(15,7), stopp(15,7), steg(15&
                                       ,7)
                                    antel(15,7) = anf7 + antel(15,6)
                                    if (.not.(antel(15,7)<=antal .and. ansats(&
                                       15,6,1)<=2 .and. ansats(15,6,0)<=2)) &
                                       cycle
                                    do plusf7 = min(anf7,16), max(anf7 - 14,0)&
                                       , -1
                                    ansats(15,7,1) = plusf7
                                    ansats(15,7,0) = anf7 - plusf7
!     15l
                                    call slug (15, 8, varmax, varupp, varned, &
                                       ansats, org, lock(15,8), dubbel, low, &
                                       start(15,8), stopp(15,8))
                                    do anf8 = start(15,8), stopp(15,8), steg(15&
                                       ,8)
                                    antel(15,8) = anf8 + antel(15,7)
                                    if (.not.(antel(15,8)<=antal .and. ansats(&
                                       15,7,1)<=2 .and. ansats(15,7,0)<=2)) &
                                       cycle
                                    do plusf8 = min(anf8,18), max(anf8 - 16,0)&
                                       , -1
                                    ansats(15,8,1) = plusf8
                                    ansats(15,8,0) = anf8 - plusf8
!     15m
                                    call slug (15, 9, varmax, varupp, varned, &
                                       ansats, org, lock(15,9), dubbel, low, &
                                       start(15,9), stopp(15,9))
                                    do anf9 = start(15,9), stopp(15,9), steg(15&
                                       ,9)
                                    antel(15,9) = anf9 + antel(15,8)
                                    if (.not.(antel(15,9)<=antal .and. ansats(&
                                       15,8,1)<=2 .and. ansats(15,8,0)<=2)) &
                                       cycle
                                    do plusf9 = min(anf9,20), max(anf9 - 18,0)&
                                       , -1
                                    ansats(15,9,1) = plusf9
                                    ansats(15,9,0) = anf9 - plusf9
!     15n
                                    call slug (15, 10, varmax, varupp, varned, &
                                       ansats, org, lock(15,10), dubbel, low, &
                                       start(15,10), stopp(15,10))
                                    do anfa = start(15,10), stopp(15,10), steg(&
                                       15,10)
                                    antel(15,10) = anfa + antel(15,9)
                                    !if (.not.(antel(15,10)==antal .and. ansats(&
                                    if (.not.(antel(15,10)<=antal .and. ansats(&
                                       15,9,1)<=2 .and. ansats(15,9,0)<=2)) &
                                       cycle
                                    do plusfa = min(anfa,22), max(anfa - 20,0)&
                                       , -1
                                    ansats(15,10,1) = plusfa
                                    ansats(15,10,0) = anfa - plusfa
                                    !if (ansats(15,10,1)>2 .or. ansats(15,10,0)>&
                                    !   2) cycle
! n = 16
                                    ng16 = 16
!     16s
                                    call slug (ng16, 0, varmax, varupp, varned, &
                                       ansats, org, lock(ng16,0), dubbel, low, &
                                       start(ng16,0), stopp(ng16,0))
                                    do ang0 = start(ng16,0), stopp(ng16,0), steg(ng16,0)
                                    antel(ng16,0) = ang0 + antel(ng16-1,10)
                                    if (.not.(antel(ng16,0)<=antal .and. ansats(&
                                       ng16-1,10,1)<=2 .and. ansats(ng16-1,10,0)<=2)) &
                                       cycle
                                    ansats(ng16,0,0) = ang0
!     16p
                                    call slug (ng16, 1, varmax, varupp, varned, &
                                       ansats, org, lock(ng16,1), dubbel, low, &
                                       start(ng16,1), stopp(ng16,1))
                                    do ang1 = start(ng16,1), stopp(ng16,1), steg(ng16,1)
                                    antel(ng16,1) = ang1 + antel(ng16,0)
                                    if (antel(ng16,1) > antal) cycle
                                    do plusg1 = min(ang1,4), max(ang1 - 2,0),-1
                                    ansats(ng16,1,1) = plusg1
                                    ansats(ng16,1,0) = ang1 - plusg1
!     16d
                                    call slug (ng16, 2, varmax, varupp, varned, &
                                       ansats, org, lock(ng16,2), dubbel, low, &
                                       start(ng16,2), stopp(ng16,2))
                                    do ang2 = start(ng16,2), stopp(ng16,2), steg(ng16,2)
                                    antel(ng16,2) = ang2 + antel(ng16,1)
                                    if (antel(ng16,2) > antal) cycle
                                    do plusg2 = min(ang2,6), max(ang2 - 4,0),-1
                                    ansats(ng16,2,1) = plusg2
                                    ansats(ng16,2,0) = ang2 - plusg2
!     16f
                                    call slug (ng16, 3, varmax, varupp, varned, &
                                       ansats, org, lock(ng16,3), dubbel, low, &
                                       start(ng16,3), stopp(ng16,3))
                                    do ang3 = start(ng16,3), stopp(ng16,3), steg(ng16,3)
                                    antel(ng16,3) = ang3 + antel(ng16,2)
                                    if (antel(ng16,3) > antal) cycle
                                    do plusg3 = min(ang3,8), max(ang3 - 6,0),-1
                                    ansats(ng16,3,1) = plusg3
                                    ansats(ng16,3,0) = ang3 - plusg3
!     16g
                                    call slug (ng16, 4, varmax, varupp, varned, &
                                       ansats, org, lock(ng16,4), dubbel, low, &
                                       start(ng16,4), stopp(ng16,4))
                                    do ang4 = start(ng16,4), stopp(ng16,4), steg(ng16,4)
                                    antel(ng16,4) = ang4 + antel(ng16,3)
                                    if (antel(ng16,4) > antal) cycle
                                    do plusg4 = min(ang4,10), max(ang4 - 8,0),-1
                                    ansats(ng16,4,1) = plusg4
                                    ansats(ng16,4,0) = ang4 - plusg4
!     16h
                                    call slug (ng16, 5, varmax, varupp, varned, &
                                       ansats, org, lock(ng16,5), dubbel, low, &
                                       start(ng16,5), stopp(ng16,5))
                                    do ang5 = start(ng16,5), stopp(ng16,5), steg(ng16,5)
                                    antel(ng16,5) = ang5 + antel(ng16,4)
                                    if (antel(ng16,5)>antal .or. ansats(ng16,4,1)>2&
                                       ) cycle
                                    do plusg5 = min(ang5,12), max(ang5 - 10,0),-1
                                    ansats(ng16,5,1) = plusg5
                                    ansats(ng16,5,0) = ang5 - plusg5
!     16i
                                    call slug (ng16, 6, varmax, varupp, varned, &
                                       ansats, org, lock(ng16,6), dubbel, low, &
                                       start(ng16,6), stopp(ng16,6))
                                    do ang6 = start(ng16,6), stopp(ng16,6), steg(ng16,6)
                                    antel(ng16,6) = ang6 + antel(ng16,5)
                                    if (.not.(antel(ng16,6)<=antal .and. ansats(&
                                       ng16,5,1)<=2 .and. ansats(ng16,5,0)<=2)) &
                                       cycle
                                    do plusg6 = min(ang6,14), max(ang6 - 12,0),-1
                                    ansats(ng16,6,1) = plusg6
                                    ansats(ng16,6,0) = ang6 - plusg6
!     16k
                                    call slug (ng16, 7, varmax, varupp, varned, &
                                       ansats, org, lock(ng16,7), dubbel, low, &
                                       start(ng16,7), stopp(ng16,7))
                                    do ang7 = start(ng16,7), stopp(ng16,7), steg(ng16,7)
                                    antel(ng16,7) = ang7 + antel(ng16,6)
                                    if (.not.(antel(ng16,7)<=antal .and. ansats(&
                                       ng16,6,1)<=2 .and. ansats(ng16,6,0)<=2)) &
                                       cycle
                                    do plusg7 = min(ang7,16), max(ang7 - 14,0),-1
                                    ansats(ng16,7,1) = plusg7
                                    ansats(ng16,7,0) = ang7 - plusg7
!     16l
                                    call slug (ng16, 8, varmax, varupp, varned, &
                                       ansats, org, lock(ng16,8), dubbel, low, &
                                       start(ng16,8), stopp(ng16,8))
                                    do ang8 = start(ng16,8), stopp(ng16,8), steg(ng16,8)
                                    antel(ng16,8) = ang8 + antel(ng16,7)
                                    if (.not.(antel(ng16,8)<=antal .and. ansats(&
                                       ng16,7,1)<=2 .and. ansats(ng16,7,0)<=2)) &
                                       cycle
                                    do plusg8 = min(ang8,18), max(ang8 - 16,0),-1
                                    ansats(ng16,8,1) = plusg8
                                    ansats(ng16,8,0) = ang8 - plusg8
!     16m
                                    call slug (ng16, 9, varmax, varupp, varned, &
                                       ansats, org, lock(ng16,9), dubbel, low, &
                                       start(ng16,9), stopp(ng16,9))
                                    do ang9 = start(ng16,9), stopp(ng16,9), steg(ng16,9)
                                    antel(ng16,9) = ang9 + antel(ng16,8)
                                    if (.not.(antel(ng16,9)<=antal .and. ansats(&
                                       ng16,8,1)<=2 .and. ansats(ng16,8,0)<=2)) &
                                       cycle
                                    do plusg9 = min(ang9,20), max(ang9 - 18,0),-1
                                    ansats(ng16,9,1) = plusg9
                                    ansats(ng16,9,0) = ang9 - plusg9
!     16n
                                    call slug (ng16, 10, varmax, varupp, varned, &
                                       ansats, org, lock(ng16,10), dubbel, low, &
                                       start(ng16,10), stopp(ng16,10))
                                    do anga = start(ng16,10), stopp(ng16,10), steg(ng16,10)
                                    antel(ng16,10) = anga + antel(ng16,9)
                                    !if (.not.(antel(ng16,10)==antal .and. ansats(&
                                    if (.not.(antel(ng16,10)<=antal .and. ansats(&
                                       ng16,9,1)<=2 .and. ansats(ng16,9,0)<=2)) &
                                       cycle
                                    do plusga = min(anga,22), max(anga - 20,0),-1
                                    ansats(ng16,10,1) = plusga
                                    ansats(ng16,10,0) = anga - plusga
                                    !if (ansats(ng16,10,1)>2 .or. ansats(ng16,10,0)>&
                                    !   2) cycle
! n = 17
                                    nh17 = 17
!     17s
                                    call slug (nh17, 0, varmax, varupp, varned, &
                                       ansats, org, lock(nh17,0), dubbel, low, &
                                       start(nh17,0), stopp(nh17,0))
                                    do anh0 = start(nh17,0), stopp(nh17,0), steg(nh17,0)
                                    antel(nh17,0) = anh0 + antel(nh17-1,10)
                                    if (.not.(antel(nh17,0)<=antal .and. ansats(&
                                       nh17-1,10,1)<=2 .and. ansats(nh17-1,10,0)<=2)) &
                                       cycle
                                    ansats(nh17,0,0) = anh0
!     17p
                                    call slug (nh17, 1, varmax, varupp, varned, &
                                       ansats, org, lock(nh17,1), dubbel, low, &
                                       start(nh17,1), stopp(nh17,1))
                                    do anh1 = start(nh17,1), stopp(nh17,1), steg(nh17,1)
                                    antel(nh17,1) = anh1 + antel(nh17,0)
                                    if (antel(nh17,1) > antal) cycle
                                    do plush1 = min(anh1,4), max(anh1 - 2,0),-1
                                    ansats(nh17,1,1) = plush1
                                    ansats(nh17,1,0) = anh1 - plush1
!     17d
                                    call slug (nh17, 2, varmax, varupp, varned, &
                                       ansats, org, lock(nh17,2), dubbel, low, &
                                       start(nh17,2), stopp(nh17,2))
                                    do anh2 = start(nh17,2), stopp(nh17,2), steg(nh17,2)
                                    antel(nh17,2) = anh2 + antel(nh17,1)
                                    if (antel(nh17,2) > antal) cycle
                                    do plush2 = min(anh2,6), max(anh2 - 4,0),-1
                                    ansats(nh17,2,1) = plush2
                                    ansats(nh17,2,0) = anh2 - plush2
!     17f
                                    call slug (nh17, 3, varmax, varupp, varned, &
                                       ansats, org, lock(nh17,3), dubbel, low, &
                                       start(nh17,3), stopp(nh17,3))
                                    do anh3 = start(nh17,3), stopp(nh17,3), steg(nh17,3)
                                    antel(nh17,3) = anh3 + antel(nh17,2)
                                    if (antel(nh17,3) > antal) cycle
                                    do plush3 = min(anh3,8), max(anh3 - 6,0),-1
                                    ansats(nh17,3,1) = plush3
                                    ansats(nh17,3,0) = anh3 - plush3
!     17g
                                    call slug (nh17, 4, varmax, varupp, varned, &
                                       ansats, org, lock(nh17,4), dubbel, low, &
                                       start(nh17,4), stopp(nh17,4))
                                    do anh4 = start(nh17,4), stopp(nh17,4), steg(nh17,4)
                                    antel(nh17,4) = anh4 + antel(nh17,3)
                                    if (antel(nh17,4) > antal) cycle
                                    do plush4 = min(anh4,10), max(anh4 - 8,0),-1
                                    ansats(nh17,4,1) = plush4
                                    ansats(nh17,4,0) = anh4 - plush4
!     17h
                                    call slug (nh17, 5, varmax, varupp, varned, &
                                       ansats, org, lock(nh17,5), dubbel, low, &
                                       start(nh17,5), stopp(nh17,5))
                                    do anh5 = start(nh17,5), stopp(nh17,5), steg(nh17,5)
                                    antel(nh17,5) = anh5 + antel(nh17,4)
                                    if (antel(nh17,5)>antal .or. ansats(nh17,4,1)>2&
                                       ) cycle
                                    do plush5 = min(anh5,12), max(anh5 - 10,0),-1
                                    ansats(nh17,5,1) = plush5
                                    ansats(nh17,5,0) = anh5 - plush5
!     17i
                                    call slug (nh17, 6, varmax, varupp, varned, &
                                       ansats, org, lock(nh17,6), dubbel, low, &
                                       start(nh17,6), stopp(nh17,6))
                                    do anh6 = start(nh17,6), stopp(nh17,6), steg(nh17,6)
                                    antel(nh17,6) = anh6 + antel(nh17,5)
                                    if (.not.(antel(nh17,6)<=antal .and. ansats(&
                                       nh17,5,1)<=2 .and. ansats(nh17,5,0)<=2)) &
                                       cycle
                                    do plush6 = min(anh6,14), max(anh6 - 12,0),-1
                                    ansats(nh17,6,1) = plush6
                                    ansats(nh17,6,0) = anh6 - plush6
!     17k
                                    call slug (nh17, 7, varmax, varupp, varned, &
                                       ansats, org, lock(nh17,7), dubbel, low, &
                                       start(nh17,7), stopp(nh17,7))
                                    do anh7 = start(nh17,7), stopp(nh17,7), steg(nh17,7)
                                    antel(nh17,7) = anh7 + antel(nh17,6)
                                    if (.not.(antel(nh17,7)<=antal .and. ansats(&
                                       nh17,6,1)<=2 .and. ansats(nh17,6,0)<=2)) &
                                       cycle
                                    do plush7 = min(anh7,16), max(anh7 - 14,0),-1
                                    ansats(nh17,7,1) = plush7
                                    ansats(nh17,7,0) = anh7 - plush7
!     17l
                                    call slug (nh17, 8, varmax, varupp, varned, &
                                       ansats, org, lock(nh17,8), dubbel, low, &
                                       start(nh17,8), stopp(nh17,8))
                                    do anh8 = start(nh17,8), stopp(nh17,8), steg(nh17,8)
                                    antel(nh17,8) = anh8 + antel(nh17,7)
                                    if (.not.(antel(nh17,8)<=antal .and. ansats(&
                                       nh17,7,1)<=2 .and. ansats(nh17,7,0)<=2)) &
                                       cycle
                                    do plush8 = min(anh8,18), max(anh8 - 16,0),-1
                                    ansats(nh17,8,1) = plush8
                                    ansats(nh17,8,0) = anh8 - plush8
!     17m
                                    call slug (nh17, 9, varmax, varupp, varned, &
                                       ansats, org, lock(nh17,9), dubbel, low, &
                                       start(nh17,9), stopp(nh17,9))
                                    do anh9 = start(nh17,9), stopp(nh17,9), steg(nh17,9)
                                    antel(nh17,9) = anh9 + antel(nh17,8)
                                    if (.not.(antel(nh17,9)<=antal .and. ansats(&
                                       nh17,8,1)<=2 .and. ansats(nh17,8,0)<=2)) &
                                       cycle
                                    do plush9 = min(anh9,20), max(anh9 - 18,0),-1
                                    ansats(nh17,9,1) = plush9
                                    ansats(nh17,9,0) = anh9 - plush9
!     17n
                                    call slug (nh17, 10, varmax, varupp, varned, &
                                       ansats, org, lock(nh17,10), dubbel, low, &
                                       start(nh17,10), stopp(nh17,10))
                                    do anha = start(nh17,10), stopp(nh17,10), steg(nh17,10)
                                    antel(nh17,10) = anha + antel(nh17,9)
                                    !if (.not.(antel(nh17,10)==antal .and. ansats(&
                                    if (.not.(antel(nh17,10)<=antal .and. ansats(&
                                       nh17,9,1)<=2 .and. ansats(nh17,9,0)<=2)) &
                                       cycle
                                    do plusha = min(anha,22), max(anha - 20,0),-1
                                    ansats(nh17,10,1) = plusha
                                    ansats(nh17,10,0) = anha - plusha
                                    !if (ansats(nh17,10,1)>2 .or. ansats(nh17,10,0)>&
                                    !   2) cycle
! n = 18
                                    ni18 = 18
!     18s
                                    call slug (ni18, 0, varmax, varupp, varned, &
                                       ansats, org, lock(ni18,0), dubbel, low, &
                                       start(ni18,0), stopp(ni18,0))
                                    do ani0 = start(ni18,0), stopp(ni18,0), steg(ni18,0)
                                    antel(ni18,0) = ani0 + antel(ni18-1,10)
                                    if (.not.(antel(ni18,0)<=antal .and. ansats(&
                                       ni18-1,10,1)<=2 .and. ansats(ni18-1,10,0)<=2)) &
                                       cycle
                                    ansats(ni18,0,0) = ani0
!     18p
                                    call slug (ni18, 1, varmax, varupp, varned, &
                                       ansats, org, lock(ni18,1), dubbel, low, &
                                       start(ni18,1), stopp(ni18,1))
                                    do ani1 = start(ni18,1), stopp(ni18,1), steg(ni18,1)
                                    antel(ni18,1) = ani1 + antel(ni18,0)
                                    if (antel(ni18,1) > antal) cycle
                                    do plusi1 = min(ani1,4), max(ani1 - 2,0),-1
                                    ansats(ni18,1,1) = plusi1
                                    ansats(ni18,1,0) = ani1 - plusi1
!     18d
                                    call slug (ni18, 2, varmax, varupp, varned, &
                                       ansats, org, lock(ni18,2), dubbel, low, &
                                       start(ni18,2), stopp(ni18,2))
                                    do ani2 = start(ni18,2), stopp(ni18,2), steg(ni18,2)
                                    antel(ni18,2) = ani2 + antel(ni18,1)
                                    if (antel(ni18,2) > antal) cycle
                                    do plusi2 = min(ani2,6), max(ani2 - 4,0),-1
                                    ansats(ni18,2,1) = plusi2
                                    ansats(ni18,2,0) = ani2 - plusi2
!     18f
                                    call slug (ni18, 3, varmax, varupp, varned, &
                                       ansats, org, lock(ni18,3), dubbel, low, &
                                       start(ni18,3), stopp(ni18,3))
                                    do ani3 = start(ni18,3), stopp(ni18,3), steg(ni18,3)
                                    antel(ni18,3) = ani3 + antel(ni18,2)
                                    if (antel(ni18,3) > antal) cycle
                                    do plusi3 = min(ani3,8), max(ani3 - 6,0),-1
                                    ansats(ni18,3,1) = plusi3
                                    ansats(ni18,3,0) = ani3 - plusi3
!     18g
                                    call slug (ni18, 4, varmax, varupp, varned, &
                                       ansats, org, lock(ni18,4), dubbel, low, &
                                       start(ni18,4), stopp(ni18,4))
                                    do ani4 = start(ni18,4), stopp(ni18,4), steg(ni18,4)
                                    antel(ni18,4) = ani4 + antel(ni18,3)
                                    if (antel(ni18,4) > antal) cycle
                                    do plusi4 = min(ani4,10), max(ani4 - 8,0),-1
                                    ansats(ni18,4,1) = plusi4
                                    ansats(ni18,4,0) = ani4 - plusi4
!     18h
                                    call slug (ni18, 5, varmax, varupp, varned, &
                                       ansats, org, lock(ni18,5), dubbel, low, &
                                       start(ni18,5), stopp(ni18,5))
                                    do ani5 = start(ni18,5), stopp(ni18,5), steg(ni18,5)
                                    antel(ni18,5) = ani5 + antel(ni18,4)
                                    if (antel(ni18,5)>antal .or. ansats(ni18,4,1)>2&
                                       ) cycle
                                    do plusi5 = min(ani5,12), max(ani5 - 10,0),-1
                                    ansats(ni18,5,1) = plusi5
                                    ansats(ni18,5,0) = ani5 - plusi5
!     18i
                                    call slug (ni18, 6, varmax, varupp, varned, &
                                       ansats, org, lock(ni18,6), dubbel, low, &
                                       start(ni18,6), stopp(ni18,6))
                                    do ani6 = start(ni18,6), stopp(ni18,6), steg(ni18,6)
                                    antel(ni18,6) = ani6 + antel(ni18,5)
                                    if (.not.(antel(ni18,6)<=antal .and. ansats(&
                                       ni18,5,1)<=2 .and. ansats(ni18,5,0)<=2)) &
                                       cycle
                                    do plusi6 = min(ani6,14), max(ani6 - 12,0),-1
                                    ansats(ni18,6,1) = plusi6
                                    ansats(ni18,6,0) = ani6 - plusi6
!     18k
                                    call slug (ni18, 7, varmax, varupp, varned, &
                                       ansats, org, lock(ni18,7), dubbel, low, &
                                       start(ni18,7), stopp(ni18,7))
                                    do ani7 = start(ni18,7), stopp(ni18,7), steg(ni18,7)
                                    antel(ni18,7) = ani7 + antel(ni18,6)
                                    if (.not.(antel(ni18,7)<=antal .and. ansats(&
                                       ni18,6,1)<=2 .and. ansats(ni18,6,0)<=2)) &
                                       cycle
                                    do plusi7 = min(ani7,16), max(ani7 - 14,0),-1
                                    ansats(ni18,7,1) = plusi7
                                    ansats(ni18,7,0) = ani7 - plusi7
!     18l
                                    call slug (ni18, 8, varmax, varupp, varned, &
                                       ansats, org, lock(ni18,8), dubbel, low, &
                                       start(ni18,8), stopp(ni18,8))
                                    do ani8 = start(ni18,8), stopp(ni18,8), steg(ni18,8)
                                    antel(ni18,8) = ani8 + antel(ni18,7)
                                    if (.not.(antel(ni18,8)<=antal .and. ansats(&
                                       ni18,7,1)<=2 .and. ansats(ni18,7,0)<=2)) &
                                       cycle
                                    do plusi8 = min(ani8,18), max(ani8 - 16,0),-1
                                    ansats(ni18,8,1) = plusi8
                                    ansats(ni18,8,0) = ani8 - plusi8
!     18m
                                    call slug (ni18, 9, varmax, varupp, varned, &
                                       ansats, org, lock(ni18,9), dubbel, low, &
                                       start(ni18,9), stopp(ni18,9))
                                    do ani9 = start(ni18,9), stopp(ni18,9), steg(ni18,9)
                                    antel(ni18,9) = ani9 + antel(ni18,8)
                                    if (.not.(antel(ni18,9)<=antal .and. ansats(&
                                       ni18,8,1)<=2 .and. ansats(ni18,8,0)<=2)) &
                                       cycle
                                    do plusi9 = min(ani9,20), max(ani9 - 18,0),-1
                                    ansats(ni18,9,1) = plusi9
                                    ansats(ni18,9,0) = ani9 - plusi9
!     18n
                                    call slug (ni18, 10, varmax, varupp, varned, &
                                       ansats, org, lock(ni18,10), dubbel, low, &
                                       start(ni18,10), stopp(ni18,10))
                                    do ania = start(ni18,10), stopp(ni18,10), steg(ni18,10)
                                    antel(ni18,10) = ania + antel(ni18,9)
                                    !if (.not.(antel(ni18,10)==antal .and. ansats(&
                                    if (.not.(antel(ni18,10)<=antal .and. ansats(&
                                       ni18,9,1)<=2 .and. ansats(ni18,9,0)<=2)) &
                                       cycle
                                    do plusia = min(ania,22), max(ania - 20,0),-1
                                    ansats(ni18,10,1) = plusia
                                    ansats(ni18,10,0) = ania - plusia
                                    !if (ansats(ni18,10,1)>2 .or. ansats(ni18,10,0)>&
                                    !   2) cycle
! n = 19
                                    nj19 = 19
!     19s
                                    call slug (nj19, 0, varmax, varupp, varned, &
                                       ansats, org, lock(nj19,0), dubbel, low, &
                                       start(nj19,0), stopp(nj19,0))
                                    do anj0 = start(nj19,0), stopp(nj19,0), steg(nj19,0)
                                    antel(nj19,0) = anj0 + antel(nj19-1,10)
                                    if (.not.(antel(nj19,0)<=antal .and. ansats(&
                                       nj19-1,10,1)<=2 .and. ansats(nj19-1,10,0)<=2)) &
                                       cycle
                                    ansats(nj19,0,0) = anj0
!     19p
                                    call slug (nj19, 1, varmax, varupp, varned, &
                                       ansats, org, lock(nj19,1), dubbel, low, &
                                       start(nj19,1), stopp(nj19,1))
                                    do anj1 = start(nj19,1), stopp(nj19,1), steg(nj19,1)
                                    antel(nj19,1) = anj1 + antel(nj19,0)
                                    if (antel(nj19,1) > antal) cycle
                                    do plusj1 = min(anj1,4), max(anj1 - 2,0),-1
                                    ansats(nj19,1,1) = plusj1
                                    ansats(nj19,1,0) = anj1 - plusj1
!     19d
                                    call slug (nj19, 2, varmax, varupp, varned, &
                                       ansats, org, lock(nj19,2), dubbel, low, &
                                       start(nj19,2), stopp(nj19,2))
                                    do anj2 = start(nj19,2), stopp(nj19,2), steg(nj19,2)
                                    antel(nj19,2) = anj2 + antel(nj19,1)
                                    if (antel(nj19,2) > antal) cycle
                                    do plusj2 = min(anj2,6), max(anj2 - 4,0),-1
                                    ansats(nj19,2,1) = plusj2
                                    ansats(nj19,2,0) = anj2 - plusj2
!     19f
                                    call slug (nj19, 3, varmax, varupp, varned, &
                                       ansats, org, lock(nj19,3), dubbel, low, &
                                       start(nj19,3), stopp(nj19,3))
                                    do anj3 = start(nj19,3), stopp(nj19,3), steg(nj19,3)
                                    antel(nj19,3) = anj3 + antel(nj19,2)
                                    if (antel(nj19,3) > antal) cycle
                                    do plusj3 = min(anj3,8), max(anj3 - 6,0),-1
                                    ansats(nj19,3,1) = plusj3
                                    ansats(nj19,3,0) = anj3 - plusj3
!     19g
                                    call slug (nj19, 4, varmax, varupp, varned, &
                                       ansats, org, lock(nj19,4), dubbel, low, &
                                       start(nj19,4), stopp(nj19,4))
                                    do anj4 = start(nj19,4), stopp(nj19,4), steg(nj19,4)
                                    antel(nj19,4) = anj4 + antel(nj19,3)
                                    if (antel(nj19,4) > antal) cycle
                                    do plusj4 = min(anj4,10), max(anj4 - 8,0),-1
                                    ansats(nj19,4,1) = plusj4
                                    ansats(nj19,4,0) = anj4 - plusj4
!     19h
                                    call slug (nj19, 5, varmax, varupp, varned, &
                                       ansats, org, lock(nj19,5), dubbel, low, &
                                       start(nj19,5), stopp(nj19,5))
                                    do anj5 = start(nj19,5), stopp(nj19,5), steg(nj19,5)
                                    antel(nj19,5) = anj5 + antel(nj19,4)
                                    if (antel(nj19,5)>antal .or. ansats(nj19,4,1)>2&
                                       ) cycle
                                    do plusj5 = min(anj5,12), max(anj5 - 10,0),-1
                                    ansats(nj19,5,1) = plusj5
                                    ansats(nj19,5,0) = anj5 - plusj5
!     19i
                                    call slug (nj19, 6, varmax, varupp, varned, &
                                       ansats, org, lock(nj19,6), dubbel, low, &
                                       start(nj19,6), stopp(nj19,6))
                                    do anj6 = start(nj19,6), stopp(nj19,6), steg(nj19,6)
                                    antel(nj19,6) = anj6 + antel(nj19,5)
                                    if (.not.(antel(nj19,6)<=antal .and. ansats(&
                                       nj19,5,1)<=2 .and. ansats(nj19,5,0)<=2)) &
                                       cycle
                                    do plusj6 = min(anj6,14), max(anj6 - 12,0),-1
                                    ansats(nj19,6,1) = plusj6
                                    ansats(nj19,6,0) = anj6 - plusj6
!     19k
                                    call slug (nj19, 7, varmax, varupp, varned, &
                                       ansats, org, lock(nj19,7), dubbel, low, &
                                       start(nj19,7), stopp(nj19,7))
                                    do anj7 = start(nj19,7), stopp(nj19,7), steg(nj19,7)
                                    antel(nj19,7) = anj7 + antel(nj19,6)
                                    if (.not.(antel(nj19,7)<=antal .and. ansats(&
                                       nj19,6,1)<=2 .and. ansats(nj19,6,0)<=2)) &
                                       cycle
                                    do plusj7 = min(anj7,16), max(anj7 - 14,0),-1
                                    ansats(nj19,7,1) = plusj7
                                    ansats(nj19,7,0) = anj7 - plusj7
!     19l
                                    call slug (nj19, 8, varmax, varupp, varned, &
                                       ansats, org, lock(nj19,8), dubbel, low, &
                                       start(nj19,8), stopp(nj19,8))
                                    do anj8 = start(nj19,8), stopp(nj19,8), steg(nj19,8)
                                    antel(nj19,8) = anj8 + antel(nj19,7)
                                    if (.not.(antel(nj19,8)<=antal .and. ansats(&
                                       nj19,7,1)<=2 .and. ansats(nj19,7,0)<=2)) &
                                       cycle
                                    do plusj8 = min(anj8,18), max(anj8 - 16,0),-1
                                    ansats(nj19,8,1) = plusj8
                                    ansats(nj19,8,0) = anj8 - plusj8
!     19m
                                    call slug (nj19, 9, varmax, varupp, varned, &
                                       ansats, org, lock(nj19,9), dubbel, low, &
                                       start(nj19,9), stopp(nj19,9))
                                    do anj9 = start(nj19,9), stopp(nj19,9), steg(nj19,9)
                                    antel(nj19,9) = anj9 + antel(nj19,8)
                                    if (.not.(antel(nj19,9)<=antal .and. ansats(&
                                       nj19,8,1)<=2 .and. ansats(nj19,8,0)<=2)) &
                                       cycle
                                    do plusj9 = min(anj9,20), max(anj9 - 18,0),-1
                                    ansats(nj19,9,1) = plusj9
                                    ansats(nj19,9,0) = anj9 - plusj9
!     19n
                                    call slug (nj19, 10, varmax, varupp, varned, &
                                       ansats, org, lock(nj19,10), dubbel, low, &
                                       start(nj19,10), stopp(nj19,10))
                                    do anja = start(nj19,10), stopp(nj19,10), steg(nj19,10)
                                    antel(nj19,10) = anja + antel(nj19,9)
                                    !if (.not.(antel(nj19,10)==antal .and. ansats(&
                                    if (.not.(antel(nj19,10)<=antal .and. ansats(&
                                       nj19,9,1)<=2 .and. ansats(nj19,9,0)<=2)) &
                                       cycle
                                    do plusja = min(anja,22), max(anja - 20,0),-1
                                    ansats(nj19,10,1) = plusja
                                    ansats(nj19,10,0) = anja - plusja
                                    !if (ansats(nj19,10,1)>2 .or. ansats(nj19,10,0)>&
                                    !   2) cycle
! n = 20
                                    nk20 = 20
!     20s
                                    call slug (nk20, 0, varmax, varupp, varned, &
                                       ansats, org, lock(nk20,0), dubbel, low, &
                                       start(nk20,0), stopp(nk20,0))
                                    do ank0 = start(nk20,0), stopp(nk20,0), steg(nk20,0)
                                    antel(nk20,0) = ank0 + antel(nk20-1,10)
                                    if (.not.(antel(nk20,0)<=antal .and. ansats(&
                                       nk20-1,10,1)<=2 .and. ansats(nk20-1,10,0)<=2)) &
                                       cycle
                                    ansats(nk20,0,0) = ank0
!     20p
                                    call slug (nk20, 1, varmax, varupp, varned, &
                                       ansats, org, lock(nk20,1), dubbel, low, &
                                       start(nk20,1), stopp(nk20,1))
                                    do ank1 = start(nk20,1), stopp(nk20,1), steg(nk20,1)
                                    antel(nk20,1) = ank1 + antel(nk20,0)
                                    if (antel(nk20,1) > antal) cycle
                                    do plusk1 = min(ank1,4), max(ank1 - 2,0),-1
                                    ansats(nk20,1,1) = plusk1
                                    ansats(nk20,1,0) = ank1 - plusk1
!     20d
                                    call slug (nk20, 2, varmax, varupp, varned, &
                                       ansats, org, lock(nk20,2), dubbel, low, &
                                       start(nk20,2), stopp(nk20,2))
                                    do ank2 = start(nk20,2), stopp(nk20,2), steg(nk20,2)
                                    antel(nk20,2) = ank2 + antel(nk20,1)
                                    if (antel(nk20,2) > antal) cycle
                                    do plusk2 = min(ank2,6), max(ank2 - 4,0),-1
                                    ansats(nk20,2,1) = plusk2
                                    ansats(nk20,2,0) = ank2 - plusk2
!     20f
                                    call slug (nk20, 3, varmax, varupp, varned, &
                                       ansats, org, lock(nk20,3), dubbel, low, &
                                       start(nk20,3), stopp(nk20,3))
                                    do ank3 = start(nk20,3), stopp(nk20,3), steg(nk20,3)
                                    antel(nk20,3) = ank3 + antel(nk20,2)
                                    if (antel(nk20,3) > antal) cycle
                                    do plusk3 = min(ank3,8), max(ank3 - 6,0),-1
                                    ansats(nk20,3,1) = plusk3
                                    ansats(nk20,3,0) = ank3 - plusk3
!     20g
                                    call slug (nk20, 4, varmax, varupp, varned, &
                                       ansats, org, lock(nk20,4), dubbel, low, &
                                       start(nk20,4), stopp(nk20,4))
                                    do ank4 = start(nk20,4), stopp(nk20,4), steg(nk20,4)
                                    antel(nk20,4) = ank4 + antel(nk20,3)
                                    if (antel(nk20,4) > antal) cycle
                                    do plusk4 = min(ank4,10), max(ank4 - 8,0),-1
                                    ansats(nk20,4,1) = plusk4
                                    ansats(nk20,4,0) = ank4 - plusk4
!     20h
                                    call slug (nk20, 5, varmax, varupp, varned, &
                                       ansats, org, lock(nk20,5), dubbel, low, &
                                       start(nk20,5), stopp(nk20,5))
                                    do ank5 = start(nk20,5), stopp(nk20,5), steg(nk20,5)
                                    antel(nk20,5) = ank5 + antel(nk20,4)
                                    if (antel(nk20,5)>antal .or. ansats(nk20,4,1)>2&
                                       ) cycle
                                    do plusk5 = min(ank5,12), max(ank5 - 10,0),-1
                                    ansats(nk20,5,1) = plusk5
                                    ansats(nk20,5,0) = ank5 - plusk5
!     20i
                                    call slug (nk20, 6, varmax, varupp, varned, &
                                       ansats, org, lock(nk20,6), dubbel, low, &
                                       start(nk20,6), stopp(nk20,6))
                                    do ank6 = start(nk20,6), stopp(nk20,6), steg(nk20,6)
                                    antel(nk20,6) = ank6 + antel(nk20,5)
                                    if (.not.(antel(nk20,6)<=antal .and. ansats(&
                                       nk20,5,1)<=2 .and. ansats(nk20,5,0)<=2)) &
                                       cycle
                                    do plusk6 = min(ank6,14), max(ank6 - 12,0),-1
                                    ansats(nk20,6,1) = plusk6
                                    ansats(nk20,6,0) = ank6 - plusk6
!     20k
                                    call slug (nk20, 7, varmax, varupp, varned, &
                                       ansats, org, lock(nk20,7), dubbel, low, &
                                       start(nk20,7), stopp(nk20,7))
                                    do ank7 = start(nk20,7), stopp(nk20,7), steg(nk20,7)
                                    antel(nk20,7) = ank7 + antel(nk20,6)
                                    if (.not.(antel(nk20,7)<=antal .and. ansats(&
                                       nk20,6,1)<=2 .and. ansats(nk20,6,0)<=2)) &
                                       cycle
                                    do plusk7 = min(ank7,16), max(ank7 - 14,0),-1
                                    ansats(nk20,7,1) = plusk7
                                    ansats(nk20,7,0) = ank7 - plusk7
!     20l
                                    call slug (nk20, 8, varmax, varupp, varned, &
                                       ansats, org, lock(nk20,8), dubbel, low, &
                                       start(nk20,8), stopp(nk20,8))
                                    do ank8 = start(nk20,8), stopp(nk20,8), steg(nk20,8)
                                    antel(nk20,8) = ank8 + antel(nk20,7)
                                    if (.not.(antel(nk20,8)<=antal .and. ansats(&
                                       nk20,7,1)<=2 .and. ansats(nk20,7,0)<=2)) &
                                       cycle
                                    do plusk8 = min(ank8,18), max(ank8 - 16,0),-1
                                    ansats(nk20,8,1) = plusk8
                                    ansats(nk20,8,0) = ank8 - plusk8
!     20m
                                    call slug (nk20, 9, varmax, varupp, varned, &
                                       ansats, org, lock(nk20,9), dubbel, low, &
                                       start(nk20,9), stopp(nk20,9))
                                    do ank9 = start(nk20,9), stopp(nk20,9), steg(nk20,9)
                                    antel(nk20,9) = ank9 + antel(nk20,8)
                                    if (.not.(antel(nk20,9)<=antal .and. ansats(&
                                       nk20,8,1)<=2 .and. ansats(nk20,8,0)<=2)) &
                                       cycle
                                    do plusk9 = min(ank9,20), max(ank9 - 18,0),-1
                                    ansats(nk20,9,1) = plusk9
                                    ansats(nk20,9,0) = ank9 - plusk9
!     20n
                                    call slug (nk20, 10, varmax, varupp, varned, &
                                       ansats, org, lock(nk20,10), dubbel, low, &
                                       start(nk20,10), stopp(nk20,10))
                                    do anka = start(nk20,10), stopp(nk20,10), steg(nk20,10)
                                    antel(nk20,10) = anka + antel(nk20,9)
                                    !if (.not.(antel(nk20,10)==antal .and. ansats(&
                                    if (.not.(antel(nk20,10)<=antal .and. ansats(&
                                       nk20,9,1)<=2 .and. ansats(nk20,9,0)<=2)) &
                                       cycle
                                    do pluska = min(anka,22), max(anka - 20,0),-1
                                    ansats(nk20,10,1) = pluska
                                    ansats(nk20,10,0) = anka - pluska
                                    !if (ansats(nk20,10,1)>2 .or. ansats(nk20,10,0)>&
                                    !   2) cycle
! n = 21
                                    nl21 = 21
!     21s
                                    call slug (nl21, 0, varmax, varupp, varned, &
                                       ansats, org, lock(nl21,0), dubbel, low, &
                                       start(nl21,0), stopp(nl21,0))
                                    do anl0 = start(nl21,0), stopp(nl21,0), steg(nl21,0)
                                    antel(nl21,0) = anl0 + antel(nl21-1,10)
                                    if (.not.(antel(nl21,0)<=antal .and. ansats(&
                                       nl21-1,10,1)<=2 .and. ansats(nl21-1,10,0)<=2)) &
                                       cycle
                                    ansats(nl21,0,0) = anl0
!     21p
                                    call slug (nl21, 1, varmax, varupp, varned, &
                                       ansats, org, lock(nl21,1), dubbel, low, &
                                       start(nl21,1), stopp(nl21,1))
                                    do anl1 = start(nl21,1), stopp(nl21,1), steg(nl21,1)
                                    antel(nl21,1) = anl1 + antel(nl21,0)
                                    if (antel(nl21,1) > antal) cycle
                                    do plusl1 = min(anl1,4), max(anl1 - 2,0),-1
                                    ansats(nl21,1,1) = plusl1
                                    ansats(nl21,1,0) = anl1 - plusl1
!     21d
                                    call slug (nl21, 2, varmax, varupp, varned, &
                                       ansats, org, lock(nl21,2), dubbel, low, &
                                       start(nl21,2), stopp(nl21,2))
                                    do anl2 = start(nl21,2), stopp(nl21,2), steg(nl21,2)
                                    antel(nl21,2) = anl2 + antel(nl21,1)
                                    if (antel(nl21,2) > antal) cycle
                                    do plusl2 = min(anl2,6), max(anl2 - 4,0),-1
                                    ansats(nl21,2,1) = plusl2
                                    ansats(nl21,2,0) = anl2 - plusl2
!     21f
                                    call slug (nl21, 3, varmax, varupp, varned, &
                                       ansats, org, lock(nl21,3), dubbel, low, &
                                       start(nl21,3), stopp(nl21,3))
                                    do anl3 = start(nl21,3), stopp(nl21,3), steg(nl21,3)
                                    antel(nl21,3) = anl3 + antel(nl21,2)
                                    if (antel(nl21,3) > antal) cycle
                                    do plusl3 = min(anl3,8), max(anl3 - 6,0),-1
                                    ansats(nl21,3,1) = plusl3
                                    ansats(nl21,3,0) = anl3 - plusl3
!     21g
                                    call slug (nl21, 4, varmax, varupp, varned, &
                                       ansats, org, lock(nl21,4), dubbel, low, &
                                       start(nl21,4), stopp(nl21,4))
                                    do anl4 = start(nl21,4), stopp(nl21,4), steg(nl21,4)
                                    antel(nl21,4) = anl4 + antel(nl21,3)
                                    if (antel(nl21,4) > antal) cycle
                                    do plusl4 = min(anl4,10), max(anl4 - 8,0),-1
                                    ansats(nl21,4,1) = plusl4
                                    ansats(nl21,4,0) = anl4 - plusl4
!     21h
                                    call slug (nl21, 5, varmax, varupp, varned, &
                                       ansats, org, lock(nl21,5), dubbel, low, &
                                       start(nl21,5), stopp(nl21,5))
                                    do anl5 = start(nl21,5), stopp(nl21,5), steg(nl21,5)
                                    antel(nl21,5) = anl5 + antel(nl21,4)
                                    if (antel(nl21,5)>antal .or. ansats(nl21,4,1)>2&
                                       ) cycle
                                    do plusl5 = min(anl5,12), max(anl5 - 10,0),-1
                                    ansats(nl21,5,1) = plusl5
                                    ansats(nl21,5,0) = anl5 - plusl5
!     21i
                                    call slug (nl21, 6, varmax, varupp, varned, &
                                       ansats, org, lock(nl21,6), dubbel, low, &
                                       start(nl21,6), stopp(nl21,6))
                                    do anl6 = start(nl21,6), stopp(nl21,6), steg(nl21,6)
                                    antel(nl21,6) = anl6 + antel(nl21,5)
                                    if (.not.(antel(nl21,6)<=antal .and. ansats(&
                                       nl21,5,1)<=2 .and. ansats(nl21,5,0)<=2)) &
                                       cycle
                                    do plusl6 = min(anl6,14), max(anl6 - 12,0),-1
                                    ansats(nl21,6,1) = plusl6
                                    ansats(nl21,6,0) = anl6 - plusl6
!     21k
                                    call slug (nl21, 7, varmax, varupp, varned, &
                                       ansats, org, lock(nl21,7), dubbel, low, &
                                       start(nl21,7), stopp(nl21,7))
                                    do anl7 = start(nl21,7), stopp(nl21,7), steg(nl21,7)
                                    antel(nl21,7) = anl7 + antel(nl21,6)
                                    if (.not.(antel(nl21,7)<=antal .and. ansats(&
                                       nl21,6,1)<=2 .and. ansats(nl21,6,0)<=2)) &
                                       cycle
                                    do plusl7 = min(anl7,16), max(anl7 - 14,0),-1
                                    ansats(nl21,7,1) = plusl7
                                    ansats(nl21,7,0) = anl7 - plusl7
!     21l
                                    call slug (nl21, 8, varmax, varupp, varned, &
                                       ansats, org, lock(nl21,8), dubbel, low, &
                                       start(nl21,8), stopp(nl21,8))
                                    do anl8 = start(nl21,8), stopp(nl21,8), steg(nl21,8)
                                    antel(nl21,8) = anl8 + antel(nl21,7)
                                    if (.not.(antel(nl21,8)<=antal .and. ansats(&
                                       nl21,7,1)<=2 .and. ansats(nl21,7,0)<=2)) &
                                       cycle
                                    do plusl8 = min(anl8,18), max(anl8 - 16,0),-1
                                    ansats(nl21,8,1) = plusl8
                                    ansats(nl21,8,0) = anl8 - plusl8
!     21m
                                    call slug (nl21, 9, varmax, varupp, varned, &
                                       ansats, org, lock(nl21,9), dubbel, low, &
                                       start(nl21,9), stopp(nl21,9))
                                    do anl9 = start(nl21,9), stopp(nl21,9), steg(nl21,9)
                                    antel(nl21,9) = anl9 + antel(nl21,8)
                                    if (.not.(antel(nl21,9)<=antal .and. ansats(&
                                       nl21,8,1)<=2 .and. ansats(nl21,8,0)<=2)) &
                                       cycle
                                    do plusl9 = min(anl9,20), max(anl9 - 18,0),-1
                                    ansats(nl21,9,1) = plusl9
                                    ansats(nl21,9,0) = anl9 - plusl9
!     21n
                                    call slug (nl21, 10, varmax, varupp, varned, &
                                       ansats, org, lock(nl21,10), dubbel, low, &
                                       start(nl21,10), stopp(nl21,10))
                                    do anla = start(nl21,10), stopp(nl21,10), steg(nl21,10)
                                    antel(nl21,10) = anla + antel(nl21,9)
                                    !if (.not.(antel(nl21,10)==antal .and. ansats(&
                                    if (.not.(antel(nl21,10)<=antal .and. ansats(&
                                       nl21,9,1)<=2 .and. ansats(nl21,9,0)<=2)) &
                                       cycle
                                    do plusla = min(anla,22), max(anla - 20,0),-1
                                    ansats(nl21,10,1) = plusla
                                    ansats(nl21,10,0) = anla - plusla
                                    if (ansats(nl21,10,1)>2 .or. ansats(nl21,10,0)>&
                                       2) cycle
! n = 22
                                    nm22 = 22
!     22s
                                    call slug (nm22, 0, varmax, varupp, varned, &
                                       ansats, org, lock(nm22,0), dubbel, low, &
                                       start(nm22,0), stopp(nm22,0))
                                    do anm0 = start(nm22,0), stopp(nm22,0), steg(nm22,0)
                                    antel(nm22,0) = anm0 + antel(nm22-1,10)
                                    if (.not.(antel(nm22,0)<=antal .and. ansats(&
                                       nm22-1,10,1)<=2 .and. ansats(nm22-1,10,0)<=2)) &
                                       cycle
                                    ansats(nm22,0,0) = anm0
!     22p
                                    call slug (nm22, 1, varmax, varupp, varned, &
                                       ansats, org, lock(nm22,1), dubbel, low, &
                                       start(nm22,1), stopp(nm22,1))
                                    do anm1 = start(nm22,1), stopp(nm22,1), steg(nm22,1)
                                    antel(nm22,1) = anm1 + antel(nm22,0)
                                    if (antel(nm22,1) > antal) cycle
                                    do plusm1 = min(anm1,4), max(anm1 - 2,0),-1
                                    ansats(nm22,1,1) = plusm1
                                    ansats(nm22,1,0) = anm1 - plusm1
!     22d
                                    call slug (nm22, 2, varmax, varupp, varned, &
                                       ansats, org, lock(nm22,2), dubbel, low, &
                                       start(nm22,2), stopp(nm22,2))
                                    do anm2 = start(nm22,2), stopp(nm22,2), steg(nm22,2)
                                    antel(nm22,2) = anm2 + antel(nm22,1)
                                    if (antel(nm22,2) > antal) cycle
                                    do plusm2 = min(anm2,6), max(anm2 - 4,0),-1
                                    ansats(nm22,2,1) = plusm2
                                    ansats(nm22,2,0) = anm2 - plusm2
!     22f
                                    call slug (nm22, 3, varmax, varupp, varned, &
                                       ansats, org, lock(nm22,3), dubbel, low, &
                                       start(nm22,3), stopp(nm22,3))
                                    do anm3 = start(nm22,3), stopp(nm22,3), steg(nm22,3)
                                    antel(nm22,3) = anm3 + antel(nm22,2)
                                    if (antel(nm22,3) > antal) cycle
                                    do plusm3 = min(anm3,8), max(anm3 - 6,0),-1
                                    ansats(nm22,3,1) = plusm3
                                    ansats(nm22,3,0) = anm3 - plusm3
!     22g
                                    call slug (nm22, 4, varmax, varupp, varned, &
                                       ansats, org, lock(nm22,4), dubbel, low, &
                                       start(nm22,4), stopp(nm22,4))
                                    do anm4 = start(nm22,4), stopp(nm22,4), steg(nm22,4)
                                    antel(nm22,4) = anm4 + antel(nm22,3)
                                    if (antel(nm22,4) > antal) cycle
                                    do plusm4 = min(anm4,10), max(anm4 - 8,0),-1
                                    ansats(nm22,4,1) = plusm4
                                    ansats(nm22,4,0) = anm4 - plusm4
!     22h
                                    call slug (nm22, 5, varmax, varupp, varned, &
                                       ansats, org, lock(nm22,5), dubbel, low, &
                                       start(nm22,5), stopp(nm22,5))
                                    do anm5 = start(nm22,5), stopp(nm22,5), steg(nm22,5)
                                    antel(nm22,5) = anm5 + antel(nm22,4)
                                    if (antel(nm22,5)>antal .or. ansats(nm22,4,1)>2&
                                       ) cycle
                                    do plusm5 = min(anm5,12), max(anm5 - 10,0),-1
                                    ansats(nm22,5,1) = plusm5
                                    ansats(nm22,5,0) = anm5 - plusm5
!     22i
                                    call slug (nm22, 6, varmax, varupp, varned, &
                                       ansats, org, lock(nm22,6), dubbel, low, &
                                       start(nm22,6), stopp(nm22,6))
                                    do anm6 = start(nm22,6), stopp(nm22,6), steg(nm22,6)
                                    antel(nm22,6) = anm6 + antel(nm22,5)
                                    if (.not.(antel(nm22,6)<=antal .and. ansats(&
                                       nm22,5,1)<=2 .and. ansats(nm22,5,0)<=2)) &
                                       cycle
                                    do plusm6 = min(anm6,14), max(anm6 - 12,0),-1
                                    ansats(nm22,6,1) = plusm6
                                    ansats(nm22,6,0) = anm6 - plusm6
!     22k
                                    call slug (nm22, 7, varmax, varupp, varned, &
                                       ansats, org, lock(nm22,7), dubbel, low, &
                                       start(nm22,7), stopp(nm22,7))
                                    do anm7 = start(nm22,7), stopp(nm22,7), steg(nm22,7)
                                    antel(nm22,7) = anm7 + antel(nm22,6)
                                    if (.not.(antel(nm22,7)<=antal .and. ansats(&
                                       nm22,6,1)<=2 .and. ansats(nm22,6,0)<=2)) &
                                       cycle
                                    do plusm7 = min(anm7,16), max(anm7 - 14,0),-1
                                    ansats(nm22,7,1) = plusm7
                                    ansats(nm22,7,0) = anm7 - plusm7
!     22l
                                    call slug (nm22, 8, varmax, varupp, varned, &
                                       ansats, org, lock(nm22,8), dubbel, low, &
                                       start(nm22,8), stopp(nm22,8))
                                    do anm8 = start(nm22,8), stopp(nm22,8), steg(nm22,8)
                                    antel(nm22,8) = anm8 + antel(nm22,7)
                                    if (.not.(antel(nm22,8)<=antal .and. ansats(&
                                       nm22,7,1)<=2 .and. ansats(nm22,7,0)<=2)) &
                                       cycle
                                    do plusm8 = min(anm8,18), max(anm8 - 16,0),-1
                                    ansats(nm22,8,1) = plusm8
                                    ansats(nm22,8,0) = anm8 - plusm8
!     22m
                                    call slug (nm22, 9, varmax, varupp, varned, &
                                       ansats, org, lock(nm22,9), dubbel, low, &
                                       start(nm22,9), stopp(nm22,9))
                                    do anm9 = start(nm22,9), stopp(nm22,9), steg(nm22,9)
                                    antel(nm22,9) = anm9 + antel(nm22,8)
                                    if (.not.(antel(nm22,9)<=antal .and. ansats(&
                                       nm22,8,1)<=2 .and. ansats(nm22,8,0)<=2)) &
                                       cycle
                                    do plusm9 = min(anm9,20), max(anm9 - 18,0),-1
                                    ansats(nm22,9,1) = plusm9
                                    ansats(nm22,9,0) = anm9 - plusm9
!     22n
                                    call slug (nm22, 10, varmax, varupp, varned, &
                                       ansats, org, lock(nm22,10), dubbel, low, &
                                       start(nm22,10), stopp(nm22,10))
                                    do anma = start(nm22,10), stopp(nm22,10), steg(nm22,10)
                                    antel(nm22,10) = anma + antel(nm22,9)
                                    !if (.not.(antel(nm22,10)==antal .and. ansats(&
                                    if (.not.(antel(nm22,10)<=antal .and. ansats(&
                                       nm22,9,1)<=2 .and. ansats(nm22,9,0)<=2)) &
                                       cycle
                                    do plusma = min(anma,22), max(anma - 20,0),-1
                                    ansats(nm22,10,1) = plusma
                                    ansats(nm22,10,0) = anma - plusma
                                    if (ansats(nm22,10,1)>2 .or. ansats(nm22,10,0)>&
                                       2) cycle
! n = 23
                                    nn23 = 23
!     23s
                                    call slug (nn23, 0, varmax, varupp, varned, &
                                       ansats, org, lock(nn23,0), dubbel, low, &
                                       start(nn23,0), stopp(nn23,0))
                                    do ann0 = start(nn23,0), stopp(nn23,0), steg(nn23,0)
                                    antel(nn23,0) = ann0 + antel(nn23-1,10)
                                    if (.not.(antel(nn23,0)<=antal .and. ansats(&
                                       nn23-1,10,1)<=2 .and. ansats(nn23-1,10,0)<=2)) &
                                       cycle
                                    ansats(nn23,0,0) = ann0
!     23p
                                    call slug (nn23, 1, varmax, varupp, varned, &
                                       ansats, org, lock(nn23,1), dubbel, low, &
                                       start(nn23,1), stopp(nn23,1))
                                    do ann1 = start(nn23,1), stopp(nn23,1), steg(nn23,1)
                                    antel(nn23,1) = ann1 + antel(nn23,0)
                                    if (antel(nn23,1) > antal) cycle
                                    do plusn1 = min(ann1,4), max(ann1 - 2,0),-1
                                    ansats(nn23,1,1) = plusn1
                                    ansats(nn23,1,0) = ann1 - plusn1
!     23d
                                    call slug (nn23, 2, varmax, varupp, varned, &
                                       ansats, org, lock(nn23,2), dubbel, low, &
                                       start(nn23,2), stopp(nn23,2))
                                    do ann2 = start(nn23,2), stopp(nn23,2), steg(nn23,2)
                                    antel(nn23,2) = ann2 + antel(nn23,1)
                                    if (antel(nn23,2) > antal) cycle
                                    do plusn2 = min(ann2,6), max(ann2 - 4,0),-1
                                    ansats(nn23,2,1) = plusn2
                                    ansats(nn23,2,0) = ann2 - plusn2
!     23f
                                    call slug (nn23, 3, varmax, varupp, varned, &
                                       ansats, org, lock(nn23,3), dubbel, low, &
                                       start(nn23,3), stopp(nn23,3))
                                    do ann3 = start(nn23,3), stopp(nn23,3), steg(nn23,3)
                                    antel(nn23,3) = ann3 + antel(nn23,2)
                                    if (antel(nn23,3) > antal) cycle
                                    do plusn3 = min(ann3,8), max(ann3 - 6,0),-1
                                    ansats(nn23,3,1) = plusn3
                                    ansats(nn23,3,0) = ann3 - plusn3
!     23g
                                    call slug (nn23, 4, varmax, varupp, varned, &
                                       ansats, org, lock(nn23,4), dubbel, low, &
                                       start(nn23,4), stopp(nn23,4))
                                    do ann4 = start(nn23,4), stopp(nn23,4), steg(nn23,4)
                                    antel(nn23,4) = ann4 + antel(nn23,3)
                                    if (antel(nn23,4) > antal) cycle
                                    do plusn4 = min(ann4,10), max(ann4 - 8,0),-1
                                    ansats(nn23,4,1) = plusn4
                                    ansats(nn23,4,0) = ann4 - plusn4
!     23h
                                    call slug (nn23, 5, varmax, varupp, varned, &
                                       ansats, org, lock(nn23,5), dubbel, low, &
                                       start(nn23,5), stopp(nn23,5))
                                    do ann5 = start(nn23,5), stopp(nn23,5), steg(nn23,5)
                                    antel(nn23,5) = ann5 + antel(nn23,4)
                                    if (antel(nn23,5)>antal .or. ansats(nn23,4,1)>2&
                                       ) cycle
                                    do plusn5 = min(ann5,12), max(ann5 - 10,0),-1
                                    ansats(nn23,5,1) = plusn5
                                    ansats(nn23,5,0) = ann5 - plusn5
!     23i
                                    call slug (nn23, 6, varmax, varupp, varned, &
                                       ansats, org, lock(nn23,6), dubbel, low, &
                                       start(nn23,6), stopp(nn23,6))
                                    do ann6 = start(nn23,6), stopp(nn23,6), steg(nn23,6)
                                    antel(nn23,6) = ann6 + antel(nn23,5)
                                    if (.not.(antel(nn23,6)<=antal .and. ansats(&
                                       nn23,5,1)<=2 .and. ansats(nn23,5,0)<=2)) &
                                       cycle
                                    do plusn6 = min(ann6,14), max(ann6 - 12,0),-1
                                    ansats(nn23,6,1) = plusn6
                                    ansats(nn23,6,0) = ann6 - plusn6
!     23k
                                    call slug (nn23, 7, varmax, varupp, varned, &
                                       ansats, org, lock(nn23,7), dubbel, low, &
                                       start(nn23,7), stopp(nn23,7))
                                    do ann7 = start(nn23,7), stopp(nn23,7), steg(nn23,7)
                                    antel(nn23,7) = ann7 + antel(nn23,6)
                                    if (.not.(antel(nn23,7)<=antal .and. ansats(&
                                       nn23,6,1)<=2 .and. ansats(nn23,6,0)<=2)) &
                                       cycle
                                    do plusn7 = min(ann7,16), max(ann7 - 14,0),-1
                                    ansats(nn23,7,1) = plusn7
                                    ansats(nn23,7,0) = ann7 - plusn7
!     23l
                                    call slug (nn23, 8, varmax, varupp, varned, &
                                       ansats, org, lock(nn23,8), dubbel, low, &
                                       start(nn23,8), stopp(nn23,8))
                                    do ann8 = start(nn23,8), stopp(nn23,8), steg(nn23,8)
                                    antel(nn23,8) = ann8 + antel(nn23,7)
                                    if (.not.(antel(nn23,8)<=antal .and. ansats(&
                                       nn23,7,1)<=2 .and. ansats(nn23,7,0)<=2)) &
                                       cycle
                                    do plusn8 = min(ann8,18), max(ann8 - 16,0),-1
                                    ansats(nn23,8,1) = plusn8
                                    ansats(nn23,8,0) = ann8 - plusn8
!     23m
                                    call slug (nn23, 9, varmax, varupp, varned, &
                                       ansats, org, lock(nn23,9), dubbel, low, &
                                       start(nn23,9), stopp(nn23,9))
                                    do ann9 = start(nn23,9), stopp(nn23,9), steg(nn23,9)
                                    antel(nn23,9) = ann9 + antel(nn23,8)
                                    if (.not.(antel(nn23,9)<=antal .and. ansats(&
                                       nn23,8,1)<=2 .and. ansats(nn23,8,0)<=2)) &
                                       cycle
                                    do plusn9 = min(ann9,20), max(ann9 - 18,0),-1
                                    ansats(nn23,9,1) = plusn9
                                    ansats(nn23,9,0) = ann9 - plusn9
!     23n
                                    call slug (nn23, 10, varmax, varupp, varned, &
                                       ansats, org, lock(nn23,10), dubbel, low, &
                                       start(nn23,10), stopp(nn23,10))
                                    do anna = start(nn23,10), stopp(nn23,10), steg(nn23,10)
                                    antel(nn23,10) = anna + antel(nn23,9)
                                    !if (.not.(antel(nn23,10)==antal .and. ansats(&
                                    if (.not.(antel(nn23,10)<=antal .and. ansats(&
                                       nn23,9,1)<=2 .and. ansats(nn23,9,0)<=2)) &
                                       cycle
                                    do plusna = min(anna,22), max(anna - 20,0),-1
                                    ansats(nn23,10,1) = plusna
                                    ansats(nn23,10,0) = anna - plusna
                                    if (ansats(nn23,10,1)>2 .or. ansats(nn23,10,0)>&
                                       2) cycle
! n = 24
                                    no24 = 24
!     24s
                                    call slug (no24, 0, varmax, varupp, varned, &
                                       ansats, org, lock(no24,0), dubbel, low, &
                                       start(no24,0), stopp(no24,0))
                                    do ano0 = start(no24,0), stopp(no24,0), steg(no24,0)
                                    antel(no24,0) = ano0 + antel(no24-1,10)
                                    if (.not.(antel(no24,0)<=antal .and. ansats(&
                                       no24-1,10,1)<=2 .and. ansats(no24-1,10,0)<=2)) &
                                       cycle
                                    ansats(no24,0,0) = ano0
!     24p
                                    call slug (no24, 1, varmax, varupp, varned, &
                                       ansats, org, lock(no24,1), dubbel, low, &
                                       start(no24,1), stopp(no24,1))
                                    do ano1 = start(no24,1), stopp(no24,1), steg(no24,1)
                                    antel(no24,1) = ano1 + antel(no24,0)
                                    if (antel(no24,1) > antal) cycle
                                    do pluso1 = min(ano1,4), max(ano1 - 2,0),-1
                                    ansats(no24,1,1) = pluso1
                                    ansats(no24,1,0) = ano1 - pluso1
!     24d
                                    call slug (no24, 2, varmax, varupp, varned, &
                                       ansats, org, lock(no24,2), dubbel, low, &
                                       start(no24,2), stopp(no24,2))
                                    do ano2 = start(no24,2), stopp(no24,2), steg(no24,2)
                                    antel(no24,2) = ano2 + antel(no24,1)
                                    if (antel(no24,2) > antal) cycle
                                    do pluso2 = min(ano2,6), max(ano2 - 4,0),-1
                                    ansats(no24,2,1) = pluso2
                                    ansats(no24,2,0) = ano2 - pluso2
!     24f
                                    call slug (no24, 3, varmax, varupp, varned, &
                                       ansats, org, lock(no24,3), dubbel, low, &
                                       start(no24,3), stopp(no24,3))
                                    do ano3 = start(no24,3), stopp(no24,3), steg(no24,3)
                                    antel(no24,3) = ano3 + antel(no24,2)
                                    if (antel(no24,3) > antal) cycle
                                    do pluso3 = min(ano3,8), max(ano3 - 6,0),-1
                                    ansats(no24,3,1) = pluso3
                                    ansats(no24,3,0) = ano3 - pluso3
!     24g
                                    call slug (no24, 4, varmax, varupp, varned, &
                                       ansats, org, lock(no24,4), dubbel, low, &
                                       start(no24,4), stopp(no24,4))
                                    do ano4 = start(no24,4), stopp(no24,4), steg(no24,4)
                                    antel(no24,4) = ano4 + antel(no24,3)
                                    if (antel(no24,4) > antal) cycle
                                    do pluso4 = min(ano4,10), max(ano4 - 8,0),-1
                                    ansats(no24,4,1) = pluso4
                                    ansats(no24,4,0) = ano4 - pluso4
!     24h
                                    call slug (no24, 5, varmax, varupp, varned, &
                                       ansats, org, lock(no24,5), dubbel, low, &
                                       start(no24,5), stopp(no24,5))
                                    do ano5 = start(no24,5), stopp(no24,5), steg(no24,5)
                                    antel(no24,5) = ano5 + antel(no24,4)
                                    if (antel(no24,5)>antal .or. ansats(no24,4,1)>2&
                                       ) cycle
                                    do pluso5 = min(ano5,12), max(ano5 - 10,0),-1
                                    ansats(no24,5,1) = pluso5
                                    ansats(no24,5,0) = ano5 - pluso5
!     24i
                                    call slug (no24, 6, varmax, varupp, varned, &
                                       ansats, org, lock(no24,6), dubbel, low, &
                                       start(no24,6), stopp(no24,6))
                                    do ano6 = start(no24,6), stopp(no24,6), steg(no24,6)
                                    antel(no24,6) = ano6 + antel(no24,5)
                                    if (.not.(antel(no24,6)<=antal .and. ansats(&
                                       no24,5,1)<=2 .and. ansats(no24,5,0)<=2)) &
                                       cycle
                                    do pluso6 = min(ano6,14), max(ano6 - 12,0),-1
                                    ansats(no24,6,1) = pluso6
                                    ansats(no24,6,0) = ano6 - pluso6
!     24k
                                    call slug (no24, 7, varmax, varupp, varned, &
                                       ansats, org, lock(no24,7), dubbel, low, &
                                       start(no24,7), stopp(no24,7))
                                    do ano7 = start(no24,7), stopp(no24,7), steg(no24,7)
                                    antel(no24,7) = ano7 + antel(no24,6)
                                    if (.not.(antel(no24,7)<=antal .and. ansats(&
                                       no24,6,1)<=2 .and. ansats(no24,6,0)<=2)) &
                                       cycle
                                    do pluso7 = min(ano7,16), max(ano7 - 14,0),-1
                                    ansats(no24,7,1) = pluso7
                                    ansats(no24,7,0) = ano7 - pluso7
!     24l
                                    call slug (no24, 8, varmax, varupp, varned, &
                                       ansats, org, lock(no24,8), dubbel, low, &
                                       start(no24,8), stopp(no24,8))
                                    do ano8 = start(no24,8), stopp(no24,8), steg(no24,8)
                                    antel(no24,8) = ano8 + antel(no24,7)
                                    if (.not.(antel(no24,8)<=antal .and. ansats(&
                                       no24,7,1)<=2 .and. ansats(no24,7,0)<=2)) &
                                       cycle
                                    do pluso8 = min(ano8,18), max(ano8 - 16,0),-1
                                    ansats(no24,8,1) = pluso8
                                    ansats(no24,8,0) = ano8 - pluso8
!     24m
                                    call slug (no24, 9, varmax, varupp, varned, &
                                       ansats, org, lock(no24,9), dubbel, low, &
                                       start(no24,9), stopp(no24,9))
                                    do ano9 = start(no24,9), stopp(no24,9), steg(no24,9)
                                    antel(no24,9) = ano9 + antel(no24,8)
                                    if (.not.(antel(no24,9)<=antal .and. ansats(&
                                       no24,8,1)<=2 .and. ansats(no24,8,0)<=2)) &
                                       cycle
                                    do pluso9 = min(ano9,20), max(ano9 - 18,0),-1
                                    ansats(no24,9,1) = pluso9
                                    ansats(no24,9,0) = ano9 - pluso9
!     24n
                                    call slug (no24, 10, varmax, varupp, varned, &
                                       ansats, org, lock(no24,10), dubbel, low, &
                                       start(no24,10), stopp(no24,10))
                                    do anoa = start(no24,10), stopp(no24,10), steg(no24,10)
                                    antel(no24,10) = anoa + antel(no24,9)
                                    !if (.not.(antel(no24,10)==antal .and. ansats(&
                                    if (.not.(antel(no24,10)<=antal .and. ansats(&
                                       no24,9,1)<=2 .and. ansats(no24,9,0)<=2)) &
                                       cycle
                                    do plusoa = min(anoa,22), max(anoa - 20,0),-1
                                    ansats(no24,10,1) = plusoa
                                    ansats(no24,10,0) = anoa - plusoa
                                    if (ansats(no24,10,1)>2 .or. ansats(no24,10,0)>&
                                       2) cycle
! n = 25
                                    np25 = 25
!     25s
                                    call slug (np25, 0, varmax, varupp, varned, &
                                       ansats, org, lock(np25,0), dubbel, low, &
                                       start(np25,0), stopp(np25,0))
                                    do anp0 = start(np25,0), stopp(np25,0), steg(np25,0)
                                    antel(np25,0) = anp0 + antel(np25-1,10)
                                    if (.not.(antel(np25,0)<=antal .and. ansats(&
                                       np25-1,10,1)<=2 .and. ansats(np25-1,10,0)<=2)) &
                                       cycle
                                    ansats(np25,0,0) = anp0
!     25p
                                    call slug (np25, 1, varmax, varupp, varned, &
                                       ansats, org, lock(np25,1), dubbel, low, &
                                       start(np25,1), stopp(np25,1))
                                    do anp1 = start(np25,1), stopp(np25,1), steg(np25,1)
                                    antel(np25,1) = anp1 + antel(np25,0)
                                    if (antel(np25,1) > antal) cycle
                                    do plusp1 = min(anp1,4), max(anp1 - 2,0),-1
                                    ansats(np25,1,1) = plusp1
                                    ansats(np25,1,0) = anp1 - plusp1
!     25d
                                    call slug (np25, 2, varmax, varupp, varned, &
                                       ansats, org, lock(np25,2), dubbel, low, &
                                       start(np25,2), stopp(np25,2))
                                    do anp2 = start(np25,2), stopp(np25,2), steg(np25,2)
                                    antel(np25,2) = anp2 + antel(np25,1)
                                    if (antel(np25,2) > antal) cycle
                                    do plusp2 = min(anp2,6), max(anp2 - 4,0),-1
                                    ansats(np25,2,1) = plusp2
                                    ansats(np25,2,0) = anp2 - plusp2
!     25f
                                    call slug (np25, 3, varmax, varupp, varned, &
                                       ansats, org, lock(np25,3), dubbel, low, &
                                       start(np25,3), stopp(np25,3))
                                    do anp3 = start(np25,3), stopp(np25,3), steg(np25,3)
                                    antel(np25,3) = anp3 + antel(np25,2)
                                    if (antel(np25,3) > antal) cycle
                                    do plusp3 = min(anp3,8), max(anp3 - 6,0),-1
                                    ansats(np25,3,1) = plusp3
                                    ansats(np25,3,0) = anp3 - plusp3
!     25g
                                    call slug (np25, 4, varmax, varupp, varned, &
                                       ansats, org, lock(np25,4), dubbel, low, &
                                       start(np25,4), stopp(np25,4))
                                    do anp4 = start(np25,4), stopp(np25,4), steg(np25,4)
                                    antel(np25,4) = anp4 + antel(np25,3)
                                    if (antel(np25,4) > antal) cycle
                                    do plusp4 = min(anp4,10), max(anp4 - 8,0),-1
                                    ansats(np25,4,1) = plusp4
                                    ansats(np25,4,0) = anp4 - plusp4
!     25h
                                    call slug (np25, 5, varmax, varupp, varned, &
                                       ansats, org, lock(np25,5), dubbel, low, &
                                       start(np25,5), stopp(np25,5))
                                    do anp5 = start(np25,5), stopp(np25,5), steg(np25,5)
                                    antel(np25,5) = anp5 + antel(np25,4)
                                    if (antel(np25,5)>antal .or. ansats(np25,4,1)>2&
                                       ) cycle
                                    do plusp5 = min(anp5,12), max(anp5 - 10,0),-1
                                    ansats(np25,5,1) = plusp5
                                    ansats(np25,5,0) = anp5 - plusp5
!     25i
                                    call slug (np25, 6, varmax, varupp, varned, &
                                       ansats, org, lock(np25,6), dubbel, low, &
                                       start(np25,6), stopp(np25,6))
                                    do anp6 = start(np25,6), stopp(np25,6), steg(np25,6)
                                    antel(np25,6) = anp6 + antel(np25,5)
                                    if (.not.(antel(np25,6)<=antal .and. ansats(&
                                       np25,5,1)<=2 .and. ansats(np25,5,0)<=2)) &
                                       cycle
                                    do plusp6 = min(anp6,14), max(anp6 - 12,0),-1
                                    ansats(np25,6,1) = plusp6
                                    ansats(np25,6,0) = anp6 - plusp6
!     25k
                                    call slug (np25, 7, varmax, varupp, varned, &
                                       ansats, org, lock(np25,7), dubbel, low, &
                                       start(np25,7), stopp(np25,7))
                                    do anp7 = start(np25,7), stopp(np25,7), steg(np25,7)
                                    antel(np25,7) = anp7 + antel(np25,6)
                                    if (.not.(antel(np25,7)<=antal .and. ansats(&
                                       np25,6,1)<=2 .and. ansats(np25,6,0)<=2)) &
                                       cycle
                                    do plusp7 = min(anp7,16), max(anp7 - 14,0),-1
                                    ansats(np25,7,1) = plusp7
                                    ansats(np25,7,0) = anp7 - plusp7
!     25l
                                    call slug (np25, 8, varmax, varupp, varned, &
                                       ansats, org, lock(np25,8), dubbel, low, &
                                       start(np25,8), stopp(np25,8))
                                    do anp8 = start(np25,8), stopp(np25,8), steg(np25,8)
                                    antel(np25,8) = anp8 + antel(np25,7)
                                    if (.not.(antel(np25,8)<=antal .and. ansats(&
                                       np25,7,1)<=2 .and. ansats(np25,7,0)<=2)) &
                                       cycle
                                    do plusp8 = min(anp8,18), max(anp8 - 16,0),-1
                                    ansats(np25,8,1) = plusp8
                                    ansats(np25,8,0) = anp8 - plusp8
!     25m
                                    call slug (np25, 9, varmax, varupp, varned, &
                                       ansats, org, lock(np25,9), dubbel, low, &
                                       start(np25,9), stopp(np25,9))
                                    do anp9 = start(np25,9), stopp(np25,9), steg(np25,9)
                                    antel(np25,9) = anp9 + antel(np25,8)
                                    if (.not.(antel(np25,9)<=antal .and. ansats(&
                                       np25,8,1)<=2 .and. ansats(np25,8,0)<=2)) &
                                       cycle
                                    do plusp9 = min(anp9,20), max(anp9 - 18,0),-1
                                    ansats(np25,9,1) = plusp9
                                    ansats(np25,9,0) = anp9 - plusp9
!     25n
                                    call slug (np25, 10, varmax, varupp, varned, &
                                       ansats, org, lock(np25,10), dubbel, low, &
                                       start(np25,10), stopp(np25,10))
                                    do anpa = start(np25,10), stopp(np25,10), steg(np25,10)
                                    antel(np25,10) = anpa + antel(np25,9)
                                    if (.not.(antel(np25,10)==antal .and. ansats(&
                                    !if (.not.(antel(np25,10)<=antal .and. ansats(&
                                       np25,9,1)<=2 .and. ansats(np25,9,0)<=2)) &
                                       cycle
                                    do pluspa = min(anpa,22), max(anpa - 20,0),-1
                                    ansats(np25,10,1) = pluspa
                                    ansats(np25,10,0) = anpa - pluspa
                                    if (ansats(np25,10,1)>2 .or. ansats(np25,10,0)>&
                                       2) cycle

                                    par = 0
                                    elar = 0
                                    do i = 1, 25
                                    do j = 0, min(10,i - 1)
                                    do k = 0, min(j,1)
                                    elar = elar + ansats(i,j,k)
                                    par = mod(par + j*ansats(i,j,k),2)
                                    end do
                                    end do
                                    end do
                                    if (par /= par0) cycle

                                    if (elar == antal) then
                                    call gen (ansats, posn, posl, skal, cf, &
                                       first, minj, maxj, par0)
                                    else
                                    write (*, *) 'FEL'
                                    endif
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      if (first) then
         rewind (fil_1)
      else
         rewind (fil_2)
      endif
      if (cf == 0) then
         write (*, 1005) 'No configuration state has been generated.'
      else if (cf == 1) then
         write (*, 1005) 'One configuration state has been generated.'
      else if (cf < 10) then
         write (*, 1001) cf, ' configuration states have been generated.'
      else if (cf < 100) then
         write (*, 1002) cf, ' configuration states have been generated.'
      else if (cf < 1000) then
         write (*, 1003) cf, ' configuration states have been generated.'
      else if (cf < 10000) then
         write (*, 1004) cf, ' configuration states have been generated.'
      else if (cf < 100000) then
         write (*, 1006) cf, ' configuration states have been generated.'
      else
         write (*, *) cf, ' configuration states have been generated.'
      endif
! 1000 format(A)
 1001 format(' ',i1,a)
 1002 format(' ',i2,a)
 1003 format(' ',i3,a)
 1004 format(' ',i4,a)
 1005 format(' ',a)
 1006 format(' ',i5,a)
 5000 format(11i2)
      return
      end subroutine blanda
