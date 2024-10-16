!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     RCSFSYM
!
!     This subroutine removes un-needed labelling CSFs
!      
!     Let n = 7 be the max principal quantum number for the
!     symmetry-ordered orbitals then 1s(2)6s7s should be kept 
!     wheras 1s(2)6s7d, 1s(2)7s6d 
!     should be removed and only 1s(2)7s7d kept      
!     
!     WRitten by Per Jonsson, Malmo University, Sweden
!     May 2017
! 
!     Modified by Chongyang Chen, Fudan University, China
!     Feb 2021
!
!     Last Modification by CYC,  Dec. 2023
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!...Translated by  Gediminas Gaigalas  04/21/2021


      subroutine rcsfsym
      implicit none
      character*1500 :: string(5),stringappend,line4
      character*1500 :: stringnonsym  ! cyc
      character*300 :: line1,line2,line3
      character*2 :: string1,string2,string3,string4,nmax(21),sym(21) 
      character*2 :: nmaxnon(21),nstring
      integer :: i,j,nlength,nonsym,n,keep,nc,m,iunit
      integer :: nlengthcore           ! cyc
      integer :: nmaxi(21),nmaxnoni(21),n0(21)

!cyc==================================================================== 
      character*4        :: char4a, char4b
      integer            :: nmax_spec, lmax_spec,nmaxcan, nonsym_free
      integer            :: orb_spec_nl(0:25),n1,n2,ndiscard0
      integer            :: nmax_orbset(100),nmrcon, ntime
      logical            :: freeorb
      common/orb_spec/nmax_spec,lmax_spec,nmaxcan,nonsym_free,       &
            orb_spec_nl,nmax_orbset,nmrcon,freeorb,ntime
!cyc==================================================================== 

!      write(*,*) 'nonsym_free=',nonsym_free
      write(*,*)
      write(*,*) 'Write the leading five lines'
      write(*,*)

      stringappend = repeat(' ',1500)  ! Initialize stringappend
      nmax = '  '
      nmaxnon = '  '

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

      n0(1)  = 1
      n0(2)  = 2
      n0(3)  = 2
      n0(4)  = 3
      n0(5)  = 3
      n0(6)  = 4
      n0(7)  = 4
      n0(8)  = 5
      n0(9)  = 5
      n0(10) = 6
      n0(11) = 6
      n0(12) = 7
      n0(13) = 7
      n0(14) = 8
      n0(15) = 8
      n0(16) = 9
      n0(17) = 9
      n0(18) = 10
      n0(19) = 10
      n0(20) = 11
      n0(21) = 11

      open (19, file = 'rcsfg.out', status = 'old', form = 'formatted',&
                action='readwrite')
      open (1312, status = 'scratch', form = 'formatted')
      open (1313, status = 'scratch', form = 'formatted')
      open (2425, file = 'rlabel.out', status = 'unknown',        &
                  form = 'formatted', action = 'write')

! Read and analyze header to determine the highest principal quantum
! number for each symmetry. Also add orbitals to the header

      do i = 1, 5
         read (19,'(a)') string(i)
      end do

      nlength = len_trim(string(4)) + 1
! cyc, nlengthcore is used to determine the number of core-orbital. 
      nlengthcore = len_trim(string(2)) + 1
      if (nlengthcore.lt.5) nlengthcore = 0
! cyc

!cyc: Flexible labelling orbital set, using nonsym_free 
      nonsym = nonsym_free
!      write(*,*)'nonsym=',nonsym
      write(2425,*)nonsym, '  = nonsym'

! Analyze the labelling orbitals
! cyc =================================================================
! Some labelling orbitals are stored in string(2)
      stringnonsym = string(4)(1:nonsym*5)
      if (nlengthcore.gt.0)                                           &
       stringnonsym = trim(string(2)) // ' ' //                      &
                      string(4)(1:(nonsym-nlengthcore/5)*5)
      write(*,*)'Labeling orbitals:'
      write(*,'(a)')trim(stringnonsym)
      write(2425,'(a)')trim(stringnonsym)
      do i = 1,nonsym
         string1 = stringnonsym(5*(i-1)+4:5*i)   ! symmmetry string
         string2 = stringnonsym(5*(i-1)+2:5*i-2) ! n string
         do j = 1,21
            if (string1.eq.sym(j)) then
               nmaxnon(j) = string2
            end if
         end do
      end do
      write(*,*)'Highest Labelling orbitals:'
      do i = 1,21
        if (nmaxnon(i).ne.'  ') then
          write(*,*) nmaxnon(i),sym(i)
        end if
      end do
! cyc =================================================================

! Analyze the orbtitals and determine the maximum n for each symmetry
      do i = nonsym-nlengthcore/5+1,nlength/5
         string1 = string(4)(5*(i-1)+4:5*i)   ! symmmetry string
         string2 = string(4)(5*(i-1)+2:5*i-2) ! n string
         do j = 1,21
            if (string1.eq.sym(j)) then
               nmax(j) = string2
            end if
         end do
      end do
      write(*,*)'Highest symmetry-ordered orbitals:'
      do i = 1,21
       if (nmax(i).ne.'  ') then
         write(*,*) nmax(i),sym(i)
       end if
      end do

! Obtain line4
! Do some manipulations to write out all orbitals on line 4
      do i = 1,21
         read(nmax(i),'(i2)') nmaxi(i) ! Convert nmax string to integer
         read(nmaxnon(i),'(i2)') nmaxnoni(i)
      end do

      do i = 1,21
        if (nmaxi(i).gt.0) then
          do j = max0(n0(i),nmaxnoni(i)+1),nmaxi(i)
            write(nstring,'(i2)') j  ! Convert integer j to string
            if (sym(i)(2:2).ne.'-') sym(i)(2:2) = '+'
            stringappend = trim(stringappend)//' '//nstring//sym(i)
          end do
        end if
      end do
      line4 = trim(stringappend)
      !write(*,'(A)')trim(line4)
      !write(*,*)"nmrcon=",nmrcon

! Reading the THREE lines for each CSF, Counting the discarded CSFs
      ndiscard0 = 0
   10 continue
      read (19, '(a)', end = 99) line1
      read (19, '(a)') line2
      read (19, '(a)') line3

! Count the number of symmetry-ordered orbitals in the CSF 
      n = 0
      do i = nonsym-nlengthcore/5+1,nlength/5
         if (index(line1,string(4)(5*(i-1)+1:5*i)).gt.0) n = n + 1
      end do

! All CSFs kept, processed later in rcsfsym_cyc.f
      keep = 1
      m = len_trim(line1)
! Except for the following situations, like 6p+ 7p-, 7s 10s
! In all cases:
      if (n.eq.2) then
       char4a = line1(m-7:m-4)
       if (char4a(4:4).ne.'-') char4a(4:4) = '+'
       char4b = line1(m-16:m-13) 
       if (char4b(4:4).ne.'-') char4b(4:4) = '+'
       n1 = 0
       n2 = 0
       n1 = index(line4, trim(adjustl(char4a)))
       n2 = index(line4, trim(adjustl(char4b)))
       if (n1.eq.0 .or. n2.eq.0) then
         write(*,*)"Sth error in rcsfsym.f Line 225 ..."
         write(*,'(A)')trim(line1)
         write(*,'(A)')trim(line2)
         write(*,'(A)')trim(line3)
         write(*,'(A)')trim(line4)
         stop "Error, Unexpected n1.eq.0 .or. n2.eq.0 in rcsfsym.f ..."
       endif 
       if (n2.gt.n1) keep = 0
      endif 

      if (n.eq.2 .and. keep.eq.1) then
        string1 = line1(m-5:m-4)
        string3 = line1(m-14:m-13)
        if (string1.eq.string3) then
          string2 = line1(m-7:m-6)
          string4 = line1(m-16:m-15)
          read(string2, *)n1
          read(string4, *)n2
          ! TYPE 4 CSFs shoud meet n1 = n2 + 1
          if (n2+1 .ne. n1) keep = 0
        endif
      endif 
! End of all cases

! If Ntimes (nmrcon) = 2, keep TYPE 2, 3, 5 CSFs with nmax(i)/sym(i)
      if (nmrcon.eq.2) then
       if (n.eq.1) then
        ! TYPE 2 and 5
        nc = 0
        string1 = line1(m-5:m-4)         ! sym of symmetry-ordered orbital 1
        if (string1(2:2).ne.'-') string1(2:2) = '+'
        string2 = line1(m-7:m-6)         ! n of symmetry-ordered orbital 1
        do i = 1,21
          if (string1.eq.sym(i)) then
           if (string2.eq.nmax(i)) then
             nc = nc + 1
             exit
           endif
          end if
        end do
        if (nc.ne.1) keep = 0 

       elseif (n.eq.2 .and. keep.eq.1) then
        string1 = line1(m-5:m-4)
        if (string1(2:2).ne.'-') string1(2:2) = '+'
        string2 = line1(m-7:m-6)
        string3 = line1(m-14:m-13)
        if (string3(2:2).ne.'-') string3(2:2) = '+'
        string4 = line1(m-16:m-15)
        nc = 0
        do i = 1,21
         if (string1.eq.sym(i)) then
           if (string2.eq.nmax(i)) nc = nc + 1
         end if
         if (string3.eq.sym(i)) then
           if (string4.eq.nmax(i)) nc = nc + 1
         end if
        end do
        if (nc.ne.2 .and. string1.ne.string3) keep = 0 
       endif
      endif
! End of Ntimes = 2

! Write the labelling- and symmetry-ordered-CSFs into file 1312 and 1313,
! respectively. Count the discarded CSFs.
      if (keep.eq.0) ndiscard0 = ndiscard0 + 1
      if (keep.eq.1) then
         if (n.eq.0) then 
            iunit = 1312
         else
            iunit = 1313
         end if
         write(iunit,'(a)') trim(line1)
         write(iunit,'(a)') trim(line2)
         write(iunit,'(a)') trim(line3)
      end if

      goto 10
! End of reading the original CSFs.

  99  continue
      write(*,*)"Number of the discarded CSFs in rcsfsym:", ndiscard0

      rewind(19)
      rewind(1312)
      rewind(1313)

! Write the first THREE lines
      do i = 1,3
         write(19,'(a)') trim(string(i))
      end do

! Remove '+' from the string
      do i = 1,1500
        if (stringappend(i:i).eq.'+') stringappend(i:i) = ' '
      end do

! Write lines 4 and 5
      write(19,'(a)')                                                  &
              string(4)(1:5*(nonsym-nlengthcore/5))//trim(stringappend)
      write(19,'(a)') string(5)

! Start by writing all the labelling CSFs
   11 continue
      read (1312,'(a)', end = 998) line1
      read (1312,'(a)') line2
      read (1312,'(a)') line3
      write(19,'(a)') trim(line1)
      write(19,'(a)') trim(line2)
      write(19,'(a)') trim(line3)
      go to 11 
  998 continue 

! Here comes the symmetry-ordered CSFs, i.e., CSFG
   12 continue
      read (1313,'(a)', end = 999) line1
      read (1313,'(a)') line2
      read (1313,'(a)') line3
      write(19,'(a)') trim(line1)
      write(19,'(a)') trim(line2)
      write(19,'(a)') trim(line3)
      go to 12  
  999 continue     
      close(19)
      close(1312)
      close(1313)
      close(2425)

      end
    
