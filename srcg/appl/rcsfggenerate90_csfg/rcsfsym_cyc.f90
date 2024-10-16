!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     RCSFSYM
!
!     This subroutine removes un-needed CSFs built from only labelling
!     orbitals.
!      
!     Let n = 7 be the max principal quantum number for the
!     symmetry-ordered 
!     orbitals then 1s(2)6s7s should be kept wheras 1s(2)6s7d, 1s(2)7s6d 
!     should be removed and only 1s(2)7s7d kept      
!     
!     WRitten by Per Jonsson, Malmo University, Sweden
!     May 2017
!     
!     Modified by Chongyang Chen, Fudan University, China
!     May 2020
!     
!     Modified for the flexible rcsf-generation, such as the nmax is
!     different for different kind of correlation.
!     
!     Last Modification by CYC,  Dec. 2023
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!...Translated by  Gediminas Gaigalas  04/21/2021
! Last modification, Chongyang Chen, Fudan university,  May 2021

      subroutine rcsfsym_cyc
!
!
      implicit none
      character*1500 :: string(5),stringorb
      character*1500 :: stringnonsym  ! cyc
      character*300 :: line1,line2,line3
      character*2 :: string1,string2,string3,string4,nmax(21),sym(21) 
      character*2 :: nmaxnon(21),nstring
      integer :: i,j,nlength,nonsym,n,keep,nc,m,iunit
      integer :: nlengthcore           ! cyc
      integer :: nmaxi(21),nmaxnoni(21),n0(21)

! cyc 
      integer :: nblock,icsf,iflag1,iflag2,nkeep,ndiscard
      integer :: iblk,ncfmin,ncfmax,n1,n2
      integer :: TotCSFs_perblock(3,50),ntotcsf0
      integer, dimension(:), allocatable :: keepc
      integer, dimension(:), allocatable :: num_symorb
      logical, dimension(:), allocatable :: flag
      character*300, dimension(:), allocatable :: conf1
      character*300, dimension(:), allocatable :: conf2
      character*300, dimension(:), allocatable :: conf3

      character*300, dimension(:), allocatable :: conf1a
      integer, dimension(:), allocatable :: lconf1,lconf1a
      character*7  , dimension(:), allocatable :: conf1b
      character*6  , dimension(:), allocatable :: conf1c
      character*100, dimension(:), allocatable :: conf2a,conf3a
      integer, dimension(:), allocatable :: lconf2a,lconf3a
  
      integer keep1, discard1, TIME_BEGIN, TIME_END, TIME_RATE

!cyc==================================================================== 
      integer            :: nmax_spec, lmax_spec,nmaxcan, nonsym_free
      integer            :: orb_spec_nl(0:25)
      integer            :: nmax_orbset(100),nmrcon, ntime
      logical            :: freeorb
      common/orb_spec/nmax_spec,lmax_spec,nmaxcan,nonsym_free,       &
            orb_spec_nl,nmax_orbset,nmrcon,freeorb, ntime
!cyc==================================================================== 

      write(*,*)
      write(*,*) 'Weed out unneeded CSFs to obtain CSFGs in rcsfsym_cyc ...'
!      write(*,*) 'Freeorb=',freeorb

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
! the first n-value of the above each symmetry.
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

! Input 
      open (19, file = 'rcsfg.out', status = 'old', form = 'formatted', &
                action='readwrite')
! Output
! Temp files for each blocks
! CSFs built from the spectroscopy orbitals
!      open (1312, status = 'scratch', form = 'formatted')
! CSFs involving at least one of the correlation orbitals
!      open (1313, status = 'scratch', form = 'formatted')

! Read and analyze header to determine the highest principal quantum
! number for each symmetry. Also add orbitals to the header

      do i = 1, 5
         read (19,'(a)') string(i)
         ! First five lines have been generated correctly in rcsfsym
!         write(20,'(a)') trim(string(i))
      end do
      stringorb = trim(string(4))
      m = len_trim(stringorb)
      if (stringorb(m:m).ne.'-') stringorb=stringorb(1:m)//'+'
      m = m +1 
      if (mod(m,5).ne.0) then
       write(*,*)'m=',m
       stop "Unexpected mod(m,5).ne.0 ..."
      endif  
      m = m/5
      do i=1,m
        if (stringorb(i*5:i*5).eq.' ') stringorb(i*5:i*5) = '+'
      enddo

      nlength = len_trim(string(4)) + 1
! cyc, nlengthcore is used to determin the number of core-orbital. 
      nlengthcore = len_trim(string(2)) + 1
      if (nlengthcore.lt.5) nlengthcore = 0
! cyc

!DEBUG      write(*,*) trim(string(5))

!      write(*,*) ' Give the number of labelling orbitals'

!cyc: Flexible labelling set, using nonsym_free 
      nonsym = nonsym_free

! Analyze the labelling orbitals
! cyc =================================================================
! Some labelling orbitals are stored in string(2)
      stringnonsym = string(4)(1:nonsym*5)
      if (nlengthcore.gt.0)                                           &
       stringnonsym = trim(string(2)) // ' ' //                      &
                      string(4)(1:(nonsym-nlengthcore/5)*5)
      ! max n value for each labelling orbital
      do i = 1,nonsym
         string1 = stringnonsym(5*(i-1)+4:5*i)   ! symmetry string
         string2 = stringnonsym(5*(i-1)+2:5*i-2) ! n string
         do j = 1,21
            if (string1.eq.sym(j)) then
               nmaxnon(j) = string2
            end if
         end do
      end do
! cyc =================================================================
!      write(*,*)"Highest Labelling orbitals: "
!      write(*,'(A)')trim(stringnonsym)
!      do i = 1,21
!        if (nmaxnon(i).ne.'  ') then
!          write(*,*) nmaxnon(i),sym(i)
!        end if
!      end do

! Analyze the orbitals and determine the maximum n for each symmetry
! cyc =================================================================
! However, in the flexible CSF generation, the maximum n for VV may be
! different from the CV and CC ones. Just keep the following original
! lines, but keep in mind that nmax(j) should be used carefully.
!      write(*,*)"Highest symmetry-ordered orbitals: "
!!      do i = nonsym+1,nlength/5
!      do i = nonsym-nlengthcore/5+1,nlength/5
!! cyc =================================================================
!         string1 = string(4)(5*(i-1)+4:5*i)   ! symmmetry string
!         string2 = string(4)(5*(i-1)+2:5*i-2) ! n string
!         do j = 1,21
!            if (string1.eq.sym(j)) then
!               nmax(j) = string2
!            end if
!         end do
!      end do
!      do i = 1,21
!        if (nmax(i).ne.'  ') then
!          write(*,*) nmax(i),sym(i)
!        end if
!      end do

! Reading the lines in file 19, write the needed ones into file 20,
! dealt with Block by Block. 

! Section 1: Read through file 19, obtain the number of blocks, and the 
! number of the CSFs in each block.
! 
      TotCSFs_perblock = 0 
      nblock = 1
      do
        read(19,'(a)',end=7) line1
        if (line1(1:2).eq.' *') then
          nblock = nblock + 1
          read(19,'(a)') line1
        end if
        read(19,'(a)') line2
        read(19,'(a)') line3
        TotCSFs_perblock(1,nblock) = TotCSFs_perblock(1,nblock) + 1
      end do
7     continue
!      write(*,*)'nblock = ', nblock
!      write(*,*)TotCSFs_perblock(1,1:nblock)

      ntotcsf0 = 0
      do i = 1, nblock
        ntotcsf0 = ntotcsf0 + TotCSFs_perblock(1,i)
      enddo
      !write(*,*)'ntotcsf0 =', ntotcsf0

      ! Allocate the working arrays
      allocate(conf1(ntotcsf0))
      allocate(conf2(ntotcsf0))
      allocate(conf3(ntotcsf0))
      allocate(keepc(ntotcsf0))
      allocate(flag(ntotcsf0))
      ! Number of symmetry-ordered orbitals of each CSF 
      allocate(num_symorb(ntotcsf0))
      
      allocate(conf1a(ntotcsf0))
      allocate(conf1b(ntotcsf0))
      allocate(conf1c(ntotcsf0))
      allocate(conf2a(ntotcsf0))
      allocate(conf3a(ntotcsf0))
      conf1a = ' '
      conf1b = ' '
      conf1c = ' '
      conf2a = ' '
      conf3a = ' '

      allocate(lconf1(ntotcsf0))
      allocate(lconf1a(ntotcsf0))
      allocate(lconf2a(ntotcsf0))
      allocate(lconf3a(ntotcsf0))
      lconf1   = 0
      lconf1a  = 0
      lconf2a  = 0
      lconf3a  = 0

! Section 2: 
! Record the configurations (line 1 of each CSF)
! Count the number of symmetry-ordered orbitals in the CSF
! Set keepc value (partially)
      nkeep = 0
      ndiscard = 0
      icsf = 0

!      write(*,*)'Now reading CSFs ...'
!      CALL SYSTEM_CLOCK(TIME_BEGIN,TIME_RATE)
      rewind (19)
      do i = 1, 5 
        read(19,*)
      enddo

      do 100
        read(19,'(a)',end=8) line1
        if (line1(1:2).eq.' *') read(19,'(a)') line1
        read(19,'(a)') line2
        read(19,'(a)') line3
        icsf = icsf + 1 

        conf1(icsf) = trim(line1)
        conf2(icsf) = trim(line2)
        conf3(icsf) = trim(line3)
        keepc(icsf) = 0
        flag(icsf) = .false.
! Count the number of symmetry-ordered orbitals in the CSF 
        n = 0
        do i = nonsym-nlengthcore/5+1,nlength/5
           if (index(line1,string(4)(5*(i-1)+1:5*i)).gt.0) n = n + 1
        end do
        num_symorb(icsf) = n

        if (n.eq.0) then 
          ! Labelling CSFs, kept
          keepc(icsf) = 1
          flag(icsf) = .true.
          nkeep = nkeep + 1
        elseif (n.eq.1) then
          ! TYPE 2 or TYPE 5 CSFGs
          nc = 0
          m = len_trim(line1)
          string1 = line1(m-5:m-4)         ! symmetry of symmetry-ordered orbital 1
          string2 = line1(m-7:m-6)         ! n of symmetry-ordered orbital 1
          do i = 1,21
            if (string1.eq.sym(i)) then
             if (string2.eq.nmax(i)) then
               nc = nc + 1
               exit
             endif
            end if   
          end do   
          ! symmetry-ordered orbital with maxest n value of the sym should be kept
          if (nc.eq.1) then 
            ! [core] 10p-( 1), TYPE 2
            ! [core] 10p-( 2), TYPE 5
            keepc(icsf) = 1
            flag(icsf) = .true.
            nkeep = nkeep + 1
          !elseif (nc.eq.0) then
           ! [core] 7p-( 1), kept possibly
           ! [core] 7p-( 2), kept possibly
           ! to be dealt with later by the comparison between different CONF1(ICSF)
          endif
        elseif (n.eq.2) then 
          ! TYPE 3 or TYPE 4 CSFs
          m = len_trim(line1)
          string1 = line1(m-5:m-4)         ! sym of symmetry-ordered orbital 1
          if (string1(2:2).ne.'-') string1(2:2)='+'
          string2 = line1(m-7:m-6)         ! n of symmetry-ordered orbital 1
          string3 = line1(m-14:m-13)       ! sym of symmetry-ordered orbital 2
          if (string3(2:2).ne.'-') string3(2:2)='+'
          string4 = line1(m-16:m-15)       ! n of symmetry-ordered orbital 2
          n1 = 0
          n2 = 0
          n1 = index(stringorb,trim(adjustl(string2//string1)))
          n2 = index(stringorb,trim(adjustl(string4//string3)))
          if (n2.gt.n1) then
           ! CSFs like 6p 7p- should be discarded
           keepc(icsf) = 0
           flag(icsf) = .true.
           ndiscard = ndiscard + 1 
           cycle
          endif

          string1 = line1(m-5:m-4)
          string3 = line1(m-14:m-13)
          nc = 0 
          do i = 1,21
           if (string1.eq.sym(i)) then
             if (string2.eq.nmax(i)) nc = nc + 1
           end if   
           if (string3.eq.sym(i)) then
             if (string4.eq.nmax(i)) nc = nc + 1
           end if   
          end do

          if (nc.eq.2) then
           ! 10p- 10p  kept, TYPE 3 CSF   
           keepc(icsf) = 1
           flag(icsf) = .true.
           nkeep = nkeep + 1
          elseif (nc.eq.1) then
           if (string1.eq.string3) then
            ! TYPE 4 CSF
            ! 9p- 10p- kept, here nmax("p-") = "10"
            ! FLEXIBLE orbital set, this could be discarded. 
            ! For example: 
            !            Labelling set        : 5s,5p,4d,4f
            !            symmetry-ordered set : 7s,7p,7d,7f,7g,6h
            ! CSF: [core] 5f- 7f- should be discarded 
            ! Compare n1 and n2 to resolve this issue
            read(string2,*)n1
            read(string4,*)n2
            if (n2+1.eq.n1) then
              ! 6f- 7f- kept
              keepc(icsf) = 1
              flag(icsf) = .true.
              nkeep = nkeep + 1
            else
              ! 5f- 7f- discarded
              keepc(icsf) = 0
              flag(icsf) = .true.
              ndiscard = ndiscard + 1
            endif
           endif
          !elseif (nc.eq.0) then
          ! 7p- 7p could be possibly kept, to be dealt with later by
          ! the comparison between different CONF1(ICSF)
          endif
        else
          write(*,'(A)')trim(line1)
          write(*,'(A)')trim(line2)
          write(*,'(A)')trim(line3)
          write(*,*)'Error, Too many symmetry-ordered orbitals, n=', n
          stop &
      'Error, Too many symmetry-ordered orbitals in rcsfsym_cyc.f ...'
        endif
100   end do
8     continue
      CALL SYSTEM_CLOCK(TIME_END,TIME_RATE)
!      WRITE(*,*)'TIME USED A (s):', (TIME_END-TIME_BEGIN)/TIME_RATE

! Section 3: 
! In the cases of the nmax values are different for VV, CV, and CC
! electron correlation, the symmetry-ordered CSFs are determined by the
! comparison between different CSFs. The comparisons one by one is
! very time-consuming, improved by the followings.

! Section 3a:
! Improvement 1: 
! Removing the blanks in Conf1, Conf2, Conf3
! Obtain the strings Conf1a,Conf1b,Conf1c, and Conf2a, Conf3a
! Also the length for each string

!      write(*,*)'CSFs obtained, now pre-processing them 3A ...'
      CALL SYSTEM_CLOCK(TIME_BEGIN,TIME_RATE)

      do icsf = 1, ntotcsf0
       if (num_symorb(icsf).eq.0) then
         ! Kept
         cycle
       else
         line1 = trim(conf2(icsf))
         m = len_trim(line1)
         ! The following line cann't be used within OPENMP
         ! implementation. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !call removeblank(line1,m,conf2a(icsf),lconf2a(icsf))
         call removeblank(line1,m,line2,n)
         conf2a(icsf) = line2(1:n)
         lconf2a(icsf)= n

         line1 = trim(conf3(icsf))
         m = len_trim(line1)
         call removeblank(line1,m,line2,n)
         conf3a(icsf) = line2(1:n)
         lconf3a(icsf)= n

         line1 = conf1(icsf)
         m = len_trim(line1)
         lconf1(icsf) = m
         conf1c(icsf) = line1(m-5:m)
         if (num_symorb(icsf).eq.1) then
           call removeblank(line1,m-8,line2,n)
           conf1a(icsf) = line2(1:n)
           lconf1a(icsf)= n
         else
          conf1b(icsf) = line1(m-14:m-8)
          call removeblank(line1,m-17,line2,n)
          conf1a(icsf) = line2(1:n)
          lconf1a(icsf)=n
         endif 
       endif 
      enddo
      CALL SYSTEM_CLOCK(TIME_END,TIME_RATE)
!      WRITE(*,*)'TIME USED 3A (s):', (TIME_END-TIME_BEGIN)/TIME_RATE

! Section 3b:
! Remove unneeded CSFs by the detail comparisons

      CALL SYSTEM_CLOCK(TIME_BEGIN,TIME_RATE)
      ncfmin = 1
      do iblk = 1, nblock
       ncfmax = ncfmin + TotCSFs_perblock(1,iblk) - 1
       discard1 = 0
       keep1    = 0
       do icsf = ncfmin, ncfmax
        if (mod(icsf,50000).eq.0)                                     &
      write(*,*)'nblock, ntotcsf0, iblock, ncfmin, ncfmax, icsf=',    &
                nblock, ntotcsf0, iblk, ncfmin, ncfmax, icsf
        if (flag(icsf)) cycle
        flag(icsf) = .true.
        keepc(icsf) = 1
        n = num_symorb(icsf) 
        !Keep or discard
        do i = icsf + 1, ncfmax
!===================================================================== 
!          ! The CSFs with the same conf1a string are grouped together
!          ! in the original rcsf.out, hence,  
           if (.not.freeorb) then
            if (lconf1a(icsf).ne.lconf1a(i)  .or.                     &
             conf1a(icsf)(1:lconf1a(icsf)) .ne.                       &
             conf1a(i)(1:lconf1a(i))) exit
           else 
          ! The above assumption should be checked carefully in case of
          ! the calculation performed for one new isoelectronic sequence
          ! ion. Otherwise, use the following three lines: (very
          ! time-consuming) 
            if (lconf1a(icsf).ne.lconf1a(i))  cycle
            if (conf1a(icsf)(1:lconf1a(icsf))                         &
                .ne.conf1a(i)(1:lconf1a(i))) cycle
           endif  
!===================================================================== 
          ! Compare the other integers
          if (num_symorb(i).ne.n) cycle
          if (lconf1(icsf) .ne. lconf1(i)) cycle
          if (lconf2a(icsf).ne.lconf2a(i)) cycle
          if (lconf3a(icsf).ne.lconf3a(i)) cycle 

          ! Compare the other strings
          if (conf1c(icsf) .ne. conf1c(i)) cycle
          if (n.eq.2 .and. conf1b(icsf).ne.conf1b(i)) cycle
          if (conf3a(icsf)(1:lconf3a(icsf)) .ne.                       &
              conf3a(i)(1:lconf3a(i))) cycle
          if (conf2a(icsf)(1:lconf2a(icsf)) .ne.                       &
              conf2a(i)(1:lconf2a(i))) cycle

          ! Here, find another CSFG (the i^{th} CSF) which could
          ! reproduce the icsf-th one
          keepc(icsf) = 0
          exit
        enddo
        if (keepc(icsf).eq.1) then
         ! No similar CSFG found, keep it
         keep1 = keep1 + 1
        else 
         discard1 = discard1 + 1
        endif
       enddo 
       ncfmin = ncfmax + 1
       nkeep = nkeep + keep1
       ndiscard = ndiscard + discard1
      enddo
      CALL SYSTEM_CLOCK(TIME_END,TIME_RATE)
!      WRITE(*,*)'TIME USED 3B (s):', (TIME_END-TIME_BEGIN)/TIME_RATE
!      write(*,*)
      write(*,*)"Total CSFs entering  rcsfsym_cyc:", ntotcsf0
      write(*,*)"Weed out CSFs within rcsfsym_cyc:", ndiscard
      write(*,*)"CSF(G)s kept:", nkeep
      if (ndiscard + nkeep /= ntotcsf0) &
        write(*,*)"Warning!!! ndiscard + nkeep /= ntotcsf0 ..."

      rewind(19)
      do i=1,5
        read(19,*)
      enddo
      
      keep1 = 0
      ncfmin = 1
      do iblk = 1, nblock
        ncfmax = ncfmin + TotCSFs_perblock(1,iblk) - 1
        do icsf = ncfmin, ncfmax
          if (keepc(icsf).eq.1) then
            keep1 = keep1 + 1
            write(19,'(a)')trim(conf1(icsf))
            write(19,'(a)')trim(conf2(icsf))
            write(19,'(a)')trim(conf3(icsf))
          endif
        enddo
        ncfmin = ncfmax + 1
        if (iblk.lt.nblock) write(19, '(a2)')' *'
      enddo
      close(19)
      if (keep1 /= nkeep) &
         write(*,*)"Warning !!! keep1 /= nkeep within rcsfsym_cyc.f90"

      return
      end
