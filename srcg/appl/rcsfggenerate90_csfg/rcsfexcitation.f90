!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! RCSFEXCITATION
!
! This subroutine generates excitation input to jjgen
!
! Per Jonsson, Malmo University, October 2013
!                                Corrected Feb 2017
! 
! Modified for FLEXIBLE orbital set,
! Chongyang Chen, Fudan university,  Jan 2021
!
!     
!     Last Modification by CYC,  Dec. 2023
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!...The CSFG version was modified by  Gediminas Gaigalas  04/21/2021
!
! Last modification, Chongyang Chen, Fudan university,  May 2023
 
subroutine rcsfexcitation
implicit none
integer            :: i,k,j,nmax,lmax,mr,n(25),l(25),occ(25),nr,lr,occr,nl(0:25)
integer            :: jmin,jmax,nexc,occupied,ncore,nstart(0:6),lstart(0:6)
integer            :: nc(30),lc(30),ncoreloop,nclose,norbitalstring
integer            :: conflength,norb,orblength(25),orbstart(0:25),flag(25),nfound
integer            :: number1,number2,norbstrings,orbstring_l(11),orbstring_n(11)
integer            :: ncoreorbitals,lmaxcore,ntimes,ndouble,jl,jr,pos,ncore5to4
integer            :: ncore6to5,nfoundf,nfoundd

character(len=1)   :: lstring,ans
character(len=2)   :: sel(25),nstring
character(len=135) :: config
character(len=135) :: configvect(1000)
character(len=44)  :: orbitalstring
character(len=11)  :: decoding
character(len=3)   :: orbital(11)
character(len=22)  :: corestring5
character(len=30)  :: corestring6

!cyc====================================================================
integer            :: nmax_spec, lmax_spec,nmaxcan, nonsym_free
integer            :: orb_spec_nl(0:25), ntime
integer            :: nmax_orbset(100),nmrcon
logical            :: freeorb
common/orb_spec/nmax_spec,lmax_spec,nmaxcan,nonsym_free,              &
                orb_spec_nl,nmax_orbset,nmrcon,freeorb,ntime
integer            :: nminus1_orbs
!cyc====================================================================

! Including the modifications provided by Jon to allow for the generation of
! CSF's from MR sets including configurations with n>9 
integer            :: stat, n_len

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call starttime (ntime,'RCSFGGenerate_csfg')

decoding = 'spdfghiklmn'

open(11, file='excitationdata',status='unknown',form='formatted')
open(12, file='rcsfggenerate.log',status='unknown',form='formatted')

! For CSFG, generate clist.ref

open(113, file='clist.ref',status='unknown',form='formatted')

! Generate initial input to jjgen

! New list
write(11,'(a)') '*'

nc = 0
lc = 0
ncoreorbitals = 0
lmaxcore = 0
freeorb = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*)
write(*,*) 'RCSFGGENERATE_CSFG'
write(*,*) 'This program generates a list of CSFG CSFs'
write(*,*)
write(*,*) 'Configurations should be entered in spectroscopic notation'
write(*,*) 'with occupation numbers and indications if orbitals are'
write(*,*) 'closed (c), inactive (i), active (*) or has a minimal'
write(*,*) 'occupation e.g. 1s(2,1)2s(2,*)'
write(*,*) 'Outputfiles: rcsfg.out rlabel.out rcsfggenerate.log rcsfg.log clist.ref excitation'
write(*,*)
!write(*,*) 'Confirm user specified ordering available in clist.ref (u)'
!read(*,*) ans
!if (ans.ne.'u') stop
write(11,'(a)') 'u'
!write(12,*) 'u', ' ! Orbital order'
write(*,*)
write(*,*) 'Select core '
write(*,*) '       0: No core'
write(*,*) '       1: He (       1s(2)                  =  2 electrons)'
write(*,*) '       2: Ne ([He] + 2s(2)2p(6)             = 10 electrons)'
write(*,*) '       3: Ar ([Ne] + 3s(2)3p(6)             = 18 electrons)'
write(*,*) '       4: Kr ([Ar] + 3d(10)4s(2)4p(6)       = 36 electrons)'
write(*,*) '       5: Xe ([Kr] + 4d(10)5s(2)5p(6)       = 54 electrons)'
write(*,*) '       6: Rn ([Xe] + 4f(14)5d(10)6s(2)6p(6) = 86 electrons)'
read(*,*) ncore
write(12,*) ncore, ' ! Selected core'

! Due to restrictions in the underlying jjgen we have do use the following "trick"
! if ncore == 5 then we redefine the core to ncore = 4 and add
! 4d(10,c)5s(2,c)5p(6,c) to each of the configuration strings
! if ncore == 6 then we redefine the core to ncore = 5 and add
! 4f(14,c)5d(10,c)6s(2,c)6p(6,c) to each of the configuration strings
! Per J Feb 2017

ncore5to4 = 0
ncore6to5 = 0
if (ncore.eq.5) then
   ncore = 4
   ncore5to4 = 1
end if
if (ncore.eq.6) then
   ncore = 5
   ncore6to5 = 1
end if

corestring5 = '4d(10,c)5s(2,c)5p(6,c)'
corestring6 = '4f(14,c)5d(10,c)6s(2,c)6p(6,c)'
write(*,*)

! End correction 2017

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Analyze core and define the orbitals in the core

select case(ncore)
   case(1)
      ncoreorbitals = 1   ! 1 core orbital
      lmaxcore = 0        !
      nc(1) = 1           ! 1s
      lc(1) = 0           ! 1s
   case(2)
      ncoreorbitals = 3   ! 3  core orbitals
      lmaxcore = 1        !
      nc(1) = 1           ! 1s
      lc(1) = 0           ! 1s
      nc(2) = 2           ! 2s
      lc(2) = 0           ! 2s
      nc(3) = 2           ! 2p
      lc(3) = 1           ! 2p
   case(3)
      ncoreorbitals = 5   ! 5 core orbitals
      lmaxcore = 1        !
      nc(1) = 1           ! 1s
      lc(1) = 0           ! 1s
      nc(2) = 2           ! 2s
      lc(2) = 0           ! 2s
      nc(3) = 2           ! 2p
      lc(3) = 1           ! 2p
      nc(4) = 3           ! 3s
      lc(4) = 0           ! 3s
      nc(5) = 3           ! 3p
      lc(5) = 1           ! 3p
   case(4)
      ncoreorbitals = 8   ! 8 core orbitals
      lmaxcore = 2        !
      nc(1) = 1           ! 1s
      lc(1) = 0           ! 1s
      nc(2) = 2           ! 2s
      lc(2) = 0           ! 2s
      nc(3) = 2           ! 2p
      lc(3) = 1           ! 2p
      nc(4) = 3           ! 3s
      lc(4) = 0           ! 3s
      nc(5) = 3           ! 3p
      lc(5) = 1           ! 3p
      nc(6) = 3           ! 3d
      lc(6) = 2           ! 3d
      nc(7) = 4           ! 4s
      lc(7) = 0           ! 4s
      nc(8) = 4           ! 4p
      lc(8) = 1           ! 4p
   case(5)
      ncoreorbitals = 11  ! 11 core orbitals
      lmaxcore = 2        !
      nc(1) = 1           ! 1s
      lc(1) = 0           ! 1s
      nc(2) = 2           ! 2s
      lc(2) = 0           ! 2s
      nc(3) = 2           ! 2p
      lc(3) = 1           ! 2p
      nc(4) = 3           ! 3s
      lc(4) = 0           ! 3s
      nc(5) = 3           ! 3p
      lc(5) = 1           ! 3p
      nc(6) = 3           ! 3d
      lc(6) = 2           ! 3d
      nc(7) = 4           ! 4s
      lc(7) = 0           ! 4s
      nc(8) = 4           ! 4p
      lc(8) = 1           ! 4p
      nc(9) = 4           ! 4d
      lc(9) = 2           ! 4d
      nc(10) = 5          ! 5s
      lc(10) = 0          ! 5s
      nc(11) = 5          ! 5p
      lc(11) = 1          ! 5p
   case(6)
      ncoreorbitals = 15  ! 15 core orbitals
      lmaxcore = 3        !
      nc(1) = 1           ! 1s
      lc(1) = 0           ! 1s
      nc(2) = 2           ! 2s
      lc(2) = 0           ! 2s
      nc(3) = 2           ! 2p
      lc(3) = 1           ! 2p
      nc(4) = 3           ! 3s
      lc(4) = 0           ! 3s
      nc(5) = 3           ! 3p
      lc(5) = 1           ! 3p
      nc(6) = 3           ! 3d
      lc(6) = 2           ! 3d
      nc(7) = 4           ! 4s
      lc(7) = 0           ! 4s
      nc(8) = 4           ! 4p
      lc(8) = 1           ! 4p
      nc(9) = 4           ! 4d
      lc(9) = 2           ! 4d
      nc(10) = 5          ! 5s
      lc(10) = 0          ! 5s
      nc(11) = 5          ! 5p
      lc(11) = 1          ! 5p
      nc(12) = 4          ! 4f
      lc(12) = 3          ! 4f
      nc(13) = 5          ! 5d
      lc(13) = 2          ! 5d
      nc(14) = 6          ! 6s
      lc(14) = 0          ! 6s
      nc(15) = 6          ! 6p
      lc(15) = 1          ! 6p
end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ntimes = 1
110 continue
n = 0
l = 0
occ = 0

! Input the mutireference for later processing

if (ntimes.eq.1) then
  write(*,*) 'Generation of CSFs in labeling space'
  write(*,*) 'Entered orbitals define the labeling orbital set'
else
  write(*,*) 'Generation of CSFs in labeling space and, if entered orbitals extend '
  write(*,*) 'the labeling orbital set, also defining CSFs in correlation space'
  write(*,*) 'Orbitals extending the labeling orbital set appear symmetry ordered'
end if

write(*,*) 'Enter list of (maximum 1000) configurations.'
write(*,*) 'End list with a blank line or an asterisk (*)'
write(*,*)

do j = 1,1000

10 continue
   write(*,*) 'Give configuration', j
   read(*,'(a)') configvect(j)
   if (trim(configvect(j)).eq.''.or.trim(configvect(j)).eq.'*') then
      mr = j-1
      write(12,'(a)') '*'
      goto 99
   end if

!  Initial check, each orbital need to be closed, inactive, or minimal

   conflength = len(trim(configvect(j)))
   number1 = 0
   number2 = 0
   do k = 1,conflength
      if (configvect(j)(k:k).eq.'(') number1 = number1 + 1
      if (configvect(j)(k:k).eq.',') number2 = number2 + 1
   end do
   if (number1.ne.number2) then
      write(*,*)  'Each orbital must be closed (c), inactive (i), active (*)'
      write(*,*)  'or have a minimal occupation; redo!'
      goto 10
   end if
   write(12,'(a)') trim(configvect(j))

! Here we append string if we have redefined core Per J Feb 2017

   if (ncore5to4.eq.1) then
      configvect(j) = trim(corestring5//configvect(j))
      if (conflength + 22.gt.135) then
         write(*,*) 'In rcsfexcitation, increase length of config and configvect'
         stop
      end if
   end if
   if (ncore6to5.eq.1) then
      configvect(j) = trim(corestring6//configvect(j))
      if (conflength + 30.gt.135) then
         write(*,*) 'In rcsfexcitation, increase length of config and configvect'
         stop
      end if
   end if

! End correction Feb 2017

end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Input and process orbital set

99 continue
if (ntimes.eq.1) then
  write(*,*) 'Give set of active orbitals, as defined by the highest principal quantum number'
  write(*,*) 'in a comma delimited list in s,p,d etc order, e.g. 3s,3p,3d or 4s,4p,4d,4f or'
  write(*,*) '5s,5p,5d,5f,5g etc'
else
  write(*,*) 'Give set of active orbitals, as defined by the highest principal quantum number'
  write(*,*) 'per l-symmetry, in a comma delimited list in s,p,d etc order, e.g.  5s,4p,3d'
end if

read(*,'(a)') orbitalstring

! If ncore5to4 = 1 and the orbitalstring only contains 4f then we need to add 5s

if (ncore5to4.eq.1) then
   if (trim(orbitalstring).eq.'4f') then
      orbitalstring = '5s,4f'
!      write(*,'(a)') trim(orbitalstring)
   end if
end if

! If ncore6to5 = 1 we again rely on a trick: the orbitals string must contain an f orbital.
! If not we just insert 3f. In addition, if all orbitals in the string has n = 5 then we
! have to add 6s

if (ncore6to5.eq.1) then
   if ( (trim(orbitalstring).eq.'5f').or.(trim(orbitalstring).eq.'5f,5g').or. &
        (trim(orbitalstring).eq.'5g') ) then
      orbitalstring = '6s,'//trim(orbitalstring)
   end if

   nfoundf = 0
   do i = 1,len_trim(orbitalstring)
      if (orbitalstring(i:i).eq.'f') then
         nfoundf = 1
      end if
   end do
   if (nfoundf.eq.0) then
   ! Did not find f. Insert 3f either after d symmetry or before g symmetry
      nfoundd = 0
      do i = 1,len_trim(orbitalstring)
         if (orbitalstring(i:i).eq.'d') then
            nfoundd = 1
            if (orbitalstring(i+1:i+1).eq.',') then
               orbitalstring = orbitalstring(1:i)//',3f'//orbitalstring(i+1:len_trim(orbitalstring))
            else
               orbitalstring = orbitalstring(1:i)//',3f'
            end if
            exit
         end if
      end do
      if (nfoundd.eq.0) then
         orbitalstring = '3f,'//trim(orbitalstring)
      end if
   end if
!   write(*,'(a)') trim(orbitalstring)
!   pause
end if

! De-code orbital string

! Checks on orbital string

jl=1;
jr=0;
k=1;
do while ( index(orbitalstring(jl:len_trim(orbitalstring)),',').gt.0 )
   pos = index(orbitalstring(jl:len_trim(orbitalstring)),',')
   jr = jr + pos
   if (len_trim(orbitalstring(jl:jr-1)).eq.3) then
      ! "10s"
      orbital(k) = orbitalstring(jl:jr-1)
   else if(len_trim(orbitalstring(jl:jr-1)).eq.2) then
      ! " 9s"
      orbital(k) = " " // orbitalstring(jl:jr-1)
   else
      write(*,*) 'Orbitals should be given in comma delimited list, redo!'
      goto 99
   end if
   jl = jr + 1
   k = k +1
end do
! The last one of the orbitalstring
jr = len_trim(orbitalstring)
if (len_trim(orbitalstring(jl:jr)).eq.3) then
   orbital(k) = orbitalstring(jl:jr)
else if (len_trim(orbitalstring(jl:jr)).eq.2) then
   orbital(k) = " " // orbitalstring(jl:jr)
else
   write(*,*) 'Orbitals should be given in comma delimited list, redo!'
   goto 99
end if
write(12,'(a)') trim(orbitalstring)

norbstrings = k

! Determine highest n for orbital set and highest n for each symmetry

nmax = 0
do i = 1,norbstrings
!   write(*,'(a)') orbital(i)
  select case(orbital(i)(3:3))
   case('s')
      lmax = 0
   case('p')
      lmax = 1
   case('d')
      lmax = 2
   case('f')
      lmax = 3
   case('g')
      lmax = 4
    case('h')
      lmax = 5
   case('i')
      lmax = 6
   case('k')
      lmax = 7
   case('l')
      lmax = 8
   case('m')
      lmax = 9
   case('n')
      lmax = 10
   case default
      write(*,*) 'Orbital quantum numbers should be in range s to n'
      stop
  end select
  orbstring_l(i) = lmax

  select case(orbital(i)(1:2))
   case(' 1')
      nmax = 1
   case(' 2')
      nmax = 2
   case(' 3')
      nmax = 3
   case(' 4')
      nmax = 4
   case(' 5')
      nmax = 5
   case(' 6')
      nmax = 6
   case(' 7')
      nmax = 7
   case(' 8')
      nmax = 8
   case(' 9')
      nmax = 9
   case('10')
      nmax = 10
   case('11')
      nmax = 11
   case('12')
      nmax = 12
   case('13')
      nmax = 13
   case('14')
      nmax = 14
   case('15')
      nmax = 15
   case('16')
      nmax = 16
   case('17')
      nmax = 17
   case('18')
      nmax = 18
   case('19')
      nmax = 19
   case('20')
      nmax = 20
   case('21')
      nmax = 21
   case('22')
      nmax = 22
   case('23')
      nmax = 23
   case('24')
      nmax = 24
   case('25')
      nmax = 25
   case default
      write(*,*) 'Principal quantum numbers should be <= 25'
      stop
   end select
   orbstring_n(i) = nmax
end do

! Determine highest n and l for orbital set as well as spectrocopic notation

nmax = maxval(orbstring_n(1:norbstrings))
if (lmax.lt.lmaxcore) lmax = lmaxcore
lstring = decoding(lmax+1:lmax+1)

!cyc, determin the n and l (maximum) values of the labelling orbitals
if (ntimes.eq.1) then
  nmax_spec = nmax
  lmax_spec = lmax
endif
nmrcon = ntimes
nmax_orbset(nmrcon) = nmax

! Initialize nl
nl=0
!cyc

do i = 0,lmax
   nfound = 0
   do j = 1,norbstrings
      if (orbstring_l(j).eq.i) then
        nl(i) = orbstring_n(j)
        nfound = 1
      end if
   end do
   if (nfound.eq.0) nl(i) = 0
end do

!cyc
if (ntimes.eq.1) then
  orb_spec_nl=0
  orb_spec_nl(0:lmax)=nl(0:lmax)
  freeorb = .false.
  do i = 0,lmax-1
    if (orb_spec_nl(i) < orb_spec_nl(i+1)) freeorb = .true.
  enddo
endif
!cyc

if (lmax.ge.nmax) then
   write(*,*) 'Orbital quantum number should be less than n'
   stop
end if

! Write the canonical orbitals to clist.ref
! cyc: write only the orbitals which are not higher than those of orbitalstring,
! such as 5s, 4p, 4d, and 4f, ONLY 1s-5s, 2p-4p, 3d-4d, and 4f are written into
! file 113.
if (ntimes.eq.1) nonsym_free=0
if (ntimes.eq.1) then
  nmaxcan = nmax
  nminus1_orbs = 0  ! Record the total orbitals lower than nmaxcan
!PERJ  do i = 0, 10
!PERJ    write(*,*)'i=',i, ' nl(i)=',nl(i)
!PERJ  enddo
  do i = 1,nmax
    write (nstring,'(i2)') i
    do j = 0,i-1
! cyc
      if (i.gt.nl(j)) cycle
      ! nonsym_free: relativistic orbitals
      if (j.eq.0) then
        nonsym_free=nonsym_free+1
      else
        nonsym_free=nonsym_free+2
      endif
!PERJ      write(*,*)'nmax=',nmax,' n=',i, ' l=',j, ' nrelorbs=',nonsym_free
! cyc
      if (j.eq.0) then
        write(113,'(a)') trim(adjustl(nstring))//'s'
      elseif (j.eq.1) then
        write(113,'(a)') trim(adjustl(nstring))//'p'
      elseif (j.eq.2) then
        write(113,'(a)') trim(adjustl(nstring))//'d'
      elseif (j.eq.3) then
        write(113,'(a)') trim(adjustl(nstring))//'f'
      elseif (j.eq.4) then
        write(113,'(a)') trim(adjustl(nstring))//'g'
      elseif (j.eq.5) then
        write(113,'(a)') trim(adjustl(nstring))//'h'
      elseif (j.eq.6) then
        write(113,'(a)') trim(adjustl(nstring))//'i'
      elseif (j.eq.7) then
        write(113,'(a)') trim(adjustl(nstring))//'k'
      elseif (j.eq.8) then
        write(113,'(a)') trim(adjustl(nstring))//'l'
      elseif (j.eq.9) then
        write(113,'(a)') trim(adjustl(nstring))//'m'
      elseif (j.eq.10) then
        write(113,'(a)') trim(adjustl(nstring))//'n'
      end if
    end do
    if (i.eq.nmax-1) nminus1_orbs = nonsym_free
  end do
!PERJ  write(*,*)"nminus1_orbs=",nminus1_orbs
  ! freeorb used in rcsfsym_cyc.f90, determining how to remove the unneeded CSFs
  if (nminus1_orbs.eq.(nmax-1)*(nmax-1)) then
    if (.not.freeorb) freeorb = .false.
  else
    freeorb = .true.
  endif
!PERJ  write(*,*)"Freeorb=",freeorb
end if

write(*,*) 'Resulting 2*J-number? lower, higher (J=1 -> 2*J=2 etc.)'
read(*,*) jmin,jmax
write(12,*) jmin,jmax, ' ! Lower and higher 2*J'
write(*,*) 'Number of excitations (if negative number e.g. -2, correlation '
write(*,*) 'orbitals will always be doubly occupied)                        '
read(*,*) nexc
ndouble = 0
if (nexc.lt.0) then
   ndouble = 1
end if
write(12,*) nexc, ' ! Number of excitations '

!write(*,*)
!write(*,*) 'TESTING ORBITAL SET'
!do i = 0,lmax
!   write(*,*) 'orbital set, l, nmax',i,nl(i)
!end do
!write(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do j = 1,mr

! Highest principal quantum number
   write(11,*) nmax

! Highest orbital quantum number
   write(11,'(a)') lstring

! All these orbitals active
   write(11,'(a)') 'n'

! Limitations of population
   write(11,'(a)') '*'

   config = configvect(j)

! Analyze config and extract n, l and occupation

   ! Determine length of orbit substrings. Depends on number of positions
   ! that occupation and type indicator occupies.
   conflength = len(trim(config))
   k = 0
   n = 0 ! CYC: initialization
   orbstart(0) = 0
   do i = 1, conflength
      if(config(i:i).eq.')') then
         k = k + 1
         orblength(k) = i - orbstart(k-1)
         orbstart(k) = i
      end if
   end do
   norb = k  ! Number of orbit substrings in configuration

   flag(:) = 0
   do k = 1,norb
      i = 1+orbstart(k-1)
      if(config(i+4:i+4).eq.',') flag(k) = 1 ! Flag to facilitate determination if
                                             ! occupation/type indicator occupies one/two positions

      ! Find length of the n integer
      if (config(i + 2:i + 2) .eq. '(') then
         n_len = 1
      else
         n_len = 2
      end if

      ! get principal qn
      read (config(i:i + n_len - 1), *, iostat=stat) nr
      !select case(config(i:i))
      !   case('1')
      !      nr = 1
      !   case('2')
      !      nr = 2
      !   case('3')
      !      nr = 3
      !   case('4')
      !      nr = 4
      !   case('5')
      !      nr = 5
      !   case('6')
      !      nr = 6
      !   case('7')
      !      nr = 7
      !   case('8')
      !      nr = 8
      !   case('9')
      !      nr = 9
      !end select

      !select case(config(i+1:i+1))
      select case(config(i + n_len : i + n_len))
         case('s')
            lr = 0
         case('p')
            lr = 1
         case('d')
            lr = 2
         case('f')
            lr = 3
         case('g')
            lr = 4
         case('h')
            lr = 5
         case('i')
            lr = 6
         case('k')
            lr = 7
         case('l')
            lr = 8
         ! CYC:
         case('m')
            lr = 9
      end select

      ! If occupation number "occupies" one position
      if (flag(k) .eq. 1) then
         read (config(i + 2 + n_len:i + 2 + n_len), *, iostat=stat) occr
      else
         read (config(i + 2 + n_len:i + 3 + n_len), *, iostat=stat) occr
      end if
      !if(flag(k).eq.1) then
      !   select case(config(i+3:i+3))
      !     case('1')
      !        occr = 1
      !     case('2')
      !        occr = 2
      !     case('3')
      !        occr = 3
      !     case('4')
      !        occr = 4
      !     case('5')
      !        occr = 5
      !     case('6')
      !        occr = 6
      !     case('7')
      !        occr = 7
      !     case('8')
      !        occr = 8
      !     case('9')
      !        occr = 9
      !    end select
      !else
      !   select case(config(i+3:i+4))
      !     case(' 1')
      !        occr = 1
      !     case(' 2')
      !        occr = 2
      !     case(' 3')
      !        occr = 3
      !     case(' 4')
      !        occr = 4
      !     case(' 5')
      !        occr = 5
      !     case(' 6')
      !        occr = 6
      !     case(' 7')
      !        occr = 7
      !     case(' 8')
      !        occr = 8
      !     case(' 9')
      !        occr = 9
      !     case('10')
      !        occr = 10
      !     case('11')
      !        occr = 11
      !     case('12')
      !        occr = 12
      !     case('13')
      !        occr = 13
      !     case('14')
      !        occr = 14
      !    end select
      ! end if
      if (nr.gt.nmax) then
         write(*,*) 'n,nmax',nr,nmax
         write(*,*) 'n in config greater than nmax'
         stop
      end if

      n(k) = nr
      l(k) = lr
      occ(k) = occr

      ! Determine start position and number of positions for type indicator
      if(orblength(k).eq.7) then
         sel(k) = config(i+5:i+5)
      elseif((orblength(k).eq.8).and.(flag(k).eq.0)) then
         sel(k) = config(i+6:i+6)
      elseif((orblength(k).eq.8).and.(flag(k).eq.1)) then
         sel(k) = config(i+5:i+6)
      else
         sel(k) = config(i+6:i+7)
      end if

!CPJ  Orbitalerna i configurationen med occupation and indicators
!PERJ      write(*,*) 'n,l,occ,sel',n(k),l(k),occ(k), sel(k)

   end do

! Generate input to jjgen from this configuration

! Highest n for reference configuration
   write(11,*) maxval(n)

! Predefined core, write only once

   if ((j.eq.1).and.(ntimes.eq.1)) then
      if (ncore.eq.0) then
         write(11,'(a)') '*'
      else
         write(11,'(a)') 'c'
         write(11,*) ncore
      end if
   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   write(*,*) 'TESTING CONFIG',j

!   do k = 1, ncoreorbitals
!      write(*,*) 'coreorbitals',k,nc(k),lc(k)
!   end do

!   do k = 1, norb
!      write(*,*) 'orbitals with occupation ',k,n(k),l(k),occ(k),sel(k)
!   end do

!   do nr = 1,maxval(n)
!      do lr = 0,min(nr-1,lmax)
!         write(*,*) 'orbitals upp to nmax of configuration',nr,lr
!      end do
!   end do

!   do nr = maxval(n)+1,nmax
!      do lr = 0,min(nr-1,lmax)
!         write(*,*) 'orbitals upp to total nmax',nr,lr
!      end do
!   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do nr = 1,maxval(n)
      do lr = 0,min(nr-1,lmax)

! Check if this orbital was previously defined as closed. If so do not write anything for this

          nclose = 0
          do k = 1,ncoreorbitals
             if ((nr.eq.nc(k)).and.(lr.eq.lc(k))) then
                nclose = 1
             end if
          end do
          if (nclose.eq.1) cycle

! Find out if orbital is closed then add it to core orbitals
! In for the next reference configuration there will be no question for this

          do k = 1, norb
             if ((nr.eq.n(k)).and.(lr.eq.l(k))) then
                if (trim(adjustl(sel(k))).eq.'c') then
                   ncoreorbitals = ncoreorbitals + 1
                   nc(ncoreorbitals) = n(k)
                   lc(ncoreorbitals) = l(k)
                end if
             end if
          end do

! Find out occupation and options for this orbital

          occupied = 0
          do k = 1, norb
             if ((nr.eq.n(k)).and.(lr.eq.l(k))) then
                write(11,*) occ(k), '  # nr, lr=', nr, lr
                write(11,'(a)') trim(adjustl(sel(k)))
                occupied = 1
             end if
          end do
          if (occupied.eq.0) then
             write(11,*) 0, '  # nr, lr=', nr, lr
!CPJ             write(11,'(a)') '*'
!Correction
             if (nr.le.nl(lr)) then
                if (ndouble.eq.1) then
                   write(11,'(a)') 'd'
                else
                   write(11,'(a)') '*'
                end if
             else
                write(11,'(a)') 'i'
             end if
!Correction
          end if

      end do
   end do

! Loop over remaining orbitals in the orbital set

   if (nmax.gt.maxval(n)) then
      do nr = maxval(n)+1,nmax
         do lr = 0,min(nr-1,lmax)
            if (ntimes.gt.1) then                              ! symmetry-ordered Orbitals part
              if (((nl(lr)-1).le.nr).and.(nr.le.nl(lr))) then  ! Only highest and second highest n
                if (ndouble.eq.1) then
                   write(11,'(a)') 'd'
                else
                   write(11,'(a1,14x,a9,2i12)') '*', '# nr, lr=', nr, lr
                end if
              else
! cyc: make all labelling orbitals active 
                 if (nr.le.nmax_spec.and.nr.le.nl(lr)) then
                     write(11,'(a1,14x,a9,2i12)') '*', '# nr, lr=', nr, lr
                 else
                     write(11,'(a1,14x,a9,2i12)') 'i', '# nr, lr=', nr, lr
                 endif
              end if
            else
              if (nr.le.nl(lr)) then
                if (ndouble.eq.1) then
                  write(11,'(a)') 'd'
                else
                  write(11,'(a1,14x,a9,2i12)') '*', '# nr, lr=', nr, lr
                end if
              else
                write(11,'(a1,14x,a9,2i12)') 'i', '# nr, lr=', nr, lr
! cyc
              end if
            end if
         end do
      end do
   end if

! Write JMIN, JMAX and number of excitations

   write(11,*) jmin,jmax
   write(11,*) iabs(nexc)

   if (j.ne.mr) then
      write(11,'(a)') 'y'
   else
      write(*,'(a)') ' Generate more lists ? (y/n)'
      read(*,*) ans
      write(12,'(a)') ans
      if ((ans.eq.'y').or.(ans.eq.'Y')) then
         write(11,'(a)') 'y'
         ntimes = ntimes + 1
         goto 110
      else
         write(11,'(a)') '*'
      end if
   end if

end do

close(11)

do j = 0,10
! cyc  do i = max(nmaxcan+1,j+1),15
    do i = max(orb_spec_nl(j)+1, j+1), 25
    write (nstring,'(i2)') i
    if (j.eq.0) then
      write(113,'(a)') trim(adjustl(nstring))//'s'
    elseif (j.eq.1) then
      write(113,'(a)') trim(adjustl(nstring))//'p'
    elseif (j.eq.2) then
      write(113,'(a)') trim(adjustl(nstring))//'d'
    elseif (j.eq.3) then
      write(113,'(a)') trim(adjustl(nstring))//'f'
    elseif (j.eq.4) then
      write(113,'(a)') trim(adjustl(nstring))//'g'
    elseif (j.eq.5) then
      write(113,'(a)') trim(adjustl(nstring))//'h'
    elseif (j.eq.6) then
      write(113,'(a)') trim(adjustl(nstring))//'i'
    elseif (j.eq.7) then
      write(113,'(a)') trim(adjustl(nstring))//'k'
    elseif (j.eq.8) then
      write(113,'(a)') trim(adjustl(nstring))//'l'
    elseif (j.eq.9) then
      write(113,'(a)') trim(adjustl(nstring))//'m'
    elseif (j.eq.10) then
      write(113,'(a)') trim(adjustl(nstring))//'n'
    end if
  end do
end do

close(113)


end subroutine rcsfexcitation
