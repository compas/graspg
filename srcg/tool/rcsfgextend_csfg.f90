!                                                                      *
! Extend the low-n orbital set to high-n ones.                         *
! Written by Chongyang Chen, Fudan University, Shanghai                *
!                                              Sep 2022                *
!                                                                      *
!                                                                      *
!     Last modifications by CYC, Aug 2024                              *
!                                                                      *
!                                                                      *
!   Input  files: <name>.l, <name>.g (in CSFG-format)                  *
!   Output files: <name?>.g (in CSFG-format)                           *
!                                                                      *
program rcsfextend_csfg
implicit none
integer :: i, j, k, m, n, nlayer, norb, norblayer, nsymmetrymatch, norbcomp, ncsf, nwrite, ncount
integer :: jr, jl, pos, ncsflist(50), nmin, nmax
character(len=200) :: string1, string2, string3, name
character(len=1500) :: orbitalstring
character(len=3) :: charn
character(len=3) :: orb(300),orbital(25),orbcomp(300)
character(len=4) :: orbrel(300),char4
character(len=200) :: orbitallayer,label(50)

!CYC====================================================================
character*3    :: char3a, char3b
logical        :: flagsym,yinf
integer        :: nonsym,nonsyminf,nonsymcore
!character*11   :: string11_nrorb    ! 'spdfghiklmn'
character*1    :: char_l1, char_l2
character*2    :: char_rl1, char_rl2
character*1500 :: symorbset_Enmin,symorbset_Enn,stringnonsym,peelorbsmall
character*300  :: line1  ! The first line of CSF
integer        :: iwrite

!string11_nrorb='spdfghiklmn'
!write(*,*)"p =", index(string11_nrorb, 'p')


!CYC====================================================================

write(*,*) 'RCSFEXTEND_CSFG'
write(*,*) 'Extend a list name.g of CSFGs formed by larger set of active orbitals. '
write(*,*) 'Orbital sets are specified by giving the highest principal quantum number'
write(*,*) 'per l-symmetry, in a comma delimited list in s,p,d etc order, e.g. 5s,4p,3d'
write(*,*) 'Input file:   name.g, name.l'
write(*,*) 'Output files: namel<abel1>.g, name<label2>.g ...'
write(*,*) 'Output files: name<label1>.l, name<label2>.l ...'
write(*,*)
write(*,*) 'Written by Chongyang Chen, Fudan, Sep 2022 ...'
write(*,*)

write(*,*) 'Name of state'
read(*,'(a)') name
! Check <name>.l
yinf = .false.
inquire(file=trim(name)//'.l', exist=yinf)
if (yinf) then
  open (22,file=trim(name)//'.l', status='old',form='formatted')
  read (22, *)nonsyminf
  read (22, '(a)')stringnonsym
  close (22)
else
  write(*,*)"Error!", trim(name)//'.l', " is not found ..."
  STOP
endif
! Open file

open(unit=36,file=trim(name)//'.g',status='old')

! Read five first line save line four containing the string of orbitals

read(36,'(a)') string1
read(36,'(a)') string1
i = len_trim(string1) + 1
if (i.lt.5) i = 0
nonsymcore = i / 5
! string "Peel subshells:"
read(36,'(a)') string1
! Peel subshells
read(36,'(a)') orbitalstring
peelorbsmall = trim(orbitalstring)
! string "CSF(s):"
read(36,'(a)') string1

! Number of orbitals
norb = (len_trim(orbitalstring)+1)/5

! Get individual orbital strings non-relativistic / relativistic notation
read(unit=orbitalstring,fmt='(300(x,a3,x))') orb(1:300)
read(unit=orbitalstring,fmt='(300(x,a4))') orbrel(1:300)

! nonsym, flagsym
if (yinf) then
  flagsym = .true.
  ! number of labelling orbitals locating within peer-string
  nonsym = nonsyminf - nonsymcore
else
  STOP "Sth Error!"
endif
! string for all the symmetry-ordered (symoblic) Dirac orbitals:
symorbset_Enmin=orbitalstring(nonsym*5+1:len_trim(orbitalstring))//' '

! output for check
!write(*,*)"Orbitals set: "
!write(*,'(A)')trim(orbitalstring)
write(*,'(A)')"Labelling symmetries:"
write(*,'(300(x,a3,x))')orb(1 : nonsym)
!write(*,'(300(x,a4))')orbrel(1 : nonsym)
!write(*,*)"Non-relativistic symmetry-ordered symmetries:"
!write(*,'(300(x,a3,x))')orb(nonsym+1 : norb)
write(*,*)"Relativistic symmetry-ordered symmetries:"
!write(*,'(300(x,a4))')orbrel(nonsym+1 : norb)
write(*,'(A)')trim(symorbset_Enmin)

! Loop over numberof layers

write(*,*) 'Number of orbital sets'
read(*,*) nlayer
do k = 1,nlayer
991 continue
   write(*,*) 'Orbital set',k
   write(*,*) 'Give set of active orbitals, as defined by the highest principal quantum number'
   write(*,*) 'per l-symmetry, in a comma delimited list in s,p,d etc order, e.g. 5s,4p,3d'
   read(*,'(a)') orbitalstring
   write(*,*) 'Give file label'
   read(*,'(a)') label(k)
   open(unit=48+k,file=trim(name)//trim(label(k))//'.g',status='unknown')

! Number of orbitals in orbital set

   jl=1
   jr=0
   i=1
   do while ( index(orbitalstring(jl:len_trim(orbitalstring)),',').gt.0 )
      pos = index(orbitalstring(jl:len_trim(orbitalstring)),',')
      jr = jr + pos
      if (len_trim(orbitalstring(jl:jr-1)).eq.3) then
         orbital(i) = orbitalstring(jl:jr-1)
      else if (len_trim(orbitalstring(jl:jr-1)).eq.2) then
         orbital(i) = " " // orbitalstring(jl:jr-1)
      else
         write(*,*) 'Orbitals should be given in comma delimited list, redo!'
         goto 991
      end if
      jl = jr + 1
      i = i +1
   end do

   jr = len_trim(orbitalstring)
   if (len_trim(orbitalstring(jl:jr)).eq.3) then
      orbital(i) = orbitalstring(jl:jr)
   else if (len_trim(orbitalstring(jl:jr)).eq.2) then
      orbital(i) = " " // orbitalstring(jl:jr)
   else
      write(*,*) 'Orbitals should be given in comma delimited list, redo!'
      goto 991
   end if

   norblayer = i

! For current orbital layer find out the additional orbitals
! Checking
   do i = 1,norb
      nsymmetrymatch = 0
      do j = 1,norblayer
         if (orb(i)(3:3).eq.orbital(j)(3:3)) then
            nsymmetrymatch = 1
            if (lgt(orb(i)(1:2),orbital(j)(1:2))) then
              write(*,*)"Error, larger orbitals is smaller than original set: "
              write(*,*)orb(i), orbital(j)
              stop 
            end if
         end if
      end do
      if (nsymmetrymatch.eq.0) then
         write(*,*)"Error, larger orbitals has not found orbital: ", orb(i)
         stop 
         !norbcomp = norbcomp + 1
         !orbcomp(norbcomp) =  orb(i)
      end if
   end do
! Add the needed Dirac orbitals (correlation orbitals)
   norbcomp = 0
   orbitalstring = peelorbsmall(1:nonsym*5) 
   jl = nonsym*5 + 1
   orbrel(norb+1) = 'abcd'
   do i = nonsym+1, norb
     if (orbrel(i)(3:4) /= orbrel(i+1)(3:4)) then
       jr = 5*i
       n = len_trim(orbitalstring)
       if (orbitalstring(n:n).eq.'-') then
         orbitalstring = trim(orbitalstring) // peelorbsmall(jl:jr)
       else
         orbitalstring = trim(orbitalstring) // ' ' // peelorbsmall(jl:jr)
       endif
       jl = jr + 1
       ! Adding
       read(orbrel(i)(1:2), '(i2)')nmin
       do j = 1, norblayer
         if (orb(i)(3:3).eq.orbital(j)(3:3)) then
           read(orbital(j)(1:2),'(i2)')nmax
           exit
         endif
       enddo
       do n = nmin + 1, nmax
         write(charn,'(i3)')n
         if (orbrel(i)(4:4).eq.'-') then
           orbitalstring = trim(orbitalstring) // charn //orbrel(i)(3:4)
         else
           orbitalstring = trim(orbitalstring) //' '// charn //orbrel(i)(3:4)
         endif
       enddo
     endif
   enddo
   write(*,*)"Peer Subshells for layer k", k
   write(*,'(A)')trim(orbitalstring)
   rewind(36)
   read(36,'(a)') string1
   write(48+k,'(a)') trim(string1)
   read(36,'(a)') string1
   write(48+k,'(a)') trim(string1)
   read(36,'(a)') string1
   write(48+k,'(a)') trim(string1)
   read(36,'(a)') string1
   write(48+k,'(a)') trim(orbitalstring)

   write(*,*)"Symmetry-ordered symmetry set for En"//trim(label(k))//":"
   m = len_trim(orbitalstring)
   symorbset_Enn=orbitalstring(nonsym*5+1:m)//' '
   do i = 1, (len_trim(symorbset_Enn)+1)/5
      if (symorbset_Enn(i*5:i*5).ne.'-') symorbset_Enn(i*5:i*5) = '+'
   enddo
   write(*,'(A)')trim(symorbset_Enn)

   ! string "CSF(s):"
   read(36,'(a)') string1
   write(48+k,'(a)') trim(string1)

   ncsf = 0
   do
      read(36,'(a)',end=99) string1
      if (string1(2:2).eq.'*') then
         write(48+k,'(a)') trim(string1)
         read(36,'(a)') string1
      end if
      read(36,'(a)') string2
      read(36,'(a)') string3

! CYC===================================================================
! Labelling orbital set: 1s 2s 2p 3s 3p 3d 4s 4p 4d 4f 5s 5p
!
! Small CSFG orbital set: 8s,8p,8d,8f,8g,8h,8i,8k
! 
! Large CSFG orbital set: 11s,11p,11d,10f,9g,9h,9i,8k
! Symmetry-ordered orbital set: [6-11]s, [6-11]p, [5-11]d, [5-10]f, [5-9]g, [6-9]h, [7,8,9]i, 8k

! Orbital sets:
! En8 : 8s,8p,8d,8f,8g,8h,8i,8k
! En9 : 9s,9p,9d,9f,9g,9h,9i,8k
! En10: 10s,10p,10d,10f,9g,9h,9i,8k
! En11: 11s,11p,11d,10f,9g,9h,9i,8k

! Extend En8 CSFG to obtain En11 CSFG expansion.
!
! A: Labelling CSF, output it without any changes
  !  TYPE 1 CSF kept.

! Symmetry-ordered CSF (CSFG):
! C: Change the n-number of the symmetry-ordered orbitals to produce the
  !    needed symmetry-ordered CSFs
  !    "11s,11p,11d,10f,9g,9h,9i,8k" <== "8s,8p,8d,8f,8g,8h,8i,8k"
  !    TYPE 2: "[core A] 11s ( 1)" <== " 8s ( 1)"
  !    TYPE 3: "[core B] 11s ( 1) 11p-( 1)" <== " 8s ( 1)  8p-( 1)"
  !    TYPE 5: "[core C] 11s ( 2)" <== " 8s ( 2)"
  !
  !    TYPE 4:  
  !      En11:  "10s ( 1) 11s ( 1)"  <== " 7s ( 1)  8s ( 1)"
! CYC===================================================================
      iwrite = 0
      m = len_trim(string1)
      line1 = trim(string1)

      ! The " nl" of last subshell
      char3a = string1(m-7:m-5)

      ! "l" character
      char_l1 = string1(m-5:m-5)

      ! Relativistic notation of the orbital
      char_rl1 = string1(m-5:m-4)
      if (char_rl1(2:2).ne.'-') char_rl1(2:2) = '+'

      if (index(trim(symorbset_Enmin),char3a).eq.0) then
      ! Labelling CSF (TYPE 1)
        iwrite = 1

      ! symmetry-ordered CSF (CSFG):
      elseif (index(trim(symorbset_Enn), char_l1).gt.0) then
        ! index(trim(symorbset_Enn), char_l1).gt.0 means there is char_l1
        ! symmetry within symorbset_Enn
        n = len_trim(symorbset_Enn)

        ! For the last second subshell
        char3b  = string1(m-16:m-14)
        char_l2 = string1(m-14:m-14)
        char_rl2= string1(m-14:m-13)
        if (char_rl2(2:2).ne.'-') char_rl2(2:2) = '+'

        if (index(trim(symorbset_Enmin),char3b).eq.0) then
        ! The last second orbital is labelling symmetry
        ! TYPE 2 / TYPE 5
        ! symorbset_Enn contains the char_l1 symmetry, output
          iwrite = 1
          ! Replacement (if needed) for the last subshell
          do i = n, 1, -1
            if (symorbset_Enn(i:i).eq.char_l1) then
              if (lgt(symorbset_Enn(i-2:i),char3a))             &
                line1(m-7:m-5)=symorbset_Enn(i-2:i)
              exit
            endif
          enddo

        ! TYPE 3 / TYPE 4
        elseif (index(trim(symorbset_Enn),char_l2).gt.0) then
          ! The last second orbital is symmetry-ordered symmetry (TYPE 3 / 4
          ! CSFGs) AND symorbset_Enn contains the char_l2
          ! symmetry Then DO ...
          if (char_rl1.eq.char_rl2) then
            ! TYPE 4
            do i = n, 1, -1
              if (symorbset_Enn(i-1:i).eq.char_rl1) then
                ! There should be at least TWO same symmetries
                if (symorbset_Enn(i-6:i-5).eq.char_rl1) then
                  ! There are at least TWO same char_rl1 symmetry within
                  ! symorbset_Enn, output
                  iwrite = 1
                  line1(m-7:m-6)   = symorbset_Enn(i-3:i-2)
                  line1(m-16:m-15) = symorbset_Enn(i-8:i-7)
                  exit
                else
                  Stop "Something Error!..."
                endif
              endif
            enddo
          else
            ! TYPE 3
            ! char_l1 and char_l2 symmetries, output
            iwrite = 1

            ! Replacement (if needed) for the last subshell
            do i = n, 1, -1
              if (symorbset_Enn(i:i).eq.char_l1) then
                if (lgt(symorbset_Enn(i-2:i),char3a))      &
                  line1(m-7:m-5)=symorbset_Enn(i-2:i)
                exit
              endif
            enddo

            ! Replacement (if needed) for the last second subshell
            do i = n, 1, -1
              if (symorbset_Enn(i:i).eq.char_l2) then
                if (lgt(symorbset_Enn(i-2:i),char3b))      &
                  line1(m-16:m-14)=symorbset_Enn(i-2:i)
                exit
              endif
            enddo
          endif
        endif
      endif

      ! Output the needed CSFs
      if (iwrite.eq.1) then
        write(48 + k, '(a)') trim(line1)
        write(48 + k, '(a)') trim(string2)
        write(48 + k, '(a)') trim(string3)
        ncsf = ncsf + 1
      !else
      !  write(*,'(A)')trim(string1)
      !  write(*,'(A)')trim(string2)
      endif
   end do

99 continue

   ncsflist(k) = ncsf
end do

write(*,*)

do k = 1,nlayer
   write(*,*) 'List ',trim(name)//trim(label(k))//'.g',' based on orbital set',k
   write(*,*) 'contains',ncsflist(k),'CSFGs'
   write(*,*)
end do

! The followings are modified from rcsfgblocksplit_csfg.f90
write(*,*)
do k = 1,nlayer
   string1 = 'cp '// trim(name)//'.l ' // trim(name)//trim(label(k))//'.l'
   j = system (trim(string1))
!cjb  use EXECUTE_COMMAND_LINE if function 'system' is not supported
!cjb  call execute_command_line (trim(string1), exitstat=j)
   write(*,*) 'Exit status of the name.l file copying for layer',k,'was',j
enddo

end program rcsfextend_csfg
