program rcsfgsplit_csfg

! Written by Per Jonsson, August 2014                                  *
!                                                                      *
! Modified from rcsfsplit.f90 for CSFG                                 *
!             by Chongyang Chen, Fudan University, China      Jan 2021 *
!                                                                      *
!...The CSFG version was modified by  Gediminas Gaigalas  04/28/2021   *
!                                                                      *
!   Last modifications by CYC, Aug 2024                                *
!                                                                      *
!   Input  files: rlabel.inp, <name>.g (in CSFG-format)                *
!   Output files: <name?>.g (in CSFG-format)                           *
! 
implicit none
integer :: i, j, k, m, n, nlayer, norb, norblayer, nsymmetrymatch, norbcomp, ncsf, nwrite, ncount
integer :: jr, jl, pos, ncsflist(50)
character(len=200) :: string1, string2, string3, name, fnonsym
character(len=1500) :: orbitalstring
character(len=3) :: orb(300),orbital(25),orbcomp(300)
character(len=4) :: orbrel(300)
character(len=200) :: orbitallayer,label(50)

!CYC====================================================================
character*3    :: char3a, char3b
logical        :: flagsym,yinf
integer        :: nonsym,nonsyminf,nonsymcore
!character*11   :: string11_nrorb    ! 'spdfghiklmn'
character*1    :: char_l1, char_l2
character*2    :: char_rl1, char_rl2
character*1500 :: symorbset_Enmax,symorbset_Enn,stringnonsym
character*300  :: line1  ! The first line of CSF
integer        :: iwrite, ierr
integer        :: system

!string11_nrorb='spdfghiklmn'
!write(*,*)"p =", index(string11_nrorb, 'p')

!CYC====================================================================

write(*,*) 'RCSFGSPLIT_CSFG'
write(*,*) 'Splits a list name.g of CSFGs into a number of lists with CSFGs that '
write(*,*) 'can be formed from different sets of active orbitals. '
write(*,*) 'Orbital sets are specified by giving the highest principal quantum number'
write(*,*) 'per l-symmetry, in a comma delimited list in s,p,d etc order, e.g. 5s,4p,3d'
write(*,*) 'Input file: <name>.g rlabel.inp'
write(*,*) 'Output files: <name><label1>.g, <name><label2>.g, ...'
write(*,*) '              <name><label1>.l, <name><label2>.l, ...'
write(*,*)

write(*,*) 'Name of state'
read(*,'(a)') name
! CYC-2024/06/21
! Always use rlabel.inp
  yinf = .true.
! CYC-2023/11/28/
! Reading the information for the Labelling-orbitals
      fnonsym = 'rlabel.inp' 
      OPEN (2023,FILE=TRIM(fnonsym), STATUS='OLD', &
               FORM='FORMATTED',IOSTAT=IERR)
      IF (IERR == 0) THEN
        READ(2023, *) nonsyminf
        READ(2023, '(A)')stringnonsym
        CLOSE(2023)
      ELSE
        WRITE (*, *) TRIM(fnonsym) // &
             ' not exists, rcsfgsplit_csfg Stop Error!'
        STOP
      ENDIF
! CYC-2023/11/28/ end

! Open file
open(unit=36,file=trim(name)//'.g',status='old', &
               FORM='FORMATTED',IOSTAT=IERR)

! Read five first line save line four containing the string of orbitals

read(36,'(a)') string1
read(36,'(a)') string1
i = len_trim(string1) + 1
if (i.lt.5) i = 0
nonsymcore = i / 5
read(36,'(a)') string1
read(36,'(a)') orbitalstring
read(36,'(a)') string1

! Number of orbitals

norb = (len_trim(orbitalstring)+1)/5

! Get individual orbital strings non-relativistic / relativistic notation

read(unit=orbitalstring,fmt='(300(x,a3,x))') orb(1:300)
read(unit=orbitalstring,fmt='(300(x,a4))') orbrel(1:300)

! nonsym, flagsym
! CYC 2024/08/01 : ! Always use rlabel.inp
!if (yinf) then
  flagsym = .true.
  nonsym = nonsyminf - nonsymcore
!else
!  ! In case of there is no symmetry-ordered-orbitals, initialize nonsym
!  nonsym = (len_trim(orbitalstring) + 1)/5
!  !nonsym = 0
!  flagsym = .false.
!  ! skip 1s = 2s
!  if (index(trim(orbitalstring), '1s') .gt. 0) then
!    n = 3
!  else
!    n = 2
!  endif
!  do i = n, norb
!     if (orbrel(i-1)(3:4).eq.orbrel(i)(3:4)) then
!       nonsym = i - 2
!       flagsym = .true.
!       exit
!     endif
!  enddo
!  !if (.not. flagsym) then
!  !  write(*,*)"Warning, it is not symmetry-ordered-CSFs file ..."
!  !  write(*,*)"Program exits without any spliting ..."
!  !  !stop "Warning, it is not symmetry-ordered-CSFs file ..."
!  !endif
!endif
! CYC 2024/08/01 END
symorbset_Enmax=orbitalstring(nonsym*5+1:len_trim(orbitalstring))//' '

! output for check
!write(*,*)"Orbitals set: "
!write(*,'(A)')trim(orbitalstring)
write(*,*)"Labelling symmetries:"
write(*,'(300(x,a3,x))')orb(1 : nonsym)
!write(*,'(300(x,a4))')orbrel(1 : nonsym)
!write(*,*)"Non-relativistic symmetry-ordered symmetries:"
!write(*,'(300(x,a3,x))')orb(nonsym+1 : norb)
write(*,*)"Relativistic symmetry-ordered symmetries:"
!write(*,'(300(x,a4))')orbrel(nonsym+1 : norb)
write(*,'(A)')trim(symorbset_Enmax)

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
   open(unit=48+k,file=trim(name)//trim(label(k))//'.g',status='unknown', &
               FORM='FORMATTED',IOSTAT=IERR)

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

! For current orbital layer find out the compliment orbitals

   norbcomp = 0
   do i = 1,norb
      nsymmetrymatch = 0
      do j = 1,norblayer
         if (orb(i)(3:3).eq.orbital(j)(3:3)) then
            nsymmetrymatch = 1
            if (lgt(orb(i)(1:2),orbital(j)(1:2))) then
               norbcomp = norbcomp + 1
               orbcomp(norbcomp) =  orb(i)
            end if
         end if
      end do
      if (nsymmetrymatch.eq.0) then
         norbcomp = norbcomp + 1
         orbcomp(norbcomp) =  orb(i)
      end if
   end do

   rewind(36)
   read(36,'(a)') string1
   write(48+k,'(a)') trim(string1)
   read(36,'(a)') string1
   write(48+k,'(a)') trim(string1)
   read(36,'(a)') string1
   write(48+k,'(a)') trim(string1)
   read(36,'(a)') orbitalstring

!  Find out and write the orbitals for this layer in relativistic notation

   read(unit=orbitalstring,fmt='(300(x,a4))') orbrel(1:300)
   orbitalstring = ''
   ncount = 0
   do i = 1,norb
      nwrite = 0
      do j = 1,norbcomp
         if (orbrel(i)(1:3).eq.orbcomp(j)) then
            nwrite = 1
! CYC
            exit
! CYC
         end if
      end do

!  Write out the ones not in the complement orbital set

      if (nwrite.eq.0) then
         orbitalstring = orbitalstring(1:ncount*5)//' '//orbrel(i)
         ncount = ncount + 1
      end if
   end do
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

! "11s,11p,11d,10f,9g,9h,9i,8k"
! Symmetry-ordered orbital set: [6-11]s, [6-11]p, [5-11]d, [5-10]f, [5-9]g, [6-9]h, [7,8,9]i, 8k

! Orbital sets:
! En5 : 5s,5p,5d,5f,5g
! En6 : 6s,6p,6d,6f,6g,6h
! En7 : 7s,7p,7d,7f,7g,7h,7i
! En8 : 8s,8p,8d,8f,8g,8h,8i,8k
! En9 : 9s,9p,9d,9f,9g,9h,9i,8k
! En10: 10s,10p,10d,10f,9g,9h,9i,8k
! En11: 11s,11p,11d,10f,9g,9h,9i,8k
!
! A: Labelling CSF, output it without any changes
  !  TYPE 1 CSF kept.

! Symmetry-ordered CSF (CSFG):
! B: If there is any nonrelativistic orbital (such as 'k') which is not one of the
  !    orbitalstring (for En6, En7), discard it.
  !    "11s,11p,11d,10f,9g,9h,9i,8k" / "6s,6p,6d,6f,6g,6h"
  !    Any CSF with "i" and "k" symmetries should be discarded for En6.
  !
  !    Otherwise:
! C: Change the n-number of the symmetry-ordered orbitals to produce the
  !    needed CSFGs
  !    "11s,11p,11d,10f,9g,9h,9i,8k" / "6s,6p,6d,6f,6g,6h"
  !    TYPE 2: "[core A] 11s ( 1)" ==> " 6s ( 1)"
  !    TYPE 3: "[core B] 11s ( 1) 11p-( 1)" ==> " 6s ( 1)  6p-( 1)"
  !    TYPE 5: "[core C] 11s ( 2)" ==> " 6s ( 2)"
  !
  !    TYPE 4:  "10s ( 1) 11s ( 1)"  ==>
  !      En10:  "10s ( 1) 11s ( 1)"  ==> " 9s ( 1) 10s ( 1)"
  !      En9 :  "10s ( 1) 11s ( 1)"  ==> " 8s ( 1)  9s ( 1)"
  !      En8 :  "10s ( 1) 11s ( 1)"  ==> " 7s ( 1)  8s ( 1)"
  !      En7 :  "10s ( 1) 11s ( 1)"  ==> " 6s ( 1)  7s ( 1)"
  !      En6 :  5s is one labelling orbital, discard it,
  !             The CSF with " 5s ( 1)  6s ( 1)" generated
  !             by " 5s ( 1) 11s ( 1)" ==> " 5s ( 1)  6s ( 1)" (TYPE 2 step)
  !             Of course, their cores are different, as there are different
  !             electrons (1 / 2 ) in the symmetry-ordered orbitals.
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

      if (index(trim(symorbset_Enmax),char3a).eq.0) then
      ! Labelling CSF (TYPE 1)
      ! In the calculations for the n = 2 levels of N-like ions,
      ! the MRs are enlarged by including the n = 3 and 4 complexes,
      ! hence, the CSFs with the n = 4 orbitals should be discarded
      ! for the En3 calculation
        if (index(trim(orbitalstring),char3a).eq.0) then
         iwrite = 0
        else
         iwrite = 1
        endif

      ! symmetry-ordered CSF (CSFG):
      elseif (index(trim(symorbset_Enn), char_l1).gt.0) then
        ! index(trim(symorbset_Enn), char_l1).gt.0 means there is char_l1
        ! symmetry within symorbset_Enn
        n = len_trim(symorbset_Enn)

        ! For the last second subshell
        char3b   = string1(m-16:m-14)
        char_l2 = string1(m-14:m-14)
        char_rl2= string1(m-14:m-13)
        if (char_rl2(2:2).ne.'-') char_rl2(2:2) = '+'

        if (index(trim(symorbset_Enmax),char3b).eq.0) then
        ! The last second orbital is labelling symmetry
        ! TYPE 2 / TYPE 5
        ! symorbset_Enn contains the char_l1 symmetry (Line 313), output
          iwrite = 1
          ! Replacement (if needed) for the last subshell
          do i = n, 1, -1
            if (symorbset_Enn(i:i).eq.char_l1) then
              if (lgt(char3a,symorbset_Enn(i-2:i)))             &
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
                ! Is there another the same symmetry-ordered orbital ($l\pm$)?
                if (symorbset_Enn(i-6:i-5).eq.char_rl1) then
                  ! There are at least TWO same char_rl1 symmetry within
                  ! symorbset_Enn, output
                  iwrite = 1
                  ! Replacements (if needed) for the last two subshells
                  ! with the same char_rl1 symmetry.
                  if (lgt(char3a,symorbset_Enn(i-3:i-1))) then
                    line1(m-7:m-6)   = symorbset_Enn(i-3:i-2)
                    line1(m-16:m-15) = symorbset_Enn(i-8:i-7)
                  endif
                  exit
                endif
              endif
            enddo
          else
            ! TYPE 3
            ! Lines 328 and 354 ensure that symorbset_Enn contains both
            ! char_l1 and char_l2 symmetries, output
            iwrite = 1

            ! Replacement (if needed) for the last subshell
            do i = n, 1, -1
              if (symorbset_Enn(i:i).eq.char_l1) then
                if (lgt(char3a,symorbset_Enn(i-2:i)))      &
                  line1(m-7:m-5)=symorbset_Enn(i-2:i)
                exit
              endif
            enddo

            ! Replacement (if needed) for the last second subshell
            do i = n, 1, -1
              if (symorbset_Enn(i:i).eq.char_l2) then
                if (lgt(char3b,symorbset_Enn(i-2:i)))      &
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
   write(*,*) 'contains',ncsflist(k),'CSFs'
   write(*,*)
end do

! The followings are modified from rcsfgblocksplit_csfg.f90
write(*,*)
do k = 1,nlayer
   string1 = 'cp '//' rlabel.inp '// trim(name)//trim(label(k))//'.l'
   j = system (trim(string1))
!cjb  use EXECUTE_COMMAND_LINE if function 'system' is not supported
!cjb  call execute_command_line (trim(string1), exitstat=j)
   write(*,*) 'Exit status of the name.l file copying for layer',k,'was',j
enddo

end program rcsfgsplit_csfg
