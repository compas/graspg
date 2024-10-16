program rmixaccumulate_csfg
  Use symexpand_mod

  Use rcsfsymexpand_I

  implicit none
! CYC
  integer, parameter :: nblevmax = 100
! CYC
  integer :: i,ii,j,jj,k,kk,l,ios,err,m
  integer :: nelec, ncftot, nw, nvectot, nvecsize, nblock, NDUM, ICSF
  integer :: nb, nevblk(200), iatjp, iaspa
  integer :: ivec
  integer :: nrelorb, indexans, co(200), ierr
  integer :: jr, jl, pos, norblayer, nlengthcore, nonsymcore 
  integer :: norbcomp, norb, nsymmetrymatch, norbmax, norbold, iadd, ncount
  integer, allocatable :: checkcsf(:,:), ind(:,:), ind2(:,:), ncfblk(:)
  integer, allocatable :: ncfblockout(:)

  double precision :: eav,eval
  double precision :: dc2evec, totc2evecblk, c2evecblklim
  double precision, allocatable :: evec(:,:,:),c2evec(:,:),c2evecsym(:,:)

  character(len=256) :: fnonsym
  character(len=200) :: state, file1, file2, file_rcsl, flsym, flsym2
  character(len=200) :: string2, string3
  character(len=600) :: header(5), string
  character(len=300) :: line1  ! The first line of CSF
  character(len=1500) :: string1, orbitalstring, orbitalstringnew, symorbset_Enn
  character(len=3) :: char3a, char3b
  character(len=1) :: char_l1, char_l2
  character(len=2) :: char_rl1, char_rl2
  character(len=5) :: relorbitals(200)
  character(len=6) :: G92MIX
  character(len=1) :: ciflag,sortflag,redeflag
  character(len=2) :: orbadd
  character(len=3) :: orbital(25)
  character(len=4) :: orbrel(300),orbcomp(300)
  character(len=200), allocatable  :: conf(:,:), coupling(:,:), spin(:,:)

  logical :: lsort

  write(*,*) '***************************************************************************'
  write(*,*) 'Welcome to program rmixaccumulate_csfg'
  write(*,*)
  write(*,*) 'The program accumulates the squared weights of the CSFGs from smaller CI.'
  write(*,*) 'The CSFGs in the output list can be sorted by mixing coefficents.'
  write(*,*) 'A priori condensation expansions can be obtained by redefining the'
  write(*,*) 'highest orbitals for each symmetry.'
  write(*,*)
  write(*,*) 'Input files: <state>.(c)m, <state>.g, and <state>.l'
  write(*,*) 'Output file: rcsfg.out, rcsfg.rem, rcsf.out, '
  write(*,*) '             rcsfg.nmax(optional)'
  write(*,*)
  write(*,*) ' Modified from rmixaccumulate.f90 written by J. Ekman & P. Jonsson Feb 2016'
  write(*,*) '                               CSFG-version written by Yanting Li, Apr 2022'
  write(*,*) '                                         Modify by  Chongyang Chen     2023'
  write(*,*) '***************************************************************************'

  write(*,*)
  write(*,*) 'Give name of the state: '
  read(*,*) state
  write(*,*) 'Expansion coefficients resulting from CI calculation (y/n)? '
  read(*,*) ciflag
  if(ciflag.ne.'y'.and.ciflag.ne.'n') then
     write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
     call exit()
  end if
  write(*,*) 'Fraction of total wave function [0-1] to be included in reduced list: '
  read(*,*) c2evecblklim
  write(*,*) 'Defining CSFs in output file sorted by mixing coefficients (y/n)? '
  read(*,*) sortflag
  if(sortflag.ne.'y'.and.sortflag.ne.'n') then
     write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
     call exit()
  end if
  lsort = .false.
  if (sortflag.eq.'y') lsort = .true.

!CYC  l = index(state,' ')
  l = len_trim(state)
  if(ciflag.eq.'y') then
     file1 = state(1:l)//'.cm'
  else
     file1 = state(1:l)//'.m'
  end if
  file2 = state(1:l)//'.c'
  flsym = state(1:l)//'.g'
  fnonsym = state(1:l)//'.l'
  call rcsfsymexpand(trim(flsym), trim(fnonsym), trim(file2), 0)
!
! Copy <name>.l to rlabel.inp, it will be used in last expansions.
!
  ierr = system ('cp ' // trim(fnonsym) // '  rlabel.inp')

  NDUM = MAXVAL(TotCSFs_perblock(2,1:NBLOCKSYM))

  ! Open and read data from mix file <state>.(c)m
  open (7,file=file1,status='old',form='unformatted')
  read (7,iostat=ios) G92MIX
  !write(*,*) G92MIX
  read (7) nelec, ncftot, nw, nvectot, nvecsize, nblock
  !write(*,*) '  nelec   = ', nelec
  !write(*,*) '  ncftot  = ', ncftot
  !write(*,*) '  nw      = ', nw
  !write(*,*) '  nblock  = ', nblock
  !write(*,*)

  ! Allocate various arrays
  allocate( conf(nblock,ncftot) )
  allocate( coupling(nblock,ncftot) )
  allocate( spin(nblock,ncftot) )
  allocate( ncfblockout(nblock) )
  allocate( ncfblk(ncftot) )
  allocate( ind(nblock,ncftot), ind2(nblock,ncftot) )
  allocate( checkcsf(nblock,ncftot) )
  allocate( evec(nblock,nblevmax,ncftot), c2evec(nblock,ncftot) )

  ! Continue to read data from mixing file <state>.(c)m
  ! Mixing coefficients are stored in 3D array: evec(block,eig,csf)
  evec(:,:,:) = 0.d0
  write(*,*)
!PERJ  write(*,*) 'Block data read from mixing file'
!PERJ  write(*,*) '        block        ncf         nev        2j+1          parity'
  do i=1, nblock
     READ (7,end=98) nb, ncfblk(i), nevblk(i), iatjp, iaspa
     !ALLOCATE(MAPCS(ncfblk(i)))
!PERJ     write(*,*) nb, ncfblk(i), nevblk(i), iatjp, iaspa
!CYC-2023
     if (nevblk(i).gt.nblevmax) then
       write(*,'(A,i5,A,i5)')"Error, nevblk(i).gt.nblevmax, ", nevblk(i), ">", nblevmax
       write(*,*)"Pls enlarge enough nblevmax, re-compile, re-run ..."
       stop 
     endif
!CYC-2023 end
     if(nevblk(i).gt.0) then
!        DO j = nvecpat + 1, nvecpat + nevblk(i)
!           ! ivec(i)   = ivec(i) + ncfpat ! serial # of the state
!           iatjpo(j) = iatjp
!           iaspar(j) = iaspa
!        ENDDO

        read (7) (ivec, j = 1,nevblk(i))
        !write(*,*) ivec(nvecpat+j)

        read (7) eav, (eval, j = 1, nevblk(i))
        !write(*,*) eav

        read (7) ((evec(i,k,j), j = 1, ncfblk(i)), k = 1, nevblk(i))
        !write(*,*) evec(nvecsizpat + ncfpat+j + (k-1)*ncftot)
     end if
     !nvecpat = nvecpat + nevblk(i)
     !ncfpat = ncfpat + ncfblk(i)
     !nvecsizpat = nvecsizpat + nevblk(i)*ncftot
     !DEALLOCATE(MAPCS)
  end do
98 continue
  close(7)

  ! For each CSF, calculate the sum of square of expansion coefficents for eigenvalues
  ! Block divided
  c2evec(:,:) = 0.0
  !write(*,*)
  do i=1,nblock
     !do j=1,2
        do k=1,ncfblk(i)
           do j=1,nevblk(i)
              c2evec(i,k) = c2evec(i,k)+evec(i,j,k)**2.d0
           end do
           c2evec(i,k) = c2evec(i,k)/dble(nevblk(i))
        end do
     !end do
  end do

!---------------------------------------------------------------------------
!Yanting for CSFG map
! 
!c2evecsym  -- Accumulation of squared weights 
!              of all CSFs belonging to an angular template
!ICSF       -- Position in name.c
!MAP2(1,j,i)-- The number of CSFs expanded by angular template 'j'
!              and block 'i'
! 
!---------------------------------------------------------------------------

  ALLOCATE(c2evecsym(nblock,NDUM))
  c2evecsym(:,:) = 0
  do i = 1,nblock 
    !ALLOCATE(MAPCS(ncfblk(i)))
    ICSF = 0
    DO j = 1, TotCSFs_perblock(2,i)
       DO k = 1, MAP2(1,j,i)
          ICSF = ICSF + 1
          !MAPCS(ICSF) = j
          c2evecsym(i,j) = c2evecsym(i,j) + c2evec(i,ICSF)
       ENDDO
       !c2evecsym(i,j) = c2evecsym(i,j)/dble(MAP2(1,j,i))
    ENDDO
  enddo

  ! Sort square of expansions coefficients (c2evecsym) and index table (ind)
  ! Compute sum of square of expansion coefficients for each block
  ! and "flag" CSF:s that contribute to user defined fraction of total wave functions
  ! write(*,*)
  checkcsf(:,:) = 0
  ind2(:,:) = 0
  do i = 1, nblock

     if (TotCSFs_perblock(2,i).gt.1) then
        call HPSORT(TotCSFs_perblock(2,i),c2evecsym(i,1:TotCSFs_perblock(2,i)),ind(i,1:TotCSFs_perblock(2,i)))
     end if
!CYC-2023
     !do k=TotCSFs_perblock(2,i),TotCSFs_perblock(2,i)-nevblk(i)+1,-1
     !   kk = TotCSFs_perblock(2,i) + 1 - k
     !   ind2(i,kk) = ind(i,k)
     !   if((MAP2(1,ind2(i,kk),i) /= 1) .and. (totc2evecblk.le.c2evecblklim)) then
     !     if(sortflag.eq.'n') then
     !        checkcsf(i,ind2(i,kk)) = 1
     !     else if(sortflag.eq.'y') then
     !        checkcsf(i,kk) = 1
     !     end if
     !   endif
     !end do
!CYC-2023 end

     totc2evecblk = 0.d0
     do k=TotCSFs_perblock(2,i),1,-1
        totc2evecblk = totc2evecblk + c2evecsym(i,k)
!CYC-2023
        ! ind(i,k) ==> k: the ind(i,k)-th CSFG within <state>.g corresponds to the k-th c_i^2 one.
        ind2(i,ind(i,k)) = k
        if (MAP2(1,ind(i,k),i) == 1 .and. c2evecsym(i,k).gt.1.0d-20) then
          ! Labelling CSFs / single CSFG, included nomatter the value of totc2evecblk
          checkcsf(i,ind(i,k)) = 1 
        elseif (totc2evecblk .le.c2evecblklim) then 
          checkcsf(i,ind(i,k)) = 1
        endif
!        kk = TotCSFs_perblock(2,i) + 1 - k 
!        if ((MAP2(1,ind(i,k),i) == 1) .and. (sortflag.eq.'y')) then
!            checkcsf(i,kk) = 1
!!Keep the CSFs built from the canonical orbitals and Only remove angular templates
!        elseif((MAP2(1,ind(i,k),i) /= 1) .and. (totc2evecblk.le.c2evecblklim)) then
!           if(sortflag.eq.'n') then
!              checkcsf(i,ind2(i,kk)) = 1
!           else if(sortflag.eq.'y') then
!              checkcsf(i,kk) = 1
!           end if
!        end if
!CYC-2023 end
     end do

  end do

  ! Open and read data from input file <state>.g
  open (8,file=flsym,status='old',form='formatted')
  do j=1,5
     read(8,'(a)') header(j)
  end do

  string = header(4)
  ii = 1
  do
     if(string((5*ii-3):(5*ii)).eq.'   ') then
       exit
     endif
     relorbitals(ii) = string((5*ii-3):(5*ii))
     ii = ii +1
  end do
  nrelorb = ii - 1
!PERJ  write(*,*) 'nrelorb: ', nrelorb

  i = 1
  j = 1
  do
     read(8,'(a)',end=99) string
     if(string(2:2).eq.'*') then
        i = i + 1
        j = 1
        read(8,'(a)') conf(i,j)
     else
        conf(i,j) = string
     end if
     read(8,'(a)') coupling(i,j)
     read(8,'(a)') spin(i,j)
     j = j + 1
  end do
99 continue
  close(8)

  ! Open output file rcsfg.out and write reduced CSFG list
  ! Open output file rcsfg.rem and write removed CSFG list
  co(:) = 0
  ncfblockout(:) = 0
  open (10,file='rcsfg.out',status='unknown',form='formatted')
  open (17,file='rcsfg.rem',status='unknown',form='formatted')
  do j=1,5
     write(10,'(a)') trim(header(j))
     write(17,'(a)') trim(header(j))
  end do
  do i = 1, nblock
!CYC-2023
     do k = TotCSFs_perblock(2,i), 1, -1
       if (lsort) then
         ! output the CSFGs according to the weigths
         j = ind(i,k)  ! the original position in <state>.g
         l = k         ! Its c^2
       else
         j = TotCSFs_perblock(2,i) + 1 - k ! the original position in <state>.g
         l = ind2(i,j) ! Its c^2
       endif

       if (checkcsf(i,j).eq.1) then
         ncfblockout(i) = ncfblockout(i) + 1
         write(10,'(a)') trim(conf(i,j))
         write(10,'(a)') trim(coupling(i,j))
         write(10,'(a)') trim(spin(i,j))
         !write(10,*)c2evecsym(i,l)
!     do j=1, TotCSFs_perblock(2,i)
!        if(MAP2(1,j,i) == 1) then
!            ncfblockout(i) = ncfblockout(i) + 1
!            write(10,'(a)') trim(conf(i,j))
!            write(10,'(a)') trim(coupling(i,j))
!            write(10,'(a)') trim(spin(i,j))
!        else if ((checkcsf(i,j).eq.1).or.(TotCSFs_perblock(2,i).eq.1)) then ! PJ Treat also special case with 1 CSF
!           ncfblockout(i) = ncfblockout(i) + 1
!           if(sortflag.eq.'n') then
!              write(10,'(a)') trim(conf(i,j))
!
!              if(sum(co(1:nrelorb)).lt.nrelorb) then
!                 do ii=1,nrelorb
!                    indexans = index(conf(i,j),relorbitals(ii)(1:4))
!                    if(indexans > 0) then
!                       co(ii) = 1
!                    end if
!                 end do
!              end if
!
!              write(10,'(a)') trim(coupling(i,j))
!              write(10,'(a)') trim(spin(i,j))
!           else if(sortflag.eq.'y') then
!              write(10,'(a)') trim(conf(i,ind2(i,j)))
!
!              if(sum(co(1:nrelorb)).lt.nrelorb) then
!                 do ii=1,nrelorb
!                    indexans = index(conf(i,ind2(i,j)),relorbitals(ii)(1:4))
!                    if(indexans > 0) then
!                       co(ii) = 1
!                    end if
!                 end do
!              end if
!
!              write(10,'(a)') trim(coupling(i,ind2(i,j)))
!              write(10,'(a)') trim(spin(i,ind2(i,j)))
!           end if
       else
           write(17,'(a)') trim(conf(i,j))
           write(17,'(a)') trim(coupling(i,j))
           write(17,'(a)') trim(spin(i,j))
!           write(17,*) "c2evecsym = ", c2evecsym(i,j)
           write(17,*) "c2evecsym = ", c2evecsym(i,l)
!CYC-2023 end
           write(17,*) 
       end if
     end do
     if(i.lt.nblock) then
        write(10,'(a2)') ' *'
        write(17,'(a2)') ' *'
     endif
  end do
  close(10)
  close(17)


  ! Print information of reduced list
  write(*,*)
  write(*,*) 'Number of CSFGs written to rcsfg.out'
  write(*,*) '        block        NCSF(G)'
  do i = 1, nblock
     write(*,*) i, ncfblockout(i)
  end do

! Deallocate the arrays
  DEALLOCATE(c2evecsym)
  DEALLOCATE(conf)
  DEALLOCATE(coupling)
  DEALLOCATE(spin)
  DEALLOCATE(ncfblockout)
  DEALLOCATE(ind)
  DEALLOCATE(ind2)
  DEALLOCATE(checkcsf)
  DEALLOCATE(evec)

! Allocated in rcsfsymexpand
  DEALLOCATE(map1)
  DEALLOCATE(map2)

!-------Redefine the orbitals of all the CSFGs-------------------------
  write(*,*) 'Redefine the CSFGs (y/n)? '
  read(*,*) redeflag
  if(redeflag.ne.'y'.and.redeflag.ne.'n') then
     write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
     call exit()
  end if
  if (redeflag.eq.'y') then
    open (78,file='rcsfg.out', status = 'old',  form = 'formatted')
    open (79,file='rcsfg.nmax',status='unknown',form='formatted')
    read(78,'(a)') string1
    write(79,'(a)') trim(string1)
    read(78,'(a)') string1
    write(79,'(a)') trim(string1)
    nlengthcore = len_trim(string1) + 1
    if (nlengthcore.lt.5) nlengthcore = 0
    nonsymcore = nlengthcore/5
    read(78,'(a)') string1
    write(79,'(a)') trim(string1)
    read(78,'(a)') string1
    write(*,*) trim(string1)
991 write(*,*) 'Give set of active orbitals, as defined by the highest principal&
                 quantum number'
    write(*,*) 'per l-symmetry, in a comma delimited list in s,p,d etc order,&
                  e.g. 5s,4p,3d'
    read(*,'(a)') orbitalstring

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

    m = len_trim(string1)
    symorbset_Enn = string1((nonsym-nonsymcore)*5+1:m)//' '

!   Number of orbitals
    norb = (len_trim(header(4))+1)/5
    read(unit=header(4),fmt='(300(x,a4))') orbrel(1:300)

!    Find out and write the orbitals for this layer in relativistic notation

    norbcomp = nonsym-nonsymcore
    do i = nonsym-nonsymcore+1,norb
       nsymmetrymatch = 0
       do j = 1,norblayer
          if (orbrel(i)(3:3).eq.orbital(j)(3:3)) then
             nsymmetrymatch = 1
             norbcomp = norbcomp + 1
             orbcomp(norbcomp) = orbrel(i)
             if ( i == norb .or.((i .LT. norb) .AND. &
                (orbrel(i)(3:4) .NE. orbrel(i+1)(3:4)))) then
                read(orbital(j)(1:2),'(I2)') norbmax
                read(orbrel(i)(1:2),'(I2)') norbold
                if (norbmax .GT. norbold) then
                   do iadd = norbold+1 , norbmax
                     write(orbadd,'(I2)') iadd
                     norbcomp = norbcomp + 1
                     orbcomp(norbcomp) = orbadd // orbrel(i)(3:4)
                   enddo
                endif
             endif
          end if
       end do
       if (nsymmetrymatch.eq.0) then
          norbcomp = norbcomp + 1
          orbcomp(norbcomp) =  orbrel(i)
       end if
    end do
     
    orbitalstringnew = ''
    ncount = 0
    write(*,*)nonsym, nonsymcore
    do i = 1, norbcomp
       if (i .le. nonsym-nonsymcore) then
         orbitalstringnew = orbitalstringnew(1:ncount*5)//' '//orbrel(i)
         ncount = ncount + 1
       else
         orbitalstringnew = orbitalstringnew(1:ncount*5)//' '//orbcomp(i)
         ncount = ncount + 1
       endif
    enddo
    write(*,*)'Redefined symmetry-ordered set:'
    write(*,*)trim(orbitalstringnew)
    write(79,'(a)') trim(orbitalstringnew)

! string "CSF(s):"
    read(78,'(a)') string1
    write(79,'(a)') trim(string1)

! Loop over CSFG in rcsfg.out
    do
       read(78,'(a)',end=109) string1
       if (string1(2:2).eq.'*') then
          write(79,'(a)') trim(string1)
          read(78,'(a)') string1
       end if
       read(78,'(a)') string2
       read(78,'(a)') string3
      
       m = len_trim(string1)
       line1 = trim(string1)

       ! The rightmost subshell
       ! The " nl" of last subshell
       char3a = string1(m-7:m-5)

       ! "l" character
       char_l1 = string1(m-5:m-5)

       ! Relativistic notation of the orbital
       char_rl1 = string1(m-5:m-4)
       if (char_rl1(2:2).ne.'-') char_rl1(2:2) = '+'

       if (index(trim(symorbset_Enn),char3a).gt.0) then
       ! Correlation CSFG:

         ! For the second rightmost subshell
         char3b   = string1(m-16:m-14)
         char_l2 = string1(m-14:m-14)
         char_rl2= string1(m-14:m-13)
         if (char_rl2(2:2).ne.'-') char_rl2(2:2) = '+'
         if (index(trim(symorbset_Enn),char3b).eq.0) then
        ! The second rightmost orbital is labelling symmetry
        ! TYPE 1 / TYPE 4
          do i = 1,norblayer
           if (orbital(i)(3:3).eq.char3a(3:3)) then
             if (lgt(orbital(i)(1:2),char3a(1:2))) then
               line1(m-7:m-6)=orbital(i)(1:2)
             endif 
           endif 
          enddo
         else
         ! TYPE 2 / TYPE 3
          if (char_rl1.eq.char_rl2) then
            ! TYPE 3
            do i = 1,norblayer
             if (orbital(i)(3:3).eq.char_l1) then
               if (lgt(orbital(i)(1:2),char3a(1:2))) then
                 read(orbital(i)(1:2),'(I2)') norbmax
                 write(orbadd,'(I2)') norbmax-1
                 line1(m-7:m-6)=orbital(i)(1:2) 
                 line1(m-16:m-15) = orbadd(1:2)
               endif
             endif
            enddo
          else
          ! TYPE 2
            do i = 1,norblayer
             if (orbital(i)(3:3).eq.char_l1) then
               if (lgt(orbital(i)(1:2),char3a(1:2))) then
                 line1(m-7:m-6)=orbital(i)(1:2)
               endif
             endif
             ! the second rightmost subshell
             if (orbital(i)(3:3).eq.char_l2) then
               if (lgt(orbital(i)(1:2),char3a(1:2))) then
                 line1(m-16:m-15)=orbital(i)(1:2)
               endif
             endif
            enddo 
          endif
         endif

       endif
       write(79,'(a)') trim(line1)
       write(79,'(a)') trim(string2)
       write(79,'(a)') trim(string3)
       
    enddo
109 continue
    close(79)
    close(78)
    flsym2 = 'rcsfg.nmax' 
  else
    ! Not redefined
    flsym2 = 'rcsfg.out'
  endif
 
! Expanding CSFGs
  ierr = system('rm -f rcsfg.inp')
  ierr = system ('cp ' // trim(flsym2) // '  rcsfg.inp')
  ierr = system ('cp ' // trim(fnonsym) // ' rlabel.inp')
  call rcsfsymexpand('rcsfg.inp', 'rlabel.inp', 'rcsf.out', 1)
  ierr = system ('rm -f rcsfg.inp rlabel.inp')

  stop "Normal Exit"
end program rmixaccumulate_csfg

SUBROUTINE HPSORT(N,RA,IND)
  integer N,IND(N),L,IR,I,J
  double precision RA(N), RRA

  L=N/2+1
  IR=N
  do I=1,N
     IND(I) = I
  end do

  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
    L=L-1
    RRA=RA(L)
    IRRA = IND(L)  ! je
  else
    RRA=RA(IR)
    IRRA = IND(IR) ! je
    RA(IR)=RA(1)
    IND(IR)=IND(1)  !je
    IR=IR-1
    if(IR.eq.1)then
      RA(1)=RRA
      IND(1)=IRRA !je
      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
     if(J < IR)then
        if(RA(J) < RA(J+1))  J=J+1
     end if
     if(RRA < RA(J))then
        RA(I)=RA(J)
        IND(I)=IND(J) !je
        I=J; J=J+J
     else
        J=IR+1
     end if

     goto 20
  end if
  RA(I)=RRA
  IND(I)=IRRA
  goto 10
END
