!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      *
!     rcsfsymexpand                                                    *
!                                                                      *
!     This subroutine reads the CSFG list rcsfg.inp and adds            * 
!     CSFs to obtain a full normal list for checking and validation    *
!                                                                      *
!     Also obtained by this subroutine are the parameters (variables)  *
!     in symexpand_mod                                                 *
!                                                                      *
!     WRitten by Per Jonsson, Malmo University, Sweden                 *
!     May 2017                                                         *
!                                                                      *
!     Modify by Chongyang Chen, Fudan University, PRC                  *
!     May 2020                                                         *
!     Last Modification by CYC, Aug 2024                               *
!                                                                      *
!cyc====================================================================
!CYC-2023/11/27                                                        *
!     Moved into lib9290_csfg, called (possibly) by rci_mpi_csfg /     *
!     rci_mpi_csfg_block / rmixaccumulate_csfg, ...                    *
!                                                                      *
!     Add interface file rcsfsymexpand_I.f90                           *
!                                                                      *
!  Dummy arguments:                                                    *
!        <name>    : the string of state                               *
!        iflag = 1 : generate   <name>.c in normal CSF format.         *
!        iflag = 0 : not output <name>.c in normal CSF format.         *
!                                                                      *
!For convience, generate the full <name>.c to obtain the full CSFs list*
!Called by rci_mpi_csfg /rci_mpi_csfg_block: iflag = 1,                *
!         by others:                         iflag = 0                 *
!                                                                      *
!CYC-2023/11/27 end                                                    *
!cyc====================================================================

      subroutine rcsfsymexpand(NAME, FILE_RL, FILE_OUT, IFLAG)
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use symexpand_mod
      use symexpand_mod, nblock=>nblocksym
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setlaborb_I
      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(*), INTENT(IN)   ::  NAME
      CHARACTER(*), INTENT(IN)   ::  FILE_RL
      CHARACTER(*), INTENT(IN)   ::  FILE_OUT
      INTEGER,      INTENT(IN)   ::  IFLAG
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character*1500 :: string(5)
      character*300 :: line1,line2,line3,lineadd
      character*2 :: string1,string2,string3,string4,sym(21)
      integer :: i,j,nlength,n,keep,nc,m,iunit,norb(21)
      integer :: npos1,npos2,npos3,npos4,ntime,n1,n2,nmin_sym(21)
      integer :: maxcsf,ncount1,ncount2,nlengthc,nend,icsf
!cyc====================================================================
!count for Type 1 CSFs
      integer :: ncount0,ncount3 
      integer :: nonsymcore
      logical :: yout, lcsfg
      character*1500 :: stringnonsym
!cyc====================================================================
!-----------------------------------------------------------------------
      write(*,*)
      write(*,*) 'Now perform rcsfgexpand ...'
!      write(*,*)

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

      norb = 0

!CYC-2023/11/26
!   Set the number of labelling orbitals. 
      CALL SETLABORB(TRIM(FILE_RL))
      stringnonsym = trim(strlaborb)
!   Set output <name>.c
      if (iflag == 0) then
        yout = .false.
      else
        yout = .true.
      endif
!CYC-2023/11/26 end

! Open files
      i = len_trim(name)
      open(19,file=trim(name), status = 'old',form ='formatted')
      if (yout) then
        open(20, file   = trim(FILE_OUT), &
                 status ='unknown', form = 'formatted')
      endif

! Start by looping through the file to count the number of blocks and
! CSFs in each block

      do i = 1, 5
         read(19,'(a)') string(i)
      end do

      TotCSFs_perblock=0
      nblock = 1
      do 
        read(19,'(a)',end=9) line1
        if (line1(1:2).eq.' *') then
          nblock = nblock + 1 
          read(19,'(a)') line1
        end if
        read(19,'(a)') line2  
        read(19,'(a)') line3  
        TotCSFs_perblock(2,nblock) = TotCSFs_perblock(2,nblock) + 1
      end do
9     continue

! Allaocate space for map
!PERJ      do i = 1,nblock
         !write(*,*) 'Block',i,'number of CSFs',ncsfsymperblock(i)
!PERJ         write(*,*) 'Block',i,'number of CSFs',TotCSFs_perblock(2,i)
!PERJ      end do
      maxcsf = maxval(TotCSFs_perblock(2,:))
!PERJ      write(*,*) 'maxcsf',maxcsf
      allocate(map1(maxcsf,nblock))
      allocate(map2(4,maxcsf,nblock))
      map1 = 0
      map2 = 0

! Rewind
      rewind(19)

! Read and analyze header

      do i = 1, 5
         read(19,'(a)') string(i)
         if (yout) write(20,'(a)') trim(string(i))
      end do

      nlength = len_trim(string(4)) + 1
      nlengthc = len_trim(string(2)) + 1
      if (nlengthc.lt.5) nlengthc = 0
      nonsymcore = nlengthc / 5 
!CYC-2023/11/26
      ! number of labelling orbitals within Peer shells
      nonsym = nonsyminf - nonsymcore
!CYC-2023/11/26 end
!cyc====================================================================
!      if (.not.iflag) then
!         write(*,*)"Waring, there is no symmetry-ordered orbitals...."
!         stop "No symmetry-ordered orbitals, exiting..."
!      endif
!cyc====================================================================
!PERJ      write(*,*) 'Number of labelling orbitals'
!PERJ      write(*,*) 'nonsym = ', nonsym + nlengthc/5
      i=nonsym
      read(string(4)(5*(i-1)+1:5*i-2),*)nmaxgen
!PERJ      write(*,*)'nmaxgen=',nmaxgen

! Loop through the symmetry-ordered orbitals and deterine the number of orbitals for
! each symmetry
      nsym_orb(:,:)=0
      nmin_sym = 0 
      do i = nonsym+1,nlength/5
        do j = 1,21
          if (string(4)(5*(i-1)+4:5*i).eq.sym(j)) then
             norb(j) = norb(j) + 1
             if (norb(j).eq.1) then
               read(string(4)(5*(i-1)+1:5*(i-1)+3),'(i3)')nmin_sym(j)
             endif
!cyc====================================================================
             nsym_orb(j,1)=nsym_orb(j,1)+1
             if (nsym_orb(j,2).eq.0) nsym_orb(j,2)=i+nlengthc/5 
             nsym_orb(j,3)=i+nlengthc/5 
          end if
        end do
      end do
!PERJ      write(*,*)'Symmetry-ordered orbitals: '
!cyc====================================================================
!PERJ      do i = 1,21
!PERJ        write(*,*) nsym_orb(i,:),norb(i),sym(i),nmin_sym(i) 
!PERJ      end do

! Loop through the CSFs and determine the type
      icsf = 0
      nblock = 1
      ncount0 = 0  ! TYPE-1 (Labelling) CSF counter   
      ncount1 = 0  ! Labelling CSFs and Correlation CSFG counter
      ncount2 = 0  ! Total number normal CSF in the expanded <name>.c

! Counting the leading labeling CSFs, used by jj2lsj_csfg, CYC, Aug 2024
      lcsfg   = .false. 
      scfDF1  = 0
!cyc====================================================================
      do 
        read(19,'(a)',end=99) line1
        if (line1(1:2).eq.' *') then
!cyc====================================================================
          write(*,*)'block',nblock, ' NCSF(L) =', ncount0,            &
           ' NCSF(G) =', ncount1, ' NCSF =', ncount0+ncount2 
          if (map1(ncount1,nblock).ne.map2(3,ncount1,nblock)) then
            write(*,*)'1: ncount1, nblock, map1( , ), map2(3, )=',     &
            ncount1,nblock, map1(ncount1,nblock), map2(3,ncount1,nblock)
            stop &
             'Unexpected map1(ncount1,nblock).ne.map2(3,ncount1,nblock)'
          endif

          TotCSFs_perblock(1,nblock)=ncount0
          if (TotCSFs_perblock(2,nblock).ne.ncount1)                 &
             stop 'Error in rcsfsymexpand ===1=== ...'
          TotCSFs_perblock(3,nblock)=ncount0 + ncount2
          icsf  = 0
          lcsfg = .false.
!cyc====================================================================
          if (yout) write(20,'(a)')trim(line1)
          read(19,'(a)') line1
          nblock = nblock + 1
!cyc====================================================================
          ncount0 = 0
          ncount1 = 0
          ncount2 = 0
!          write(98,'(a)')'-1'
        end if
        read(19,'(a)') line2
        read(19,'(a)') line3
!cyc====================================================================
!Count of the original Labelling and symmetry-ordered-CSF is move here, from the
!bottom of this file.
        ncount1 = ncount1 + 1
        ! Total number CSFs expanded by the ncount1-th CSF in rcsfg.inp
        map2(1,ncount1,nblock)=0
        ! The index in the fully expanded CSF (XXX.c) of the first CSF 
        ! expanded by the ncount1-th CSFG in XXX.g
        map2(2,ncount1,nblock)=0
        ! The index in the fully expanded CSF (XXX.c) of the final CSF 
        ! expanded by the ncount1-th CSFG in XXX.g
        map2(3,ncount1,nblock)=0
        ! Flags for Labelling CSFs and Correlation CSFGs:
        ! 0: the former; 1: the latter (CSFG)
        map2(4,ncount1,nblock) = 0
!cyc====================================================================

! Loop through the symmetry-ordered orbitals

        npos1 = 0
        npos2 = 0
        npos3 = 0
        npos4 = 0
        ntime = 0
        do i = 1,len_trim(line1)/9
          do j = nonsym+1,nlength/5
            if (line1(9*(i-1)+1:9*i-4).eq.string(4)(5*(j-1)+1:5*j)) then
              if (ntime.eq.0) then
                npos1 = i
                npos2 = j 
                ntime = ntime + 1
              else if (ntime.eq.1) then  
                npos3 = i
                npos4 = j 
                ntime = ntime + 1
              end if
            end if
          end do
        end do 
!        write(*,*) npos1,npos2,npos3,npos4,ntime
!cyc====================================================================
!        write(*,*)'ncount1=',ncount1, ' ntime=',ntime,  
!     & ' npos1=', npos1, ' npos2=',npos2, ' npos3=',npos3, 
!     & ' npos4=', npos4, ' ',trim(line1)
!cyc====================================================================
!  Write CSFs of type 1
        if (ntime.eq.0) then
          if (yout) then
            write(20,'(a)') trim(line1)
            write(20,'(a)') trim(line2)
            write(20,'(a)') trim(line3)
          endif
!cyc====================================================================
          ncount0 = ncount0 + 1
          icsf = icsf + 1
          map2(1,ncount1,nblock) = 1
          map2(2,ncount1,nblock) = icsf
          map2(3,ncount1,nblock) = icsf
          map2(4,ncount1,nblock) = 0          
!cyc====================================================================
! Counting the leading labeling CSFs, used by jj2lsj_csfg, CYC, Aug 2024
          if (.not.lcsfg) scfDF1(nblock) = scfDF1(nblock) + 1
        else
          lcsfg = .true.
        end if
      
        if (ntime.eq.1) then  ! CSF of type 2 or 5
          map2(4,ncount1,nblock) = 1
!  Check the number of symmetry-ordered orbitals with this symmetry
          do j = 1,21
            if (string(4)(5*(npos2-1)+4:5*npos2).eq.sym(j)) then
              read(line1(9*(npos1-1)+1:9*(npos1-1)+3),'(i3)')nend
              n = nend - nmin_sym(j) + 1
              exit
            end if
          end do

!  Replacements 7s    --> 4s,5s,6s,7s           
!               7s(2) --> 4s(2),5s(2),6s(2),7s(2)            

!  Singe loop over symmetry-ordered orbitals to do the replacement

          lineadd = line1 
          do j = 1,n
            lineadd(9*(npos1-1)+1:9*npos1-4) =                       &
             string(4)(5*(npos2-(n-j)-1)+1:5*(npos2-(n-j))) 
            !7s-6s, 7s(2)-6s(2)
            ncount2 = ncount2 + 1
            icsf = icsf + 1
            if (j.eq.1) map2(2,ncount1,nblock) = icsf
            if (j.eq.n) then
              map2(3,ncount1,nblock) = icsf
              map2(1,ncount1,nblock) = n
            endif
 
            if (yout) then
              write(20,'(a)') trim(lineadd)
              write(20,'(a)') trim(line2)
              write(20,'(a)') trim(line3)
            endif
          end do
        end if

        if (ntime.eq.2) then ! Type 3,4
          map2(4,ncount1,nblock) = 1
          lineadd = line1 
!cyc====================================================================
!The following commented out if-statement has a bug for some ORBs, such as 
!"  1s   2s   2p-  2p |  3s  4s   3p-  4p-  3p   4p   3d-  4d-  3d   4d   4f-  4f" 
!          if (npos4-1.ne.npos2) then ! Two symmetries e.g. 7s7p
!changed as: 
          if(string(4)(5*(npos4-1)+4:5*npos4).ne.                    &
             string(4)(5*(npos2-1)+4:5*npos2)) then
!cyc====================================================================
!  Check the number of symmetry-ordered orbitals with these two symmetries
            do j = 1,21
              if (string(4)(5*(npos2-1)+4:5*npos2).eq.sym(j)) then
                !n1 = norb(j)
                read(line1(9*(npos1-1)+1:9*(npos1-1)+3),'(i3)')nend
                n1 = nend - nmin_sym(j) + 1
                exit
              end if
            end do

            do j = 1,21
              if (string(4)(5*(npos4-1)+4:5*npos4).eq.sym(j)) then
                !n2 = norb(j)
                read(line1(9*(npos3-1)+1:9*(npos3-1)+3),'(i3)')nend
                n2 = nend - nmin_sym(j) + 1
                exit
              end if
            end do

!  Replacements 7s7p --> 4s4p,4s5p,4s6p,4s7p,5s4p,5s5p,5s6p,5s7,          
!                        ..... 7s6p,7s7p
!  Double loop over symmetry-ordered orbitals to do the replacement
            lineadd = line1 
            do i = 1,n1
              lineadd(9*(npos1-1)+1:9*npos1-4) =                     &
                string(4)(5*(npos2-(n1-i)-1)+1:5*(npos2-(n1-i))) 
              do j = 1,n2 
                lineadd(9*(npos3-1)+1:9*npos3-4) =                   &
                 string(4)(5*(npos4-(n2-j)-1)+1:5*(npos4-(n2-j))) 
                if (yout) then 
                  write(20,'(a)') trim(lineadd)
                  write(20,'(a)') trim(line2)
                  write(20,'(a)') trim(line3)
                endif
                ncount2 = ncount2 + 1

                icsf = icsf + 1
                if (i.eq.1 .and. j.eq.1) map2(2,ncount1,nblock) = icsf
                if (i.eq.n1.and. j.eq.n2) then
                  map2(3,ncount1,nblock) = icsf
                  map2(1,ncount1,nblock) = icsf-map2(2,ncount1,nblock)+1
                endif
              end do
            end do
 
          else        ! Same symmetry e.g. 6s7s

!  Check the number of symmetry-ordered orbitals with this symmetry

            do j = 1,21
              if (string(4)(5*(npos2-1)+4:5*npos2).eq.sym(j)) then
                !n = norb(j)
                read(line1(9*(npos3-1)+1:9*(npos3-1)+3),'(i3)')nend
                n = nend - nmin_sym(j) + 1
                exit
              end if
            end do

!  Replacements 6s7s --> 4s5s,4s6s,4s7s,5s6s,5s7s,6s7s

!  Double loop over symmetry-ordered orbitals to do the replacement

            lineadd = line1 

!  Please note that npos4 is the postion of 7s and that npos2 that
!  of 6s. We do the distribution in relation to npos4           
            do i = 1,n-1
              lineadd(9*(npos1-1)+1:9*npos1-4) =                     &
                string(4)(5*(npos4-(n-i)-1)+1:5*(npos4-(n-i))) 
              do j = i+1,n 
                lineadd(9*(npos3-1)+1:9*npos3-4) =                   &
                 string(4)(5*(npos4-(n-j)-1)+1:5*(npos4-(n-j))) 
                if (yout) then
                  write(20,'(a)') trim(lineadd)
                  write(20,'(a)') trim(line2)
                  write(20,'(a)') trim(line3)
                endif
                ncount2 = ncount2 + 1

                icsf = icsf + 1
                if (i.eq.1 .and. j.eq.2) map2(2,ncount1,nblock) = icsf
                if (i.eq.n-1.and. j.eq.n) then
                  map2(3,ncount1,nblock) = icsf
                  map2(1,ncount1,nblock) = icsf-map2(2,ncount1,nblock)+1
                endif
              end do
            end do

          end if 
        end if
        map1(ncount1,nblock) = ncount0+ncount2
        if (map1(ncount1,nblock).ne.map2(3,ncount1,nblock)) then
          write(*,*)'2: ncount1, nblock, map1( , ), map2(3, )=',      &
      ncount1,nblock, map1(ncount1,nblock), map2(3,ncount1,nblock)
        stop 'Unexpected map1(ncount1,nblock).ne.map2(3,ncount1,nblock)'
        endif

!cyc====================================================================
!      write(98,*)'nblock=',nblock,' ncount1=',ncount1, ' ntime=',ntime, 
!     & ' ncount2=', ncount2, ' irow=', ncount0+ncount2, ' ', trim(line1)
!cyc====================================================================

      end do

99    continue
      close(19)
      if (yout) close(20)
!cyc====================================================================
      !Include the core orbitals, then it is used by findtype ...
!CYC-2023/11/26
      nonsym = nonsym + nonsymcore
      if (nonsym /= nonsyminf) then
        write(*,*)'Error, nonsym /= nonsyminf ...'
        STOP "Error, nonsym /= nonsyminf ..."
      endif
!CYC-2023/11/26 end
      TotCSFs_perblock(1,nblock)=ncount0
      if (TotCSFs_perblock(2,nblock).ne.ncount1)                     &
        stop 'Error in rcsfsymexpand ...'
      TotCSFs_perblock(3,nblock)=ncount0 + ncount2
      write(*,*)                                                     &
        'block',nblock,' NCSF(L) =',ncount0,' NCSF(G) =',ncount1,    &
        ' NCSF =', ncount0+ncount2
!cyc====================================================================

! This is for testing the module      
!      call testmodule

      return
      end subroutine rcsfsymexpand
