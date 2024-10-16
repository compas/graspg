!***********************************************************************
!                                                                      *
      program rcsfgexpand_csfg
!                                                                      *
!     rcsfexpand_csfg                                                  *
!                                                                      *
!     This program reads the CSFG list rcsf.inp and adds               *
!     CSFs to obtain a full normal list for checking and validation    *
!                                                                      *
!     WRitten by Per Jonsson, Malmo University, Sweden                 *
!     May 2017                                                         *
!                                                                      *
!     Modified by Chongyang Chen, Fudan University, Shanghai, China    *
!     May 2020                                                         *
!                                                                      *
!     Last modifications by CYC, Dec 2023                              *
!                                                                      *
!   Input  files: rcsfg.inp, rlabel.inp                                *
!   Output files: rcsf.out (in normal CSF-format)                      *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas  July 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character*1500 :: string(5)
      character*300 :: line1,line2,line3,lineadd
      character*2 :: string1,string2,string3,string4,sym(21)
      integer :: i,j,nlength,n,keep,nc,m,iunit,norb(21),nmin_sym(21)
      integer :: npos1,npos2,npos3,npos4,ntime,n1,n2,nend
      integer :: maxcsf,ncount1,ncount2
! cyc: ncount0 for Type 1 CSFs (CSF consisting of by labelling orbitals)
      integer :: ncount0,np1,np2
      logical :: iflag_symorb,yinf 
      integer :: nonsyminf,nonsymcore
      character*300 stringnonsym, FNAME
      integer :: nlengthcore,ierr

      integer :: nonsym,nblock,ncsfsymperblock(50)
!-----------------------------------------------
! cyc
      write(*,*)
      write(*,*) 'RCSFGEXPAND_CSFG'
      write(*,*) 'This program expands a list of CSFGs'
      write(*,*) 'Input files: rcsfg.inp, rlabel.inp'
      write(*,*) 'Output file: rcsf.out'
      write(*,*)

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

! Reading the information for the labelling-orbitals
! nonsymorb.inf needed for very flexible orbital-set.
!
! CYC-2024 
      FNAME = 'rlabel.inp'
      OPEN (2023,FILE=TRIM(FNAME), STATUS='OLD', &
               FORM='FORMATTED',IOSTAT=IERR)
      IF (IERR == 0) THEN
        READ(2023, *) nonsyminf
        READ(2023, '(A)')stringnonsym
        CLOSE(2023)
      ELSE
        WRITE (*, *) TRIM(FNAME) // &
             ' not exists, rcsfgexpand_csfg Stop Error!'
        STOP
      ENDIF
! CYC-2023/11/28/ end
      
! Open files    
      open (19,file = 'rcsfg.inp', status = 'old', form = 'formatted')
      open (20,file = 'rcsf.out', status = 'unknown',form = 'formatted')
      
! Start by looping through the file to count the number of blocks and
! CSFs in each block
      do i = 1, 5
         read(19,'(a)') string(i)
      end do
! cyc, nlengthcore is used to determine the number of core-orbital. 
      nlengthcore = len_trim(string(2)) + 1
      if (nlengthcore.lt.5) nlengthcore = 0
      nonsymcore=nlengthcore/5
! cyc
      ncsfsymperblock = 0
      nblock = 1
      do 
        read(19,'(a)',end=9) line1
        if (line1(1:2).eq.' *') then
          nblock = nblock + 1 
          read(19,'(a)') line1
        end if
        read(19,'(a)') line2  
        read(19,'(a)') line3  
        ncsfsymperblock(nblock) = ncsfsymperblock(nblock) + 1
      end do

9     continue

!! Allaocate space for map
!      do i = 1,nblock
!         write(*,*) 'Block',i,'number of CSFs',ncsfsymperblock(i)
!      end do
!      maxcsf = maxval(ncsfsymperblock)
!      write(*,*) 'maxcsf',maxcsf
!      allocate(map(maxcsf,nblock))

! Rewind
      rewind(19)

! Read and analyze header
      do i = 1, 5
         read(19,'(a)') string(i)
         write(20,'(a)') trim(string(i))
      end do

      nlength = len_trim(string(4)) + 1
! 
! cyc=================================================================== 
! Use nonsym = nonsyminf - nonsymcore
      nonsym = nonsyminf - nonsymcore
!! cyc=================================================================== 
!      if (.not.iflag_symorb) then
!         write(*,*)"Waring, there is no symmetry-ordered orbitals...."
!         nonsym=nlength/5
!         !stop "No symmetry-ordered orbitals, exiting..."
!      endif
      write(*,'(a)')trim(string(2))
      write(*,*)'nonsymcore=',nonsymcore
      write(*,'(a)')trim(string(4))
      write(*,*) 'Number of labelling orbitals within string4 : '
      write(*,*) nonsym
      write(*,'(a)')string(4)(1:nonsym*5)
! cyc=================================================================== 

! Loop through the symmetry-ordered orbitals and deterine the number of orbitals for
! each symmetry
      do i = nonsym+1,nlength/5
        do j = 1,21
          if (string(4)(5*(i-1)+4:5*i).eq.sym(j)) then
            norb(j) = norb(j) + 1
            if (norb(j).eq.1) then
              read(string(4)(5*(i-1)+1:5*(i-1)+3),'(i3)')nmin_sym(j)
              !write(*,'(a5,2x,i4)')string(4)(5*(i-1)+1:5*i), nmin_sym(j)
            endif
          end if
        end do
      end do
! cyc=================================================================== 
!PERJ      write(*,*)'Symmetry-ordered orbitals: '
! cyc=================================================================== 
!PERJ      do i = 1,21
!PERJ        write(*,*) norb(i),sym(i),nmin_sym(i) 
!PERJ      end do

! Loop through the CSFs and determine the type
      nblock = 1
      ncount1 = 0
      ncount2 = 0
!cyc
      ncount0 = 0
!cyc
      do 
        read(19,'(a)',end=99) line1
        if (line1(1:2).eq.' *') then
!cyc
        write(*,*)                                                   &
        'block',nblock, ' NCSF(L) =', ncount0, ' NCSF(G) =', ncount1,&
        ' NCSF =', ncount0+ncount2
!cyc
          write(20,'(a)') trim(line1)
          read(19,'(a)') line1
          nblock = nblock + 1
          ncount0 = 0
          ncount1 = 0
          ncount2 = 0
        end if
        read(19,'(a)') line2  
        read(19,'(a)') line3  
!cyc
        ncount1 = ncount1 + 1

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
!cyc
!        write(*,*)'ncount1=',ncount1, ' ntime=',ntime,  
!     & ' npos1=', npos1, ' npos2=',npos2, ' npos3=',npos3, 
!     & ' npos4=', npos4, ' ',trim(line1)
!cyc

!  Write CSFs of type 1
        if (ntime.eq.0) then
          ncount0 = ncount0 + 1
          write(20,'(a)') trim(line1)
          write(20,'(a)') trim(line2)
          write(20,'(a)') trim(line3)
        end if
      
        if (ntime.eq.1) then  ! CSF of type 2 or 5

!  Check the number of symmetry-ordered orbitals with this symmetry
          do j = 1,21
            if (string(4)(5*(npos2-1)+4:5*npos2).eq.sym(j)) then
              !n = norb(j)
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
!              string(4)(5*(npos2-2)+1:5*(npos2-1)) 
            write(20,'(a)') trim(lineadd)
            write(20,'(a)') trim(line2)
            write(20,'(a)') trim(line3)
            ncount2 = ncount2 + 1
          end do
        end if

        if (ntime.eq.2) then ! Type 3,4
          lineadd = line1 
! cyc: the following statement has a bug for some ORBs, 
! such as "  1s   2s   2p-  2p |  3s  4s   3p-  4p-  3p   4p   3d-  4d-  3d   4d   4f-  4f" 
!          if (npos4-1.ne.npos2) then ! Two symmetries e.g. 7s7p
! changed as: 
          if (string(4)(5*(npos4-1)+4:5*npos4).ne.                   &
                string(4)(5*(npos2-1)+4:5*npos2)) then

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
            ! TYPE - 3 
            do i = 1,n1
              lineadd(9*(npos1-1)+1:9*npos1-4) =                     &
                string(4)(5*(npos2-(n1-i)-1)+1:5*(npos2-(n1-i))) 
              do j = 1,n2 
                lineadd(9*(npos3-1)+1:9*npos3-4) =                   &
                  string(4)(5*(npos4-(n2-j)-1)+1:5*(npos4-(n2-j))) 
                write(20,'(a)') trim(lineadd)
                write(20,'(a)') trim(line2)
                write(20,'(a)') trim(line3)
                ncount2 = ncount2 + 1
              end do
            end do
 
          else        ! Same symmetry e.g. 6s7s

!  Check the number of symmetry-ordered orbitals with this symmetry
            ! TPYE - 4
            do j = 1,21
              if (string(4)(5*(npos4-1)+4:5*npos4).eq.sym(j)) then
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
                write(20,'(a)') trim(lineadd)
                write(20,'(a)') trim(line2)
                write(20,'(a)') trim(line3)
                ncount2 = ncount2 + 1
              end do
            end do

          end if 
        end if
      end do

99    continue
!cyc
      write(*,*)                                                     &
      'block',nblock, ' NCSF(L) =', ncount0, ' NCSF(G) =', ncount1,  &
      ' NCSF =', ncount0+ncount2
!cyc
      stop 'Normal Exit ...'

      end program rcsfgexpand_csfg
