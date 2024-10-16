!***********************************************************************
!     ------------------------------------------------------------------
!     JJGEN  -- program to generate
!
!     Written by: Lennart Sturesson
!
!     2:nd version
!
!     Last edited Januar 2, 1997
!
!     ------------------------------------------------------------------
!
      subroutine jjgen15
!...The CSFG version was modified by  Gediminas Gaigalas  04/21/2021
!...Switches:
!
! Last modification for CSFG, Chongyang Chen, Fudan university, May 2021 
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use reffa_I
      use adder_I
      use matcin_I
      use blandc_I
      use matain_I
      use fivelines_I
      use blanda_I
      use merge_I
      use open79_I
      use matbin_I
      use copy7t9_I
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: logfil = 31
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(25,0:10) :: org
      integer :: varmax, skal, anel, par
      integer , dimension(25,0:10) :: low
      integer , dimension(220) :: posn, posl
      integer :: nmax
      integer , dimension(25) :: lim
      integer :: minj, maxj, cfmax, ii
      logical , dimension(25,0:10) :: lock, closed
      logical :: slut
      logical , dimension(25,0:10) :: med, dubbel
      logical :: advexp, second
      character :: x
!-----------------------------------------------
!cyc====================================================================
      integer            :: nmax_spec, lmax_spec,nmaxcan, nonsym_free
      integer            :: orb_spec_nl(0:25), ntime
      integer            :: nmax_orbset(100),nmrcon
      logical            :: freeorb
      common/orb_spec/nmax_spec,lmax_spec,nmaxcan,nonsym_free,        &
                      orb_spec_nl,nmax_orbset,nmrcon,freeorb,ntime
!cyc====================================================================

      open(unit=193,file='excitationdata',status='unknown')
      write(*,*) ' Excitationdata file opened'


      open(unit=logfil,file='rcsfg.log',status='unknown')
      write(*,*)
      write(*,*) 'THE LINES BELOW ARE DEBUG OUTPUT FROM THE'
      write(*,*) 'CSF GENERATOR:   PLEASE IGNORE !!!       '
      write(*,*)
   10 write(*,200) '    *  : new list'
!      write(*,200) '    a  : add to existing list'
      write(*,200) '    e  : expand existing list'
      write(*,200) '    q  : quit'
      read(193,100) X
      write(logfil,200) ' Option : ',X
      call Reffa(posn,posl)
      advexp = .FALSE.
! PJ GG
      closed = .false.
      if (X.EQ.'a' .OR. X.EQ.'A') then
         call Adder(closed,med,slut,anel,par,.FALSE.)
         if (slut) then
            write(*,200)
            write(*,200) 'The clist.inp-file is not readable! '
            stop
         endif
         write(logfil,200) ' New reference set.'
         second = .TRUE.
      elseif (X.EQ.'e' .OR. X.EQ.'E') then
         advexp = .TRUE.
         call Adder(closed,med,slut,anel,par,.TRUE.)
         if (slut) then
            write(*,200)
            write(*,200) 'The clist.inp-file is not readable! '
            stop
         endif
         write(logfil,200) ' File as reference sets.'
         call  &
            Matcin(lock,closed,med,varmax,cfmax,nmax,minJ,maxJ,lim)
         call  &
            Blandc(varmax,cfmax,lock,med,minJ,maxJ,nmax,posn,posl,lim)
         second = .FALSE.
      else
         call Matain(org,lock,closed,varmax,skal,nmax,anel,par, &
                                    low,minJ,maxJ,lim,dubbel)
         call Fivelines(org,lock,closed,.TRUE.,posn,posl)
         call Blanda(org,varmax,lock,minJ,maxJ,skal,nmax,low,   &
                                  posn,posl,lim,dubbel,.TRUE.)
         second = .FALSE.
      endif
           ii=0
      if(.not.second) then
         call Merge(.TRUE.,posn,posl,ii)
         if(advexp) ii=ii+1
         call open79(ii)
      endif

      do
         call Matbin(org,lock,closed,varmax,skal,second,anel,   &
                            par,low,nmax,lim,dubbel,minJ,maxJ)
         if(.not.second) exit
         call Fivelines(org,lock,closed,.FALSE.,posn,posl)
         call Blanda(org,varmax,lock,minJ,maxJ,skal,nmax,low,   &
                                posn,posl,lim,dubbel,.FALSE.)
         call Merge(.false.,posn,posl,ii)
         ii=ii+1
         call open79(ii)
         second = .false.
      enddo
!      write(*,200) 'The merged file is called rcsf.out.'
      if(mod(ii,2).eq.0.and.ii.ne.0) call copy7t9
      if(mod(ii,2).ne.0) then
        open(unit=93,file='fil1.dat',status='unknown')
        close(unit=93,status='delete')
      end if
      if (ii.eq.0) then
         close(unit=7,status='delete')
      end if

      close(unit=7)
      close(unit=9)

      ! Here all CSFs are kept by setting keep=1 in rcsfsym.
      !close(193,status='delete') ! PJ
      close(193) ! CYC: keep excitationdata file for check

      call rcsfsym
      call rcsfblock

      ! Here the unneeded CSFs are safely weeped out.
      call rcsfsym_cyc
      call stoptime (ntime, 'RCSFGGenerate_CSFG') 
      !stop "Normal Exit ......"
  100 format(A)
  200 format(' ',10A)
  300 format(' ',A,I2,A)
  400 format(' ',A,I3,A)
      end
