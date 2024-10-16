!***********************************************************************
!                                                                      *
      SUBROUTINE FINDTYPE(NCFSMALL,NCFGEN,NCFTOT)
!                                                                      *
!   Reads the array with occupation numbers and deterimines the type   *
!   of CSF and the ordernumber of the orbitals in ORB_SD-MR            *
!   In the examplification of the span we assume that ORB_GEN          *
!   is 1s,2s,2p-,2p                                                    *
!                                                                      *
!   Type 1: Ordinary CSF, with no orbital in ORB_SD-MR, e.g. 2s2p      *
!   Type 2: CSFG with one orbital in ORB_SD-MR, e.g. 2s7p              *
!           spanning 2s3p, 2s4p, ..., 2s7p                             *
!   Type 3: CSFG with two orbitals of different symmetry               *
!           in ORB_SD-MR, e.g. 7s7p                                    *
!           spanning 3s3p, 3s4p,..., 3s7p, 4s3p, 4s4p, ...,4s7p,       *
!           7s3p,7s4p,...,7s7p                                         *
!   Type 4: CSFG with two orbitals of the same symmetry                *
!           in ORB_SD-MR, e.g. 6s7s                                    *
!           spanning 3s4s,3s5s,..,3s7s,4s5s,4s6s,4s7s,5s6s,5s7s,6s7s   *
!   Type 5: CSFG with a doubly occupied orbital in ORB_SD_MR           *
!           e.g. 7s(2)                                                 *
!           spanning 3s(2),4s(2),...,7s(2)                             *
!                                                                      *
!   NCF    = number of CSFs in the CSFG list                           *
!   NCFGEN = number of CSFs of type 1 in the CSFG list                 *
!   NCFTOT = total number of CSFs spanned by the CSFs in the           *
!            CSFG list                                                 *
!                                                                      *
!   NTYPE(6,NCF)                                                       *
!                                                                      *
!   NTYPE(1,J) = type of CSF J                                         *
!   NTYPE(2,J) = number of spanned CSFs from CSFG J                    *
!   NTYPE(3,J) = lower position of first  symmetry-ordered orb in CSF J*
!   NTYPE(4,J) = upper position of first  symmetry-ordered orb in CSF J*
!   NTYPE(5,J) = lower position of second symmetry-ordered orb in CSF J*
!   NTYPE(6,J) = upper position of second symmetry-ordered orb in CSF J*
!                                                                      *
!   Written by Per Jonsson, Malmo University April 2017                *
!                                                                      *
!   Last modifications by CYC, Fudan university, Dec 2023              *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use symexpand_mod 
      use symmatrix_mod 
      USE parameter_def
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      USE orb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE iq_I
      USE indexsym_I
      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: NCFSMALL, NCFGEN, NCFTOT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(LEN=2) :: char2
      INTEGER          :: indx, nsym, I, J, K, N, NDIFF, NMAX
      INTEGER          :: IFLAGSYM
!-----------------------------------------------
!
      NTYPE = 0
      ncsfDF1 = 0

! Define the highest principal quantum number in ORB_GEN. 

!      WRITE(*,*) 'In findtype, give highest n in ORB_GEN'
!      READ(*,*) NMAX

! Work out the number of orbitals in ORB_GEN assuming that all
! orbitals with n <= nmax are included

! NMAX = 2 --> 1s,2s,2p-2p                          --> NORBGEN = 4      
! NMAX = 3 --> additional: 3s,3p-,3p,3d-,3d         --> NORBGEN = 9      
! NMAX = 4 --> additional: 4s,4p-,4p,4d-,4d,4f-,4f, --> NORBGEN = 16      
! etc 
!cyc====================================================================
!      if(nonsym.eq.4)then
!        NMAX = 2
!      else if(nonsym.eq.9)then
!        NMAX = 3
!      else
!       write(*,*)'stop in findtype.f'
!       stop
!      end if
!      SELECT CASE(NMAX)
!        CASE(1)
!          NORBGEN = 1
!        CASE(2)
!          NORBGEN = 4
!        CASE(3)
!          NORBGEN = 9
!        CASE(4)
!          NORBGEN = 16
!        CASE(5)
!          NORBGEN = 25
!        CASE(6)
!          NORBGEN = 36
!        CASE(7)
!          NORBGEN = 49
!        CASE(8)
!          NORBGEN = 64
!        CASE(9)
!          NORBGEN = 81
!      END SELECT    
!      if (nonsym.ne.1.and.nonsym.ne.4.and.nonsym.ne.9.and.        &
!         nonsym.ne.16.and.nonsym.ne.25.and.nonsym.ne.36.and.      &
!         nonsym.ne.49.and.nonsym.ne.64.and.nonsym.ne.81) then
!         write(*,*)'Error, Stop in findtype.f...'
!         stop 'Abnormal Exit...'
!      endif
      !NMAX=int(sqrt(nonsym))
      NMAX=nmaxgen
      NORBGEN=nonsym
!      write(*,*)'NMAX, NORBGEN, nw=', NMAX, NORBGEN, nw
!cyc====================================================================
      IFLAGSYM = 0
! Loop over CSFs (including added ones) for current block 
      DO J = 1,NCF
        N = 0

! Scan occupation numbers. Put upper orbital positions in element 4 and 6
! We only have to loop over the orbitals not in ORB_GEN
!        WRITE(*,*)'J=',J,' IQ=',(IQ(I,J),I=1,NW)

        DO I = NORBGEN+1,NW
          IF (IQ(I,J).EQ.2) THEN  ! Doubly occupied, We know that we have type 5
            N = 5                 ! Set N = 5 to flag this
            NTYPE(4,J) = I
            EXIT                  ! We are done with this CSF       
          ELSEIF (IQ(I,J).EQ.1) THEN
            N = N + 1
            NTYPE(2*N+2,J) = I    ! Element 4 for n = 1 and 6 for n = 2
            IF (N.EQ.2) EXIT      ! We are done with this CSF
          END IF
        END DO
        IF (N.NE.0) IFLAGSYM = 1
! Now compute the type and lower positions. The latter will be in element 3 and 5

        IF (N.EQ.0) THEN                
          NTYPE(1,J) = 1
          ! ncsfDF1 used in maneig for initialization (eig-pairs)
          IF (IFLAGSYM .EQ. 0) ncsfDF1 = ncsfDF1 + 1
        ELSEIF (N.EQ.1) THEN
          NTYPE(1,J) = 2
!cyc====================================================================
!          NTYPE(3,J) = NTYPE(4,J) - NP(NTYPE(4,J)) + NMAX + 1  
          indx = indexsym(NH(NTYPE(4,J)))
          NTYPE(3,J) = nsym_orb(indx,2)  
!cyc====================================================================
        ELSEIF (N.EQ.2) THEN      ! Two cases: equal or unequal symmetry
          IF (NH(NTYPE(4,J)).EQ.NH(NTYPE(6,J))) THEN  ! Equal
            NTYPE(1,J) = 4
!cyc====================================================================
!            NTYPE(3,J) = NTYPE(4,J) - NP(NTYPE(4,J)) + NMAX + 1 
!            NTYPE(5,J) = NTYPE(6,J) - NP(NTYPE(6,J)) + NMAX + 2
            indx = indexsym(NH(NTYPE(4,J)))
            NTYPE(3,J) = nsym_orb(indx,2) 
            NTYPE(5,J) = nsym_orb(indx,2) + 1
!cyc====================================================================
          ELSE                                        ! Unequal
            NTYPE(1,J) = 3
!cyc====================================================================
!            NTYPE(3,J) = NTYPE(4,J) - NP(NTYPE(4,J)) + NMAX + 1 
!            NTYPE(5,J) = NTYPE(6,J) - NP(NTYPE(6,J)) + NMAX + 1
            indx = indexsym(NH(NTYPE(4,J)))
            NTYPE(3,J) = nsym_orb(indx,2)
            indx=indexsym(NH(NTYPE(6,J)))
            NTYPE(5,J) = nsym_orb(indx,2)
!cyc====================================================================

          END IF                                        
        ELSEIF (N.EQ.5) THEN
          NTYPE(1,J) = 5
!cyc====================================================================
!          NTYPE(3,J) = NTYPE(4,J) - NP(NTYPE(4,J)) + NMAX + 1
          indx = indexsym(NH(NTYPE(4,J)))
          NTYPE(3,J) = nsym_orb(indx,2)
!cyc====================================================================
        END IF
      END DO        

! Count number of CSFs of type 1 as well as the total number of
! CSFs spanned by the symmetry-ordered CSFs      

      NCFTOT = 0
      NCFGEN = 0
      MAXSPAN = 0
      DO K = 1,NCFSMALL
        !J = MAP(K)
        ! CYC: Smallest list is used, everyone should be counted.
        J = K
        SELECT CASE (NTYPE(1,J))
          CASE(1)
            NDIFF = 1 
            NCFGEN = NCFGEN + NDIFF
          CASE(2)
            NDIFF = NTYPE(4,J) - NTYPE(3,J) + 1 
          CASE(3)  
            NDIFF = (NTYPE(4,J) - NTYPE(3,J) + 1)*                   &
                    (NTYPE(6,J) - NTYPE(5,J) + 1)       
          CASE(4) 
            NDIFF = ((NTYPE(6,J) - NTYPE(5,J) + 2)*                  &
                               (NTYPE(6,J) - NTYPE(5,J)+1 ))/2
          CASE(5)
            NDIFF = NTYPE(4,J) - NTYPE(3,J) + 1
        END SELECT
        NTYPE(2,J) = NDIFF
        NCFTOT = NCFTOT + NDIFF
        IF (NDIFF.GT.MAXSPAN) THEN
          MAXSPAN = NDIFF
        END IF

      END DO
!      IF(NTYPE(1,J) .EQ. 4)Then
!        NTYPE(5,J)=NTYPE(5,J)+1
!      END IF
!      WRITE(*,*)'ncsfDF1 =', ncsfDF1
      RETURN
      END SUBROUTINE FINDTYPE
