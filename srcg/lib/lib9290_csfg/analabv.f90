!***********************************************************************
!                                                                      *
      SUBROUTINE ANALABV(IC, IR, IV)
!                                                                      *
! Determine the positions of symmetry-ordered-orbitals within LABV,    *
! which are obtained by using Dirac-Coulomb Hamiltonian.               *
!                                                                      *
! Written by Chongyang Chen, Fudan university, Shanghai,  2020         *
!     Last modification by CYC                        Nov 2023         *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE symmatrix_mod,  ONLY: NSAMESYMCR, NSYMCOL, NSYMCR, &
                                NSYMROW, LABV, NORBGEN, IPSYM
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER IC, IR, IV
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER I, J, K, IERROR
!-----------------------------------------------
!
      NSYMCR = 0
      NSYMCOL = 0
      NSYMROW = 0
      NSAMESYMCR = 0
      IERROR = 0
      
      DO I = 1, 4 
         IF (LABV(I) .GT. NORBGEN) THEN
            NSYMCR = NSYMCR + 1 
            IF (I .LE. 2 ) THEN
            ! Two-eletron electronic interaction 
            ! symmetry-ordered orbitals of column basis: IC-CSF
               NSYMCOL = NSYMCOL + 1 
            ELSE
            ! Two-eletron electronic interaction 
            ! symmetry-ordered orbitals of row    basis: IR-CSF
               NSYMROW = NSYMROW + 1
            ENDIF
            IPSYM(NSYMCR) = I
         ENDIF
      ENDDO

! Checking, Types 2 
       
!      IF (NSYMCR .EQ. 1) THEN
!         IF (IPSYM(1) .NE. 2 ) THEN 
!           WRITE(*,'(3HIC=,I7,2X,3HIR=,I7,2X,3HIV=,I3,A)')
!     :          IC, IR, IV, 
!     :          "  The symmetry-ordered-orb doesn't located at position 2 ...."
!           STOP 'Error in analabv ***1***'
!         ENDIF 
!      ENDIF
! Checking, Types 3, 4, and 5

!      write(*,*)'IC,IR,IV,NSYMCR,NSYMCOL,NSYMROW=',
!     &     IC,IR,IV,NSYMCR,NSYMCOL,NSYMROW

      RETURN
      END SUBROUTINE ANALABV
