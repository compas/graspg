!***********************************************************************
!                                                                      *
      SUBROUTINE genintbreit1 (myid, nprocs, NB, j2max)
!                                                                      *
!  Input:                                                              *
!     myid, nprocs                                                     *
!  Output:                                                             *
!   N - Number of integrals                                            *
!   j2max - max of 2*J                                                 *
!                                                                      *
!       Generate the list of Breit Integrals of type 1                 *
!       that could arise from a set of orbitals                        *
!       of orbitals.                                                   *
!                                                                      *
!       This routine is similar to genintrk                            *
!                                                                      *
!     Written by Per Jonsson            October 2014                   *
!     Last modification by Chongyang Chen   Nov 2023                   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE memory_man
      USE bilst_C
      USE orb_C, ONLY: nkj, nw, NP, NKL
      USE kkstartbreit_C
      USE symmatrix_mod, ONLY: KSTARTBREIT1IA
      use csfg_decide_C, ONLY: NDISCARDBR, LDISCARDBR, &
                               NDISBR, LDISBR
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE brintf_I
      USE triangbreit1_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: myid, nprocs
      INTEGER             :: NB, j2max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: GEN
      INTEGER :: I, KEY, K, IA, IB, IC, ID
!-----------------------------------------------
!
      KEY = NW + 1
      KSTARTBREIT1(0) = 1
      KSTARTBREIT1IA(0) = 0
!
!   Find 2*JMAX; it should not be greater than PARAMETER KMAX
!
      J2MAX = NKJ(1)
      DO I = 2, NW
         IF (NKJ(I) .GT. J2MAX) J2MAX = NKJ(I)
      ENDDO

      IF (J2MAX .GT. KMAX) THEN
         STOP 'genintbreit1: KMAX too small'
      ENDIF
!
!   Make the breit integrals. Note that there are no symmetry relations that
!   can be utilized. See paper by Grasnt.
!   When GEN is false, sweep through to find dimension
!
      GEN = .FALSE.

  999 NB = 0
      DO K = 0, J2MAX
         DO IA = 1, NW
! CYC NOTE:
! The contributions involving only s and p orbitals are limited by only
! LDISCARDBR and LDISBR, whereas the others are limited also by
! NDISCARDBR and NDISBR. 
           IF (NKL(IA).GT.1.AND.NDISCARDBR.AND.NP(IA).GT.NDISBR) THEN
             IF (GEN) KSTARTBREIT1IA(K*NW+IA) = NB
             CYCLE
           ENDIF
           IF (LDISCARDBR .AND. NKL(IA).GT.LDISBR) THEN
             IF (GEN) KSTARTBREIT1IA(K*NW+IA) = NB
             CYCLE
           ENDIF
           DO IB = 1, NW
             IF (NKL(IB).GT.1.AND.NDISCARDBR.AND.NP(IB).GT.NDISBR) CYCLE
             IF (LDISCARDBR.AND.NKL(IB).GT.LDISBR) CYCLE
             DO IC = 1, NW
               IF (NKL(IC).GT.1.AND.NDISCARDBR.AND.NP(IC).GT.NDISBR) CYCLE
               IF (LDISCARDBR.AND.NKL(IC).GT.LDISBR) CYCLE
               DO ID = 1, NW
                 IF (NKL(ID).GT.1.AND.NDISCARDBR.AND.NP(ID).GT.NDISBR) CYCLE
                 IF (LDISCARDBR.AND.NKL(ID).GT.LDISBR) CYCLE

                 IF (TRIANGBREIT1(IA,IB,IC,ID,K)) THEN
                    NB = NB + 1
                    IF (GEN .AND.(MOD(NB,nprocs) .EQ. myid)) THEN
                       INDTP1(NB) = ((IA*KEY+IB)*KEY+IC)*KEY+ID
                       VALTP1(NB) = BRINTF (1,IA,IB,IC,ID,K)
                    ENDIF
                 ENDIF
               ENDDO
             ENDDO
           ENDDO
           IF (GEN) KSTARTBREIT1IA(K*NW+IA) = NB
         ENDDO
         KSTARTBREIT1(K+1) = NB + 1
      ENDDO
!
!     Allocate memory for integral book keeping
!
      IF (.NOT. GEN) THEN
         CALL ALLOC (INDTP1,NB,'INDT1','genintbreit1')
         CALL ALLOC (VALTP1,NB,'VALT1','genintbreit1')

! Initialization is necessary in the mpi version

         DO i = 1, NB
            INDTP1(i) = 0
            VALTP1(i) = 0.d0
         ENDDO

         IF (myid .EQ. 0)                                         &
            PRINT *, 'Computing',NB,' Breit integrals of type 1'

         GEN = .TRUE.
         GOTO 999
      ENDIF

      RETURN
      END SUBROUTINE genintbreit1
