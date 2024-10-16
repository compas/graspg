!***********************************************************************
!                                                                      *
      SUBROUTINE BRINT1 (IA,IB,IC,ID,K,TEGRAL,ICC,IRR,IVV)
!                                                                      *
!   Returns integrals for the transverse photon interaction.           *
!                                                                      *
!   Written by Per Jonsson                   Octaober 2014             *
!   Modify  by G Gaigalas                         May 2021             *
!                                                                      *
!   Speed up the search by using KSTARTBREIT1IA array                  *
!   Modify  by Chongyang Chen                         2020             *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE bilst_C
      USE orb_C,           ONLY: NW
      USE kkstartbreit_C
      use symmatrix_mod,   ONLY: KSTARTBREIT1IA
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ia, ib, ic, id, k
      INTEGER, INTENT(IN) :: ICC, IRR, IVV
      REAL(DOUBLE), INTENT(out) :: tegral
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER ::  index, loc, ju, jl, jm, key
!-----------------------------------------------
!
      KEY = NW + 1
      INDEX = ((IA*KEY+IB)*KEY+IC)*KEY+ID
!
!CYC
      !JL = KSTARTBREIT1(K)
      !JU = KSTARTBREIT1(K+1) - 1
      JL = KSTARTBREIT1IA(K*NW+IA-1) + 1
      JU = KSTARTBREIT1IA(K*NW+IA)
!
      IF (INDEX < INDTP1(JL).OR.INDEX > INDTP1(JU)) THEN
        WRITE(*,*) 'Something wrong in brint1'
        STOP
      ENDIF
!
!   The index is within the range of the indices stored; search
!   for it in the list of indices
!
    1 IF (JU-JL .GT. 1) THEN
         JM = (JU+JL)/2
         IF (INDTP1(JM) > INDEX) THEN
            JU = JM
         ELSE
            JL = JM
         ENDIF
         GOTO 1
      ENDIF
!
!   The range is bracketed to the extent possible
!
      IF (INDEX .EQ. INDTP1(JU)) THEN
         LOC = JU
      ELSEIF (INDEX .EQ. INDTP1(JL)) THEN
         LOC = JL
      ELSE
         WRITE(*,*) K,IA,IB,IC,ID,INDEX,'IC=',ICC,' IR=',IRR,' IV=',IVV
         WRITE(*,*) 'Brint1 Integral not found'
         STOP
      ENDIF
!
!   Return the value of the integral
!   from storage

      TEGRAL = VALTP1(LOC)
!
      RETURN
      END SUBROUTINE BRINT1
