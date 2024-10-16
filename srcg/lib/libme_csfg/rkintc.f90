!***********************************************************************
!                                                                      *
      SUBROUTINE RKINTC(IA0, IB0, IC0, ID0, K, TEGRAL)
!                                                                      *
!                         k                                            *
!   This routine returns R (abcd) integrals.                           *
!                                                                      *
!   Written by Per Jonsson                                             *
!   Modify  by G Gaigalas                         May 2021             *
!                                                                      *
!   Speed up the search by using KSTARTBREIT1IA array                  *
!   Modify  by Chongyang Chen                         2020             *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE cteilsrk_C
      USE orb_C
      USE kkstart_C
      use symmatrix_mod,   ONLY: KSTARTIA
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! CYC
!      INTEGER, INTENT(INOUT) :: IA0, IB0, IC0, ID0
      INTEGER, INTENT(IN) :: IA0, IB0, IC0, ID0
      INTEGER, INTENT(IN) :: K
      REAL(DOUBLE), INTENT(OUT) :: TEGRAL
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA, IB, IC, ID
      INTEGER :: KEY, ISWAP, INDEX, JL, JU, JM, LOC
      LOGICAL :: FOUND, FIRST
!-----------------------------------------------
!
      KEY = NW + 1
!
!   Ensure that the indices are in `canonical' order
!   Compute the composite (packed) index
!

!cychen
! CYC 2023/11/15: ensure that IA0, IB0, IC0 and ID0 won't be changed in
! this subroutine.
      IA=IA0
      IB=IB0
      IC=IC0
      ID=ID0

!cychen
      IF (IA > IC) THEN
         ISWAP = IC
         IC = IA
         IA = ISWAP
      ENDIF
      IF (IB > ID) THEN
         ISWAP = ID
         ID = IB
         IB = ISWAP
      ENDIF
      IF (IA > IB) THEN
         ISWAP = IB
         IB = IA
         IA = ISWAP
         ISWAP = ID
         ID = IC
         IC = ISWAP
      ENDIF

! CYC:
      IF (IA == IB .AND. IC > ID) THEN
        ISWAP = IC
        IC = ID
        ID = ISWAP
      ENDIF
     
!      if (IA0 == 2 .AND. IB0 == 3 .AND. IC0 == 4 .AND. ID0 == 2) THEN
!        WRITE(2095, *)'rintk : IA0, IB0, IC0, ID0 =', IA0, IB0, IC0, ID0
!        WRITE(2095, *)'rintk : IA,   IB,  IC,  ID =', IA, IB, IC, ID
!      ENDIF

      INDEX = ((IA*KEY + IB)*KEY + IC)*KEY + ID
!
! CYC
      !JL = KSTART(K)
      !JU = KSTART(K+1) - 1
      JL = KSTARTIA(K*NW+IA-1) + 1
      JU = KSTARTIA(K*NW+IA) 
      IF (INDEX<INDTEIRK(JL) .OR. INDEX>INDTEIRK(JU)) THEN
         WRITE (*, *)'JL = ', JL, '  JU=',JU
         WRITE (*, *) 'INDEX = ', INDEX
         WRITE (*, *)'INDTEIRK(JL)=', INDTEIRK(JL)
         WRITE (*, *)'INDTEIRK(JU)=', INDTEIRK(JU)
         WRITE (*, *) 'Something wrong in rkintc'
         STOP
      ENDIF
!
!   The index is within the range of the indices stored; search
!   for it in the list of indices
!
    1 CONTINUE
      IF (JU - JL > 1) THEN
         JM = (JU + JL)/2
         IF (INDTEIRK(JM) > INDEX) THEN
            JU = JM
         ELSE
            JL = JM
         ENDIF
         GO TO 1
      ENDIF
!
!   The range is bracketed to the extent possible
!
      IF (INDEX == INDTEIRK(JU)) THEN
         LOC = JU
      ELSE IF (INDEX == INDTEIRK(JL)) THEN
         LOC = JL
      ELSE
         WRITE (*, *) K, IA, IB, IC, ID, INDEX
         WRITE (*, *)'JL = ', JL, '  JU=',JU
         WRITE (*, *) 'INDEX = ', INDEX
         WRITE (*, *)'INDTEIRK(JL)=', INDTEIRK(JL)
         WRITE (*, *)'INDTEIRK(JU)=', INDTEIRK(JU)
         WRITE (*, *) 'Something wrong in rkintc'
         WRITE (*, *) 'Rkintc Integral not found'
         STOP
      ENDIF
!
!   Return the value of the integral
!   from storage

      TEGRAL = VALTEIRK(LOC)
!      if (IA0 == 2 .AND. IB0 == 3 .AND. IC0 == 4 .AND. ID0 == 2) THEN
!        WRITE(2095, *)'rintk : IA, IB, IC, ID, INDEX, LOC =', &
!                      IA, IB, IC, ID, INDEX, LOC, TEGRAL
!      ENDIF
!
      RETURN
      END SUBROUTINE RKINTC
