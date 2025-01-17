!***********************************************************************
!                                                                      *
      SUBROUTINE TALK(JA, JB, NU, IA, IB, IC, ID, ITYPE, COEF)
!                                                                      *
!   Print  coefficients  and  integral  parameters  if IBUG1 > 0 and   *
!   write to disk.                                                     *
!                                                                      *
!                                           Last update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:16:29   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: KEYORB
      USE BUFFER_C
      USE debug_C
      USE orb_C
      USE cons_C
!CYC:
      USE csfg_decide_C,        ONLY:  NDISCARDBR, LDISCARDBR, NDISBR, LDISBR

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JA, JB, NU, IA, IB, IC, ID, ITYPE
      REAL(DOUBLE), INTENT(IN) :: COEF
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   Print coefficient if requested
!
      IF (IBUG1 /= 0) WRITE (99, 300)JA,JB,NP(IA),NH(IA),NP(IB),NH(IB),&
         NP(IC),NH(IC),NP(ID),NH(ID),NU,ITYPE,COEF

! CYC: Discarding the Breit interaction according to the NDISCARDBR,
! LDISCARDBR, NDISBR, and LDISBR. 
! NOTE:
! The contributions involving only s and p orbitals are limited by only
! LDISCARDBR and LDISBR, whereas the others are limited also by
! NDISCARDBR and NDISBR. 
!
! The followings worked only for the original RCI_MPI code,  they are
! wrong for RCI_MPI_CSFG because the values of n would be changed within
! matrixblockXX.f90.
      !IF (NKL(IA).GT.1 .OR. NKL(IB).GT.1 .OR.
      !&
      !    NKL(IC).GT.1 .OR. NKL(ID).GT.1) THEN
      !  IF (NDISCARDBR.AND.(NP(ID).GT.NDISBR .OR. NP(IC).GT.NDISBR
      !  .OR.&
      !                      NP(IB).GT.NDISBR .OR. NP(IA).GT.NDISBR))
      !                      &
      !     RETURN
      !ENDIF

! These work for RCI_MPI_CSFG, as the l-values are not changed
! throughout matrixblockXX.f90
      IF (LDISCARDBR.AND.(NKL(ID).GT.LDISBR .OR. NKL(IC).GT.LDISBR .OR.&
                          NKL(IB).GT.LDISBR .OR. NKL(IA).GT.LDISBR))   &
         RETURN
!
!   Increment counter
!
      IF(DABS(COEF) > EPS) THEN
         NVCOEF = NVCOEF + 1
!
!   Ensure that arrays are of adequate size; reallocate if necessary
!
         IF (NVCOEF > NBDIM) CALL ALCBUF (2)
!
!   Store integral indices and coefficient in COMMON/BUFFER/
!
         LABEL(1,NVCOEF) = IA
         LABEL(2,NVCOEF) = IB
         LABEL(3,NVCOEF) = IC
         LABEL(4,NVCOEF) = ID
         LABEL(5,NVCOEF) = NU
         LABEL(6,NVCOEF) = ITYPE
         COEFF(NVCOEF)   = COEF
      END IF
!
      RETURN
!
  300 FORMAT(2(1X,1I2),4(1X,I2,A2),1X,1I2,1X,1I2,1X,1P,D19.12)
      RETURN
!
      END SUBROUTINE TALK
