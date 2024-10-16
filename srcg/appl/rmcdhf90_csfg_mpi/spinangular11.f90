!***********************************************************************
!                                                                      *
      SUBROUTINE SPINANGULAR11 (JASAV, JBSAV, IUNITF, NPOS, JCSFASAV, JCSFBSAV)
!                                                                      *
!   This routine controls the computation  and storage of the values   *
!   and all indices of the angular coefficients                        *
!                                                                      *
!                                       k                              *
!                   T  (ab)            V  (abcd)                       *
!                    rs                 rs                             *
!                                                                      *
!   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
!   c and d are orbital sequence numbers.  r and s are configuration   *
!   state function indices.                                            *
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, RKCO_GG,       *
!                        TNSRJJ.                                       *
!               [GENMCP]: FNDBEG, SETSDA, SORT.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 28 Sep 1993   *
!   Modified by C. Froese Fischer for block computation.               *
!   Modified by G. Gaigalas and J. Bieron for new spin-angular         *
!   integration                                        01 April 2012   *
!                                                                      *
!***********************************************************************
!                                                                      *
!   Taken from mcpmpi_gg.f90                                           *
!                                                                      *
!   Modify for Spin-Angular integration: (TYPE 1 and 1 normal CSFs)    *
!              <CSF JA      | H_DC | CSF JB       >                    *
!                                                                      *
!   Written by Chongyang Chen, Fudan University, Shanghai, Oct  2023   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,   ONLY:  NNNW, KEYORB
      USE csfg_memory_man
      USE MPI_C     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I
      USE rkco_GG_I
      USE RINTI_I
      USE labpack_I
      USE SLATER_I
!-----------------------------------------------
!   C O M M O N  B L O C K S
!-----------------------------------------------
      USE BUFFER_C,  ONLY: NVCOEF, LABEL, COEFF
      USE DEBUG_C,   ONLY: LDBPA
      USE MCP_C,     ONLY: KMAX, DIAG, LFORDR
      USE ORB_C     
  
      Use symmatrix_mod, KMAXTmp=>KMAX 
      Use csfg_tv_C

      IMPLICIT NONE
      EXTERNAL CORD
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  JASAV, JBSAV, IUNITF, NPOS, JA, JB, JCSFASAV, JCSFBSAV
      !INTEGER  LLISTV(0:KMAX)
      LOGICAL  LINCR
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: KEYSQ = KEYORB*KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KK,IA,IB,LAB,IC,JCSFA,JCSFB,  &
         ID, NSWAP, ISWAP, LAC, LBD, NTGI
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      REAL(DOUBLE) :: VCOEFF, TSHELL1, TCOEFF, TEGRAL, ENONSYM

      LOGICAL :: F0INT, LOSCAL, LRKCO 
      INTEGER :: NDIFF1,NDIFF2,NORBCOLL,NORBCOLU,NORBROWL,NORBROWU
      INTEGER :: IV, I, J, K, L, M, N

!-----------------------------------------------
!
      JA = JASAV +  NCFGPAST
      JB = JBSAV +  NCFGPAST
      JCSFA = JCSFASAV
      JCSFB = JCSFBSAV

! IF (.NOT. LCHM), Calculate Potentials:
! There is no T-contributions for diagonal matrixelement for labeling
! CSFs.
!
! Accumulate the one-body contributions, build the Da LAB for the
! potentials involving the T-coefficients. 
      IA = 0
      CALL read_TV(1, JASAV, JBSAV, IUNITF, 1, IA, IB, TCOEFF)
      IF (IA /= 0 .AND. ABS(TCOEFF) > CUTOFF) THEN
        IF (LABTVFRST) NTPT = NTPT + 1
        IF (LABTVFRST.OR.LCPOT) THEN
          CALL SETCOF_NDA_CSFG(1, 1, IA, IB, TCOEFF)
        ELSE
          CALL IABINT(IA, IB, TEGRAL)
          EMTBLOCK(1,1) = EMTBLOCK(1,1) + TCOEFF * TEGRAL 
        ENDIF
      ENDIF

!
! Read the V-coefficients, accumulate the contributions of two-body
! interactions,  build the LAB for calculations of X-pot, Y-pot,
! Lagrange multipliers.
!
      NVCOEF = 0
      CALL read_TV(2, JASAV, JBSAV, IUNITF, 1, IA, IB, TSHELL1)
      DO 2 IV = 1, NVCOEF
         VCOEFF = COEFF(IV)
         LABV(1:5) = LABEL(1:5,IV)
         IF (LABTVFRST) NTPV = NTPV + 1
         IF (LABTVFRST.OR.LCPOT) THEN
           CALL SETCOF_NXY_CSFG(1, 1, VCOEFF)
         ELSE
           CALL twobody_DC(TEGRAL)
           EMTBLOCK(1,1) = EMTBLOCK(1,1) + VCOEFF * TEGRAL 
         ENDIF
    2 CONTINUE

      IF (LCHM.AND.LTRANSFER) CALL TRANSFER_CSFG(JASAV+NCFGPAST, JBSAV+NCFGPAST)

      RETURN
      END SUBROUTINE SPINANGULAR11
