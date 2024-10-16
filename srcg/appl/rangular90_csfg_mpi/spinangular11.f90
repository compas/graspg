!***********************************************************************
!                                                                      *
      SUBROUTINE SPINANGULAR11 (JASAV, JBSAV, IUNITF, NPOS)
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
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I
      USE rkco_GG_I
!-----------------------------------------------
!   C O M M O N  B L O C K S
!-----------------------------------------------
      USE BUFFER_C,  ONLY: NVCOEF, LABEL, COEFF
      USE DEBUG_C,   ONLY: LDBPA
      USE MCP_C,     ONLY: KMAX, DIAG, LFORDR

      Use symmatrix_mod, KMAXTmp=>KMAX 
      Use csfg_tv_C

      IMPLICIT NONE
      EXTERNAL CORD
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  JASAV, JBSAV, IUNITF, NPOS, JA, JB
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
      INTEGER :: KK,IA,IB,LAB,IC, &
         ID, NSWAP, ISWAP, LAC, LBD, NTGI
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      REAL(DOUBLE) :: VCOEFF, TSHELL1

      LOGICAL :: F0INT, LOSCAL, LRKCO 
      INTEGER :: NDIFF1,NDIFF2,NORBCOLL,NORBCOLU,NORBROWL,NORBROWU
      INTEGER :: IV, I, J, K, L, M, N

!-----------------------------------------------
!
!   LINCR is .TRUE. if NPOS is to be incremented by 1; there
!   is always a diagonal element in each column
!
      JA = JASAV
      JB = JBSAV
      LOSCAL = .FALSE.
      LRKCO = .FALSE.

      IF (JB /= JA) THEN
         LINCR = .TRUE.
      ELSE
         NPOS = NPOS + 1
         LINCR = .FALSE.

         ! The diagonla matrix element is always nonzero.
         LTRANSFER = .TRUE. 
         LNonZero(1,1) = .TRUE.
      ENDIF

      IA = 0
      IB = 0
      IF (JB /= JA) THEN
        ! Compute T coefficients
        CALL ONESCALAR(JA,JB,IA,IB,TSHELL)
        CALL SAVETV(1, JASAV, JBSAV, IUNITF, 1, IA, IB, TSHELL(1))
        IF (IA /= 0 .AND. IA /= IB .AND.       &
               ABS (TSHELL(1)) > CUTOFF) THEN
           IF (.NOT. LOSCAL) THEN
             LOSCAL = .TRUE.
             ! NONESCALAR = NONESCALAR + 1
             LTRANSFER = .TRUE.
             LNonZero(1,1) = .TRUE.
           ENDIF

           !NPOS = NPOS + 1
           !LINCR = .FALSE.
           !LLISTT = LLISTT + 1
           !LAB = MIN (IA,IB) * KEY + MAX (IA,IB)
            
        ENDIF
      ENDIF
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs. NVCOEF is initialized here but changed
!   in /rkco/cor[d] via COMMON.
!
      NVCOEF = 0
      CALL RKCO_GG (JA, JB, CORD, 1, 1)
      CALL SAVETV(2, JASAV, JBSAV, IUNITF, 1, IA, IB, TSHELL1)
      DO 2 IV = 1, NVCOEF
         VCOEFF = COEFF(IV)
         IF (ABS (VCOEFF) > CUTOFF) THEN
           !IA = LABEL(1,IV)
           !IB = LABEL(2,IV)
           !IC = LABEL(3,IV)
           !ID = LABEL(4,IV)
           !KK = LABEL(5,IV)
           !F0INT=(KK.EQ.0).AND.(IA.EQ.IC).AND.(IB.EQ.ID)
           !IF (.NOT. F0INT) THEN
!
! There are V-COEFF contributing to Potentials and H matrixelements.
! Record this call of RKCO_GG.
             IF (.NOT. LRKCO) THEN
               LRKCO = .TRUE.
               ! NRKCO = NRKCO + 1
               LTRANSFER = .TRUE.
               LNonZero(1,1) = .TRUE.
             ENDIF

! Swap index to make sure IA <= IC, IB <= ID and record the number
! of swaps
             !NSWAP = 0
             !IF (IA > IC) THEN
             !   ISWAP = IC
             !   IC = IA
             !   IA = ISWAP
             !   NSWAP = NSWAP + 1
             !ENDIF
             !IF (IB > ID) THEN
             !   ISWAP = ID
             !   ID = IB
             !   IB = ISWAP
             !   NSWAP = NSWAP + 1
             !ENDIF
             !IF (LINCR) THEN
             !   NPOS = NPOS + 1
             !   LINCR = .FALSE.
             !ENDIF
             !LLISTV(KK) = LLISTV(KK) + 1
             !LAC = IA * KEY + IC
             !LBD = IB * KEY + ID
             !IF (LAC .LT. LBD) THEN
             !   LAB = LAC * KEYSQ + LBD
             !ELSE
             !   LAB = LBD * KEYSQ + LAC
             !ENDIF
             !WRITE (32+KK) JA, NPOS, LAB, VCOEFF
           !ENDIF
         ENDIF
    2 CONTINUE

!   All angular coefficients for this pair of CSFs have been
!   generated; update file 30

!      IF ((JB .EQ. JA) .OR. (.NOT. LINCR))   &
 

      IF (LTRANSFER) CALL TRANSFER_CSFG(JASAV, JBSAV)

      RETURN
      END SUBROUTINE SPINANGULAR11
