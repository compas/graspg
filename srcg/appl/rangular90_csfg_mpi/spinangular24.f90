!***********************************************************************
!                                                                      *
      SUBROUTINE SPINANGULAR24 (JASAV, JBSAV, IUNITF, NPOS)
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
!   Modify for Spin-Angular integration: (TYPE 2 and TYPE 4 CSFGs)     *
!              <[core A] nk | H_DC | [core B] (n'-1)k' n'k'>           *
!   [core] is built by labelling symmetries,                           *
!   nk, (n'-1)k', n'k' are kappa-ordered orbitals                      *
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

      LOGICAL :: F0INT, LOSCAL, LRKCO, FLAGL, FLAGU, LSYM1, LSYM3
      INTEGER :: NDIFF1,NDIFF2,NORBCOL,NORBROW
      INTEGER :: IV, I, J, K, L, M, N, IRW, ICW, MTYPE
!-----------------------------------------------
!
!   LINCR is .TRUE. if NPOS is to be incremented by 1; there
!   is always a diagonal element in each column
!
      IF (LTRANSPOSE) THEN
        JA = JBSAV
        JB = JASAV
      ELSE
        JA = JASAV
        JB = JBSAV
      ENDIF

      LOSCAL = .FALSE.
      LRKCO = .FALSE.

      IF (JB /= JA) THEN
         LINCR = .TRUE.
      ELSE
         NPOS = NPOS + 1
         LINCR = .FALSE.
      ENDIF

! CSFs generated by the JB-th CSFG, Row demension
      M = NTYPE(2,JB)
! CSFs generated by the JA-th CSFG, Column demension
      N = NTYPE(2,JA)
! Number of symmetry-ordered orbs of JA-CSFG
      NORBCOL = NTYPE(6,JA) - NTYPE(3,JA) + 1
! Number of symmetry-ordered orbs of JB-CSFG 
      NORBROW = NTYPE(4,JB) - NTYPE(3,JB) + 1
!
! To use the lowest CSFs pair by default, then if the JB CSFG has one of
! the TWO symmetry-ordered orbitals of JA CSFG, the IRW and ICW meet the
! requirement that they have the same orbital by default. 
!
      IRW = IRFICT
      ICW = ICFICT
      ! < IR [CoreA] 7s | ==> < IR [CoreA] 4s|
      CALL FICTIOUS_CSF(2, IRFICT,   JB, NTYPE(3,JB), NTYPE(4,JB), 0, 0)
      ! | IC [CoreB] 6s 7s> | ==> | IC [CoreB] 4s 5s>
      CALL FICTIOUS_CSF(4, ICFICT, JA, NTYPE(3,JA), NTYPE(4,JA), &
                                     NTYPE(3,JA)+1, NTYPE(6,JA))
!
! No onebody contributions.
!
      FLAGL = .FALSE.
      FLAGU = .FALSE.
      IF (NTYPE(3,JB) /= NTYPE(3,JA)) GOTO 101

      FLAGL = .TRUE.
      IF (NORBROW > 1) FLAGU = .TRUE.
      ! < IR [CoreA] 7s | ==> < IR [CoreA] 5s|
      CALL FICTIOUS_CSF(2, IRFICT+1, JB, NTYPE(3,JB)+1, NTYPE(4,JB), 0, 0)
!
! There is possible onebody contributions between TYPE 2 - 4 CSFGs.
!
! Onebody contributions.
!
      DO MTYPE = 1, 2
        LOSCAL = .FALSE.
!
! MTYPE = 1: < [Core A] K | [Core B] I J>, K = I
! MTYPE = 2: < [Core A] K | [Core B] I J>, K = J
! 
        IF (MTYPE == 1) THEN
          IRW = IRFICT
          ICW = ICFICT
        ELSE
          IF (.NOT. FLAGU) CYCLE
          IRW = IRFICT + 1
          ICW = ICFICT
        ENDIF     
          ! Compute T coefficients
         
        IA = 0
        IB = 0
        CALL ONESCALAR(ICW, IRW, IA, IB, TSHELL)
        CALL SAVETV(1, JASAV, JBSAV, IUNITF, MTYPE, IA, IB, TSHELL(1))
        IF (IA /= 0 .AND. IA /= IB .AND.       &
               ABS (TSHELL(1)) > CUTOFF) THEN
           IF (.NOT. LOSCAL) THEN
             LOSCAL = .TRUE.
             ! NONESCALAR = NONESCALAR + 1
             LTRANSFER = .TRUE.

             DO K = 1, NORBROW
               N = 0
               DO I = 1, NORBCOL-1
                 DO J = I + 1, NORBCOL
                   N = N + 1
                   IF (MTYPE == 1 .AND. K == I) LNonZero(K,N) = .TRUE.
                   IF (MTYPE == 2 .AND. K == J) LNonZero(K,N) = .TRUE.
                 ENDDO
               ENDDO
             ENDDO
           ENDIF
!
! TSHELL(1) should be used for the semi-diagonal matrix elements,
! for (FLAGL) K == I. 
!           
! TSHELL(1) should be used for the semi-diagonal matrix elements,
! for (FLAGU) K == J. The value of TSHELL(1) may differ from the above
! one.  
!           
           !NPOS = NPOS + 1
           !LINCR = .FALSE.
           !LLISTT = LLISTT + 1
           !LAB = MIN (IA,IB) * KEY + MAX (IA,IB)
            
        ENDIF
      ENDDO

101   CONTINUE
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs. NVCOEF is initialized here but changed
!   in /rkco/cor[d] via COMMON.
!
      DO MTYPE = 1, 2
        IF (MTYPE == 1) THEN
          IRW = IRFICT
          ICW = ICFICT
        ELSE
          IF (.NOT. FLAGU) CYCLE
          IRW = IRFICT + 1
          ICW = ICFICT
        ENDIF     

        LSYM1 = .FALSE.
        LSYM3 = .FALSE.
        LRKCO = .FALSE.
        NVCOEF = 0
        CALL RKCO_GG (ICW, IRW, CORD, 1, 1)
        CALL SAVETV(2, JASAV, JBSAV, IUNITF, MTYPE, IA, IB, TSHELL1)
        DO 2 IV = 1, NVCOEF
           VCOEFF = COEFF(IV)
           IF (ABS (VCOEFF) > CUTOFF) THEN
             !IA = LABEL(1,IV)
             !IB = LABEL(2,IV)
             !IC = LABEL(3,IV)
             !ID = LABEL(4,IV)
             !KK = LABEL(5,IV)
             NSYMCR = COUNT(LABEL(1:4,IV) > NORBGEN)
             !F0INT=(KK.EQ.0).AND.(IA.EQ.IC).AND.(IB.EQ.ID)
             !IF (.NOT. F0INT) THEN
               IF (.NOT.LRKCO) THEN
                 LRKCO = .TRUE.
                 ! NRKCO = NRKCO + 1
                 LTRANSFER = .TRUE.
               ENDIF

               ! The following loops should be further checked by using
               ! ANALABV, if (NSYMCR .EQ. 3) then no modifications
               ! needed.
               IF (NSYMCR == 1 .AND. (.NOT. LSYM1)) THEN
                 LSYM1 = .TRUE.
                 N = 0
                 DO I = 1, NORBCOL - 1
                   DO J = I + 1, NORBCOL
                     N = N + 1   
                     DO K = 1, NORBROW
                       IF (MTYPE == 1 .AND. K == I) LNonZero(K,N)=.TRUE.
                       IF (MTYPE == 2 .AND. K == J) LNonZero(K,N)=.TRUE.
                     ENDDO
                   ENDDO 
                 ENDDO
               ELSEIF (NSYMCR == 3 .AND. (.NOT. LSYM3)) THEN
                 ! Full sub-matrix needed
                 LSYM3 = .TRUE.
                 LNonZero(1:NTYPE(2,JB),1:NTYPE(2,JA)) = .TRUE. 
               ENDIF

! Swap index to make sure IA <= IC, IB <= ID and record the number
! of swaps
             !  NSWAP = 0
             !  IF (IA > IC) THEN
             !     ISWAP = IC
             !     IC = IA
             !     IA = ISWAP
             !     NSWAP = NSWAP + 1
             !  ENDIF
             !  IF (IB > ID) THEN
             !     ISWAP = ID
             !     ID = IB
             !     IB = ISWAP
             !     NSWAP = NSWAP + 1
             !  ENDIF
             !  IF (LINCR) THEN
             !     NPOS = NPOS + 1
             !     LINCR = .FALSE.
             !  ENDIF
             !  LLISTV(KK) = LLISTV(KK) + 1
             !  LAC = IA * KEY + IC
             !  LBD = IB * KEY + ID
             !  IF (LAC .LT. LBD) THEN
             !     LAB = LAC * KEYSQ + LBD
             !  ELSE
             !     LAB = LBD * KEYSQ + LAC
             !  ENDIF
             !  WRITE (32+KK) JA, NPOS, LAB, VCOEFF
             !ENDIF
           ENDIF
    2   CONTINUE
!
!   All angular coefficients for this pair of CSFs have been
!   generated; update file 30
!

!        IF ((JB .EQ. JA) .OR. (.NOT. LINCR))   &
 
      ENDDO

      IF (LTRANSFER) THEN
        IF (LTRANSPOSE) LNonZero = TRANSPOSE(LNonZero)
        CALL TRANSFER_CSFG(JASAV, JBSAV)
      ENDIF

      RETURN
      END SUBROUTINE SPINANGULAR24
