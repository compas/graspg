!***********************************************************************
!                                                                      *
      SUBROUTINE SPINANGULAR22 (JASAV, JBSAV, IUNITF, NPOS, JCSFASAV, JCSFBSAV)
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
!   Modify for Spin-Angular integration: (TYPE 2 and TYPE 2 CSFGs)     *
!              <[core A] nk | H_DC | [core B] n'k'>                    *
!   [core] is built by labelling symmetries,                           *
!   nk, n'k' are kappa-ordered orbitals                                *
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

      LOGICAL :: F0INT, LOSCAL, LRKCO, LDIAG, LOFFD
      INTEGER :: NDIFF1,NDIFF2,NORBCOL,NORBROW
      INTEGER :: IV, I, J, K, L, M, N, IRW, ICW
      LOGICAL :: FLAGU
!-----------------------------------------------
!
      IF (LTRANSPOSE) THEN
        JA = JBSAV +  NCFGPAST
        JB = JASAV +  NCFGPAST
        JCSFA = JCSFBSAV
        JCSFB = JCSFASAV
      ELSE
        JA = JASAV +  NCFGPAST
        JB = JBSAV +  NCFGPAST
        JCSFA = JCSFASAV
        JCSFB = JCSFBSAV
      ENDIF

      LOSCAL = .FALSE.
      LRKCO = .FALSE.

      IF (JB /= JA) THEN
         LINCR = .TRUE.
      ELSE
         NPOS = NPOS + 1
         LINCR = .FALSE.
      ENDIF
      LDIAG = .FALSE.
      LOFFD = .FALSE.
! CSFs generated by the JB-th CSFG, Row demension
      M = NTYPE(2,JB)
! CSFs generated by the JA-th CSFG, Column demension
      N = NTYPE(2,JA)
! Number of symmetry-ordered orbs of JA-CSFG
      NORBCOL = NTYPE(4,JA) - NTYPE(3,JA) + 1
! Number of symmetry-ordered orbs of JB-CSFG 
      NORBROW = NTYPE(4,JB) - NTYPE(3,JB) + 1
! Check
      IF (NORBCOL /= N .OR. NORBROW /= M) THEN
        WRITE(*,*)"NORBCOL /= N .OR. NORBROW /= M, ", &
                  NORBCOL, N, NORBROW, M
        STOP "NORBCOL /= N .OR. NORBROW /= M IN SpinAngular22.f90 ..."
      ENDIF
!
! There is possible onebody contributions between TYPE 2 - 2 CSFGs.
!
      IF (NTYPE(3,JA) .EQ. NTYPE(3,JB)) THEN
        FLAGU = .TRUE.
        IRW = IRFICT
        ICW = ICFICT
      ELSE
        FLAGU = .FALSE.
        IRW = JB
        ICW = JA
      ENDIF
      IF (.NOT. FLAGU) GOTO 101

! IF (.NOT. LCHM), Calculate Potentials, T-contributions calculated
! in SETCOF_T_CSFG
!      IF (.NOT. LCHM) GOTO 101

      IA = 0
      TEGRAL = 0.D0
      CALL read_TV(1, JASAV, JBSAV, IUNITF, 1, IA, IB, TCOEFF)
      IF (IA /= 0 .AND. IA /= IB .AND.       &
            ABS (TCOEFF) > CUTOFF) THEN
        IF (IA > IB) THEN
          ISWAP = IB
          IB = IA
          IA = ISWAP
        ENDIF
        IF (LABTVFRST) NTPT = NTPT + 1
!
! TSHELL(1) should be used for the min(M,N) diagonal matrix elements
! without any changes, as IA and IB are both labelling orbitals. 
!        
        ! Diagonal matrixelement needed
        IF (LCHM) THEN
          CALL IABINT(IA, IB, TEGRAL)
          ENONSYM = TCOEFF * TEGRAL
        ENDIF

        DO I = 1, min(NORBCOL, NORBROW)
          IF (LABTVFRST.OR.LCPOT) THEN
            CALL SETCOF_NDA_CSFG(I, I, IA, IB, TCOEFF)
          ELSE
            EMTBLOCK(I,I) = ENONSYM
          ENDIF
        ENDDO
      ENDIF

101   CONTINUE
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs. NVCOEF is initialized here but changed
!   in /rkco/cor[d] via COMMON.
!
      ENONSYM = 0.D0
      NVCOEF = 0
      CALL read_TV(2, JASAV, JBSAV, IUNITF, 1, IA, IB, TSHELL1)
      !IF (NVCOEF.GT.0) &
      !  CALL PRINTLABELV(2, JA, JB, 22, 1, 0, 0, 0.0d0)
      DO 2 IV = 1, NVCOEF
        VCOEFF = COEFF(IV)

! Determine the position of the symmetry-ordered orbitals
        LABV(1:5) = LABEL(1:5, IV)
        IF (LABTVFRST) NTPV = NTPV + 1
        CALL ANALABV(JASAV, JBSAV, IV)

        IF (NSYMCR .EQ. 0) THEN
          IF (LABTVFRST.OR.LCPOT) THEN
            DO I = 1, MIN(NTYPE(2,JB), NTYPE(2,JA))
              CALL SETCOF_NXY_CSFG(I, I, VCOEFF) 
            ENDDO
          ELSE
            CALL twobody_DC(TEGRAL)
            ENONSYM = ENONSYM + VCOEFF * TEGRAL
          ENDIF 

! Matrixelement between Type 1 - 2, 5, Loop for
! symmetry-ordered-orbitals
        ELSEIF (NSYMCR .EQ. 1) THEN
          WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
          STOP 'Warning!!! NSYMCR ERROR ***1*** 22'

        ELSEIF (NSYMCR .EQ. 2) THEN
! Replace the symmetry-ordered-orbital within LABV, the position is
! determied by IPSym(1), IPSym(2)
          DO N = 1, NORBCOL
            LABV(IPSym(1)) = NTYPE(3, JA) + N - 1
            DO M = 1, NORBROW
              LABV(IPSym(2)) = NTYPE(3, JB) + M - 1
              IF (LABTVFRST.OR.LCPOT) THEN
                CALL SETCOF_NXY_CSFG(M, N, VCOEFF)
              ELSE
                CALL twobody_DC(TEGRAL)
                EMTBLOCK(M,N) = EMTBLOCK(M,N) + VCOEFF * TEGRAL
              ENDIF
            ENDDO
          ENDDO

        ELSEIF (NSYMCR .EQ. 3) THEN
          WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
               'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
          STOP 'Warning!!! NSYMCR ERROR ***3*** 22'

        ELSEIF (NSYMCR .EQ. 4) THEN
          WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
               'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
          STOP 'Warning!!! NSYMCR ERROR ***4*** 22'
        ENDIF
    2 CONTINUE
!
! Add the common parts for the semi-diagonal matrixelement
!
      IF (LCHM.AND.DABS(ENONSYM) .GT. CUTOFF) THEN
        DO I = 1, MIN(NORBCOL, NORBROW)
           EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM
        ENDDO
      ENDIF

      IF (LCHM.AND.LTRANSFER) THEN
        IF (LTRANSPOSE) EMTBLOCK=TRANSPOSE(EMTBLOCK)
        CALL TRANSFER_CSFG(JASAV+NCFGPAST, JBSAV+NCFGPAST)
      ENDIF

      RETURN
      END SUBROUTINE SPINANGULAR22