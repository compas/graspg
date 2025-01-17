!***********************************************************************
!                                                                      *
      SUBROUTINE SPINANGULAR13 (JASAV, JBSAV, IUNITF, NPOS, JCSFASAV, JCSFBSAV)
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
!   Modify for Spin-Angular integration: (Normal CSF and TYPE 3 CSFG)  *
!              <CSF         | H_DC | [core] nk n'k'>                   *
!   nk, n'k' are symmetry-ordered (by kappa).                          *
!   [core] is built by labelling symmetries.                           *
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

      LOGICAL :: F0INT, LOSCAL, LRKCO
      INTEGER :: NDIFF1,NDIFF2,NORBCOLL,NORBCOLU,NORBROWL,NORBROWU
      INTEGER :: IV, I, J, K, L, M, N
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

!   CSFs generated by the JA-th CSFG.        
      N = NTYPE(2,JA)

! Number of symmetry-ordered orbs of JA-CSFG
      NORBCOLL = NTYPE(4,JA) - NTYPE(3,JA) + 1
      NORBCOLU = NTYPE(6,JA) - NTYPE(5,JA) + 1
     
!
! There is no one-body contributions between TYPE 1 and 3 CSFGs
! 
      IA = 0
      IB = 0
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs. NVCOEF is initialized here but changed
!   in /rkco/cor[d] via COMMON.
!
      NVCOEF = 0
      CALL read_TV(2, JASAV, JBSAV, IUNITF, 1, IA, IB, TSHELL1)
      DO 2 IV = 1, NVCOEF
        VCOEFF = COEFF(IV)

! Determine the position of the symmetry-ordered orbitals
        LABV(1:5) = LABEL(1:5, IV)
        IF (LABTVFRST) NTPV = NTPV + 1
        !IF (LABTVFRST) CALL SAVEP_VK(VCOEFF)
        CALL ANALABV(JASAV, JBSAV, IV)

        IF (NSYMCR .EQ. 0) THEN
          WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
          STOP 'Warning!!! NSYMCR ERROR ***1*** 13'

! Matrixelement between Type 1 - 2, 5, Loop for
! symmetry-ordered-orbitals
        ELSEIF (NSYMCR .EQ. 1) THEN
          WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
          STOP 'Warning!!! NSYMCR ERROR ***1*** 13'

        ELSEIF (NSYMCR .EQ. 2) THEN
! Replace the symmetry-ordered-orbital within LABV, the position is
! determied by IPSym(1), IPSym(2)
            N = 0  ! Column index
            DO I = 1, NORBCOLL
              LABV(IPSym(1)) = NTYPE(3, JA) + I - 1
              DO J = 1, NORBCOLU
                N = N + 1
                LABV(IPSym(2)) = NTYPE(5, JA) + J - 1
                IF (LABTVFRST.OR.LCPOT) THEN
                  CALL SETCOF_NXY_CSFG(1, N, VCOEFF)
                ELSE
                  CALL twobody_DC(TEGRAL)
                  EMTBLOCK(1,N) = EMTBLOCK(1,N) + VCOEFF * TEGRAL
                ENDIF
              ENDDO
            ENDDO

        ELSEIF (NSYMCR .EQ. 3) THEN
          WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
               'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
          STOP 'Warning!!! NSYMCR ERROR ***3*** 13'

        ELSEIF (NSYMCR .EQ. 4) THEN
          WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
               'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
          STOP 'Warning!!! NSYMCR ERROR ***4*** 13'
        ENDIF
    2 CONTINUE

      IF (LCHM.AND.LTRANSFER) THEN
        IF (LTRANSPOSE) EMTBLOCK=TRANSPOSE(EMTBLOCK)
        CALL TRANSFER_CSFG(JASAV+NCFGPAST, JBSAV+NCFGPAST)
      ENDIF

      RETURN
      END SUBROUTINE SPINANGULAR13
