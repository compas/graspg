!***********************************************************************
!                                                                      *
      SUBROUTINE SPINANGULAR44 (JASAV, JBSAV, IUNITF, NPOS, JCSFASAV, JCSFBSAV)
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
!   Modify for Spin-Angular integration: (TYPE 4 and 4 CSFGS)          *
!              <[core A] (n-1)k nk | H_DC | [core B] (n'-1)k' n'k'>    *
!   [core] is built by labelling symmetries,                           *
!   nks are kappa-ordered orbitals                                     *
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

      LOGICAL :: F0INT, LRKCO, FLAGL, FLAGU, LSYM2, LSYM4, LSYM0
      INTEGER :: NDIFF1,NDIFF2,NORBCOL,NORBROW,LABTMP
      INTEGER :: IV, I, J, K, L, M, N, IRW, ICW, IORBMIN, MTYPE
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

! CSFs generated by the JB-th CSFG, Row demension
      M = NTYPE(2,JB)
! CSFs generated by the JA-th CSFG, Column demension
      N = NTYPE(2,JA)
! Number of symmetry-ordered orbs of JA-CSFG
      NORBCOL = NTYPE(6,JA) - NTYPE(3,JA) + 1
! Number of symmetry-ordered orbs of JB-CSFG 
      NORBROW = NTYPE(6,JB) - NTYPE(3,JB) + 1
! The lowest orbital of the involving symmetry.
      IORBMIN = NTYPE(3,JA) 
!
! To use the lowest CSFs pair by default, then if the JB CSFG has one of
! the TWO symmetry-ordered orbitals of JA CSFG, the IRW and ICW meet the
! requirement that they have the same orbital by default. 
!
      IRW = IRFICT
      ICW = ICFICT
      ! <IR [Core A] 6s 7s|  ==> <IR [Core A] 4s 5s|
      !CALL FICTIOUS_CSF(4, IRFICT, JB, NTYPE(3,JB), NTYPE(4,JB), &
      !                               NTYPE(3,JB)+1, NTYPE(6,JB))
      ! |IC [Core B] 6s 7s>| ==> |IC [Core B] 4s 5s>
      !CALL FICTIOUS_CSF(4, ICFICT, JA, NTYPE(3,JA), NTYPE(4,JA), &
      !                               NTYPE(3,JA)+1, NTYPE(6,JA))

      FLAGL = .FALSE.
      FLAGU = .FALSE.
      IF (NTYPE(3,JB) /= NTYPE(3,JA)) GOTO 101

! Both JA and JB have the same Kappa-ordered smmetry orbitals.
! There might be onebody contributions. 

      FLAGL = .TRUE.
      IF (NORBCOL >= 3 .OR. NORBROW >= 3) THEN
        FLAGU = .TRUE.
        ! |IC [Core B] 6s 7s>| ==> |IC [Core B] 5s 6s>
        !CALL FICTIOUS_CSF(4, ICFICT+1, JA, NTYPE(3,JA)+1, NTYPE(4,JA), &
        !                                   NTYPE(3,JA)+2, NTYPE(6,JA))
      ! <IR [Core A] 6s 7s|  ==> <IR [Core A] 5s 6s|
        !CALL FICTIOUS_CSF(4, IRFICT+1, JB, NTYPE(3,JB)+1, NTYPE(4,JB), &
        !                                   NTYPE(3,JB)+2, NTYPE(6,JB))
      ENDIF 
!
! There are onebody contributions for semi-diagonal and some off-diagonal
! elements of the block generated by the JA - JB CSFGs.
!
! Onebody contributions.
!
! Semi diagonal elements
! MTYPE = 1: < [Core] K L | h1 | [Core] I J>
! < (Core A) K L (IRFICT) | H1 | (Core B) I J (ICFICT) >
! K = I, L  = J         : 4s 5s (IRFICT) -- 4s 5s (ICFICT)
!
! Off-diagonal elements
! < (Core A) K L (IRFICT) | H1 | (Core B) I J (ICFICT+1) >
! L > K = J > I,  5s 6s -- 4s 5s
! K < L = I < J,  4s 5s -- 5s 6s
!
      ENONSYM = 0.D0
      DO MTYPE = 1, 2
!
! MTYPE = 2: < [Core] K L | h1 | [Core] I J>, K < L = I < J
! 
        IF (MTYPE == 1) THEN
          IRW = IRFICT
          ICW = ICFICT
        ELSE
          IF (.NOT. FLAGU) CYCLE
          IRW = IRFICT
          ICW = ICFICT + 1
        ENDIF     
        ! Read T coefficients
        IA = 0
        IB = 0
        CALL read_TV(1, JASAV, JBSAV, IUNITF, MTYPE, IA, IB, TCOEFF)
        IF (IA /= 0 .AND. IA /= IB .AND.       &
               ABS (TCOEFF) > CUTOFF) THEN

          IF (LABTVFRST) NTPT = NTPT + 1

          IF (MTYPE == 1) THEN
          ! Semi-diagonal elements, the one-body contribution arise from
          ! the labelling orbitals.
            IF (IA > NORBGEN .OR. IB > NORBGEN) THEN
              WRITE(*,*)'JASAV, JBSAV, IA, IB =', JASAV, JBSAV, IA, IB
              STOP 'IA > NORBGEN .OR. IB > NORBGEN in 44.f90 ...'
            ENDIF
            IF (IA > IB) THEN
              ISWAP = IA
              IA = IB
              IB = ISWAP 
            ENDIF

            IF (LCHM) THEN
              CALL IABINT(IA, IB, TEGRAL)
              ENONSYM = TCOEFF * TEGRAL
            ENDIF

            ! ENONSYM is common to all Semi-diagonal matrixelements
            IF (NTYPE(2,JB).EQ.NTYPE(2,JA)) THEN
              DO I = 1, NTYPE(2,JA)
                IF (LABTVFRST.OR.LCPOT) THEN
                  CALL SETCOF_NDA_CSFG(I, I, IA, IB, TCOEFF)
                ELSE
                  EMTBLOCK(I,I) = ENONSYM 
                ENDIF
              ENDDO
            ELSE
              ! NTYPE(2,IC).NE.NTYPE(2,IR). There might be more effective
              ! method to find the position of the semi-diagonal
              ! matrixelements. 
              N = 0 ! Column index
              DO I = 1, NORBCOL - 1
               DO J = I + 1, NORBCOL
                N = N + 1
                M = 0 ! Row index
                DO K = 1, NORBROW - 1
                 DO L = K + 1, NORBROW
                  M = M + 1
                  IF (K.EQ.I .AND. L.EQ.J) THEN
                    IF (LABTVFRST.OR.LCPOT) THEN
                      CALL SETCOF_NDA_CSFG(M, N, IA, IB, TCOEFF)
                    ELSE
                      EMTBLOCK(M,N) = ENONSYM 
                    ENDIF
                  ENDIF
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
            ENDIF
          
          ELSEIF(MTYPE == 2) THEN
            IF (IA.LE.NORBGEN .OR. IB.LE.NORBGEN) THEN 
              WRITE(*,*)'JASAV, JBSAV, IA, IB =', JASAV, JBSAV, IA, IB
              STOP 'IA.LE.NORBGEN .OR. IB.LE.NORBGEN in 44.f90 ...'
            ENDIF
            M = 0
            DO K = 1, NORBROW - 1
              DO L = K + 1, NORBROW
                M = M + 1 
                N = 0
                DO I = 1, NORBCOL - 1
                  DO J = I + 1, NORBCOL
                    N = N + 1
                    ! Full matrix needed.
                    IF (K == J) THEN
                    ! L > K = J > I
                      IA = I + NTYPE(3,JA) - 1
                      IB = L + NTYPE(3,JB) - 1
                      IF (LABTVFRST.OR.LCPOT) THEN
                        CALL SETCOF_NDA_CSFG(M, N, IA, IB, TCOEFF)
                      ELSE
                        CALL IABINT(IA, IB, TEGRAL)
                        EMTBLOCK(M,N) = TCOEFF * TEGRAL 
                      ENDIF

                    ELSEIF (L == I) THEN
                    ! K < L = I < J
                      IA = K + NTYPE(3,JB) - 1
                      IB = J + NTYPE(3,JA) - 1
                      IF (LABTVFRST.OR.LCPOT) THEN
                        CALL SETCOF_NDA_CSFG(M, N, IA, IB, TCOEFF)
                      ELSE
                        CALL IABINT(IA, IB, TEGRAL)
                        EMTBLOCK(M,N) = TCOEFF * TEGRAL 
                      ENDIF
                    ENDIF 
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF 
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
      ENONSYM = 0.D0
      DO MTYPE = 1, 3
        IF (MTYPE == 1) THEN
        ! MTYPE = 1 :
        ! <IR K L | h1 | IC  I J> for K = I, L = J
        ! <IR K L | h1 | IC  I J> for K = I, L /= J
        ! < 4s 5s | 4s 5s> 
          IRW = IRFICT
          ICW = ICFICT
        ELSEIF (.NOT. FLAGL) THEN
        ! Different Kappa-ordered orbitals, only MTYPE = 1 is needed.
          EXIT
        ELSEIF (MTYPE == 2) THEN
        ! <IR K L | h1 | IC  I J> for K < L = I < J
        ! |IRFICT    , [Core A] 4s 5s>
        ! |ICFICT + 1, [Core B] 5s 6s>
        ! < 4s 5s | 5s 6s>
          IF (NORBCOL .LT. 3) CYCLE
          IRW = IRFICT
          ICW = ICFICT + 1
        ELSEIF(MTYPE == 3) THEN
        ! <IR K L | h1 | IC  I J> for L > K = J > I
        ! |IRFICT + 1, [Core A] 5s 6s>
        ! |ICFICT    , [Core B] 4s 5s>
        ! < 5s 6s | 4s 5s>
        ! Could this kind of match be removed? 
          IF (NORBROW .LT.3) CYCLE
          IRW = IRFICT + 1
          ICW = ICFICT
        ENDIF

        NVCOEF = 0
        CALL read_TV(2, JASAV, JBSAV, IUNITF, MTYPE, IA, IB, TSHELL1)
        DO 2 IV = 1, NVCOEF
           VCOEFF = COEFF(IV)

! Determine the position of the symmetry-ordered orbitals
          LABV(1:5) = LABEL(1:5, IV)
          IF (LABTVFRST) NTPV = NTPV + 1
          CALL ANALABV(JASAV, JBSAV, IV)

          IF (NSYMCR .EQ. 0) THEN
            IF (LABTVFRST.OR.LCPOT) THEN
              N = 0
              DO I = 1, NORBCOL - 1
                DO J = I + 1, NORBCOL
                  N = N + 1
                  M = 0
                  DO K = 1, NORBROW - 1
                    DO L = K + 1, NORBROW
                      M = M + 1
                      IF (K == I .AND. L == J) &
                        CALL SETCOF_NXY_CSFG(M, N, VCOEFF)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              CALL twobody_DC(TEGRAL)
              ENONSYM = ENONSYM + VCOEFF * TEGRAL
            ENDIF

          ELSEIF (NSYMCR .EQ. 1) THEN
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                 'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***1*** 44'

          ELSEIF (NSYMCR .EQ. 2) THEN
! Replace the symmetry-ordered-orbital within LABV, the position is
! determied by IPSym(1), IPSym(2)
!    IF (MTYPE == 1) THEN
!      IF (IORBCOL == NTYPE(3,JA) .AND. L == J) &
!                                 LNonZero(M,N) = .TRUE.
!      IF (IORBCOL == NTYPE(5,JA) .AND. K == I) &
!                                 LNonZero(M,N) = .TRUE.
!    ENDIF
!    IF (MTYPE == 2 .AND. L == I) LNonZero(M,N) = .TRUE.
!    IF (MTYPE == 3 .AND. K == J) LNonZero(M,N) = .TRUE.

            N = 0
            DO I = 1, NORBCOL - 1
              DO J = I + 1, NORBCOL
                N = N + 1
                M = 0
                DO K = 1, NORBROW - 1
                  DO L = K + 1, NORBROW 
                    M = M + 1
                    ! full matrix needed
                    IF (MTYPE == 1) THEN
        ! <IR K L | h12 | IC  I J> for K = I, L = J
        ! <IR K L | h12 | IC  I J> for K = I, L /= J
                      IF (L == I) CYCLE ! MTYPE 2
                      IF (K == J) CYCLE ! MTYPE 3

                ! K = I, L = J (Semi-diagonal) / L!= J (Off-diagonal)
                IF (K.EQ.I.AND.LABEL(IPSYM(1),IV).EQ.IORBMIN+1) THEN
                ! < (Core A) K L | 1/r12 | (Core B) I J >: K and I
                ! orbitals are orthogonal normalized, leaving TWO
                ! labelling orbitals and TWO symmetry-ordered orbitals
                ! involving L and J. 
                  LABV(IPSYM(1)) = NTYPE(3, JA) + J - 1
                  LABV(IPSYM(2)) = NTYPE(3, JB) + L - 1
                  IF (LABTVFRST.OR.LCPOT) THEN
                    CALL SETCOF_NXY_CSFG(M, N, VCOEFF)
                  ELSE
                    CALL twobody_DC(TEGRAL)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + VCOEFF * TEGRAL
                  ENDIF 
                ENDIF
                ! L = J, K = I (Semi-diagonal) / K!= I (Off-diagonal)
                IF (L.EQ.J.AND.LABEL(IPSYM(1),IV).EQ.IORBMIN) THEN
                ! < (Core A) K L | 1/r12 | (Core B) I J >: L and J
                ! orbitals are orthogonal normalized, leaving TWO
                ! labelling orbitals and TWO symmetry-ordered orbitals
                ! involving K and I. 
                  LABV(IPSYM(1)) = NTYPE(3, JA) + I - 1
                  LABV(IPSYM(2)) = NTYPE(3, JB) + K - 1
                  IF (LABTVFRST.OR.LCPOT) THEN
                    CALL SETCOF_NXY_CSFG(M, N, VCOEFF)
                  ELSE
                    CALL twobody_DC(TEGRAL)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + VCOEFF * TEGRAL
                  ENDIF
                ENDIF

                    ELSEIF (MTYPE == 2 .AND. L == I) THEN
                ! <IR K L | h12 | IC  I J> for K < L = I < J
                ! L = I, K!= J;  4s 5s -- 5s 6s
                ! < (Core A) K L | 1/r12 | (Core B) I J >: L and I
                ! orbitals are orthogonal normalized, leaving TWO
                ! labelling orbitals and TWO symmetry-ordered orbitals
                ! involving K and J. 

                      !IF (LABEL(IPSYM(1),IV).EQ.NTYPE(3,JA)+2) THEN
                        LABV(IPSYM(1)) = NTYPE(3, JA) + J - 1
                        LABV(IPSYM(2)) = NTYPE(3, JB) + K - 1
                        IF (LABTVFRST.OR.LCPOT) THEN
                          CALL SETCOF_NXY_CSFG(M, N, VCOEFF)
                        ELSE
                          CALL twobody_DC(TEGRAL)
                          EMTBLOCK(M,N) = EMTBLOCK(M,N) + VCOEFF * TEGRAL
                        ENDIF
                      !ENDIF
                    ELSEIF (MTYPE == 3 .AND. K == J) THEN
                ! <IR K L | h12 | IC  I J> for L > K = J > I and K > I
                ! K = J, L!= I;  5s 6s -- 4s 5s
                ! < (Core A) K L | 1/r12 | (Core B) I J >: K and J
                ! orbitals are orthogonal normalized, leaving TWO
                ! labelling orbitals and TWO symmetry-ordered orbitals
                ! involving L and I. 
                      !IF (LABEL(IPSYM(1),IV).EQ.NTYPE(3,JA)) THEN
                        LABV(IPSYM(1)) = NTYPE(3, JA) + I - 1
                        LABV(IPSYM(2)) = NTYPE(3, JB) + L - 1
                        IF (LABTVFRST.OR.LCPOT) THEN
                          CALL SETCOF_NXY_CSFG(M, N, VCOEFF)
                        ELSE
                          CALL twobody_DC(TEGRAL)
                          EMTBLOCK(M,N) = EMTBLOCK(M,N) + VCOEFF * TEGRAL
                        ENDIF
                      !ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO  
            ENDDO

          ELSEIF (NSYMCR .EQ. 3) THEN
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                 'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***3*** 44'

          ELSEIF (NSYMCR .EQ. 4) THEN
            IF (MTYPE /= 1) THEN
            ! In rangular_mpi_csfg: savetv.f90, V-Coefficients with 
            ! NSYMCR = 4 and MTYPE = 2, 3 are discarded.
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                   'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
              STOP 'Warning!!! NSYMCR ERROR ***4*** 44'
            ENDIF
            N = 0
            DO I = 1, NORBCOL - 1
              LABV(1) = NTYPE(3, JA) + I - 1
              DO J = I + 1, NORBCOL
                N = N + 1
                LABV(2) = NTYPE(3, JA) + J - 1
                M = 0
                DO K = 1, NORBROW - 1
                  DO L = K + 1, NORBROW
                    M = M + 1
                    ! Full matrix needed
                    LABV(3) = NTYPE(3, JB) + K - 1
                    LABV(4) = NTYPE(3, JB) + L - 1
                    IF (LABEL(3,IV).GT.LABEL(4,IV)) THEN
                      LABTMP = LABV(3)
                      LABV(3) = LABV(4)
                      LABV(4) = LABTMP
                    ENDIF
                    IF (LABTVFRST.OR.LCPOT) THEN
                      CALL SETCOF_NXY_CSFG(M, N, VCOEFF)
                    ELSE
                      CALL twobody_DC(TEGRAL)
                      EMTBLOCK(M,N) = EMTBLOCK(M,N) + VCOEFF * TEGRAL
                    ENDIF 
                  ENDDO 
                ENDDO
              ENDDO
            ENDDO
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
    2   CONTINUE
      ENDDO
!
! Add the common parts for the diagonal matrixelement
!
      IF (LCHM.AND.DABS(ENONSYM) .GT. CUTOFF) THEN
        N = 0 
        DO I = 1, NORBCOL - 1
          DO J = I + 1, NORBCOL
            N = N + 1
            M = 0
            DO K = 1, NORBROW - 1
              DO L = K + 1, NORBROW
                M = M + 1
                IF (K == I .AND. L == J) &
                EMTBLOCK(M,N) = EMTBLOCK(M,N) + ENONSYM
              ENDDO
            ENDDO
          ENDDO  
        ENDDO
      ENDIF

      IF (LCHM.AND.LTRANSFER) THEN
        IF (LTRANSPOSE) EMTBLOCK=TRANSPOSE(EMTBLOCK)
        CALL TRANSFER_CSFG(JASAV+NCFGPAST, JBSAV+NCFGPAST)
      ENDIF

      RETURN
      END SUBROUTINE