!***********************************************************************
!                                                                      *
      SUBROUTINE SPINANGULAR4  (JASAV, JBSAV, IUNITF, NPOS, JCSFASAV, JCSFBSAV)
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
!   Modify for Spin-Angular integration: (Same TYPE 4 CSFGS)           *
!              <[core A] (n-1)k nk | H_DC | [core A] (n-1)k nk>        *
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
      REAL(DOUBLE) :: ENONSYMA, ENONSYMB
      
      LOGICAL :: F0INT, LOSCAL, LRKCO, FLAGLI
      INTEGER :: NDIFF1,NDIFF2,NORBCOL,NORBROW
      INTEGER :: IV, I, J, K, L, M, N, IRW, ICW, IORBMIN, MTYPE, LABTMP
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
! Check, it should be diagnoal here, i.e., IRTOT.eq.ICTOT
      IF (IRTOT /= ICTOT) THEN
        WRITE(*,*)"IRTOT /= ICTOT, ",IRTOT, ICTOT
        STOP "IRTOT /= ICTOT in spinangular3.f90 ..."
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
! V-coefficients involving only labelling orbitals are not output, if
! fictious CSFs are used for diagonal elements, same as spinangular5.f90
      ! 6s 7s
      IRW = JB
      ! 6s 7s
      ICW = JA

      FLAGLI = .FALSE.
      IF (NORBCOL >= 3) THEN
        FLAGLI = .TRUE.
        !! 6s 7s (JB, Here JB =  JA) ==> 5s 7s 
        !CALL FICTIOUS_CSF(4, IRFICT+1, JB, NTYPE(4,JB)-1, NTYPE(4,JB), &
        !                                             0,           0)
        !! 6s 7s (JB, Here JB =  JA) ==> 5s 6s 
        !CALL FICTIOUS_CSF(4, IRFICT+2, JB, NTYPE(4,JB)-1, NTYPE(4,JB), &
        !                                   NTYPE(6,JB)-1, NTYPE(6,JB))
      ENDIF 
!
! There are onebody contributions for diagonal and some off-diagonal
! elements between the diagonal block generated by the JA - JB ( = JA) CSFGs.
!
! Onebody contributions.
!
! Diagonal elements  
! MTYPE = 1: < [Core] K L | h1 | [Core] I J>
! K = I, L = J, TSHELL = nqa, occupation number
!
      IF (LCHM) THEN
        ENONSYM = 0.D0
        DO IA = 1, NORBGEN
          TCOEFF = IQA(IA, JCSFA)
          IF (ABS(TCOEFF) .GT. CUTOFF) THEN
            CALL IABINT(IA, IA, TEGRAL)
            ENONSYM = ENONSYM + TCOEFF * TEGRAL
          ENDIF
        ENDDO
        TCOEFF = 1.D0
        N = 0
        DO I = 1, NORBCOL - 1
          IA = I + NTYPE(3,JA) - 1
          CALL IABINT(IA, IA, ENONSYMA)
          DO J = I + 1, NORBCOL
            N = N + 1  
            IB = J + NTYPE(3,JA) - 1
            CALL IABINT(IB, IB, ENONSYMB)
            EMTBLOCK(N, N) = ENONSYM + ENONSYMA + ENONSYMB
          ENDDO
        ENDDO
      ENDIF
!
! Off-diagonal elements, Part 1:
! TSHELL(1) = 1.0d0 for K = I and L < J, or K < I and L = J
! Off-diagonal, < K L | h1 | I J> : A) K = I, L < J; B) K < I,  L = J
! For above two cases, TCOEFF = 1.0d0

! In case of NTYPE(2,JA) == NTYPE(2,JB) == 1, there is no off-diagonal
! elements, hence there is no T-Coef.
      IF (NTYPE(2,JA) > 1) THEN
      ! i.e., NORBCOL >= 3  
        TCOEFF = 1.D0 
        IA = NTYPE(3,JA)
        IB = NTYPE(4,JA)
        IF (LABTVFRST) NTPT = NTPT + 1
      ENDIF

      M = 0
      DO K = 1, NORBROW - 1
        DO L = K + 1, NORBROW
          M = M + 1
          N = 0
          DO I = 1, NORBCOL - 1
            DO J = I + 1, NORBCOL
              N = N + 1
              ! Upper-triangle matrix needed
              IF (K > I) CYCLE
              IF (K == I .AND. L > J) CYCLE

              ! i.e., NORBCOL >= 3  
              IF (K == I .AND. L < J) THEN
                IA = L + NTYPE(3, JB) - 1
                IB = J + NTYPE(3, JA) - 1
                IF (LABTVFRST.OR.LCPOT) THEN
                  CALL SETCOF_NDA_CSFG(M, N, IA, IB, TCOEFF)
                ELSE
                  CALL IABINT(IA, IB, TEGRAL)
                  EMTBLOCK(M,N) = TEGRAL 
                ENDIF
              ENDIF

              ! i.e., NORBCOL >= 3  
              IF (K < I .AND. L == J) THEN
                IA = K + NTYPE(3, JB) - 1
                IB = I + NTYPE(3, JA) - 1
                IF (LABTVFRST.OR.LCPOT) THEN
                  CALL SETCOF_NDA_CSFG(M, N, IA, IB, TCOEFF)
                ELSE
                  CALL IABINT(IA, IB, TEGRAL)
                  EMTBLOCK(M,N) = TEGRAL 
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

! Off-diagonal elements, Part 2:
!     K < L = I <  J  , calculated explicitly.
      DO MTYPE = 2, 2
!
! MTYPE = 2: < [Core] K L | h1 | [Core] I J>, K < L = I < J
! 
        IF (MTYPE == 1) THEN
          IRW = JB
          ICW = JA
        ELSE
          !IF (.NOT. FLAGLI) CYCLE
          IF (NORBCOL < 3) CYCLE
          !IRW = IRFICT + 2
          !ICW = JA 
        ENDIF     
        ! Read T coefficients
        IA = 0
        IB = 0
        CALL read_TV(1, JASAV, JBSAV, IUNITF, MTYPE, IA, IB, TCOEFF)
        IF (IA /= 0 .AND. IA /= IB .AND.       &
               ABS (TCOEFF) > CUTOFF) THEN
          IF (LABTVFRST) NTPT = NTPT + 1
          M = 0
          DO K = 1, NORBROW - 1
            DO L = K + 1, NORBROW
              M = M + 1 
              N = 0
              DO I = 1, NORBCOL-1
                DO J = I + 1, NORBCOL
                  N = N + 1
                  ! Upper-triangle matrix needed.
                  IF (K > I) CYCLE
                  IF (K == I .AND. L > J) CYCLE 
                  IF (MTYPE == 2 .AND. L == I) THEN
                  ! K < L = I < J
                    IA = K + NTYPE(3, JB) - 1
                    IB = J + NTYPE(3, JA) - 1
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
      ENDDO
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs. NVCOEF is initialized here but changed
!   in /rkco/cor[d] via COMMON.
!
200   CONTINUE
      ENONSYM = 0.D0
      DO MTYPE = 1, 3
        IF (MTYPE == 1) THEN
        ! MTYPE = 1 :
        ! <IR K L | h12 | IC  I J> for K = I, L = J
        ! <IR K L | h12 | IC  I J> for K = I, L /= J
          IRW = JB
          ICW = JA
        ELSEIF(MTYPE == 2) THEN
        ! <IR K L | h12 | IC  I J> for K < L = I < J
        !  K < L = I < J: <5s 6s | h12 | 6s 7s>
          !IF (.NOT. FLAGLI) CYCLE
          IF (NORBCOL < 3) CYCLE
          IRW = IRFICT + 2
          ICW = JA
        ELSEIF(MTYPE == 3) THEN
        ! <IR K L | h12 | IC  I J> for K < L = J > I and K < I
        !  K < L = I < J: <5s 7s | h12 | 6s 7s>
          !IF (.NOT. FLAGLI) CYCLE
          IF (NORBCOL < 3) CYCLE
          IRW = IRFICT + 1
          ICW = JA
        ENDIF

        NVCOEF = 0
        CALL read_TV(2, JASAV, JBSAV, IUNITF, MTYPE, IA, IB, TSHELL1)
        !IF (NVCOEF > 0) &
        !  CALL PRINTLABELV(2, JA, JB, 4, MTYPE, 0, 0, 0.0d0)
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

          ELSEIF (NSYMCR .EQ. 1) THEN
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                 'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***1*** 4'

          ELSEIF (NSYMCR .EQ. 2) THEN
! Replace the symmetry-ordered-orbital within LABV, the position is
! determied by IPSym(1), IPSym(2)
            N = 0
            DO I = 1, NORBCOL - 1
              DO J = I + 1, NORBCOL
                N = N + 1
                M = 0
                DO K = 1, NORBROW - 1
                  DO L = K + 1, NORBROW 
                    M = M + 1
                    ! Upper-triangle matrix needed
                    IF (K.GT.I) CYCLE
                    IF (K == I .AND. L .GT. J) CYCLE
                    ! Note, here K is always less than J, because J > I >= K
                    IF (MTYPE == 1) THEN
        ! <IR K L | h12 | IC  I J> for K = I, L = J
        ! <IR K L | h12 | IC  I J> for K = I, L /= J
                      IF (K == I) THEN
                        IF ( L == J) THEN
                        ! Diagonal elements 
                         IF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,JA)-1) THEN
                           LABV(IPSYM(1)) = NTYPE(3, JA) + I - 1
                           LABV(IPSYM(2)) = NTYPE(3, JB) + K - 1
                         ELSEIF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,JA)) THEN
                           LABV(IPSYM(1)) = NTYPE(3, JA) + J - 1
                           LABV(IPSYM(2)) = NTYPE(3, JB) + L - 1
                         ENDIF
                         IF (LABTVFRST.OR.LCPOT) THEN
                           CALL SETCOF_NXY_CSFG(M, N, VCOEFF)
                         ELSE
                           CALL twobody_DC(TEGRAL)
                           EMTBLOCK(M,N) = EMTBLOCK(M,N) + VCOEFF * TEGRAL
                         ENDIF
                        ELSEIF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,JA)) THEN
                        ! K = I, L < J 
                         LABV(IPSYM(1)) = NTYPE(3, JA) + J - 1
                         LABV(IPSYM(2)) = NTYPE(3, JB) + L - 1
                         IF (LABTVFRST.OR.LCPOT) THEN
                           CALL SETCOF_NXY_CSFG(M, N, VCOEFF)
                         ELSE
                           CALL twobody_DC(TEGRAL)
                           EMTBLOCK(M,N) = EMTBLOCK(M,N) + VCOEFF * TEGRAL
                         ENDIF
                        ENDIF
                      ENDIF
                    ELSEIF (MTYPE == 2 .AND. L == I) THEN
        ! <IR K L | h12 | IC  I J> for K < L = I < J
                      IF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,JA)) THEN
                        LABV(IPSYM(1)) = NTYPE(3, JA) + J - 1
                        LABV(IPSYM(2)) = NTYPE(3, JB) + K - 1
                        IF (LABTVFRST.OR.LCPOT) THEN
                          CALL SETCOF_NXY_CSFG(M, N, VCOEFF)
                        ELSE
                          CALL twobody_DC(TEGRAL)
                          EMTBLOCK(M,N) = EMTBLOCK(M,N) + VCOEFF * TEGRAL
                        ENDIF
                      ENDIF
                    ELSEIF (MTYPE == 3 .AND. L == J) THEN
        ! <IR K L | h12 | IC  I J> for K < L = J > I and K < I
                      IF (K >= I) CYCLE 
                      IF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,JA)-1) THEN
                        LABV(IPSYM(1)) = NTYPE(3, JA) + I - 1
                        LABV(IPSYM(2)) = NTYPE(3, JB) + K - 1
                        IF (LABTVFRST.OR.LCPOT) THEN
                          CALL SETCOF_NXY_CSFG(M, N, VCOEFF)
                        ELSE
                          CALL twobody_DC(TEGRAL)
                          EMTBLOCK(M,N) = EMTBLOCK(M,N) + VCOEFF * TEGRAL
                        ENDIF 
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO  
            ENDDO

          ELSEIF (NSYMCR .EQ. 3) THEN
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                 'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***3*** 4'

          ELSEIF (NSYMCR .EQ. 4) THEN
            IF (MTYPE /= 1) THEN
            ! In rangular_mpi_csfg: savetv.f90, V-Coefficients with 
            ! NSYMCR = 4 and MTYPE = 2, 3 are discarded.
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                   'IC=', JASAV, 'IR=', JBSAV, 'LABV=', LABV(1:5)
              STOP 'Warning!!! NSYMCR ERROR ***4*** 4'
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
                    ! Upper-triangle matrix needed
                    IF (K.GT.I) CYCLE
                    IF (K == I .AND. L .GT. J) CYCLE

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
    2   CONTINUE
      ENDDO
!
! Add the common parts for the diagonal matrixelement
!
      IF (LCHM.AND.DABS(ENONSYM) .GT. CUTOFF) THEN
        DO I = 1, NTYPE(2,JA)
          EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM
        ENDDO
      ENDIF

      IF (LCHM.AND.LTRANSFER) THEN
        IF (LTRANSPOSE) EMTBLOCK=TRANSPOSE(EMTBLOCK)
        CALL TRANSFER_CSFG(JASAV+NCFGPAST, JBSAV+NCFGPAST)
      ENDIF

      RETURN
      END SUBROUTINE
