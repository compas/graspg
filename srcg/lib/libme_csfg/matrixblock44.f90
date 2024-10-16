!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXBLOCK44(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,     &
                 NCORE,ELSTO)
!                                                                      *
!   This subroutine calls onescalar and computes one electron          *
!   matrix elements when IC = IR are of type 4.                        *
!                                                                      *
!   Written by CHONG-YANG CHEN                               JUNE 2020 *
!   Last modification by C. Y. Chen                          Dec  2023 *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP
      USE symmatrix_mod
      USE buffer_C,   ONLY: NVCOEF, LABEL, COEFF
      USE debug_C,    ONLY: IBUG1
      USE decide_C
      USE def_C
      USE orb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
      EXTERNAL BREID,CORD
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,NCORE
      REAL(DOUBLE)        :: ELSTO
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-20
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: FLAGU
      REAL(DOUBLE) :: EMTTMP,ATWINV,ENONSYM,TCOEFF,VCOEFF,TCOEF
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER :: I,IA,IB,IV,ISWAP,IC0,ICPTI,ICPTJ,ICW,IR0,IRPTK,IRPTL,K
      INTEGER :: IRW,J,L,LABVTMP,M,MTYPE,N,NORB,NORBCOL,NORBROW
      INTEGER :: IORBMAX, IORBMIN, NUMORB
      INTEGER :: ICORBIF(3,4), IRORBIF(3,4)
!-----------------------------------------------------------------------
!     PRINT THE ONE-body (T-) coefficients, two-body electron and Breit
!     interaction (V-) coefficients. For clean, output only if FLAGU =.TRUE.
      !LPRINT = .TRUE.

      ATWINV = 1.D0/EMN
      IBUG1 = 0
      IC0 = IC
      IR0 = IR
!   Set all matrix elements to zero
      !EMTBLOCK = 0.0D0

! Number of symmetry-ordered orbs for the IC- and IR-th symmetry-ordered-CSFs
      NORBCOL = NTYPE(6,IC) - NTYPE(3,IC) + 1
      NORBROW = NTYPE(6,IR) - NTYPE(3,IR) + 1

      !NO ONEBODY CONTRIBUTION
      FLAGU = .FALSE.

      IF (NTYPE(3,IR).NE.NTYPE(3,IC)) GOTO 101
      FLAGU = .TRUE.

      NUMORB  = MAX(NORBCOL,NORBROW)
      IORBMAX = MAX(NTYPE(6,IR),NTYPE(6,IC))
      IORBMIN = NTYPE(3,IC)
       
! The smallest list version:
! Step 1: make IR / IC having the same n1 l n2 l pair
!         label them as IR0 / IC0 (Origianl the largest list version).
      ! Here, if IORBMIN = NTYPE(4,IR), then IORBMIN+1 = NTYPE(6,IR)
      ! Here, if IORBMIN = NTYPE(4,IC), then IORBMIN+1 = NTYPE(6,IC)
      ! Then, IR / IC are just copyed into IRFICT / ICFICT.
      ! fictious_csf.f90 can handled this kind of case. 
      !IF (IC.EQ.2756 .AND. IR.EQ.163) CALL PRINTFICTCSF(IC, IR, IC)
      !IF (IC.EQ.2756 .AND. IR.EQ.163) CALL PRINTFICTCSF(IC, IR, IR)
      CALL FICTIOUS_CSF(4, IRFICT, IR, IORBMIN, NTYPE(4,IR), &
                                     IORBMIN+1, NTYPE(6,IR))
      CALL FICTIOUS_CSF(4, ICFICT, IC, IORBMIN, NTYPE(4,IC), &
                                     IORBMIN+1, NTYPE(6,IC))
      IRORBIF(1,1) = IRFICT
      IRORBIF(2,1) = IORBMIN
      IRORBIF(3,1) = IORBMIN+1

      ICORBIF(1,1) = ICFICT
      ICORBIF(2,1) = IORBMIN
      ICORBIF(3,1) = IORBMIN+1
      !IF (IC.EQ.2756 .AND. IR.EQ.163)CALL PRINTFICTCSF(IC, IR, IRFICT)
      !IF (IC.EQ.2756 .AND. IR.EQ.163)CALL PRINTFICTCSF(IC, IR, ICFICT)
      
! Step 2: From IRFICT / ICFICT, construct other needed FICTIOUS CSFs.
! FLEXIBLE: VV-, CV-, CC-CSFs use different orbitals: (Breit interaction)
!  MTYPE = 2 : IR0 + 1 -- IC0                 :    K = I, L!= J;     4s 6s (IR1) -- 4s 5s (IC0)
!  MTYPE = 3 : IR0 + NORBROW - 1 -- IC0       :    K = J, L!= I;     5s 6s (IR2) -- 4s 5s (IC0)
!  MTYPE = 4 : IR0  -- IC0 + NORBCOL - 1      :    L = I, K!= J;     4s 5s (IR0) -- 5s 6s (IC2)
!  MTYPE = 5 : IR0 + NORBROW - 1 -- IC0 + 1   :    L = J, K!= I;     5s 6s (IR2) -- 4s 6s (IC1)
!  MTYPE = 6 : IR0    --  IC0                 :    K = I, L = J;     4s 5s (IR0) -- 4s 5s (IC0)
      IF (NUMORB.GE.3) THEN
        ! IR1 / IC1 4s5s ==> 4s6s
        CALL FICTIOUS_CSF(4, IRFICT+1, &
                             IRFICT, IORBMIN+2, IORBMIN+1, 0, 0)
        CALL FICTIOUS_CSF(4, ICFICT+1, &
                             ICFICT, IORBMIN+2, IORBMIN+1, 0, 0)
        !IF (IC.EQ.2756 .AND. IR.EQ.163)CALL PRINTFICTCSF(IC, IR, IRFICT+1)
        !IF (IC.EQ.2756 .AND. IR.EQ.163)CALL PRINTFICTCSF(IC, IR, ICFICT+1)
        IRORBIF(1,2) = IRFICT+1
        IRORBIF(2,2) = IORBMIN
        IRORBIF(3,2) = IORBMIN+2

        ICORBIF(1,2) = ICFICT+1
        ICORBIF(2,2) = IORBMIN
        ICORBIF(3,2) = IORBMIN+2

        ! IR2 / IC2 4s5s ==> 5s6s
        CALL FICTIOUS_CSF(4, IRFICT+2, IRFICT, IORBMIN+1, IORBMIN , &
                                               IORBMIN+2, IORBMIN+1)
        CALL FICTIOUS_CSF(4, ICFICT+2, ICFICT, IORBMIN+1, IORBMIN , &
                                               IORBMIN+2, IORBMIN+1)
        !IF (IC.EQ.2756 .AND. IR.EQ.163)CALL PRINTFICTCSF(IC, IR, IRFICT+2)
        !IF (IC.EQ.2756 .AND. IR.EQ.163)CALL PRINTFICTCSF(IC, IR, ICFICT+2)
        IRORBIF(1,3) = IRFICT+2
        IRORBIF(2,3) = IORBMIN+1
        IRORBIF(3,3) = IORBMIN+2

        ICORBIF(1,3) = ICFICT+2
        ICORBIF(2,3) = IORBMIN+1
        ICORBIF(3,3) = IORBMIN+2
      ENDIF

      IF (NUMORB.GE.4) THEN
!  MTYPE = 1 : K!=I,K!=J,L!=I,L!=J; 6s 7s (IR3) -- 4s 5s (IC0)
!  IR3  4s5s ==> 6s7s
        CALL FICTIOUS_CSF(4, IRFICT+3, IRFICT, IORBMIN+2, IORBMIN , &
                                               IORBMIN+3, IORBMIN+1)
        !IF (IC.EQ.2756 .AND. IR.EQ.163)CALL PRINTFICTCSF(IC, IR, IRFICT+3)
        IRORBIF(1,4) = IRFICT+3
        IRORBIF(2,4) = IORBMIN+2
        IRORBIF(3,4) = IORBMIN+3
      ENDIF

      IF (LPRINT) & 
        CALL PRINTLABELV(0, IC, IR, 44, 0, IC0, IR0, 0.d0)

!   Call onescalar
      ENONSYM = 0.0D0
      DO MTYPE = 1, 2 
       TSHELL = 0.D0
       TCOEFF = 0.0D0
       IA = 0
       IB = 0
!   < K L (IR) | h1 | I J (IC) >
       IF (MTYPE.EQ.1) THEN
        ! IR0    --  IC0 :  
        ! K = I, L = J;  4s 5s -- 4s 5s 
        ICW = ICORBIF(1,1)
        IRW = IRORBIF(1,1)
        CALL ONESCALAR(ICW, IRW, IA,IB,TSHELL)
        IF (IA == 0 .OR. ABS(TSHELL(1)).LT.CUTOFF) CYCLE
        TCOEFF = DBLE(TSHELL(1))
        IF (LPRINT) & 
          CALL PRINTLABELV(1, IC0, IR0, 44, MTYPE, IA, IB, TSHELL(1))

       ELSEIF (MTYPE.EQ.2) THEN
        ! L > K = J > I,  5s 6s -- 4s 5s
        ! K < L = I < J,  4s 5s -- 5s 6s
        IF (NORBROW.LT.3) CYCLE 
        ICW = ICORBIF(1,1)
        IRW = IRORBIF(1,3)
        CALL ONESCALAR(ICW, IRW, IA,IB,TSHELL)
        IF (IA == 0 .OR. ABS(TSHELL(1)).LT.CUTOFF) CYCLE
        TCOEFF = DBLE(TSHELL(1))
        IF (LPRINT) & 
          CALL PRINTLABELV(1, IC0, IR0, 44, MTYPE, IA, IB, TSHELL(1))
       ENDIF

!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
!
!   Type 4 - 4 
!
! Part1: semi-diagonal, Core orbital contribution
       IF (MTYPE.EQ.1) THEN
         ! K = I, L = J;  (Core A) 4s 5s -- (Core B) 4s 5s
         ! Here, IA /= 0 .AND. TCOEFF .GT.CUTOFF:
         ENONSYM = 0.0D0
         ! This contribution should arise from the labelling 
         ! orbital-pair.
         IF (IA.GT.NORBGEN.OR.IB.GT.NORBGEN) THEN ! IA > NORBGEB --> IA symmetry-ordered orb
           WRITE(*,*)'IC=',IC,IC*3+5,' IR=',IR,IR*3+5,' BLOCK44 ...'
           STOP 'Unexpected IA.GT.NORBGEN.OR.IB.GT.NORBGEN'
         ENDIF 
         IF (IA.GT.IB) THEN
           ISWAP = IA
           IA = IB
           IB = ISWAP
         ENDIF
         ! Core orbital contribution
         CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
         ENONSYM = EMTTMP
         
         ! ENONSYM is common to all Semi-diagonal matrixelements
         IF (NTYPE(2,IC).EQ.NTYPE(2,IR)) THEN
           DO I = 1, NTYPE(2,IC)
             EMTBLOCK(I,I) = ENONSYM            
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
               IF (K.EQ.I .AND. L.EQ.J) EMTBLOCK(M,N) = ENONSYM
              ENDDO
             ENDDO
            ENDDO
           ENDDO
         ENDIF 
       ELSEIF (MTYPE.EQ.2) THEN
         ! Here, IA /= 0 .AND. TCOEFF .GT.CUTOFF:
         ! Check:
         IF (IA.LT.NORBGEN .OR. IB.LT.NORBGEN) THEN
          ! This kind of contribution should arise from the orbital-pair
          ! above NORBGEN.
           WRITE(*,*)'IC=',IC,IC*3+5,' IR=',IR,IR*3+5,' BLOCK44 ...'
           STOP 'Unexpected IA.LT.NORBGEN .OR. IB.LT.NORBGEN'
         ENDIF
! Obtain onebody contributions for off-semi-diagonal matrixments.
         NCOEC = NCOEC + NORB
         N = 0
         DO I = 1, NORBCOL-1
          DO J = I + 1, NORBCOL
           N = N + 1  ! COLUMN INDEX
           M = 0
           DO K = 1, NORBROW-1
             DO L = K + 1, NORBROW
               M = M + 1 ! ROW INDEX

               IF (K.EQ.J) THEN
                ! K = J, L!= I; 5s 6s (K L) -- 4s 5s (I J)
                ! L > K = J > I  
                IA = I + NTYPE(3,IC) - 1
                IB = L + NTYPE(3,IR) - 1
                CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
                EMTBLOCK(M,N) = EMTTMP
               ENDIF

               IF (L.EQ.I) THEN
                ! L = I, K!= J; 4s 5s (K L) -- 5s 6s (I J)
                ! K < L = I < J
                IA = K + NTYPE(3,IR) - 1
                IB = J + NTYPE(3,IC) - 1 
                CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
                EMTBLOCK(M,N) = EMTTMP
               ENDIF

             ENDDO
           ENDDO
          ENDDO
         END DO
       ENDIF
      ENDDO  ! END MTYPE=1, 2

101   CONTINUE

!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation 

!   ENONSYM that does not involve the symmetry-ordered orbital. 
      ENONSYM = 0.D0
      ! Note: MTYPE = 1, 3, 4 are used. 
      DO MTYPE = 1, 5
       NVCOEF = 0

! FLEXIBLE: VV-, CV-, CC-CSFs use different orbitals: 
! Final implementation uses the followings: 
!  MTYPE = 1 : IRFICT   -- ICFICT   :  K = I, L = J;  4s 5s -- 4s 5s 
!  MTYPE = 2 : IRFICT+1 -- ICFICT   :  K = I, L!= J;  4s 6s -- 4s 5s 
!  MTYPE = 3 : IRFICT+2 -- ICFICT   :  K = J, L!= I;  5s 6s -- 4s 5s
!  MTYPE = 4 : IRFICT   -- ICFICT+2 :  L = I, K!= J;  4s 5s -- 5s 6s
!  MTYPE = 5 : IRFICT+2 -- ICFICT+1 :  L = J, K!= I;  5s 6s -- 4s 6s

!  The following case is treated within MTYPE 1.
!  IRFICT+3 -- ICFICT :   K!=I,K!=J,L!=I,L!=J; 6s 7s -- 4s 5s
!  After checks, MTYPE = 2 and 5 use also the same coefficients obtained
!  within MTYPE = 1 run.

       IF (MTYPE.EQ.1) THEN
!  MTYPE = 1 : IR0    --  IC0               :  K = I, L = J;  4s 5s -- 4s 5s 
         IF (.NOT.FLAGU) THEN
           ICW = IC
           IRW = IR
         ELSE
           ICW = ICORBIF(1,1)
           IRW = IRORBIF(1,1)
         ENDIF  
         CALL RKCO_GG (ICW, IRW, CORD, INCOR, 1)

       ELSEIF (MTYPE.EQ.2) THEN
!  MTYPE = 2 : IR0 + 1 -- IC0               :  K = I, L!= J;  4s 6s -- 4s 5s 
         IF (.NOT.FLAGU) EXIT
         ! Calculated within MTYPE = 1 loop
         CYCLE 

       ELSEIF (MTYPE.EQ.3) THEN
         ! L > K = J > I
         IF (NORBROW.LT.3) CYCLE
!  MTYPE = 3 : IR0 + NORBROW - 1 -- IC0     :  K = J, L!= I;  5s 6s -- 4s 5s
         ICW = ICORBIF(1,1)
         IRW = IRORBIF(1,3)
         CALL RKCO_GG (ICW, IRW, CORD, INCOR, 1)

       ELSEIF (MTYPE.EQ.4) THEN
!  MTYPE = 4 : IR0  -- IC0 + NORBCOL - 1    :  L = I, K!= J;  4s 5s -- 5s 6s
         ! K < L = I < J
         IF (NORBCOL.LT.3) CYCLE
         ICW = ICORBIF(1,3)
         IRW = IRORBIF(1,1)
         CALL RKCO_GG (ICW, IRW, CORD, INCOR, 1)

       ELSEIF (MTYPE.EQ.5) THEN
!  MTYPE = 5 : IR0 + NORBROW - 1 -- IC0 + 1 :  L = J, K!= I;  5s 6s -- 4s 6s
         CYCLE
         ! Calculated within MTYPE = 1 loop
       ENDIF

       IF (LPRINT .AND. FLAGU .AND. NVCOEF.GT.0) & 
          CALL PRINTLABELV(2, ICW, IRW, 44, MTYPE, IA, IB, TSHELL(1))
       IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
       DO IV = 1, NVCOEF
        VCOEFF = COEFF(IV)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC) 

! Determine the position of the symmetry-ordered orbitals
          LABV=LABEL(1:6, IV)
          CALL ANALABV(IC, IR, IV)

          IF (NSYMCR .EQ. 0) THEN 
            ! Common parts for SEMI-diagonal matrixelement
            IF (MTYPE.NE.1) THEN
             WRITE(*,*)'IC=',IC,IC*3+5,' IR=',IR,IR*3+5,' BLOCK 44 ...'
             STOP 'Unexpected NSYMCR .EQ. 0' 
            ENDIF
            CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
            ENONSYM = ENONSYM + EMTTMP
          
          ELSEIF (NSYMCR .EQ. 1 ) THEN 
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***1***'

          ELSEIF (NSYMCR .EQ. 2) THEN
           IF (NTYPE(3,IC).NE.NTYPE(3,IR)) THEN
            WRITE(*,*)'IC=',IC,IC*3+5,' IR=',IR,IR*3+5, 'Block 44'
            STOP 'Unexpected NTYPE(3,IC).NE.NTYPE(3,IR) ...'
           ENDIF
           N = 0
           DO I = 1, NORBCOL-1
            DO J = I + 1, NORBCOL
             N = N + 1
             M = 0
             DO K = 1, NORBROW-1
              DO L = K + 1, NORBROW
               M = M + 1 
               ! FULL MATRIX NEEDED
               ! IF (K .GT. I ) CYCLE
               ! IF (K .EQ. I .AND. L.GT.J) CYCLE

               IF (MTYPE.EQ.1) THEN
                IF (K.EQ.J) CYCLE ! MTYPE 3
                IF (L.EQ.I) CYCLE ! MTYPE 4

                ! K = I, L = J (Semi-diagonal) / L!= J (Off-diagonal, MTYPE 2)
                !IF (K.EQ.I.AND.LABEL(IPSYM(1),IV).EQ.NTYPE(6,ICW)) THEN
                IF (K.EQ.I.AND.LABEL(IPSYM(1),IV).EQ.IORBMIN+1) THEN
                ! < (Core A) K L | 1/r12 | (Core B) I J >: K and I
                ! orbitals are orthogonal normalized, leaving TWO
                ! non-symmetry-ordered orbitals and TWO symmetry-ordered orbitals
                ! involving L and J. 
                  LABV(IPSYM(1)) = NTYPE(3, IC) + J - 1
                  LABV(IPSYM(2)) = NTYPE(3, IR) + L - 1
                  CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                ENDIF
                ! L = J, K = I (Semi-diagonal) / K!= I (Off-diagonal, MTYPE 5)
                IF (L.EQ.J.AND.LABEL(IPSYM(1),IV).EQ.IORBMIN) THEN
                ! < (Core A) K L | 1/r12 | (Core B) I J >: L and J
                ! orbitals are orthogonal normalized, leaving TWO
                ! non-symmetry-ordered orbitals and TWO symmetry-ordered orbitals
                ! involving K and I. 
                  LABV(IPSYM(1)) = NTYPE(3, IC) + I - 1
                  LABV(IPSYM(2)) = NTYPE(3, IR) + K - 1
                  CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                ENDIF

               ELSEIF (MTYPE.EQ.3) THEN
                 IF (K.NE.J) CYCLE
               ! L > K = J > I
!  MTYPE = 3 : IR0 + NORBROW - 1 -- IC0     :  K = J, L!= I;  5s 6s -- 4s 5s
               ! IF  (LABEL(IPSYM(1),IV).EQ.NTYPE(4,ICW)) THEN
                ! < (Core A) K L | 1/r12 | (Core B) I J >: K and J
                ! orbitals are orthogonal normalized, leaving TWO
                ! non-symmetry-ordered orbitals and TWO symmetry-ordered orbitals
                ! involving L and I. 
                 LABV(IPSYM(1)) = NTYPE(3, IC) + I - 1
                 LABV(IPSYM(2)) = NTYPE(3, IR) + L - 1
                 CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
               ! ENDIF

               ELSEIF (MTYPE.EQ.4) THEN
                 IF (L.NE.I) CYCLE
               ! K < L = I < J
!  MTYPE = 4 : IR0  -- IC0 + NORBCOL - 1    :  L = I, K!= J;  4s 5s -- 5s 6s
               ! IF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,ICW)) THEN
                ! < (Core A) K L | 1/r12 | (Core B) I J >: L and I
                ! orbitals are orthogonal normalized, leaving TWO
                ! non-symmetry-ordered orbitals and TWO symmetry-ordered orbitals
                ! involving K and J. 
                 LABV(IPSYM(1)) = NTYPE(3, IC) + J - 1
                 LABV(IPSYM(2)) = NTYPE(3, IR) + K - 1
                 CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
               ! ENDIF
               ENDIF

               !IF (K .EQ. I) THEN 
               ! IF (MTYPE .NE. 1) CYCLE  ! Using <IR |h12| IC>
               !  ! OUTER symmetry-ordered ORB - CORE PAIR
               ! IF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,ICW)) THEN
               !   LABV(IPSYM(1)) = NTYPE(3, IC) + J - 1
               !   LABV(IPSYM(2)) = NTYPE(3, IR) + L - 1
               !   CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
               !   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
               ! ENDIF

               ! IF (L .EQ. J) THEN
               !  ! DIAGONAL MATRIXELEMENT IN EMTBLOCK
               !  ! INNER symmetry-ordered ORB - CORE PAIR
               !  IF (LABEL(IPSYM(1),IV).EQ.NTYPE(4,ICW)) THEN
               !   LABV(IPSYM(1)) = NTYPE(3, IC) + I - 1
               !   LABV(IPSYM(2)) = NTYPE(3, IR) + K - 1
               !   CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
               !   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
               !  ENDIF 
               ! ENDIF 

               !ELSE ! K.NE.I, DIFFERENT INNER ORB-PAIR
               ! ! Here K can be equal to J 
               ! IF (MTYPE.EQ.3 .AND. K.EQ.J) THEN
               !   ! 10s11s (IR) - 9s10s (IC)
               !   IF  (LABEL(IPSYM(1),IV).EQ.NTYPE(4,ICW)) THEN
               !    LABV(IPSYM(1)) = NTYPE(3, IC) + I - 1
               !    LABV(IPSYM(2)) = NTYPE(3, IR) + L - 1
               !    CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
               !    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP     
               !   ENDIF  
               ! ENDIF

               ! ! AND THE OUTER ORB (L) OF IR IS SAME AS THE INNER 
               ! ! (I) or OUTER (J) orbitals OF IC 
               ! IF (MTYPE.EQ.4. AND. L.EQ.I) THEN 
               !  ! 9s10s (IR) - 10s11s (IC)
               !  IF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,ICW)) THEN
               !   LABV(IPSYM(1)) = NTYPE(3, IC) + J - 1
               !   LABV(IPSYM(2)) = NTYPE(3, IR) + K - 1
               !   CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
               !   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
               !  ENDIF
               ! ENDIF  

               ! IF (MTYPE.EQ.1 .AND. L.EQ.J) THEN 
               !  ! PAIR AS 9s11s (IR) -- 10s 11s (IC)
               !  IF (LABEL(IPSYM(1),IV).EQ.NTYPE(4,ICW)) THEN
               !   LABV(IPSYM(1)) = NTYPE(3, IC) + I - 1
               !   LABV(IPSYM(2)) = NTYPE(3, IR) + K - 1
               !   CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
               !   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
               !  ENDIF
               ! ENDIF
               !ENDIF

              ENDDO
             ENDDO
            ENDDO
           ENDDO

          ELSEIF (NSYMCR .EQ. 3) THEN
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***3***'

          ELSEIF (NSYMCR.EQ.4 .AND. MTYPE.EQ.1) THEN 
           ! Contribution arsing from the integral among FOUR
           ! symmetry-ordered-orbitals. 
           ! VCOEFF ARE SAME FOR ALL possible nl n'l - n''l n'''l pairs
           ! It means the two-electron interaction involving four 
           ! symmetry-ordered-orbitals haveing the same coeffients no matter 
           ! MTYPE=2,3,4,5,...
           N=0 
           DO I = 1, NORBCOL-1
            LABV(1) = NTYPE(3, IC) + I - 1
            DO J = I + 1, NORBCOL 
             LABV(2) = NTYPE(3, IC) + J - 1
             N = N + 1
             M = 0
             DO K = 1, NORBROW-1
              DO L = K + 1, NORBROW
               M = M + 1
               ! Full matrix needed
               !IF (K .GT. I) CYCLE ! IR INNER ALWAYS <= IC INNER ORB
               !IF (K .EQ. I .AND. L.GT.J) CYCLE  ! IF INNER IS SAME,
               LABV(3) = NTYPE(3, IR) + K - 1
               LABV(4) = NTYPE(3, IR) + L - 1
               IF (LABEL(3,IV).GT.LABEL(4,IV)) THEN
                  LABVTMP = LABV(3)
                  LABV(3) = LABV(4)
                  LABV(4) = LABVTMP
               ENDIF
               CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
               EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
              ENDDO
             ENDDO
            ENDDO
           ENDDO

          ENDIF 
        ENDIF
       ENDDO

      ENDDO 
! common parts for the diagonal matrixelement 
      IF (DABS(ENONSYM) .GT. CUTOFF) THEN
       IF (NTYPE(3,IC).NE.NTYPE(3,IR)) THEN
        WRITE(*,*)'IC=',IC,IC*3+5,' IR=',IR,IR*3+5, ' Block 44 '
        STOP 'Unexpected NTYPE(3,IC).NE.NTYPE(3,IR) ...'
       ENDIF 

       IF (NTYPE(2,IC).EQ.NTYPE(2,IR)) THEN      
        DO I = 1, NTYPE(2, IC)
          EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM
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
            IF (K.EQ.I .AND. L.EQ.J) EMTBLOCK(M,N)=EMTBLOCK(M,N)+ENONSYM
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF
      ENDIF

      IBUG1 = 0
!   Accumulate the contribution from the two-electron
!   transverse interaction operator

! FLEXIBLE: VV-, CV-, CC-CSFs use different orbitals:
!  MTYPE = 1 : IR0 + 2*(NORBROW-2) + 1 -- IC0 : K!=I,K!=J,L!=I,L!=J; 6s 7s (IR3) -- 4s 5s (IC0)
!  MTYPE = 2 : IR0 + 1 -- IC0                 :    K = I, L!= J;     4s 6s (IR1) -- 4s 5s (IC0)
!  MTYPE = 3 : IR0 + NORBROW - 1 -- IC0       :    K = J, L!= I;     5s 6s (IR2) -- 4s 5s (IC0)
!  MTYPE = 4 : IR0  -- IC0 + NORBCOL - 1      :    L = I, K!= J;     4s 5s (IR0) -- 5s 6s (IC2)
!  MTYPE = 5 : IR0 + NORBROW - 1 -- IC0 + 1   :    L = J, K!= I;     5s 6s (IR2) -- 4s 6s (IC1)
!  MTYPE = 6 : IR0    --  IC0                 :    K = I, L = J;     4s 5s (IR0) -- 4s 5s (IC0)

! iorbmax = max(NTYPE(6,IC),NTYPE(6,IR)
! iorbmin = NTYPE(3,IC) = NTYPE(3,IR)

      IF (LTRANS .AND. (INC2.EQ.1)) THEN 
       ENONSYM = 0.D0
       DO MTYPE = 1, 7
        NVCOEF = 0
        COEFF = 0
        IF (MTYPE.EQ.1) THEN
         IF (NTYPE(3,IC).NE.NTYPE(3,IR)) CYCLE
         ! < (Core A) 9s 10s | h_br | (Core B) 7s 8s > 
         ! < K L | h_br | I J >: K != I, K != J, L != I, and, L != J
         IF (NORBROW.LT.4) CYCLE ! At least 4 orbitals in the same sym.  
         !CALL RKCO_GG (IC, IR-5, BREID, 1, 2)
!  MTYPE = 1 : IR0 + 2*(NORBROW-2) + 1 -- IC0 : K!=I,K!=J,L!=I,L!=J; 6s 7s -- 4s 5s
!  Smallest version: 
!  6s7s :   IR0 + 2*(NORBROW-2) + 1 ==> IRORBIF(1,4), IRW = 4
!  4s5s :   IC0                     ==> ICORBIF(1,1), ICW = 1
         !ICW = IC0
         !IRW = IR0 + 2*(NORBROW-2) + 1
         !CALL RKCO_GG (ICW, IRW, BREID, 1, 2)
         ICW = 1
         IRW = 4
         CALL RKCO_GG (ICORBIF(1,ICW), IRORBIF(1,IRW), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.2) THEN
         IF (NTYPE(3,IC).NE.NTYPE(3,IR)) CYCLE
         ! < (Core A) 8s 10s | h_br | (Core B) 8s 9s >
         ! < K L | h_br | I J >: K = I, L != J
         IF (NORBROW.LT.3) CYCLE
         !CALL RKCO_GG (IC-1, IR-2, BREID, 1, 2)
!  MTYPE = 2 : IR0 + 1 -- IC0                 :  K = I, L!= J;  4s 6s --  4s 5s 
!  Smallest version: 
!  4s6s :   IR0 + 1                 ==> IRORBIF(1,2), IRW = 2
!  4s5s :   IC0                     ==> ICORBIF(1,1), ICW = 1
         !ICW = IC0
         !IRW = IR0 + 1
         !CALL RKCO_GG (ICW, IRW, BREID, 1, 2)
         ICW = 1
         IRW = 2
         CALL RKCO_GG (ICORBIF(1,ICW), IRORBIF(1,IRW), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.3) THEN
         IF (NTYPE(3,IC).NE.NTYPE(3,IR)) CYCLE
         ! < (Core A) 8s 9s | h_br | (Core B) 9s 10s >
         ! < K L | h_br | I J >: K = J, L != I
         IF (NORBROW.LT.3) CYCLE
         !CALL RKCO_GG (IC-2, IR, BREID, 1, 2)
!  MTYPE = 3 : IR0 + NORBROW - 1 -- IC0       :  K = J, L!= I;  5s 6s --  4s 5s
!  Smallest version: 
!  5s6s :   IR0 + NORBROW - 1       ==> IRORBIF(1,3), IRW = 3
!  4s5s :   IC0                     ==> ICORBIF(1,1), ICW = 1
         !ICW = IC0
         !IRW = IR0 + NORBROW - 1
         !CALL RKCO_GG (ICW, IRW, BREID, 1, 2)
         ICW = 1
         IRW = 3
         CALL RKCO_GG (ICORBIF(1,ICW), IRORBIF(1,IRW), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.4) THEN
         IF (NTYPE(3,IC).NE.NTYPE(3,IR)) CYCLE
         ! < (Core) 9s 10s | h_br | (Core ) 8s 9s >
         ! < K L | h_br | I J >: K < L = I < J
         IF (NORBCOL.LT.3) CYCLE
         !CALL RKCO_GG (IC, IR-2, BREID, 1, 2)
!  MTYPE = 4 : IR0  -- IC0 + NORBCOL - 1      :  L = I, K!= J;  4s 5s --  5s 6s
!  Smallest version: 
!  4s5s :   IR0                     ==> IRORBIF(1,1), IRW = 1
!  5s6s :   IC0 + NORBCOL - 1       ==> ICORBIF(1,3), ICW = 3
         !ICW = IC0 + NORBCOL - 1
         !IRW = IR0
         !CALL RKCO_GG (ICW, IRW, BREID, 1, 2)
         ICW = 3
         IRW = 1
         CALL RKCO_GG (ICORBIF(1,ICW), IRORBIF(1,IRW), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.5) THEN
         IF (NTYPE(3,IC).NE.NTYPE(3,IR)) CYCLE
         ! < (Core) 9s 10s | h_br | (Core) 8s 10s >
         ! < K L | h_br | I J >: K != I, L = J
         IF (NORBROW.LT.3) CYCLE
         IF (NORBCOL.LT.3) CYCLE
         !CALL RKCO_GG (IC, IR-1, BREID, 1, 2)
!  MTYPE = 5 : IR0 + NORBROW - 1 -- IC0 + 1   :  L = J, K!= I;  5s 6s --  4s 6s
!  Smallest version: 
!  5s6s :   IR0 + NORBROW - 1       ==> IRORBIF(1,3), IRW = 3
!  4s6s :   IC0 + 1                 ==> ICORBIF(1,2), ICW = 2
         !ICW = IC0 + 1
         !IRW = IR0 + NORBROW - 1 
         !CALL RKCO_GG (ICW, IRW, BREID, 1, 2)
         ICW = 2
         IRW = 3
         CALL RKCO_GG (ICORBIF(1,ICW), IRORBIF(1,IRW), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.6) THEN
         IF (NTYPE(3,IC).NE.NTYPE(3,IR)) CYCLE
         ! < (Core) 9s 10s | h_br | (Core) 9s 10s >
         ! < K L | h_br | I J >: K = I, and L = J
         !CALL RKCO_GG (IC, IR, BREID, 1, 2)
!  MTYPE = 6 : IR0    --  IC0                 :  K = I, L = J;  4s 5s --  4s 5s 
!  Smallest version: 
!  4s5s :   IR0                     ==> IRORBIF(1,1), IRW = 1
!  4s5s :   IC0                     ==> ICORBIF(1,1), ICW = 1
         !ICW = IC0
         !IRW = IR0
         !CALL RKCO_GG (ICW, IRW, BREID, 1, 2)
         ICW = 1
         IRW = 1
         CALL RKCO_GG (ICORBIF(1,ICW), IRORBIF(1,IRW), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.7) THEN
         IF (NTYPE(3,IC).EQ.NTYPE(3,IR)) CYCLE
         ! < (Core ) 9s 10s | h_br | (Core ) 9p 10p >
         ICW = IC
         IRW = IR
         CALL RKCO_GG (ICW, IRW, BREID, 1, 2)
        ENDIF
        IF (LPRINT .AND. FLAGU .AND. NVCOEF.GT.0) & 
          CALL PRINTLABELV(3, ICW, IRW, 44, MTYPE, IA, IB, TSHELL(1))
        IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
        DO 10 IV = 1, NVCOEF
          VCOEFF = COEFF(IV)
          IF (DABS (VCOEFF) .GT. CUTOFF) THEN
            NMCBP = NMCBP + NTYPE(2,IC)
! Determine the position of the symmetry-ordered orbitals
            LABV = LABEL(1:6,IV)
            CALL ANALABV(IC, IR, IV)

            IF (NSYMCR .EQ. 0) THEN
! Common for all semi-diagonal matrixelement
! All the TWO pair integrals < K L (IR) | h_br | I J > involving the four
! symmetry-ordered-orbitals are orthogonal normalized, leaving four non-symmetry-ordered
! orbitals, thus this kind of contributions are common to all
! semi-diagonal matrixelement.
             IF (MTYPE.NE.6) THEN
               WRITE(*,*)'IC,IR,IV,MTYPE=',IC,IR,IV,MTYPE
               STOP 'Unexpected MTYPE.NE.6 in matrixblock44.f ...'
             ENDIF
             CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
             ENONSYM = ENONSYM + EMTTMP

            ELSEIF (NSYMCR .EQ. 1) THEN
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')           &
                   'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
              STOP 'Warning!!! NSYMCR ERROR ***1***'

            ELSEIF (NSYMCR .EQ. 2) THEN
             IF (MTYPE.EQ.1 .OR. MTYPE.EQ.7) THEN
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')           &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
              STOP 'Unexpected MTYPE.EQ.1 .OR. MTYPE.EQ.7 in B44.f ...'
             ENDIF
             N = 0  ! Column index
             DO I = 1, NORBCOL-1
              DO J = I + 1, NORBCOL
                N = N + 1
                M = 0
                DO K = 1, NORBROW - 1
                 DO L = K + 1, NORBROW
                  M = M + 1  ! Row index
                  ! Full matrix needed
                  ! IF (K .GT. I) CYCLE ! IR INNER ALWAYS <= IC INNER ORB
                  ! IF (K .EQ. I .AND. L.GT.J) CYCLE 

                  IF (MTYPE.EQ.2 .AND. K.EQ.I .AND. L.NE.J) THEN
                   ! < K L | h_br | I J >: K = I,  L != J
                   ! < (Core A) K L | h_br | (Core B) I J >: K and I
                   ! orbitals are orthogonal normalized, leaving TWO
                   ! non-symmetry-ordered orbitals and TWO symmetry-ordered orbitals
                   ! involving L and J. 
                   !IF (LABEL(IPSym(1),IV).EQ.NTYPE(6,IRW)) THEN 
                   IF (LABEL(IPSym(1),IV).EQ.IRORBIF(3,IRW)) THEN 
                    LABV(IPSym(1)) = NTYPE(3, IR) + L - 1
                    LABV(IPSym(2)) = NTYPE(3, IC) + J - 1
                   ELSE
                    LABV(IPSym(2)) = NTYPE(3, IR) + L - 1
                    LABV(IPSym(1)) = NTYPE(3, IC) + J - 1
                   ENDIF  
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP

                  ELSEIF (MTYPE.EQ.3 .AND. K.EQ.J) THEN 
                   ! < K L | h_br | I J >: L > K  = J > I
                   ! < (Core A) K L | h_br | (Core B) I J >: K and J
                   ! orbitals are orthogonal normalized, leaving TWO
                   ! non-symmetry-ordered orbitals and TWO symmetry-ordered orbitals
                   ! involving L and I. 
                   !IF (LABEL(IPSym(1),IV).EQ.NTYPE(6,IRW)) THEN
                   IF (LABEL(IPSym(1),IV).EQ.IRORBIF(3,IRW)) THEN
                    LABV(IPSym(1)) = NTYPE(3, IR) + L - 1
                    LABV(IPSym(2)) = NTYPE(3, IC) + I - 1
                   ELSE
                    LABV(IPSym(2)) = NTYPE(3, IR) + L - 1
                    LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
                   ENDIF
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP

                  ELSEIF (MTYPE.EQ.4 .AND. L.EQ.I) THEN 
                   ! < K L | h_br | I J >: K < L = I < J
                   ! < (Core A) K L | h_br | (Core B) I J >: L and I
                   ! orbitals are orthogonal normalized, leaving TWO
                   ! non-symmetry-ordered orbitals and TWO symmetry-ordered orbitals
                   ! involving K and J. 
                   !IF (LABEL(IPSym(1),IV).EQ.NTYPE(4,IRW)) THEN
                   IF (LABEL(IPSym(1),IV).EQ.IRORBIF(2,IRW)) THEN
                    LABV(IPSym(1)) = NTYPE(3, IR) + K - 1
                    LABV(IPSym(2)) = NTYPE(3, IC) + J - 1
                   ELSE
                    LABV(IPSym(2)) = NTYPE(3, IR) + K - 1
                    LABV(IPSym(1)) = NTYPE(3, IC) + J - 1
                   ENDIF
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP

                  ELSEIF (MTYPE.EQ.5 .AND. K.NE.I .AND. L.EQ.J) THEN  
                   ! < K L | h_br | I J >: K != I, L = J
                   ! < (Core A) K L | h_br | (Core B) I J >: L and J
                   ! orbitals are orthogonal normalized, leaving TWO
                   ! non-symmetry-ordered orbitals and TWO symmetry-ordered orbitals
                   ! involving K and I. 
                   !IF (LABEL(IPSym(1),IV).EQ.NTYPE(4,IRW)) THEN
                   IF (LABEL(IPSym(1),IV).EQ.IRORBIF(2,IRW)) THEN
                    LABV(IPSym(1)) = NTYPE(3, IR) + K - 1
                    LABV(IPSym(2)) = NTYPE(3, IC) + I - 1
                   ELSE
                    LABV(IPSym(2)) = NTYPE(3, IR) + K - 1
                    LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
                   ENDIF
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP

                  ELSEIF (MTYPE.EQ.6 .AND. K.EQ.I .AND. L.EQ.J) THEN
                    ! < K L | h_br | I J >: K = I, and L = J
                   !IF (LABEL(IPSym(1),IV).EQ.NTYPE(4,ICW)) THEN
                   IF (LABEL(IPSym(1),IV).EQ.ICORBIF(2,ICW)) THEN
                   ! < (Core A) K L | h_br | (Core B) I J >: L and J
                   ! orbitals are orthogonal normalized, leaving TWO
                   ! non-symmetry-ordered orbitals and TWO symmetry-ordered orbitals
                   ! involving K and I. 
                    LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
                    LABV(IPSym(2)) = NTYPE(3, IR) + K - 1
                   ELSE
                   ! < (Core A) K L | h_br | (Core B) I J >: K and I
                   ! orbitals are orthogonal normalized, leaving TWO
                   ! non-symmetry-ordered orbitals and TWO symmetry-ordered orbitals
                   ! involving L and J. 
                    LABV(IPSym(1)) = NTYPE(3, IC) + J - 1
                    LABV(IPSym(2)) = NTYPE(3, IR) + L - 1
                   ENDIF
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ENDIF 
                 ENDDO
                ENDDO
              ENDDO
             ENDDO

            ELSEIF (NSYMCR .EQ. 3) THEN
             WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')            &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
             STOP 'Warning!!! NSYMCR ERROR ***3***'

            ELSEIF (NSYMCR .EQ. 4) THEN
! Determine the positions of the symmetry-ordered orbitals 
             IRPTK = 0
             IRPTL = 0
             ICPTI = 0
             ICPTJ = 0
             IF (MTYPE.EQ.1) THEN
              ! K != I, K != J, L != I, and, L != J
              DO I = 1, 4
               !IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) IRPTK = I
               !IF (LABEL(I,IV).EQ.NTYPE(6,IRW)) IRPTL = I
               !IF (LABEL(I,IV).EQ.NTYPE(4,ICW)) ICPTI = I
               !IF (LABEL(I,IV).EQ.NTYPE(6,ICW)) ICPTJ = I
               IF (LABEL(I,IV).EQ.IRORBIF(2,IRW)) IRPTK = I
               IF (LABEL(I,IV).EQ.IRORBIF(3,IRW)) IRPTL = I
               IF (LABEL(I,IV).EQ.ICORBIF(2,ICW)) ICPTI = I
               IF (LABEL(I,IV).EQ.ICORBIF(3,ICW)) ICPTJ = I
              ENDDO

             ELSEIF (MTYPE.EQ.2) THEN
              !  K = I, L != J
              DO I = 1, 4
               !IF (LABEL(I,IV).EQ.NTYPE(6,IRW)) IRPTL = I
               !IF (LABEL(I,IV).EQ.NTYPE(6,ICW)) ICPTJ = I
               !IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) THEN
               IF (LABEL(I,IV).EQ.IRORBIF(3,IRW)) IRPTL = I
               IF (LABEL(I,IV).EQ.ICORBIF(3,ICW)) ICPTJ = I
               IF (LABEL(I,IV).EQ.IRORBIF(2,IRW)) THEN
                 IF (IRPTK.EQ.0) THEN
                   IRPTK = I
                 ELSE
                   ICPTI = I
                 ENDIF
               ENDIF
              ENDDO

             ELSEIF (MTYPE.EQ.3) THEN
              ! K = J, L > K = J > I, hence L != I
              DO I = 1, 4
               !IF (LABEL(I,IV).EQ.NTYPE(6,IRW)) IRPTL = I
               !IF (LABEL(I,IV).EQ.NTYPE(4,ICW)) ICPTI = I
               !IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) THEN
               IF (LABEL(I,IV).EQ.IRORBIF(3,IRW)) IRPTL = I
               IF (LABEL(I,IV).EQ.ICORBIF(2,ICW)) ICPTI = I
               IF (LABEL(I,IV).EQ.IRORBIF(2,IRW)) THEN
                 IF (IRPTK.EQ.0) THEN
                   IRPTK = I
                 ELSE
                   ICPTJ = I
                 ENDIF
               ENDIF
              ENDDO
               
             ELSEIF (MTYPE.EQ.4) THEN
              ! L = I, K != J
              DO I = 1, 4
               !IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) IRPTK = I
               !IF (LABEL(I,IV).EQ.NTYPE(6,ICW)) ICPTJ = I
               !IF (LABEL(I,IV).EQ.NTYPE(6,IRW)) THEN
               IF (LABEL(I,IV).EQ.IRORBIF(2,IRW)) IRPTK = I
               IF (LABEL(I,IV).EQ.ICORBIF(3,ICW)) ICPTJ = I
               IF (LABEL(I,IV).EQ.IRORBIF(3,IRW)) THEN
                 IF (IRPTL.EQ.0) THEN
                   IRPTL = I
                 ELSE
                   ICPTI = I
                 ENDIF
               ENDIF
              ENDDO

             ELSEIF (MTYPE.EQ.5) THEN
              ! K != I, L = J
              DO I = 1, 4
               !IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) IRPTK = I
               !IF (LABEL(I,IV).EQ.NTYPE(4,ICW)) ICPTI = I
               !IF (LABEL(I,IV).EQ.NTYPE(6,IRW)) THEN
               IF (LABEL(I,IV).EQ.IRORBIF(2,IRW)) IRPTK = I
               IF (LABEL(I,IV).EQ.ICORBIF(2,ICW)) ICPTI = I
               IF (LABEL(I,IV).EQ.IRORBIF(3,IRW)) THEN
                 IF (IRPTL.EQ.0) THEN
                   IRPTL = I
                 ELSE
                   ICPTJ = I
                 ENDIF
               ENDIF
              ENDDO

             ELSEIF (MTYPE.EQ.6) THEN
              ! K = I, L = J
              DO I = 1, 4
               !IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) THEN
               IF (LABEL(I,IV).EQ.IRORBIF(2,IRW)) THEN
                 IF (IRPTK.EQ.0) THEN
                   IRPTK = I
                 ELSE
                   ICPTI = I
                 ENDIF 
               ENDIF
               !IF (LABEL(I,IV).EQ.NTYPE(6,IRW)) THEN
               IF (LABEL(I,IV).EQ.IRORBIF(3,IRW)) THEN
                 IF (IRPTL.EQ.0) THEN
                   IRPTL = I
                 ELSE
                   ICPTJ = I
                 ENDIF
               ENDIF 
              ENDDO

             ELSEIF (MTYPE.EQ.7) THEN
              DO I = 1, 4
               IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) IRPTK = I
               IF (LABEL(I,IV).EQ.NTYPE(6,IRW)) IRPTL = I
               IF (LABEL(I,IV).EQ.NTYPE(4,ICW)) ICPTI = I
               IF (LABEL(I,IV).EQ.NTYPE(6,ICW)) ICPTJ = I
              ENDDO
             ENDIF
! Check
             IF(IRPTK*IRPTL*ICPTI*ICPTJ.EQ.0) THEN
               WRITE(*,*)                                            &
                  'IC=',IC,' IR=',IR,' IV=',IV,' LABV=',LABV(1:5),   &
                  ' MTYPE=', MTYPE,                                  &
                  'IRPTK,IRPTL,ICPTI,ICPTJ=',IRPTK,IRPTL,ICPTI,ICPTJ
               STOP 'Unexpected IRPTK*IRPTL*ICPTI*ICPTJ.EQ.0 in b44.f .'
             ENDIF
 
! Loop over the symmetry-ordered orbitals
             N = 0  ! Column index
             DO I = 1, NORBCOL-1
              DO J = I + 1, NORBCOL
                N = N + 1
                M = 0
                DO K = 1, NORBROW - 1
                 DO L = K + 1, NORBROW
                  M = M + 1  ! Row index
                  ! Full matrix needed
                  ! IF (K .GT. I) CYCLE 
                  ! IF (K .EQ. I .AND. L.GT.J) CYCLE
!  MTYPE = 1 : IR0 + 2*(NORBROW-2) + 1 -- IC0 : K!=I,K!=J,L!=I,L!=J; 6s 7s -- 4s 5s
                  IF (MTYPE.EQ.1 .AND.                               &
                     (K.EQ.I .OR. K.EQ.J .OR. L.EQ.I .OR. L.EQ.J)) CYCLE

!  MTYPE = 2 : IR0 + 1 -- IC0                 :  K = I, L!= J;  4s 6s --  4s 5s 
                  IF (MTYPE.EQ.2 .AND. (K.NE.I .OR. L.EQ.J)) CYCLE

!  MTYPE = 3 : IR0 + NORBROW - 1 -- IC0       :  K = J, L!= I;  5s 6s --  4s 5s
!  If K.EQ.J, then L!=I because L > K = J > I
                  IF (MTYPE.EQ.3 .AND. (K.NE.J)) CYCLE

!  MTYPE = 4 : IR0  -- IC0 + NORBCOL - 1      :  L = I, K!= J;  4s 5s --  5s 6s
!  If L.EQ.I, then K!=J because K < L = I < J
                  IF (MTYPE.EQ.4 .AND. (L.NE.I)) CYCLE

!  MTYPE = 5 : IR0 + NORBROW - 1 -- IC0 + 1   :  L = J, K!= I;  5s 6s --  4s 6s
                  IF (MTYPE.EQ.5 .AND. (K.EQ.I .OR. L.NE.J)) CYCLE

!  MTYPE = 6 : IR0    --  IC0                 :  K = I, L = J;  4s 5s --  4s 5s 
                  IF (MTYPE.EQ.6 .AND. (K.NE.I .OR. L.NE.J)) CYCLE
                  LABV(IRPTK) = NTYPE(3,IR) + K - 1
                  LABV(IRPTL) = NTYPE(3,IR) + L - 1
                  LABV(ICPTI) = NTYPE(3,IC) + I - 1
                  LABV(ICPTJ) = NTYPE(3,IC) + J - 1
                  CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                 ENDDO
                ENDDO
              ENDDO
             ENDDO
            ENDIF
          ENDIF
   10   CONTINUE
       ENDDO

! Add the common part of the semi-diagonal matrixelements
       IF (ABS(ENONSYM) .GT. CUTOFF) THEN
        IF (NTYPE(3,IC).NE.NTYPE(3,IR)) THEN
          WRITE(*,*)'IC,IR,ENONSYM=',IC,IR,ENONSYM
          STOP 'Unexpected NTYPE(3,IC).NE.NTYPE(3,IR) in block44.f ...'
        ENDIF
        IF (NTYPE(2,IC).EQ.NTYPE(2,IR)) THEN
         DO I = 1, NTYPE(2,IC)
          EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM
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
             IF (K.EQ.I.AND.L.EQ.J) EMTBLOCK(M,N)=EMTBLOCK(M,N)+ENONSYM
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ENDIF 

      ENDIF

!cyc  Transfer EMTBLOCK to EMTSYM 
      IF (LTRANSFER) call transfer_cyc(IC, IR)

      RETURN
      END SUBROUTINE MATRIXBLOCK44
