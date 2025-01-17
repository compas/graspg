!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXBLOCK45(ICSAV,IRSAV,NCOEC,INCOR,NCTEC,INC2,     &
                 NMCBP,NCORE,ELSTO)
!                                                                      *
!   This subroutine calls onescalar and computes one electron          *
!   matrix elements when IC and IR are of type 4 and 5                 * 
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
      use symmatrix_mod
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
      INTEGER :: ICSAV,IRSAV,NCOEC,INCOR,NCTEC,INC2,NMCBP,NCORE
      REAL(DOUBLE)        :: ELSTO
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!   Matrix elements smaller than CUTOFF are not accumulated
!
!      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-20
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAl :: FLAGU
      REAL(DOUBLE) :: EMTTMP,ATWINV,ENONSYM,TCOEFF,VCOEFF
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER :: IA,IB,IV,ISWAP,I,IA0,IC,IC0,ICPTI,ICPTIJ,ICPTJ
      INTEGER :: ICW,IFLAG12,IR,IR0,IRPTK,IRPTKL,IRPTL,IRW,J,K,L,M,N
      INTEGER :: MTYPE,NORBCOL,NORBROW, ICORBIF(2,3)
!-----------------------------------------------------------------------
! A: 
!  1s ( 1)  2s ( 1) 10s ( 1) 11s ( 1)
!      1/2      1/2      1/2      1/2
!                    0      1/2      0+
! B:
!  1s ( 1)  2s ( 1) 10s ( 1) 11s ( 1)
!      1/2      1/2      1/2      1/2
!                    1      1/2      0+
! C:
!  1s ( 2) 11s ( 2)
!
!                  0+
! <A|h|C> differs from <B|h|C>
!

! Exchange IC and IR if needed, to make IC is TYPE 5, IR TYPE 4
      IF (LTRANSPOSE) THEN
         IC = IRSAV
         IR = ICSAV
      ! The first COLUMN generated by symmetry-ordered-CSF ICC
      ! The first ROW generated by symmetry-ordered-CSF IRR
      !   IC0 = IRTOT + 1
      !   IR0 = ICTOT + 1
      ELSE
         IC = ICSAV
         IR = IRSAV
      !   IC0 = ICTOT + 1
      !   IR0 = IRTOT + 1
      ENDIF
      IR0 = IR
 
      ATWINV = 1.D0/EMN
      IBUG1 = 0

!   Set all matrix elements to zero
      !EMTBLOCK = 0.0d0
      ENONSYM = 0.0D0
      FLAGU = .FALSE.
! Number of symmetry-ordered orbs for the IC- and IR-th symmetry-ordered-CSFs
      NORBCOL = NTYPE(4,IC) - NTYPE(3,IC) + 1 
      NORBROW = NTYPE(6,IR) - NTYPE(3,IR) + 1 

!   For type 4 and 5, there is onebody contributions when the
!   symmetry-ordered
!   orbitails are all same.
      IF (NTYPE(3,IR).NE.NTYPE(3,IC)) GOTO 101
      FLAGU = .TRUE.

      ! Buiding the FICTIOUTS CSF as needed.
      ! L > K = I  
      CALL FICTIOUS_CSF(5, ICFICT, IC, NTYPE(4,IR), NTYPE(4,IC), 0, 0)
      ICORBIF(1,1) = ICFICT
      ICORBIF(2,1) = NTYPE(4,IR)

      ! K < L = I  
      CALL FICTIOUS_CSF(5, ICFICT+1, IC, NTYPE(6,IR), NTYPE(4,IC), 0, 0)
      ICORBIF(1,2) = ICFICT+1
      ICORBIF(2,2) = NTYPE(6,IR)

      IF (NORBCOL.GE.3 .OR. NORBROW.GE.3) THEN
      ! K /= I AND L /= I   
        IF (NORBCOL.GT.NORBROW) THEN
          ICORBIF(1,3) = IC
          ICORBIF(2,3) = NTYPE(4,IC)
        ELSE
          CALL FICTIOUS_CSF(5, ICFICT+2, IC, NTYPE(3,IC),     & 
                                             NTYPE(4,IC), 0, 0)
          ICORBIF(1,3) = ICFICT+2
          ICORBIF(2,3) = NTYPE(3,IC)
        ENDIF
      ENDIF     

!   Call onescalar 
!   < K L (IR) | h1 | I J (IC) >
!   IR0: (Core A) 4s 5s -- IC0 (Core B) 4s2
      DO MTYPE = 1, 2
       TSHELL = 0.D0
       IF (MTYPE.EQ.1) THEN
        IF (NORBCOL.LT.2) CYCLE
        
        ! 4s 5s / 5s2, K != L = I = J
        ICW = ICORBIF(1,2)
        IRW = IR0  
        CALL ONESCALAR(ICW,IRW,IA,IB,TSHELL)
       ELSE
        ! 4s 5s / 4s2, L != K = I = J
        ICW = ICORBIF(1,1)
        IRW = IR0
        CALL ONESCALAR(ICW,IRW,IA,IB,TSHELL)
       ENDIF
!   
!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
!
       IF (IA .NE. 0) THEN
        TCOEFF = DBLE(TSHELL(1))
        LTRANSFER = .TRUE.
        IF (DABS(TCOEFF) .GT. CUTOFF) THEN
          NCOEC = NCOEC + NTYPE(2,IC)  
          N = 0
          DO I = 1, NORBCOL
            IA0 = I + NTYPE(3,IC) - 1
            N = N + 1
            M = 0
            DO K = 1,NORBROW-1
             DO L = K + 1, NORBROW
              M = M + 1 
              IF (K.NE.I .AND. L.NE.I) CYCLE

              IF (K.EQ.I .AND. MTYPE.EQ.2) THEN
                IB = L + NTYPE(3,IR) - 1
                IA = IA0  ! IA might be changed by the following statements.
                IF (IA.GT.IB) THEN
                 ISWAP = IA
                 IA = IB
                 IB = ISWAP
                ENDIF
                CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
                EMTBLOCK(M,N) = EMTTMP
              ENDIF

              IF (L.EQ.I .AND. MTYPE.EQ.1) THEN
                IB = K + NTYPE(3,IR) - 1
                IA = IA0  ! IA might be changed by the following statements.
                IF (IA.GT.IB) THEN
                 ISWAP = IA
                 IA = IB
                 IB = ISWAP
                ENDIF
                CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
                EMTBLOCK(M,N) = EMTTMP
              ENDIF

             ENDDO
            ENDDO
          END DO
        ENDIF
       ENDIF
      ENDDO

101   CONTINUE
      IBUG1 = 0
!
!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation; the latter
!   is computed first because the orbital indices may be
!   permuted by RKINTC
!
! One more RKCO_GG call to deal with the possible 'PHASE-FACTOR'.
      DO MTYPE = 1, 2 
       ENONSYM = 0.0D0
       NVCOEF = 0
       IF (MTYPE.EQ.1) THEN 
         ! 4s 5s -- 5s2
         ! K < L = I = J
         IF (FLAGU) THEN
           IF (NORBCOL.LT.2) CYCLE
           ICW = ICORBIF(1,2)
           IRW = IR0 
         ELSE
           ICW = IC
           IRW = IR
         ENDIF 
         CALL RKCO_GG (ICW, IRW,CORD, INCOR, 1)
       ELSE
         IF (.NOT.FLAGU) CYCLE
         ! 4s 5s -- 4s2
         ! L > K = I = J
         ICW = ICORBIF(1,1)
         IRW = IR
         CALL RKCO_GG (ICW, IRW,CORD, INCOR, 1)
       ENDIF
!
       IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
       DO IV = 1, NVCOEF
        VCOEFF = COEFF(IV)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC) 

! Determine the position of the symmetry-ordered orbitals
          LABV=LABEL(1:6, IV)
          CALL ANALABV(IC, IR, IV)

! THERE ARE AT MOST THREE symmetry-ordered-orbs FOR TYPE 2 - 4 MATRIXELEMENTS.
          IF (NSYMCR .EQ. 0) THEN 
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')             &
                 'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***0***'

          ELSEIF (NSYMCR .EQ. 1 ) THEN 
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')             &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***1***'

          ELSEIF (NSYMCR .EQ. 2 ) THEN
           IF (.NOT.FLAGU) THEN
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')             &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***2***'
           ENDIF

           N = 0
           DO I = 1, NORBCOL
            LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
            N = N + 1
            M = 0
            DO K = 1, NORBROW-1
              DO L = K + 1, NORBROW
                M = M + 1
                IF (K.NE.I .AND. L.NE.I) CYCLE

                ! < (Core A) K L | 1/r_12 | (Core B) I J >: K and I
                ! orbitals are orthogonal normalized, leaving TWO
                ! non-symmetry-ordered orbitals and TWO symmetry-ordered orbitals
                ! involving L and J (= I, TYPE 5)
                IF (K.EQ.I .AND. MTYPE.EQ.2) THEN
                 LABV(IPSym(2)) = NTYPE(3, IR) + L - 1
                 CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                ENDIF

                IF (L.EQ.I .AND. MTYPE.EQ.1) THEN
                ! < (Core A) K L | 1/r_12 | (Core B) I J >: L and J
                ! orbitals are orthogonal normalized, leaving TWO
                ! non-symmetry-ordered orbitals and TWO symmetry-ordered orbitals
                ! involving K and I (= J, TYPE 5)
                 LABV(IPSym(2)) = NTYPE(3, IR) + K - 1
                 CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                ENDIF

              ENDDO
            ENDDO
           ENDDO

          ELSEIF (NSYMCR .EQ. 3 ) THEN
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***3***'
          
          ELSEIF (NSYMCR .EQ. 4 ) THEN
           IF (MTYPE.NE.1) CYCLE

           N = 0
           DO I = 1, NORBCOL
            N = N + 1
            LABV(1) = NTYPE(3, IC) + I - 1
            LABV(2) = NTYPE(3, IC) + I - 1
            M = 0
            DO K = 1, NORBROW-1
              DO L = K + 1, NORBROW
                M = M + 1
                LABV(3) = NTYPE(3, IR) + K - 1
                LABV(4) = NTYPE(3, IR) + L - 1
                CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
              ENDDO
            ENDDO
           ENDDO

          ENDIF 
        ENDIF
       ENDDO
      ENDDO

!
      IBUG1 = 0
!   Accumulate the contribution from the two-electron
!   transverse interaction operator
      IF (LTRANS .AND. (INC2.EQ.1)) THEN 
       ENONSYM = 0.0d0
       DO MTYPE = 1, 3
        NVCOEF = 0
        IF (MTYPE.EQ.1) THEN
         IF (.NOT.FLAGU) CYCLE
         ! < (Core A) 9s 10s | h_br | (Core B) 8s2 >
         ! < K L | h_br | I J >: K < L != I = J 
         IF (NORBROW.LT.3 .AND. NORBCOL.LT.3) CYCLE ! At least 3 orbitals in the same sym.  
         ! K /= I AND L /= I
         IRW = IR
         ICW = 3
         CALL RKCO_GG (ICORBIF(1,ICW), IRW, BREID, 1, 2)

        ELSEIF (MTYPE.EQ.2) THEN
         IF (.NOT.FLAGU) CYCLE 
         ! L > K = I = J  
         ! < (Core)  9s 10s | h_br | (Core) 9s2 >
         ! 4s 5s -- 4s2
         ICW = 1
         IRW = IR 
         CALL RKCO_GG (ICORBIF(1,ICW), IRW, BREID, 1, 2)

        ELSEIF (MTYPE.EQ.3) THEN
         ! FLAGU: K < L = I = J  
         ! < (Core)  9s 10s | h_br | (Core) 10s2 >
         ! .NOT. FLAGU: 
         ! < (Core)  9s 10s | h_br | (Core) 10p2 >
         IF (FLAGU) THEN
           IF (NORBCOL.LT.2) CYCLE
           ! 4s 5s -- 5s2
           ! K < L = I = J 
           ICW = 2
           IRW = IR0
           CALL RKCO_GG (ICORBIF(1,ICW), IRW, BREID, 1, 2)
         ELSE
           CALL RKCO_GG (IC, IR, BREID, 1, 2)
         ENDIF
        ENDIF
        IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
        DO 10 IV = 1, NVCOEF
          VCOEFF = COEFF(IV)
          IF (DABS (VCOEFF) .GT. CUTOFF) THEN
            NMCBP = NMCBP + NTYPE(2,IC)
! Determine the position of the symmetry-ordered orbitals
            LABV = LABEL(1:6,IV)
            CALL ANALABV(IC, IR, IV)

            IF (NSYMCR .EQ. 0) THEN
! No common part for 4 - 5 matrixelement
             !IF (MTYPE.NE.6) THEN
               WRITE(*,*)'IC,IR,IV,MTYPE=',IC,IR,IV,MTYPE
               STOP 'Unexpected NSYMCR .EQ. 0 in matrixblock45.f ...'
             !ENDIF
             !CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
             !ENONSYM = ENONSYM + EMTTMP

            ELSEIF (NSYMCR .EQ. 1) THEN
             WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')            &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
             STOP 'Warning!!! NSYMCR ERROR ***1***'
             
! Matrixelement between Type 4 - 5, Loop for symmetry-ordered-orbitals
            ELSEIF (NSYMCR .EQ. 2) THEN
             IF (MTYPE.EQ.1 .OR. (.NOT.FLAGU)) THEN
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')           &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
              STOP 'Unexpected MTYPE.EQ.1 .OR. (.NOT.FLAGU) in B45.f '
             ENDIF

             IRPTKL = 0
             ICPTIJ = 0
             IFLAG12 = 0
             IF (MTYPE.EQ.2 .AND.                                    &
                   LABEL(IPSym(1),IV).NE.NTYPE(6,IRW)) IFLAG12 =1
             IF (MTYPE.EQ.3 .AND.                                    &
                   LABEL(IPSym(1),IV).NE.NTYPE(4,IRW)) IFLAG12 = 1
             IF (IFLAG12.EQ.0) THEN
               IRPTKL = IPSym(1)
               ICPTIJ = IPSym(2)
             ELSE
               IRPTKL = IPSym(2)
               ICPTIJ = IPSym(1)
             ENDIF

             N = 0  ! Column index
             DO I = 1, NORBCOL
                J = I
                N = I
                M = 0
                DO K = 1, NORBROW - 1
                 DO L = K + 1, NORBROW
                  M = M + 1  ! Row index
                  ! Full matrix needed
                  ! IF (K .GT. I) CYCLE ! IR INNER ALWAYS <= IC INNER ORB
                  ! IF (K .EQ. I .AND. L.GT.J) CYCLE 
                  IF (MTYPE.EQ.2 .AND. K.EQ.I) THEN
                    ! < K L | h_br | I J >: K = I
                    LABV(IRPTKL) = NTYPE(3, IR) + L - 1
                    LABV(ICPTIJ) = NTYPE(3, IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ELSEIF (MTYPE.EQ.3 .AND. L.EQ.I) THEN 
                    ! < K L | h_br | I J >: L = I, FLAGU
                    LABV(IRPTKL) = NTYPE(3, IR) + K - 1
                    LABV(ICPTIJ) = NTYPE(3, IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ENDIF 
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
         ! 4s 5s -- 6s2, or 5s 6s -- 4s2
         ! K /= I AND L /= I 
              DO I = 1, 4
               IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) IRPTK = I
               IF (LABEL(I,IV).EQ.NTYPE(6,IRW)) IRPTL = I
               IF (LABEL(I,IV).EQ.ICORBIF(2,ICW)) THEN
                 IF (ICPTI.EQ.0) THEN
                   ICPTI = I
                 ELSE
                   ICPTJ = I
                 ENDIF
               ENDIF
              ENDDO

             ELSEIF (MTYPE.EQ.2) THEN
         ! 4s 5s  -- 4s2
         ! L > K = I
              DO I = 1, 4 
                IF (LABEL(I,IV).EQ.NTYPE(6,IRW)) IRPTL = I
                IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) THEN
                  IF (IRPTK.EQ.0) THEN
                    IRPTK = I
                  ELSEIF (ICPTI.EQ.0) THEN
                    ICPTI = I
                  ELSE
                    ICPTJ = I
                  ENDIF
                ENDIF
              ENDDO

             ELSEIF (MTYPE.EQ.3) THEN
         ! FLAGU: K < L = I = J  
         ! < (Core)  9s 10s | h_br | (Core) 10s2 >
         ! .NOT. FLAGU: 
         ! < (Core)  9s 10s | h_br | (Core) 10p2 >
              IF (FLAGU) THEN
               DO I = 1, 4
                IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) IRPTK = I
                IF (LABEL(I,IV).EQ.NTYPE(6,IRW)) THEN
                  IF (IRPTL.EQ.0) THEN
                    IRPTL = I
                  ELSEIF (ICPTI.EQ.0) THEN
                    ICPTI = I
                  ELSE
                    ICPTJ = I
                  ENDIF
                ENDIF
               ENDDO
              ELSE
               DO I = 1, 4 
                IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) IRPTK = I
                IF (LABEL(I,IV).EQ.NTYPE(6,IRW)) IRPTL = I
                IF (LABEL(I,IV).EQ.NTYPE(4,ICW)) THEN
                  IF (ICPTI.EQ.0) THEN
                    ICPTI = I
                  ELSE
                    ICPTJ = I
                  ENDIF
                ENDIF
               ENDDO
              ENDIF
             ENDIF 
! Check
      IF (IRPTK*IRPTL*ICPTI*ICPTJ.EQ.0) THEN
        WRITE(*,*)'IC,IR,IV,IRPTK*IRPTL*ICPTI*ICPTJ= ',              &
        IC,IR,IV,IRPTK,IRPTL,ICPTI,ICPTJ
        STOP 'Unexpected IRPTK*IRPTL*ICPTI*ICPTJ.EQ.0...'
      ENDIF       
! Loop over the symmetry-ordered orbitals
             N = 0  ! Column index
             DO I = 1, NORBCOL
                J = I 
                N = I
                M = 0
                DO K = 1, NORBROW - 1
                 DO L = K + 1, NORBROW
                  M = M + 1  ! Row index

                  IF (MTYPE.EQ.1.AND.(K.EQ.I.OR.L.EQ.I)) CYCLE
                  IF (MTYPE.EQ.2.AND.K.NE.I) CYCLE
                  IF (MTYPE.EQ.3.AND.FLAGU.AND.L.NE.I) CYCLE
                  LABV(IRPTK) = NTYPE(3,IR) + K - 1
                  LABV(IRPTL) = NTYPE(3,IR) + L - 1
                  LABV(ICPTI) = NTYPE(3,IC) + I - 1
                  LABV(ICPTJ) = NTYPE(3,IC) + J - 1
                  CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                 ENDDO
                ENDDO
             ENDDO
            ENDIF
          ENDIF
   10   CONTINUE
       ENDDO
      ENDIF

!cyc  Transfer EMTBLOCK to EMTSYM 
      IF (LTRANSPOSE) EMTBLOCK=TRANSPOSE(EMTBLOCK)
      IF (LTRANSFER) Call transfer_cyc(ICSAV, IRSAV)

      RETURN
      END SUBROUTINE MATRIXBLOCK45

