!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXBLOCK5(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,      &
                 NCORE,ELSTO)
!                                                                      *
!   This subroutine calls onescalar and computes one electron          *
!   matrix elements when IC = IR are of type 5                         *
!                                                                      *
!   Written by CHONG-YANG CHEN                               JUNE 2020 *
!   Last modification by C. Y. Chen                          Dec  2022 *
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
!
!   Matrix elements smaller than CUTOFF are not accumulated
!
!      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-20
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: EMTTMP,ATWINV,ENONSYM,TCOEFF,VCOEFF,TCOEFFSAVE
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER :: IA,IB,IV,I,ICPTI,ICPTJ,IRPTK,IRPTL,J,K,L,M,MTYPE,N
      INTEGER :: NORBCOL,NORBROW
!-----------------------------------------------------------------------
! Check, it should be diagnoal here, i.e., IRTOT.eq.ICTOT
      IF (ICTOT .NE. IRTOT) THEN
        WRITE(*,*)'IC=',IC,' IR=',IR,' ICTOT=',ICTOT,' IRTOT=',IRTOT
        STOP 'Error, ICTOT .NE. IRTOT in matrixblock5 ...'
      ENDIF

      ATWINV = 1.D0/EMN
      IBUG1 = 0

!   Set all matrix elements to zero
      !EMTBLOCK = 0.0D0
      TSHELL = 0.D0
      ENONSYM = 0.0D0

      NORBCOL = NTYPE(4,IC) - NTYPE(3,IC) + 1
      NORBROW = NTYPE(4,IR) - NTYPE(3,IR) + 1

      ! Buiding the FICTIOUTS CSF as needed.
      IF (NORBCOL.GE.2) THEN
        CALL FICTIOUS_CSF(5, IRFICT, IR, NTYPE(3,IR),     &
                                         NTYPE(4,IR), 0, 0)
      ENDIF
!   Call onescalar        
      CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
!   
!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
!
      DO IA = 1,NW
        TCOEFF = DBLE(TSHELL(IA))
        IF (DABS (TCOEFF) .GT. CUTOFF) THEN
          IF (IA.GT.NORBGEN) THEN ! IA > NORBGEB --> IA symmetry-ordered orb
            TCOEFFSAVE = TCOEFF
          ELSE
            CALL onebody_DC_MS_VP(IA,IA,ATWINV,TCOEFF,EMTTMP)
            ENONSYM = ENONSYM + EMTTMP
          END IF
        END IF
      ENDDO

!  To get the diagonal matrix elements, add contributions from each of
!  the symmetry-ordered orbitals to ENONSYM 
      NCOEC = NCOEC + NORBCOL
      DO I = 1,NORBCOL
        IA = I + NTYPE(3,IC) - 1
        CALL onebody_DC_MS_VP(IA,IA,ATWINV,TCOEFFSAVE,EMTTMP)
        EMTBLOCK(I,I) = ENONSYM + EMTTMP
      END DO

!   No off-diagonal part within the diagonal symmetry-ordered block.

!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation 

!   ENONSYM that does not include the symmetry-ordered orbital. 
      NVCOEF = 0
      ENONSYM = 0.D0

      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
      DO I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC) 

! Determine the position of the symmetry-ordered orbitals
          LABV=LABEL(1:6, I)
          CALL ANALABV(IC, IR, I)

          IF (NSYMCR .EQ. 0) THEN 
            CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
            ENONSYM = ENONSYM + EMTTMP

          ELSEIF (NSYMCR .EQ. 1) THEN 
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***1***'

          ELSEIF (NSYMCR .EQ. 2) THEN
           DO J = NORBCOL, 1, -1
            LABV(IPSYM(1)) = NTYPE(3, IC) + J - 1
            LABV(IPSYM(2)) = NTYPE(3, IR) + J - 1
            CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
            EMTBLOCK(J,J) = EMTBLOCK(J,J) + EMTTMP
           ENDDO

          ELSEIF (NSYMCR .EQ. 3 ) THEN
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***3***'

          ELSEIF (NSYMCR .EQ. 4 ) THEN
           DO J = 1, NORBCOL
            DO K = 1, J ! TRIANGLE MATRIX NEEDED
              LABV(IPSYM(1)) = NTYPE(3, IC) + J - 1
              LABV(IPSYM(2)) = NTYPE(3, IC) + J - 1
              LABV(IPSYM(3)) = NTYPE(3, IC) + K - 1
              LABV(IPSYM(4)) = NTYPE(3, IC) + K - 1
              CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
              EMTBLOCK(K,J) = EMTBLOCK(K,J) + EMTTMP
            ENDDO 
           ENDDO 
          ENDIF 
        ENDIF
      ENDDO
! Add the common parts for the diagonal matrixelement 
      IF (DABS(ENONSYM) .GT. CUTOFF ) THEN
        DO I = 1, NORBCOL
           EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM
        ENDDO
      ENDIF

      IBUG1 = 0
!   Accumulate the contribution from the two-electron
!   transverse interaction operator
      IF (LTRANS .AND. (INC2.EQ.1)) THEN 
       ENONSYM = 0.0d0
       DO MTYPE = 1, 2
        NVCOEF = 0
        COEFF = 0
        IF (MTYPE.EQ.1) THEN
         ! < (Core A) 9s2 | h_br | (Core B) 10s2 >
         ! < K L | h_br | I J >: (K = L) != (I = J)
         IF (NORBROW.LT.2) CYCLE ! At least 2 orbitals in the same sym.  
         !CALL RKCO_GG (IC, IR-1, BREID, 1, 2)
         CALL RKCO_GG (IC, IRFICT, BREID, 1, 2)

        ELSEIF (MTYPE.EQ.2) THEN
         ! < (Core)  10s2   | h_br | (Core) 10s2 >
         CALL RKCO_GG (IC, IR, BREID, 1, 2)

        ENDIF

        DO 10 IV = 1, NVCOEF
          VCOEFF = COEFF(IV)
          IF (DABS (VCOEFF) .GT. CUTOFF) THEN
            NMCBP = NMCBP + NTYPE(2,IC)
! Determine the position of the symmetry-ordered orbitals
            LABV = LABEL(1:6,IV)
            CALL ANALABV(IC, IR, IV)

            IF (NSYMCR .EQ. 0) THEN
! No common part for 4 - 5 matrixelement
             IF (MTYPE.NE.2) THEN
               WRITE(*,*)'IC,IR,IV,MTYPE=',IC,IR,IV,MTYPE
               STOP 'Unexpected NSYMCR .EQ. 0 in matrixblock5.f ...'
             ENDIF
             CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
             ENONSYM = ENONSYM + EMTTMP

! Matrixelement between Type 3 - 4, Loop for symmetry-ordered-orbitals
            ELSEIF (NSYMCR .EQ. 1) THEN
             WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')            &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
             STOP 'Warning!!! NSYMCR ERROR ***1***'
             
            ELSEIF (NSYMCR .EQ. 2) THEN
             IF (MTYPE.EQ.1) THEN
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')           &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
              STOP 'Unexpected MTYPE.EQ.1 in B5.f '
             ENDIF

             N = 0  ! Column index
             DO I = 1, NORBCOL
                J = I
                N = I
                M = 0
                ! Diagonal matrixelements
                DO K = I, I 
                  L = K 
                  M = K  ! Row index
                  ! IF (K .GT. I) CYCLE ! IR INNER ALWAYS <= IC INNER ORB
                  ! IF (K .EQ. I .AND. L.GT.J) CYCLE 
                  ! IF (K.EQ.I) THEN
                    ! < K L | h_br | I J >: K = I
                    LABV(IPSYM(1)) = NTYPE(3, IC) + K - 1
                    LABV(IPSYM(2)) = NTYPE(3, IC) + I - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ! ENDIF 
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
         ! < (Core A) 9s2    | h_br | (Core B) 10s2 >
              DO I = 1, 4
               !IF (LABEL(I,IV).EQ.NTYPE(4,IR)-1) THEN
               IF (LABEL(I,IV).EQ.NTYPE(3,IR)) THEN
                 IF (IRPTK.EQ.0) THEN
                   IRPTK = I
                 ELSE
                   IRPTL = I
                 ENDIF
               ENDIF
               IF (LABEL(I,IV).EQ.NTYPE(4,IC)-0) THEN
                 IF (ICPTI.EQ.0) THEN
                   ICPTI = I
                 ELSE
                   ICPTJ = I
                 ENDIF
               ENDIF
              ENDDO

             ELSEIF (MTYPE.EQ.2) THEN
         ! < (Core)  10s2   | h_br | (Core)10s2 >
         ! I, J, K, AND L have the same values 
              IRPTK = 1
              IRPTL = 2
              ICPTI = 3
              ICPTJ = 4 
             ENDIF 
! Check
!      IF (IRPTK*IRPTL*ICPTI*ICPTJ.EQ.0) THEN
!        WRITE(*,*)'IC,IR,IV,IRPTK*IRPTL*ICPTI*ICPTJ= ',              &
!        IC,IR,IV,IRPTK,IRPTL,ICPTI,ICPTJ
!        STOP 'Unexpected IRPTK*IRPTL*ICPTI*ICPTJ.EQ.0...'
!      ENDIF

! Loop over the symmetry-ordered orbitals
             N = 0  ! Column index
             DO I = 1, NORBCOL
                J = I 
                N = I
                M = 0
                ! Up-triangle matrix needed
                DO K = 1, I
                 L = K
                  M = K      ! Row index
                  IF (MTYPE.EQ.1.AND.K.EQ.I) CYCLE
                  IF (MTYPE.EQ.2.AND.K.NE.I) CYCLE
                  LABV(IRPTK) = NTYPE(3,IR) + K - 1
                  LABV(IRPTL) = NTYPE(3,IR) + L - 1
                  LABV(ICPTI) = NTYPE(3,IC) + I - 1
                  LABV(ICPTJ) = NTYPE(3,IC) + J - 1
                  CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                ENDDO
             ENDDO
            ENDIF
          ENDIF
   10   CONTINUE
       ENDDO
! Add the common parts for diagonal matrixelement
       IF (DABS(ENONSYM) .GT. CUTOFF) THEN
        DO I = 1, NTYPE(2,IC)
          EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM  
        ENDDO
       ENDIF
      ENDIF


!cyc  Transfer EMTBLOCK to EMTSYM 
      call transfer_cyc(IC, IR)

      RETURN
      END SUBROUTINE MATRIXBLOCK5
