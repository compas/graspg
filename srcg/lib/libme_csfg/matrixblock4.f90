!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXBLOCK4(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,      &
                 NCORE,ELSTO)
!                                                                      *
!   This subroutine calls onescalar and computes one electron          *
!   matrix elements when IC = IR are of type 4.                        *
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
      USE stat_C
      USE mpi_C
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
      REAL(DOUBLE) :: EMTTMP,ATWINV,ENONSYM,TCOEFF,VCOEFF,TCOEF
      REAL(DOUBLE) :: TSAVE1,TSAVE2
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER :: IA,IB,IV,I,ICPTIJ,ICPTI,ICPTJ,IRPTK,IRPTKL,IRPTL,J,K,L
      INTEGER :: LABVTMP,M,MTYPE,N,NORB,NORBCOL,NORBROW,NSYMCOL0
      INTEGER :: IORBIF(3,4), IRW
!-----------------------------------------------------------------------
!      WRITE(*,*) 'ONESCALAR In matrixblock12 IC,IR',IC,IR

! Check, it should be diagnoal here, i.e., IRTOT.eq.ICTOT
      IF (ICTOT .NE. IRTOT) THEN
        WRITE(*,*)'IC=',IC,' IR=',IR,' ICTOT=',ICTOT,' IRTOT=',IRTOT
        STOP 'Error, ICTOT .NE. IRTOT in matrixblock4 ...'
      ENDIF

      ATWINV = 1.D0/EMN
      IBUG1 = 0

!   Set all matrix elements to zero
      !EMTBLOCK = 0.0D0

! Number of symmetry-ordered orbs for the IC- and IR-th symmetry-ordered-CSFs
      NORBCOL = NTYPE(6,IC) - NTYPE(3,IC) + 1
      NORBROW = NORBCOL

      !IF (IR.eq.121) THEN
      !  WRITE(*,*)'IR, MAP(IR)=', IR, MAP(IR)
      !  WRITE(*,'(250(I3))')(I, I=1,NW)
      !  WRITE(*,'(250(I3))')(IQA(I,IR),I=1,NW)
      !  WRITE(*,'(250(I3))')(JQSA(I,1,IR),I=1,NW)
      !  WRITE(*,'(250(I3))')(JQSA(I,2,IR),I=1,NW)
      !  WRITE(*,'(250(I3))')(JQSA(I,3,IR),I=1,NW)
      !  WRITE(*,'(250(I3))')(JCUPA(I,IR),I=1,NW)
      !ENDIF
! CALL RKCO_GG (IC, IR-5, BREID, 1, 2)
! < K L | h_br | L J >: K = I, L < J
! CALL RKCO_GG (IR-1, IR-2, BREID, 1, 2)
! < K L | h_br | L J >: K < L = I < J
! CALL RKCO_GG (IC, IR-2, BREID, 1, 2)
! < K L | h_br | L J >: K != I, L = J
! CALL RKCO_GG (IC, IR-1, BREID, 1, 2)
! < K L | h_br | L J >: K = I, and L = J
! CALL RKCO_GG (IC, IR, BREID, 1, 2)

! Needing IR-5, IR-2, IR-1 FICTIOUS CSFs:
      IF (NORBROW.GE.4) THEN
        ! 6s 7s (IR/IC, Here IR =  IC) ==> 4s 5s (IR-5)
        IORBIF(1,1) = IRFICT
        IORBIF(2,1) = NTYPE(3,IR)
        IORBIF(3,1) = NTYPE(5,IR) 
        CALL FICTIOUS_CSF(4, IRFICT, IR, NTYPE(4,IR)-2, NTYPE(4,IR),  &
                                         NTYPE(6,IR)-2, NTYPE(6,IR))
        ! CALL PRINT
      ENDIF
      IF (NORBROW.GE.3) THEN
        ! 6s 7s (IR, Here IR =  IC) ==> 5s 7s (IR-1)
        IORBIF(1,2) = IRFICT + 1
        IORBIF(2,2) = NTYPE(4,IR) - 1
        IORBIF(3,2) = NTYPE(6,IR)
        CALL FICTIOUS_CSF(4, IRFICT+1, IR, NTYPE(4,IR)-1, NTYPE(4,IR), &
                                                     0,           0)
        ! 6s 7s (IR, Here IR =  IC) ==> 5s 6s (IR-2)
        IORBIF(1,3) = IRFICT + 2
        IORBIF(2,3) = NTYPE(4,IR) - 1
        IORBIF(3,3) = NTYPE(6,IR) - 1
        CALL FICTIOUS_CSF(4, IRFICT+2, IR, NTYPE(4,IR)-1, NTYPE(4,IR), &
                                           NTYPE(6,IR)-1, NTYPE(6,IR))
      ENDIF
      ! The original IR/IC symmetry-ordered-CSF obtained from <state>.g
      IORBIF(1,4) = IR
      IORBIF(2,4) = NTYPE(4,IR)
      IORBIF(3,4) = NTYPE(6,IR)
!
!   Call onescalar
      DO MTYPE = 1, 2
       IA = 0
       TSHELL = 0.D0
       TCOEFF = 0.0D0
       IF (MTYPE.EQ.1) THEN
         ! < (Core A) K L (IR) | H1 | (Core A) I J (IC) >
         ! K = I, L  = J         : 5s 7s (IR) -- 5s 7s (IC)
         ! K = I, L != J         : 5s 6s (IR) -- 5s 7s (IC)
         ! K < I, L  = J         : 5s 7s (IR) -- 6s 7s (IC)
         CALL ONESCALAR(IC,IR,IA,IB,TSHELL)

       ELSEIF (MTYPE.EQ.2) THEN
         ! <IR-2| h1 | IC  > for K < L = I < J 
         !  5s 6s (IR) -- 6s 7s (IC)
         IF (NORBROW.LT.3) CYCLE  
         CALL ONESCALAR(IC, IRFICT+2, IA, IB, TSHELL)
       ENDIF

! Check 
!   
!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
!
!   Type 4,
!
! Part1: diagonal, Core orbital contribution
       IF (MTYPE.EQ.1) THEN
         ENONSYM = 0.0D0
         DO IA = 1,NW
           TCOEFF = DBLE(TSHELL(IA))
           IF (DABS (TCOEFF) .GT. CUTOFF) THEN
             IF (IA.GT.NORBGEN) THEN ! IA > NORBGEB --> IA symmetry-ordered orb
               EXIT ! TCOEFF recorded
             ELSE 
               ! Core orbital contribution
               CALL onebody_DC_MS_VP(IA,IA,ATWINV,TCOEFF,EMTTMP)
               ENONSYM = ENONSYM + EMTTMP
             END IF
           END IF
         ENDDO
       ELSE 
         ! MTYPE 2, record TCOEFF
         TCOEFF = DBLE(TSHELL(1))
         IF (IA.EQ.0 .OR. DABS(TCOEFF).LT.CUTOFF) CYCLE
       ENDIF 

! Obtain onebody contributions, add those from the symmetry-ordered orbitals to
! ENONSYM. Triangle matrix needed.  
       NCOEC = NCOEC + NORB
       N = 0
       DO I = 1, NORBCOL-1
        DO J = I + 1, NORBCOL
         N = N + 1  ! COLUMN INDEX
         M = 0
         DO K = 1, NORBROW-1
           DO L = K + 1, NORBROW
            M = M + 1 ! ROW INDEX
            ! TRIANGLE MATRIX NEEDED
            IF (K .GT. I) CYCLE
            IF (K .EQ. I .AND. L.GT.J) CYCLE

            IF (MTYPE.EQ.1) THEN
              IF (K.EQ.I) THEN
                ! For K = I,  L <= J
                IA = L + NTYPE(3,IR) - 1
                IB = J + NTYPE(3,IC) - 1
                CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
                EMTBLOCK(M,N) = EMTTMP

                !DIAGONAL MATRIXELEMENT IN EMTBLOCK
                IF (L.EQ.J) THEN
                  ! K = I, L = J
                  IA = K + NTYPE(3,IR) - 1
                  IB = I + NTYPE(3,IC) - 1
                  CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP + ENONSYM
                ENDIF
              ELSEIF (L.EQ.J) THEN
                 ! K.NE.I .AND. L.EQ.I
                 ! K < I, L = J
                 IA = K + NTYPE(3,IR) - 1
                 IB = I + NTYPE(3,IC) - 1 
                 CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
                 EMTBLOCK(M,N) = EMTTMP
              ENDIF
            ! MTYPE 2
            ELSEIF (L.EQ.I) THEN
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
      ENDDO

!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation 

!   ENONSYM that does not involve the symmetry-ordered orbital. 
      ENONSYM = 0.D0
      ! In fact, call RKCO_GG three times.
      DO MTYPE = 1, 5 
       NVCOEF = 0
       COEFF = 0.D0

! To label <K L| h1 | I J> of the CSFs generated by (core) 10s 11s
! MTYPE = 1 : <IR  | h1 | IC  > for K = I, L = J: 10s 11s -- 10s 11s
! MTYPE = 2 : <IR-2| h1 | IC-1> for K = I, L!= J:  9s 10s --  9s 11s
! MTYPE = 3 : <IR  | h1 | IC-2> for K = J, L!= I: 10s 11s --  9s 10s
! MTYPE = 4 : <IR-2| h1 | IC  > for L = I, K!= J:  9s 10s -- 10s 11s
! MTYPE = 5 : <IR-1| h1 | IC  > for L = J, K!= I:  9s 11s -- 10s 11s

       IF (MTYPE.EQ.1) THEN
         CALL RKCO_GG (IC, IR, CORD, INCOR, 1)       ! CALL 1
       ELSEIF (MTYPE.EQ.2) THEN
         ! K = I, L!= J
         ! Calculated within MTYPE = 1.
         CYCLE !Use the coefficients from MTYPE=1.
       ELSEIF (MTYPE.EQ.3) THEN
         ! Triangle matrix needed.  
         ! HERE K is never equal to J, because K <= I < J 
         CYCLE
       ELSEIF (MTYPE.EQ.4) THEN
         ! K < L = I < J 
         IF (NORBROW.LT.3) CYCLE
         !CALL RKCO_GG (IC, IR-2, CORD, INCOR, 1)  
         CALL RKCO_GG (IC, IRFICT+2, CORD, INCOR, 1) ! CALL 2 
       ELSEIF (MTYPE.EQ.5) THEN
         ! K < L = J > I, K /= I 
         IF (NORBROW.LT.3) CYCLE
         !CALL RKCO_GG (IC, IR-1, CORD, INCOR, 1) 
         CALL RKCO_GG (IC, IRFICT+1, CORD, INCOR, 1) ! CALL 3
       ENDIF

       DO IV = 1, NVCOEF
        VCOEFF = COEFF(IV)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC) 

! Determine the position of the symmetry-ordered orbitals
          LABV=LABEL(1:6, IV)
          CALL ANALABV(IC, IR, IV)

          IF (NSYMCR .EQ. 0) THEN 
            ! Common parts for diagonal matrixelement
            IF (MTYPE.NE.1) THEN
             WRITE(*,*)'IC=',IC,IC*3+5,' IR=',IR,IR*3+5,' BLOCK 4 ...'
             STOP 'Unexpected MTYPE.NE.1' 
            ENDIF
            CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
            ENONSYM = ENONSYM + EMTTMP

          ELSEIF (NSYMCR .EQ. 1 ) THEN 
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***1***'

          ELSEIF (NSYMCR .EQ. 2) THEN
           N = 0
           DO I = 1, NORBCOL-1
            DO J = I + 1, NORBCOL
             N = N + 1
             M = 0
             DO K = 1, NORBROW-1
              DO L = K + 1, NORBROW
               M = M + 1 
               ! TRIANGLE MATRIX NEEDED
               IF (K .GT. I ) CYCLE
               IF (K .EQ. I .AND. L.GT.J) CYCLE

               IF (K .EQ. I) THEN 
                IF (MTYPE .NE. 1) CYCLE  ! Using <IR |h12| IC>
                IF (L .EQ. J) THEN
                 ! DIAGONAL MATRIXELEMENT IN EMTBLOCK
                 ! INNER symmetry-ordered orb - CORE PAIR
                 IF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,IC) - 1) THEN
                  LABV(IPSYM(1)) = NTYPE(3, IC) + I - 1
                  LABV(IPSYM(2)) = NTYPE(3, IR) + K - 1

                 ! OUTER symmetry-ordered orb - CORE PAIR
                 ELSEIF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,IC)) THEN
                  LABV(IPSYM(1)) = NTYPE(3, IC) + J - 1
                  LABV(IPSYM(2)) = NTYPE(3, IR) + L - 1
                 ENDIF

                 CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP

                ! K.EQ.I .AND. L.NE.J
                ELSEIF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,IC)) THEN
                 LABV(IPSYM(1)) = NTYPE(3, IC) + J - 1
                 LABV(IPSYM(2)) = NTYPE(3, IR) + L - 1
                 CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                ENDIF

               ELSE ! K.NE.I, DIFFERENT INNER ORB-PAIR
                ! BUT THE OUTER ORB (L) OF IR IS SAME AS THE INNER 
                ! (I) or OUTER (J) orbitals OF IC 
                IF (MTYPE.EQ.4 .AND. L.EQ.I) THEN 
                 ! K < L = I < J
                 ! 9s10s (IR) - 10s11s (IC)
                 IF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,IC)) THEN
                  LABV(IPSYM(1)) = NTYPE(3, IC) + J - 1
                  LABV(IPSYM(2)) = NTYPE(3, IR) + K - 1
                  CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                 ENDIF
                ELSEIF (MTYPE.EQ.5 .AND. L.EQ.J) THEN 
                 ! K < L = J > I, K /= I 
                 ! PAIR AS 9s11s (IR) -- 10s 11s (IC)
                 IF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,IC) - 1) THEN
                  LABV(IPSYM(1)) = NTYPE(3, IC) + I - 1
                  LABV(IPSYM(2)) = NTYPE(3, IR) + K - 1
                  CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                 ENDIF
                ENDIF
               ENDIF

               ! Note, here K is always less than J, because J > I >= K. 

              ENDDO
             ENDDO
            ENDDO
           ENDDO

          ELSEIF (NSYMCR .EQ. 3) THEN
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***3***'

          ELSEIF (NSYMCR.EQ.4 .AND. MTYPE.EQ.1) THEN 
           ! VCOEFF ARE SAME FOR ALL possible nl n'l - n''l n'''l pairs
           N=0 
           DO I = 1, NORBCOL-1
            DO J = I + 1, NORBCOL 
             N = N + 1
             M = 0
             DO K = 1, NORBROW-1
              DO L = K + 1, NORBROW
               M = M + 1
               ! Triangle matrix needed
               IF (K .GT. I) CYCLE ! IR INNER ALWAYS <= IC INNER ORB
               IF (K .EQ. I .AND. L.GT.J) CYCLE  ! IF INNER IS SAME,
               LABV(1) = NTYPE(3, IC) + I - 1
               LABV(2) = NTYPE(3, IC) + J - 1
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
! Add the common parts for the diagonal matrixelement 
      IF (ABS(ENONSYM) .GT. CUTOFF ) THEN
        DO I = 1, NTYPE(2, IC)
           EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM
        ENDDO
      ENDIF

! Off-diagonal matrixelements have been obtained as above, assuming
! their V-Coefficients are identical to those of NSYMCR=4 in IC-IC pair.

      IBUG1 = 0
!
!   Accumulate the contribution from the two-electron
!   transverse interaction operator
      IF (LTRANS .AND. (INC2.EQ.1)) THEN !Kai:What do these arugments mean?
       ENONSYM = 0.D0
       DO MTYPE = 1, 5
        NVCOEF = 0
        COEFF = 0
        IF (MTYPE.EQ.1) THEN
         ! < (Core) 9s 10s (IC) | h_br | (Core) 7s 8s > 
         ! < K L (IR) | h_br | L J >: K != I, L != I, and, L != J
         IF (NORBCOL.LT.4) CYCLE ! At least 4 orbitals in the same sym.  
         !CALL RKCO_GG (IC, IR-5, BREID, 1, 2)
         CALL RKCO_GG (IC, IORBIF(1,1), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.2) THEN
         ! < (Core) 8s 10s (IC) | h_br | (Core ) 8s 9s >
         ! < K L (IR) | h_br | L J >: K = I, L < J
         IF (NORBCOL.LT.3) CYCLE
         !CALL RKCO_GG (IC-1, IR-2, BREID, 1, 2)
         CALL RKCO_GG (IORBIF(1,2), IORBIF(1,3), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.3) THEN
         ! < (Core) 9s 10s (IC) | h_br | (Core ) 8s 9s >
         ! < K L (IR) | h_br | L J >: K < L = I < J
         IF (NORBCOL.LT.3) CYCLE
         !CALL RKCO_GG (IC, IR-2, BREID, 1, 2)
         CALL RKCO_GG (IC, IORBIF(1,3), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.4) THEN
         ! < (Core) 9s 10s (IC) | h_br | (Core) 8s 10s >
         ! < K L (IR) | h_br | L J >: K != I, L = J
         IF (NORBCOL.LT.3) CYCLE
         !CALL RKCO_GG (IC, IR-1, BREID, 1, 2)
         CALL RKCO_GG (IC, IORBIF(1,2), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.5) THEN
         ! < (Core) 9s 10s (IC) | h_br | (Core) 9s 10s >
         ! < K L (IR) | h_br | L J >: K = I, and L = J
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
! Common for all diagonal matrixelement
             IF (MTYPE.NE.5) THEN
               WRITE(*,*)'IC,IR,IV,MTYPE=',IC,IR,IV,MTYPE
               STOP 'Unexpected MTYPE.NE.5 in matrixblock4.f ...'
             ENDIF
             CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
             ENONSYM = ENONSYM + EMTTMP

            ELSEIF (NSYMCR .EQ. 1) THEN
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')           &
                   'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
              STOP 'Warning!!! NSYMCR ERROR ***1***'

            ELSEIF (NSYMCR .EQ. 2) THEN
             IF (MTYPE.EQ.1) THEN
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')           &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
              STOP 'Unexpected MTYPE.EQ.1 in matrixblock4.f ...'
             ENDIF
! Determine the positions of the symmetry-ordered orbitals
             IRPTKL = 0
             ICPTIJ = 0
             IF (MTYPE.EQ.2 .OR. MTYPE.EQ.3) THEN
               ! MTYPE 2:  K = I, L < J,  <IC-1 : IR-2 >
               ! MTYPE 3: K < L = I < J,  <IC   : IR-2 >
               IF (LABEL(IPSym(1),IV).EQ.NTYPE(6,IC)) THEN
                 ICPTIJ = IPSym(1)
                 IRPTKL = IPSym(2)
               ELSE
                 ICPTIJ = IPSym(2)
                 IRPTKL = IPSym(1)
               ENDIF
             ELSEIF (MTYPE.EQ.4) THEN
               ! < K L (IR) | h_br | L J >: K != I, L = J 
               IF (LABEL(IPSym(1),IV).EQ.NTYPE(6,IC)-1) THEN
                 ICPTIJ = IPSym(1)
                 IRPTKL = IPSym(2)
               ELSE
                 ICPTIJ = IPSym(2)
                 IRPTKL = IPSym(1)
               ENDIF
             ENDIF
! Loop over symmetry-ordered orbitals 
             N = 0  ! Column index
             DO I = 1, NORBCOL-1
              DO J = I + 1, NORBCOL
                N = N + 1
                M = 0
                DO K = 1, NORBROW - 1
                 DO L = K + 1, NORBROW
                  M = M + 1  ! Row index
                  ! Up-triangle matrix needed
                  IF (K .GT. I) CYCLE ! IR INNER ALWAYS <= IC INNER ORB
                  IF (K .EQ. I .AND. L.GT.J) CYCLE 

                  IF (MTYPE.EQ.2 .AND. K.EQ.I .AND. L.LT.J) THEN
                    ! < K L | h_br | L J >: K = I,  L < J
                    LABV(IRPTKL) = NTYPE(3, IR) + L - 1
                    LABV(ICPTIJ) = NTYPE(3, IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ELSEIF (MTYPE.EQ.3 .AND. L.EQ.I) THEN ! K < L = I < J
                    ! < K L | h_br | L J >: K < L = I < J
                    LABV(IRPTKL) = NTYPE(3, IR) + K - 1
                    LABV(ICPTIJ) = NTYPE(3, IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ELSEIF (MTYPE.EQ.4 .AND. L.EQ.J .AND. K.NE.I) THEN  
                    ! < K L | h_br | L J >: K != I, L = J
                    LABV(IRPTKL) = NTYPE(3, IR) + K - 1
                    LABV(ICPTIJ) = NTYPE(3, IC) + I - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ELSEIF (MTYPE.EQ.5 .AND. K.EQ.I .AND. L.EQ.J) THEN
                    ! < K L | h_br | L J >: K = I, and L = J
                    IF (LABEL(IPSym(1),IV).EQ.NTYPE(6,IC)-1) THEN
                     LABV(IPSym(1)) = NTYPE(3, IR) + K - 1
                     LABV(IPSym(2)) = NTYPE(3, IC) + I - 1
                    ELSE
                     LABV(IPSym(1)) = NTYPE(3, IR) + L - 1
                     LABV(IPSym(2)) = NTYPE(3, IC) + J - 1
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
             IF (MTYPE.EQ.1) THEN
               DO I = 1, 4
                 IF (LABEL(I,IV).EQ.NTYPE(6,IR)-3) IRPTK = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IR)-2) IRPTL = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IC)-1) ICPTI = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IC)-0) ICPTJ = I
               ENDDO
             ENDIF
 
             IF (MTYPE.EQ.2) THEN
               ICPTI = 0
               DO I = 1, 4
                 IF (LABEL(I,IV).EQ.NTYPE(6,IR)-1) IRPTL = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IC)) ICPTJ = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IC)-2) THEN
                   IF (ICPTI.EQ.0) THEN
                     ICPTI = I
                   ELSE
                     IRPTK = I
                   ENDIF 
                 ENDIF
               ENDDO 
             ENDIF

             IF (MTYPE.EQ.3) THEN
               ICPTI = 0
               DO I = 1, 4
                 IF (LABEL(I,IV).EQ.NTYPE(6,IR)-2) IRPTK = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IC)) ICPTJ = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IC)-1) THEN
                   IF (ICPTI.EQ.0) THEN
                     ICPTI = I
                   ELSE
                     IRPTL = I
                   ENDIF
                 ENDIF
               ENDDO
             ENDIF

             IF (MTYPE.EQ.4) THEN
               IRPTL = 0
               DO I = 1, 4
                 IF (LABEL(I,IV).EQ.NTYPE(6,IR)-2) IRPTK = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IC)-1) ICPTI = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IC)) THEN
                   IF (IRPTL.EQ.0) THEN
                     IRPTL = I
                   ELSE
                     ICPTJ = I
                   ENDIF
                 ENDIF
               ENDDO
             ENDIF

             IF (MTYPE.EQ.5) THEN
              IRPTL = 0
              ICPTI = 0
              DO I = 1, 4
               IF (LABEL(I,IV).EQ.NTYPE(6,IC)-1) THEN
                 IF (ICPTI.EQ.0) THEN
                   ICPTI = I
                 ELSE
                   IRPTK = I
                 ENDIF
               ENDIF
               IF (LABEL(I,IV).EQ.NTYPE(6,IC)) THEN
                 IF (IRPTL.EQ.0) THEN
                   IRPTL = I
                 ELSE
                   ICPTJ = I
                 ENDIF
               ENDIF
              ENDDO
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
                  ! Up-triangle matrix needed
                  IF (K .GT. I) CYCLE ! IR INNER ALWAYS <= IC INNER ORB
                  IF (K .EQ. I .AND. L.GT.J) CYCLE

                  IF (MTYPE.EQ.1.AND.K.NE.I.AND.L.NE.I.AND.L.NE.J)   &
                  THEN 
                    ! < K L | h_br | L J >: K != I, L != I, L != J
                    LABV(IRPTK) = NTYPE(3,IR) + K - 1
                    LABV(IRPTL) = NTYPE(3,IR) + L - 1
                    LABV(ICPTI) = NTYPE(3,IC) + I - 1
                    LABV(ICPTJ) = NTYPE(3,IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ELSEIF (MTYPE.EQ.2 .AND. K.EQ.I .AND. L.NE.J) THEN
                    ! < K L | h_br | L J >: K = I, L != J
                    LABV(IRPTK) = NTYPE(3,IR) + K - 1
                    LABV(IRPTL) = NTYPE(3,IR) + L - 1
                    LABV(ICPTI) = NTYPE(3,IC) + I - 1
                    LABV(ICPTJ) = NTYPE(3,IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ELSEIF (MTYPE.EQ.3 .AND. L.EQ.I) THEN
                    ! < K L | h_br | L J >: K < L = I < J
                    LABV(IRPTK) = NTYPE(3,IR) + K - 1
                    LABV(IRPTL) = NTYPE(3,IR) + L - 1
                    LABV(ICPTI) = NTYPE(3,IC) + I - 1
                    LABV(ICPTJ) = NTYPE(3,IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ELSEIF (MTYPE.EQ.4 .AND. K.NE.I .AND. L.EQ.J) THEN
                    ! < K L | h_br | L J >: K != I, L = J
                    LABV(IRPTK) = NTYPE(3,IR) + K - 1
                    LABV(IRPTL) = NTYPE(3,IR) + L - 1
                    LABV(ICPTI) = NTYPE(3,IC) + I - 1
                    LABV(ICPTJ) = NTYPE(3,IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ELSEIF (MTYPE.EQ.5 .AND. K.EQ.I .AND. L.EQ.J) THEN
                    ! < K L | h_br | L J >: K = I, and L = J
                    LABV(IRPTK) = NTYPE(3,IR) + K - 1
                    LABV(IRPTL) = NTYPE(3,IR) + L - 1
                    LABV(ICPTI) = NTYPE(3,IC) + I - 1
                    LABV(ICPTJ) = NTYPE(3,IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ENDIF
                 ENDDO
                ENDDO
              ENDDO
             ENDDO
            ENDIF
          ENDIF
   10   CONTINUE
       ENDDO

! Add the common part of the diagonal matrixelements
       IF (DABS(ENONSYM) .GT. CUTOFF) THEN
        DO I = 1, NTYPE(2,IC)
          EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM
        ENDDO
       ENDIF 

      ENDIF

!cyc  Transfer EMTBLOCK to EMTSYM 
      call transfer_cyc(IC, IR)

      RETURN
      END SUBROUTINE MATRIXBLOCK4
