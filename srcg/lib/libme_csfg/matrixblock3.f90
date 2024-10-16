!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXBLOCK3(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,      &
                 NCORE,ELSTO)
!                                                                      *
!   This subroutine calls onescalar and computes one electron          *
!   matrix elements when IC and IR are of type 1 and 2                 *         
!                                                                      *
!   Written by Per Jönsson & Kai Wang                        May 2017  *
!                                                                      *
!   Modified by CHONG-YANG CHEN                              JUNE 2020 *
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
      REAL(DOUBLE) :: EMTTMP,ATWINV,ENONSYM,TCOEFF,VCOEFF
      REAL(DOUBLE) :: TCOEFFSAVE,TCOEFFSAVE1,TCOEFFSAVE2
      REAL(DOUBLE) :: TSAVE1,TSAVE2
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER :: IA,IB,IV,I,ICPTI,ICPTJ,IRPTK,IRPTL,J,K,L,LABVTMP,M
      INTEGER :: MTYPE,NORB,NORBCOLL,NORBCOLU,NORBROWL,NORBROWU,N
      INTEGER :: NSYMCOL0, IRORBIF(3,4)
!-----------------------------------------------------------------------
! Check, it should be diagnoal here, i.e., IRTOT.eq.ICTOT
      IF (ICTOT .NE. IRTOT) THEN
        WRITE(*,*)'IC=',IC,' IR=',IR,' ICTOT=',ICTOT,' IRTOT=',IRTOT
        STOP 'Error, ICTOT .NE. IRTOT in matrixblock3 ...'
      ENDIF

      ATWINV = 1.D0/EMN
      IBUG1 = 0

!   Set all matrix elements to zero
      !EMTBLOCK = 0.0D0
      TSHELL = 0.D0

! Number of symmetry-ordered orbs for the IC- and IR-th symmetry-ordered-CSFs
      NORBCOLL = NTYPE(4,IC) - NTYPE(3,IC) + 1
      NORBCOLU = NTYPE(6,IC) - NTYPE(5,IC) + 1
      NORBROWL = NTYPE(4,IR) - NTYPE(3,IR) + 1 
      NORBROWU = NTYPE(6,IR) - NTYPE(5,IR) + 1
!
!   Call onescalar
      TSHELL = 0.D0
      CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
! There must be onebody contributions.       
      IF (IA.EQ.0) THEN
        WRITE(*,*)'IC=',IC,' IR=',IR,' ICTOT=',ICTOT,' IRTOT=',IRTOT
        STOP 'Error, IA = 0 in matrixblock3.f ...'
      ENDIF
!   
!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
!
!   Type 3
!
! Part1: diagonal, common
      NSYMCOL0 = 0
      ENONSYM = 0.0D0
      TSAVE1 = 0.0D0
      TSAVE2 = 0.0D0
      DO IA = 1, NW
        TCOEFF = DBLE(TSHELL(IA))
        IF (DABS (TCOEFF) .GT. CUTOFF) THEN
          IF (IA.GT.NORBGEN) THEN ! IA > NORBGEB --> IA symmetry-ordered orb
            NSYMCOL0 = NSYMCOL0 + 1
            IF (NSYMCOL0 .EQ. 1) TSAVE1 = TCOEFF
            IF (NSYMCOL0 .EQ. 2) TSAVE2 = TCOEFF
          ELSE
            CALL onebody_DC_MS_VP(IA,IA,ATWINV,TCOEFF,EMTTMP)
            ENONSYM = ENONSYM + EMTTMP
          END IF
        END IF
      ENDDO
! CHECK, IT SHOULD MEET TSAVE1=TSAVE2=1.0D0
      IF(DABS(TSAVE1-TSAVE2).GT.1.0D-8.OR.DABS(TSAVE1-1.0D0).GT.1.0D-8)&
      THEN
        WRITE(*,*)'IC=',IC,' IR=',IR,' ICTOT=',ICTOT,' IRTOT=',IRTOT
        WRITE(*,*)'TCOEFFSAVE1=',TCOEFFSAVE1,' TCOEFFSAVE2=',TCOEFFSAVE2
        STOP 'Error, Something unexpected in matrixblock3 ...'
      ENDIF
      TCOEFFSAVE=TSAVE1

! To get the diagonal matrix elements, add contributions from each of
! the symmetry-ordered orbitals to ENONSYM 
! The off-diagonal elements are also calculated
      NCOEC = NCOEC + NORB
      N = 0
      DO I = 1, NORBCOLL
       DO J = 1, NORBCOLU
          N = N + 1  ! COLUMN INDEX
          M = 0
          DO K = 1, NORBROWL
           DO L = 1, NORBROWU
            M = M + 1 ! ROW INDEX
            ! TRIANGLE MATRIX NEEDED
            IF (K .GT. I) CYCLE
            IF (K .EQ. I .AND. L.GT.J) CYCLE 
            IF (K .EQ. I) THEN
             IA = L + NTYPE(5,IR) - 1
             IB = J + NTYPE(5,IC) - 1  
             CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFFSAVE,EMTTMP)
             EMTBLOCK(M,N) = EMTTMP

             IF (L .EQ. J) THEN
              !DIAGONAL MATRIXELEMENT IN EMTBLOCK
              IA = K + NTYPE(3,IR) - 1
              IB = I + NTYPE(3,IC) - 1  
              CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFFSAVE,EMTTMP)
              EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP + ENONSYM
             ENDIF
            ! K.NE.I, but L.EQ.J 
            ELSEIF (L .EQ. J) THEN
              IA = K + NTYPE(3,IR) - 1
              IB = I + NTYPE(3,IC) - 1
              CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFFSAVE,EMTTMP)
              EMTBLOCK(M,N) = EMTTMP
            ENDIF
           ENDDO
          ENDDO
       ENDDO
      END DO

!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation 

!   ENONSYM that does not include the symmetry-ordered orbital. 
      NVCOEF = 0
      ENONSYM = 0.D0
      COEFF = 0.0D0

      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
      DO IV = 1, NVCOEF
        VCOEFF = COEFF(IV)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC) 

          LABV=LABEL(1:6, IV)
          CALL ANALABV(IC, IR, IV)

          IF (NSYMCR .EQ. 0) THEN 
            CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
            ENONSYM = ENONSYM + EMTTMP

          ELSEIF (NSYMCR .EQ. 1 ) THEN 
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***1***'

          ELSEIF (NSYMCR .EQ. 2 ) THEN
           N = 0   ! Column index
           DO I = 1, NORBCOLL
            DO J = 1, NORBCOLU
             N = N + 1
             M = 0  ! Row index
             DO K = 1, NORBROWL
              DO L = 1, NORBROWU
               M = M + 1 
               !TRIANGLE MATRIX NEEDED
               IF (K .GT. I) CYCLE
               IF (K .EQ. I .AND. L.GT.J) CYCLE 

               IF (K .EQ. I) THEN 
                IF (L .EQ. J) THEN
                 !DIAGONAL MATRIXELEMENT IN EMTBLOCK
                 IF (LABEL(IPSYM(1),IV).EQ.NTYPE(4,IC)) THEN ! INNER symmetry-ordered ORB - CORE PAIR
                  LABV(IPSYM(1)) = NTYPE(3, IC) + I - 1
                  LABV(IPSYM(2)) = NTYPE(3, IR) + K - 1
                 ELSEIF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,IC)) THEN ! OUTER symmetry-ordered ORB - CORE PAIR
                  LABV(IPSYM(1)) = NTYPE(5, IC) + J - 1
                  LABV(IPSYM(2)) = NTYPE(5, IR) + L - 1
                 ELSE
      STOP 'UNEXPECTED CASES FOR DIAGONAL MATRIXELEMENTS WITHIN TYPE 3!'
                 ENDIF
                 CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP

                ELSE  ! K.EQ.I AND L.NE.J
                   IF (LABEL(IPSYM(1),IV).EQ.NTYPE(6,IC)) THEN
                    LABV(IPSYM(1)) = NTYPE(5, IC) + J - 1
                    LABV(IPSYM(2)) = NTYPE(5, IR) + L - 1
                    CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                   ENDIF                  
                ENDIF

               ! K.NE.I
               ! BUT THE OUTER ORB (L) OF IR IS SAME AS THE OUTER ORB (J) OF IC
               ELSEIF (L .EQ. J) THEN 
                IF (LABEL(IPSYM(1),IV).EQ.NTYPE(4,IC)) THEN
                  LABV(IPSYM(1)) = NTYPE(3, IC) + I - 1
                  LABV(IPSYM(2)) = NTYPE(3, IR) + K - 1
                  CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                ENDIF
               ENDIF
              
              ENDDO
             ENDDO
            ENDDO
           ENDDO

          ELSEIF (NSYMCR .EQ. 3 ) THEN
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***3***'

          ELSEIF (NSYMCR .EQ. 4 ) THEN 
           N=0 
           DO I = 1, NORBCOLL
            DO J = 1, NORBCOLU
             N = N + 1
             M = 0
             DO K = 1, NORBROWL
              DO L = 1, NORBROWU
               M = M + 1
               ! TRIANGLE MATRIX NEEDED
               IF (K .GT. I) CYCLE 
               IF (K .EQ. I .AND. L.GT.J) CYCLE  
               LABV(1) = NTYPE(3, IC) + I - 1
               LABV(2) = NTYPE(5, IC) + J - 1
               LABV(3) = NTYPE(3, IR) + K - 1
               LABV(4) = NTYPE(5, IR) + L - 1
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
! Add the common parts for the diagonal matrixelement 
       IF (DABS(ENONSYM) .GT. CUTOFF) THEN
        DO I = 1, NTYPE(2, IC)
           EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM
        ENDDO
      ENDIF

      IBUG1 = 0
!
!   Accumulate the contribution from the two-electron
!   transverse interaction operator
      IF (LTRANS .AND. (INC2.EQ.1)) THEN 
       ! Buiding the FICTIOUTS CSF as needed
         
       IF (NORBROWU.GE.2) THEN
         CALL FICTIOUS_CSF(2, IRFICT, IR, NTYPE(6,IR)-1, NTYPE(6,IR), 0, 0)
         IRORBIF(1,2) = IRFICT
       ENDIF
       IF (NORBROWL.GE.2) THEN 
         CALL FICTIOUS_CSF(2, IRFICT+1, IR, NTYPE(4,IR)-1, NTYPE(4,IR), 0, 0)
         IRORBIF(1,3) = IRFICT+1
       ENDIF

       ! Call twices "FICTIOUS_CSF(2,..)" to shift  (Core) 10s10d to
       ! (Core) 9s9d
       IF (NORBROWU.GE.2 .AND. NORBROWL.GE.2) THEN
         CALL FICTIOUS_CSF(2, IRFICT+2, IRFICT+1, &
                              NTYPE(6,IR)-1, NTYPE(6,IR), 0, 0)
         IRORBIF(1,1) = IRFICT+2
       ENDIF
 
       ENONSYM = 0.D0
       DO MTYPE = 1, 4
        NVCOEF = 0
        COEFF = 0
        IF (MTYPE.EQ.1) THEN
         ! < (Core) 10s10d | h_br | (Core) 9s 9d > 
         ! < K L | h_br | L J >: K != I, L != J
         IF (NORBCOLL.LT.2) CYCLE ! At least 2 orbitals in the same sym.  
         IF (NORBCOLU.LT.2) CYCLE ! At least 2 orbitals in the same sym. 
         
         !CALL RKCO_GG (IC, IR-NORBCOLU-1, BREID, 1, 2)
         CALL RKCO_GG (IC, IRORBIF(1,1), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.2) THEN
         ! < (Core) 10s10d | h_br | (Core ) 10s 9d >
         ! < K L | h_br | L J >: K = I, L != J
         IF (NORBCOLU.LT.2) CYCLE
         !CALL RKCO_GG (IC, IR-1, BREID, 1, 2)
         CALL RKCO_GG (IC, IRORBIF(1,2), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.3) THEN
         ! < (Core) 10s10d | h_br | (Core ) 9s10d >
         ! < K L | h_br | L J >: K != I,  L = J
         IF (NORBCOLL.LT.2) CYCLE
         !CALL RKCO_GG (IC, IR-NORBCOLU, BREID, 1, 2)
         CALL RKCO_GG (IC, IRORBIF(1,3), BREID, 1, 2)

        ELSEIF (MTYPE.EQ.4) THEN
         ! < (Core) 10s10d | h_br | (Core) 10s10d >
         ! < K L | h_br | L J >: K = I, L = J
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
             IF (MTYPE.NE.4) THEN
               WRITE(*,*)'IC,IR,IV,MTYPE=',IC,IR,IV,MTYPE
               STOP 'Unexpected MTYPE.NE.4 in matrixblock4.f ...'
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
             IRPTK = 0
             IRPTL = 0
             ICPTI = 0
             ICPTJ = 0
             IF (MTYPE.EQ.2) THEN
               ! K = I, L != J
               IF (LABEL(IPSym(1),IV).EQ.NTYPE(6,IR)-1) THEN
                 IRPTL = IPSym(1)
                 ICPTJ = IPSym(2)
               ELSE
                 IRPTL = IPSym(2)
                 ICPTJ = IPSym(1)
               ENDIF   
             ELSEIF (MTYPE.EQ.3) THEN
               ! K != I, L = J
               IF (LABEL(IPSym(1),IV).EQ.NTYPE(4,IR)-1) THEN
                 IRPTK = IPSym(1)
                 ICPTI = IPSym(2)
               ELSE
                 IRPTK = IPSym(2)
                 ICPTI = IPSym(1)
               ENDIF
             !ELSEIF (MTYPE.EQ.4) THEN
               ! Treated below 
             ENDIF
! Loop over symmetry-ordered orbitals 
             N = 0  ! Column index
             DO I = 1, NORBCOLL
              DO J = 1, NORBCOLU
                N = N + 1
                M = 0
                DO K = 1, NORBROWL
                 DO L = 1, NORBROWU
                  M = M + 1  ! Row index
                  ! Up-triangle matrix needed
                  IF (K .GT. I) CYCLE ! IR INNER ALWAYS <= IC INNER ORB
                  IF (K .EQ. I .AND. L.GT.J) CYCLE 

                  IF (MTYPE.EQ.2 .AND. K.EQ.I .AND. L.NE.J) THEN
                    ! < K L | h_br | L J >: K = I,  L != J
                    LABV(IRPTL) = NTYPE(5, IR) + L - 1
                    LABV(ICPTJ) = NTYPE(5, IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP

                  ELSEIF (MTYPE.EQ.3 .AND. K.NE.I .AND. L.EQ.J) THEN 
                    ! < K L | h_br | L J >: K != I,  L =  J
                    LABV(IRPTK) = NTYPE(3, IR) + K - 1
                    LABV(ICPTI) = NTYPE(3, IC) + I - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP

                  ELSEIF (MTYPE.EQ.4 .AND. K.EQ.I .AND. L.EQ.J) THEN
                    ! < K L | h_br | L J >: K = I, and L = J
                    IF (LABEL(IPSym(1),IV).EQ.NTYPE(4,IR)) THEN
                     LABV(IPSym(1)) = NTYPE(3, IR) + K - 1
                     LABV(IPSym(2)) = NTYPE(3, IC) + I - 1
                     CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                     EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                    ELSEIF (LABEL(IPSym(1),IV).EQ.NTYPE(6,IC)) THEN
                     LABV(IPSym(1)) = NTYPE(5, IR) + L - 1
                     LABV(IPSym(2)) = NTYPE(5, IC) + J - 1
                     CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                     EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                    ENDIF
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
               ! < K L | h_br | L J >: K != I,  L != J
               DO I = 1, 4
                 IF (LABEL(I,IV).EQ.NTYPE(4,IR)-1) IRPTK = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IR)-1) IRPTL = I
                 IF (LABEL(I,IV).EQ.NTYPE(4,IC)-0) ICPTI = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IC)-0) ICPTJ = I
               ENDDO
             ENDIF
 
             IF (MTYPE.EQ.2) THEN
               ! < K L | h_br | L J >: K = I,  L != J
               IRPTK = 0
               DO I = 1, 4
                 IF (LABEL(I,IV).EQ.NTYPE(6,IR)-1) IRPTL = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IC)) ICPTJ = I
                 IF (LABEL(I,IV).EQ.NTYPE(4,IR)) THEN
                   IF (IRPTK.EQ.0) THEN
                     IRPTK = I
                   ELSE
                     ICPTI = I
                   ENDIF 
                 ENDIF
               ENDDO 
             ENDIF

             IF (MTYPE.EQ.3) THEN
               ! < K L | h_br | L J >: K != I,  L = J
               IRPTL = 0
               DO I = 1, 4
                 IF (LABEL(I,IV).EQ.NTYPE(4,IR)-1) IRPTK = I
                 IF (LABEL(I,IV).EQ.NTYPE(4,IC)) ICPTI = I
                 IF (LABEL(I,IV).EQ.NTYPE(6,IR)) THEN
                   IF (IRPTL.EQ.0) THEN
                     IRPTL = I
                   ELSE
                     ICPTJ = I
                   ENDIF
                 ENDIF
               ENDDO
             ENDIF

             IF (MTYPE.EQ.4) THEN
              IRPTK = 0
              IRPTL = 0
              DO I = 1, 4
               IF (LABEL(I,IV).EQ.NTYPE(4,IC)) THEN
                 IF (IRPTK.EQ.0) THEN
                   IRPTK = I
                 ELSE
                   ICPTI = I
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
             DO I = 1, NORBCOLL
              DO J = 1, NORBCOLU
                N = N + 1
                M = 0
                DO K = 1, NORBROWL
                 DO L = 1, NORBROWU
                  M = M + 1  ! Row index
                  ! Up-triangle matrix needed
                  IF (K .GT. I) CYCLE ! IR INNER ALWAYS <= IC INNER ORB
                  IF (K .EQ. I .AND. L.GT.J) CYCLE

                  IF (MTYPE.EQ.1 .AND. K.NE.I .AND. L.NE.J)          &
                  THEN 
                    ! < K L | h_br | L J >: K != I, L != I, L != J
                    LABV(IRPTK) = NTYPE(3,IR) + K - 1
                    LABV(IRPTL) = NTYPE(5,IR) + L - 1
                    LABV(ICPTI) = NTYPE(3,IC) + I - 1
                    LABV(ICPTJ) = NTYPE(5,IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ELSEIF (MTYPE.EQ.2 .AND. K.EQ.I .AND. L.NE.J) THEN
                    ! < K L | h_br | L J >: K = I, L != J
                    LABV(IRPTK) = NTYPE(3,IR) + K - 1
                    LABV(IRPTL) = NTYPE(5,IR) + L - 1
                    LABV(ICPTI) = NTYPE(3,IC) + I - 1
                    LABV(ICPTJ) = NTYPE(5,IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ELSEIF (MTYPE.EQ.3 .AND. K.NE.I .AND. L.EQ.J) THEN
                    ! < K L | h_br | L J >: K != I, L = J
                    LABV(IRPTK) = NTYPE(3,IR) + K - 1
                    LABV(IRPTL) = NTYPE(5,IR) + L - 1
                    LABV(ICPTI) = NTYPE(3,IC) + I - 1
                    LABV(ICPTJ) = NTYPE(5,IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ELSEIF (MTYPE.EQ.4 .AND. K.EQ.I .AND. L.EQ.J) THEN
                    ! < K L | h_br | L J >: K = I, and L = J
                    LABV(IRPTK) = NTYPE(3,IR) + K - 1
                    LABV(IRPTL) = NTYPE(5,IR) + L - 1
                    LABV(ICPTI) = NTYPE(3,IC) + I - 1
                    LABV(ICPTJ) = NTYPE(5,IC) + J - 1
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
      END SUBROUTINE MATRIXBLOCK3
