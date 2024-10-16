!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXBLOCK24(ICSAV,IRSAV,NCOEC,INCOR,NCTEC,INC2,     &
                 NMCBP,NCORE,ELSTO)
!                                                                      *
!   This subroutine calls onescalar and computes one electron          *
!   matrix elements when IC and IR are of type 2 and 4                 *
!                                                                      *
!   Written by Per Jönsson & Kai Wang                        May 2017 *
!                                                                      *
!***********************************************************************
!   Modified by CHONG-YANG CHEN                              JUNE 2020 *
!   Last modification by C. Y. Chen                          Dec  2023 *
!***********************************************************************
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP
      use symmatrix_mod
!-----------------------------------------------
!   C O M M O N  B L O C K S
!-----------------------------------------------
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
!      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-20
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: FLAG_IC1_IR1,FLAG_IC2_IR2,FLAGU
      REAL(DOUBLE) :: EMTTMP,ATWINV,ENONSYM,TCOEFF,VCOEFF
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER :: I,IA,IB,IV,ISWAP,IA0,IAIB0,IAIR0,IBIR0,IC,IC0,IC1,IC2
      INTEGER :: ICPTI,ICPTJ,ICW,IERROR,IR,IR0,IR1,IR2,IRPTK,IRW,J,K
      INTEGER :: M,MTYPE,N,NDIFF,NORBCOL,NORBROW
      INTEGER :: IRKEQI,IRKEQJ,IRKNEIJ
      INTEGER :: IRORB4(3)   ! The symmetry-ordered-orbital positions
!-----------------------------------------------------------------------
!
!  Set IC as TYPE 4, IR as TYPE 2.
!
      IF (LTRANSPOSE) THEN
         IC = IRSAV
         IR = ICSAV
      ELSE
         IC = ICSAV
         IR = IRSAV
      ENDIF

      ATWINV = 1.D0/EMN
      IBUG1 = 0

!   Set all matrix elements to zero
      !EMTBLOCK = 0.0d0
      ENONSYM = 0.0D0
! Number of symmetry-ordered orbs for the IC- and IR-th symmetry-ordered-CSFs
      NORBCOL = NTYPE(6,IC) - NTYPE(3,IC) + 1 
      NORBROW = NTYPE(4,IR) - NTYPE(3,IR) + 1 

!   There is possibly onebody contribution if IR has one of the two symmetry-ordered-Orbs
!   of IC symmetry-ordered-CSF.
      FLAGU = .FALSE.
      IF (NTYPE(3,IR).NE.NTYPE(3,IC)) GOTO 101
      FLAGU = .TRUE.

! Three cases are to be handled.
      ! symmetry-ordered-orbital of Row is equal to the inner one of Column
      ! IRW (4s) -- ICW (4s 5s)

      ! symmetry-ordered-orbital of Row is equal to the outer one of Column
      ! IR1 (5s) -- IC1 (4s 5s)

      ! Without same symmetry-ordered-orbital between the Row and Column
      ! IR2 (6s) -- IC2 (4s 5s), OR, IR2 (4s) -- IC2 (5s 6s)

! Call onescalar 

!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
      IF (LPRINT) &
        CALL PRINTLABELV(0, IC, IR, 24, 0, 0, 0, 0.d0)

! CYC++: Both MTYPE = 1 and 2 have to be retained:
! See file fort.1024
!BE  MTYPE=  0  ICW=     155  ICW*3+5=     470  IRW=     111  IRW*3+5=     338  MatrixblockXX= 24
!    IR0=     111  IR0*3+5=     338  IRE=     114  IRE*3+5=     347
!    IC0=     303  IC0*3+5=     914  ICE=     308  ICE*3+5=     929
!OB  MTYPE=  1  ICW=     155  ICW*3+5=     470  IRW=     492  IRW*3+5=    1481  MatrixblockXX= 24
!    IA=  13  IB=   1  TSHELL1=   1.225E+00
!OB  MTYPE=  2  ICW=     155  ICW*3+5=     470  IRW=     111  IRW*3+5=     338  MatrixblockXX= 24
!    IA=  12  IB=   1  TSHELL1=  -1.225E+00

      ICW = IC
      DO MTYPE = 1, 2
        IA = 0
        TSHELL = 0.D0
        IF (MTYPE.EQ.1) THEN
          !  < K (IR) | h1 | I J (IC) >: K = I 
          ! Construct IRFICT FICTIOUS CSF as needed.
          ! IC (1s2 9s 10s)   -- IR (1s 2s2 7s / 15s) are both replaced by:
          ! IC (1s2 9s 10s)   -- IRW (1s 2s2 9s)
          IRW = IRFICT
          IRORB4(1) = NTYPE(4,IC)
          CALL FICTIOUS_CSF(2, IRW, IR, NTYPE(4,IC), NTYPE(4,IR), 0, 0)
          ! For two-body DC and Breit
          IRKEQI = IRW   ! < IR (Core A) K -- IC (Core B) I J >: K.EQ.I 
        ELSE 
          !  < K (IR) | h1 | I J (IC) >: K = J 
          ! Construct IRFICT FICTIOUS CSF as needed.
          ! IC (1s2 9s 10s)   -- IR (1s 2s2 7s / 15s) are both replaced by:
          ! IC (1s2 9s 10s)   -- IRW (1s 2s2 10s)
          IRW = IRFICT+1
          IRORB4(2) = NTYPE(6,IC)
          CALL FICTIOUS_CSF(2, IRW, IR, NTYPE(6,IC), NTYPE(4,IR), 0, 0)
          ! For two-body DC and Breit
          IRKEQJ = IRW  ! < IR (Core A) K -- IC (Core B) I J >: K.EQ.J 
        ENDIF

        CALL ONESCALAR(ICW,IRW,IA,IB,TSHELL)
        IF (LPRINT) &
          CALL PRINTLABELV(1, ICW, IRW, 24, MTYPE, IA, IB, TSHELL(1))

        IF (IA.EQ.0.OR.DABS(TSHELL(1)).LT.CUTOFF) CYCLE
        LTRANSFER = .TRUE.
        TCOEFF = TSHELL(1)
        IF (IA.GT.IB) THEN
          ISWAP = IA
          IA = IB
          IB = ISWAP
        ENDIF
        IF (IA.GT.NORBGEN) THEN
          WRITE(*,*)'IC=',IC,IC*3+5,' IR=',IR,IR*3+5, ' BLOCK 2-4 ...'
          STOP 'Unexpected IA.GT.NORBGEN ...'
        ENDIF
        IF (IB.LE.NORBGEN) THEN
          WRITE(*,*)'IC=',IC,IC*3+5,' IR=',IR,IR*3+5, ' BLOCK 2-4 ...'
          STOP 'Unexpected IB.LE.NORBGEN ...'
        ENDIF

        NCOEC = NCOEC + NTYPE(2,IC)  
        N = 0  !Column index
        DO I = 1, NORBCOL-1
         DO J = I + 1, NORBCOL
          N = N + 1
          DO K = 1, NORBROW
            M = K
            IF (MTYPE.EQ.1 .AND. K.EQ.I) THEN
              IB = J + NTYPE(3,IC) - 1
              CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
              EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP 
            ELSEIF (MTYPE.EQ.2 .AND. K.EQ.J) THEN
            ! IC (1s2 4s 5s) -- IR (1s 2s2 5s)
              IB = I + NTYPE(3,IC) - 1
              CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
              EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP 
            ENDIF
          ENDDO
         ENDDO
        END DO
      ENDDO

101   CONTINUE
      IBUG1 = 0
!
!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation; 

! CYC++: MTYPE = 3 are removed. MTYPE = 1 / MTYPE = 2 loop also worked
! for (K /= I and K /= J, MTYPE=3) pairs. There are "Phase factor" 
! beteen the integrals invovling only ONE symmetry-ordered orbitals. The MTYPE =
! 2 loop are partially removed as below, according to LPHASET and
! TCOEF0, and LPHASE_NV.
! There is bug for the partial removement of MTYPE=2.
  
! See fort.1024, the below is one example. 

!TB  MTYPE=  1  ICW=     155  ICW*3+5=     470  IRW=     492  IRW*3+5=    1481  MatrixblockXX= 24
!    IVCOEF=          1          2          3          4
!    Label1=          1          2         12         12
!    Label2=         13         13         13         13
!    Label3=          1          2          1         12
!    Label4=          1          1         12          1
!    Label5=          0          0          0          0
!    Coeff =  1.225E+00  1.225E+00 -1.225E+00  1.225E+00
!TB  MTYPE=  2  ICW=     155  ICW*3+5=     470  IRW=     111  IRW*3+5=     338  MatrixblockXX= 24
!    IVCOEF=          1          2          3          4
!    Label1=          1          2         12         12
!    Label2=         12         12         13         13
!    Label3=          1          2          1         13
!    Label4=          1          1         13          1
!    Label5=          0          0          0          0
!    Coeff = -1.225E+00 -1.225E+00 -1.225E+00  1.225E+00
!TB  MTYPE=  3  ICW=     155  ICW*3+5=     470  IRW=     494  IRW*3+5=    1487  MatrixblockXX= 24
!    IVCOEF=          1          2
!    Label1=         12         12
!    Label2=         13         13
!    Label3=          1         10
!    Label4=         10          1
!    Label5=          0          0
!    Coeff = -1.225E+00  1.225E+00

      ICW = IC
      ! "CALL RKCO_GG" within MTYPE = 3 loop is removed.  
      ! The method to partially remove "CALL RKCO_GG" within MTYPE = 2
      ! worked for Be-like ions, but didn't work for B-like ions. 
      DO MTYPE = 1, 4
       NVCOEF = 0
       IF (MTYPE.EQ.1) THEN
         IF (.NOT.FLAGU) CYCLE 
         ! IC (1s2 4s 5s) -- IR (1s 2s2 4s)
         !  < K (IR) | h12 | I J (IC) >: K = I 
         IRW = IRKEQI
         CALL RKCO_GG (ICW, IRKEQI, CORD, INCOR, 1)

       ELSEIF (MTYPE.EQ.2) THEN
         IF (.NOT.FLAGU) CYCLE
         IF (NORBROW.LT.2) CYCLE 
         ! IC (1s2 4s 5s) -- IR (1s 2s2 5s)
         !  < K (IR) | h12 | I J (IC) >: K = J 
         IRW = IRKEQJ 
         CALL RKCO_GG (ICW, IRKEQJ, CORD, INCOR, 1)

       ELSEIF (MTYPE.EQ.3) THEN
         ! K /= I and K /= J
         IF (.NOT.(FLAGU)) CYCLE
         ! At least THREE obitals of same symmetry.
         IF (NORBROW.LT.3 .AND. NORBCOL.LT.3) CYCLE
         ! IC (1s2 6s 7s) -- IR (1s 2s2 8s)  
         ! OR:
         ! IC (1s2 6s 7s) -- IR (1s 2s2 6s)
         IF (NORBROW.LE.NORBCOL) THEN
            IRW = IRFICT + 2
            IRORB4(3) = NTYPE(3,IC)
            CALL FICTIOUS_CSF(2, IRW, IR, NTYPE(3,IC), NTYPE(4,IR), 0, 0)
         ELSE
           IRW = IR
           IRORB4(3) = NTYPE(4,IR)
         ENDIF
         ! For Breit
         IRKNEIJ = IRW  ! K /= I and K /= J
         
         ! CYCLE here, because Breit section needs the FICTIOUS CSF.
         ! Use the coefficients calculated within MTYPE = 1
         CYCLE
         CALL RKCO_GG (ICW, IRW, CORD, INCOR, 1)

       ELSEIF (MTYPE.EQ.4) THEN
        IF (FLAGU) CYCLE
        CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
       ENDIF

       IF (LPRINT .AND. FLAGU .AND. NVCOEF.GT.0) &
         CALL PRINTLABELV(2, ICW, IRW, 24, MTYPE, 0, 0, 0.d0)
       IF (NVCOEF .GT. 0) LTRANSFER = .TRUE. 
! CONTINUE THE MATRIXELEMNET CALCULATION
       DO IV = 1, NVCOEF
        VCOEFF = COEFF(IV)
        IF (DABS (VCOEFF) .LE. CUTOFF) CYCLE
        NCTEC = NCTEC + NTYPE(2,IC) 

! Determine the position of the symmetry-ordered orbitals
        LABV=LABEL(1:6, IV)
        CALL ANALABV(IC, IR, IV)

! THERE ARE AT MOST THREE symmetry-ordered-ORBS FOR TYPE 2 - 4 MATRIXELEMENTS.
        IF (NSYMCR .EQ. 0) THEN 
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
               'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***0***'

! Matrixelement between Type 2 - 4, ONLY when the IR-th CSF has ONE
! of the TWO symmetry-ordered-orbitals in the IC-th CSF.
        ELSEIF (NSYMCR .EQ. 1) THEN 
           N = 0  ! Column index
           DO I = 1, NORBCOL-1
            DO J = I + 1, NORBCOL
               N = N + 1
               DO K = 1, NORBROW
                M = K  ! Row index
                IF (MTYPE.EQ.1 .AND. K.EQ.I) THEN 
                  LABV(IPSym(1)) = NTYPE(3, IC) + J - 1
                  CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP

                ELSEIF (MTYPE.EQ.2 .AND. K.EQ.J) THEN
                  LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
                  CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                ENDIF
               ENDDO
            ENDDO
           ENDDO

        ELSEIF (NSYMCR .EQ. 2) THEN
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')             &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***2***'

        ELSEIF (NSYMCR .EQ. 3) THEN
           IF (MTYPE.NE.1 .AND. MTYPE.NE.4) CYCLE
           N = 0
           DO I = 1, NORBCOL-1
            LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
            DO J = I + 1, NORBCOL
               LABV(IPSym(2)) = NTYPE(3, IC) + J - 1
               N = N + 1
               DO K = 1, NORBROW
                M = K
                !IF (MTYPE.EQ.1) THEN
                !  IF (K.NE.I) CYCLE
                !ELSEIF (MTYPE.EQ.2) THEN
                !  IF (K.NE.J) CYCLE 
                !ELSEIF (MTYPE.EQ.3) THEN
                !  IF (K.EQ.I .OR. K.EQ.J) CYCLE
                !ENDIF
                LABV(IPSym(3)) = NTYPE(3, IR) + K - 1
                CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
               ENDDO
            ENDDO
           ENDDO
          
        ELSEIF (NSYMCR .EQ. 4) THEN
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***4***'
        ENDIF 

       ENDDO
      ENDDO

!   Accumulate the contribution from the two-electron
!   transverse interaction operator
      IBUG1 = 0
      ENONSYM = 0.0d0
      IF (LTRANS .AND. (INC2.EQ.1)) THEN
       ICW = IC
       DO MTYPE = 1, 4
        NVCOEF = 0
        IF (MTYPE.EQ.1) THEN
         ! < (Core) 9s 10s (IC) | h_br | (Core) 8s (IR) >
         ! < (Core) 5s  6s (IC) I J | h_br | (Core) 4s (IR) K >  OR 
         ! < (Core) 4s  5s (IC) I J | h_br | (Core) 6s (IR) K >  
         ! K != I, K != J
         IF (.NOT. FLAGU) CYCLE
         IF (NORBROW.LT.3 .AND. NORBCOL.LT.3) CYCLE
         CALL RKCO_GG (IC, IRKNEIJ, BREID, 1, 2)
         IRW = IRKNEIJ

        ELSEIF (MTYPE.EQ.2) THEN
         ! < (Core) 4s  5s (IC) I J | h_br | (Core) 4s (IR) K >
         ! K .EQ. I
         IF (.NOT. FLAGU) CYCLE
         CALL RKCO_GG (IC, IRKEQI, BREID, 1, 2)
         IRW = IRKEQI

        ELSEIF (MTYPE.EQ.3) THEN
         ! < (Core) 4s  5s (IC) I J | h_br | (Core) 5s (IR) K >
         ! K .EQ. J
         IF (.NOT. FLAGU) CYCLE
         CALL RKCO_GG (IC, IRKEQJ, BREID, 1, 2)
         IRW = IRKEQJ

        ELSEIF (MTYPE.EQ.4) THEN
         ! < (Core) 9s 10s | h_br | (Core) 10p >
         IF (FLAGU) CYCLE
         CALL RKCO_GG (IC, IR, BREID, 1, 2)
        ENDIF

        IF (LPRINT .AND. FLAGU .AND. NVCOEF.GT.0) &
          CALL PRINTLABELV(3, IC, IRW, 24, MTYPE, 0, 0, 0.d0)
        IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
        DO 10 IV = 1, NVCOEF
          VCOEFF = COEFF(IV)
          IF (DABS (VCOEFF) .GT. CUTOFF) THEN
            NMCBP = NMCBP + NTYPE(2,IC)
! Determine the position of the symmetry-ordered orbitals
            LABV = LABEL(1:6,IV)
            CALL ANALABV(IC, IR, IV)

            IF (NSYMCR .EQ. 0) THEN
! No contributions in matriblock24...
              WRITE(*,*)'IC,IR,IV,MTYPE=',IC,IR,IV,MTYPE
              STOP 'Unexpected NSYMCR .EQ. 0 in B24 ...'

! Matrixelement between Type 2 - 4, Loop for symmetry-ordered-orbitals
            ELSEIF (NSYMCR .EQ. 1) THEN
             IF (MTYPE.EQ.1 .OR. MTYPE.EQ.4 .OR.                     &
                 NTYPE(3,IC).NE.NTYPE(3,IR)) THEN
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')           &
                   'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
              STOP 'Warning!!! NSYMCR ERROR ***1***'
             ENDIF

             N = 0  ! Column index
             DO I = 1, NORBCOL-1
              DO J = I + 1, NORBCOL
                 N = N + 1
                 DO K = 1, NORBROW
                  M = K  ! Row index
                  IF (MTYPE.EQ.2 .AND. K.EQ.I) THEN
                   LABV(IPSym(1)) = NTYPE(3, IC) + J - 1
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ENDIF

                  IF (MTYPE.EQ.3 .AND. K.EQ.J) THEN
                   LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ENDIF
                 ENDDO
              ENDDO
             ENDDO

            ELSEIF (NSYMCR .EQ. 2) THEN
             WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')            &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
             STOP 'Warning!!! NSYMCR ERROR ***2***'

            ELSEIF (NSYMCR .EQ. 3) THEN
             IRPTK = 0
             ICPTI = 0
             ICPTJ = 0 
! Determine the positions of the symmetry-ordered orbitals 
             IF (MTYPE.EQ.1) THEN
               DO I = 1, 4
                IF (LABEL(I,IV).EQ.IRORB4(3)) IRPTK = I
                IF (LABEL(I,IV).EQ.NTYPE(4,IC)) ICPTI = I
                IF (LABEL(I,IV).EQ.NTYPE(6,IC)) ICPTJ = I
               ENDDO
             ELSEIF (MTYPE.EQ.2) THEN
               ICPTI = 0
               IF (IRORB4(1).NE.NTYPE(4,IC)) THEN
                 WRITE(*,*)"Unexpected IORB4(1).NE.NTYPE(4,IC) ..."
                 STOP "Unexpected IORB4(1).NE.NTYPE(4,IC) ..."
               ENDIF
               DO I = 1, 4
                IF (LABEL(I,IV).EQ.NTYPE(4,IC)) THEN
                  IF (ICPTI.EQ.0) THEN
                    ICPTI = I
                  ELSE
                    IRPTK = I
                  ENDIF
                ENDIF 
                IF (LABEL(I,IV).EQ.NTYPE(6,IC)) ICPTJ = I
               ENDDO
             ELSEIF (MTYPE.EQ.3) THEN
               IRPTK = 0
               IF (IRORB4(2).NE.NTYPE(6,IC)) THEN
                 WRITE(*,*)"Unexpected IORB4(2).NE.NTYPE(6,IC) ..."
                 STOP "Unexpected IORB4(2).NE.NTYPE(6,IC) ..."
               ENDIF
               DO I = 1, 4
                IF (LABEL(I,IV).EQ.NTYPE(6,IC)) THEN
                  IF (IRPTK.EQ.0) THEN
                    IRPTK = I
                  ELSE
                    ICPTJ = I
                  ENDIF
                ENDIF
                IF (LABEL(I,IV).EQ.NTYPE(4,IC)) ICPTI = I
               ENDDO
             ELSEIF (MTYPE.EQ.4) THEN
               DO I = 1, 4
                IF (LABEL(I,IV).EQ.NTYPE(4,IR)) IRPTK = I
                IF (LABEL(I,IV).EQ.NTYPE(4,IC)) ICPTI = I
                IF (LABEL(I,IV).EQ.NTYPE(6,IC)) ICPTJ = I
               ENDDO
             ENDIF

! Loop over the symmetry-ordered orbitals
             N = 0
             DO I = 1, NORBCOL-1
              LABV(ICPTI) = NTYPE(3, IC) + I - 1
              DO J = I + 1, NORBCOL
                 LABV(ICPTJ) = NTYPE(3, IC) + J - 1
                 N = N + 1 ! Column index
                 DO K = 1, NORBROW
                  M = K    ! Row index
                  IF (MTYPE.EQ.1) THEN
                    IF (K.EQ.I .OR. K.EQ.J) CYCLE
                  ELSEIF (MTYPE.EQ.2) THEN
                    IF (K.NE.I) CYCLE
                  ELSEIF (MTYPE.EQ.3) THEN
                    IF (K.NE.J) CYCLE
                  ENDIF
                  LABV(IRPTK) = NTYPE(3, IR) + K - 1
                  !LABV(ICPTI) = NTYPE(3, IC) + I - 1
                  !LABV(ICPTJ) = NTYPE(3, IC) + J - 1
                  CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                 ENDDO
              ENDDO
             ENDDO

            ELSEIF (NSYMCR .EQ. 4) THEN
             WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')            &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
             STOP 'Warning!!! NSYMCR ERROR ***4***'
            ENDIF

          ENDIF
   10   CONTINUE
       ENDDO
      ENDIF
      
!cyc  Transfer EMTBLOCK to EMTSYM 
      IF (LTRANSPOSE) EMTBLOCK=TRANSPOSE(EMTBLOCK)
      IF (LTRANSFER) Call transfer_cyc(ICSAV, IRSAV)

      RETURN
      END SUBROUTINE MATRIXBLOCK24
