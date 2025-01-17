!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXBLOCK25(ICSAV,IRSAV,NCOEC,INCOR,NCTEC,INC2,     &
                 NMCBP,NCORE,ELSTO)
!                                                                      *
!   This subroutine calls onescalar and computes one electron          *
!   matrix elements when IC and IR are of type 2 and 5                 *
!                                                                      *
!   Written by Per Jönsson & Kai Wang                         May 2017 *
!                                                                      *
!***********************************************************************
!   Modified by CHONG-YANG CHEN                              JUNE 2020 *
!   Last modification by C. Y. Chen                          Jan  2022 *
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
!      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-20
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: FLAGSAME
      REAL(DOUBLE) :: EMTTMP,ATWINV,ENONSYM,TCOEFF,VCOEFF
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER :: I,IA,IB,IC,IC0,IV,ISWAP,ICPTI,ICPTJ,ICW,IR,IR0,IRPTK
      INTEGER :: IRW,J,K,M,MTYPE,N,NORBCOL,NORBROW,IR1,IORBKNEI
!-----------------------------------------------------------------------
!
!  Set IC as TYPE 5, IR as TYPE 2.
!
      IF (LTRANSPOSE) THEN
         IC = IRSAV
         IR = ICSAV
      !! The first COLUMN generated by symmetry-ordered-CSF ICC
      !! The first ROW generated by symmetry-ordered-CSF IRR
      !   IC0 = IRTOT + 1
      !   IR0 = ICTOT + 1
      ELSE
         IC = ICSAV
         IR = IRSAV
      ! CYC++: Smallest list used, the followings are meaningless.
      !   IC0 = ICTOT + 1
      !   IR0 = IRTOT + 1
      ENDIF
      IC0 = IC
      IR0 = IR

      ATWINV = 1.D0/EMN
      IBUG1 = 0

!   Set all matrix elements to zero
      !EMTBLOCK = 0.0d0
      ENONSYM = 0.0D0
! Number of symmetry-ordered orbs for the IC- and IR-th symmetry-ordered-CSFs
      NORBCOL = NTYPE(4,IC) - NTYPE(3,IC) + 1 
      NORBROW = NTYPE(4,IR) - NTYPE(3,IR) + 1 
      FLAGSAME = .FALSE.
      
!   There is possibly onebody contribution if IR has one of the two symmetry-ordered-Orbs
!   of IC symmetry-ordered-CSF.
      IF (NTYPE(3,IR).NE.NTYPE(3,IC)) GOTO 101
      FLAGSAME = .TRUE.
      
      ! Buiding the FICTIOUTS CSF as needed.
      IF (NORBCOL.NE.NORBROW) THEN
        ! (Core A) 9s / 5s (IR) -- (Core B) 7s2 (IC)  ==> 
        ! (Core A) 7s (IRFICT)  -- (Core B) 7s2 (IC)
        CALL FICTIOUS_CSF(2, IRFICT, IR, NTYPE(4,IC), NTYPE(4,IR), 0, 0)
        IR0 = IRFICT  ! K  = I
        IR1 = IR      ! K /= I
        IORBKNEI = NTYPE(4,IR)
      ELSEIF (NORBROW.GE.2) THEN
        ! (Core A) 7s (IR)      -- (Core B) 7s2 (IC)  ==> 
        ! (Core A) 4s (IRFICT)  -- (Core B) 7s2 (IC)
        CALL FICTIOUS_CSF(2, IRFICT, IR, NTYPE(3,IR), NTYPE(4,IR), 0, 0)
        IR0 = IR      ! K  = I
        IR1 = IRFICT  ! K /= I
        IORBKNEI = NTYPE(3,IR)
      ENDIF

      IF (LPRINT) &
        CALL PRINTLABELV(0, IC, IR, 25, 0, 0, 0, 0.d0)
       
!   Call onescalar        
      TSHELL = 0.D0
      ICW = IC0
      IRW = IR0
      CALL ONESCALAR(ICW,IRW,IA,IB,TSHELL)

      IF (LPRINT) &
        CALL PRINTLABELV(1, ICW, IRW, 25, 0, IA, IB, TSHELL(1))

!   
!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
!
      IF (IA .NE. 0) THEN
!   Ensure that the indices are in `canonical' order
        LTRANSFER = .TRUE.
        IF (IA .GT. IB) THEN
          ISWAP = IB
          IB = IA
          IA = ISWAP
          IF (IB.LE.NORBGEN) then
            WRITE(*,*)'IC,IR,IA,IB,TSHELL=',IC,IR,IA,IB,TSHELL(1)
            STOP 'Error, IA.gt.IB and IB.le.NORBGEN in matrixblock12...'
          ENDIF
        ENDIF

        TCOEFF = DBLE(TSHELL(1))
        IF (DABS(TCOEFF) .GT. CUTOFF) THEN
          NCOEC = NCOEC + NTYPE(2,IC) 
          DO I = 1, MIN(NORBCOL,NORBROW)
            IB = NTYPE(3,IC) + I - 1
            CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
            EMTBLOCK(I,I) = EMTTMP 
          ENDDO
        ENDIF
      ENDIF

101   CONTINUE
      IBUG1 = 0
!
!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation
!
      NVCOEF = 0
      ICW = IC0
      IRW = IR0
      CALL RKCO_GG (ICW, IRW, CORD, INCOR, 1)
      IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
      IF (LPRINT .AND. FLAGSAME .AND. NVCOEF.GT.0) &
        CALL PRINTLABELV(2, ICW, IRW, 25, 0, 0, 0, 0.d0)

      ENONSYM = 0.0D0
      DO IV = 1, NVCOEF
        VCOEFF = COEFF(IV)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC) 

! Determine the position of the symmetry-ordered orbitals
          LABV=LABEL(1:6, IV)
          CALL ANALABV(IC, IR, IV)

! THERE ARE AT MOST THREE symmetry-ordered-orbs FOR TYPE 2 - 5 MATRIXELEMENTS.
          IF (NSYMCR .EQ. 0) THEN 
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')             &
                 'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***0***'

! Matrixelement between Type 2 - 5, ONLY when the IR-th CSF has ONE
! of the TWO symmetry-ordered-orbitals in the IC-th CSF.
          ELSEIF (NSYMCR .EQ. 1) THEN 
           IF (NSYMCOL.NE.1) THEN
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')             &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***1***'
           ENDIF

           DO I = 1, MIN(NORBCOL,NORBROW)
             LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
             CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
             EMTBLOCK(I,I) = EMTBLOCK(I,I) + EMTTMP
           ENDDO

          ELSEIF (NSYMCR .EQ. 2) THEN
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')             &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***2***'

          ELSEIF (NSYMCR .EQ. 3) THEN
           !FULL MATRIX NEEDED, NOT TRIANGLE.
           DO I = 1, NORBCOL
             LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
             LABV(IPSym(2)) = NTYPE(3, IC) + I - 1
             DO K = 1, NORBROW
               LABV(IPSym(3)) = NTYPE(3, IR) + K - 1
               CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
               EMTBLOCK(K,I) = EMTBLOCK(K,I) + EMTTMP
             ENDDO
           ENDDO
          
          ELSEIF (NSYMCR .EQ. 4) THEN
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***4***'
          ENDIF 
        ENDIF
      ENDDO

!
      IBUG1 = 0
!   Accumulate the contribution from the two-electron
!   transverse interaction operator
      IF (LTRANS .AND. (INC2.EQ.1)) THEN
       ENONSYM = 0.0d0
       DO MTYPE = 1, 2
        NVCOEF = 0
        COEFF = 0.D0
        IF (MTYPE.EQ.1) THEN
         IF (.NOT.FLAGSAME) CYCLE
         IF (NORBROW.LT.2.AND.NORBCOL.LT.2) CYCLE
         ! K /= I
         ICW = IC0
         IRW = IR1
         CALL RKCO_GG (ICW, IRW, BREID, 1, 2)

        ELSEIF (MTYPE.EQ.2) THEN
         ! < (Core)  4p | h_br | (Core )  4s2 >
         ! < (Core)  4s | h_br | (Core )  4s2 >
         ! K = I
         ICW = IC0
         IRW = IR0
         CALL RKCO_GG (ICW, IRW, BREID, 1, 2)
        ENDIF

       IF (LPRINT .AND. FLAGSAME .AND. NVCOEF.GT.0) &
         CALL PRINTLABELV(3, ICW, IRW, 25, MTYPE, 0, 0, 0.d0)
        IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
        DO 10 IV = 1, NVCOEF
          VCOEFF = COEFF(IV)
          IF (DABS (VCOEFF) .GT. CUTOFF) THEN
            NMCBP = NMCBP + NTYPE(2,IC)
! Determine the position of the symmetry-ordered orbitals
            LABV = LABEL(1:6,IV)
            CALL ANALABV(IC, IR, IV)

            IF (NSYMCR .EQ. 0) THEN
! No contributions in matriblock25...
              WRITE(*,*)'IC,IR,IV,MTYPE=',IC,IR,IV,MTYPE
              STOP 'Unexpected NSYMCR .EQ. 0 in B25 ...'

! Matrixelement between Type 2 - 5, Loop for symmetry-ordered-orbitals
            ELSEIF (NSYMCR .EQ. 1) THEN
             IF (MTYPE.EQ.1.OR.(.NOT.FLAGSAME)) THEN
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')           &
                   'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
              STOP 'Unexpected MTYPE.EQ.1.OR.(.NOT.FLAGSAME) in B25.f '
             ENDIF

             N = 0  ! Column index
             DO I = 1, NORBCOL
                 N = I
                 DO K = 1, NORBROW
                  M = K  ! Row index
                  IF (K.EQ.I) THEN
                   LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ENDIF
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
                IF (LABEL(I,IV).EQ.IORBKNEI) IRPTK = I
                IF (LABEL(I,IV).EQ.NTYPE(4,ICW)) THEN
                  IF (ICPTI.EQ.0) THEN
                    ICPTI = I
                  ELSE
                    ICPTJ = I
                  ENDIF
                ENDIF
               ENDDO
             ELSEIF (MTYPE.EQ.2) THEN 
              IF (FLAGSAME) THEN 
               ! K = I = J
               IRPTK = IPSym(1)
               ICPTI = IPSym(2)
               ICPTJ = IPSym(3)
              ELSE
               DO I = 1, 4
                 IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) IRPTK = I
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

! Loop over the symmetry-ordered orbitals
             N = 0
             DO I = 1, NORBCOL
              J = I
                 N = I
                 DO K = 1, NORBROW
                  M = K  ! Row index
                  IF (MTYPE.EQ.1) THEN
                   IF (K.EQ.I) CYCLE
                   LABV(IRPTK) = NTYPE(3, IR) + K - 1
                   LABV(ICPTI) = NTYPE(3, IC) + I - 1
                   LABV(ICPTJ) = NTYPE(3, IC) + J - 1
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ELSEIF (MTYPE.EQ.2) THEN
                   IF (FLAGSAME.AND.K.NE.I) CYCLE
                   LABV(IRPTK) = NTYPE(3, IR) + K - 1
                   LABV(ICPTI) = NTYPE(3, IC) + I - 1
                   LABV(ICPTJ) = NTYPE(3, IC) + J - 1
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ENDIF
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
      END SUBROUTINE MATRIXBLOCK25
