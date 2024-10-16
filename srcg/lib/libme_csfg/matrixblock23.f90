!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXBLOCK23(ICSAV,IRSAV,NCOEC,INCOR,NCTEC,INC2,     &
                 NMCBP,NCORE,ELSTO)
!                                                                      *
!   This subroutine calls onescalar and computes one electron          *
!   matrix elements when IC and IR are of type 2 and 3                 *
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
      LOGICAL :: FLAGL,FLAGU
      REAL(DOUBLE) :: EMTTMP,ATWINV,ENONSYM,TCOEFF,VCOEFF
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER :: I,IA,IB,IC,IC0,IV,ISWAP,ICPTI,ICPTJ,ICW,IR,IR0,IRPTK
      INTEGER :: IRW,J,K,M,MTYPE,N,NORBCOLL,NORBCOLU,NORBROW
      INTEGER :: LABVI,LABVJ,LABVK,LABVL,IR1ORB,IR1
!-----------------------------------------------------------------------
!  Set IC as TYPE 3, IR as TYPE 2.
!
      IF (LTRANSPOSE) THEN
         IC = IRSAV
         IR = ICSAV
      ! The following two lines are meanlingless within Smallest list
      ! verion.
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
      IC0 = IC
      IR0 = IR

      ATWINV = 1.D0/EMN
      IBUG1 = 0

!   Set all matrix elements to zero
      !EMTBLOCK = 0.0d0
      ENONSYM = 0.0D0
! Number of symmetry-ordered orbs for the IC- and IR-th symmetry-ordered-CSFs
      NORBCOLL = NTYPE(4,IC) - NTYPE(3,IC) + 1 
      NORBCOLU = NTYPE(6,IC) - NTYPE(5,IC) + 1 
      NORBROW  = NTYPE(4,IR) - NTYPE(3,IR) + 1 
      FLAGL = .FALSE.
      FLAGU = .FALSE.

!   There is possibly onebody contribution if IR has one of the two symmetry-ordered-Orbs
!   of IC symmetry-ordered-CSF.
      ! CSFs-Pair without same symmetry-ordered-orbital
      IF (NTYPE(3,IR).NE.NTYPE(3,IC).AND.                  &
          NTYPE(3,IR).NE.NTYPE(5,IC)) GOTO 101
      
      ! CSFs-Pair with same symmetry-ordered-orbital
      IF (NTYPE(3,IR).EQ.NTYPE(3,IC)) THEN
        FLAGL = .TRUE.
        ! In case of FLEXIBLE CSFs, possible NORBROW.NE.NORBCOLL
      
      ELSEIF (NTYPE(3,IR).EQ.NTYPE(5,IC)) THEN
        FLAGU = .TRUE.
        ! In case of FLEXIBLE CSFs, possible NORBROW.NE.NORBCOLL
      ENDIF

      ! Buiding the FICTIOUTS CSF as needed.
      IF (FLAGL) THEN
        IF (NORBROW .NE. NORBCOLL) THEN
          ! (Core A) 9s / 5s (IR) -- (Core B) 7s 7p (IC)  ==> 
          ! (Core A) 7s (IRFICT)  -- (Core B) 7s 7p (IC)
          CALL FICTIOUS_CSF(2, IRFICT, IR, NTYPE(4,IC), NTYPE(4,IR),0,0)
          IR0 = IRFICT  ! K  = I
          IR1 = IR      ! K /= I
          IR1ORB = NTYPE(4,IR)
        ELSEIF (NORBROW.GE.2) THEN  ! NORBROW.EQ.NORBCOLL
          ! (Core A) 7s      (IR) -- (Core B) 7s 7p (IC)  ==> 
          ! (Core A) 4s (IRFICT)  -- (Core B) 7s 7p (IC)
          CALL FICTIOUS_CSF(2, IRFICT, IR, NTYPE(3,IR), NTYPE(4,IR),0,0)
          IR0 = IR      ! K  = I
          IR1 = IRFICT  ! K /= I
          IR1ORB = NTYPE(3,IR)
        ENDIF
      ENDIF
      
      IF (FLAGU) THEN
        IF (NORBROW .NE. NORBCOLU) THEN
          ! (Core A) 9p / 5p (IR) -- (Core B) 7s 7p (IC)  ==> 
          ! (Core A) 7p (IRFICT)  -- (Core B) 7s 7p (IC)
          CALL FICTIOUS_CSF(2, IRFICT, IR, NTYPE(6,IC), NTYPE(4,IR),0,0)
          IR0 = IRFICT  ! K  = J
          IR1 = IR      ! K /= J
          IR1ORB = NTYPE(4,IR)
        ELSEIF (NORBROW.GE.2) THEN ! NORBROW.EQ.NORBCOLU
          ! (Core A) 7p      (IR) -- (Core B) 7s 7p (IC)  ==> 
          ! (Core A) 4p (IRFICT)  -- (Core B) 7s 7p (IC)
          CALL FICTIOUS_CSF(2, IRFICT, IR, NTYPE(3,IR), NTYPE(4,IR),0,0)
          IR0 = IR      ! K  = J
          IR1 = IRFICT  ! K /= J
          IR1ORB = NTYPE(3,IR)
        ENDIF
      ENDIF

      IF (LPRINT .AND. (FLAGL.OR.FLAGU)) &
        CALL PRINTLABELV(0, IC, IR, 23, 0, 0, 0, 0.d0)

      IA = 0
      IB = 0 
      TSHELL = 0.D0
      ! Here, FLAGL or FLAGU is .true.
      ! IR 4p- -- IC 4p- 4p+ 
      ! IR 4p+ -- IC 4p- 4p+
      ! Changing IC/IR, ensure that IR and IC have the same symmetry-ordered-orbitals
      ICW = IC0
      IRW = IR0
      CALL ONESCALAR(ICW,IRW,IA,IB,TSHELL)
!   
!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
!
!   For for type 2 and 3 there are only onebody contributions when
!   the symmetry-ordered-orbital of IR is same as one of the two symmetry-ordered ones of IC
!
      IF (LPRINT.AND.IA.NE.0.AND.DABS(TSHELL(1)).GT.CUTOFF) &
        CALL PRINTLABELV(1, IC, IR, 23, 0, IA, IB, TSHELL(1))

      IF (IA .NE. 0) THEN
         LTRANSFER = .TRUE.
!        WRITE(*,*)'IC, IR, IA, IB=', IC, IR, IA, IB
!   Ensure that the indices are in `canonical' order
        IF (IA .GT. IB) THEN
          ISWAP = IB
          IB = IA
          IA = ISWAP
          IF (IA.GT.NORBGEN.OR.IB.LE.NORBGEN) then
            WRITE(*,*)'IC,IR,IA,IB,TSHELL=',IC,IR,IA,IB,TSHELL(1)
            STOP 'Error, IA.gt.IB and IB.le.NORBGEN in matrixblock23*1*'
          ENDIF
        ENDIF
        TCOEFF = DBLE(TSHELL(1))
        IF (DABS(TCOEFF) .GT. CUTOFF) THEN
          NCOEC = NCOEC + NTYPE(2,IC)  
          N = 0  ! Column index
          DO I = 1, NORBCOLL
           LABVI = I + NTYPE(3,IC) - 1
           DO J = 1, NORBCOLU
             LABVJ = J + NTYPE(5,IC) - 1
             N = N + 1
            DO K = 1,NORBROW
             M = K  ! Row index
             IF (FLAGL .AND. K.EQ.I) THEN
              !IB = J + NTYPE(5,IC) - 1
              IB = LABVJ
             ELSEIF (FLAGU .AND. K.EQ.J) THEN
              !IB = I + NTYPE(3,IC) - 1
              IB = LABVI
             ELSE
              CYCLE
             ENDIF
             CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
             EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP 
            ENDDO
           ENDDO
          END DO
        ENDIF
      ENDIF

101   CONTINUE
      IBUG1 = 0
!
!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation; the latter
!   is computed first because the orbital indices may be
!   permuted by RKINTC
     
      NVCOEF = 0
      ENONSYM = 0.0D0
! Possible NSYMCR.EQ.1 contributions, using IC0 and IR0, insteadly IC/IR
      ICW = IC0
      IRW = IR0
      CALL RKCO_GG (ICW, IRW, CORD, INCOR, 1)
      IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
      IF (LPRINT .AND. (FLAGL.OR.FLAGU) .AND. NVCOEF.GT.0) &
        CALL PRINTLABELV(2, ICW, IRW, 23, 0, 0, 0, 0.d0)

      DO IV = 1, NVCOEF
        VCOEFF = COEFF(IV)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC) 

! Determine the position of the symmetry-ordered orbitals
          LABV=LABEL(1:6, IV)
          CALL ANALABV(IC, IR, IV)

          IF (NSYMCR .EQ. 0) THEN 
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')             &
                 'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***0***'

! Matrixelement between Type 2 - 3, ONLY when the IR-th CSF has ONE
! of the TWO symmetry-ordered-orbitals in the IC-th CSF.
          ELSEIF (NSYMCR .EQ. 1) THEN 
           N = 0  ! Column index
           DO I = 1, NORBCOLL
             LABVI = NTYPE(3, IC) + I - 1
            DO J = 1, NORBCOLU
               LABVJ = NTYPE(5, IC) + J - 1
               N = N + 1
               DO K = 1, NORBROW
                M = K  ! Row index
                IF (FLAGL .AND. K.EQ.I) THEN
                 !LABV(IPSym(1)) = NTYPE(5, IC) + J - 1
                 LABV(IPSym(1)) = LABVJ
                ELSEIF (FLAGU .AND. K.EQ.J) THEN
                 !LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
                 LABV(IPSym(1)) = LABVI
                ELSE
                 CYCLE  
                ENDIF
                CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
               ENDDO
            ENDDO
           ENDDO

          ELSEIF (NSYMCR .EQ. 2) THEN
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')             &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***2***'

          ELSEIF (NSYMCR .EQ. 3) THEN
           N = 0
           DO I = 1, NORBCOLL
            LABVI = NTYPE(3, IC) + I - 1   
            DO J = 1, NORBCOLU
               LABVJ = NTYPE(5, IC) + J - 1
               N = N + 1
               DO K = 1, NORBROW
                M = K
                !LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
                !LABV(IPSym(2)) = NTYPE(5, IC) + J - 1
                LABV(IPSym(1)) = LABVI
                LABV(IPSym(2)) = LABVJ
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
        ENDIF
      ENDDO

! There is no common part in Matrixelement 2 - 3
!
!   Accumulate the contribution from the two-electron
!   transverse interaction operator
!

! CYC++: MTYPE = 1 should be kept. See the differences between the
! MTYPE = 1 and MTYPE = 2 coefficients within fort.1023
      IBUG1 = 0
      ENONSYM = 0.0d0
      IF (LTRANS .AND. (INC2.EQ.1)) THEN
       DO MTYPE = 1, 3
        NVCOEF = 0
        COEFF = 0.D0
        IF (MTYPE.EQ.1) THEN
         ! FLAGL.OR.FLAGU, CSFs-Pairs without same symmetry-ordered-orbital
         ! < (Core) 10s10d | h_br | (Core) 9s >
         ! < (Core) 10s10d | h_br | (Core) 9d >
         IF (.NOT.(FLAGL.OR.FLAGU)) CYCLE
         IF (FLAGL .AND. NORBROW.LT.2 .AND. NORBCOLL.LT.2) CYCLE
         IF (FLAGU .AND. NORBROW.LT.2 .AND. NORBCOLU.LT.2) CYCLE

         !IF (FLAGL) THEN
         !! FLAGL: IR 5p- -- IC 4p- 4p+, OR IR 4p- -- IC 5p- 4p+
         !  IF (NORBROW.GE.2) THEN
         !    ICW = IC0
         !    IRW = IR0 + 1
         !  ELSEIF (NORBCOLL.GE.2) THEN
         !    ICW = IC0 + NORBCOLU
         !    IRW = IR0
         !  ENDIF
         !ELSEIF (FLAGU) THEN
         !! FLAGU: IR 5p+ -- IC 4p- 4p+, OR IR 4p+ -- IC 4p- 5p+
         !  IF (NORBROW.GE.2) THEN
         !    ICW = IC0
         !    IRW = IR0 + 1
         !  ELSEIF (NORBCOLU.GE.2) THEN
         !    ICW = IC0 + 1
         !    IRW = IR0
         !  ENDIF
         !ENDIF

         ! Smallest list version
         ICW = IC0
         IRW = IR1
         CALL RKCO_GG (ICW, IRW, BREID, 1, 2)

        ELSEIF (MTYPE.EQ.2) THEN
         ! FLAGL.OR.FLAGU, CSFs-Pairs with same symmetry-ordered-orbital
         IF (.NOT.(FLAGL.OR.FLAGU)) CYCLE
         ! < (Core) 10s | h_br | (Core ) 10s 10d >  IR -- IC
         ! < (Core) 10d | h_br | (Core ) 10s 10d >
         ! < (Core)  4s | h_br | (Core )  4s  4d >  IR0 -- IC0
         ! < (Core)  4d | h_br | (Core )  4s  4d >
         ICW = IC0
         IRW = IR0
         CALL RKCO_GG (ICW, IRW, BREID, 1, 2)

        ELSEIF (MTYPE.EQ.3) THEN
         ! .NOT.(FLAGL.OR.FLAGU)
         IF (FLAGL.OR.FLAGU) CYCLE
         ! < (Core) 10p  (IR)  | h_br | (Core) 10s 10d (IC) >  
         ! < (Core) 4p   (IR0) | h_br | (Core) 4s 4d (IC0) >
         ICW = IC0
         IRW = IR0
         CALL RKCO_GG (ICW, IRW, BREID, 1, 2)
        ENDIF

      IF (LPRINT .AND. (FLAGL.OR.FLAGU) .AND. NVCOEF.GT.0) &
        CALL PRINTLABELV(3, ICW, IRW, 23, MTYPE, 0, 0, 0.d0)
        IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
        DO 10 IV = 1, NVCOEF
          VCOEFF = COEFF(IV)
          IF (DABS (VCOEFF) .GT. CUTOFF) THEN
            NMCBP = NMCBP + NTYPE(2,IC)
! Determine the position of the symmetry-ordered orbitals
            LABV = LABEL(1:6,IV)
            CALL ANALABV(IC, IR, IV)

            IF (NSYMCR .EQ. 0) THEN
! No contributions in matriblock23...
              WRITE(*,*)'IC,IR,IV,MTYPE=',IC,IR,IV,MTYPE
              STOP 'Unexpected NSYMCR .EQ. 0 in B24 ...'

! Matrixelement between Type 2 - 3, Loop for symmetry-ordered-orbitals
            ELSEIF (NSYMCR .EQ. 1) THEN
             IF (MTYPE.EQ.1 .OR. MTYPE.EQ.3) THEN
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')           &
                   'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
              STOP 'Unexpected MTYPE.EQ.1 .OR. MTYPE.EQ.3 in B23.f ...'
             ENDIF

             N = 0  ! Column index
             DO I = 1, NORBCOLL
              DO J = 1, NORBCOLU
                 N = N + 1
                 DO K = 1, NORBROW
                  M = K  ! Row index

                  ! < (Core) 10s | h_br | (Core ) 10s 10d > IR -- IC
                  ! < (Core)  4s | h_br | (Core )  4s  4d > IR0-- IC0
                  IF (FLAGL .AND. K.EQ.I) THEN
                   LABV(IPSym(1)) = NTYPE(5, IC) + J - 1
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ENDIF

                  ! < (Core) 10d | h_br | (Core ) 10s 10d > IR -- IC
                  ! < (Core)  4d | h_br | (Core )  4s  4d > IR0-- IC0
                  IF (FLAGU .AND. K.EQ.J) THEN
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
                IF (LABEL(I,IV).EQ.IR1ORB) IRPTK = I
                IF (LABEL(I,IV).EQ.NTYPE(4,ICW)) ICPTI = I
                IF (LABEL(I,IV).EQ.NTYPE(6,ICW)) ICPTJ = I
               ENDDO
             ELSEIF (MTYPE.EQ.2) THEN
               IRPTK = 0
               IF (FLAGL) THEN
                DO I = 1, 4
                   IF (LABEL(I,IV).EQ.NTYPE(6,ICW)) ICPTJ = I
                   IF (LABEL(I,IV).EQ.NTYPE(4,ICW)) THEN
                     IF (IRPTK.EQ.0) THEN
                       IRPTK = I
                     ELSE
                       ICPTI = I
                     ENDIF
                   ENDIF
                ENDDO  
               ELSEIF (FLAGU) THEN
                DO I = 1, 4
                   IF (LABEL(I,IV).EQ.NTYPE(4,ICW)) ICPTI = I
                   IF (LABEL(I,IV).EQ.NTYPE(6,ICW)) THEN
                     IF (IRPTK.EQ.0) THEN
                       IRPTK = I
                     ELSE
                       ICPTJ = I
                     ENDIF
                   ENDIF
                ENDDO  
               ENDIF
             ELSEIF (MTYPE.EQ.3) THEN
               DO I = 1, 4
                IF (LABEL(I,IV).EQ.NTYPE(4,IRW)) IRPTK = I
                IF (LABEL(I,IV).EQ.NTYPE(4,ICW)) ICPTI = I
                IF (LABEL(I,IV).EQ.NTYPE(6,ICW)) ICPTJ = I
               ENDDO
             ENDIF

! Loop over the symmetry-ordered orbitals
             N = 0  ! Column index
             DO I = 1, NORBCOLL
              LABVI = NTYPE(3, IC) + I - 1
              DO J = 1, NORBCOLU
                 LABVJ = NTYPE(5, IC) + J - 1
                 N = N + 1
                 DO K = 1, NORBROW
                  M = K  ! Row index
                  IF (MTYPE.EQ.1) THEN
                   IF (FLAGL .AND. K.NE.I) THEN
                    LABV(IRPTK) = NTYPE(3, IR) + K - 1
                    LABV(ICPTI) = LABVI
                    LABV(ICPTJ) = LABVJ
                    !LABV(ICPTI) = NTYPE(3, IC) + I - 1
                    !LABV(ICPTJ) = NTYPE(5, IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                   ENDIF
                   IF (FLAGU .AND. K.NE.J) THEN
                    LABV(IRPTK) = NTYPE(3, IR) + K - 1
                    LABV(ICPTI) = LABVI
                    LABV(ICPTJ) = LABVJ
                    !LABV(ICPTI) = NTYPE(3, IC) + I - 1
                    !LABV(ICPTJ) = NTYPE(5, IC) + J - 1
                    CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                   ENDIF

                  ELSEIF (MTYPE.EQ.2) THEN
                   IF (FLAGL .AND. K.NE.I) CYCLE
                   IF (FLAGU .AND. K.NE.J) CYCLE
                   LABV(IRPTK) = NTYPE(3, IR) + K - 1
                   LABV(ICPTI) = LABVI
                   LABV(ICPTJ) = LABVJ
                   !LABV(ICPTI) = NTYPE(3, IC) + I - 1
                   !LABV(ICPTJ) = NTYPE(5, IC) + J - 1
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP

                  ELSEIF (MTYPE.EQ.3) THEN
                   LABV(IRPTK) = NTYPE(3, IR) + K - 1
                   LABV(ICPTI) = LABVI
                   LABV(ICPTJ) = LABVJ
                   !LABV(ICPTI) = NTYPE(3, IC) + I - 1
                   !LABV(ICPTJ) = NTYPE(5, IC) + J - 1
                   CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                   EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
                  ENDIF
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
      END SUBROUTINE MATRIXBLOCK23
