!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXBLOCK2(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,      &
                 NCORE,ELSTO)
!                                                                      *
!   This subroutine calls onescalar and computes one electron          *
!   matrix elements when IC and IR are of type 2, IC == IR             *         
!                                                                      *
!   Written by Per Jönsson & Kai Wang                       May 2017  *
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
      USE symmatrix_mod
      USE buffer_C,   ONLY: NVCOEF, LABEL, COEFF
      USE debug_C,    ONLY: IBUG1
      USE decide_C
      USE def_C
      USE orb_C
      USE stat_C
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
      REAL(DOUBLE) :: EMTTMP,ATWINV,ENONSYM,TCOEFF,VCOEFF,TCOEFFSAVE
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER :: I,J,M,IA,IB,IV,ISWAP,MTYPE,N,NORBCOL,NORBROW
!-----------------------------------------------------------------------
! Check, it should be diagnoal here, i.e., IRTOT.eq.ICTOT
 
      IF (ICTOT .NE. IRTOT) THEN
        WRITE(*,*)'IC=',IC,' IR=',IR,' ICTOT=',ICTOT,' IRTOT=',IRTOT
        STOP 'Error, ICTOT .NE. IRTOT in matrixblock2 ...'
      ENDIF

      ATWINV = 1.D0/EMN
      IBUG1 = 0
      ! Number of symmetry-ordered orbs of IC-CSFG
      NORBCOL = NTYPE(4,IC) - NTYPE(3,IC) + 1
      ! Number of symmetry-ordered orbs of IR-CSFG 
      NORBROW = NTYPE(4,IR) - NTYPE(3,IR) + 1 
!   Set all matrix elements to zero
      !EMTBLOCK = 0.0D0
      TSHELL = 0.D0
      ENONSYM = 0.0D0
      IA = 0
      TCOEFFSAVE = 0.D0

      IF (LPRINT) &
        CALL PRINTLABELV(0, IC, IR, 2, 0, 0, 0, 0.d0)
!   Call onescalar        

      CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
      IF (LPRINT.AND.IA.NE.0) &
        CALL PRINTLABELV(1, IC, IR, 2, 0, IA, IB, TSHELL(1))
!   
!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
!
!  Diagonal matrielements, IC = IR, K = I       
!  To get the diagonal matrix elements, add contributions from each of
!  the labelling orbitals to ENONSYM 
      DO IA = 1,NW
        TCOEFF = DBLE(TSHELL(IA))
        IF (DABS (TCOEFF) .GT. CUTOFF) THEN
          IF (IA.GT.NORBGEN) THEN ! IA > NORBGEB --> IA ! Number of symmetry-ordered orbs
            TCOEFFSAVE = TCOEFF
            EXIT
          ELSE
            CALL onebody_DC_MS_VP(IA,IA,ATWINV,TCOEFF,EMTTMP)
            ENONSYM = ENONSYM + EMTTMP
          END IF
        END IF
      ENDDO
!  Diagonal matrielements, IC = IR, K = I, add the contributions arising
!  from the symmetry-ordered orbs.
      NCOEC = NCOEC + NORBCOL
      DO I = 1,NORBCOL
        IA = I + NTYPE(3,IC) - 1
        CALL onebody_DC_MS_VP(IA,IA,ATWINV,TCOEFFSAVE,EMTTMP)
        EMTBLOCK(I,I) = ENONSYM + EMTTMP
      END DO

!  Here comes the off-diagonal part within the diagonal symmetry-ordered block.
      IF (NORBROW.GT.1) THEN
       IA = 0
       TSHELL = 0.D0
       ! < (Core A) I (IC)   | h_br | (Core A) K >: I /= K
! CYC++ Here, the T-coefficient is just set as TCOEFFSAVE, the
! FICTIOUS CSF construction is kept for Breit calculation.  
       ! Construct FICTIOUS CSF, saving it as IQA(:, IRFICT),
       ! JQSA(:, :, IRFICT), and JCUPA(:, IRFICT).
       ! [ (Core) 7s (IR) ==> (Core) 6s (IRFICT) ]
       CALL FICTIOUS_CSF(2, IRFICT, IR, NTYPE(3,IR), NTYPE(4,IR), 0, 0)
       !Call PRINTFICTCSF(IC, IR, IRFICT)
!======================================================================
       !CALL ONESCALAR(IC,IRFICT,IA,IB,TSHELL)
       !IF (LPRINT.AND.IA.NE.0.AND.ABS(TSHELL(1)).GT.CUTOFF) &
       !  CALL PRINTLABELV(1, IC, IR, 2, 2, IA, IB, TSHELL(1))

       !TCOEFF = DBLE(TSHELL(1))
       !IF (IA /= 0 .AND. DABS(TCOEFF) .GT. CUTOFF) THEN
!======================================================================
       IF (DABS(TCOEFFSAVE).GT.CUTOFF) THEN
        NCOEC = NCOEC + NTYPE(2,IC)
        ! Triangle matrix needed  
        DO I = 1, NORBROW-1
          IA = I + NTYPE(3,IC) - 1
          DO J = I + 1, NORBCOL
            IB = J + NTYPE(3,IC) - 1
            ! IA is definitely smaller than IB
            !CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
            CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFFSAVE,EMTTMP)
            EMTBLOCK(I,J) = EMTTMP
          END DO
        END DO
       ENDIF
      ENDIF

!  Accumulate the contributions from the two-electron
!  Coulomb operator and the mass polarisation 

      NVCOEF = 0
      ENONSYM = 0.D0

!  There is no "Phase problem" for the two-body electronic
!  interaction, call RKCO_GG once.
      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
      IF (LPRINT .AND. NVCOEF.GT.0) &
        CALL PRINTLABELV(2, IC, IR, 2, 0, 0, 0, 0.d0)
       
      DO IV = 1, NVCOEF
        VCOEFF = COEFF(IV)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC) 

! Determine the position of the symmetry-ordered orbitals
          LABV=LABEL(1:6, IV)
          CALL ANALABV(IC, IR, IV)

          IF (NSYMCR .EQ. 0) THEN 
!   Diagonal case IC = IR. Start with computing the contribution
!   ENONSYM that does not involve the symmetry-ordered orbital. 
!            WRITE(*,'(a3, i7, 2x, i7, 2x, a5, 5i3)')
!     :           'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
!            STOP 'Warning!!! NSYMCR ERROR ***0***'
            CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
            ENONSYM = ENONSYM + EMTTMP

          ELSEIF (NSYMCR .EQ. 1) THEN 
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***1***'

! Matrixelement between Type 2 - 2
! Loop for symmetry-ordered-orbitals
          ELSEIF (NSYMCR .EQ. 2) THEN
            DO N = 1, NORBCOL ! Column index
              LABV(IPSYM(1)) = NTYPE(3, IC) + N - 1
              ! Triangle matrix needed
              DO M = 1, N   ! Row index
                LABV(IPSYM(2)) = NTYPE(3, IR) + M - 1
                CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
              ENDDO
            ENDDO

          ELSEIF (NSYMCR .EQ. 3) THEN
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***3***'

          ELSEIF (NSYMCR .EQ. 4) THEN
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***4***'
          ENDIF 
        ENDIF
      ENDDO
! Add the common parts for the diagonal matrixelement 
       IF (DABS(ENONSYM) .GT. CUTOFF) THEN
        DO I = 1, NORBCOL
           EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM
        ENDDO
      ENDIF

!  Accumulate the contribution from the two-electron
!  transverse interaction operator

! There is "Phase problem" for the Breit interaction? call RKCO_GG twice.

! CYC++: Please see file fort.1002, the question is how to re-use the
! coefficiencts calculated by MTYPE = 2, without MTYPE = 1 loop. 
! It seems to be very difficult, the TYPEs of Breit integrals differ
! from each other for MTYPE = 1 and 2 loop, as well as the number of
! V-coefficient, and also the values.
!BR  MTYPE=  1  ICW=     111  ICW*3+5=     338  IRW=     492  IRW*3+5=    1481  MatrixblockXX=  2
!    IVCOEF=          1          2          3          4          5          6          7          8          9
!    Label1=          2          2         10          2         10          1         10          1         10
!    Label2=          2         10          2         10          2         10          1         10          1
!    Label3=         13         13          2          2         13         13          1          1         13
!    Label4=         10          2         13         13          2          1         13         13          1
!    Label5=          1          1          1          1          1          1          1          1          1
!    Label6=          3          1          1          1          1          1          1          1          1
!    Coeff =  1.333E+00  6.667E-01  6.667E-01  6.667E-01  6.667E-01  6.667E-01  6.667E-01  6.667E-01  6.667E-01
!BR  MTYPE=  2  ICW=     111  ICW*3+5=     338  IRW=     492  IRW*3+5=    1481  MatrixblockXX=  2
!    IVCOEF=          1          2          3          4          5          6          7          8          9         10         11
!    Label1=          2         13          2          2          1          1          2          2          1         13         13
!    Label2=          2          2         13         13          1          2          1          1         13          1          1
!    Label3=         13         13         13          2          1          1          1          2          1          1         13
!    Label4=         13          2          2         13          1          2          2          1         13         13          1
!    Label5=          1          1          1          1          1          1          1          1          1          1          1
!    Label6=          4          5          5          5          4          5          5          5          5          5          5
!    Coeff =  2.667E+00  6.667E-01  1.333E+00  6.667E-01  2.667E+00  6.667E-01  1.333E+00  6.667E-01  6.667E-01  1.333E+00  6.667E-01

      IBUG1 = 0
      IF (LTRANS .AND. (INC2.EQ.1)) THEN 
       ENONSYM = 0.0d0
       DO MTYPE = 1, 2
        NVCOEF = 0
        COEFF = 0.0D0  
        IF (MTYPE.EQ.1) THEN
         ! < (Core A) 10s (IC) | h_br | (Core A) 9s >
         IF (NORBROW.LT.2) CYCLE
         ! < (Core A) I (IC)   | h_br | (Core A) K >: I /= K
         ! FICTIOUS_CSF(2, IRFICT, IR, NTYPE(3,IR), NTYPE(4,IR), 0, 0) called
         CALL RKCO_GG (IC, IRFICT, BREID, 1, 2)

         IF (LPRINT .AND. NVCOEF.GT.0) &
           CALL PRINTLABELV(3, IC, IRFICT, 2, MTYPE, 0, 0, 0.d0)

        ELSEIF (MTYPE.EQ.2) THEN
         ! < (Core A) 10s | h_br | (Core A) 10s >
         ! < (Core A) I (IC)   | h_br | (Core A) K >: I = K
         CALL RKCO_GG (IC, IR, BREID, 1, 2)

         IF (LPRINT .AND. NVCOEF.GT.0) &
           CALL PRINTLABELV(3, IC, IRFICT, 2, MTYPE, 0, 0, 0.d0)
        ENDIF

        DO 10 IV = 1, NVCOEF
          VCOEFF = COEFF(IV)
          IF (DABS (VCOEFF) .GT. CUTOFF) THEN
            NMCBP = NMCBP + NTYPE(2,IC)
! Determine the position of the symmetry-ordered orbitals
            LABV = LABEL(1:6,IV)
            CALL ANALABV(IC, IR, IV)

            IF (NSYMCR .EQ. 0) THEN
! Same for the diagonal matrixelement
              IF (MTYPE.NE.2) THEN
                WRITE(*,*)'IC,IR,IV,MTYPE=',IC,IR,IV,MTYPE
                STOP                                                 &
                   'Unexpected MTYPE.NE.2 in matrixblock2.f ...'
              ENDIF
              CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
              ENONSYM = ENONSYM + EMTTMP

! Matrixelement between Type 1 - 2, 5, Loop for symmetry-ordered-orbitals
            ELSEIF (NSYMCR .EQ. 1) THEN
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')           &
                   'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
              STOP 'Warning!!! NSYMCR ERROR ***1***'

            ELSEIF (NSYMCR .EQ. 2) THEN
             DO N = 1, NORBCOL ! Column index
              ! Triangle matrix needed
              DO M = 1, N ! Row index
               IF (MTYPE.EQ.1) THEN
                 IF (M.EQ.N) CYCLE
                 IF (LABEL(IPSym(1),IV).EQ.NTYPE(4,IC)) THEN
                   LABV(IPSym(1)) = NTYPE(3, IC) + N - 1
                   LABV(IPSym(2)) = NTYPE(3, IR) + M - 1
                 ELSE
                   LABV(IPSym(1)) = NTYPE(3, IR) + M - 1
                   LABV(IPSym(2)) = NTYPE(3, IC) + N - 1
                 ENDIF
                 CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
               ELSEIF (MTYPE.EQ.2) THEN
                 IF (M.NE.N) CYCLE
                 LABV(IPSym(1)) = NTYPE(3, IC) + N - 1
                 LABV(IPSym(2)) = NTYPE(3, IR) + M - 1
                 CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP
               ENDIF
              ENDDO
             ENDDO

            ELSEIF (NSYMCR .EQ. 3) THEN
             WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')            &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
             STOP 'Warning!!! NSYMCR ERROR ***3***'

            ELSEIF (NSYMCR .EQ. 4) THEN
             WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')            &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
             STOP 'Warning!!! NSYMCR ERROR ***4***'
            ENDIF

          ENDIF
   10   CONTINUE
       ENDDO
! Adding the common BREIT part
       IF (ENONSYM .NE. 0.0D0) THEN
         !IF (NTYPE(3,IC).NE.NTYPE(3,IR).OR.NORBCOL.NE.NORBROW) THEN
         !  WRITE(*,*)'IC=',IC, IC*3+5, ' IR=',IR, IR*3+5
         !  STOP 'UNEXPECTED ERROR IN BLOCK 2 - 2 ..FF..'
         !ENDIF
         DO I = 1, NORBCOL
           EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM
         ENDDO
       ENDIF
      ENDIF

!cyc  Transfer EMTBLOCK to EMTSYM 
      call transfer_cyc(IC, IR)

      RETURN
      END SUBROUTINE MATRIXBLOCK2
