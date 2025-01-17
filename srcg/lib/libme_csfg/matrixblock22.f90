!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXBLOCK22(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,     &
                 NCORE,ELSTO)
!                                                                      *
!   This subroutine calls onescalar and computes one electron          *
!   matrix elements when IC and IR are of type 2 and 2                 *
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
      USE parameter_def,   ONLY: NNNP, NNNW
      USE symmatrix_mod
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
      LOGICAL :: FLAGU
      REAL(DOUBLE) :: EMTTMP,ATWINV,ENONSYM,TCOEFF,VCOEFF
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER :: I,IA,IB,IV,ISWAP,IC0,ICW,IL0,IR0,IRW,JU0,M,MTYPE,N
      INTEGER :: NORBCOL,NORBROW
!-----------------------------------------------------------------------
      ! The first COLUMN generated by symmetry-ordered-CSF ICC
      ! The first ROW generated by symmetry-ordered-CSF IRR
!      IC0 = ICTOT + 1
!      IR0 = IRTOT + 1

      ATWINV = 1.D0/EMN
      IBUG1 = 0

!   Set all matrix elements to zero
      !EMTBLOCK = 0.0d0

      IF (NTYPE(3,IC) .EQ. NTYPE(3,IR)) THEN
        FLAGU = .TRUE.
      ELSE
        FLAGU = .FALSE.
      ENDIF

      NORBCOL = NTYPE(4,IC) - NTYPE(3,IC) + 1
      NORBROW = NTYPE(4,IR) - NTYPE(3,IR) + 1

      ICW = IC
      IRW = IR 
      IF (.NOT.FLAGU) GOTO 101

!   Call onescalar        
! There is non-zero contribution for the CSFs with the same symmetry-ordered-Orb 
! but different core-configurations, also only for the diagonal-pair CSFs
! generated by such two type 2 symmetry-ordered-CSFs, such are (Core A) 11s --
! (Core B) 11s

! FlEXIBLE CSFs, such as: 
! IR:  1s2 2s 10s, VV up to nmax1 = 10
! IC1: 1s 2s2 8s,  CV up to nmax2 =  8
! IC2: 2s2 2s 6s,  CC up to nmax3 =  6 
! Needing modification, match the lower n value of symmetry-ordered
! orbital of IR and IC.
! For the above THREE CSFs, nsym in IR should be respectively replaced 
! by 8 in < IR |h_1, h_2, h_br| IC1 >, and 
! by 6 in < IR |h_1, h_2, h_br| IC2 >.
! In the above two cases, there is no common onebody contribution.

      IF (LPRINT) &
        CALL PRINTLABELV(0, IC, IR, 22, 0, 0, 0, 0.d0)

!   
!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
      TSHELL = 0.D0
      IA = 0
      IB = 0
      IF (NORBROW.NE.NORBCOL) THEN
! CYC++:! Construct IRFICT as the symmetry-ordered orbitals are not equal for IC/IR.
        ! < IR: (Core B) K | h1 | IC: (CoreA) I>: set K = I to accout
        ! for the possible one-body contributions. 

        ! IC: (Core A) 7s  --  IR: (Core B) 11s
        ! IR: (Core B) 11s ==> IRFICT: (Core B) 7s
        ! OR:
        ! IC: (Core A) 7s  --  IR: (Core B) 5s
        ! IR: (Core B) 5s  ==> IR: (Core B) 7s
        CALL FICTIOUS_CSF(2, IRFICT, IR, NTYPE(4,IC), NTYPE(4,IR), 0, 0)
        IRW = IRFICT
      ELSE
        IRW = IR 
      ENDIF
      CALL ONESCALAR(ICW,IRW,IA,IB,TSHELL)
      IF (LPRINT.AND.IA.NE.0.AND.DABS(TSHELL(1)).GT.CUTOFF) &
        CALL PRINTLABELV(1, IC, IR, 22, 1, IA, IB, TSHELL(1))

!
!   IA = 0, the CSFs differ by two or more orbitals and there is no 
!   onebody interaction;
!   IA != 0 and IA != IB, the CSFs differ by one orbital. 
!   For such CSF-pair, the one-body contribution arises from
!   the core-orbitals, such as < 1s2 2s 11s | h1 | 1s2 3s 11s>, the
!   symmetry-ordered-orbitals are orthogonal normalized.
!
      IF (IA .NE. 0) THEN
        LTRANSFER = .TRUE.
!   Ensure that the indices are in `canonical' order
        IF (IA .GT. IB) THEN
          ISWAP = IB
          IB = IA
          IA = ISWAP
          IF (TSHELL(1).NE.0.0d0.AND.                                &
               (IB.GT.NORBGEN.OR.NTYPE(3,IC).NE.NTYPE(3,IR))) THEN
              WRITE(*,*)'IC,IR,IA,IB,TSHELL=',IC,IR,IA,IB,TSHELL(1)
              STOP 'Error, IA.gt.IB and IB.gt.NORBGEN in matrixblock12.'
          ENDIF
        ENDIF

        TCOEFF = DBLE(TSHELL(1))
        IF (DABS(TCOEFF) .GT. CUTOFF) THEN
          NCOEC = NCOEC + NTYPE(2,IC) 
          !No symmetry-ordered-Orb contributions, no needing to change IB
          CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
!   The contribution is common for all semi-diagonal matrixelement.
          DO I = 1,MIN(NORBCOL,NORBROW)
            EMTBLOCK(I,I) = EMTBLOCK(I,I) + EMTTMP 
          END DO
        ENDIF
      ENDIF

101   CONTINUE
!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation 

      NVCOEF = 0
      ENONSYM = 0.0D0
      ! IF (FLAGU): IRW / IRFICT has already been set as needed. 
      CALL RKCO_GG (ICW, IRW, CORD, INCOR, 1)
      IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
      IF (LPRINT .AND. FLAGU .AND. NVCOEF.GT.0) &
        CALL PRINTLABELV(2, IC, IR, 22, 0, 0, 0, 0.d0)

      DO IV = 1, NVCOEF
        VCOEFF = COEFF(IV)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC) 

! Determine the position of the symmetry-ordered orbitals
          LABV=LABEL(1:6, IV)
          CALL ANALABV(IC, IR, IV)

          IF (NSYMCR .EQ. 0) THEN 
! There are core twobody contributions for the two TYPE 2 CSFs with same
! symmetry-ordered-orbital, such as [Core A] 11s -- [Core B] 11s
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
              DO M = 1, NORBROW ! Row index
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
! Add the common parts for the semi-diagonal matrixelement 
      IF (DABS(ENONSYM) .GT. CUTOFF) THEN
        IF (NTYPE(3,IC).NE.NTYPE(3,IR)) THEN
          WRITE(*,*)'IC=',IC, IC*3+5, ' IR=',IR, IR*3+5
          STOP 'UNEXPECTED ERROR IN BLOCK 2 - 2 ..FF..'
        ENDIF
        DO I = 1, MIN(NORBCOL,NORBROW)
          EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM
        ENDDO
      ENDIF

!   Accumulate the contribution from the two-electron
!   transverse interaction operator

! CYC++: Please see file fort.1022, it seems that MTYPE= 3 should be
! kept to deal with the "Phase problem". 1s,2s,...,3d-,3d+ (orbitals
! 1-9) are non-symmetry-ordered (labelling) orbitals. 

!BR  MTYPE=  2  ICW=     114  ICW*3+5=     347  IRW=     111  IRW*3+5=     338  MatrixblockXX= 22
!    IVCOEF=          1          2          3          4          5          6          7          8          9
!    Label1=          5          5         13          5         13          1          2          1          2
!    Label2=          2         13          5         13          5          2          1          2          1
!    Label3=         13         13          2          2         13          5          1          1          5
!    Label4=         13          2         13         13          2          1          5          5          1
!    Label5=          1          1          1          1          1          1          1          1          1
!    Label6=          3          1          1          1          1          1          1          1          1
!    Coeff =  1.333E+00  6.667E-01  6.667E-01  6.667E-01  6.667E-01  6.667E-01  6.667E-01  6.667E-01  6.667E-01
!BR  MTYPE=  3  ICW=     114  ICW*3+5=     347  IRW=     492  IRW*3+5=    1481  MatrixblockXX= 22
!    IVCOEF=          1          2          3          4          5          6          7          8
!    Label1=          5          2          5          2          5         12          5         12
!    Label2=          2          5          2          5         12          5         12          5
!    Label3=         13         12         12         13         13          2          2         13
!    Label4=         12         13         13         12          2         13         13          2
!    Label5=          1          1          1          1          1          1          1          1
!    Label6=          1          1          1          1          1          1          1          1
!    Coeff =  6.667E-01  6.667E-01  6.667E-01  6.667E-01  6.667E-01  6.667E-01  6.667E-01  6.667E-01

      IBUG1 = 0
      IF (LTRANS .AND. (INC2.EQ.1)) THEN
       ENONSYM = 0.0d0
       DO MTYPE = 1, 3
        NVCOEF = 0
        IF (MTYPE.EQ.1) THEN
         ! Different symmetry-ordered orbitals
         ! < (Core A) 10p-| h_br | (Core B) 10p >
         IF (FLAGU) CYCLE
         CALL RKCO_GG (ICW, IRW, BREID, 1, 2)
         IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
        ELSEIF (MTYPE.EQ.2) THEN
         ! Semi-diagonal matrixelement
         IF (.NOT.FLAGU) CYCLE
         ! IR (CORE A) 4s -- IC (CORE B) 4s
         ! < IR: (Core B) K | h_br | IC: (CoreA) I>: K = I
         ! ICW / IRW have the same symmetry-ordered orbital within ONESCALAR/DC
         CALL RKCO_GG (ICW, IRW, BREID, 1, 2)


        ELSEIF (MTYPE.EQ.3) THEN
         IF (.NOT.FLAGU) CYCLE
         IF (NORBROW.LT.2 .AND. NORBCOL.LT.2) CYCLE
         ! Non-semi-diagonal matrixelement
         ! < (Core A) 10s | h_br | (Core B) 9s >
         ! IR 4s -- IC 5s, OR, IR 5s -- IC 4s
         ! < (Core A) 5s (IR) | h_br | (Core B) 4s (IC) >
         ! 1s,2s,3s are labelling orbitals
         IF (NORBROW.NE.NORBCOL) THEN
           ICW = IC
           IRW = IR
         ELSE
         ! Here, there are at least TWO sym-orbitals same as NTYPE(4,IR)/NTYPE(4,IC)
           ICW = IC
           IRW = IRFICT
           CALL FICTIOUS_CSF(2, IRFICT, IR, NTYPE(4,IR)-1, NTYPE(4,IR), 0, 0)
         ENDIF 
         CALL RKCO_GG (ICW, IRW, BREID, 1, 2)
        ENDIF

        IF (LPRINT .AND. FLAGU .AND. NVCOEF.GT.0) &
          CALL PRINTLABELV(3, ICW, IRW, 22, MTYPE, 0, 0, 0.d0)

        DO 10 IV = 1, NVCOEF
          VCOEFF = COEFF(IV)
          IF (DABS (VCOEFF) .GT. CUTOFF) THEN
            NMCBP = NMCBP + NTYPE(2,IC)
! Determine the position of the symmetry-ordered orbitals
            LABV = LABEL(1:6,IV)
            CALL ANALABV(IC, IR, IV)

            IF (NSYMCR .EQ. 0) THEN
! Same for the semi-diagonal matrixelement
              IF (MTYPE.NE.2.OR.NTYPE(3,IC).NE.NTYPE(3,IR)) THEN
                WRITE(*,*)'IC,IR,IV,MTYPE=',IC,IR,IV,MTYPE
                STOP                                                 &
                   'Unexpected MTYPE.NE.2.OR.NTYPE(3,IC).NE.NTYPE(3,IR)'
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
              IL0 = NTYPE(3, IC) + N - 1
              DO M = 1, NORBROW ! Row index
               JU0 = NTYPE(3, IR) + M - 1
               IF (MTYPE.EQ.1) THEN
                 IF (LABEL(IPSym(1),IV).EQ.NTYPE(4,ICW)) THEN
                   LABV(IPSym(1)) = IL0
                   LABV(IPSym(2)) = JU0
                 ELSE
                   LABV(IPSym(1)) = JU0
                   LABV(IPSym(2)) = IL0
                 ENDIF
                 CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP

               ELSEIF (MTYPE.EQ.2) THEN
                 IF (M.NE.N) CYCLE
                 LABV(IPSym(1)) = IL0
                 LABV(IPSym(2)) = JU0
                 CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + EMTTMP

               ELSEIF (MTYPE.EQ.3) THEN
                 IF (M.EQ.N) CYCLE
                 IF (LABEL(IPSym(1),IV).EQ.NTYPE(4,IC)) THEN
                   LABV(IPSym(1)) = IL0
                   LABV(IPSym(2)) = JU0
                 ELSE
                   LABV(IPSym(1)) = JU0
                   LABV(IPSym(2)) = IL0
                 ENDIF
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
       IF (DABS(ENONSYM) .GT. CUTOFF) THEN
        IF (NTYPE(3,IC).NE.NTYPE(3,IR)) THEN
          WRITE(*,*)'IC=',IC, IC*3+5, ' IR=',IR, IR*3+5
          STOP 'UNEXPECTED ERROR IN BLOCK 2 - 2 ..BREIT..'
        ENDIF
        DO I = 1, MIN(NORBCOL,NORBROW)
           EMTBLOCK(I,I) = EMTBLOCK(I,I) + ENONSYM
        ENDDO
       ENDIF
      ENDIF
      
!cyc  Transfer EMTBLOCK to EMTSYM
      IF (LTRANSFER) Call transfer_cyc(IC, IR)

      RETURN
      END SUBROUTINE MATRIXBLOCK22
