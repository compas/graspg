!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXBLOCK11(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,     &
                 NCORE,ELSTO)
!                                                                      *
!   This subroutine calls onescalar and computes one electron          *
!   matrix elements when IC and IR are of type 1 and 1                 *
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
      REAL(DOUBLE) :: EMTTMP,ATWINV,ENONSYM,TCOEFF,VCOEFF
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER :: IA,IB,IIA,IV,ISWAP
!-----------------------------------------------------------------------

      ATWINV = 1.D0/EMN
      IBUG1 = 0

!   Set all matrix elements to zero
      !EMTBLOCK = 0.0d0

!   Call onescalar 
      TSHELL = 0.D0
      IA = 0
      IB = 0
      CALL ONESCALAR(IC,IR,IA,IB,TSHELL)

      IF (IA .NE. 0) THEN
        LTRANSFER = .TRUE. 
        IF (IA .EQ. IB) THEN
         !Diagonal matrixelement, IC = IR.
         DO IIA = 1, NW
           TCOEFF = TSHELL(IIA)
           IF (ABS(TCOEFF) .GT. CUTOFF) THEN
             NCOEC = NCOEC + 1
             CALL onebody_DC_MS_VP(IIA,IIA,ATWINV,TCOEFF,EMTTMP)
             EMTBLOCK(1,1) = EMTBLOCK(1,1) + EMTTMP
           ENDIF
         ENDDO
        ELSE 
!   Ensure that the indices are in `canonical' order
         IF (IA .GT. IB) THEN
          ISWAP = IB
          IB = IA
          IA = ISWAP
         ENDIF
         TCOEFF = DBLE(TSHELL(1))
         IF (ABS(TCOEFF) .GT. CUTOFF) THEN
           NCOEC = NCOEC + 1
           CALL onebody_DC_MS_VP(IA,IB,ATWINV,TCOEFF,EMTTMP)
           EMTBLOCK(1,1) = EMTBLOCK(1,1) + EMTTMP
         ENDIF
        ENDIF
      ENDIF

      IBUG1 = 0
!
!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation; the latter
!   is computed first because the orbital indices may be
!   permuted by RKINTC
!
      ENONSYM = 0.0D0
      NVCOEF = 0
!
      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
      IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
      DO IV = 1, NVCOEF
        VCOEFF = COEFF(IV)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + 1
          LABV=LABEL(1:6, IV)
          CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
          EMTBLOCK(1,1) = EMTBLOCK(1,1) + EMTTMP
        ENDIF
      ENDDO

      IBUG1 = 0
      IF (LTRANS .AND. (INC2.EQ.1)) THEN
!
!   Accumulate the contribution from the two-electron
!   transverse interaction operator
!
        NVCOEF = 0
        CALL RKCO_GG (IC, IR, BREID, 1, 2)
        IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
        DO 10 IV = 1, NVCOEF
          VCOEFF = COEFF(IV)
          IF (DABS (VCOEFF) .GT. CUTOFF) THEN
            NMCBP = NMCBP + NTYPE(2,IC)
            LABV(1:6) = LABEL(1:6,IV)
            CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
            EMTBLOCK(1,1) = EMTBLOCK(1,1) + EMTTMP
          ENDIF
   10   CONTINUE
      ENDIF
      
!cyc  Transfer EMTBLOCK to EMTSYM 
      IF (LTRANSFER) Call transfer_cyc(IC, IR)     

      RETURN
      END SUBROUTINE MATRIXBLOCK11
