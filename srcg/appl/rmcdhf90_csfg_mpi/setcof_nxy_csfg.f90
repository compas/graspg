!=======================================================================
!                                                                      *
! This subroutine a) increase NCRPV                                    *
!                 b) Construct the NXA, NYA arrays.                    * 
!                 c) Calculate XA and YA coefficients                  *
!                                                                      *
! Written by Chongyang Chen,      Fudan University,          Jan 2024  *
!                                                                      *
!=======================================================================

!!!!!!!!!
!CYC: Jan 2024
!  Needing modifications for zero-first calculation, the present version
!  is implemneted for ZF approximation with BLOCK interaction.
!
!  For the original GRASP-2018 ZF approximation, i.e., only the diagonal 
!  elements are included, the interaction within the sub-BLCOK generated 
!  by the same ICG - IRG CSFG (ICG = IRG) pair are discarded, 
!  modifications are needed.
!
!CYC: Aug 2024
!  The above issue is resolved, see calddrs.f90
!!!!!!!!!
      SUBROUTINE SETCOF_NXY_CSFG(IROFF, ICOFF, VCOEFF)

      USE vast_kind_param, ONLY: DOUBLE, LONG
      USE parameter_def,   ONLY: KEYORB
      USE mpi_C
      Use csfg_memory_man
      Use csfg_tv_C
      Use csfg_scf_C,      ONLY: LSKIPPOT, DDRS, XAWRK
      Use symmatrix_mod,   ONLY: IRTOT, ICTOT, LABV, LTRANSPOSE

      USE setcof_xy_csfg_I

      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)  :: IROFF, ICOFF
      REAL(DOUBLE), INTENT(IN)  :: VCOEFF
      
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: KEYSQ = KEYORB*KEYORB

      INTEGER        :: IRR, ICC, K, IC1, IR1
      INTEGER        :: LAB, LAC, LBD
      INTEGER(BYTE)  :: IA, IB, IC, ID, ISWAP
      LOGICAL        :: IFLAG 
      REAL(DOUBLE)   :: CONTR
!-----------------------------------------------
! The following should be checked carefully.
      IA = LABV(1); IB = LABV(2); IC = LABV(3); ID = LABV(4) 
      IF (LSKIPPOT(IA) .AND. LSKIPPOT(IB) .AND. &
          LSKIPPOT(IC) .AND. LSKIPPOT(ID)) RETURN
      K  = LABV(5) 

! Transpose back the indexes of IC - IR
      IF (LTRANSPOSE) THEN
        ICC = IROFF + ICTOT
        IRR = ICOFF + IRTOT
        IC1 = IROFF
        IR1 = ICOFF
      ELSE
        ICC = ICOFF + ICTOT
        IRR = IROFF + IRTOT
        IC1 = ICOFF
        IR1 = IROFF
      ENDIF

! CYC: The following constrcution of LABEL and de-construction in
! SETCOF_XY_CSFG could be improved. For extremely large-scale
! calculation, the total number of SLATER integrals is very large, maybe
! exceeding 2^31, all the LABELs packed and unpacked repeatedly.

! Construct LABELs
      IF (IA > IC) THEN
        ISWAP = IC
        IC = IA
        IA = ISWAP
      ENDIF
      IF (IB > ID) THEN
        ISWAP = ID
        ID = IB
        IB = ISWAP
      ENDIF
      IF (IA == IB .AND. IC > ID) THEN
        ISWAP = IC
        IC = ID
        ID = ISWAP
      ENDIF

      LAC = IA * KEY + IC
      LBD = IB * KEY + ID
      IF (LAC .LT. LBD) THEN
        LAB = LAC * KEYSQ + LBD
      ELSE
        LAB = LBD * KEYSQ + LAC
      ENDIF

      IF (LABTVFRST) THEN
! Only construct the NXA, NYA arrays.
        IFLAG = .NOT. LABTVFRST
! Increament of NCRPV, Total number of SLATER integrals involving
! potentials, it differs from those construct H-matrix.
        NCRPV = NCRPV + 1
        CALL SETCOF_XY_CSFG(IFLAG, LAB, K, CONTR) 
      ELSE
! LCPOT = .TRUE., accumulate the xa and ya coefficients.
        IFLAG = LCPOT
        !CONTR = DSUBRS(EOL,IR,IC,JBLOCK)*VCOEFF
        CONTR = DDRS(IR1,IC1)*VCOEFF
        !*** off-diagonal contributions have double the weight
        IF (IRR /= ICC) CONTR = CONTR + CONTR
        IF (ABS(CONTR).GT.1.0D-15) &
          CALL SETCOF_XY_CSFG(IFLAG, LAB, K, CONTR)
      ENDIF 
      
      !CALL SETCOF_XY_CSFG(IFLAG, IA, IB, IC, ID, K, CONTR)

      RETURN
      END SUBROUTINE SETCOF_NXY_CSFG
