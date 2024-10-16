!=======================================================================
!                                                                      *
! This subroutine a) increase NCRPT, NCRPST(NTPT)                      *
!                 b) IF (LFRST) build NDAOPT                           *
!                 b) IF (LPOT)  calculate DAOPT                        *
!                                                                      *
! Written by Chongyang Chen,      Fudan University,          Aug 2024  *
!                                                                      *
!=======================================================================
!!!!!!!!!
!  
!  Needing modifications for zero-first calculation, the present version
!  is implemneted for ZF approximation with BLOCK interaction.
!
!  For the original GRASP-2018 ZF approximation, i.e., only the diagonal 
!  elements are included, the interaction within the BLCOK generated by 
!  the same ICG - IRG CSFG pair are discarded, modifications are needed.
!
!  The above issue is resolved in calddrs.f90
!!!!!!!!!
      SUBROUTINE SETCOF_NDA_CSFG(IROFF, ICOFF, IA, IB, TCOEFF)

      USE vast_kind_param, ONLY: DOUBLE, LONG
      USE parameter_def,   ONLY: KEYORB

      Use orb_C,           ONLY: NW
      Use csfg_memory_man
      Use csfg_tv_C
      Use symmatrix_mod,   ONLY: IRTOT, ICTOT, LTRANSPOSE
      Use csfg_scf_C
      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)  :: IROFF, ICOFF, IA, IB
      REAL(DOUBLE), INTENT(IN)  :: TCOEFF
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: KEYSQ = KEYORB*KEYORB

      INTEGER :: IRR, ICC, IC1, IR1, LAB
      INTEGER :: I, J, ITHIS, IPOS
      LOGICAL :: LFOUND 
      REAL(DOUBLE) :: CONTR, UCFJ, SUM
!-----------------------------------------------
!
      IF (IA == IB .OR. IA.GT.IB) THEN
        WRITE (*,*)'IOFFC, IOFFR, IA, IB =', ICOFF, IROFF, IA, IB
        STOP 'Error, unexpected IA == IB .OR. &
              IA.GT.IB in savep_tlabcr.f90!'
      ENDIF
       
      IF (LSKIPPOT(IA) .AND. LSKIPPOT(IB)) RETURN
! Increament of NCRPT
      NCRPT = NCRPT + 1

! Transpose back the indexes of IC - IR. Here IC and IR are the CSF
! indexes.  
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

      IF (LCPOT) THEN
        CONTR = DDRS(IR1, IC1) * TCOEFF * 0.5D0 
        ! off-diagonal contributions have double the weight
        IF (IRR /= ICC) CONTR = CONTR + CONTR
        IF (ABS(CONTR) .LT. 1.0D-20) RETURN
      ENDIF

      DO I = 1, 2
        IF (I == 1) THEN
          J = IA; ITHIS = IB
        ELSE
          J = IB; ITHIS = IA
        ENDIF

        CALL LOCATE(ITHIS,NDCOFOPT(J),NDAOPT(:,J),LFOUND,IPOS)

        IF (LABTVFRST) THEN
          ! Construct the NDAOPT, Obtain NDCOFOPT
          IF (.NOT.LFOUND) THEN
            ! NDCOFOPT(:) and NDAOPT(:,:) are initialized in scfmpi.f90
            NDCOFOPT(J) = NDCOFOPT(J) + 1
            NDCOF = NDCOFOPT(J)
            IF (NDCOF > SIZE(NDAOPT,DIM=1)) THEN
             NDDIM = 2*NDDIM
             CALL RALLOC(NDAOPT, NDDIM, NW, 'NDAOPT', 'SETCOF_NDA_CSFG')
            ENDIF
            NDAOPT(IPOS+2:NDCOF,J) = NDAOPT(IPOS+1:NDCOF-1,J)
            NDAOPT(IPOS+1,J) = ITHIS
          ENDIF
        ELSE
          ! Here LCPOT = .TRUE. 
          !IF (.NOT.LCPOT) STOP "Error! .not. LCPOT in setcof_nda_csfg.."
          !IF (.NOT.LFOUND) THEN
          ! WRITE(*, *)"Error, Unexpected .NOT.LFOUND in SETCOF_NDA_CSFG"
          ! STOP "Error, Unexpected .NOT.LFOUND in SETCOF_NDA_CSFG"
          !ENDIF
          DAOPT(IPOS,J) = DAOPT(IPOS,J) +  CONTR/UCF(J)
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE SETCOF_NDA_CSFG