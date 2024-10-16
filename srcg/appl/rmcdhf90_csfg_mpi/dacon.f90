
!***********************************************************************
!                                                                      *
      SUBROUTINE DACON (J)
!                                                                      *
!   This  routine  includes  the  contribution from the off-diagonal   *
!   I(a,b) integrals in the 'exchange' term.                           *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT.                                        *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!
!   Last modification by Chongyang Chen, Fudan University, May, 2022   *
!   Output the contributions, neglecting when NP(J) > 15
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,   ONLY:  NNNP
      USE debug_C
      USE def_C
      USE grid_C
      USE npot_C
      USE orb_C
      USE pote_C
      USE csfg_scf_C
      USE tatb_C
      USE wave_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: J
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IORB, MFI, I, MFMAX 
      REAL(DOUBLE) :: TWOC, COEFF, FK, RPORII, PFI, QFI, ZBCI
!
!-----------------------------------------------
! CYC-DEBUG 
      REAL(DOUBLE) :: IabP(NNNP), IabQ(NNNP)

! CYC May, 2022
!      IF (NP(J) > 15) RETURN

! Debug printout: header
!
      IF (LDBPR(31) .OR. LDBPR(32)) THEN
        WRITE (99, 300) NP(J), NH(J)
        IabP = 0.D0
        IabQ = 0.D0
      ENDIF
! CYC May, 2022
!
!-----------------------------------------------
!
      TWOC = C + C
!
      MFMAX = 0
      DO K = 1, NDCOF
!
         IORB = NDA(K)
         CALL DPBDT (IORB)
         MFI = MF(IORB)
         COEFF = DA(K)
         FK = DBLE(NAK(IORB))
! CYC-DEBUG
         ! Debug printout: composition
         IF (MFI > MFMAX) MFMAX = MFI
         IF (LDBPR(31)) WRITE (99, 301) COEFF, NP(J), NH(J), NP(IORB), &
            NH(IORB), NP(IORB), NH(IORB)
!
         DO I = 2, MFI
            RPORII = 1.0D0/(H*RPOR(I))
            PFI = PF(I,IORB)
            QFI = QF(I,IORB)
            ZBCI = ZZ(I)/C
            XP(I) = XP(I) + COEFF*(TA(I)*RPORII+FK*PFI-(TWOC*R(I)+ZBCI)*QFI)
            XQ(I) = XQ(I) + COEFF*(TB(I)*RPORII-FK*QFI+ZBCI*PFI)
            IF (LDBPR(32)) THEN
              IabP(I) = IabP(I) + COEFF*(TA(I)*RPORII+FK*PFI-(TWOC*R(I)+ZBCI)*QFI)
              IabQ(I) = IabQ(I) + COEFF*(TB(I)*RPORII-FK*QFI+ZBCI*PFI)
            ENDIF
         END DO
!
      END DO

! CYC-DEBUG
      IF (LDBPR(32)) THEN
        WRITE (99, 302)
        DO I = 1, MFMAX
          WRITE (99, 303) R(I), IabP(I), IabQ(I), R(I), XP(I), XQ(I)
        ENDDO
        WRITE (99, 304) NP(J), NH(J)
      ENDIF
!
      RETURN
!
  300 FORMAT(/,/,' Exchange potential contributions from off-diagonal' &
          ' I(ab), COEF divided by C for ',1I2,1A2,' orbital :'/,/)
  301 FORMAT(1X,1P,D21.14,'* I (',1I2,1A2,',',1I2,1A2, ') ','* P (',1I2,1A2,')')
  302 FORMAT(/,/,31X,'(P)',19X,'(Q)',41X,'(P)',19X,'(Q)'/,&
         ' --------- r --------- ------ I  (r) -------',&
         ' ------ I  (r) -------', &
         ' --------- r --------- ------ X  (r) -------',&
         ' ------ X  (r) -------')
  303 FORMAT(1P,6(1X,1D21.14))
  304 FORMAT(/,/,2X, 'End build Direct, Exchange (V-Coef and T-Coef)' &
       ' Potential for', 2x, 1I2, 1A2/,/'=============')
!
      END SUBROUTINE DACON
