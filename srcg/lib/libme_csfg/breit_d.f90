!***********************************************************************
!                                                                      *
      SUBROUTINE BREIT_D(IC, IR, IV, VCOEFF, NCORE, ELSTO, EMTTMP)
!                                                                      *
!   This subroutine perform the BREIT-INTERACTON calculations          *
!                                                                      *
!   Modified by CHONG-YANG CHEN           JUNE 2020                    *
!                                                                      *
!   Last modification by Chongyang Chen   Nov 2023                     *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE orb_C,           ONLY: NP, NKL
      use csfg_decide_C,   ONLY: NDISCARDBR, LDISCARDBR, NDISBR, LDISBR
      use symmatrix_mod,   ONLY: LABV
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
      EXTERNAL BREID, CORD
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER      :: IC, IR, IV, NCORE
      REAL(DOUBLE) :: VCOEFF, ELSTO, EMTTMP
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(6) :: LABW
      INTEGER               :: ITYPE, I
      REAL(DOUBLE)          :: CONTR, TEGRAL
!-----------------------------------------------
!
      ! In BRITN?, LABV may be changed. 
      LABW(1:6)=LABV(1:6)
      !write(*,*)'IC,IR,IV,VCOEFF,LABV=',IC,IR,IV,VCOEFF,LABV(1:6)
      EMTTMP = 0.D0
      TEGRAL = 0.D0
! CYC: Discarding the Breit interaction according to the NDISCARDBR,
! LDISCARDBR, NDISBR, and LDISBR. 
! NOTE:
! The contributions involving only s and p orbitals are limited by only
! LDISCARDBR and LDISBR, whereas the others are limited also by
! NDISCARDBR and NDISBR. 
      IF (NDISCARDBR) THEN
        DO I = 1, 4
          IF (NKL(LABW(I)).GT.1 .AND. NP(LABW(I)).GT.NDISBR) RETURN
        ENDDO
      ENDIF
      IF (LDISCARDBR) THEN
        DO I = 1, 4
          IF (NKL(LABW(I)).GT.LDISBR) RETURN
        ENDDO
      ENDIF

      ITYPE = ABS(LABW(6))       
      IF (ITYPE .EQ. 1) THEN
        CALL BRINT1 (LABW(1), LABW(2),                                &
                     LABW(3), LABW(4),                                &
                     LABW(5), TEGRAL,IC,IR,IV)
      ELSEIF (ITYPE .EQ. 2) THEN
        CALL BRINT2 (LABW(1), LABW(2),                                &
                     LABW(3), LABW(4),                                &
                     LABW(5), TEGRAL,IC,IR,IV)
      ELSEIF (ITYPE .EQ. 3) THEN
        CALL BRINT3 (LABW(1), LABW(2),                                &
                     LABW(3), LABW(4),                                &
                     LABW(5), TEGRAL)
      ELSEIF (ITYPE .EQ. 4) THEN
        CALL BRINT4 (LABW(1), LABW(2),                                &
                     LABW(3), LABW(4),                                &
                     LABW(5), TEGRAL)
      ELSEIF (ITYPE .EQ. 5) THEN
        CALL BRINT5 (LABW(1), LABW(2),                                &
                     LABW(3), LABW(4),                                &
                     LABW(5), TEGRAL)
      ELSEIF (ITYPE .EQ. 6) THEN
        CALL BRINT6 (LABW(1), LABW(2),                                &
                     LABW(3), LABW(4),                                &
                     LABW(5), TEGRAL)
      ENDIF
      CONTR = VCOEFF*TEGRAL

      IF (LABW(6) .GT. 0) THEN
        EMTTMP = CONTR
      ELSE
!                ...It comes here only when ic=ir=1
!                   clue: rkco<-breid<-talk<-label(6,i)
        !To make sure the above statements, check:
        IF (IC.NE.1.OR.IR.NE.1) THEN
          WRITE(*,*)'IC,IR,IV,LABV=',IC,IR,IV,LABV(1:6)
          STOP 'Unexpected IC.NE.1.OR.IR.NE.1 in BREIT_D.f ...'
        ENDIF 
        NCORE = NCORE + 1
        ELSTO = ELSTO + CONTR 
      ENDIF

      RETURN
      END SUBROUTINE BREIT_D
