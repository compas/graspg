!*******************************************************************
!                                                                  *
      SUBROUTINE EL2(JJA,JJB,JA,JB,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 04  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 2        *
!                                              N'2 = N2 + 2        *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1122,ITREXG,IXJTIK,PERKO2,       *
!                        RECO,RECO2,SIXJ,SPEAK                     *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, HALF, EPS
      USE m_C,             ONLY: JLIST, NPEEL
      USE orb_C,           ONLY: NAK, NKL
      USE trk_C
!CYC 
      USE csfg_decide_C,        ONLY: LDISCARDBR, LDISBR
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE reco_I
      USE reco2_I
      USE itrexg_I
      USE itrig_I
      USE ixjtik_I
      USE snrc_I
      USE coulom_I
      USE gg1122_I
      USE sixj_I
      USE cxk_I
      USE talk_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,ICOLBREI
!      DIMENSION J(2)
!      DIMENSION COND(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IAT,IA,IB,IBRD,IBRE,IP1,IP2,IG1,IG2,IKK,II,I2,I3, &
                 IFAZ,IFAZ1,IFAZFRCS,J12,JAA,JBB, &
                 KRA,L1,L2,N,ND1,ND2,NE1,NE2,NUP1,NU,MU
      INTEGER, DIMENSION(2) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
      REAL(DOUBLE)          :: AA,A1,AB,BB,QM1,QM2,QM3,QM4,RECC,SI
      REAL(DOUBLE), DIMENSION(12)    :: S
      REAL(DOUBLE), DIMENSION(12,20) :: COND
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      IF(JA > JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      CALL RECO(JAA,JBB,JBB,JBB,1,IAT)
      IF(IAT == 0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
!CYC:
      IF (ICOLBREI == 2) THEN
        IF (LDISCARDBR .AND. (NKL(IA).GT.LDISBR .OR.NKL(IB).GT.LDISBR)) RETURN
      ENDIF

      QM1=HALF
      QM2=HALF
      QM3=-HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      J(1)=ID1(3)
      J(2)=ID2(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      IP1=ITREXG(J(1),J(1),J(2),J(2),IKK)+1
      IF(IKK <= 0)RETURN
      IG1=IP1+IKK-1
! * * *                      * * *                      * * *
!     THE CASE 1122   + + - -
!
      IP2=ITREXG(J(1),J(2),J(1),J(2),IKK)+1
      IF(IKK <= 0) RETURN
      IG2=IP2+IKK-1
      IF (ICOLBREI == 2) THEN
        IS(1)=IA
        IS(2)=IA
        IS(3)=IB
        IS(4)=IB
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD <= 0)RETURN
        DO II=1,20
          COND(1,II) =ZERO
          COND(2,II) =ZERO
          COND(3,II) =ZERO
          COND(4,II) =ZERO
          COND(5,II) =ZERO
          COND(6,II) =ZERO
          COND(7,II) =ZERO
          COND(8,II) =ZERO
          COND(9,II) =ZERO
          COND(10,II)=ZERO
          COND(11,II)=ZERO
          COND(12,II)=ZERO
        END DO
      END IF
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
        IF (ICOLBREI == 1) THEN
          CALL COULOM(L1,L1,L2,L2,ID1(5),ID1(5),ID2(5),ID2(5),KRA,A1)
          IF(DABS(A1) < EPS) CYCLE
          A1=-HALF*A1
        END IF
        AB=ZERO
          DO I3=IP1,IG1,2
            J12=(I3-1)/2
            CALL RECO2(JAA,JBB,J12*2,0,IAT,RECC)
            IF(IAT /= 0) THEN
              IF(IXJTIK(J(1),J(2),KRA*2,J(2),J(1),J12*2) /= 0) THEN
                CALL GG1122(J12,J12,QM1,QM2,QM3,QM4,AA)
                IF(DABS(AA) > EPS) THEN
                  CALL RECO2(JAA,JBB,J12*2,1,IAT,RECC)
                  AA=AA*RECC
                  CALL SIXJ(J(1),J(2),KRA*2,J(2),J(1),J12*2,0,SI)
                  AA=AA*SI*DSQRT(DBLE(I3))
                  IFAZ=IK1(3)+IK2(3)+KRA*2+J12*2
                  IF((IFAZ/4)*4 /= IFAZ)AA=-AA
                  AB=AB+AA
                END IF
              END IF
            END IF
          END DO
!
!       TRANSFORM FANO & RACAH PHASE CONVENTION
!       TO CONDON & SHORTLEY PHASE CONVENTION
!
        IFAZFRCS=1
        IFAZ1=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
        IF((IFAZ1/4)*4 /= IFAZ1)IFAZFRCS=-IFAZFRCS
!
        IF (ICOLBREI == 1) THEN
          BB=A1*AB*DBLE(IFAZFRCS)
          IF(DABS(BB) > EPS)CALL SPEAK(JJA,JJB,IA,IA,IB,IB,KRA,BB)
        ELSE IF (ICOLBREI == 2) THEN
          NU=KRA
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+1) /= 0) .AND.   &
               (ITRIG(KS(2),KS(4),NU+NU+1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(NU > 0) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA+1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU-1) /= 0) .AND.  &
               (ITRIG(KS(2),KS(4),NU+NU-1) /= 0)) THEN
              N=(NU-ND1)/2+1
              IF(N <= ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO MU = 1,4
                  COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
                END DO
              END IF
            END IF
          END IF
          NU=KRA-1
          IF(((NU-ND1)/2)*2 == (NU-ND1)) THEN
            IF((ITRIG(KS(1),KS(3),NU+NU+3) /= 0) .AND.   &
               (ITRIG(KS(2),KS(4),NU+NU+3) /= 0)) THEN
              IF(NU >= 0) THEN
                N=(NU-ND1)/2+1
                IF(N < ND2) THEN
                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                  DO MU = 1,12
                    COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      IF (ICOLBREI .EQ. 2) THEN
        DO N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N == ND2) CYCLE
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
        END DO
      END IF
      RETURN
      END SUBROUTINE EL2
