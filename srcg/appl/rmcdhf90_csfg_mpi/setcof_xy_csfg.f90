!=====================================================================**
!                                                                      *
      SUBROUTINE SETCOF_XY_CSFG(IFLAG, LAB0, K, SUMV)
!                                                                      *
!   This  subroutine  sets  up the coefficients and orbital pointers   *
!   for the direct and exchange potentials for orbitals involved in    *
!   LAB0. Copied parts of the orginal GRASP SETCOF by C. Y. Chen 2021  *
!                                                                      *
!   Call(s) to: [LIB92]: alloc, dalloc.                                *
!               [RSCF92]: alcsca, dsubrs.                              *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 21 Dec 1992   *
!   Modified by Xinghong He                 Last update: 21 Dec 1997   *
!   Modified by G. Gaigalas               Last revision: 07 Sep 2021   *
!      Spin-angular coefficiants have been put on memory               *

!=====================================================================**
!   Sorted version, modified by Chongyang Chen, Fudan University,      *
!                                               Shanghai, Dec. 2021    *
!   CSFG   version, modified by Chongyang Chen, Fudan University,      *
!                                               Shanghai, Dec. 2023    *
!   Using vectors of NXAWRK, XAWRK, NYAWRK, YAWRK, C. Y. Chen          *
!                                               Shanghai, Aug. 2024    *
!                                                                      *
!   IFLAG == F: Setup the sorted NXAOPT, NYAOPT arrays;                *
!               also the array sizes NXCOFOPT, NYCOFOPT                *
!   IFLAG == T: Obtain the coefficients for all orbitals               *
!                                                                      *
!=====================================================================**
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE, LONG, BYTE
      USE parameter_def,    ONLY: KEYORB
      USE csfg_memory_man
      USE orb_C
      USE hblock_C
      USE hmat_C
      USE iounit_C
      USE MCPA_C
      USE mpi_C
      USE pos_C

      use csfg_scf_C
      use csfg_tv_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE locateall_I
      USE locate_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL, INTENT(IN)        :: IFLAG
      INTEGER, INTENT(IN)        :: LAB0
      INTEGER, INTENT(IN)        :: K
      REAL(DOUBLE), INTENT(IN)   :: SUMV
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      REAL, PARAMETER :: EPS = 1.0D-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER       :: INDEXS(4)
      REAL(DOUBLE)  :: UCFJ, SUM 
      INTEGER       :: IPOS, ILABWORK
      LOGICAL       :: LFOUND, FIRST
!=======================================================================
!   Initializations
!=======================================================================
      REAL    :: nsecs
      INTEGER :: NDIFFJ, JARR(4), IJLOOP, J, LAB
      INTEGER :: I, IFOUND, IIND, IL, INDIND, IORB, IRANK, ITHIS, &
                 ITHIS2, IYO1, IYO2, LOC1, LOC2
!      INTEGER :: IA, IB, IC, ID
!=======================================================================
!
!   IFLAG == F: Setup the sorted NDAOPT, NXAOPT, NYAOPT arrays;        *
!               also the array sizes NDCOFOPT, NXCOFOPT, NYCOFOPT      *
!   IFLAG == T: Obtain the coefficients for all orbitals               *
!
      FIRST = .NOT. IFLAG
      LAB = LAB0
      !===                       k
      !=== Decode the labels of R (abcd)
      INDEXS(4) = MOD(LAB,KEY)  ! ID?
      LAB = LAB/KEY
      INDEXS(2) = MOD(LAB,KEY)  ! IB?
      LAB = LAB/KEY
      INDEXS(3) = MOD(LAB,KEY)  ! IC?
      INDEXS(1) = LAB/KEY       ! IA?
!      ID = INDEXS(1); IB = INDEXS(2); IC = INDEXS(3); IA = INDEXS(4)   
       
      ! Count the number of different orbitals for LAB.
      JARR = -1
      JARR(1) = INDEXS(1)
      NDIFFJ = 1
      DO I = 2, 4
        IF (COUNT(JARR == INDEXS(I)) == 0) THEN
          NDIFFJ = NDIFFJ + 1
          JARR(NDIFFJ) = INDEXS(I)
        ENDIF
      ENDDO

      ! Calculate the coefficients for the NDIFFJ orbitals
      IJLOOP0: DO IJLOOP = 1, NDIFFJ
        J = JARR(IJLOOP)

!CYC 2024/03 FIX
        IF (LSKIPPOT(J)) CYCLE

        ! NO matter FIRST is true or false, the followings are right
        NYCOF = NYCOFOPT(J)

        UCFJ = UCF(J)
        IF (IFLAG) SUM = 0.5D0*SUMV/UCFJ

        !=== Determine the number of indices that match
        IRANK = COUNT(INDEXS==J)

        SELECT CASE (IRANK)
        CASE (1)

!=======================================================================
!   One matching index: exchange potential contribution
!=======================================================================
           DO IIND = 1, 4
              IF (INDEXS(IIND) /= J) CYCLE   ! at least one
              IL = IIND + 2
              IF (IL > 4) IL = IL - 4
              IORB = INDEXS(IL)
              IL = IIND + 1
              IF (IL > 4) IL = IL - 4
              IYO1 = INDEXS(IL)
              IL = IIND + 3
              IF (IL > 4) IL = IL - 4
              IYO2 = INDEXS(IL)
              ITHIS = ((IORB*KEY + IYO2)*KEY + IYO1)*KEY + K
              EXIT
           END DO

           !=== Check ithis against the previously recorded ===
           IF (FIRST) THEN
             CALL LOCATE(ITHIS, JXIPOS(IORB,J), &
                         NXAIORB(:,IORB,J), LFOUND, IPOS)
             IF (.NOT. LFOUND) THEN
               JXIPOS(IORB,J) = JXIPOS(IORB,J) + 1
               IF (JXIPOS(IORB,J) > SIZE(NXAIORB,DIM=1)) THEN
                 NDIMX = 2*NDIMX
                 CALL RALLOC(NXAIORB, NDIMX, NW, NW, &
                             'NXAIOB','SETCOF_XY_CSFG')
               ENDIF
               IPOS = IPOS + 1
               NXAIORB(IPOS+1:JXIPOS(IORB,J),IORB,J) = &
                       NXAIORB(IPOS:JXIPOS(IORB,J)-1,IORB,J)
               NXAIORB(IPOS,IORB,J) = ITHIS
             ENDIF
           ELSE
             ! To speed up the search, locate the range by IORB
             CALL LOCATEALL(ITHIS, JXIPOS(IORB-1,J)+1, JXIPOS(IORB,J), &
                            NXAWRK, LFOUND, IPOS)
             ! Check
             !IF (.NOT.LFOUND) THEN
             !  WRITE(*, *)"Error, Unexpected .NOT.LFOUND ===04==="
             !  STOP "Error, Unexpected .NOT.LFOUND ===04==="
             !ENDIF
             XAWRK(IPOS) = XAWRK(IPOS) + SUM
           ENDIF

        CASE (2)

!=======================================================================
!   Two matching indices: either direct or exchange potential
!   contribution
!=======================================================================
           IFOUND = 0
           DO IIND = 1, 4
              IF (INDEXS(IIND) /= J) CYCLE
              IF (IFOUND == 0) THEN
                 LOC1 = IIND
                 IFOUND = IFOUND + 1
              ELSE IF (IFOUND == 1) THEN
                 LOC2 = IIND
                 EXIT
              ENDIF
           END DO

           IF (LOC2 - LOC1 == 2) THEN
!
!   Direct contribution
!
              !=== Find ithis ===
              IL = LOC1 + 3
              IF (IL > 4) IL = IL - 4
              IYO2 = INDEXS(IL)
              IL = LOC1 + 1
              IF (IL > 4) IL = IL - 4
              IYO1 = INDEXS(IL)
              ITHIS = (IYO2*KEY + IYO1)*KEY + K

              IF (FIRST) THEN
                CALL LOCATE(ITHIS,NYCOF,NYAOPT(:,J),LFOUND,IPOS)
                IF (.NOT.LFOUND) THEN
                  !=== Not found, add an item ===
                  NYCOF = NYCOF + 1
                  IF (NYCOF > SIZE(NYAOPT,DIM=1)) THEN
                    NYDIM = 2*NYDIM
                    CALL RALLOC(NYAOPT, NYDIM, NW, &
                               'NYAOPT','SETCOF_XY_CSFG')
                  ENDIF
                  NYAOPT(IPOS+2:NYCOF,J) = NYAOPT(IPOS+1:NYCOF-1,J)
                  NYAOPT(IPOS+1,J) = ITHIS

                  NYCOFOPT(J) = NYCOF
                ENDIF
              ELSE
                CALL LOCATEALL(ITHIS, NYCOFW(J-1)+1, NYCOFW(J), &
                               NYAWRK, LFOUND, IPOS)
                ! Check
                !IF (.NOT.LFOUND) THEN
                !  WRITE(*, *)"Error, Unexpected .NOT.LFOUND ===05==="
                !  STOP "Error, Unexpected .NOT.LFOUND ===05==="
                !ENDIF
                YAWRK(IPOS) = YAWRK(IPOS) + SUM + SUM
              ENDIF
           ELSE
!
!   Exchange contribution
!
              !=== Find ithis ===
              IL = LOC1 + 2
              IF (IL > 4) IL = IL - 4
              IORB = INDEXS(IL)
              IL = LOC1 + 1
              IF (IL > 4) IL = IL - 4
              IYO1 = INDEXS(IL)
              IL = LOC1 + 3
              IF (IL > 4) IL = IL - 4
              IYO2 = INDEXS(IL)
              ITHIS = ((IORB*KEY + IYO2)*KEY + IYO1)*KEY + K

              IF (FIRST) THEN
                CALL LOCATE(ITHIS, JXIPOS(IORB,J), &
                            NXAIORB(:,IORB,J), LFOUND, IPOS)
                IF (.NOT. LFOUND) THEN
                  JXIPOS(IORB,J) = JXIPOS(IORB,J) + 1
                  IF (JXIPOS(IORB,J) > SIZE(NXAIORB,DIM=1)) THEN
                    NDIMX = 2*NDIMX
                    CALL RALLOC(NXAIORB, NDIMX, NW, NW, &
                               'NXAIOB', 'SETCOF_XY_CSFG')
                  ENDIF
                  IPOS = IPOS + 1
                  NXAIORB(IPOS+1:JXIPOS(IORB,J),IORB,J) = &
                          NXAIORB(IPOS:JXIPOS(IORB,J)-1,IORB,J)
                  NXAIORB(IPOS,IORB,J) = ITHIS
                ENDIF
              ELSE
                CALL LOCATEALL(ITHIS, JXIPOS(IORB-1,J)+1, JXIPOS(IORB,J), &
                               NXAWRK, LFOUND, IPOS)
                ! Check
                !IF (.NOT.LFOUND) THEN
                !  WRITE(*, *)"Error, Unexpected .NOT.LFOUND ===02B=="
                !  STOP "Error, Unexpected .NOT.LFOUND ===02B=="
                !ENDIF
                XAWRK(IPOS) = XAWRK(IPOS) + SUM + SUM ! XAWRK has been initiallized
              ENDIF
           ENDIF

        CASE (3)

!=======================================================================
!   Three matching indices: direct and exchange potential contributions
!=======================================================================

           !=== Find ithis AND ithis2
           ITHIS = -911
           ITHIS2 = -911
           DO IIND = 1, 4
              IF (INDEXS(IIND) == J) CYCLE
              INDIND = INDEXS(IIND)
              IYO2 = INDIND
              IYO1 = J
              ITHIS = (IYO2*KEY + IYO1)*KEY + K
              IORB = INDIND
              IYO1 = J
              IYO2 = J
              ITHIS2 = ((IORB*KEY + IYO2)*KEY + IYO1)*KEY + K
           END DO

           IF (ITHIS==(-911) .OR. ITHIS2==(-911)) STOP 'ithis2'
!
! Direct potential
! 
           IF (FIRST) THEN
              CALL LOCATE(ITHIS,NYCOF,NYAOPT(:,J),LFOUND,IPOS)
              IF (.NOT.LFOUND) THEN
                ! Not found, add an item ===
                NYCOF = NYCOF + 1
                IF (NYCOF > SIZE(NYAOPT,DIM=1)) THEN
                  NYDIM = 2*NYDIM
                  CALL RALLOC (NYAOPT, NYDIM, NW, &
                              'NYAOPT', 'SETCOF_XY_CSFG')
                ENDIF
                NYAOPT(IPOS+2:NYCOF,J) = NYAOPT(IPOS+1:NYCOF-1,J)
                NYAOPT(IPOS+1,J) = ITHIS

                NYCOFOPT(J) = NYCOF
              ENDIF
           ELSE
              CALL LOCATEALL(ITHIS, NYCOFW(J-1)+1, NYCOFW(J), &
                             NYAWRK, LFOUND, IPOS)
              ! Check
              !IF (.NOT.LFOUND) THEN
              !  WRITE(*, *)"Error, Unexpected .NOT.LFOUND ===07==="
              !  STOP "Error, Unexpected .NOT.LFOUND ===07==="
              !ENDIF
              YAWRK(IPOS) = YAWRK(IPOS) + SUM + SUM ! YAWRK has been initiallized
           ENDIF

!
! Exchange potential
! 
           IF (FIRST) THEN
             CALL LOCATE(ITHIS2, JXIPOS(IORB,J), &
                         NXAIORB(:,IORB,J), LFOUND, IPOS)
             IF (.NOT. LFOUND) THEN
               JXIPOS(IORB,J) = JXIPOS(IORB,J) + 1
               IF (JXIPOS(IORB,J) > SIZE(NXAIORB,DIM=1)) THEN
                 NDIMX = 2*NDIMX
                 CALL RALLOC(NXAIORB,NDIMX,NW,NW,'NXAIOB','SETALLCOF')
               ENDIF
               IPOS = IPOS + 1
               NXAIORB(IPOS+1:JXIPOS(IORB,J),IORB,J) = &
                       NXAIORB(IPOS:JXIPOS(IORB,J)-1,IORB,J)
               NXAIORB(IPOS,IORB,J) = ITHIS2
             ENDIF
           ELSE
             CALL LOCATEALL(ITHIS2, JXIPOS(IORB-1,J)+1, JXIPOS(IORB,J), &
                            NXAWRK, LFOUND, IPOS)
             ! Check
             !IF (.NOT.LFOUND) THEN
             !  WRITE(*, *)"Error, Unexpected .NOT.LFOUND ===02B=="
             !  STOP "Error, Unexpected .NOT.LFOUND ===02B=="
             !ENDIF
             XAWRK(IPOS) = XAWRK(IPOS) + SUM  ! XAWRK has been initiallized
           ENDIF

        CASE (4)

!=======================================================================
!   Four matching indices: direct potential contribution
!=======================================================================

           !=== Find ithis
           IYO2 = J
           IYO1 = J
           ITHIS = (IYO2*KEY + IYO1)*KEY + K

           IF (FIRST) THEN
              CALL LOCATE(ITHIS,NYCOF,NYAOPT(:,J),LFOUND,IPOS)
              IF (.NOT.LFOUND) THEN
                ! Not found, add an item ===
                NYCOF = NYCOF + 1
                IF (NYCOF > SIZE(NYAOPT,DIM=1)) THEN
                  NYDIM = 2*NYDIM
                  CALL RALLOC(NYAOPT, NYDIM, NW, &
                             'NYAOPT', 'SETCOF_XY_CSFG')
                ENDIF
                NYAOPT(IPOS+2:NYCOF,J) = NYAOPT(IPOS+1:NYCOF-1,J)
                NYAOPT(IPOS+1,J) = ITHIS

                NYCOFOPT(J) = NYCOF
              ENDIF
           ELSE
              CALL LOCATEALL(ITHIS, NYCOFW(J-1)+1, NYCOFW(J), &
                             NYAWRK, LFOUND, IPOS)
              ! Check
              !IF (.NOT.LFOUND) THEN
              !  WRITE(*, *)"Error, Unexpected .NOT.LFOUND ===09==="
              !  STOP "Error, Unexpected .NOT.LFOUND ===09==="
              !ENDIF
              YAWRK(IPOS) = YAWRK(IPOS) + 4.0*SUM ! YAWRK has been initiallized
           ENDIF
        END SELECT

      ENDDO IJLOOP0  ! The FOUR orbitals involving LAB
!
      END SUBROUTINE SETCOF_XY_CSFG
!

