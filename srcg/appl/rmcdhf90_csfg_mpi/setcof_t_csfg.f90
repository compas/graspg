!***********************************************************************
!                                                                      *
      SUBROUTINE SETCOF_T_CSFG(EOL, IFLAG)
!                                                                      *
!   This subroutine prepares the data involving one-electron integrals *
!   Called by SETLAGmpi                                                *
!   Modified from SETCOF.F90, by Chongyang Chen, Fudan University,     *
!                                                Shanghai,  Dec 2021   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE, LONG
      USE parameter_def,    ONLY: KEYORB
      USE csfg_memory_man
      USE orb_C
      USE hblock_C
      USE hmat_C
      USE iounit_C
      USE MCPA_C
      USE mpi_C
      USE pos_C
      USE csfg_scf_C
      USE hmat_C,           ONLY: NELMNT
!CYC      USE rmcdhf_mem_C,     ONLY: INDX, ICLMN, COEFF, IMIN2_30, IMAX2_30
!      USE cons_C,           ONLY: EPS

      use csfg_tv_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dsubrs_I
!CYC      USE fco_I
!CYC      USE gco_I
!CYC      USE read1_mem_I
!CYC      USE read3_mem_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL  :: EOL, IFLAG
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      CHARACTER*7, PARAMETER :: MYNAME = 'SETVCOF'
      REAL, PARAMETER :: EPS = 1.0D-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: FIRST, LFOUND 
      INTEGER :: I, IA, IB, IC, IR, ITHIS, J, JBLOCK, &
                 LAB, LOC, IPOS, NTC
      INTEGER :: NBEG, NBEGCR, NUMCR, IOFF

      REAL(DOUBLE) :: UCFJ, SUMR, YKAB, XKAB, SUMT,  CONTR, TCOEFF
      REAL    :: nsecs, starttime, endtime

!=======================================================================
!   Generate DA coefficients; these arise from the one-electron
!   integrals
!=======================================================================
!   CPC 21, (1980) 207, Eq. (8): T_{rs}I_{ab}
!
!   IFLAG == F: Setup the sorted NDAOPT                                *
!               also the array sizes NDCOFOPT                          *
!   IFLAG == T: 
!        Obtain the coefficients involved in potentials calculation    *
!
      FIRST = .NOT. IFLAG
      IF (FIRST) THEN
        NDCOFOPT = 0
      ELSE
        DAOPT = 0.D0
      ENDIF
!
!=======================================================================
!   Generate DA coefficients; these arise from the one-electron
!   integrals.
!   CPC 21, (1980) 207, Eq. (8): T_{rs}I_{ab}
!=======================================================================

      IF (FIRST) GOTO 101

! TSUMG: allocated in match_labt.f90
! Initialization
      TSUMG(1:NCNTGT) = 0.D0
! Accumulate the onebody contribution
      STARTTIME = MPI_WTIME() 
      DO JBLOCK = 1, NBLOCK
        NBEG = SUM(NBPT(0:JBLOCK-1))
        NBEGCR = SUM(NBCRPT(0:JBLOCK-1))
        IOFF = 0
        DO I = 1, NBPT(JBLOCK)
          TCOEFF = TCOEFP(NBEG+I)
          NUMCR  = NCRPST(NBEG+I)
          DO J = 1, NUMCR
            LOC = NBEGCR + IOFF + J  
            IR = ICIRPT(1, LOC)
            IC = ICIRPT(2, LOC)
            CONTR = DSUBRS(EOL,IR,IC,JBLOCK)*TCOEFF
            ! off-diagonal contributions have double the weight
            IF (IR /= IC) CONTR = CONTR + CONTR

            ! / UCFJ is included within setallcof.f90
            ! SUM = 0.5D0 * SUM / UCFJ
            CONTR = 0.5d0 * CONTR
            ! Add up all the contributions from LAB integral:
            ! Locate LAB by using MATLABT
            TSUMG(MATLABT(LOC)) = TSUMG(MATLABT(LOC)) + CONTR
          ENDDO
          IOFF = IOFF + NUMCR
        ENDDO
        IF (IOFF /= NBCRPT(JBLOCK)) THEN
          WRITE(*,*)'JBLOCK, IOFF, NBCRPT(JBLOCK) =', &
                     JBLOCK, IOFF, NBCRPT(JBLOCK)
          STOP &
           'Error! Unexpected IOFF /= NBCRPT(JBLOCK) in settvcof_csfg..'
        ENDIF
      END DO 

      ENDTIME = MPI_WTIME()
      if (myid.eq.0) write(*,*) &
        "SETCOF_T_CSFG TSUMG time (s): ", ENDTIME - STARTTIME 

101   CONTINUE

      STARTTIME = MPI_WTIME()
      DO NTC = 1, NCNTGT
        LAB = LABALLT(NTC)
        IA = MOD(LAB, KEY)
        IB = LAB / KEY

        SUMT = TSUMG(NTC)
        ! Generate the da coefficients for both IA and IB orbitals.
        DO I = 1, 2
          IF (I == 1) THEN 
            J = IA
            ITHIS = IB
            UCFJ = UCF(J)
          ELSE
            J = IB
            ITHIS = IA
            UCFJ = UCF(J)
          ENDIF
          ! NDCOFOPT(J) is initialized at the beginning
          CALL LOCATE(ITHIS,NDCOFOPT(J),NDAOPT(:,J),LFOUND,IPOS)
          IF (FIRST) THEN
            ! Obtain the sorted array NDAOPT(:, J)
            IF (.NOT.LFOUND) THEN
              NDCOFOPT(J) = NDCOFOPT(J) + 1
              NDCOF = NDCOFOPT(J)
              IF (NDCOF > SIZE(NDAOPT,DIM=1)) THEN
                NDDIM = 2*NDDIM
                CALL RALLOC(NDAOPT, NDDIM, NW, 'NDAOPT', 'SETALLCOF')
              ENDIF
              NDAOPT(IPOS+2:NDCOF,J) = NDAOPT(IPOS+1:NDCOF-1,J)
              NDAOPT(IPOS+1,J) = ITHIS
            ENDIF
          ELSE
            ! Check
            !IF (.NOT.LFOUND) THEN
            !  WRITE(*, *)"Error, Unexpected .NOT.LFOUND ===03==="
            !  STOP "Error, Unexpected .NOT.LFOUND ===03==="
            !ENDIF
            ! DAOPT has been initiallized
            DAOPT(IPOS,J) = DAOPT(IPOS,J) + SUMT/UCFJ
          ENDIF
        ENDDO
      END DO

      ENDTIME = MPI_WTIME()  
      if (myid.eq.0 .and. (.not.first)) write(*,*) &
        "SETCOF_T_CSFG: Potentials (DAOPT) time (s) : ", &
        ENDTIME - STARTTIME
      if (myid.eq.0 .and. first) write(*,*) &
        "SETCOF_T_CSFG: NDAOPT, NDCOF (Lab) time (s): ", &
        ENDTIME - STARTTIME
      if (myid.eq.0) write(*,*)

      RETURN
      END SUBROUTINE SETCOF_T_CSFG
