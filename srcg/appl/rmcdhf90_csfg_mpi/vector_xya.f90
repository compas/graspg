!***********************************************************************
!                                                                      *
      SUBROUTINE VECTOR_XYA
!                                                                      *
!   This subroutine vectorizes the two-demension arrays of NXAOPT,     *
!   XAOPT, NYAOPT, and YAOPT, to save memory                           *
!                                                                      *
!   Written by Chongyang Chen, Fudan University, Shanghai, Aug, 2024   *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE, LONG
      USE parameter_def,    ONLY: KEYORB, NNNW
      USE mpi_C
      USE orb_C,            ONLY: NW

      USE csfg_memory_man
      USE csfg_scf_C,       ONLY: NXCOFOPT, NXAOPT, XAOPT,  &
                                  NYCOFOPT, NYAOPT, YAOPT,  &
                                  NXCOFW,   NXAWRK, XAWRK,  &
                                  NYCOFW,   NYAWRK, YAWRK,  &
                                  NTOTXA,   NTOTYA, JXIPOS, &
                                  TMPXYA

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, IORB, NUMB, NUME, NUMT, NXY
      INTEGER :: JXIPOST(0:NNNW, NNNW) = 0
!-----------------------------------------------
! Vectorize NXAOPT and XAOPT
      NTOTXA = SUM(NXCOFOPT(1:NW))
      IF (NTOTXA .LT. 0) STOP 'Error! NTOTXA overflow in VECTORIZE_XYA'
      IF (MYID == 0) &
        WRITE(*, *)"Total XA-Term for all orbitals =", NTOTXA

! Alloc arrays
      CALL ALLOC(NXAWRK, NTOTXA, 'NXAWRK', 'VECTORIZE_XYA')
      CALL ALLOC( XAWRK, NTOTXA, ' XAWRK', 'VECTORIZE_XYA')

! Set NXCOFW
      NXCOFW(0) = 0
      DO I = 1, NW
        NXCOFW(I) = SUM(NXCOFOPT(1:I))
      ENDDO

! XAWRK will be calculated in SETCOF_NXY_CSFG by calling SETCOF_XY_CSFG,
! during SCF procedure

! Set NXAWRK
      NUMB = 0
      DO I = 1, NW
        NUMT = NXCOFOPT(I)
        NUME = NUMB + NUMT
        NXAWRK(NUMB+1:NUME) = NXAOPT(1:NUMT, I)
        NUMB = NUME
      ENDDO

! Re-set JXIPOS(0:NNNW, NNNW), which are employed to limit the range to
! speed up the binary search procedure.
      JXIPOST(:,:) = JXIPOS(:,:)
      DO I = 1, NW
        NUMB = NXCOFW(I-1)
        JXIPOST(0, I) = NUMB
        DO IORB = 1, NW
          JXIPOST(IORB, I) = NUMB + JXIPOST(IORB, I)
        ENDDO
      ENDDO
      JXIPOS(:,:) = JXIPOST(:,:)

!-----------------------------------------------------------------------
! Allocate the temp arrays for MPI_ALLREDUCE in SETCOF_CYC
! Deallocated in SCFMPI
      NXY = MAX(MAXVAL(NXCOFOPT(1:NW)), MAXVAL(NYCOFOPT(1:NW)))
      CALL ALLOC (TMPXYA, NXY, 'TMPXYA', 'VECTOR_XYA')

!-----------------------------------------------------------------------
! Vectorize NYAOPT and YAOPT
      NTOTYA = SUM(NYCOFOPT(1:NW))
      IF (NTOTYA .LT. 0) STOP 'Error! NTOTYA overflow in VECTORIZE_XYA'
      IF (MYID == 0) &
        WRITE(*, *)"Total YA-Term for all orbitals =", NTOTYA

! Alloc arrays
      CALL ALLOC(NYAWRK, NTOTYA, 'NYAWRK', 'VECTORIZE_XYA')
      CALL ALLOC( YAWRK, NTOTYA, ' YAWRK', 'VECTORIZE_XYA')

! Set NXCOFW
      NYCOFW(0) = 0
      DO I = 1, NW
        NYCOFW(I) = SUM(NYCOFOPT(1:I))
      ENDDO

! YAWRK will be calculated in SETCOF_NXY_CSFG by calling SETCOF_XY_CSFG,
! during SCF procedure

! Set NYAWRK
      NUMB = 0
      DO I = 1, NW
        NUMT = NYCOFOPT(I)
        NUME = NUMB + NUMT
        NYAWRK(NUMB+1:NUME) = NYAOPT(1:NUMT, I)
        NUMB = NUME
      ENDDO

!! Check
!      IF (MYID == 0) THEN
!        WRITE(31855, '(10X,214I10)')NXCOFOPT(1:NW)
!        WRITE(31856, '(10X,214I10)')NXCOFW(1:NW)
!        NUMB = 0
!        DO I = 1, NW
!          NUMT = NXCOFOPT(I)
!          NUME = NUMB + NUMT
!          WRITE(31855, *)"Iorb, NUMB, NUMT, NUME =", I, NUMB, NUMT, NUME
!          WRITE(31855, '(5000I10)')NXAOPT(1:NUMT, I)
!          WRITE(31855, '(2X,A8, 215I10)')"JXIPOS =", JXIPOS (0:NW,I)
!
!          WRITE(31856, *)"Iorb, NUMB, NUMT, NUME =", I, NUMB, NUMT, NUME
!          WRITE(31856, '(5000I10)')NXAWRK(NUMB+1:NUME)
!          WRITE(31856, '(2X,A8, 215I10)')"JXIPOS =", JXIPOST(0:NW,I) - NXCOFW(I-1)
!          WRITE(31856, '(2X,A8, 215I10)')"JXIPOST=", JXIPOST(0:NW,I)
!          NUMB = NUME
!        ENDDO 
!      ENDIF

      RETURN
      END SUBROUTINE VECTOR_XYA
