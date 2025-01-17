!***********************************************************************
!                                                                      *
      SUBROUTINE SPODMV(N, M, B, C)
!                                                                      *
!   Matrix-matrix product: C = AB.  A  sparse  representation of the   *
!   lower triangle of the  (NxN)  matrix  A  is read from the disk     *
!                                                                      *
!   This is an adaptation of  Andreas Stathopulos'  routine  SPSBMV,   *
!   and is specific to GRASP2 derivatives.                             *
!                                                                      *
!   Call(s) to: [AUXBLAS]: DINIT/SINIT;                                *
!               [SPBLAS]: DAXPYI/SAXPYI, DDOTI/SDOTI.                  *
!                                                                      *
!   F A Parpia and A Stathopoulos         Last revision: 13 Oct 1992   *
!   Block Version by Xinghong He          Last revision: 18 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE eigv_C
!CYC      USE Where_C
      USE fposition_C
      USE mpi_C,           ONLY:  myid, nprocs
      use symmatrix_mod,   ONLY:  IMCDF, NODENCOLS, NODECOLS
      USE MPI_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dinit_I
      USE posfile_I
      USE dmerge_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: M
      REAL(DOUBLE), DIMENSION(N,M), INTENT(IN) :: B
      REAL(DOUBLE), DIMENSION(N,M), INTENT(INOUT) :: C
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(N) :: IROW
!cjb  INTEGER :: MYID, NPROCS, NCF, NCFDUM, ICCUTDUM, MYIDDUM, NPROCSDUM, ICOL&
      INTEGER :: I, NCF, NCFDUM, ICCUTDUM, MYIDDUM, NPROCSDUM, ICOL,   &
           NELC, IR, IV, IREADH
      REAL(DOUBLE), DIMENSION(N) :: EMT
      REAL(DOUBLE) :: ELSTO, DIAG, DL, time0, time00, time1
      DATA IREADH/0/
!-----------------------------------------------
!
!     !...nposition+1 is the current position of the .res file
!     !...It is set in matrix and used in maneig, spodmv
!
!-----------------------------------------------------------------------
      !IF (MYID == 0) WRITE (6, *) 'Calling spodmv...'
      IREADH = IREADH + 1
      !IF (MYID .EQ. 0) PRINT *, 'CALLING SPODMV, IREADH=',IREADH
      NCF = N

!   Initialise the result matrix; note that this is specific to the
!   data structure of DVDSON
      time0 = MPI_WTIME()
      time00 = time0
      CALL DINIT (N*M, 0.D0, C, 1)
      time1 = MPI_WTIME()
      !IF (MYID .EQ. 0) PRINT *, 'DINIT array C (s) ...', time1-time0
      time0 = time1

      !...moved from maneig before "CALL GDVD (SPODMV..."
      CALL POSFILE (0, IMCDF, NPOSITION)

      READ (IMCDF) NCFDUM, ICCUTDUM, MYIDDUM, NPROCSDUM
      IF (NCF/=NCFDUM .OR. MYID/=MYIDDUM .OR. NPROCSDUM/=NPROCS) STOP &
         'spodmv: ncf read wrong'

!      DO ICOL = MYID + 1, N, NPROCS
      DO I = 1, NODENCOLS
         ICOL = NODECOLS(I)
         READ (IMCDF) NELC, ELSTO, (EMT(IR),IR=1,NELC), (IROW(IR),IR=1,NELC)
         DO IV = 1, M
            DIAG = C(ICOL,IV) + (EMT(NELC)-EAV)*B(ICOL,IV)
            CALL DMERGE (NELC - 1, B(1:N,IV), C(1:N,IV), IROW(1:N), &
                         EMT(1:N), B(ICOL,IV), DL)
            C(ICOL,IV) = DIAG + DL
         END DO
      END DO
      time1 = MPI_WTIME()
      !IF (MYID .EQ. 0) PRINT *, 'ReadH, Calculation (s) ...',time1-time0

      CALL gdsummpi (C, N*M)
      time0 = MPI_WTIME()
      !IF (MYID .EQ. 0) PRINT *, 'gdsumpi time (s):', time0-time1
      IF (MYID .EQ. 0) PRINT *, 'IReadH, SPODMV (s) ...',  &
                       IREADH, time0 - time00

      RETURN
      END SUBROUTINE SPODMV
