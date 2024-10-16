      SUBROUTINE COFPOTmpi(EOL, J, NPTS)
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE pote_C
      USE mpi_C
      USE orb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setcof_cyc_I
      USE ypot_I
      USE xpot_I
      USE lagcon_I
      USE dacon_I
      USE POTE_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J, I
      INTEGER  :: NPTS
      LOGICAL  :: EOL
      REAL(DOUBLE), DIMENSION(NPTS) :: TMPMPI
!-----------------------------------------------------------------------
! CYC: Reload the saved data
!      CALL SETCOF (EOL, J)
      CALL SETCOF_CYC(EOL,J)

      CALL YPOT (J)
      CALL XPOT (J)
      CALL LAGCON (J, NPROCS)
! CYC May, 2022
!      CALL DACON
      CALL DACON (J)
! CYC May, 2022

      CALL MPI_Allreduce (YP, tmpmpi, npts, MPI_DOUBLE_PRECISION,  &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL dcopy (npts, tmpmpi, 1, YP, 1)

!      WRITE(9025,'(A, I4, 2X, I2, A2, 2X, I6,2X, A)')'J = ',  &
!            J, NP(J), NH(J), NPTS, 'YP:'
!      WRITE(9025, '(1P,5(D22.14))')(YP(I),I=1,NPTS)
!      WRITE(9025,*)
      CALL MPI_Allreduce (XP, tmpmpi, npts, MPI_DOUBLE_PRECISION,  &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL dcopy (npts, tmpmpi, 1, XP, 1)
!      WRITE(9025,'(A, I4, 2X, I2, A2, 2X, I6,2X, A)')'J = ',  &
!            J, NP(J), NH(J), NPTS, 'XP:'
!      WRITE(9025, '(1P,5(D22.14))')(XP(I),I=1,NPTS)
!      WRITE(9025,*)

      CALL MPI_Allreduce (XQ, tmpmpi, npts, MPI_DOUBLE_PRECISION,  &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL dcopy (npts, tmpmpi, 1, XQ, 1)
!      WRITE(9025,'(A, I4, 2X, I2, A2, 2X, I6,2X, A)')'J = ',  &
!            J, NP(J), NH(J), NPTS, 'XQ:'
!      WRITE(9025, '(1P,5(D22.14))')(XQ(I),I=1,NPTS)
!      WRITE(9025,*)

      RETURN
      END SUBROUTINE COFPOTmpi
