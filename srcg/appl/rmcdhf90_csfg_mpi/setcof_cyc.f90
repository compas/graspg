!***********************************************************************
!                                                                      *
      SUBROUTINE SETCOF_CYC(EOL, J)
!                                                                      *
!   This  subroutine reloads the coefficients and orbital pointers,    *
!   which are saved by setcof.f90, called by:                          *
!            SETLAGmpi / IMPROVmpi :: COFPOTmpi                        *
!                                                                      *
!   Written by Chongyang Chen, Fudan University, Shanghai, Dec, 2021   *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:21:02   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE, LONG
      USE parameter_def,    ONLY: KEYORB
      USE csfg_memory_man
      USE orb_C
      USE mpi_C
      USE csfg_scf_C

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J
      LOGICAL  :: EOL
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      TMPXYA is declared in csfg_scf_C, allocated in VECTOR_XYA,
!      deallocated in SCFMPI
!      REAL(DOUBLE), DIMENSION(:), POINTER :: TMPXYA 
!-----------------------------------------------

      NDCOF = NDCOFOPT(J)
      NXCOF = NXCOFOPT(J)
      NYCOF = NYCOFOPT(J)

      IF (SIZE(DA) < NDCOF) THEN
        CALL RALLOC(DA, NDCOF, 'DA', 'SETCOFCYC')
        CALL RALLOC(NDA,NDCOF, 'NDA', 'SETCOFCYC')
      ENDIF
      ! DA is reduced in DACON as done in GRASP2018
      DA(1:NDCOF) = DAOPT(1:NDCOF, J)
      NDA(1:NDCOF) = NDAOPT(1:NDCOF, J)

      IF (SIZE(XA) < NXCOF) THEN
        CALL RALLOC(XA, NXCOF,'XA', 'SETCOFCYC')
        CALL RALLOC(NXA, NXCOF, 'NXA', 'SETCOFCYC')
      ENDIF
!      XA(1:NXCOF) = XAOPT(1:NXCOF, J)
!      NXA(1:NXCOF) = NXAOPT(1:NXCOF, J)
      XA (1:NXCOF) = XAWRK (NXCOFW(J-1)+1 : NXCOFW(J))
      NXA(1:NXCOF) = NXAWRK(NXCOFW(J-1)+1 : NXCOFW(J))

      ! Reduce XA as NXA are same for all nodes now
      !CALL ALLOC (TMPXYA, NXCOF, 'TMPXYA', 'SETCOFCYC')
      CALL MPI_Allreduce (XA,TMPXYA,NXCOF,MPI_DOUBLE_PRECISION, &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
      XA(1:NXCOF) = TMPXYA(1:NXCOF)
      !CALL DALLOC(TMPXYA, 'TMPXYA', 'SETCOFCYC')

      IF (SIZE(YA) < NYCOF) THEN
        CALL RALLOC(YA, NYCOF, 'YA', 'SETCOFCYC')
        CALL RALLOC(NYA, NYCOF, 'NYA', 'SETCOFCYC')
      ENDIF
      !YA(1:NYCOF) = YAOPT(1:NYCOF, J)
      !NYA(1:NYCOF) = NYAOPT(1:NYCOF, J)
      YA(1:NYCOF)  = YAWRK (NYCOFW(J-1)+1 : NYCOFW(J))
      NYA(1:NYCOF) = NYAWRK(NYCOFW(J-1)+1 : NYCOFW(J))

      ! Reduce YA as NYA are same for all nodes now
      !CALL ALLOC (TMPXYA, NYCOF, 'TMPXYA', 'SETCOFCYC')
      CALL MPI_Allreduce (YA,TMPXYA,NYCOF,MPI_DOUBLE_PRECISION, &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
      YA(1:NYCOF) = TMPXYA(1:NYCOF)
      !CALL DALLOC(TMPXYA, 'TMPXYA', 'SETCOFCYC')

      RETURN
      END SUBROUTINE SETCOF_CYC
