!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      *
!     Obtain the occupitation number of all the normal CSFs.           *
!                                                                      *
!     Chongyang Chen, Fudan University,                  Dec 2023      *
!                                                                      *
!cyc====================================================================

      subroutine csfg_iqa
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def,   ONLY: NNNW
      use csfg_memory_man
      use orb_C, ONLY: NW, IQA
      use symexpand_mod, nblock=>nblocksym
      use symmatrix_mod, ONLY: NTYPE 
      use mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, K, L, ICSF, ICSFG, ICGP


      ALLOCATE(IQA_CSF(NW, SUM(TotCSFs_perblock(3,1:NBLOCK))))
      IQA_CSF = 0

      !IF (MYID == 0) WRITE(*,*)'NONSYMINF =', NONSYMINF

      ICSF = 0 
      DO ICSFG = 1, NCFGTOT
        !write(5000,'(i8,200I3)')icsfg, IQA(1:NW, ICSFG)
        SELECT CASE (NTYPE(1,ICSFG))
        CASE (1)
          ICSF = ICSF + 1
          IQA_CSF(1:NW, ICSF) = IQA(1:NW, ICSFG)
        CASE (2)
          DO I = NTYPE(3,ICSFG), NTYPE(4,ICSFG)
            ICSF = ICSF + 1
            IQA_CSF(1:nonsyminf,ICSF) = IQA(1:nonsyminf,ICSFG)
            IQA_CSF(I,ICSF) = 1
          ENDDO
        CASE (3)
          DO I = NTYPE(3,ICSFG), NTYPE(4,ICSFG)
            DO J = NTYPE(5,ICSFG), NTYPE(6,ICSFG)
              ICSF = ICSF + 1
              IQA_CSF(1:nonsyminf,ICSF) = IQA(1:nonsyminf,ICSFG)
              IQA_CSF(I,ICSF) = 1
              IQA_CSF(J,ICSF) = 1 
            ENDDO
          ENDDO 
        CASE (4)
          DO I = NTYPE(3,ICSFG), NTYPE(6,ICSFG) - 1
            DO J = I + 1, NTYPE(6,ICSFG)
              ICSF = ICSF + 1
              IQA_CSF(1:nonsyminf,ICSF) = IQA(1:nonsyminf,ICSFG)
              IQA_CSF(I,ICSF) = 1
              IQA_CSF(J,ICSF) = 1 
            ENDDO
          ENDDO
        CASE (5)
          DO I = NTYPE(3,ICSFG), NTYPE(4,ICSFG)
            ICSF = ICSF + 1
            IQA_CSF(1:nonsyminf,ICSF) = IQA(1:nonsyminf,ICSFG)
            IQA_CSF(I,ICSF) = 2
          ENDDO
        END SELECT
      ENDDO 

      RETURN
      END SUBROUTINE CSFG_IQA
