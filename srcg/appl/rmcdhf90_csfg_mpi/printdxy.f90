
      SUBROUTINE PRINTDXY(EOL)

      USE vast_kind_param, ONLY:  DOUBLE
      USE csfg_memory_man
      USE blkidx_C
      USE default_C
      USE def_C
      USE debug_C
      USE eigv_C
      USE fixd_C
      USE hblock_C
      USE iounit_C
      USE lagr_C
      USE MCPA_C
      USE mpi_C
      USE pos_c
      USE peav_C
      USE ORB_C
      USE orba_C
      USE csfg_scf_C
      USE SYMA_C
      USE STAT_C
      USE ORTHCT_C

      Use csfg_tv_C

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: EOL
      INTEGER :: I, J


      DO J = 1, NW
        CALL SETCOF_CYC(EOL, J)
        ! Potentials associated with T-coefficients:
        write(41001, '(A, 2I8)')'J, NDCOF =', J, NDCOF
        write(41001, '(100I22)')      (I, I = 1, NDCOF)
        write(41001, '(100I22)')      (NDA(I), I = 1, NDCOF)
        write(41001, '(100(1pe22.14))') (DA(I), I = 1, NDCOF)
        write(41001, *)

        write(41002, '(A, 2I8)')'J, NYCOF =', J, NYCOF
        write(41002, '(100I22)')      (I, I = 1, NYCOF)
        write(41002, '(100I22)')      (NYA(I), I = 1, NYCOF)
        write(41002, '(100(1pe22.14))') (YA(I), I = 1, NYCOF)
        write(41002, *)

        write(41003, '(A, 2I8)')'J, NXCOF =', J, NXCOF
        write(41003, '(100I22)')      (I, I = 1, NXCOF)
        write(41003, '(100I22)')      (NXA(I), I = 1, NXCOF)
        write(41003, '(100(1pe22.14))') (XA(I), I = 1, NXCOF)
        write(41003, *)
      ENDDO
      CLOSE(41001)
      CLOSE(41002)
      CLOSE(41003)
      RETURN
      END SUBROUTINE  
