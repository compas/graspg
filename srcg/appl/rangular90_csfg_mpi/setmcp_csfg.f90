!***********************************************************************
!                                                                      *
      SUBROUTINE SETMCP_CSFG (IFLAG, NCORE, IDBLK, IBLK, FILEHEAD)
!                                                                      *
!======================================================================*
!                                                                      *
! IFLAG == 1:                                                          *
!                                                                      *
! Open and check the mcpXXX.30 which stores interacting pairs          *
!                                                                      *
! Open files ${FILEHEAD}.csfg[1,2,3,4,5],                              *
!            ${FILEHEAD}.csfg[11,12,13,14,15],                         *
!            ${FILEHEAD}.csfg[22,23,24,25],                            *
!            ${FILEHEAD}.csfg[33,34,35],                               *
!            ${FILEHEAD}.csfg[44,45], and                              *
!            ${FILEHEAD}.csfg55                                        *
! which store the spin-angular coefficients for H-matrix calculations  *
!                                                                      *
!======================================================================*
! IFLAG == 2:                                                          *
!                                                                      *
! Write the heading lines of every block.                              *
!                                                                      *
!======================================================================*
! IFLAG == 3:                                                          *
!                                                                      *
! Write the end lines of every block.                                  *
!                                                                      *
!======================================================================*
! IFLAG == 4:                                                          *
!                                                                      *
! Close the files.                                                     *
!                                                                      *
!======================================================================*
!                                                                      *
! Written by Chongyang Chen, Fudan university, Shanghai,     Oct 2023  *
!                                                                      *
!======================================================================*
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT, GETYN, OPENFL.                        *
!               [GENMCP]: GETINF.                                      *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 08 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 30 Jun 1998   *
!                                                                      *
!   Used by mcpvu, mcpmpivu                                            *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE mcp_C,     ONLY: KMAX
      USE orb_C,     ONLY: NW, NKJ, NCF
      USE hblock_C,  ONLY: NBLOCK, NCFBLK
      USE iounit_C,  ONLY: ISTDE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I
      USE openfl_I
      USE mpi_C,     ONLY: MYID, NPROCS, IERR

      Use symexpand_mod
      use symmatrix_mod, ONLY: NCFTr
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!     INTEGER, INTENT(IN) :: MYID
!     INTEGER, INTENT(IN) :: NPROCS
      INTEGER, INTENT(IN) :: IFLAG, IBLK
      INTEGER, INTENT(IN) :: NCORE
      CHARACTER, INTENT(IN) :: FILEHEAD*(*)
      CHARACTER, INTENT(IN) :: IDBLK(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, LNG, LCK, I
      LOGICAL :: FOUND, FOUND1, GETYN, YES
      CHARACTER :: CK*2

      INTEGER :: FN(20)
      DATA FN/1,11,12,13,14,15,2,22,23,24,25,3,33,34,35,4,44,45,5,55/
!-----------------------------------------------

      SELECT CASE (IFLAG)
      CASE (1)
        LNG = LEN_TRIM(FILEHEAD)
        !DO I = 1, 20
        !  K = FN(I)
        !  CALL CONVRT (K, CK, LCK)
        !  !CALL OPENFL (2000+K, 'csfg_'//FILEHEAD(1:LNG)//'.'//CK(1:LCK), &
        !  !             'FORMATTED', 'UNKNOWN', IERR)
        !  !WRITE (2000+K, *) 'MCP', NBLOCK, NCFBLK(1:NBLOCK), MYID, NPROCS
        !  IF (IERR == 0) CYCLE

        !  WRITE (ISTDE, *) 'Error when opening the mcp files'
        !  STOP 'Error when opening the mcp files ...'
        !ENDDO

! Determine KMAX; this is the number of  .mcp  files for the
! two-electron integrals

        KMAX = 0
        DO K = 1, NW
           KMAX = MAX(KMAX,NKJ(K))
        END DO
!
! Open csfg_mcpXXX.30 and csfg_mcpXXX.31
!
        LNG = LEN_TRIM(FILEHEAD)
        DO K = 30, 32
          CALL CONVRT (K, CK, LCK)
          CALL OPENFL (K, 'csfg_'//FILEHEAD(1:LNG)//'.'//CK(1:2), &
                       'UNFORMATTED', 'UNKNOWN', IERR)
          IF (IERR == 0) CYCLE
          DO I = 30, K - 1
             CLOSE(I)
          END DO
          WRITE (ISTDE, *) 'Error when opening the mcp files'
          STOP
        END DO
!
! Write the heading lines of mcp.30
!
        WRITE (30) NCORE, NBLOCK, KMAX
        !if (myid == 0) WRITE (*, *)"setmcp_csfg: NCFBLK(I) =", &
        !                           (NCFBLK(I),I=1,NBLOCK)
        ! Number of CSFGs for every Jp block.
        WRITE (30) (NCFBLK(I),I=1,NBLOCK)

        ! Jp character
        WRITE (30) (IDBLK(I),I=1,NBLOCK)
        !if (myid == 0) WRITE (*, *)"setmcp_csfg: MCP,NBLOCK,MYID,NPROCS", &
        !               '  ', 'MCP', NBLOCK, MYID, NPROCS
        WRITE (30) 'MCP', MYID, NPROCS
        WRITE (31) 'MCP', MYID, NPROCS
        WRITE (32) 'MCP', MYID, NPROCS

      CASE (2)
! Write the heading lines of every block.
        !DO I = 1, 20
        !  K = FN(I)
        !  WRITE (2000+K, *) 'MCP', IBLK, NCF
        !ENDDO
      CASE (3)
! Write the end lines of every block. 
        !DO I = 1, 20
        !  K = FN(I)
        !  WRITE (2000+K, *) 0, 0, 0
        !ENDDO
        WRITE (31) 0, 0, 0
        WRITE (32) 0, 0, 0

      CASE (4)
! Close the files. 
        !DO I = 1, 20
        !  K = FN(I)
        !  CLOSE (2000+K)
        !ENDDO

        CLOSE (30)
        CLOSE (31)
        CLOSE (32)
      END SELECT

      RETURN
      END SUBROUTINE SETMCP_CSFG
