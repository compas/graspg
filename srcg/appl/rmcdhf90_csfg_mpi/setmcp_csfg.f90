!***********************************************************************
!                                                                      *
      SUBROUTINE SETMCP_CSFG (IFLAG, NCORE, NBLKIN, IDBLK, FILEHEAD)
!                                                                      *
!======================================================================*
!                                                                      *
! IFLAG == 1:                                                          *
!                                                                      *
! Open and check the mcpXXX.30 which stores interacting CSFG pairs     *
! Open and check the mcpXXX.31 which stores spin-angular coefficients  *
! Open and check the mcpXXX.32 which stores spin-angular coefficients  *
!                                                                      *
!======================================================================*
! IFLAG == 2:                                                          *
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
!   Used by rscfvu90                                                   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE vast_kind_param, ONLY: LONG, DOUBLE
      USE csfg_memory_man
      Use def_C,     ONLY: NELEC
      USE mcp_C,     ONLY: KMAX, DIAG, LFORDR
      USE orb_C,     ONLY: NW, NKJ, NCF
      !USE orb_C
      USE hblock_C,  ONLY: NBLOCK, NCFBLK
      USE iounit_C,  ONLY: ISTDE
      USE iccu_C,    ONLY: ICCUT
!CYC      USE mpi_C,     ONLY: MYID, NPROCS, IERR
      USE mpi_C

      Use csfg_decide_C, ONLY: LSYMD 
      Use symmatrix_mod, ONLY: NCFBLKTr, NCFTOTTr
      Use csfg_tv_C
!      Use csfg_tv_C,     ONLY: IGMCP30, IGMCP31, NPAIRS, &
!                               NPASTGP, ICGPEND, IRGDET
      Use symexpand_mod, ONLY: TotCSFs_perblock, NCFGTOT
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I
      USE openfl_I

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)   :: IFLAG
      INTEGER               :: NCORE
      INTEGER, INTENT(IN)   :: NBLKIN
      CHARACTER             :: IDBLK(*)*8 
      CHARACTER, INTENT(IN) :: FILEHEAD*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER       :: K, LNG, LCK, I, IBLK, NB
      LOGICAL       :: FOUND, FOUND1, GETYN, YES
      CHARACTER     :: CK*2, MCPLAB*3

      CHARACTER     :: IDBLKmcp(50)*8
      INTEGER       :: NCFBLKmcp(0:50)
      INTEGER       :: NCOREmcp, NBLOCKmcp, KMAXmcp 
      INTEGER       :: NELECmcp, NCFmcp, NWmcp, MYIDmcp, NPROCSmcp
      INTEGER       :: IERROR, IOS, NCFNB, NINTACT, NINTACTT
      INTEGER       :: NDIM, NCFTrmcp, NDIMNZ
!-----------------------------------------------

      SELECT CASE (IFLAG)
      CASE (1)
!
! Section 1: 
! Open csfg_mcpXXX.30 and csfg_mcpXXX.31
!
        LNG = LEN_TRIM(FILEHEAD)
        ! DO K = 30, 32 + KMAX
        DO K = 30, 32
          CALL CONVRT (K, CK, LCK)
          CALL OPENFL (3000+K, 'csfg_'//FILEHEAD(1:LNG)//'.'//CK(1:2), &
                       'UNFORMATTED', 'OLD', IERR)
          IF (IERR == 0) CYCLE
          DO I = 30, K - 1
             CLOSE(I)
          END DO
          WRITE (ISTDE, *) 'Error when opening the csfg_mcp files'
          STOP
        END DO
        IGMCP30 = 3030
        IGMCP31 = 3031
        IGMCP32 = 3032
        IERROR  = 0 
!
! Section 2: 
! Read the heading lines of csfg_mcpXXX.30, which are written in
! setmcp_csfg.f90
!
        NCFBLKmcp(0) = 0
        READ (IGMCP30, IOSTAT=IOS) NCOREmcp, NBLOCKmcp, KMAXmcp
        IERROR = IERROR + ABS(IOS)
!
        READ (IGMCP30, IOSTAT=IOS) (NCFBLKmcp(I), I = 1, NBLOCKmcp)
        IERROR = IERROR + ABS(IOS)

        NCORE  = NCOREmcp
        NBLOCK = NBLOCKmcp
        CALL ALLOC (NCFBLK, 0, NBLOCK , 'NCFBLK', 'SETMCP')
        NCFBLK(1:NBLOCK) = NCFBLKmcp(1:NBLOCK)
!
        DO I = 1, NBLOCKmcp
          IF (NCFBLKmcp(I) /= TotCSFs_perblock(2, I)) THEN
            WRITE(ISTDE, *)"            NCFBLKmcp =", &
                           NCFBLKmcp(1:NBLOCKmcp)
            WRITE(ISTDE, *)"TotCSFs_perblock(2,I) =", &
                           TotCSFs_perblock(2,1:NBLOCKmcp)
            STOP "SETMCP_CSFG: NCFBLKmcp(I) /= TotCSFs_perblock(2, I)!"
          ENDIF
        ENDDO

        READ (IGMCP30, IOSTAT=IOS) (IDBLKmcp(I),  I = 1, NBLOCKmcp)
        IERROR = IERROR + ABS(IOS)
!
        IDBLK(1:NBLOCKmcp) = IDBLKmcp(1:NBLOCKmcp)

        READ (IGMCP30, IOSTAT=IOS) MCPLAB, MYIDmcp, NPROCSmcp
        IERROR = IERROR + ABS(IOS)

        ! Check
        IF (MYID/=MYIDmcp .OR. NPROCS/=NPROCSmcp) THEN
          WRITE (ISTDE, *) 'mcp files were generated under different', &
                           ' processor configuration.'
          STOP
        ENDIF
        IF (MCPLAB /= 'MCP') THEN
          WRITE (ISTDE, *) 'Not a GRASP92 csfg_MCPXXX.30 File;'
          IERROR = IERROR + 1
        ENDIF
        IF (IERROR /= 0) THEN
          WRITE (ISTDE, *) 'setmcp: Error accumulated , stopping...'
          STOP  "SETMCP_CSFG ===01==="
        ENDIF
!
! Section 3:
! Read the lines written in setmcpmpi.f90
! NCFmcp = Total number of CSFGs in rcsfg.inp, all blocks.
        READ (IGMCP30, IOSTAT=IOS) NELECmcp, NCFmcp, NWmcp
        NELEC = NELECmcp
        NCF = NCFmcp
        NW  = NWmcp

        IERROR = IERROR + ABS(IOS)
        NCFGPAST = SUM(NCFBLKmcp(1:NBLOCKmcp))
        IF (NCFmcp /= NCFGTOT .OR. NCFGPAST /= NCFmcp) THEN
          WRITE(ISTDE, *)"NCFmcp, NCFGTOT, NCFGPAST(NBLOCK) =", &
                         NCFmcp, NCFGTOT, NCFGPAST
          STOP &
             "ERROR! NCFmcp /= NCFGTOT .OR. NCFGPAST /= NCFmcp ..."
        ENDIF

        READ (IGMCP30, IOSTAT=IOS) DIAG, LFORDR, LSYMD
        IERROR = IERROR + ABS(IOS)
!        WRITE(*,*)'MYID, DIAG, LFORDR, LSYMD =', &
!                   MYID, DIAG, LFORDR, LSYMD
        READ (IGMCP30, IOSTAT=IOS) ICCUT(1:NBLOCKmcp)
        IERROR = IERROR + ABS(IOS)
        ! Check
        IF (IERROR /= 0) THEN
          WRITE (ISTDE, *) 'setmcp: Error accumulated , stopping...'
          STOP  "SETMCP_CSFG ===02==="
        ENDIF
!
! Section 4:
! NBLOCK -- Loop, obtain the interacting CSFG pairs, the non-zero
! matrixelement positions, which are written in setsda_csfg.f90
!
! Read ICGSCF, IRGSCF
!
! Column index in CSFG list for interacting pairs
        CALL ALLOC(ICGPDET, NCFGTOT, 'ICGPDET', 'SETMCP_CSFG') 
        CALL ALLOC(ICGPEND, NCFGTOT, 'ICGPEND', 'SETMCP_CSFG') 
!
! Column index in normal CSF list for non-zero matrixelements
!        CALL ALLOC(IENDCTV, NCFTOTTr, 'IENDCTV', 'SETMCP_CSFG') 
!        CALL ALLOC(NODECOLSTV,NCFTOTTr, 'NODECOLSTV', 'SETMCP_CSFG') 
!
! Row index in CSFG list for interacting pairs
! needs re-alloc
        NDIM = NCFGTOT  
        CALL ALLOC(IRGPDET,  NDIM, 'IRGDET', 'SETMCP_CSFG')
!
! Row index in normal CSF list for non-zero matrixelements
! needs re-alloc
        NDIMNZ = NCFTOTTr * 2
        !CALL ALLOC(IROWTV, NDIMNZ, 'IROWTV', 'SETMCP_CSFG') 

        !NCFGPAST, NPAIRPAST, NZPAST, NCOLPAST
        NPAIRS(1:2,0) = 0; NCOLBLK(0) = 0; NZBLK(0) = 0

        DO IBLK = 1, NBLOCK
! WRITE (30) 'MCP', NB, NCFG, NCFTr
          READ (IGMCP30, IOSTAT=IOS)MCPLAB, NB, NCFNB, NCFTrmcp
          IERROR = IERROR + ABS(IOS)
          IF (NB /= IBLK .OR. NCFNB /= NCFBLKmcp(NB) .OR. &
              NCFTrmcp /= TotCSFs_perblock(3, NB)) THEN
            WRITE(ISTDE, *)"IBLK, NB, NCFNB, NCFBLKmcp(NB), NCFTrmcp, &
                            TotCSFs_perblock(3, NB) =", &  
                            IBLK, NB, NCFNB, NCFBLKmcp(NB), NCFTrmcp, &
                            TotCSFs_perblock(3, NB)
            STOP "ERROR! NB /= IBLK .OR. NCFNB /= NCFBLKmcp(NB) ..."
          ENDIF
         
!      WRITE (30) NPAIRS(1,NB), NPAIRS(2,NB)
          READ(IGMCP30, IOSTAT=IOS)NPAIRS(1:2,NB)
          IERROR = IERROR + ABS(IOS)

!      WRITE (30) ICGPDET(1:NPAIRS(1,NB)) !
          NPAIRPAST(1) = SUM(NPAIRS(1, 0:NB-1))
          READ(IGMCP30, IOSTAT=IOS)                                    &
                       (ICGPDET(NPAIRPAST(1) + I), I = 1, NPAIRS(1, NB))
          IERROR = IERROR + ABS(IOS)

!      WRITE (30) ICGPEND(1:NPAIRS(1,NB)) !
          READ(IGMCP30, IOSTAT=IOS)                                    &
                       (ICGPEND(NPAIRPAST(1) + I), I = 1, NPAIRS(1, NB))
          IERROR = IERROR + ABS(IOS)

!      WRITE (30) IRGPDET(1:NPAIRS(2,NB))
          NPAIRPAST(2) = SUM(NPAIRS(2, 0:NB-1))
          IF (NPAIRPAST(2) + NPAIRS(2, NB) > NDIM) THEN
            NDIM = NPAIRPAST(2) + NPAIRS(2, NB)
            CALL RALLOC (IRGPDET, NDIM, 'IRGPDET', 'SETMCP_CSFG')
          ENDIF
          READ(IGMCP30, IOSTAT=IOS)                                    &
              (IRGPDET(NPAIRPAST(2) + I), I = 1, NPAIRS(2, NB))

!      WRITE (30) NELMNTCC, NODENCOLS, ncsfDF1, ncsfDF1_core
          READ(IGMCP30, IOSTAT=IOS) NZBLK(NB),  NCOLBLK(NB),   &
                                    MCPDF1(NB), MCPDF1C(NB)
          IERROR = IERROR + ABS(IOS)
!
! Reproduce in setham_pot
! 
!!!      WRITE (30) (NODECOLS(I), I = 1, NODENCOLS) !
!!          NCOLPAST = SUM(NCOLBLK(0:NB-1))
!!          READ(IGMCP30, IOSTAT=IOS)                                    &
!!                        (NODECOLSTV(NCOLPAST + I), I = 1, NCOLBLK(NB))
!!          IERROR = IERROR + ABS(IOS)
!!!      WRITE (30) (IENDC(I),    I = 1, NODENCOLS) !
!!          READ(IGMCP30, IOSTAT=IOS)                                    &
!!                        (IENDCTV (NCOLPAST + I), I = 1, NCOLBLK(NB))
!!          IERROR = IERROR + ABS(IOS)
!!
!!!      WRITE (30) (IROW(I),     I = 1, NELMNTCC)  !
!!          NZPAST = SUM(NZBLK(0:NB-1))
!!          IF (NZPAST + NZBLK(NB) > NDIMNZ) THEN
!!            NDIMNZ = NZPAST + NZBLK(NB)
!!            CALL RALLOC(IROWTV, NDIMNZ, 'IROWTV', 'SETMCP_CSFG') 
!!          ENDIF
!!          READ(IGMCP30, IOSTAT=IOS)                                    &
!!                        (IROWTV(NZPAST + I), I = 1, NZBLK(NB))
!!          IERROR = IERROR + ABS(IOS)
        ENDDO  
        ! Check
        IF (IERROR /= 0) THEN
          WRITE (ISTDE, *) 'setmcp: Error accumulated , stopping...'
          STOP  "SETMCP_CSFG ===03==="
        ENDIF
!
! Condense the array sizes
        ! Column index in CSFG list
        NDIM = SUM (NPAIRS(1, 0:NBLOCK))
        CALL RALLOC(ICGPDET, NDIM, 'ICGPEND', 'SETMCP_CSFG') 
        CALL RALLOC(ICGPEND, NDIM, 'ICGPEND', 'SETMCP_CSFG') 

        ! Column index in normal CSF list
        NDIM = SUM (NCOLBLK(0:NBLOCK))
!        CALL RALLOC(NODECOLSTV,NDIM, 'NODECOLSTV', 'SETMCP_CSFG') 
!        CALL RALLOC(IENDCTV, NDIM, 'IENDCTV',  'SETMCP_CSFG') 

        IF (MYID == 0) WRITE(*,*)"csfg_mcpXXX.30 loaded ...."
        NINTACT = SUM(NPAIRS(2,1:NBLOCK))
        CALL MPI_Reduce (NINTACT, NINTACTT, 1, MPI_INTEGER, MPI_SUM, 0, &
                 MPI_COMM_WORLD, IERR)
        IF (MYID == 0) &
          WRITE(*,*)"Total Interacting ICG - IRG CSFG Pairs :", NINTACTT
        IF (MYID == 0) WRITE(*,*)

!
! Section 5: 
! Read spin-angular coefficients from csfg_mcpXXX.31
!
        READ (IGMCP31) MCPLAB, MYIDmcp, NPROCSmcp
        ! Check
        IF (MYID/=MYIDmcp .OR. NPROCS/=NPROCSmcp) THEN
          WRITE (ISTDE, *) 'mcp files were generated under different', &
                           ' processor configuration.'
          STOP
        ENDIF
        IF (MCPLAB /= 'MCP') THEN
          WRITE (ISTDE, *) 'Not a GRASP92 csfg_MCPXXX.31 File;'
          IERROR = IERROR + 1
        ENDIF
        IF (IERROR /= 0) THEN
          WRITE (ISTDE, *) 'setmcp: Error accumulated , stopping...'
          STOP  "SETMCP_CSFG ===04==="
        ENDIF
!
! Block loop is located in setham_cyc.f90
!
! Section 6: 
! Read spin-angular coefficients from csfg_mcpXXX.32
        READ (IGMCP32) MCPLAB, MYIDmcp, NPROCSmcp
        ! Check
        IF (MYID/=MYIDmcp .OR. NPROCS/=NPROCSmcp) THEN
          WRITE (ISTDE, *) 'mcp files were generated under different', &
                           ' processor configuration.'
          STOP
        ENDIF
        IF (MCPLAB /= 'MCP') THEN
          WRITE (ISTDE, *) 'Not a GRASP92 csfg_MCPXXX.32 File;'
          IERROR = IERROR + 1
        ENDIF
        IF (IERROR /= 0) THEN
          WRITE (ISTDE, *) 'setmcp: Error accumulated , stopping...'
          STOP  "SETMCP_CSFG ===05==="
        ENDIF
!
! Block loop is located in setham_cyc.f90
!
      CASE (2)
! Close the files. 
        ! close mcpXXX_csfg.30
        ! close mcpXXX_csfg.31
        CLOSE (IGMCP30)
        CLOSE (IGMCP31)
! Deallocate arrays 
        CALL DALLOC(ICGPDET, 'ICGPDET', 'SETMCP_CSFG') 
        CALL DALLOC(ICGPEND, 'ICGPEND', 'SETMCP_CSFG') 
!        CALL DALLOC(IENDCTV, 'IENDCTV', 'SETMCP_CSFG') 
!        CALL DALLOC(NODECOLSTV,'NODECOLSTV', 'SETMCP_CSFG') 
        CALL DALLOC(IRGPDET, 'IRGDET', 'SETMCP_CSFG')
        !CALL DALLOC(IROWTV,  'IROWTV', 'SETMCP_CSFG') 

      END SELECT

      RETURN
      END SUBROUTINE SETMCP_CSFG
