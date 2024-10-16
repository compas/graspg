!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**        *****   ******  **   **  **    **   *****   ******        ***
!**       **   **  **      ***  **  ***  ***  **   **  **   **       ***
!**       **       **      ***  **  ** ** **  **       **   **       ***
!**       **  ***  ****    ** ****  ** ** **  **       ******        ***
!**       **   **  **      **  ***  **    **  **       **            ***
!**       **   **  **      **   **  **    **  **   **  **            ***
!**        *****   ******  **   **  **    **   *****   **            ***
!**                                                                  ***
!**   Program for generating the energy expression for H(DC). This   ***
!**   program is a derivative of GRASP2 (F. A. Parpia, I. P. Grant,  ***
!**   and C. F. Fischer, 1990).                                      ***
!**                                                                  ***
!**                         GRASP92 Version                          ***
!**          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
      PROGRAM GENMCPMPI
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 11 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 29 Jun 1998   *
!   Updated by Charlotte F. Fischer
!                                                                      *
!   Modified by Gediminas Gaigalas for new spin-angular integration.   *
!                                                                      *
!***********************************************************************
!                                                                      *
!   Modified for CSFG list                                             *
!   Chongyang Chen, Fudan University, Shanghai           28 Dec 2020   *
!                                                                      *
!   Last modifications by Chongyang Chen, Aug 2024                     *
!                                                                      *
!   Input  files: rcsfg.inp, rlabel.inp, isodata                       *
!   Output files: rangular.log                                         *
!                 ${MPI_TMP}/idx/csfg_mcpXXX.30, interating CSFG-pairs *
!                 ${MPI_TMP}/idx/csfg_mcpXXX.31, T-coefficients        *
!                 ${MPI_TMP}/idx/csfg_mcpXXX.32, V-coefficients        *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
      USE parameter_def,   ONLY:  NNNW
      USE csfg_memory_man
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE mpi_C
      USE default_C
      Use hblock_c
      Use orb_C
      Use stat_c
      Use mcp_c
      Use iounit_c
!CYC
      Use symexpand_mod
      Use csfg_decide_C, ONLY: LSYMD
      Use symmatrix_mod, ONLY: MAP,NCFBLKTr,NCFTOTTr,NCFTr
      Use csfg_tv_C,     ONLY: ICGPEND
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setdbg_I
      USE setmc_I
      USE setsum_I
      USE cslh_I
      USE strsum_I
      USE factt_I
      USE settmp_I
      USE lodcslmpi_csfg_I

      USE setmcp_csfg_I
      USE rcsfsymexpand_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER NBLK0, NCOUNT1, NCORE, NB
      PARAMETER (NBLK0 = 50)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL ::  DEBUG, RESTRT, YES
      CHARACTER(LEN=128) :: STARTDIR, PERMDIR, TMPDIR
      CHARACTER(LEN=128) :: FILE_RCSL, FILE_RL, FILE_OUT
      CHARACTER(LEN=8), DIMENSION(NBLK0):: IDBLK
      CHARACTER(LEN=3) :: IDSTRING
      INTEGER :: LENPERM, LENTMP
      INTEGER :: I, J, ICSF, NDUM
!-----------------------------------------------
      OPEN(UNIT=739,FILE='rangular.log',STATUS='UNKNOWN')
      startdir = ' '
      permdir = ' '
      tmpdir = ' '

!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1, &
                      startdir, permdir, tmpdir, 'RANGULAR_CSFG_MPI')
      WRITE (idstring, '(I3.3)') myid
      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)
!=======================================================================
!  Get NDEF on node-0 and then send to all nodes
!=======================================================================
      IF (myid == 0) THEN
         WRITE (istde,*)
!         WRITE (istde,'(A)') ' Default settings?  (y/n) '
         WRITE (istde,'(A)') ' Full interaction?  (y/n) '
         YES = GETYN ()
         IF (YES) THEN
            NDEF = 0
            write(739,'(A)') 'y            ! Full interaction'
         ELSE
            NDEF = 1
            write(739,'(A)') 'n            ! Full interaction'
         ENDIF
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NDEF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    setdbg - debug output control parameters
!    setmc - machine- and installation- dependent constants
!    setsum - open the summary file
!    cslh - load header of the csl file
!    setmcp - open and check the  .mcp  files
!    strsum - append a summary of the inputs to the  .sum  file
!    factt - table of logarithms of factorials setup
!=======================================================================

      CALL setdbgmpi (DEBUG, permdir(1:lenperm) // '/genmcp.dbg')
      CALL SETMC
      if (myid == 0) then
        !file_rcsl = permdir(1:lenperm)//'/rcsf.inp'
        FILE_RCSL = permdir(1:lenperm)//'/rcsfg.inp'
        FILE_RL   = permdir(1:lenperm)//'/rlabel.inp'
        FILE_OUT  = permdir(1:lenperm)//'/rcsf.out'
        call rcsfsymexpand(FILE_RCSL, FILE_RL, FILE_OUT, 0)
        DEALLOCATE(MAP2)
      endif

!=======================================================================
! Broadcast the values from process 0
!
      CALL MPI_Bcast (NONSYM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NBLOCKSYM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NMAXGEN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NSYM_ORB(:,:),63,MPI_INTEGER,0,                  &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (TotCSFs_perblock(:,:),150,MPI_INTEGER,0,         &
                      MPI_COMM_WORLD,ierr)
      NDUM = MAXVAL(TotCSFs_perblock(2,1:NBLOCKSYM))
! Allocate the arrays for myid .ne. 0
      if (myid.ne.0) then
        ALLOCATE(MAP1(NDUM,NBLOCKSYM))
        !ALLOCATE(MAP2(4,NDUM,NBLOCKSYM))
      endif
      CALL MPI_Bcast (MAP1(:,:),NDUM*NBLOCKSYM,MPI_INTEGER,0,          &
                      MPI_COMM_WORLD,ierr)
      !CALL MPI_Bcast (MAP2(:,:,:),4*NDUM*NBLOCKSYM,MPI_INTEGER,0,      &
      !                MPI_COMM_WORLD,ierr)
      NCFBLKTr(1:NBLOCKSYM) = TotCSFs_perblock(3,1:NBLOCKSYM)
      NCFTOTTr = SUM(NCFBLKTr(1:NBLOCKSYM))
!=======================================================================
!      file_rcsl = trim(file_rcsl)//'.inp'
      CALL cslhmpi_csfg (FILE_RCSL, ncore, nblk0, idblk)
      RESTRT = .FALSE.
!cjb  myid, nprocs = NOT args
!cjb  CALL setmcpmpi (myid, nprocs, ncore, idblk, 'mcp' // idstring)

!CYC:
! LSYMD, ICCUT, LFORDR
! Open file:
! csfg_mcpXXX.30 stores the interacting CSFG pairs.
! csfg_mcpXXX.31 and csfg_mcpXXX.32 store the T- and V-coefficients of 
! the above interacting CSFG pairs, respectively. 
      CALL setmcpmpi (ncore, idblk, 'mcp' // idstring)
!
!      IF (NDEF/=0 .AND. MYID==0) THEN
!        OPEN(24, file=permdir(1:lenperm)//'/'//'icut.log') 
!        CALL STRSUM
!      ENDIF
      CALL FACTT

      CALL MPI_Bcast (LSYMD, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

!=======================================================================
!     For each block, generate and sort the data
!=======================================================================
      DO NB = 1, NBLOCK
         NCF = NCFBLK(NB)                        ! This ncf goes to common
         NCFTr =  NCFBLKTr(NB)
         IF (MYID == 0) THEN
            WRITE (6, *)
            WRITE (6, *) 'Block ', NB, ',  NCSF(G)s = ', NCF
         ENDIF
         ! Check, NCF is now COUNT1 in rcsfsymexpand.f90
         IF (NCF .NE. TotCSFs_perblock(2,NB)) then
           WRITE(*,*)'NCF,NB,TotCSFs_perblock(2,NB)=',                  &
               NCF,NB,TotCSFs_perblock(2,NB)
           STOP 'Unexpected NCF .NE. TotCSFs_perblock(2,NB) ...'
         ENDIF
! Record the end calls of ONESCALAR and RKCO_GG for every
! column JA (IC).
         ALLOCATE(ICGPEND(0:NCF))
         ICGPEND(0) = 0

         !*** Write the head lines 
         CALL SETMCP_CSFG (2, ncore, idblk, NB, '  ')
!=======================================================================
! Intialize the MAP(:) in symmatrix_mod
         ALLOCATE(MAP(NCF))
         DO I = 1, TotCSFs_perblock(2,NB)
           MAP(I) = MAP1(I,NB)
         ENDDO
!=======================================================================
         !*** Load current CSL block. Memories de-allocated in mcp ***
         ! Needs fictious CSF
! CYC: Allocate the arrays of NCF+10, NCF is COUNT1 CSFGs stored in
! rcsfg.inp, more 10 elements are for FICTIOUS CSFGs.
         CALL ALLOC (iqa, NNNW, NCF+10, 'IQA', 'GENMCP')
         CALL ALLOC (jqsa, NNNW, 3, NCF+10, 'JQSA', 'GENMCP')
         CALL ALLOC (jcupa, NNNW, NCF+10, 'JCUPA', 'GENMCP')
!
         CALL LODCSLmpi_csfg (21, NCORE, NB)

!         !*** Open tmp.xx files for block NB ***
! Not used, unsorted mode used.
!         CALL SETTMP (NB, KMAX, 'tmp' // idstring)

         !*** Generation of MCP coefficients ***
         CALL MCPmpi (NB, RESTRT, 'mcp' // idstring)

         !*** Write the end lines for each block
         CALL SETMCP_CSFG (3, ncore, idblk, NB, '  ')

!         IF (MYID == 0) THEN
!            WRITE (6, *)
!            WRITE (6, *) 'Finished NB ... =', NB
!         ENDIF
! Deallocate MAP and ISPENDC for next block
         DEALLOCATE(MAP)
         DEALLOCATE(ICGPEND)
      ENDDO

      CALL SETMCP_CSFG (4, ncore, idblk, NB, '  ')

!      CLOSE(24)                                  ! Summary file
      CLOSE(739)                                 ! rangular.log
      IF (DEBUG) CLOSE(99)                       ! Debug file
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================
      CALL stopmpi2 (myid,nprocs,host,lenhost,ncount1, &
                     'RANGULAR_CSFG_MPI')
      STOP
      END PROGRAM GENMCPMPI
