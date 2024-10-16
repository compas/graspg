!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**       ******    *****    *****   *******   *****    *****        ***
!**       **   **  **   **  **   **  **       **   **  **   **       ***
!**       **   **  **       **       **       **   **       **       ***
!**       ******    *****   **       ****      *****       **        ***
!**       **  **        **  **       **          **       **         ***
!**       **   **  **   **  **   **  **         **      **           ***
!**       **   **   *****    *****   **        **      *******       ***
!**                                                                  ***
!**            Relativistic Self-Consistent-Field Program            ***
!**                                                                  ***
!**   This program is a derivative of GRASP2 (F. A. Parpia, I. P.    ***
!**   Grant, and C. F. Fischer, 1990)                                ***
!**                                                                  ***
!**                            GRASP92                               ***
!**          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
!**   Modified by C.F. FIscher for block input                       ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
      PROGRAM RSCFmpiVU
!                                                                      *
!   Entry routine for RSCFVU. Controls the entire computation.      *
!                                                                      *
!   Call(s) to: [LIB92]: SETMC, SETCON.                                *
!               [RSCF92]: CHKPLT, setcsl, setdbg
!                        SETMIX, SETRES, SETSUM, STRSUM.
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 31 Dec 1992   *
!   Block version by Xinghong He          Last revision: 17 Aug 1998   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:57:54   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
! Replaced with memory-version codes (to read mcp files) provided by GG*
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!                                                                      *
!======================================================================*
!                                                                      *
! Modified by Chongyang Chen, Fudan University, Dec 2021               *
!                                                                      *
! Improve the way to calculate coefficeints used within SCF procedure  *
! Also, mpi is implemented for Matrix-Vector multipy and inner-dot     *
! calculation within dvdson subroutine.                                *
!                                                                      *
!======================================================================*
!                                                                      *
! Modification for using CSFG list, Chongyang Chen,                    *
!                        Fudan university, Shanghai,       Dec  2023   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
       USE default_C
       USE core_C
       USE iounit_C
       USE mpi_C
!GG po to isimti
       USE TATB_C

      Use fixd_C,        ONLY: LFIX
      Use csfg_scf_C,    ONLY: LSKIPPOT
      Use symexpand_mod
      Use csfg_decide_C, ONLY: LSYMD
      Use symmatrix_mod, ONLY: MAP, MAXSPAN, NCFBLKTr, NCFTOTTr, &
                               NTYPE, NORBGEN

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setdbgmpi_I
      USE setmc_I
      USE setcon_I
      USE setsum_I
      USE setcslmpi_I
      USE getscdmpi_I
      USE strsum_I
      USE setmix_I
      USE factt_I
      USE scfmpi_I

      USE setmcp_csfg_I
      USE rcsfsymexpand_I

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NBLK0 = 50  ! Maximum number of blocks
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NCORE1, NCOUNT1, lenperm, lentmp
      LOGICAL :: EOL, YES
      CHARACTER, DIMENSION(NBLK0) :: IDBLK*8
!      CHARACTER (LEN = *) :: host
      CHARACTER (LEN = 3) :: idstring
      CHARACTER (LEN = 128) :: startdir,permdir,tmpdir,FILE_RCSL,file1,file2
      CHARACTER (LEN = 128) :: FNOLDM, FILE_RL, FILE_OUT

      INTEGER :: NDUM
!-----------------------------------------------
!
      OPEN(UNIT=734,FILE='rmcdhf.log',STATUS='UNKNOWN')

!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      startdir = '  '  ;    file_rcsl = '  '
      permdir = '  '   ;    file1     = '  '
      tmpdir = '  '    ;    file2     = '  '
      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1,       &
                  startdir, permdir, tmpdir, 'RMCDHF_CSFG_MPI')
      WRITE (idstring, '(I3.3)') myid
      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)
!=======================================================================
!  Get NDEF
!=======================================================================

      IF (MYID == 0) THEN
         WRITE (istde,*)
         WRITE (ISTDE, '(A)') 'Default settings?  (y/n) '
         YES = GETYN()
         IF (YES) THEN
            NDEF = 0
!cjb fort.734
!cjb        WRITE(734,'(A)') 'y            ! Default settings'
!PERJ
            WRITE(734,'(A)') 'y            ! Default settings'
!PERJ END
         ELSE
            NDEF = 1
         ENDIF
      ENDIF
     CALL MPI_Bcast (NDEF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    CHKPLT - compatibility of plant substitutions
!    SETDBGmpi - Debug output control parameters
!    SETMC - machine- and installation- dependent constants
!    SETCON - physical constants
!    SETSUM - open the summary file
!    SETMCP - open and check the  .mcp  files
!    SETCSLmpi - open, check, load data from, and close the  .csl  file
!    STRSUM - append a summary of the inputs to the  .sum  file
!    SETMIX - mixing coefficients file setup
!    FACTT - table of logarithms of factorials setup
!=======================================================================
! 
! Read rcsfg.inp, expands it, broadcast the data.
!
      if (myid == 0) then
        ! use rcsfg.inp as input file
        FILE_RCSL = permdir(1:lenperm)//'/rcsfg.inp'
        FILE_RL   = permdir(1:lenperm)//'/rlabel.inp'
        FILE_OUT  = permdir(1:lenperm)//'/rcsf.out'
        call csfg_expand_findtype(FILE_RCSL, FILE_RL, FILE_OUT, 1)
      endif

!
! Broadcast the values from process 0
!
      CALL MPI_Bcast (NONSYM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NBLOCKSYM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NMAXGEN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NSYM_ORB(:,:),63,MPI_INTEGER,0,                  &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (TotCSFs_perblock(:,:),150,MPI_INTEGER,0,         &
                      MPI_COMM_WORLD,ierr)

      CALL MPI_Bcast (scfDF1,NBLOCKSYM,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (MAXSPAN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NORBGEN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nonsyminf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NCFGTOT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      NDUM = MAXVAL(TotCSFs_perblock(2,1:NBLOCKSYM))
!
! Allocate the arrays for myid .ne. 0
!
      if (myid.ne.0) then
        ALLOCATE(MAP1(NDUM,NBLOCKSYM))
        ALLOCATE(NTYPE(6,NCFGTOT))
      endif
      CALL MPI_Bcast (MAP1(:,:),  NDUM*NBLOCKSYM, MPI_INTEGER, 0,      &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NTYPE(:,:), 6*NCFGTOT,      MPI_INTEGER, 0,      &
                      MPI_COMM_WORLD,ierr)
      NCFBLKTr(1:NBLOCKSYM) = TotCSFs_perblock(3,1:NBLOCKSYM)
! Check
      IF (NCFGTOT /= SUM(TotCSFs_perblock(2,1:NBLOCKSYM))) THEN
        STOP "NCFGTOT /= SUM(TotCSFs_perblock(2,1:NBLOCKSYM)) ..."
      ENDIF
      NCFTOTTr = SUM(NCFBLKTr(1:NBLOCKSYM))
      IF (MYID == 0) WRITE (*,*)"NTOTCSF(G)s, NTOTCSFs =", NCFGTOT, NCFTOTTr
      IF (MYID == 0) WRITE (*,*)
!=======================================================================
!
      CALL SETDBGmpi (permdir(1:lenperm) // '/rscf92.dbg')
      CALL SETMC
      CALL SETCON

      IF (myid .EQ. 0) CALL SETSUM (permdir(1:lenperm) // '/rmcdhf.sum')

      !CALL SETMCP (NCORE, NBLK0, IDBLK, 'mcp' // idstring)
!
! Open files csfg_mcpXXX.30 (IGMCP30 = 3030) and 
!            csfg_mcpXXX.31 (IGMCP31 = 3031) (T-coefficients)
!            csfg_mcpXXX.32 (IGMCP31 = 3032) (V-coefficients)
      CALL SETMCP_CSFG (1, NCORE, NBLK0, IDBLK, 'mcp' // idstring)
!
! Now using normal CSF list, it will only use CSFG list when everything
! is OK.
!
      !if(myid == 0) file_rcsl = permdir(1:lenperm)//'/rcsf.inp'
      !if(myid == 0) file_rcsl = permdir(1:lenperm)//'/rcsfg.inp'
      CALL SETCSLmpi (file_rcsl, ncore1, idblk)
      IF (NCORE /= NCORE1) STOP 'rscfmpivu: ncore'
!
! Occupation numbers of the NCFTOTTr normal CSFs
      CALL CSFG_IQA 

!=======================================================================
!  Gather all remaining information and perform some setup. This
!  part (routine) asks for user-inputs.
!=======================================================================

      if(myid == 0) file1 = permdir(1:lenperm) // '/isodata'
      if(myid == 0) file2 = permdir(1:lenperm) // '/rwfn.inp'
      CALL GETSCDmpi (EOL, idblk, file1, file2 )
      file1 = '  '
      if(myid == 0) file1 = permdir(1:lenperm) // '/rmix.out'
!
! Use old.m to speed up the diagonalization, to be implemeneted
!
!      if(myid == 0) FNOLDM = permdir(1:lenperm) // '/old.m'
!
      IF (myid .EQ. 0) THEN
         CALL STRSUM
         IF (EOL) CALL SETMIX (file1)
      ENDIF

      CALL FACTT

      file1 = '  '
      if(myid == 0) file1 = permdir(1:lenperm) // '/rwfn.out'

!      IF (MYID == 0) CALL REPORT_MEM(1, 'RSCFmpiVU')

      CALL scfmpi (EOL, file1)

!
! Close files csfg_mcpXXX.30 (IGMCP30 = 3030) and 
!             csfg_mcpXXX.31 (IGMCP31 = 3031)
!             csfg_mcpXXX.32 (IGMCP31 = 3032)
      CALL SETMCP_CSFG (2, NCORE, NBLK0, IDBLK, 'mcp' // idstring)

      DEALLOCATE(MAP1)
      DEALLOCATE(NTYPE)

!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      CALL stopmpi2 (myid, nprocs, host, lenhost, ncount1, &
                     'RMCDHF_CSFG_MPI')

      END PROGRAM RSCFmpiVU
