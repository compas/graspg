!***********************************************************************
!                                                                      *
      SUBROUTINE SCFmpi(EOL, RWFFILE2)
!                                                                      *
!   This  subroutine  performs  the SCF iterations. The procedure is   *
!   essentially algorithm 5.1 of C Froese Fischer, Comput Phys Rep 3   *
!   (1986) 290.                                                        *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
!               [RSCF92]: improvmpi, matrixmpi, MAXARR, newcompi,      *
!                         ORBOUT, ORTHSC, setlagMpi.                   *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 05 Aug 1998   *
!   ifort -i8 version by Alexander Kramida (AK) Last rev. 22 Mar 2016  *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!                                                                      *
!   Modified by Chongyang Chen,  Dec 2021                              *
!      Improvements                                                    *
!                  1: Ways to caclulate potentials                     *
!                  2: Calculation of Lagrange multiplier               *
!                  3: Diagonalization H matrix, MPI DVDSON             *
!                                                                      *
!   Modifcations for using CSFG list, Chongyang Chen,       Dec 2023   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:16:00   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE, LONG
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
     
      Use wave_C
      USE symmatrix_mod,  ONLY: NCFBLKTr,  NCFTOTTr, MAXSPAN, NODENCOLS
      Use symexpand_mod,  ONLY: NCFGTOT, TotCSFs_perblock
      Use csfg_tv_C  
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE newcompi_I
      USE setlagmpi_I
      USE improvmpi_I
      USE maxarr_I
      USE prwf_I
      USE orthsc_I
      USE orbout_I
      USE endsum_I
      USE setham_pot_I
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL  :: EOL
      CHARACTER  :: RWFFILE2*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, I, NIT, JSEQ, KOUNT, K, L_MPI, JBLOCK
      REAL(DOUBLE) :: WTAEV, WTAEV0, DAMPMX
      LOGICAL :: CONVG, LSORT, DVDFIRST, LCOFABV
      INTEGER :: NDDIM0, NXDIM0, NYDIM0, NCSFTMP, J2MAX
      INTEGER :: NCOUNTV0, NDIMV0, NINTACT, NINTACTT
      INTEGER :: ncount0, ncount1, ncount2, ncount3, ncount_rate, ncount_max
      REAL    :: nsecs, starttime, endtime
      INTEGER(LONG) :: NSPANVK

!-----------------------------------------------------------------------
      L_MPI = MPI_LOGICAL

!=======================================================================
!   Determine Orthonomalization order --- lsort
!=======================================================================
      IF (NDEF == 0) THEN
         LSORT = .FALSE.
      ELSE
         IF (myid .EQ. 0) THEN
  123       CONTINUE
            WRITE (ISTDE, *) 'Orthonomalization order? '
            WRITE (ISTDE, *) '     1 -- Update order'
            WRITE (ISTDE, *) '     2 -- Self consistency connected'
            READ (ISTDI, *) J
            IF (J == 1) THEN
               LSORT = .FALSE.
            ELSE IF (J == 2) THEN
               LSORT = .TRUE.
            ELSE
               WRITE (ISTDE, *) 'Input is wrong, redo...'
               GO TO 123
            ENDIF
         ENDIF
         CALL MPI_Bcast (lsort, 1, L_MPI, 0, MPI_COMM_WORLD, ierr)
      ENDIF

!=======================================================================
!   Deallocate storage that will no longer be used
!=======================================================================

!GG      CALL DALLOC (JQSA, 'JQSA', 'SCFmpi')

!=======================================================================
!   Allocate and fill in auxiliary arrays
!=======================================================================

      CALL ALLOC (NCFPAST, NBLOCK, 'NCFPAST', 'SCFmpi')
      CALL ALLOC (NCMINPAST, NBLOCK, 'NCMINPAST', 'SCFmpi')
      CALL ALLOC (NEVECPAST, NBLOCK, 'NEVECPAST', 'SCFmpi')
      CALL ALLOC (EAVBLK, NBLOCK, 'EAVBLK', 'SCFmpi')

      NCFPAST(1) = 0
      NCMINPAST(1) = 0
      NEVECPAST(1) = 0
      DO I = 2, NBLOCK
         NCFPAST(I) = NCFPAST(I-1) + NCFBLKTr(I - 1)
         NCMINPAST(I) = NCMINPAST(I-1) + NEVBLK(I - 1)
         NEVECPAST(I) = NEVECPAST(I-1) + NEVBLK(I - 1)*NCFBLKTr(I - 1)
      END DO

      !*** Size of the eigenvector array for all blocks
      NVECSIZ = NEVECPAST(NBLOCK) + NEVBLK(NBLOCK)*NCFBLKTr(NBLOCK)

      !*** Total number of the CSFs
      NCSFTMP = NCFPAST(NBLOCK) + NCFBLKTr(NBLOCK)
      IF (NCFTOTTr /= NCSFTMP) THEN
        write(*,*) NCFPAST(1:NBLOCK)
        write(*,*)NCFBLKTr(1:NBLOCK)
        write(*,*)'NCFGTOT, NCSFTMP = ', NCFGTOT, NCSFTMP
        STOP 'NCFGTOTT /= NCSFTMP ...' 
      ENDIF
!
      IF (EOL) THEN
         CALL ALLOC (EVAL,   NCMIN,   'EVAL', 'SCFmpi')
         CALL ALLOC (EVEC,   NVECSIZ, 'EVEC', 'SCFmpi')
         CALL ALLOC (IATJPO, NCMIN, 'IATJPO', 'SCFmpi')
         CALL ALLOC (IASPAR, NCMIN, 'IASPAR', 'SCFmpi')
      ENDIF
!=======================================================================
! Initialize the array sizes
      NDDIM = 32          ! Enough for n<=25 calculation
      NXDIM = 800000      ! Enough for 13 [s - m] calculation
      NYDIM = 7000        ! Enough for 13 [s - m] calculation
      NDIMT = 256*256
      NDIMV = 2048*2048   ! Enough for 13 [s - m] calculation?
      NDIMX = 1024

!=======================================================================
! CYC: Allocate the arrays to save T- and V- coefficients read from 
! csfg_mcpXXX.31 and csfg_mcpXXX.32, variables defined in csfg_tv_C.f90
! These data are used to construct the H-matrix, and the potentials.
      NDHGPT = 256*256
      CALL ALLOC (LABTH,  2, NDHGPT, 'LABTH',  'SCFmpi')
      CALL ALLOC (TCOEFH,    NDHGPT, 'TCOEFH', 'SCFmpi')

      NDHGPV = 256*256
      CALL ALLOC (NVHP,      NDHGPV, 'NVHP',   'SCFmpi')
      NDHVK = 2048*2048
      CALL ALLOC (LABVKH, 5, NDHVK,  'LABVKH', 'SCFmpi')
      CALL ALLOC (VCOEFH,    NDHVK,  'VCOEFH', 'SCFmpi')
!======================================================================
! CYC: set arrays to speed up the potential calculations
! see the details presented in Atoms 11, (2023) 12.
      CALL ALLOC (NXAIORB, NDIMX, NW, NW, 'NXAIORB','SCFmpi')
      ! JXIPOS(0:NW,NW) = 0 ! Initialized within csfg_scf_C.f90
      CALL ALLOC (NDCOFopt, NW, 'NDCOFopt', 'SCFmpi')
      CALL ALLOC (NXCOFopt, NW, 'NXCOFopt', 'SCFmpi')
      CALL ALLOC (NYCOFopt, NW, 'NYCOFopt', 'SCFmpi')
!CYC: NXAopt and XAopt could be vectorized (one-demension) to save 
!     memory, multi-demension arrays need more memory. 
!     Vectorized in Aug 2024
!CYC: NYAopt and YAopt is better to be also vectorized, though
!     their demensions are not large.
!     Vectorized in Aug 2024
      CALL ALLOC (NDAopt, NDDIM, NW, 'NDAopt','SCFmpi')
      CALL ALLOC (NXAopt, NXDIM, NW, 'NXAopt','SCFmpi')
      CALL ALLOC (NYAopt, NYDIM, NW, 'NYAopt','SCFmpi')
      CALL ALLOC (DAopt,  NDDIM, NW, 'DAopt', 'SCFmpi')
      CALL ALLOC (XAopt,  NXDIM, NW, 'XAopt', 'SCFmpi')
      CALL ALLOC (YAopt,  NYDIM, NW, 'YAopt', 'SCFmpi')
!=======================================================================
!CYC: Modified from GRASP2018::RCI_MPI module, speeding up SETHAM_POT. 
!   Calculate all the needed Rk integrals - genintrkwrap calls genintrk
!   to compute the integrals (each processor allocates the same amount
!   of memory, but only computes part of the integrals); and then do the
!   "MPI_Allreduce" stuff. The original codes are slightly modified.
      starttime = MPI_Wtime()
      CALL GENINTRKwrap_CSFG (myid, nprocs, j2max)
      !CALL GENINTRKwrap (myid, nprocs, j2max)
      endtime = MPI_Wtime()
      if (myid.eq.0)                                                &
        print *, 'Time for intrk calculation (s) = ',               &
                  endtime-starttime
!=======================================================================
! In the first calling MATRIX_CSFG, contstruct H matrix and diagonalize
! them, recored the spin-angular coefficients to construct H matrix.

! Initialization to record T and V data for H-matrix construction
! read from csfg_mcpXXX.31, 32
      NUMHT = 0;  NUMHP = 0;  NUMHVK = 0

! Columns index involved for MYID process, used in NEWCOMPI.f90
      CALL ALLOC (ICOLIDX, SUM(NCOLBLK(1:NBLOCK)), 'ICOLIDX', 'SCFmpi')
!
!CYC: Call MATRIXmpi and NEWCOmpi either EOL is .TRUE. or .FALSE.
!     For (E)OL calculations, determine the level energies and
!     mixing coefficients

!CFF   .. set the logical variable DVDFIRST
! Using the old mixing coefficients?
      DVDFIRST = .TRUE.
! CYC:
! Read spin-angular data from mcp files or from arrays?
      LREADMCP = .TRUE.
! Calculate H-Matrix and then diagonalize it?
      LCHM = .TRUE.
! Calculate potentials?
      LCPOT = .FALSE.
! Build the arrays involving potential calculations?
      LABTVFRST = .FALSE.

      starttime = MPI_Wtime()
      CALL MATRIX_CSFG(EOL, DVDFIRST)
      endtime = MPI_Wtime()
      if (myid.eq.0)                                                &
        print *, 'Time for the first call MATRIX_CSFG (s) = ',      &
                  endtime-starttime

! Spin-angular data to construct H-matrix using CSFG list
! MPI_Reduce to report the total T, and V-coefficients
      CALL MPI_Reduce (NUMHT, NDHGPT, 1, MPI_INTEGER, MPI_SUM, 0,   &
                       MPI_COMM_WORLD, IERR)
      CALL MPI_Reduce (NUMHP, NDHGPV, 1, MPI_INTEGER, MPI_SUM, 0,   &
                       MPI_COMM_WORLD, IERR)
      CALL MPI_Reduce (NUMHVK, NDHVK, 1, MPI_INTEGER, MPI_SUM, 0,   &
                       MPI_COMM_WORLD, IERR)
      if (myid == 0) print *
      if (myid == 0) print *, "Sping-angular coefficients involving &
                              HMatrix (CSFG, not spanned):"
      if (myid == 0) print *, "Number of T(ab)                =", NDHGPT
      if (myid == 0) print *, "Number of two-body CSFG-pairs  =", NDHGPV
      if (myid == 0) print *, "Number of obtained Slater (VK) =", NDHVK
!      NDHGPT = NUMHT
!      NDHGPV = NUMHP
!      NDHVK  = NUMHVK
      CALL RALLOC (TCOEFH,    NUMHT,  'TCOEFH', 'SCFmpi')
      CALL RALLOC (LABTH,  2, NUMHT,  'LABTH',  'SCFmpi')
      CALL RALLOC (NVHP,      NUMHP,  'NVHP',   'SCFmpi')
      CALL RALLOC (LABVKH, 5, NUMHVK, 'LABVKH', 'SCFmpi')
      CALL RALLOC (VCOEFH,    NUMHVK, 'VCOEFH', 'SCFmpi')

!      IF (MYID == 0) CALL REPORT_MEM(2, 'SCFmpi 1')

! For (E)OL calculations, determine the generalised occupation numbers.
      CALL NEWCOmpi (WTAEV)

! Set NEC, IECC, and LSKIPPOT, return immediately.
! LSKIPPOT are needed in SETHAM_POT to record the LABELs building 
! D-, X-, and Y- potentials.
! SETLAGmpi needs the generalized occupation numbers to determine
! LSKIPPOT, IECC.
     CALL SETLAGmpi (EOL)

!=======================================================================
! Now LSKIPTPOT obtained, loop over blocks to set:
!     NDCOFOPT, NYCOFOPT, NXCOFOPT, (sizes of the followings:)
! and NDAOPT,   NYAOPT,   NDAOPT

! Initialize the counters for T- and V-coefficients involving potentials
! NTPT:  different T-coefficients for DA potentials
! NCRPT: Spanned by the NTPT ones for DA potentials, NACRPT should be 
!        in general much larger than NATPT for large-scale calulation 
!        using many symmetry-ordered obtials.
! NTPV:  different V-coefficients for Y- and X-potentials
! NCRPV: Spanned by the NTPV ones for R_{rs}(abcd) integrals to 
!        calculate the Y- and X-potentials, NCRPV (in general) >> NTPV
      NTPT = 0;   NCRPT = 0; NTPV = 0; NCRPV = 0; 

! SETHAM_POT constructs only the arrays needed by potential
! calculations.
      LREADMCP  = .FALSE.;  LCHM = .FALSE.;  LCPOT = .FALSE.
      LABTVFRST = .TRUE.
      IF (MYID == 0) WRITE (*, *)
      IF (MYID == 0) WRITE (*, *) &
        "Setting arrays to calculate D(T)-, X-, and Y-potentials ..."

      starttime = MPI_Wtime()
      DO JBLOCK = 1, NBLOCK
        CALL SETHAM_POT(EOL, JBLOCK)
      ENDDO
      endtime = MPI_Wtime()
      if (myid.eq.0)                                                &
        print *, 'Time for setting arrays (s) = ',                  &
                  endtime-starttime

! Spin-angular data to construct potentials
! MPI_Reduce to report the total T, and V-coefficients involving
! potentials
      CALL MPI_Reduce (NTPT , NATPT , 1, MPI_INTEGER, MPI_SUM, 0,   &
                       MPI_COMM_WORLD, IERR)
      CALL MPI_Reduce (NCRPT, NACRPT, 1, MPI_INTEGER, MPI_SUM, 0,   &
                       MPI_COMM_WORLD, IERR)
      CALL MPI_Reduce (NTPV  , NATPV, 1, MPI_INTEGER, MPI_SUM, 0,   &
                       MPI_COMM_WORLD, IERR)
      CALL MPI_Reduce (int8(NCRPV), NSPANVK, 1, MPI_INTEGER8,       &
                       MPI_SUM, 0, MPI_COMM_WORLD, IERR)
      IF (MYID == 0) THEN
        PRINT *
        PRINT *, "Sping-angular coefficients invovling potentials:"
        PRINT *, 'Total number of one-body T (CSFG)         =',  NATPT
        PRINT *, 'Total number of one-body T (CSF, Spanned) =', NACRPT
        PRINT *, 'Total number of two-body VK(CSFG)         =',  NATPV
        PRINT '(1X,A,I12)', &
                'Total number of two-body VK(CSF, Spanned) =', NSPANVK
        PRINT *
      ENDIF
! NXAIORB is not needed anymore, it is used within the first calling
! SETHAM_POT with LABTVFRST = .TRUE.
      CALL DALLOC (NXAIORB, 'NXAIORB', 'SCFmpi')

!=======================================================================
! CYC, July 2022
! Additional modification to make sure that every node has the same
! NXCOFopt, NYCOFopt, NXAopt and NYAopt, then XPOT and YPOT subroutines
! could be parallelized. 
! NDCOFopt, NDAOPT, DAOPT, then DA are calculated serially.
      IF (MYID == 0) WRITE(*,*)"Gathering NXAOPT and NYAOPT arrays ..."
      call system_clock (ncount1, ncount_rate, ncount_max)
      ncount0 = ncount1
!
      CALL NXANYAmpi
!
      call system_clock (ncount2, ncount_rate, ncount_max)
      nsecs = (ncount2-ncount1) / real(ncount_rate)
      if (myid.eq.0) write(*,*) "Gathering NXAOPT and NYAOPT (s): ",nsecs
      if (myid.eq.0) write(*,*)
! CYC, July 2022

! Obtain the largest size, alloc / ralloc the corresponding arrays
      NDDIM0 = MAXVAL(NDCOFopt)
      NXDIM0 = MAXVAL(NXCOFopt)
      NYDIM0 = MAXVAL(NYCOFopt)
      CALL RALLOC (NDAopt, NDDIM0, NW, 'NDAopt','SCFmpi')
      CALL RALLOC (DAopt,  NDDIM0, NW, 'DAopt ','SCFmpi')
      CALL RALLOC (NXAopt, NXDIM0, NW, 'NXAopt','SCFmpi')
      CALL RALLOC (XAopt,  NXDIM0, NW, 'XAopt ','SCFmpi')
      CALL RALLOC (NYAopt, NYDIM0, NW, 'NYAopt','SCFmpi')
      CALL RALLOC (YAopt,  NYDIM0, NW, 'YAopt ','SCFmpi')

      IF (MYID.EQ.0) THEN
        WRITE(*,*)"Size of NDAOPT, DAOPT(:,NW) =", size(NDAOPT,DIM=1)
        WRITE(*,*)"Size of NXAOPT, XAOPT(:,NW) =", size(NXAOPT,DIM=1)
        WRITE(*,*)"Size of NYAOPT, YAOPT(:,NW) =", size(NYAOPT,DIM=1)
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call system_clock (ncount1, ncount_rate, ncount_max)
      nsecs =  (ncount1-ncount0) / real(ncount_rate) 
      if (myid.eq.0) write(*,*)
      if (myid.eq.0) write(*,*) "Prepare NXAOPT and NYAOPT (s): ",nsecs
!-----------------------------------------------------------------------
!CYC Here the maxmum sizes of NDA, DA, NXA, XA, NYA, and YA are
!    obtained, allocate them here, then no re-allocation needed.
!!AK Must initialize DA,XA,YA,... here since MPI routines rely on 
!!   their existence in all processes
      CALL ALLOC (NDA, NDDIM0, 'NDA', 'SCFmpi')
      CALL ALLOC (DA,  NDDIM0, 'DA',  'SCFmpi')
      CALL ALLOC (NXA, NXDIM0, 'NXA', 'SCFmpi')
      CALL ALLOC (XA,  NXDIM0, 'XA',  'SCFmpi')
      CALL ALLOC (NYA, NYDIM0, 'NYA', 'SCFmpi')
      CALL ALLOC (YA,  NYDIM0, 'YA',  'SCFmpi')

!-----------------------------------------------------------------------
!CYC: NXAopt and XAopt, NYAopt and YAopt, are vectorized (ONE-demension) here.
      CALL VECTOR_XYA

!      IF (MYID == 0) CALL REPORT_MEM(3, 'SCFmpi 2')

      ! Vectorized as NXAWRK and XAWRK    
      CALL DALLOC (NXAopt,  'NXAopt',   'SCFmpi')
      CALL DALLOC (XAopt,   'XAopt',    'SCFmpi')
      ! Vectorized as NYAWRK and YAWRK    
      CALL DALLOC (NYAopt,  'NYAopt',   'SCFmpi')
      CALL DALLOC (YAopt,   'YAopt',    'SCFmpi')

!=======================================================================
! DDRS save the generalized weights d_{r,s),  i.e., Eq. (13) in Atom 11,
! (2023) 12, CSFs indexes r and s are those generated by the
! corresponding CSFG pair, they are same for all possible SLATER
! integrals R(abcd; k) among the interacting CSFs pair. 

! Also for the I(ab) integral CSFs pair

! Also for the original GRASP2018 zero-first calculations.
      IF (MYID == 0) WRITE(*, *)'MAXSPAN =', MAXSPAN
      CALL ALLOC (DDRS, MAXSPAN, MAXSPAN, 'DDRS', 'SCFmpi')
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!      IF (MYID == 0) CALL REPORT_MEM(5, 'SCFmpi 3')
!=======================================================================
! The followings are always .FALSE. during the SCF procedure
      LREADMCP  = .false. 
      DVDFIRST  = .false.  
      LABTVFRST = .false.

! Debug
      !LDBPR(25:30) = .true.

      WTAEV0 = 0.0
      DO NIT = 1, NSCF
         call system_clock (ncount1, ncount_rate, ncount_max)
         ncount0 = ncount1
         IF (MYID == 0) WRITE (*, 301) NIT

!----------------------------------------------------------------------
! Section 1: Potentials
! CYC: Obtain the coefficients DAOPT, XAOPT, YAOPT for all orbitals,
!    they will be used directly within SETLAGmpi.f90 :: SETCOF_CYC, or
!    IMPROVmpi :: COFPOTmpi :: SETCOF_CYC
         IF (MYID == 0) WRITE (*, *) &
           "Setting Coefficients for, D(T)-, X-, and Y-potentials ..."
! Set LCPOT and LCHM alternately to perform SCF calculations.
         LCPOT = .TRUE.
         LCHM  = .FALSE.
         starttime = MPI_Wtime()
         DO JBLOCK = 1, NBLOCK
           CALL SETHAM_POT(EOL, JBLOCK)
         ENDDO
         endtime = MPI_Wtime()
         if (myid.eq.0)                                                &
           print *, 'Time for D-, Y-, and X-potentials (s) = ',        &
                  endtime-starttime
!
! Print da, ya, xa for check
!         IF (NIT == 1) CALL PRINTDXY(EOL)
!         STOP "Coding Stop ..."
!
!----------------------------------------------------------------------
! Section 2: Lagrange multipliers
! For all pairs constrained through a Lagrange multiplier, compute
! the Lagrange multiplier
         call system_clock (ncount2, ncount_rate, ncount_max)

         CALL SETLAGmpi (EOL)

         call system_clock (ncount3, ncount_rate, ncount_max)
         nsecs = (ncount3-ncount2) / real(ncount_rate)
         if (myid.eq.0) write(*,*)
         if (myid.eq.0) write(*,*) "Only SETLAGmpi (s): ",nsecs
         nsecs = (ncount3-ncount1) / real(ncount_rate)
         if (myid.eq.0) write(*,*) "SETCOF and SETLAGmpi (s): ",nsecs
         ncount1 = ncount3

!----------------------------------------------------------------------
! Section 3: Improve all orbitals in turn
         DAMPMX = 0.0d0
         IF (MYID == 0) WRITE (*, 302)
         DO J = 1, NW
            JSEQ = IORDER(J)
            IF (LFIX(JSEQ)) CYCLE
            CALL IMPROVmpi (EOL, JSEQ, LSORT, DAMPMX)
         END DO

         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call system_clock (ncount2, ncount_rate, ncount_max)
         nsecs = (ncount2-ncount1) / real(ncount_rate)
         if (myid.eq.0) write(*,*)
         if (myid.eq.0) write(*,*) "IMPROVmpi (s): ",nsecs
         nsecs = (ncount2-ncount0) / real(ncount_rate)
         if (myid.eq.0) write(*,*) &
           "Total time used to update orbitals (s): ",nsecs
         ncount1 = ncount2
!
!   For KOUNT = 1 to NSIC: find the least self-consistent orbital;
!   improve it
!
         DO KOUNT = 1, NSIC
            CALL MAXARR (K)
            IF (K == 0) THEN
               CONVG = .TRUE.
               GO TO 3
            ELSE
               IF (SCNSTY(K) <= ACCY) THEN
                  CONVG = .TRUE.
                  GO TO 3
               ENDIF
            ENDIF
            CALL IMPROVmpi (EOL, K, LSORT, DAMPMX)
         END DO

         CALL MAXARR (K)

         IF (K == 0) THEN
            CONVG = .TRUE.
         ELSE
            IF (SCNSTY(K) <= ACCY) THEN
               CONVG = .TRUE.
            ELSE
               CONVG = .FALSE.
            ENDIF
         ENDIF


    3    CONTINUE
         IF (LDBPR(24) .AND. MYID==0) CALL PRWF (0)

!   Perform Gram-Schmidt process
!   For OL calculation, orthst is true and orbitals are orthonormalized
!   in subroutine improv. For AL calculation, orthst is false.
         IF (.NOT.ORTHST) CALL ORTHSC

!   Write the subshell radial wavefunctions to the .rwf file

         IF (MYID == 0) CALL ORBOUT (RWFFILE2)
!----------------------------------------------------------------------
! Section 4: Construt and diagonalize the H matrix using the updated
! orbital wavefucntions. Check the covergence.
         IF (EOL) THEN
            call system_clock (ncount1, ncount_rate, ncount_max)

            starttime = MPI_Wtime()
            CALL GENINTRKwrap_CSFG (myid, nprocs, j2max)
            !CALL GENINTRKwrap (myid, nprocs, j2max)
            endtime = MPI_Wtime()
            if (myid.eq.0)                                             &
              print *, 'Time for intrk calculation (s) = ',            &
                        endtime-starttime

            LCHM = .TRUE.
            LCPOT = .FALSE.
            CALL MATRIX_CSFG(EOL, DVDFIRST)
            call system_clock (ncount2, ncount_rate, ncount_max)
            nsecs = (ncount2-ncount1) / real(ncount_rate)
            if (myid.eq.0) write(*,*) &
      "All blocks, SETHAM and Matrix Diagonalization (Total, s): ",nsecs
            ncount1 = ncount2
            CALL NEWCOmpi (WTAEV)
         ENDIF

!        Make this a relative convergence test
!        IF(ABS(WTAEV-WTAEV0).LT.1.0D-9.and.
!    &   DAMPMX.LT.1.0D-4) CONVG=.true.
!        PRINT *, 'WTAEV, WTAEV0', WTAEV, WTAEV0,
!    &             ABS((WTAEV-WTAEV0)/WTAEV)
!GG         IF (ABS((WTAEV - WTAEV0)/WTAEV) < 0.001*ACCY) CONVG = .TRUE.
!cjb unified convergence criteria in RMCDHF and RMCDHF_MPI
!cjb     IF(DABS(WTAEV-WTAEV0).LT.1.0D-8.and.              &
!cjb                   DAMPMX.LT.1.0D-2) CONVG=.true.
         IF (ABS((WTAEV - WTAEV0)/WTAEV) < 0.001*ACCY) CONVG = .TRUE.
         !IF (ABS((WTAEV - WTAEV0)/WTAEV) < 0.00005d0*ACCY) CONVG = .TRUE.
         !IF (ABS((WTAEV - WTAEV0)/WTAEV) < 0.000025d0*ACCY) CONVG = .TRUE.
         WTAEV0 = WTAEV

         call system_clock (ncount2, ncount_rate, ncount_max)
         nsecs = (ncount2-ncount0) / real(ncount_rate)
         if (myid.eq.0) write(*,*)
         if (myid.eq.0) write(*,*) &
            "This loop spends time (s): ", nsecs

! R. Si: record memory usage
!      if (myid == 0) write(*,*)'Sleep to record memory usage     & 
!      using shell command top (: RES) ...'
!      CALL SYSTEM('sleep 60')
!      STOP "Stop for monitoring Memory usage (RES) ..."

         IF (.NOT.CONVG .OR. NIT.LT.2) CYCLE

         IF (LDBPR(25) .AND. .NOT.LDBPR(24) .AND. MYID==0) CALL PRWF (0)

         GO TO 5
      END DO

      IF (MYID==0) WRITE (ISTDE,*)' Maximum iterations in SCF Exceeded.'

    5 CONTINUE

!CYC
!      DO I = 31, 32 + KMAXF
!         CLOSE(I)                                ! The MCP coefficient files
!      END DO

      IF (MYID == 0) THEN
         !CLOSE (23)     ! The .rwf file
         CLOSE(25)                               ! The .mix file
      ENDIF
!
!   Complete the summary - moved from rscf92 for easier alloc/dalloc
!
      IF (myid .EQ. 0) CALL ENDSUM

!
!CYC: Output the information of da,xa,ya
!
      IF (MYID == 0) THEN
        WRITE (*,'(A)')
        WRITE (*,'(A)')"Number of DA, XA, YA terms involving potentials"
        WRITE (*,'(A6, 3A10)')"DXY:", "nd", "nx", "ny"
        DO I = 1, NW
           WRITE (*,305)  NP(I), NH(I), NDCOFOPT(I), &
                          NXCOFOPT(I),  NYCOFOPT(I)
        END DO
        WRITE (*,*) 
      ENDIF
!
!   Deallocate storage
!
      CALL DALLOC (WT, 'WT', 'SCFmpi')                       !Either getold or getald

      IF (NEC > 0) THEN
         CALL DALLOC (IECC, 'IECC', 'SCFmpi')
         CALL DALLOC (ECV,  'ECV',  'SCFmpi')
         CALL DALLOC (IQA,  'IQA',  'SCFmpi')
      ENDIF

      IF (NDDIM > 0) THEN
         CALL DALLOC (DA,   'DA',   'SCFmpi')
         CALL DALLOC (NDA,  'NDA',  'SCFmpi')
         NDDIM = 0
      ENDIF

      IF (NXDIM > 0) THEN
         CALL DALLOC (XA,   'XA',   'SCFmpi')
         CALL DALLOC (NXA,  'NXA',  'SCFmpi')
         NXDIM = 0
      ENDIF

      IF (NYDIM > 0) THEN
         CALL DALLOC (YA,   'YA',   'SCFmpi')
         CALL DALLOC (NYA,  'NYA',  'SCFmpi')
         NYDIM = 0
      ENDIF

!=======================================================================
!CYC: DEAllocate the additional arrays, which are not GRASP2018 ones.
      CALL DALLOC (ICOLIDX,'ICOLIDX','SCFmpi')
      CALL DALLOC (DDRS,   'DDRS',   'SCFmpi')

      CALL DALLOC (TCOEFH, 'TCOEFH', 'SCFmpi')
      CALL DALLOC (LABTH,  'LABTH ', 'SCFmpi')
      CALL DALLOC (NVHP,   'NVHP',   'SCFmpi')
      CALL DALLOC (LABVKH, 'LABVKH', 'SCFmpi')
      CALL DALLOC (VCOEFH, 'VCOEFH', 'SCFmpi')

      CALL DALLOC (NDCOFopt,'NDCOFopt', 'SCFmpi')
      CALL DALLOC (NXCOFopt,'NXCOFopt', 'SCFmpi')
      CALL DALLOC (NYCOFopt,'NYCOFopt', 'SCFmpi')
      CALL DALLOC (NDAopt,  'NDAopt',   'SCFmpi')
      CALL DALLOC (DAopt,   'DAopt',    'SCFmpi')
!CYC 2024/08
      CALL DALLOC (NXAWRK,  'NXAWRK',   'SCFmpi')
      CALL DALLOC (XAWRK,   'XAWRK',    'SCFmpi')
      CALL DALLOC (NYAWRK,  'NYAWRK',   'SCFmpi')
      CALL DALLOC (YAWRK,   'YAWRK',    'SCFmpi')
      CALL DALLOC (TMPXYA,  'TMPXYA',   'SCFmpi')
!=======================================================================
      IF (EOL) THEN
         CALL DALLOC (EVAL,    'EVAL',    'SCFmpi')
         CALL DALLOC (EVEC,    'EvEC',    'SCFmpi')
         CALL DALLOC (IATJPO,  'IATJPO',  'SCFmpi')
         CALL DALLOC (IASPAR,  'IASPAR',  'SCFmpi')
         CALL DALLOC (NCMAXBLK,'NCMAXBLK','SCFmpi') ! getold.f
         CALL DALLOC (EAVBLK,  'EAVBLK',  'SCFmpi') ! getold.f
         CALL DALLOC (IDXBLK,  'IDXBLK',  'SCFmpi') ! Allocated in getold.f
         CALL DALLOC (ICCMIN,  'ICCMIN',  'SCFmpi') ! Allocated in items.f<-getold.f
      ENDIF
!
      CALL DALLOC (NCFPAST,   'NCFPAST',   'SCFmpi')
      CALL DALLOC (NCMINPAST, 'NCMINPAST', 'SCFmpi')
      CALL DALLOC (NEVECPAST, 'NEVECPAST', 'SCFmpi')

  301 FORMAT(/,' Iteration number ',1I3,/,' --------------------')
  302 FORMAT(41X,'Self-            Damping'/,&
         'Subshell    Energy    Method   P0    ',&
         'consistency  Norm-1  factor  JP',' MTP INV NNP'/)
  305 FORMAT (1X,1I2,1A2,1X,3I10)

      RETURN
      END SUBROUTINE SCFmpi
