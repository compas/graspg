!***********************************************************************
!                                                                      *
      SUBROUTINE LODRES
!                                                                      *
!   Loads the data from the  .res  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!   This subroutine is called by every node and the only differences   *
!   from the non-mpi version are:                                      *
!     1) STOP and stopmpi (not necessary)                              *
!     2) LSE needs coomunication - can be moved out of this routine.   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, GETYN, SETQIC.                         *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
!   Block version by Xinghong He          Last revision:  1 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP
      USE memory_man
      USE decide_C
      USE def_C, ONLY: nelecr, nelec, z, emn, c, accy
      USE grid_C
      USE npar_C
      USE npot_C, ONLY: zz, nnuc
      USE orb_C
      USE wave_C
      USE wfac_C
!CYC      USE where_C
      USE hblock_C
      USE iccu_C,    ONLY: iccutblk
      USE foparm_C
      USE iounit_C
      USE mpi_C
      use symmatrix_mod, ONLY: IMCDF, NCFTOTTr, NODECOLS, NODENCOLS
      use csfg_decide_C, ONLY: LCRITE, MaxMemPerProcs
      use symexpand_mod, ONLY: ncsfDF1
      use symmatrix_restart_C, ONLY: IFNCOL, ncsfDF1_core, NELMNTres
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setqic_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NCFRES, NWRES, NBLOCKRES, I, NP10, J, IOS
      LOGICAL :: YES
      CHARACTER(LEN=1) :: CHAR1
!-----------------------------------------------
!
!   Read the basic parameters of the electron cloud; check these
!   against those deduced from the  .csl file
!
      READ (IMCDF) NELECR, NCFRES, NWRES, NBLOCKRES

! Here, NCF is the total number of CSFs of <name>.g file
!      IF (NELECR/=NELEC .OR. NCFRES/=NCF .OR. NWRES/=NW .OR. NBLOCKRES/=NBLOCK&
!         ) CALL STOPMPI ('lodres: NELEC/NCF/NW does not match',myid)
      IF (NELECR /= NELEC .OR. NCFRES    /= NCFTOTTr .OR.             &
          NWRES  /=NW     .OR. NBLOCKRES /= NBLOCK)                   &
          CALL STOPMPI ('lodres: NELEC/NCF/NW does not match',myid)

!   Only ONE block calculation is allowed in the present implementaion 
      IF (nblockres .NE. 1) THEN
        CALL stopmpi ('lodres: nblockres .Ne. 1', myid)
      ENDIF
!
!   Read the nuclear parameters
!
      READ (IMCDF) Z, EMN
      READ (IMCDF) NPARM, (PARM(I),I=1,NPARM)
      READ (IMCDF) N, (ZZ(I),I=1,N), NNUC

      IF (N > NNNP) CALL STOPMPI ('lodres: N greater than NNNP',myid)
!
!   Read the physical effects specifications
!   iccutblk() is now an array of length nblock.
!
      READ (IMCDF) C, LFORDR, (ICCUTBLK(I),I=1,NBLOCK), LTRANS, WFACT, LVP, &
         LNMS, LSMS
!
!   Read the remaining parameters controlling the radial grid and the
!   grid arrays
!
      NP10 = N + 10
      READ (IMCDF) RNT, H, HP, (R(I),I=1,NP10), (RP(I),I=1,NP10), (RPOR(I),I=1,&
         NP10)
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      ACCY = H**6
!
!   Set up the coefficients for the numerical procedures
!
      CALL SETQIC
!
!   Allocate storage for the radial wavefunction arrays
!
      CALL ALLOC (PF, NNNP,NW, 'PF', 'LODMIX')
      CALL ALLOC (QF, NNNP,NW, 'QF', 'LODMIX')
!
!   Read the orbital wavefunctions and the associated arrays
!
      DO J = 1, NW
         READ (IMCDF) E(J), GAMA(J), PZ(J), MF(J)
         READ (IMCDF) (PF(I,J),I=1,MF(J)), (QF(I,J),I=1,MF(J))
      END DO
!
!   Determine if the self-energy contribution is to be estimated
!
      IF (MYID == 0) THEN
         WRITE (ISTDE, *) 'Estimate contributions from the self-energy?'
         LSE = GETYN()
         IF (LSE) THEN
           WRITE(734,'(a)') 'y            ! Estimate contributions from the self-energy?'
         ELSE
           WRITE(734,'(a)') 'n            ! Estimate contributions from the self-energy?'
         END IF

      WRITE(istde,'(/,A)') &
          "Reduce the accuracy of eigenvalue to speedup the calculation"
      WRITE(istde,'(A)') &
          "in case of the calculation performed for rmixacculate_csfg"
      WRITE(istde,'(A)')&
          "If yes, CRITE of maneigmpi.f90 would be set as 1.0d-12."
          READ(*,'(A1)')CHAR1
          IF (CHAR1.EQ.'Y' .OR. CHAR1.EQ.'y') LCRITE = .TRUE.
          WRITE(734,'(A,A)')CHAR1, ' ::   LCRITE '

         WRITE (istde,*)"Input MaxMemPerProcs (in GB) ..."
         READ *, MaxMemPerProcs
      ENDIF
      CALL MPI_BCAST (LSE, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERR)
      CALL MPI_BCAST (LCRITE, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERR)
      CALL MPI_BCAST (MaxMemPerProcs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
!
!     Load data from ncolsXXX.txt
! 
      READ (IFNCOL, *)NELECR, NCFRES, NWRES, nblockres
      IF (NELECR .NE. NELEC .OR. NCFRES .NE. NCFTOTTr .OR. &
          NWRES .NE. NW .OR. nblockres .NE. nblock) THEN
         CALL stopmpi ('lodncol: NELEC/NCF/NW does not match ', myid)
      ENDIF
      READ (IFNCOL, *)ncsfDF1, NELMNTres, ncsfDF1_core
      READ (IFNCOL, *)NODENCOLS
      ALLOCATE(NODECOLS(NODENCOLS))
      DO I = 1, NODENCOLS
        READ(IFNCOL, *)NODECOLS(I)
      ENDDO
      CLOSE(IFNCOL)

      RETURN
      END SUBROUTINE LODRES
