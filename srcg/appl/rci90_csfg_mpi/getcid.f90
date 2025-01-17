!***********************************************************************
!                                                                      *
      SUBROUTINE GETCID (isofile, rwffile, idblk)
!                                                                      *
!   Interactively determines the data governing the CI problem.        *
!   iccut is replaced by an array iccutblk(1:nblock)                   *
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, NUCPOT, RADGRD, SETQIC.                *
!               [RCI92]: SETISO, SETRWF.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!   Block version by Xinghong He          Last revision: 15 Jun 1998   *
!                                                                      *
!   Modified for CSFG-version by Chongyang Chen,  Dec 2021             *
!                                                                      *
!   Last modification for CSFG by Chongyang Chen  Dec 2023             *
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
      USE decide_C
      USE def_C,           ONLY: cvac, c, z, accy, nelec, emn
      USE default_C
      USE grid_C
      USE hblock_C
      USE iccu_C,          ONLY: iccutblk
      USE foparm_C
      USE iounit_C
      USE npar_C
      USE npot_C,          ONLY: zz, nnuc
      USE nsmdat_C
      USE orb_C
      USE wave_C
      USE wfac_C
      USE blim_C
!CYC      USE where_C
      USE qedcut_C
      USE mpi_C
      USE symexpand_mod,   ONLY: nmaxgen
      use symmatrix_mod,   ONLY: IMCDF, NCFTOTTr, NCFBLKTr
      use csfg_decide_C,   ONLY: NDISCARDBR, LDISCARDBR,               &
                                 NDISBR, LDISBR, LSYMD, MaxMemPerProcs
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setisompi_I
      USE setrwfmpi_I
      IMPLICIT NONE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      CHARACTER(LEN=*) :: isofile, rwffile
      CHARACTER(LEN=8), DIMENSION(*) ::  idblk
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       LOGICAL :: getyn, yes
       INTEGER :: i, jblock, ntmp, np10, j
!-----------------------------------------------------------------------
!
! Open, check, load data from, and close the  .iso  file
! Data loaded are: EMN,Z,/NPAR/,/NSMDAT/ where /.../ means whole
! common block.
!
!      IF (MYID == 0) PRINT *, 'Calling SETISO ...'
      CALL SETISOmpi (isofile)
!
! The speed of light, if non-default then spread from node-0
! Quantities to be obtained: C
!
      IF (NDEF .NE. 0) THEN
!cjb node0 input
        IF (myid .EQ. 0) THEN
         WRITE (istde,*) 'Revise the physical speed of light (',CVAC,  &
                           ' in a.u.) ?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE (istde,*) 'Enter the revised value:'
            READ *,C
!PERJ   
            WRITE(734,'(A)') 'y   ! revise speed of light'
            WRITE(734,*) C
!PERJ END
         ELSE
!PERJ   
            WRITE(734,'(A)') 'n   ! revise speed of light'
!PERJ END
            C = CVAC
         ENDIF
        ENDIF
        CALL MPI_Bcast (C, 1, MPI_DOUBLE_PRECISION, 0, &
                               MPI_COMM_WORLD, ierr)
      ELSE
         C = CVAC
      ENDIF
!
! Treat some CSF's as perturbation ? Broadcasting is used only in
! non-default mode, as the above case for speed of light.
! Quantities to be obtained: LFORDR
!
      IF (NDEF .NE. 0) THEN
         IF (myid .EQ. 0) THEN
            WRITE (istde,*) 'Treat contributions of some CSFs',       &
                            ' as first-order perturbations?'
            LFORDR = GETYN ()
!PERJ
            IF (LFORDR) THEN
               WRITE(734,'(A)') 'y  ! some CSFs as first-order pert' 
            ELSE
               WRITE(734,'(A)') 'n  ! some CSFs as first-order pert' 
            ENDIF
!PERJ END
         ENDIF
         CALL MPI_Bcast (LFORDR,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      ELSE
         LFORDR = .FALSE.
      ENDIF
!cyc================================================================
      LSYMD = .FALSE.
      IF (LFORDR) THEN
       IF (myid .EQ. 0) THEN
        WRITE (istde,*) 'Include symmetry-ordered-diagonal block? (y/n)'
        LSYMD = GETYN ()
!PERJ
        IF (LSYMD) THEN
           WRITE(734,'(A)') 'y  ! sym-ordered diag. block' 
        ELSE
           WRITE(734,'(A)') 'n  ! sym-ordered diag. block' 
        ENDIF
!PERJ END
       ENDIF
       CALL MPI_Bcast (LSYMD,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      ENDIF
!cyc================================================================
! Get iccutblk() from the user-input

      IF (.NOT. LFORDR) THEN
         !...Default first
         DO i = 1, nblock
           iccutblk(i) = ncfblk(i)
         ENDDO
      ELSE

      ! Let master do the i/o, then broadcast
         IF (myid .EQ. 0) THEN
!PERJ            WRITE (istde,*) 'There are ', nblock, 'blocks. They are:'
!PERJ            WRITE (istde,*) '  block     J Parity   No of CSFGs'// &
!PERJ                            '    TrueCSFs'
!PERJ            DO i = 1, nblock
!PERJ               WRITE (istde,*) i, idblk(i)(1:5), ncfblk(i), NCFBLKTr(I)
!PERJ            ENDDO

            WRITE (istde,*)
            WRITE (istde,*) 'Enter iccut (in CSFG list) for each block:'
            DO jblock = 1, nblock
               WRITE (istde,*) 'Block ', jblock, '  NCSF(G) = ',       &
                              ncfblk(jblock), ' id = ', idblk(jblock)(1:5)
  123          READ (istdi,*) ntmp
               IF (ntmp .GE. 0 .AND. ntmp .LE. ncfblk(jblock)) THEN
                  iccutblk(jblock) = ntmp
               ELSE
                  WRITE (istde,*) 'ICCUT out of range, re-enter:'
                  GOTO 123
               ENDIF
               write(734,*) ntmp,'! ICCUT for block',jblock
            ENDDO
         ENDIF

         CALL MPI_Bcast (iccutblk,nblock,MPI_INTEGER,0,            &
                                               MPI_COMM_WORLD,ierr)
      ENDIF ! .NOT. LFORDR

!*****************************************************************
!
! Pre-run ?
!
!     IF (IPRERUN .EQ. 0) THEN

!        WRITE (istde,*) ' Prerun with limited interaction?'
!        YES = GETYN ()

!        IF (YES) THEN
!           IPRERUN = 1
!           LTRANS = .FALSE.
!           LVP = .FALSE.
!           LNMS = .FALSE.
!           LSMS = .FALSE.
!           LSE = .FALSE.

!           WRITE (istde,*)  ' Give CSL cut'
!           READ *, NCSFPRE
!           WRITE (istde,*)  ' Give coefficient cut for H_0'
!           READ *, COEFFCUT1
!           WRITE (istde,*) ' Give coefficient cut for the transvers'
!    &,                  ' interaction'
!           GOTO 99
!        ENDIF
!     ENDIF
!*****************************************************************
!
! Include transverse ?
! Quantities to be obtained: LTRANS, WFACT
!
      IF (myid .EQ. 0) THEN
        WRITE (istde,*) 'Include contribution of H (Transverse)?'
        LTRANS = GETYN ()
        WRITE (istde,*) 'Modify all transverse photon frequencies?'
        YES = GETYN ()
        IF (YES) THEN
           WRITE (istde,*) 'Enter the scale factor:'
           READ *, WFACT
        ELSE
           WFACT = 1.0D00
        ENDIF
!        IF (LTRANS) THEN   ! To keep the rci.inp same, no matter
!        include Breit interaction or not. 
          ! Discard the Breit contributions invovling very high n orbital?
          WRITE (istde,*) 'Limit the Breit contributions by n?'
          NDISCARDBR = GETYN () 
          ! Discard Breit interaction for n above (with l > 1):   
          NDISBR = 15
101       WRITE (istde,*) 'Discard Breit interaction for above n: '
          READ *, NDISBR
          IF (NDISBR.LT.1) THEN
            WRITE (istde,*) 'Input error,  n should be larger than 0 !'
            GOTO 101
          ENDIF

          ! Discard the Breit contributions invovling very high l orbital?
          WRITE (istde,*) 'Limit the Breit contributions by l?' 
          LDISCARDBR = GETYN ()
          ! Discard Breit interaction for l above (affects also s and p orbitals):
102          WRITE (istde,*) 'Discard Breit interaction for above l: '
          LDISBR = 14
          READ *, LDISBR 
          IF (LDISBR.LT.0) THEN
            WRITE (istde,*) 'Input error, l should be no less than 0 !'
            GOTO 102
          ENDIF
        ENDIF
!      ENDIF ! myid .EQ. 0
      CALL MPI_Bcast (LTRANS, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (WFACT, 1, MPI_DOUBLE_PRECISION, 0,            &
                           MPI_COMM_WORLD, ierr)

      CALL MPI_Bcast (NDISCARDBR, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (LDISCARDBR, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (NDISBR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (LDISBR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
!
! Other interactions ? One logical for each case. Done altogether
!
      IF (myid .EQ. 0) THEN
         WRITE (istde,*) 'Include H (Vacuum Polarisation)?'
         LVP = GETYN ()

         WRITE (istde,*) 'Include H (Normal Mass Shift)?'
         LNMS = GETYN ()

         WRITE (istde,*) 'Include H (Specific Mass Shift)?'
         LSMS = GETYN ()

         WRITE (istde,*) 'Estimate self-energy?'
         LSE = GETYN ()
         IF (LSE.EQV..TRUE.) THEN
             NQEDCUT = 1
!CYC            WRITE (istde,*)                                           &
!CYC        'Largest n quantum number for including self-energy for orbital'
!CYC            WRITE (istde,*) 'n should be less or equal 8'
!CYC            READ *, NQEDMAX
         ELSE
            NQEDCUT = 0
         END IF
!CYC     Needing NQEDMAX / the nmax value of spectroscopy orbital for
!the frequent-dependent Breit contributions, nomatter including
!Self-energy or not. 
         WRITE (istde,*)                                               &
             'Largest n (NQEDMAX) quantum number for including '//     &
                'self-energy for orbital'
         WRITE (istde,*)                                               &
        'It is also label the max n-value for spectroscopy orbitals'
         WRITE(istde,*) &
           "In addition 1: it should be not larger than nmaxgen"
202      WRITE (istde,*) 'In addition 2: it should be less or equal 8'
         READ *, NQEDMAX
         IF (NQEDMAX.GT.nmaxgen) THEN
           WRITE(*,*) &
           "Error! NQEDMAX should be not larger than nmaxgen"
           WRITE(*,*)"Input again: "
           GOTO 202
         ENDIF

         WRITE (istde,*)"Input MaxMemPerProcs (in GB) ..."
         READ *, MaxMemPerProcs

         IF (LTRANS) THEN
          WRITE(734,'(a)')'y            ! Contribution of H Transverse?'
         ELSE
          WRITE(734,'(a)')'n            ! Contribution of H Transverse?'
         END IF
         IF (YES) THEN
           WRITE(734,'(a)') 'y            ! Modify photon frequencies?'
           WRITE(734,*) WFACT,'! Scale factor'
         ELSE
           WRITE(734,'(a)') 'n            ! Modify photon frequencies?'
         END IF

         ! PERJ Added
         IF (NDISCARDBR) THEN
          WRITE(734,'(a)')'y            ! Discard Breit for high n?'
         ELSE
          WRITE(734,'(a)')'n            ! Discard Breit for high n?'
         END IF
         WRITE(734,*) NDISBR, '! Discard Breit interaction above n'
         IF (LDISCARDBR) THEN
          WRITE(734,'(a)')'y            ! Discard Breit for high l?'
         ELSE
          WRITE(734,'(a)')'n            ! Discard Breit for high l?'
         END IF
         WRITE(734,*) LDISBR, '! Discard Breit interaction above l'
         ! PERJ END
          
         IF (LVP) THEN
           WRITE(734,'(a)') 'y            ! Vacuum polarization?'
         ELSE
           WRITE(734,'(a)') 'n            ! Vacuum polarization?'
         END IF
         IF (LNMS) THEN
           WRITE(734,'(a)') 'y            ! Normal mass shift?'
         ELSE
           WRITE(734,'(a)') 'n            ! Normal mass shift?'
         END IF
         IF (LSMS) THEN
           WRITE(734,'(a)') 'y            ! Specific mass shift?'
         ELSE
           WRITE(734,'(a)') 'n            ! Specific mass shift?'
         END IF
         IF (LSE) THEN
           WRITE(734,'(a)') 'y            ! Self energy?'
           WRITE(734,*) NQEDMAX, '! Max n for including self energy'
         ELSE
           WRITE(734,'(a)') 'n            ! Self energy?'
         END IF


         ! PERJ Added
         WRITE(734,*) MaxMemPerProcs, '! MaxMemPer Procs (in GB)'
         ! PERJ END
      ENDIF ! myid .EQ. 0

      CALL MPI_Bcast (LVP,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (LNMS, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (LSMS, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (LSE,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (MaxMemPerProcs,  1, MPI_INTEGER, 0, &
                      MPI_COMM_WORLD, ierr)

!CYC  NQEDMAX will be also used within bessel.f90
      CALL MPI_Bcast (NQEDMAX, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
!
! Parameters controlling the radial grid
!
!cff-cjb grid
!cff use default grid
!
! 99  IF (NPARM .EQ. 0) THEN
!        RNT = EXP (-65.D0/16.D0) / Z
!        H = 0.5D0**4
!        N = MIN (220,NNNP)
!     ELSE
!CFF     .. should be Z-dependent
         RNT = 2.0D-06/Z
!         RNT = 2.0D-06
         H = 5.D-2
         N = NNNP
!     ENDIF
      HP = 0.D0
!     IF ( NDEF.NE.0) THEN
!        IF (myid .EQ. 0) THEN
!           WRITE (istde,*) 'The default radial grid parameters',      &
!                          ' for this case are:'
!           WRITE (istde,*) ' RNT = ',RNT,';'
!           WRITE (istde,*) ' H = ',H,';'
!           WRITE (istde,*) ' HP = ',HP,';'
!           WRITE (istde,*) ' N = ',N,';'
!           WRITE (istde,*) ' revise these values?'
!           YES = GETYN ()
!        ENDIF
!
!         ...To prevent subsequent BCAST when YES is false
!        CALL MPI_Bcast (YES, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

!        IF (YES) THEN
!           IF (myid .EQ. 0) THEN
!              WRITE (istde,*) 'Enter RNT:'
!              READ *, RNT
!              WRITE (istde,*) 'Enter H:'
!              READ *, H
!              WRITE (istde,*) 'Enter HP:'
!              READ *, HP
!              WRITE (istde,*) 'Enter N:'
!              READ *, N
!           ENDIF
!           CALL MPI_Bcast (RNT, 1, MPI_DOUBLE_PRECISION, 0,    &
!                               MPI_COMM_WORLD, ierr)
!           CALL MPI_Bcast (H,   1, MPI_DOUBLE_PRECISION, 0,    &
!                               MPI_COMM_WORLD, ierr)
!           CALL MPI_Bcast (HP,  1, MPI_DOUBLE_PRECISION, 0,    &
!                               MPI_COMM_WORLD, ierr)
!           CALL MPI_Bcast (N,   1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!        ENDIF
!     ENDIF ! NDEF.NE.0
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      ACCY = H**6
!
!   Set up the coefficients for the numerical procedures
!
      CALL SETQIC
!
!   Generate the radial grid and all associated arrays
!
      CALL RADGRD
!
!   Generate $- r \times V_ (r)$
!
      CALL NUCPOT
!
!   Load the radial wavefunctions
!     
      CALL SETRWFmpi (rwffile)
!
!   Write the basic parameters of the model electron cloud to the
!   .res  file; this is the second record on the file --- the
!   first is the header (see SUBROUTINE SETRES)
!
! CYC:
      !WRITE (imcdf) NELEC, NCF, NW, nblock ! ncf is ncftot, from setcsll
      WRITE (imcdf) NELEC, NCFTOTTr, NW, nblock
!
!   Write the nuclear parameters and $- r \times V_ (r)$
!   to the  .res  file; these are the third, fourth, and fifth
!   records on the file
!
      WRITE (imcdf) Z, EMN
      WRITE (imcdf) NPARM,(PARM(I), I = 1, NPARM)
      WRITE (imcdf) N, (ZZ(I), I = 1, N), NNUC
!
!   Write the physical effects specification to the  .res  file.
!   iccutblk() is now an array of length nblock.
!   This is the sixth record on the file
!
      WRITE (imcdf) C, LFORDR, (ICCUTblk(i), i = 1, nblock),           &
                    LTRANS, WFACT, LVP, LNMS, LSMS
!
!   Write the grid data to the  .res  file; this is the seventh
!   record on the file
!
      NP10 = N + 10
      WRITE (imcdf) RNT, H, HP, (R(I), I = 1, NP10),                   &
                 (RP(I), I = 1, NP10), (RPOR(I), I = 1, NP10) ! (imcdf)
!
!   Write out the interpolated radial wavefunctions; there are
!   2*NW such records; thus the total number of records written
!   at this point is 7+2*NW
!
      DO J = 1, NW
         WRITE (imcdf) E(J), GAMA(J), PZ(J), MF(J)
         WRITE (imcdf) (PF(I,J), I = 1, MF(J)), (QF(I,J), I = 1, MF(J))
      ENDDO
!
      RETURN
      END SUBROUTINE GETCID
