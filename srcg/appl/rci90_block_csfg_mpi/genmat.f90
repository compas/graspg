!***********************************************************************
!                                                                      *
      SUBROUTINE GENMAT(ATWINV, JBLOCK, MYID, NPROCS, ELSTO, IRESTART, &
                        SLF_EN)
!                                                                      *
!        Generate Hamiltonian matrix for all blocks                    *
!        This routine calls setham to do the computation. It makes     *
!        sure that the hamiltonian matrix is complete.                 *
!  Only node-0 has the correct, non-zero elsto. And in any case, elsto *
!  is not added to the matrix elements in files *.res                  *
!  This routine has been re-written to work in both single- and        *
!  multi- processors.                                                  *
!                                                                      *
! Xinghong He 1998-06-23                                               *
!   Modify  by G Gaigalas                         May 2021             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17

!                                                                      *
!   Last modification for CSFG by Chongyang Chen  Dec 2023             *
!                                                                      *
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE eigv_C
!GG      USE iccu_C
      USE foparm_C
      USE orb_C
!CYC      USE where_C
      USE hblock_C
      USE hmat_C
      USE iounit_C
      use symmatrix_mod, ONLY: IMCDF, NCFTr, NODENCOLS, NODECOLS
      use symmatrix_restart_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE posfile_I
      USE setham_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: JBLOCK
      INTEGER  :: MYID
      INTEGER  :: NPROCS
      INTEGER, INTENT(OUT) :: IRESTART
      REAL(DOUBLE)  :: ATWINV
      REAL(DOUBLE)  :: ELSTO
      REAL(DOUBLE), DIMENSION(*) :: SLF_EN
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IREAD, IOS, NCFDUM, ICCUTDUM, MYIDDUM, NPROCSDUM, I, &
         IOS2, NELC, IR, IROWDUM, IPOS, NROWS, J, ICSTRT
      REAL(DOUBLE) :: STOEL, DUM, EAV0

      character idstring*3
      character msg*80, cnelemnt*30
!-----------------------------------------------
!
      NELMNT = 0                      ! Counting continues in setham
      EAV    = 0.D0
      ELSTO  = 0.D0

! See how much had been done (Hamiltonian matrix)
! irestart is set;
! iread accumulated;
! nelmnt, eav, elsto obtained (to be further modified in setham)

      IREAD = 0                       ! # of rows read, initialization necessary
!
! CYC: By now, the RESTART mode has been implemented for ONE block. 
!

      READ (IMCDF, IOSTAT=IOS) NCFDUM, ICCUTDUM, MYIDDUM, NPROCSDUM
      IRESTART = IOS

      IF (IOS == 0) THEN

         IF (NCFTr/=NCFDUM .OR. ICCUT/=ICCUTDUM .OR. MYID/=MYIDDUM &
            .OR. NPROCS /= NPROCSDUM) THEN
            WRITE (ISTDE, *) NCF, NCFDUM, ICCUT, ICCUTDUM, MYID, &
                   MYIDDUM, NPROCS, NPROCSDUM, 'check'
            STOP 'genmat:1'
         ENDIF

!CYC: CSFGs-version implementation
         IF (.NOT.allocated(NODECOLS)) THEN
           WRITE (*,*) &
      "Error! Values of NODENCOLS, NODECOLS not obtained ..."
           WRITE (*,*)"Not one ReRun case? ..........."
           STOP "Error! Values of NODENCOLS, NODECOLS not obtained ..."
         ENDIF
         !IF (MYID .EQ. 0) WRITE(*,*)"GENMAT: NODENCOLS = ", NODENCOLS
         !DO I = MYID + 1, NCF, NPROCS
         DO I = 1, NODENCOLS
            READ (IMCDF, IOSTAT=IOS2, ERR=180) &
               NELC, STOEL, (DUM,IR=2,NELC), EAV0, (IROWDUM,IR=1,NELC)
                            ! Lower triangle row-mode, diagonal last
            IF (IOS2 == 0) THEN
               !IF (MYID .EQ. 0) WRITE(251,*)"GENMAT: I, ELSTO = ", I, STOEL
               IREAD = IREAD + 1
               NELMNT = NELMNT + NELC
               EAV = EAV + EAV0
! CYC correct, 20220730
               ! ELSTO = STOEL
               IF (I.EQ.1 .AND. MYID.EQ.0) ELSTO = STOEL
!Check:
               IF (IROWdum .NE. NODECOLS(I)) THEN
       WRITE(*,*)'Wrong column index....', &
      'I,NODECOLS(I), IROW(NELC)=', I,NODECOLS(I),IROWdum
      STOP 'Unexpected IROWdum.NE.NODECOLS(I) in genmat.f Line 89 ...'
               ENDIF
            ELSE
               EXIT
            ENDIF
         END DO
         IPOS = 7 + NW + NW + IREAD + 1
!Check:
         IF (IREAD .NE. NODENCOLS) THEN
           WRITE(*,*)"Error to obtain the whole H matrixelements ...."
           STOP "Error to obtain the whole H matrixelements ...."
         ELSE
           IF (MYID == 0) WRITE(*,*) &
      "Obtain succesfully the whole H matrixelements myid ...", myid
         ENDIF
      ELSE
        IPOS = 7 + NW + NW
      ENDIF

! Find the maximum number of rows

!      NROWS = (NCFTr - MYID - 1 + NPROCS)/NPROCS
!      IF (NCFTr < NPROCS) NROWS = NCFTr/(MYID + 1)

! Report the number of rows read.
! A more suitable report on all nodes can be done here, but this will
! set a synchronization point.

!      IF (MYID == 0) THEN
!         WRITE (ISTDE, *)                                              &
!                       IREAD, ' (total ', NROWS, ') rows read from .res'
!         WRITE (24, *) IREAD, ' (total ', NROWS, ') rows read from .res'
!      END IF

! Position the file for the next record from setham
! CYC: Only ONE block, commented
!      DO I = 1, JBLOCK - 1
!         J = (NCFBLK(I) - MYID - 1 + NPROCS)/NPROCS
!         IF (NCFBLK(I) < NPROCS) J = NCFBLK(I)/(MYID + 1)
!         IPOS = IPOS + J + 1
!      END DO
      CALL POSFILE (0, IMCDF, IPOS)

      IF (IOS /= 0) THEN
         ! New calculation for this block
         WRITE (IMCDF) NCFTr, ICCUT, MYID, NPROCS
         ! ncf, iccut are specific to the current block
         ICSTRT = MYID + 1
!     ...Generate the full Hamiltonian matrix
         CALL SETHAM (MYID, NPROCS, JBLOCK, ELSTO, ICSTRT, NELMNT, ATWINV, &
            SLF_EN)
      ELSE
         NELMNTTMP = NELMNT
         NCFTMP = NCFTr
         WRITE (idstring, '(I3.3)') myid
         WRITE (cnelemnt, '(I20)')NELMNTtmp
         msg = ' ' // 'ID: ' // idstring // ' nelemnt= ' // &
              adjustl(adjustr(cnelemnt))
         CALL mpix_printmsg (msg, myid, nprocs)
         IF (NELMNT .NE. NELMNTres) THEN
           WRITE (*,*)"Unexpected NELMNT .NE. NELMNTres, myid ...", &
                      NELMNT, NELMNTres, myid
           STOP "Unexpected NELMNT .NE. NELMNTres ..."
         ENDIF
      ENDIF

      RETURN

180   WRITE (*,*)"Error read RES file, myid, NodeNCols, read ", &
                 myid, NODENCOLS, i-1
      STOP "Error read RES file .........."

      END SUBROUTINE GENMAT
