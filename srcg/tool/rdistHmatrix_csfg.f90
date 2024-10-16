MODULE PUBLICDATA
  IMPLICIT NONE
  INTEGER, DIMENSION(:), ALLOCATABLE :: NODENCOLS_P, ncsfDF1_P
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: NODECOLS_P
  
  INTEGER :: NELC
  DOUBLE PRECISION :: ELSTO 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EMT
  INTEGER, DIMENSION(:), ALLOCATABLE :: IROW
  INTEGER :: NODENCOLS
END MODULE

!=======================================================================
!                                                                      *
  program rdistHmatrix_csfg
!                                                                      *
! Program to redistribute the H-matrix using CSFGs,                    *
! Generally, this program is uesed to reduce the number of             *
! MPI processes ( CPU cores) within MANEIG, as the MPI performance     *
! is excellent to calculate H matrixelements, whereas it is inefficient*
! to search the eigen-pairs.                                           * 
!                                                                      *
! In addition, it will also reduce the memory needed to use less MPI   *
! processes to search the eigen-pairs.                                 *
!                                                                      *
! This program is accompanied with rci_csfg_block_mpi.                 *
!                                                                      *
! Written by Chongyang Chen, Fudan University, December 2021           *
!                                                                      *
!                                                                      *
!   Input  files: <Source>/XXX/rciXXX.res,                             *
!                 <Source>/XXX/ncolsXXX.txt                            *
!   Output files: <Target>/XXX/rciXXX.res,                             *
!                 <Target>/XXX/ncolsXXX.txt                            *
!                                                                      *
!   "XXX" is ID of MPI process, such as "000", "001", etc.             *
!                                                                      *
!=======================================================================
  use PUBLICDATA
  IMPLICIT NONE 
  INTEGER :: I, J, JJ, K, L, M, N, IR, NUMFS, NUMFT, IFLAG, III, IERR
  CHARACTER(LEN=512) :: DSOURCE, DTARGET, FNSOURCE, FNTARGET
  CHARACTER(LEN=512) :: FSRES, FSCOL, FTRES, FTCOL
  CHARACTER(LEN=3)   :: CIDS, CIDT
  LOGICAL            :: FRST, LDF1
  
  INTEGER :: NELEC, NCF, NW, NBLOCK
  INTEGER :: ncsfDF1, ncsfDF1_core, ncsfDF1tot_core
  INTEGER :: NTOTCOLS
  INTEGER*8 :: NELMNTres, NTOTEMNT, NELMNT
  INTEGER :: NELECR, NCFRES, NWRES, nblockres
  INTEGER :: IMCDF, M0, IJ, ICOLTMP
  
  IFLAG = 0
  !DSOURCE = "/home/cychen/HMatrix_tmpSave"
  !WRITE (*,'(A)')"Input the directory of the source HM files, such as: "
  !WRITE (*,'(A)')"HMatrix_tmpSave"
  WRITE (*,'(A)')"Input the directory of the source HM files: "
  READ (*,'(A)')Dsource
  Dsource=adjustl(adjustr(Dsource))
  CALL TESTFILES(DSOURCE, IFLAG)
  
  !DTARGET = "/home/cychen/tmp"
  !WRITE (*,'(A)')"Input the directory of the target HM files, such as:"
  !WRITE (*,'(A)')"MPI_TMP"
  WRITE (*,'(A)')"Input the directory of the target HM files:"
  READ (*,'(A)')Dtarget
  Dtarget=adjustl(adjustr(Dtarget))
  CALL TESTFILES(DTARGET, IFLAG)
  
  WRITE (*,'(A)')"Input the number of source HM files:"
  READ (*,*)NUMFS
  
  WRITE (*,'(A)')"Input the number of target HM files:"
  READ (*,*)NUMFT
  
  IF (MOD(NUMFS, NUMFT) .NE. 0) THEN
    WRITE (*,'(A)')"Error, numfs should be dividable by numft ..."
    STOP "Error, numfs should be dividable by numft ..."
  ENDIF

  N = NUMFS/NUMFT
  ALLOCATE(NODENCOLS_P(N))
  AlLOCATE(ncsfDF1_P(N))
  
  DO I = 1, NUMFT
     WRITE (CIDT, '(I3.3)') I-1
     FNTARGET = TRIM(DTARGET) // '/' // CIDT
     IFLAG = 1
     CALL TESTFILES(FNTARGET, IFLAG)
      
     ! Create the needed directory
     IF (IFLAG .NE. 0) THEN
       III = LEN_TRIM(FNTARGET)
       IERR = SYSTEM ('mkdir -p -m 775 ' // FNTARGET(1:III))
       IF (IERR .NE. 0) THEN
         WRITE (*,*)"Directory can not be created: ", FNTARGET
         STOP "Error! Directory can not be created ... "
       ENDIF
     ENDIF
  
     ! Open the target files to write
     FTRES = TRIM(FNTARGET) // '/' // 'rci' // CIDT // '.res'
     FTCOL = TRIM(FNTARGET) // '/' // 'ncols' // CIDT // '.txt' 
     !WRITE (*, '(A)') TRIM(FTRES)
     !WRITE (*, '(A)') TRIM(FTCOL)
     OPEN (20, FILE=TRIM(FTRES), STATUS='UNKNOWN', FORM='UNFORMATTED') 
     OPEN (21, FILE=TRIM(FTCOL), STATUS='UNKNOWN', FORM='FORMATTED')
  
     NTOTCOLS = 0
     NTOTEMNT = 0
     ncsfDF1tot_core = 0

     DO J = 1, N
       K = (J-1) * NUMFT + I
       WRITE (CIDS, '(I3.3)') K-1
       FNSOURCE = TRIM(DSOURCE) // '/' // CIDS
       IFLAG = 0
       CALL TESTFILES(FNSOURCE, IFLAG)
       FSCOL = TRIM(FNSOURCE) // '/' // 'ncols' // CIDS // '.txt'
       CALL TESTFILES(FSCOL, IFLAG)
  
       ! Open the source files to read, obtain: NCF, maximum value of NODENCOLS
       OPEN (16, FILE=TRIM(FSCOL), STATUS='OLD', FORM='FORMATTED') 
       
       ! Read ncolsXXX.txt
       READ (16, *) NELEC, NCF, NW, NBLOCK
       READ (16, *) ncsfDF1, NELMNT, ncsfDF1_core
       READ (16, *) NODENCOLS
  
       NCSFDF1_P(J) = ncsfDF1_core
       NODENCOLS_P(J) = NODENCOLS
       NTOTCOLS = NTOTCOLS + NODENCOLS
       NTOTEMNT = NTOTEMNT + NELMNT
       ncsfDF1tot_core = ncsfDF1tot_core + ncsfDF1_core
  
       IF (.NOT.ALLOCATED(NODECOLS_P)) ALLOCATE(NODECOLS_P(NCF,N))
       DO L = 1, NODENCOLS
         READ (16, *) NODECOLS_P(L,J)
       ENDDO
  
       CLOSE(16)
     ENDDO
     M=MAXVAL(NODENCOLS_P(1:N))

     ! Output the head lines of ncolsXXX.txt
     WRITE (21, *) NELEC, NCF, NW, NBLOCK, "  = NELEC, NCF, NW, NBLOCK" 
     WRITE (21, *) ncsfDF1, NTOTEMNT, ncsfDF1tot_core, "  = ncsfDF1, NELMNT, ncsfDF1_core " 
     WRITE (21, *) NTOTCOLS, "  = NODENCOLS"
     
     DO J = 1, N
       FRST = .FALSE.
       IF (J.EQ.1) FRST = .TRUE.
  
       K = (J-1) * NUMFT + I
       WRITE (CIDS, '(I3.3)') K-1
       FNSOURCE = TRIM(DSOURCE) // '/' // CIDS
       IFLAG = 0
       CALL TESTFILES(FNSOURCE, IFLAG)
       FSRES = TRIM(FNSOURCE) // '/' // 'rci' // CIDS // '.res'
       CALL TESTFILES(FSRES, IFLAG)
  
       ! Open the source files to read
       IMCDF = 30 + J
       OPEN (IMCDF, FILE=TRIM(FSRES), STATUS='OLD', FORM='UNFORMATTED')
       IF (FRST) &
         WRITE (*,*) "Write the head lines of rciXXX.res ...", CIDT 
       CALL READRES(FRST, NUMFS, NUMFT, I-1, K-1, IMCDF)
     ENDDO

     ! Read sequently the details of source ncolsXXX.txt and rciXXX.res, write
     ! them into the target ncolsXXX.txt and rciXXX.res. 
     IF (.NOT.(ALLOCATED(EMT))) ALLOCATE(EMT(NCF))
     IF (.NOT.(ALLOCATED(IROW))) ALLOCATE(IROW(NCF))
     DO L = 1, M
       DO J = 1, N
        IMCDF = 30 + J
        IF (L.GT.NODENCOLS_P(J)) CYCLE
        READ (IMCDF) NELC, ELSTO, (EMT(IR), IR=1,NELC), (IROW(IR), IR=1,NELC)
        IF (NODECOLS_P(L,J).NE. IROW(NELC)) THEN
          WRITE(*,*)"Unexpected NODECOLS_P(L,J).NE. IROW(NELC) ... ", &
                    NODECOLS_P(L,J), IROW(NELC)
          STOP "Unexpected NODECOLS_P(L,J).NE. IROW(NELC) ..."
        ENDIF
        WRITE (20) NELC, ELSTO, (EMT(IR), IR=1,NELC), (IROW(IR), IR=1,NELC) 
        WRITE (21, *) NODECOLS_P(L, J)
       ENDDO
     ENDDO
 
     CLOSE(20)
     CLOSE(21)
  ENDDO

  DEALLOCATE(NODENCOLS_P)
  DEALLOCATE(ncsfDF1_P)
  DEALLOCATE(NODECOLS_P)
  DEALLOCATE(EMT)
  DEALLOCATE(IROW)
  
  STOP "Normal Exit ..."

end program rdistHmatrix_csfg

SUBROUTINE TESTFILES(FDNAME, IFLAG)
  IMPLICIT NONE
  CHARACTER*(*) FDNAME
  INTEGER :: ISTATE, IFLAG
  
  ISTATE = ACCESS (TRIM(FDNAME), ' ')
  IF (ISTATE .NE. 0 .AND. IFLAG .EQ. 0) THEN
    WRITE(*, '(A31, A)')"Directory/File does not exist: ", TRIM(FDNAME)
    STOP "ERROR! Directory/File does not exist ..."
  ENDIF
  IFLAG = ISTATE
  
  RETURN
END

SUBROUTINE READRES(FRST, NSPROCS, NTPROCS, MYIDT, MYIDS, IMCDF)
  use PUBLICDATA
  IMPLICIT NONE
  INTEGER, PARAMETER :: NNN1=2010, NNNP=2000, NNNW=214
  INTEGER :: NTPROCS, NSPROCS, MYIDT, MYIDS, IMCDF
  INTEGER :: I, J, K, L, N, ICCUTblk(50)
  INTEGER :: NELECR, NCFRES, NWRES, nblockres, NELEC, NCF, NW, NBLOCK
  REAL*8  :: Z, EMN
  INTEGER :: NPARM, NP10, NNUC
  REAL*8  :: PARM(2), ZZ(NNNP)
  REAL*8  :: C,WFACT
  LOGICAL :: FRST
  LOGICAL :: LFORDR, LTRANS, LVP, LNMS, LSMS
  REAL*8  :: RNT, H, HP, R(NNN1), RP(NNN1), RPOR(NNN1)
  REAL*8  :: E(NNNW), GAMA(NNNW), PZ(NNNW)
  INTEGER :: MF(NNNW)
  REAL*8  :: PF(NNNP,NNNW), QF(NNNP,NNNW)
  INTEGER :: IR
  INTEGER :: ncfdum, iccutdum, myiddum, nprocsdum
  INTEGER :: IOS 
  CHARACTER(LEN=6) R92RES
  
  READ (imcdf, IOSTAT = IOS) R92RES
  IF ((IOS .NE. 0) .OR. (R92RES .NE. 'R92RES')) THEN
     CLOSE (imcdf)
     STOP 'ERROR: Not RCI92 .res file ....'
  ENDIF
  IF (FRST) WRITE (20) R92RES
  
  READ (imcdf) NELECR, NCFRES, NWRES, nblockres
  IF (FRST) WRITE (20) NELECR, NCFRES, NWRES, nblockres
  NELEC = NELECR
  NCF   = NCFRES
  NW    = NWRES
  NBLOCK= NBLOCKRES
  
  READ (imcdf) Z,EMN
  IF (FRST) WRITE (20)  Z,EMN
  
  READ (imcdf) NPARM,(PARM(I),I = 1,NPARM)
  IF (FRST) WRITE (20) NPARM,(PARM(I),I = 1,NPARM)
  
  READ (imcdf) N,(ZZ(I),I = 1,N),NNUC
  IF (FRST) WRITE(20) N,(ZZ(I),I = 1,N),NNUC
  
  IF (N .GT. NNNP) STOP 'ERROR! N greater than NNNP ...'
  READ (imcdf) C, LFORDR, (ICCUTblk(i), i = 1, nblock),  &
               LTRANS, WFACT, LVP, LNMS, LSMS
  IF (FRST) WRITE (20)  C, LFORDR, (ICCUTblk(i), i = 1, nblock),  &
               LTRANS, WFACT, LVP, LNMS, LSMS
  
  NP10 = N+10
  READ (imcdf) RNT,H,HP, &
       (R(I),I = 1,NP10), (RP(I),I = 1,NP10), (RPOR(I),I = 1,NP10)
  IF (FRST) WRITE (20) RNT,H,HP, &
       (R(I),I = 1,NP10), (RP(I),I = 1,NP10), (RPOR(I),I = 1,NP10)
  
  DO J = 1, NW
     READ (imcdf) E(J), GAMA(J), PZ(J), MF(J)
     READ (imcdf) (PF(I,J), I = 1, MF(J)), (QF(I,J), I = 1, MF(J))
  ENDDO
  IF (FRST) THEN
    DO J = 1, NW
       WRITE (20) E(J), GAMA(J), PZ(J), MF(J)
       WRITE (20) (PF(I,J), I = 1, MF(J)), (QF(I,J), I = 1, MF(J))
    ENDDO
  ENDIF
  
  READ (IMCDF) ncfdum, iccutdum, myiddum, nprocsdum
  IF (MYIDDUM .NE. MYIDS .OR. NPROCSDUM .NE. NSPROCS) THEN
    WRITE(*,'(A)')"ERROR! MYIDDUM .NE. MYIDS .OR. NPROCSDUM .NE. NSPROCS ..."
    WRITE(*,*)MYIDDUM, MYIDS, NPROCSDUM, NSPROCS
    STOP "ERROR! MYIDDUM .NE. MYIDS .OR. NPROCSDUM .NE. NSPROCS ..."
  ENDIF
  IF (FRST) WRITE (20) ncfdum, iccutdum, MYIDT, NTPROCS  ! myiddum ==> MYIDT; nprocsdum==>NTPROCS
  
  RETURN
END
