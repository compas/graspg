!***********************************************************************
!                                                                      *
      SUBROUTINE SETSDA_CSFG (OUTSDAdummy, IOUTT, NELMNTCC, LPRINT30,  &
                                          NB, MYID, NPROCS, FHEAD)
!                                                                      *
!***********************************************************************
!                                                                      *
!   Modified from setsda.f90 for CSFG list by Chongyang Chen, Oct 2023 *
!                                                                      *
!   Output: csfg_mcpIDX.30                                             *
!   Structure:                                                         *
!                                                                      *
!      WRITE (29) 'MCP', NB, NCFTr                                     *
!      WRITE (29) NELMNTCC, NODENCOLS, ncsfDF1, ncsfDF1_core           *
!      WRITE (29) (NODECOLS(I), I = 1, NODENCOLS), &                   *
!                 (IENDC(I),    I = 1, NODENCOLS), &                   *
!                 (IROW(I),     I = 1, NELMNTCC)                       *
!                                                                      *
!                                                                      *
!                                                                      *
!***********************************************************************
!   This routine examines lists                                        *
!                               (IC,IR,npos)                           *
!   to set up the array  IENDC  required by the Davidson eigensolver   *
!   of Stathopoulos and Fischer.                                       *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 10 Dec 1992   *
!   Modified for block interactions by C. F. Fischer        May 1997   *
!   Modified for multi-processors   by Xinghong He       03 Jul 1998   *
!
!  Currently shared by mcpblk, mcpmpi
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:52:18   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE csfg_memory_man
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE orb_C
      USE iounit_C

      Use symexpand_mod, ONLY : TotCSFs_perblock, ncsfDF1
      Use symmatrix_mod, KMAXTmp=>KMAX
      Use symmatrix_restart_C, ONLY : ncsfDF1_core
      USe csfg_tv_C

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE outsdampi_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      EXTERNAL OUTSDAdummy
      INTEGER, INTENT(IN) :: IOUTT
      INTEGER, INTENT(IN) :: NELMNTCC
      INTEGER , INTENT(IN) :: NB
      INTEGER , INTENT(IN) :: MYID
      INTEGER , INTENT(IN) :: NPROCS
      LOGICAL  :: LPRINT30
      CHARACTER , INTENT(IN) :: FHEAD*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MB, IEND, ICLAST, I, IC, NPOS, IERR
      INTEGER, DIMENSION(:), pointer :: IENDC
      INTEGER, DIMENSION(:), pointer :: IROW
      CHARACTER :: MCPLAB*3
      INTEGER :: IOFFSET, ICDUMMY, NELCDUMMY, NNONZ, NCFG
!-----------------------------------------------
      IF (MYID == 0) WRITE (6, *) &
         'Analysing sparse matrix array definition file ...', 30
!
! Allocate storage for IENDC(0:NCF)
! Alloacte IROW
      CALL ALLOC (IENDC, 0, NODENCOLS, 'IENDC', 'MCPMPI_GG')
      CALL ALLOC (IROW, NELMNTCC, 'IROW', 'MCPMPI_GG')
      IENDC(0) = 0
      IOFFSET = 0 
      REWIND(IOUTT)
      DO IC = 1, NODENCOLS
        READ(IOUTT)ICDUMMY, NELCDUMMY, (IROW(IOFFSET+I),I=1,NELCDUMMY)
        IOFFSET = IOFFSET + NELCDUMMY 
        IENDC(IC) = IOFFSET
        ! ======Check=========
        IF (NODECOLS(IC) .NE. ICDUMMY) THEN
          WRITE(*,*)"NODECOLS(IC) .NE. ICDUMMY : ", NODECOLS(IC),ICDUMMY
          STOP "Unexpected NODECOLS(IC) .NE. ICDUMMY ..."
        ENDIF
        IF (IC == NODENCOLS .AND. IOFFSET .NE. NELMNTCC) THEN
          WRITE(*,*)"IC == NODENCOLS .AND. IOFFSET .NE. NELMNTCC ..."
          WRITE(*,*)"IOFFSET, NELMNTCC =", IOFFSET, NELMNTCC
          STOP "IC == NODENCOLS .AND. IOFFSET .NE. NELMNTCC ..."
        ENDIF
        ! ======Check=========
      ENDDO
!
!    write to csfg_mcpXXX.30 file
!
      !OPEN(29, FILE='csfg_'//trim(FHEAD)//'.30', STATUS='OLD',         &
      !         FORM='UNFORMATTED', IOSTAT=IERR, POSITION='APPEND')
      !IF (IERR /= 0) THEN
      !   WRITE (ISTDE, *) ' Error when opening the file mcp.30'
      !   STOP
      !ENDIF

! Modifications needed for CSFG list
      !WRITE (29) 'MCP', NB, NCF
      !WRITE (29) NNONZ
      !WRITE (29) (IENDC(I),I=MYID + 1,NCF,NPROCS), (IROW(I),I=1,NNONZ)
!
! csfg_mcpXXX.30
!
      NCFG = TotCSFs_perblock(2,NB) 
      WRITE (30) 'MCP', NB, NCFG, NCFTr
      WRITE (30) NPAIRS(1,NB), NPAIRS(2,NB)
      WRITE (30) ICGPDET(1:NPAIRS(1,NB))
      WRITE (30) ICGPEND(1:NPAIRS(1,NB))
      WRITE (30) IRGPDET(1:NPAIRS(2,NB))
      WRITE (30) NELMNTCC, NODENCOLS, ncsfDF1, ncsfDF1_core
!
! The followings are reproduced by rmcdhf_csfg_mpi
!
      !WRITE (30) (NODECOLS(I), I = 1, NODENCOLS)
      !WRITE (30) (IENDC(I),    I = 1, NODENCOLS)
      !WRITE (30) (IROW(I),     I = 1, NELMNTCC)

      !WRITE (9030, *) 'MCP', NB, NCFG, NCFTr
      !WRITE (9030, *) NPAIRS(1,NB), NPAIRS(2,NB)
      !WRITE (9030, '(A)')"ICGPDET(1:NPAIRS(1,NB)) :"
      !WRITE (9030, '(11I10)') ICGPDET(1:NPAIRS(1,NB))

      !WRITE (9030, '(A)')"ICGPEND(1:NPAIRS(1,NB)) :"
      !WRITE (9030, '(11I10)') ICGPEND(1:NPAIRS(1,NB))

      !WRITE (9030, '(A)')"IRGPDET(1:NPAIRS(2,NB)) :"
      !WRITE (9030, '(11I10)') IRGPDET(1:NPAIRS(2,NB))

      !WRITE (9030, '(A)')"NELMNTCC, NODENCOLS, ncsfDF1, ncsfDF1_core :"
      !WRITE (9030, '(11I10)') NELMNTCC, NODENCOLS, ncsfDF1, ncsfDF1_core

      !WRITE (9030, '(A)')"(NODECOLS(I), I = 1, NODENCOLS) :"
      !WRITE (9030, '(11I10)') (NODECOLS(I), I = 1, NODENCOLS)

      !WRITE (9030, '(A)')"(IENDC(I),    I = 1, NODENCOLS) :"
      !WRITE (9030, '(11I10)') (IENDC(I),    I = 1, NODENCOLS)

      !WRITE (9030, '(A)')"(IROW(I),     I = 1, NELMNTCC) :"
      !WRITE (9030, '(11I10)') (IROW(I),     I = 1, NELMNTCC)
!
!
!   Deallocate storage
!
      CALL DALLOC (IENDC, 'IENDC', 'SETSDA')
      CALL DALLOC (IROW, 'IROW', 'SETSDA')
!
      CALL OUTSDAdummy (LPRINT30, NELMNTCC, NCFTr)

      RETURN
      END SUBROUTINE SETSDA_CSFG
