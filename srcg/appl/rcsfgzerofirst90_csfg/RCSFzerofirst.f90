!***********************************************************************
!                                                                      *
      PROGRAM RCSFGzerofirst_CSFG
!-----------------------------------------------
!                                                                      *
!   From a set of CSLs this program identifies the ones that           *
!   interact with a given multireference                               *
!                                                                      *
!   This program is a slight modification of the GENMCP program        *
!                                                                      *
!   Written by  G. Gaigalas                     Vilnius, May 2016      *
!                                                                      *
!   CSFG version by Chongyang Chen                                     *
!                            Fudan University, Shanghai, Sep 2021      *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE symzf_mod

      USE BLK_C,            only: NBLOCK,NCFBLK
      USE rang_Int_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE set_CSF_ZFlist_I
      USE lodcsl_Zero_I
      USE lodcsl_Part_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: NEXT_BLOCKZ, NEXT_BLOCKF
      INTEGER :: CSF_Number, ncount1, I, J, K, IBLOCK, NFOUND
      CHARACTER(LEN=256) :: LINEIN, LINEOUT
      INTEGER :: M, N
!-----------------------------------------------
      call starttime (ncount1, 'RCSFGzerofirst')
!
      print *, ""
      print *, "RCSFGzerofirst_CSFG:"
      print *, "            Takes a list of CSFs and partitions each symmetry"
      print *, "            block into a zero- and first-order CSF space from"
      print *, "            a zero-order list."
      print *, "            (C)   Copyright by G. Gaigalas and Ch. F. Fischer"
      print *, "            (Fortran 95 version)                NIST  (2021)."
      print *, "                                                                "
      print *, "            CSFG version, Chongyang Chen,       Fudan (2023)."
      print *, ""
      print *, "            Input files:    list with CSF(G)s to be partitioned"
      print *, "                            list with CSF(G)s defining"
      print *, "                                 the zero-order space"
      print *, "            Input files:    together with their labeling files"
      print *, "                                                       "
      print *, "            Output file:    rcsfg.out"
      print *, "            Output file:    icut"
      print *, "                                                        "
      NBLOCK = 0

      CALL SET_CSF_ZFlist
      open(301,file='icut',form='formatted',status='unknown')
      WRITE (6, *) "         Block       Zero-order Space        &
                   Full space    CHECK"
      REWIND(20)
      REWIND(21)
      DO I = 1, 5
        READ(20,*)
        READ(21,*)
      ENDDO
      NBLOCK = 0
      DO IBLOCK = 1, NBlockZ
         !write(301,*)NCFBLK(NBLOCK)
         write(301,*)NUM_in_BLKZ(IBLOCK)
         NUMZ = NUM_in_BLKZ(IBLOCK)
         NUMF = NUM_in_BLKF(IBLOCK)
         !WRITE (6,'(3X,I2,6X,I14,3X,I17)')IBLOCK,NUMZ, NUMF

! Loading the Zero-space CSFs
! Allocate the working arrays
         ALLOCATE(Found(NUMZ))
         ALLOCATE(ZCONF1(NUMZ))
         ALLOCATE(ZCONF2(NUMZ))
         ALLOCATE(ZCONF3(NUMZ))
         ALLOCATE(LZCONF1(NUMZ))
         ALLOCATE(LZCONF2(NUMZ))
         ALLOCATE(LZCONF3(NUMZ))
         CALL LODCSL_Zero (NEXT_BLOCKZ)

         IF (NBLOCK.NE.IBLOCK .OR. NOTFOUND.NE.NUMZ) THEN
          WRITE(*,*)'Unexpected NBLOCK.NE.IBLOCK .OR. NOTFOUND.NE.NUMZ'
          STOP 'Unexpected NBLOCK.NE.IBLOCK .OR. NOTFOUND.NE.NUMZ ...'
         ENDIF

! Preprocessing Zero-space CSFs
         ALLOCATE(ZCONF1B(NUMZ))
         ALLOCATE(ZCONF2B(NUMZ))
         ALLOCATE(ZCONF3B(NUMZ))
         ALLOCATE(LZCONF1B(NUMZ))
         ALLOCATE(LZCONF2B(NUMZ))
         ALLOCATE(LZCONF3B(NUMZ))
         DO I = 1, NUMZ
           LINEIN = TRIM(ZCONF1(I))
           M = LEN_TRIM(LINEIN)
           LZCONF1(I) = M
           CALL REMOVEBLANK(LINEIN, M, LINEOUT, N, I)
           ZCONF1B(I) = LINEOUT(1:N)
           LZCONF1B(I) = N

           LINEIN = TRIM(ZCONF2(I))
           M = LEN_TRIM(LINEIN)
           LZCONF2(I) = M
           CALL REMOVEBLANK(LINEIN, M, LINEOUT, N, I)
           ZCONF2B(I) = LINEOUT(1:N)
           LZCONF2B(I) = N

           LINEIN = TRIM(ZCONF3(I))
           M = LEN_TRIM(LINEIN)
           LZCONF3(I) = M
           CALL REMOVEBLANK(LINEIN, M, LINEOUT, N, I)
           ZCONF3B(I) = LINEOUT(1:N)
           LZCONF3B(I) = N
         ENDDO

! Loading the Full-space CSFs
         ALLOCATE(FLAGWRITE(NUMF))
         ALLOCATE(FCONF1(NUMF))
         ALLOCATE(FCONF2(NUMF))
         ALLOCATE(FCONF3(NUMF))
         ALLOCATE(LFCONF1(NUMF))
         ALLOCATE(LFCONF2(NUMF))
         ALLOCATE(LFCONF3(NUMF))
         CALL LODCSL_Full (NEXT_BLOCKF)

! Preprocessing Full-space CSFs
         ALLOCATE(FCONF1B(NUMF))
         ALLOCATE(FCONF2B(NUMF))
         ALLOCATE(FCONF3B(NUMF))
         ALLOCATE(LFCONF1B(NUMF))
         ALLOCATE(LFCONF2B(NUMF))
         ALLOCATE(LFCONF3B(NUMF))
         DO I = 1, NUMF
           LINEIN = TRIM(FCONF1(I))
           M = LEN_TRIM(LINEIN)
           LFCONF1(I) = M
           CALL REMOVEBLANK(LINEIN, M, LINEOUT, N, I)
           FCONF1B(I) = LINEOUT(1:N)
           LFCONF1B(I) = N

           LINEIN = TRIM(FCONF2(I))
           M = LEN_TRIM(LINEIN)
           LFCONF2(I) = M
           CALL REMOVEBLANK(LINEIN, M, LINEOUT, N, I)
           FCONF2B(I) = LINEOUT(1:N)
           LFCONF2B(I) = N

           LINEIN = TRIM(FCONF3(I))
           M = LEN_TRIM(LINEIN)
           LFCONF3(I) = M
           CALL REMOVEBLANK(LINEIN, M, LINEOUT, N, I)
           FCONF3B(I) = LINEOUT(1:N)
           LFCONF3B(I) = N
         ENDDO

         Found = 0
         NFOUND = 0
         FLAGWRITE = .true.

! Compare the Zero- and Full-CSF, determine IFLAGWRITE
         ! First loop for the Full-space, much more faster
         DO J = 1, NUMF
          DO I = 1, NUMZ
           ! The comparisons implemented within the original code, very
           ! time consuming, Changed as the follows.
           !IF (ZCONF1(I)/=FCONF1(J)) CYCLE
           !IF (ZCONF2(I)/=FCONF2(J)) CYCLE
           !IF (ZCONF3(I)/=FCONF3(J)) CYCLE

           ! Compare firstly the integers
           IF (LZCONF1(I)/=LFCONF1(J)) CYCLE
           IF (LZCONF1B(I)/=LFCONF1B(J)) CYCLE
           IF (LZCONF2(I)/=LFCONF2(J)) CYCLE
           IF (LZCONF2B(I)/=LFCONF2B(J)) CYCLE
           IF (LZCONF3(I)/=LFCONF3(J)) CYCLE
           IF (LZCONF3B(I)/=LFCONF3B(J)) CYCLE

           ! Compare the strings
           !IF (ZCONF1B(I)(1:LZCONF1B(I))/=FCONF1B(J)(1:LFCONF1B(J)))
           !CYCLE
           !IF (ZCONF2B(I)(1:LZCONF2B(I))/=FCONF2B(J)(1:LFCONF2B(J)))
           !CYCLE
           !IF (ZCONF3B(I)(1:LZCONF3B(I))/=FCONF3B(J)(1:LFCONF3B(J)))
           !CYCLE
           IF (ZCONF1B(I)/=FCONF1B(J)) CYCLE
           IF (ZCONF2B(I)/=FCONF2B(J)) CYCLE
           IF (ZCONF3B(I)/=FCONF3B(J)) CYCLE

           ! Here found the Zero-CSF within the full-space
           FOUND(I) = 1
           NFOUND = NFOUND + 1
           FLAGWRITE(J) = .FALSE.
          ENDDO
         ENDDO
         WRITE(*,*)"block",  IBLOCK, " NUMZ=", NUMZ,    &
                   " NUMF=", NUMF," NFOUND=", NFOUND
         ! Check:
         IF (NFOUND/=NUMZ) THEN
          WRITE(*,*)"Warning !!!!!!!!!!!!!!!!!!!!!!!"
          WRITE(*,*)"Some Zero-CSFs are not found in the Full-space"
          WRITE(*,*)"The first FIVE unfound CSFs are:"
          J = 0
          DO I = 1, NUMZ
            IF (FOUND(I).EQ.0) THEN
              J = J +1
      WRITE(*,*)"J=", J, " ICSF=", I, " NOT FOUND IN FULL SPACE"
      WRITE(*,'(A)')TRIM(ZCONF1(I))
      WRITE(*,'(A)')TRIM(ZCONF2(I))
      WRITE(*,'(A)')TRIM(ZCONF3(I))
              IF (J.EQ.5) EXIT
            ENDIF
          ENDDO
          ! For CSFGs, program stopped with error, 5 not found
          ! CSFs output
          IF (flag_sym_csf) then
       WRITE(*,*)"ERROR, Some Zero-CSFs are not found in the Full-space"
       STOP "ERROR, Some Zero-CSFs are not found in the Full-space"
          endif
         ENDIF

         ! Write the first-order CSFs into rcsf.out
         DO I = 1, NUMF
          IF (FLAGWRITE(I)) THEN
            WRITE(22,'(A)')TRIM(FCONF1(I))
            WRITE(22,'(A)')TRIM(FCONF2(I))
            WRITE(22,'(A)')TRIM(FCONF3(I))
          ENDIF
         ENDDO

         IF (IBLOCK.NE.NBlockZ) WRITE(22,'(A2)') ' *'

! Deallocate the working arrays
         DEALLOCATE(Found)
         DEALLOCATE(ZCONF1)
         DEALLOCATE(ZCONF2)
         DEALLOCATE(ZCONF3)
         DEALLOCATE(LZCONF1)
         DEALLOCATE(LZCONF2)
         DEALLOCATE(LZCONF3)
         DEALLOCATE(ZCONF1B)
         DEALLOCATE(ZCONF2B)
         DEALLOCATE(ZCONF3B)
         DEALLOCATE(LZCONF1B)
         DEALLOCATE(LZCONF2B)
         DEALLOCATE(LZCONF3B)
         DEALLOCATE(FLAGWRITE)
         DEALLOCATE(FCONF1)
         DEALLOCATE(FCONF2)
         DEALLOCATE(FCONF3)
         DEALLOCATE(LFCONF1)
         DEALLOCATE(LFCONF2)
         DEALLOCATE(LFCONF3)
         DEALLOCATE(FCONF1B)
         DEALLOCATE(FCONF2B)
         DEALLOCATE(FCONF3B)
         DEALLOCATE(LFCONF1B)
         DEALLOCATE(LFCONF2B)
         DEALLOCATE(LFCONF3B)
      END DO
      close(301)

      call stoptime (ncount1, 'RCSFGzerofirst_CSFG')
      STOP
      END PROGRAM RCSFGzerofirst_CSFG
