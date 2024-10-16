!***********************************************************************
!                                                                      *
      SUBROUTINE LODCSL_CSF(NCFD,CSF_Number,NCORE,NPEEL,NEXT_CSF)
!                                                                      *
!   Loads the data from the  .csl  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, IQ, ISPAR, ITJPO, JCUP, JQS,   *
!                        PACK, PARSJL, PRSRCN, PRSRSL                  *
!                                                                      *
!   Written by  G. Gaigalas                       NIST, December 2015  *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def,   ONLY:  NNNW
      USE DEBUG_C
      USE DEF_C
      USE ORB_C
      USE STAT_C
      USE TERMS_C,          only: jtab, ntab
      USE IOUNIT_C
      USE BLK_C,            only: NBLOCK,NCFBLK
      USE memory_man
      USE rang_Int_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE prsrsl_I
      USE convrt_I
      USE prsrcn_I
      USE parsjl_I
      USE pack_I
      USE iq_I
      USE jqs_I
      USE jcup_I
      USE itjpo_I
      USE ispar_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,  INTENT(IN)    :: NCORE
      INTEGER,  INTENT(INOUT) :: CSF_Number,NPEEL,NCFD
      LOGICAL,  INTENT(OUT)   :: NEXT_CSF
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NW2 = 2*NNNW
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(NNNW) :: IOCC
      INTEGER , DIMENSION(NW2)  :: IQSUB
      INTEGER , DIMENSION(NNNW) :: JX
      INTEGER :: I, NCF0
      INTEGER :: NCORP1, J, NPJ, NAKJ, LENTH, NREC &
         , IOS, IERR, LOC, NQS, NEWSIZ, ISPARC, NJX, IOC, IPTY, NQSN   &
         , NJXN, NPEELN, NOPEN, JLAST, ILAST, IOCCI, NKJI, IFULLI, NU  &
         , JSUB, IQT, NBEG, NEND, JXN, JPI, II, ITEMP, NCOREL
      LOGICAL :: EMPTY, FULL
      CHARACTER          :: RECL
      CHARACTER(LEN=256) :: RECORD
!-----------------------------------------------
!
!   Initial allocation for arrays with a dimension dependent
!   on the number of CSFs; the initial allocation must be
!   greater than 1
!
      Found_CSF = 0
      NEXT_CSF = .TRUE.
! CYC:      
      !IF(CSF_Number .EQ. 1) NCF = NCFBLK(NBLOCK) + 1
      IF(CSF_Number .EQ. 0) NCF = NCFBLK(NBLOCK) + 1
      NREC = 5
!
    3 CONTINUE
!
      C_shell(NCF)=''
      C_quant(NCF)=''
      C_coupl(NCF)=''
      READ (20, '(A)', IOSTAT=IOS) RECORD
      IF (RECORD(1:2) == ' *') THEN
         ! Check
         IF (CSF_Number /= NUM_in_BLK_Full(NBLOCK)- NUM_in_BLK(NBLOCK)) THEN
      write(*,*)"CSF_Number/=NUM_in_BLK_Full(NBLOCK)-NUM_in_BLK(NBLOCK)"
           WRITE(*, *)"CSF_Number =", CSF_Number
           WRITE(*, *)"NBLOCK =", NBLOCK
           WRITE(*, *)"NUM_in_BLK(NBLOCK) =", NUM_in_BLK(NBLOCK)
           WRITE(*, *)"NUM_in_BLK_Full(NBLOCK) =", & 
                           NUM_in_BLK_Full(NBLOCK)
           CLOSE(1312)
           CLOSE(20)
           STOP "Sth Error !!!"
         ENDIF
         ! End of block   
         NEXT_CSF = .FALSE.
         RETURN
      ENDIF

      C_shell(NCF) = RECORD
      IF (IOS == 0) THEN
!
!   Read in the occupations (q) of the peel shells; stop with a
!   message if an error occurs
!
         CALL PRSRCN (RECORD, NCORE, IOCC, IERR)
         IF (IERR /= 0) GO TO 26
!
!   Read the J_sub and v quantum numbers
!
         READ (20, '(A)', IOSTAT=IOS) RECORD
         IF (IOS /= 0) THEN
            WRITE (ISTDE, *) 'LODCSL_CSF: Expecting subshell quantum', &
               ' number specification;'
            GO TO 26
         ENDIF
         C_quant(NCF) = RECORD
         LOC = LEN_TRIM(RECORD)
         CALL PARSJL (1, NCORE, RECORD, LOC, IQSUB, NQS, IERR)
         IF (IERR /= 0) GO TO 26
!
!   Read the X, J, and (sign of) P quantum numbers
!
         READ (20, '(A)', IOSTAT=IOS) RECORD
         IF (IOS /= 0) THEN
            WRITE (ISTDE, *) 'LODCSL_CSF: Expecting intermediate ', &
               'and final angular momentum'
            WRITE (ISTDE, *) 'quantum number and final parity ', &
               'specification;'
            GO TO 26
         ENDIF
         C_coupl(NCF) = RECORD
! CYC:
! Too time-consuming
!
         IF(NotFound >= 1) THEN
            DO I =1,NCF-1
               IF(Found(I) == 0) THEN
                  IF(C_shell(I) == C_shell(NCF)) THEN
                     IF(C_quant(I) == C_quant(NCF)) THEN
                        IF(C_coupl(I) == C_coupl(NCF)) THEN
                           Found(I) = 1
                           NotFound = NotFound - 1
                           Found_CSF = 1
                           ! RETURN
                           GOTO 3
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            END DO
         ENDIF
! CYC: Output C_shell, C_quant, C_coupl into temp-file
      WRITE (1312, '(A)') TRIM(C_shell(NCF))
      WRITE (1312, '(A)') TRIM(C_quant(NCF))
      WRITE (1312, '(A)') TRIM(C_coupl(NCF))
! CYC: Accounting and storing 
         CSF_Number = CSF_Number + 1
         NCF0 = NCFD + CSF_Number
!
!   Zero out the arrays that store packed integers
!
         DO I = 1,NNNW
            IQA(I,NCF0)    = 0
            JQSA(I,1,NCF0) = 0
            JQSA(I,2,NCF0) = 0
            JQSA(I,3,NCF0) = 0
            JCUPA(I,NCF0)  = 0
         END DO
!
!   Determine the parity and all intermediate and the final
!   angular momentum quantum numbers
!
         DO I = 256, 1, -1
            IF (RECORD(I:I) == ' ') CYCLE
            LOC = I
            EXIT
         END DO
         RECL = RECORD(LOC:LOC)
         IF (RECL == '+') THEN
            ISPARC = 1
         ELSE IF (RECL == '-') THEN
            ISPARC = -1
         ELSE
            WRITE (ISTDE, *) 'LODCSL_CSF: Incorrect parity ', &
                             'specification;'
            GO TO 26
         ENDIF
         LOC = LOC - 1
!
         CALL PARSJL (2, NCORE, RECORD, LOC, JX, NJX, IERR)
         IF (IERR /= 0) GO TO 26
!
!   Set the occupation and subshell quantum number array elements
!   in IQ, JQS for the core subshells
!
         DO I = 1, NCORE
            CALL PACK (NKJ(I) + 1, I, IQA(1:NNNW,NCF0))
            CALL PACK (0, I, JQSA(1:NNNW,1,NCF0))
            CALL PACK (0, I, JQSA(1:NNNW,2,NCF0))
            CALL PACK (1, I, JQSA(1:NNNW,3,NCF0))
         END DO
!
!   Check all subshell, intermediate and final angular momentum
!   quantum numbers; set the array elements in IQ, JQS for the peel
!   subshells; set the coupling array element in JCUP and the total
!   angular momentum array element in ITJPO
!
         IOC = 0
         IPTY = 0
         NQSN = 0
         NJXN = 0
         NPEELN = 0
         NOPEN = 0
         JLAST = 0
         ILAST = 0
         NCORP1 = NCORE + 1
         DO I = NCORP1, NW
            IOCCI = IOCC(I)
            NPEELN = NPEELN + IOCCI
            NKJI = NKJ(I)
            IFULLI = NKJI + 1
            EMPTY = IOCCI == 0
            IF (.NOT.EMPTY) IOC = IOC + 1
            FULL = IOCCI == IFULLI
            IF (EMPTY .OR. FULL) THEN
               NU = 0
               JSUB = 0
            ELSE
               IPTY = IPTY + NKL(I)*IOCCI
               IF (NKJI /= 7) THEN
                  NQSN = NQSN + 1
                  IF (NQSN > NQS) THEN
                     WRITE (ISTDE, *) &
                       'LODCSL_CSF: Too few subshell quantum', &
                        ' numbers specified;'
                     GO TO 26
                  ENDIF
                  NU = 0
                  JSUB = IQSUB(NQSN)
               ELSE
                  IF (IOCCI /= 4) THEN
                     NQSN = NQSN + 1
                     IF (NQSN > NQS) THEN
                        WRITE (ISTDE, *) 'LODCSL_CSF: Too few subshell ', &
                           'quantum numbers specified;'
                        GO TO 26
                     ENDIF
                     NU = 0
                     JSUB = IQSUB(NQSN)
                  ELSE
                     NQSN = NQSN + 1
                     IF (NQSN > NQS) THEN
                        WRITE (ISTDE, *) 'LODCSL_CSF: Too few subshell ', &
                           'quantum numbers specified;'
                        GO TO 26
                     ENDIF
                     JSUB = IQSUB(NQSN)
                     IF (JSUB==4 .OR. JSUB==8) THEN
                        NU = JSUB/2
                        NQSN = NQSN + 1
                        IF (NQSN > NQS) THEN
                           WRITE (ISTDE, *) &
                             'LODCSL_CSF: Too few subshell', &
                              ' quantum numbers specified;'
                           GO TO 26
                        ENDIF
                        JSUB = IQSUB(NQSN)
                     ELSE
                        NU = 0
                     ENDIF
                  ENDIF
               ENDIF
               IQT = MIN(IOCCI,IFULLI - IOCCI)
               LOC = (IFULLI - 2)/2
               LOC = (LOC*(LOC + 1))/2 + IQT
               NBEG = JTAB(LOC+1) + 1
               NEND = JTAB(LOC+2)
! CYC
!      CALL CONVRT (NCF, RECORD, LENTH)
!      WRITE (ISTDE, *) ' CSF sequence number: '//RECORD(1:LENTH)//':'
!      WRITE (ISTDE, *) TRIM(C_shell(NCF))
!      WRITE (ISTDE, *) TRIM(C_quant(NCF))
!      WRITE (ISTDE, *) TRIM(C_coupl(NCF))
!               WRITE (ISTDE,'(5I4,A)')I, IQT, LOC,NBEG,NEND, &
!                      ' =I, IQT, LOC,NBEG,NEND'

               DO J = NBEG, NEND, 3
                  IF (NTAB(J+2) /= JSUB + 1) CYCLE
                  IF (NU == 0) THEN
                     NU = NTAB(J)
                     GO TO 9
                  ELSE
                     IF (NTAB(J) == NU) GO TO 9
                  ENDIF
               END DO
               CALL CONVRT (NP(I), RECORD, LENTH)
               WRITE (ISTDE, *) 'LODCSL_CSF: Subshell quantum numbers ', &
                  'specified incorrectly for '//RECORD(1:LENTH)//NH(I)//&
                  ' subshell.'
               GO TO 26
            ENDIF
    9       CONTINUE
            IF (.NOT.EMPTY .AND. .NOT.FULL) THEN
               NOPEN = NOPEN + 1
               IF (NOPEN > 1) THEN
                  IF (JSUB == 0) THEN
                     JXN = JLAST
                  ELSE
                     ILAST = IOC
                     NJXN = NJXN + 1
                     IF (NJXN > NJX) THEN
                        WRITE (ISTDE, *) &
                        'LODCSL_CSF: Too few intermediate', &
                           ' and final angular momentum', &
                           ' quantum numbers specified;'
                        GO TO 26
                     ENDIF
                     JXN = JX(NJXN)
                     DO J = ABS(JLAST - JSUB), JLAST + JSUB, 2
                        IF (JXN == J) GO TO 11
                     END DO
                     CALL CONVRT (NP(I), RECORD, LENTH)
                     WRITE (ISTDE, *) &
      'LODCSL_CSF: subshell to previous subshells is incorrect.'
                     GO TO 26
                  ENDIF
   11             CONTINUE
                  CALL PACK (JXN + 1, NOPEN - 1, JCUPA(1:NNNW,NCF0))
                  JLAST = JXN
               ELSE
                  JLAST = JSUB
               ENDIF
            ENDIF
            CALL PACK (IOCCI, I, IQA(1:NNNW,NCF0))
            CALL PACK (NU, I, JQSA(1:NNNW,1,NCF0))
            CALL PACK (0, I, JQSA(1:NNNW,2,NCF0))
            CALL PACK (JSUB + 1, I, JQSA(1:NNNW,3,NCF0))
         END DO
!
         DO I = MAX(1,NOPEN), NW
            CALL PACK (0, I, JCUPA(1:NNNW,NCF0))
         END DO
!
         IF (NQSN /= NQS) THEN
            WRITE (ISTDE, *) 'LODCSL_CSF: Too many subshell', &
               ' quantum numbers specified;'
            GO TO 26
         ENDIF
!
         IF (ILAST /= IOC) NJXN = NJXN + 1
         IF (NJXN /= NJX) THEN
            WRITE (ISTDE, *) 'LODCSL_CSF: Too many intermediate', &
            ' and final angular momentum', ' quantum numbers specified;'
            GO TO 26
         ENDIF
!
         IF (JX(NJXN) /= JLAST) THEN
            WRITE (ISTDE, *) 'LODCSL_CSF: Final angular momentum', &
               ' incorrectly specified;'
            GO TO 26
         ENDIF
!
         IPTY = (-1)**IPTY
         IF (IPTY /= ISPARC) THEN
            WRITE (ISTDE, *) 'LODCSL_CSF: Parity specified incorrectly;'
            GO TO 26
         ENDIF
!
         JPI = (JLAST + 1)*IPTY
         CALL PACK (JPI, NNNW, JCUPA(1:NNNW,NCF0))
!
         !IF (NCF > 1) THEN
         IF (CSF_Number > 1) THEN
            IF (NPEELN /= NPEEL) THEN
               WRITE (ISTDE, *) &
                 'LODCSL_CSF: Inconsistency in the number of electrons.'
               GO TO 26
            ENDIF
         ELSE
            NPEEL = NPEELN
         ENDIF
!
!   Successfully read a CSF; update NREC and read another CSF
!
         NREC = NREC + 3
         GO TO 3
         ! RETURN
      ELSE ! End of file
         ! Check
         IF (CSF_Number /= NUM_in_BLK_Full(NBLOCK)- NUM_in_BLK(NBLOCK)) THEN
      write(*,*)"CSF_Number/=NUM_in_BLK_Full(NBLOCK)-NUM_in_BLK(NBLOCK)"
           WRITE(*, *)"CSF_Number =", CSF_Number
           WRITE(*, *)"NBLOCK =", NBLOCK
           WRITE(*, *)"NUM_in_BLK(NBLOCK) =", NUM_in_BLK(NBLOCK)
           WRITE(*, *)"NUM_in_BLK_Full(NBLOCK) =", & 
                           NUM_in_BLK_Full(NBLOCK)
           CLOSE(1312)
           CLOSE(20)
           STOP "Sth Error !!!"
         ENDIF
         NEXT_CSF = .FALSE.
         RETURN
      ENDIF
!!
!!   Store the number of electrons in the COMMON variable
!!
!      NCOREL = 0
!      NCOREL = SUM(NKJ(:NCORE)+1)
!      NELEC = NCOREL + NPEEL
!
!      IF (LDBPA(1)) THEN
!         WRITE (99, *) 'From LODCSL:'
!         DO I = 1, NCF
!            WRITE (99, *) 'CSF ', I
!            WRITE (99, *) 'ITJPO: ', ITJPO(I)
!            WRITE (99, *) 'ISPAR: ', ISPAR(I)
!            WRITE (99, *) 'IQ: ', (IQ(J,I),J=1,NW)
!            WRITE (99, *) 'JQS(1): ', (JQS(1,J,I),J=1,NW)
!            WRITE (99, *) 'JQS(2): ', (JQS(2,J,I),J=1,NW)
!            WRITE (99, *) 'JQS(3): ', (JQS(3,J,I),J=1,NW)
!            WRITE (99, *) 'JCUP: ', (JCUP(J,I),J=1,NW - 1)
!         END DO
!      ENDIF
!!      NEXT_CSF = .FALSE.
!!
!      RETURN
!
   26 CONTINUE
      CALL CONVRT (NCF, RECORD, LENTH)
      WRITE (ISTDE, *) 'Sth Error!!!'
      WRITE (ISTDE, *) 'CSF sequence number: '//RECORD(1:LENTH)//':'
      WRITE (ISTDE, *) TRIM(C_shell(NCF0))
      WRITE (ISTDE, *) TRIM(C_quant(NCF0))
      WRITE (ISTDE, *) TRIM(C_coupl(NCF0))
   29 CLOSE(20)
      CLOSE(1312)
      STOP
!
      END SUBROUTINE LODCSL_CSF
