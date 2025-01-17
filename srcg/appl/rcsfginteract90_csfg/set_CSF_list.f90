!***********************************************************************
!                                                                      *
      SUBROUTINE SET_CSF_list(NCORE,NPEEL)
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [RCI92]: LENGTH, OPENFL.                               *
!                                                                      *
!   Written by  G. Gaigalas                   NIST, December 2015      *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!    M O D U L E S
!-----------------------------------------------
      USE IOUNIT_C
      USE ORB_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      USE prsrsl_I
      USE convrt_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,  INTENT(OUT) :: NCORE, NPEEL
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IERR, IOS, I, J, LENTH, NAKJ, NCORP1, NPJ
      LOGICAL :: FOUND
      CHARACTER :: NAME*24, FILNAM*256, RECORD*15, DEFNAM*11, FORM*11, STATUS*3
      CHARACTER :: RECORD_1*15, RECORD_2*15
      CHARACTER :: S_closed_1*181, S_closed_2*181
      CHARACTER :: S_orbitals_1*2070, S_orbitals_2*2070
!-----------------------------------------------
!
      WRITE (6, *) 'Loading Configuration Symmetry List File ...'
      FILNAM = 'rcsfgmr.inp'
      FORM = 'FORMATTED'
      STATUS = 'OLD'
      CALL OPENFL (21, FILNAM, FORM, STATUS, IERR)
      IF (IERR == 1) THEN
         WRITE (6, *) 'Error when opening', FILNAM
         STOP
      ENDIF
      ! STRING "Core subshells: "
      READ (21, '(1A15)', IOSTAT=IOS) RECORD_1
      IF (IOS/=0 .OR. RECORD_1(1:15)/='Core subshells:') THEN
         WRITE (6, *) 'Not a Configuration Symmetry List File;'
         WRITE (6, *) 'RECORD_1 =', RECORD_1
         CLOSE(21)
      ENDIF
!
      FILNAM = 'rcsfg.inp'
      FORM = 'FORMATTED'
      STATUS = 'OLD'
      CALL OPENFL (20, FILNAM, FORM, STATUS, IERR)
      IF (IERR == 1) THEN
         WRITE (6, *) 'Error when opening', FILNAM
         STOP
      ENDIF
      ! String "Core subshells: "
      READ (20, '(1A15)', IOSTAT=IOS) RECORD_2
      IF (IOS/=0 .OR. RECORD_2(1:15)/='Core subshells:') THEN
         WRITE (6, *) 'Not a Configuration Symmetry List File;'
         WRITE (6, *) 'RECORD_2 =', RECORD_2
         CLOSE(20)
      ENDIF
!
      ! Informations of core subshells 
      READ (21, '(A)') S_closed_1
      ! Informations of core subshells 
      READ (20, '(A)') S_closed_2
      IF(S_closed_1 /= S_closed_2) then
         STOP "Different close shells"
      end if
!
      FILNAM = 'rcsfg.out'
      FORM = 'FORMATTED'
      STATUS = 'OLD'
      CALL OPENFL (22, FILNAM, FORM, STATUS, IERR)
      WRITE (22, '(1A15)') RECORD_1
      I = LEN_TRIM(S_closed_1)
      WRITE (22,'(A)') S_closed_1(1:I)
!
      ! String "Peel subshells:"
      READ (21, '(1A15)', IOSTAT=IOS) RECORD_1
      ! String "Peel subshells:"
      READ (20, '(1A15)', IOSTAT=IOS) RECORD_2
      WRITE (22,'(1A15)') RECORD_1
!
      ! Informations of peel subshells 
      READ (21, '(A)') S_orbitals_1
      ! Informations of peel subshells 
      READ (20, '(A)') S_orbitals_2
!      I = LEN_TRIM(S_orbitals_1)
!      IF(S_orbitals_1(1:I) /= S_orbitals_2(1:I)) then
!         STOP "Different order of Peel orbitals"
!      end if
      ! Informations of peel subshells 
      WRITE (22, '(A)') TRIM(S_orbitals_2)
!
      READ (21, '(1A7)', IOSTAT=IOS) RECORD_1
      READ (20, '(1A7)', IOSTAT=IOS) RECORD_2
      ! String "CSF(s):"
      WRITE (22,'(1A7)') RECORD_1
!---------------------------------------------------------------------------------
      ! NUM_in_BLK for MR-CSFs
      CALL Set_CSF_number

      ! NUM_in_BLK_full for Full-CSFs
      CALL Set_CSF_number_full
!      write(*,*)"Exiting from Set_CSF_number ..."

      REWIND(20)
      READ (20, '(1A15)', IOSTAT=IOS) RECORD_2
!
!   Get the list of subshells
!
      NW = 0
!
!   Read the list of core subshells; set up the arrays NP, NAK,
!   NKL, NKJ, NH for these subshells
!
      CALL PRSRSL (20, 1)
      NCORE = NW
      NCORP1 = NW + 1
!
!   Skip the peel subshell identification header; read the list of
!   peel subshells; set up the arrays NP, NAK, NKL, NKJ, NH for
!   these subshells
!
      READ (20, *)
      CALL PRSRSL (20, 2)
      NPEEL = NW - NCORE
!
!   Ensure that the sets of core and peel subshell are disjoint
!
      DO J = NCORE + 1, NW
         NPJ = NP(J)
         NAKJ = NAK(J)
         DO I = 1, NCORE
            IF (NP(I)/=NPJ .OR. NAK(I)/=NAKJ) CYCLE
            WRITE (ISTDE, *) 'SET_CSF_list: The lists of core and', &
               ' peel subshells must form disjoint sets.'
            STOP
         END DO
      END DO
!
!   Print the number of relativistic subshells
!
      IF (NW > 1) THEN
         CALL CONVRT (NW, RECORD, LENTH)
         WRITE (6, *) 'There are '//RECORD(1:LENTH)// &
                      ' relativistic subshells;'
      ELSE
         WRITE (6, *) 'There is 1 relativistic subshell;'
      ENDIF
      READ (20, *)
!
      RETURN
      END SUBROUTINE SET_CSF_list
