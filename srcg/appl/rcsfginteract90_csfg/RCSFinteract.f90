!***********************************************************************
!                                                                      *
      PROGRAM RCSFGinteract_CSFG
!-----------------------------------------------
!                                                                      *
!   From a set of CSLs this program identifies the ones that           *
!   interact with a given multireference                               *
!                                                                      *
!   This program is a slight modification of the GENMCP program        *
!                                                                      *
!   Written by  G. Gaigalas                   NIST, December 2015      *
!                                                                      *
!   CSFG version by Chongyang Chen                                     *
!                           Fudan University, Shanghai, July 2023      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE memory_man
      USE default_C
      USE BLK_C,            only:  NBLOCK,NCFBLK
      USE orb_C
      USE STAT_C
      USE rang_Int_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE Interact_MR_I
      USE set_CSF_list_I
      USE lodcsl_MR_I
      USE lodcsl_CSF_I
      USE factt_I
      USE Interact_CSF_I

      USE testfiles_I
      
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: NEXT_BLOCK, NEXT_CSF
      INTEGER :: ICOLBREI,NCORE_MR,NPEEL_MR,NCORE,NPEEL,CSF_Number,NCFD
      INTEGER :: I_Count, ncount1

      INTEGER :: I 
      INTEGER :: JA, JB, int_CSF
      logical, dimension(:), allocatable  :: flagwrite
      CHARACTER(LEN=256) :: line1, line2, line3
!-----------------------------------------------
      call starttime (ncount1, 'RCSFGinteract_CSFG')
!
      print *, ""
      print *, &
        "RCSFinteract: Determines all the CSFs (rcsf.inp) that interact"
      print *, &
        "              with the CSFs in the multireference (rcsfmr.inp)"
      print *, &
        "              (C)  Copyright by G. Gaigalas and Ch. F. Fischer"
      print *, &
        "              (Fortran 95 version)               NIST  (2017)."

      print *, ""
      print *, &
        "              CSFG version, Chongyang Chen,      Fudan (2023)."
      print *, "              Input files: rcsfgmr.inp, rcsfg.inp"
      print *, "                           rcsfgmr.l,  rlabel.inp"
      print *, "              Output file: rcsfg.out"
      print *, ""

!CYC-2024
      print *, "Checking rcsfgmr.l and rlabel.inp ..."
      CALL TESTFILES('rcsfgmr.l')
      CALL TESTFILES('rlabel.inp')
      I = system('diff rcsfgmr.l rlabel.inp')
      IF (I /= 0) THEN
        WRITE(*, *)'Error! ' // 'rcsfgmr.l'  // ' and ' // &
                   'rlabel.inp' // ' are different ...'
        STOP &
        'Error! Labeling orbitals are different ...'
      ENDIF
!CYC-2024 End

      print *, 'Reduction based on Dirac-Coulomb (1) or'
      print *, 'Dirac-Coulomb-Breit (2) Hamiltonian?'
      READ(*,*) ICOLBREI
!
      NBLOCK = 0
      CALL FACTT
      CALL SET_CSF_list (NCORE,NPEEL)
      WRITE (6, *) " Block  MR NCSF(L) Before NCSF(G) After NCSF(G)"
      DO
         I_Count = 0
         CALL LODCSL_MR (NCORE,NPEEL,NCFD,NEXT_BLOCK)
         CSF_Number = 0
         open (1312, status = 'scratch', form = 'formatted')
         !open (1312, status = 'unknown', form = 'formatted')
         !WRITE(*, *)"Calling from LODCSL_CSF ......."
         CALL LODCSL_CSF (NCFD,CSF_Number,NCORE,NPEEL,NEXT_CSF)
         !WRITE(*, *)"Exiting from LODCSL_CSF ......."
         !IF(.NOT. NEXT_CSF) EXIT
         !IF(Found_CSF == 1) CYCLE
         IF (CSF_Number == 0 )THEN
           WRITE (6, *)NBLOCK, "   CSF_Number = 0"
           deallocate (Found)
           deallocate (C_shell)
           deallocate (C_quant)
           deallocate (C_coupl)
           deallocate (iqa)
           deallocate (jqsa)
           deallocate (jcupa)
           close(1312)
           IF(.NOT. NEXT_BLOCK) EXIT
           WRITE(22,'(A2)') ' *'
           CYCLE
         ENDIF

         NCF = NCFD + CSF_Number
         ALLOCATE(flagwrite(CSF_Number))
         FLAGWRITE = .FALSE.
         DO I = 1, CSF_Number
           JB = NCFD + I
           DO JA = 1, NCFD
             CALL Interact_CSF(JA,JB,ICOLBREI,int_CSF)
             IF (int_CSF /= 0) THEN
               I_Count = I_Count + 1
               flagwrite(I) = .TRUE.
               EXIT
             ENDIF
           ENDDO
         ENDDO
         WRITE (6,'(3X,I2,6X,I7,3X,I10,3X,I10)')                       &
                NBLOCK,NCFBLK(NBLOCK),NUM_in_BLK_Full(NBLOCK),         &
                I_count+NCFBLK(NBLOCK)
               !NBLOCK,NCFBLK(NBLOCK),CSF_Number,I_count+NCFBLK(NBLOCK)
!
         ! Output
         REWIND(1312)
         ! CSFs in rcsfmr have been written into 22
         DO I = 1, CSF_Number
           read(1312,'(A)')line1
           read(1312,'(A)')line2
           read(1312,'(A)')line3
           if (flagwrite(I)) then
             write(22, '(A)')trim(line1)
             write(22, '(A)')trim(line2)
             write(22, '(A)')trim(line3)
           endif 
         ENDDO

         deallocate (flagwrite)
         deallocate (Found)
         deallocate (C_shell)
         deallocate (C_quant)
         deallocate (C_coupl)
         deallocate (iqa)
         deallocate (jqsa)
         deallocate (jcupa)

         close(1312)
         IF(.NOT. NEXT_BLOCK) EXIT
         WRITE(22,'(A2)') ' *'
      END DO
      CLOSE(22)
      call stoptime (ncount1, 'RCSFGinteract_CSFG')
      STOP
      END PROGRAM RCSFGinteract_CSFG
