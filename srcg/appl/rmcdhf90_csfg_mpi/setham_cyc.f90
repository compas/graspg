!***********************************************************************
      SUBROUTINE setham_cyc (ICC, IRR, IC, IR)
!***********************************************************************
!                                                                      *
!   Modified for CSFG list                                             *
!   Chongyang Chen,  Fudan University, Shanghai,           Oct  2023   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,   ONLY:  NNNW, KEYORB
      USE csfg_memory_man
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I
!-----------------------------------------------
!   C O M M O N  B L O C K S
!-----------------------------------------------
      USE BUFFER_C,  ONLY: NVCOEF, LABEL, COEFF
      USE iccu_C,    ONLY: ICCUT
!CYC
      Use csfg_decide_C, ONLY : LSYMD
      Use symexpand_mod, ONLY : TotCSFs_perblock, ncsfDF1
      Use symmatrix_mod, KMAXTmp=>KMAX
      Use symmatrix_restart_C, ONLY : ncsfDF1_core 
      Use csfg_tv_C

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ICC, IRR, IC, IR
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: KEYSQ = KEYORB*KEYORB
      REAL(DOUBLE), PARAMETER :: CUTOFF0 = 1.0D-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      REAL(DOUBLE) :: VCOEFF

      INTEGER :: JA, JB, NPOS, JAA, JBB, JCSFA, JCSFB
      INTEGER :: NELMNTCC, IUNITF, NUMOS, NUMRKCO
      INTEGER :: NDIM, IPOS, NUMP, NCOLG, IOS
      
!=======================================================================
!
! CSFG list
      JA = ICC  ! Block index
      JB = IRR  ! Block index
      JAA = NCFGPAST + ICC  ! Serial number in the whole CSFG list
      JBB = NCFGPAST + IRR  ! Serial number in the whole CSFG list
! CSFs have the QA, and the coupling scheme, using CSFG list
      JCSFA = JAA ! Still CSFG
      JCSFB = JBB ! Still CSFG  
! Normal CSF list if used
      !JCSFA = NCSFPAST + IC
      !JCSFB = NCSFPAST + IR

      LTRANSPOSE = .FALSE.
      LTRANSFER = .TRUE.
      EMTBLOCK = 0.D0

      IF (NTYPE(1,JAA) == 1 .AND. NTYPE(1,JBB) == 1) THEN
        IF (JA == JB) THEN
          IUNITF = 2001
          CALL SPINANGULAR1(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
        ELSE
          IUNITF = 2011
          CALL SPINANGULAR11(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
        ENDIF

      ELSEIF (NTYPE(1,JAA) == 2 .AND. NTYPE(1,JBB) == 1) THEN
        IUNITF = 2012
        CALL SPINANGULAR12(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
      ELSEIF (NTYPE(1,JAA) == 1 .AND. NTYPE(1,JBB) == 2) THEN
        IUNITF = 2012
        LTRANSPOSE = .TRUE.
        CALL SPINANGULAR12(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)

      ELSEIF (NTYPE(1,JAA) == 3 .AND. NTYPE(1,JBB) == 1) THEN
        IUNITF = 2013
        CALL SPINANGULAR13(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
      ELSEIF (NTYPE(1,JAA) == 1 .AND. NTYPE(1,JBB) == 3) THEN
        IUNITF = 2013
        LTRANSPOSE = .TRUE.
        CALL SPINANGULAR13(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)

      ELSEIF (NTYPE(1,JAA) == 4 .AND. NTYPE(1,JBB) == 1) THEN
        IUNITF = 2014
        CALL SPINANGULAR14(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
      ELSEIF (NTYPE(1,JAA) == 1 .AND. NTYPE(1,JBB) == 4) THEN
        IUNITF = 2014
        LTRANSPOSE = .TRUE.
        CALL SPINANGULAR14(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)

      ELSEIF (NTYPE(1,JAA) == 5 .AND. NTYPE(1,JBB) == 1) THEN
        IUNITF = 2015
        CALL SPINANGULAR15(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
      ELSEIF (NTYPE(1,JAA) == 1 .AND. NTYPE(1,JBB) == 5) THEN 
        IUNITF = 2015
        LTRANSPOSE = .TRUE.
        CALL SPINANGULAR15(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)

      ELSEIF (NTYPE(1,JAA) == 2 .AND. NTYPE(1,JBB) == 2) THEN
        IF (JA == JB) THEN
          IUNITF = 2002
          CALL SPINANGULAR2(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
        ELSE
          IUNITF = 2022
          CALL SPINANGULAR22(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
        ENDIF 

      ELSEIF (NTYPE(1,JAA) == 3 .AND. NTYPE(1,JBB) == 2) THEN
        IUNITF = 2023
        CALL SPINANGULAR23(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
      ELSEIF (NTYPE(1,JAA) == 2 .AND. NTYPE(1,JBB) == 3) THEN
        IUNITF = 2023
        LTRANSPOSE = .TRUE.
        CALL SPINANGULAR23(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)

      ELSEIF (NTYPE(1,JAA) == 4 .AND. NTYPE(1,JBB) == 2) THEN
        IUNITF = 2024
        CALL SPINANGULAR24(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
      ELSEIF (NTYPE(1,JAA) == 2 .AND. NTYPE(1,JBB) == 4) THEN
        IUNITF = 2024
        LTRANSPOSE = .TRUE.
        CALL SPINANGULAR24(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)

      ELSEIF (NTYPE(1,JAA) == 5 .AND. NTYPE(1,JBB) == 2) THEN
        IUNITF = 2025
        CALL SPINANGULAR25(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
      ELSEIF (NTYPE(1,JAA) == 2 .AND. NTYPE(1,JBB) == 5) THEN
        IUNITF = 2025
        LTRANSPOSE = .TRUE.
        CALL SPINANGULAR25(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)

      ELSEIF (NTYPE(1,JAA) == 3 .AND. NTYPE(1,JBB) == 3) THEN
        IF (JA == JB) THEN
          IUNITF = 2003
          CALL SPINANGULAR3(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
        ELSE
          IUNITF = 2033
          CALL SPINANGULAR33(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
        ENDIF

      ELSEIF (NTYPE(1,JAA) == 4 .AND. NTYPE(1,JBB) == 3) THEN
        IUNITF = 2034
        CALL SPINANGULAR34(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
      ELSEIF (NTYPE(1,JAA) == 3 .AND. NTYPE(1,JBB) == 4) THEN 
        IUNITF = 2034
        LTRANSPOSE = .TRUE.
        CALL SPINANGULAR34(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)

      ELSEIF (NTYPE(1,JAA) == 5 .AND. NTYPE(1,JBB) == 3) THEN
        IUNITF = 2035
        CALL SPINANGULAR35(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
      ELSEIF (NTYPE(1,JAA) == 3 .AND. NTYPE(1,JBB) == 5) THEN 
        IUNITF = 2035
        LTRANSPOSE = .TRUE.
        CALL SPINANGULAR35(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)

      ELSEIF (NTYPE(1,JAA) == 4 .AND. NTYPE(1,JBB) == 4) THEN
        IF (JA == JB) THEN
          IUNITF = 2004
          CALL SPINANGULAR4(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
        ELSE
          IUNITF = 2044
          CALL SPINANGULAR44(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
        ENDIF

      ELSEIF (NTYPE(1,JAA) == 5 .AND. NTYPE(1,JBB) == 4) THEN
        IUNITF = 2045
        CALL SPINANGULAR45(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
      ELSEIF (NTYPE(1,JAA) == 4 .AND. NTYPE(1,JBB) == 5) THEN 
        IUNITF = 2045
        LTRANSPOSE = .TRUE.
        CALL SPINANGULAR45(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)

      ELSEIF (NTYPE(1,JAA) == 5 .AND. NTYPE(1,JBB) == 5) THEN
        IF (JA == JB) THEN
          IUNITF = 2005
          CALL SPINANGULAR5(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
        ELSE
          IUNITF = 2055
          CALL SPINANGULAR55(JA, JB, IUNITF, NPOS, JCSFA, JCSFB)
        ENDIF
      ELSE
        STOP "Sth Error within NTYPE(1,:) in mcpmpi_csfg ..." 
      ENDIF

      RETURN
      END SUBROUTINE setham_cyc
