!=======================================================================
! Save T coefficients into IUNITF: IFLAG == 1                          *
!      MTYPE = 1: mcpXXX_csfg.31                                       *
!      MTYPE = 2: mcpXXX_csfg.31                                       * 
! Save buffer_C into IUNITF:       IFLAG == 2                          *
!      MTYPE = 1: mcpXXX_csfg.32                                       *
!      MTYPE = 2: mcpXXX_csfg.32                                       *
!      MTYPE = 3: mcpXXX_csfg.32                                       *
! These coefficients are used to calculate H matrixelements.           *
! They are also used to build the X- and Y-potentials,  and            *
! calculate the Lagrange multipliers.                                  *
!                                                                      *
! Written by Chongyang Chen,      Fudan University,          Oct 2023  *
!                                                                      *
!=======================================================================
      SUBROUTINE SAVETV(IFLAG, ICW, IRW, IUNITF, MTYPE, IAA, IBB, TSHELL1)

      USE vast_kind_param, ONLY: DOUBLE, LONG, BYTE
      USE parameter_def, ONLY: KEYORB
      USE symmatrix_mod, ONLY: CUTOFF, NORBGEN, NSYMCR
      USE buffer_C,      ONLY: NVCOEF, LABEL, COEFF
      Use csfg_tv_C
      Use csfg_memory_man

      USE labpack_I

      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER    :: IFLAG, ICW, IRW, IUNITF, MTYPE
      INTEGER    :: IAA, IBB
      REAL*8     :: TSHELL1
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: KEYSQ = KEYORB*KEYORB

      INTEGER            :: I, J, NV
      INTEGER            :: ITYPE, NDISC, LPTMP

      INTEGER(BYTE)      :: IA, IB, IC, ID, LABTMP(5)
      INTEGER(BYTE), DIMENSION(:,:), pointer :: LABCDK

      INTEGER, DIMENSION(:),   pointer :: LABP
      REAL(DOUBLE)   :: VTMP
      REAL(DOUBLE),  DIMENSION(:),   pointer :: VCOEFF

!-----------------------------------------------
      SELECT CASE (IFLAG)

! T-coefficients
      CASE (1)
        IF (IAA == 0 .OR. ABS(TSHELL1) .LT. CUTOFF) RETURN

! Write to mcpXXX_csfg.31
        NONESCALAR = NONESCALAR + 1
        !IA = IAA; IB = IBB
        WRITE(31) NONESCALAR, ICW, IRW, IUNITF-2000,          &
                  MTYPE, IAA, IBB, TSHELL1
! Debug
        !WRITE(IUNITF,'(I8, 2X,4HICW=,I8, 2X,4HIRW=,I8,        &
        !              2X,5HTYPE=,I4, 2X,6HMTYPE=,I4,          &
        !              2X,3HIA=,I4, 2X,3HIB=,I4, 2X,1PE12.4)') &
        !       NONESCALAR, ICW, IRW, IUNITF-2000,             &
        !       MTYPE, IAA, IBB, TSHELL1
        RETURN 

! V-coefficients
      CASE (2)
        IF (NVCOEF == 0) RETURN

        NV = COUNT(ABS(COEFF(1:NVCOEF)) > CUTOFF)
        IF (NV == 0) RETURN

        NRKCO = NRKCO + 1
        !ALLOCATE(LABCDK(5,NV))
        !ALLOCATE(VCOEFF(NV))
        !ALLOCATE(LABP(NV))
        CALL ALLOC(LABCDK, 5, NV, 'LABCDK', 'SAVETV')
        CALL ALLOC(VCOEFF,    NV, 'VCOEFF', 'SAVETV')
        CALL ALLOC(LABP,      NV, 'LABP',   'SAVETV')

        ITYPE = IUNITF - 2000
     
        IF (MTYPE == 1) THEN
          NV = 0
          DO I = 1, NVCOEF
            IF (ABS(COEFF(I)) > CUTOFF) THEN
              LLISTV(LABEL(5,I)) = LLISTV(LABEL(5,I)) + 1
              NV = NV + 1
              VCOEFF(NV) = COEFF(I)
              LABCDK(1:5,NV) = LABEL(1:5,I)

              LABP(NV) = LABPACK(LABEL(1,I), LABEL(2,I), &
                                 LABEL(3,I), LABEL(4,I))
            ENDIF
          ENDDO

          ! Sorted Vk according to the original GRASP2018 orders
          IF (ICW /= IRW) THEN
            CALL SORTOFFDV(NV, VCOEFF, LABCDK, LABP)
          ELSE
            CALL SORTDIAGV(NV, VCOEFF, LABCDK, LABP)  
            !CALL SORTOFFDV(NV, VCOEFF, LABCDK, LABP)
          ENDIF

          WRITE(32) NRKCO, ICW, IRW, IUNITF-2000, MTYPE, NV
          WRITE(32) (LABCDK(1:5,I),I = 1, NV), VCOEFF(1:NV)
        ELSE
          ! MTYPE /= 1 for the following cases
          SELECT CASE (ITYPE)
          CASE (24)
            NDISC = 3
          CASE (34)
            NDISC = 4
          CASE (44)
            NDISC = 4
          CASE (4)
            NDISC = 4
          CASE (45)
            NDISC = 4
          CASE DEFAULT
            WRITE(*,*)"ITYPE Error in savetv.f90, ITYPE =", ITYPE
            STOP 
          END SELECT

          NV = 0 
          DO I = 1, NVCOEF
! Discard the coefficients which could be obtained by MTYPE = 1
            NSYMCR = COUNT(LABEL(1:4,I) > NORBGEN)
            IF (NSYMCR == NDISC) CYCLE

            IF (ABS(COEFF(I)) > CUTOFF) THEN
              LLISTV(LABEL(5,I)) = LLISTV(LABEL(5,I)) + 1
              NV = NV + 1
              VCOEFF(NV) = COEFF(I)
              LABCDK(1:5,NV) = LABEL(1:5,I)

              LABP(NV) = LABPACK(LABEL(1,I), LABEL(2,I), &
                                 LABEL(3,I), LABEL(4,I))
            ENDIF
          ENDDO 

          IF (NV > 0) CALL SORTOFFDV(NV, VCOEFF, LABCDK, LABP)

! Output
          WRITE(32) NRKCO, ICW, IRW, IUNITF-2000, MTYPE, NV
          IF (NV > 0) THEN
            WRITE(32) (LABCDK(1:5,I),I = 1, NV), VCOEFF(1:NV)
          ENDIF
        ENDIF 
        
! Debug
!        WRITE(IUNITF, '(I8, 2X,4HICW=,I8, 2X,4HIRW=,I8, &
!                       2X,5HTYPE=,I4, 2X,6HMTYPE=,I4, 2X,3HNV=,I4)')&
!                       NRKCO, ICW, IRW, IUNITF-2000, MTYPE, NV
!        write(IUNITF,'(4x,7hIVCOEF=,1000i11)')(J, J=1,NV)
!        write(IUNITF,'(4x,7hLabel1=,1000i11)')(LABCDK(1,J), J=1,NV)
!        write(IUNITF,'(4x,7hLabel2=,1000i11)')(LABCDK(2,J), J=1,NV)
!        write(IUNITF,'(4x,7hLabel3=,1000i11)')(LABCDK(3,J), J=1,NV)
!        write(IUNITF,'(4x,7hLabel4=,1000i11)')(LABCDK(4,J), J=1,NV)
!        write(IUNITF,'(4x,7hLabel5=,1000i11)')(LABCDK(5,J), J=1,NV)
!        write(IUNITF,'(4x,7hCoeff =,1000(1pe11.3))')(VCoeff(J), J=1,NV)
!
!        WRITE(9032, '(I8, 2X,4HICW=,I8, 2X,4HIRW=,I8, &
!                       2X,5HTYPE=,I4, 2X,6HMTYPE=,I4, 2X,3HNV=,I4)')&
!                       NRKCO, ICW, IRW, IUNITF-2000, MTYPE, NV
!        write(9032,'(4x,7hIVCOEF=,1000i11)')(J, J=1,NV)
!        write(9032,'(4x,7hLabel1=,1000i11)')(LABCDK(1,J), J=1,NV)
!        write(9032,'(4x,7hLabel2=,1000i11)')(LABCDK(2,J), J=1,NV)
!        write(9032,'(4x,7hLabel3=,1000i11)')(LABCDK(3,J), J=1,NV)
!        write(9032,'(4x,7hLabel4=,1000i11)')(LABCDK(4,J), J=1,NV)
!        write(9032,'(4x,7hLabel5=,1000i11)')(LABCDK(5,J), J=1,NV)
!        write(9032,'(4x,7hCoeff =,1000(1pe11.3))')(VCoeff(J), J=1,NV)
 
        CALL DALLOC(LABCDK, 'LABCDK', 'SAVETV')
        CALL DALLOC(VCOEFF, 'VCOEFF', 'SAVETV')
        CALL DALLOC(LABP,   'LABP',   'SAVETV')

        RETURN

      CASE DEFAULT
        WRITE (*,*)"Error IFLAG value in savetv.f90 ..." 
        STOP
      END SELECT

      RETURN
      END SUBROUTINE SAVETV
