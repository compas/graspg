!=======================================================================
! IF (LREADMCP) Read data from IGMCP31 and IGMCP32                     *
! And save the data into arrays.                                       *
!                                                                      *
! Read T coefficients into arrays:  IFLAG == 1                         *
! Read V-coefficients into arrays:  IFLAG == 2                         *
!                                                                      *
! IF (.NOT. LREADMCP)   &                                              *
!   Read T and V-coefficients from arrays                              *
!                                                                      *
! These coefficients are  used to calculate H matrixelements.          *
!                    also used to calculate potentials.                *
!                                                                      *
! Written by Chongyang Chen,      Fudan University,          Oct 2023  *
!                                                                      *
!=======================================================================
      SUBROUTINE READ_TV(IFLAG, ICW, IRW, IUNITF, MTYPE, IAA, IBB, TSHELL1)

      USE vast_kind_param, ONLY: DOUBLE, LONG
      USE parameter_def, ONLY: KEYORB
      USE csfg_memory_man
      USE symmatrix_mod, ONLY: CUTOFF, NORBGEN, NSYMCR
      USE buffer_C,      ONLY: NBDIM, NVCOEF, LABEL, COEFF
      USE csfg_tv_C
      USE mpi_C
 
      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)    :: IFLAG, ICW, IRW, IUNITF, MTYPE
      INTEGER                :: IAA, IBB
      REAL(DOUBLE)           :: TSHELL1
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: KEYSQ = KEYORB*KEYORB

      INTEGER            :: ICWTMP, IRWTMP, IUNITFTMP, MTYPETMP, NOSTMP
      INTEGER            :: I, J, NV, IA, IB, IC, ID, NRKCOTMP
      INTEGER            :: ITYPE, NDISC, IOS, ISWAP, IOFF

      INTEGER(BYTE), DIMENSION(:,:), pointer :: LABCDK

!-----------------------------------------------
      SELECT CASE (IFLAG)
!
!=======================================================================
! T-coefficients
!
      CASE (1)
       IF (LREADMCP) THEN
!
! Read T-coefficients into arrays from IGMCP31 (disk file csfg_mcpXXX.31)
!
        READ (IGMCP31, IOSTAT=IOS) NOSTMP, ICWTMP, IRWTMP, IUNITFTMP,  &
                  MTYPETMP, IAA, IBB, TSHELL1
        IF (IOS < 0 .OR. &     ! File End
                               ! or Block End 
           (NOSTMP == 0 .AND. ICWTMP == 0 .AND. IRWTMP == 0)) THEN
          IAA = 0
          IBB = 0
          TSHELL1 = 0.0D0
          BACKSPACE(IGMCP31)
          GOTO 100
        ENDIF

        IF (ICWTMP /= ICW .OR. IRWTMP /= IRW .OR. &
          IUNITFTMP+2000 /= IUNITF .OR. MTYPETMP /= MTYPE) THEN
! There is no one-body contribution betwee ICW-IRW CSFG pair for MTYPE.
! This kind of case occurs for CSFG pairs 2-4, 3-4, 4-4, 4-5 
          IAA = 0
          IBB = 0
          TSHELL1 = 0.0D0
          BACKSPACE(IGMCP31)
          GOTO 100
        ELSE
!
! Read successfully T coefficients
!
! Debug
          NONESCALAR = NONESCALAR + 1
          !WRITE(IUNITF,'(I8, 2X,4HICW=,I8, 2X,4HIRW=,I8,      &
          !            2X,5HTYPE=,I4, 2X,6HMTYPE=,I4,          &
          !            2X,3HIA=,I4, 2X,3HIB=,I4, 2X,1PE12.4)') &
          !     NONESCALAR, ICW, IRW, IUNITF-2000,             &
          !     MTYPE, IAA, IBB, TSHELL1
! Check
          IF (NOSTMP /= NONESCALAR) THEN
            WRITE (*,*) "MYID, ICW, IRW, NOSTMP, NONESCALAR =",     &
                        MYID, ICW, IRW, NOSTMP, NONESCALAR
            STOP "Error!  NOSTMP /= NONESCALAR ..."
          ENDIF
          GOTO 100
        ENDIF

! Record data into arrays TCOEFG, and LABTGP, then return
100     CONTINUE
        NUMHT = NUMHT + 1
        IF (NUMHT > NDHGPT) THEN
          NDHGPT = NDHGPT + NDHGPT
          CALL RALLOC (LABTH , 2, NDHGPT, 'LABTH ', 'READ_TV')
          CALL RALLOC (TCOEFH, NDHGPT, 'TCOEFH', 'READ_TV')
        ENDIF
        IF (IAA > IBB) THEN
          ISWAP = IAA
          IAA = IBB
          IBB = ISWAP
        ENDIF
        LABTH(1, NUMHT) = IAA
        LABTH(2, NUMHT) = IBB
        TCOEFH(NUMHT) = TSHELL1
        RETURN
!
! Read T-coefficients from array
!
       ELSE
        NUMHT = NUMHT + 1
        !IF (NUMHT > NDHGPT) STOP "Error in read_TV, NUMHT > NDHGPT ..."
        IAA = LABTH(1, NUMHT)
        IBB = LABTH(2, NUMHT)
        TSHELL1 = TCOEFH(NUMHT)

       ENDIF
!
!=======================================================================
!
! V-coefficients
!
      CASE (2)
!
! There must be V-coefficients for the recorded ICW-IRW CSFG pair
!
       IF (LREADMCP) THEN
!
! Read V-coefficients into arrays from IGMCP32 (disk file csfg_mcpXXX.32)
!
        READ (IGMCP32, IOSTAT=IOS) NRKCOTMP, ICWTMP, IRWTMP, IUNITFTMP,&
                                   MTYPETMP, NV
        IF (IOS /= 0) THEN
          WRITE(*,*)'IFLAG, ICW, IRW, IUNITF, MTYPE =', &
                    IFLAG, ICW, IRW, IUNITF, MTYPE
          STOP "Error when reading V-coefficient ..." 
        ENDIF

        IF (ICWTMP /= ICW .OR. IRWTMP /= IRW .OR.             &
            IUNITFTMP+2000 /= IUNITF .OR. MTYPETMP /= MTYPE) THEN
! There is no two-body contribution betwee ICW-IRW CSFG pair for MTYPE.
! This kind of case occurs for CSFG pairs 2-4, 3-4, 4-4, 4-5 
           NVCOEF = 0
           BACKSPACE(IGMCP32)
           GOTO 200
        ENDIF
!
! Read successfully V-coefficients for CSFG pair ICW - IRW with MTYPE
!
        NRKCO = NRKCO + 1
        IF (NV > MAXNV) MAXNV = NV
        IF (NV > NBDIM) THEN
          NBDIM = NV
          CALL ALCBUF(2)
        ENDIF

        CALL ALLOC(LABCDK, 5, NV, 'LABCDK', 'READ_TV')
        NVCOEF = NV
         
        IF (NV > 0) THEN
! LABEL is INTEGER*4 in buffer_C.f90
          ! READ (IGMCP32) (LABEL(1:5,I),I = 1, NV), COEFF(1:NV)
          READ (IGMCP32) (LABCDK(1:5,I),I = 1, NV), COEFF(1:NV)
          LABEL(1:5,1:NV) = LABCDK(1:5,1:NV)
        ENDIF
! Debug
        !WRITE(IUNITF, '(I8, 2X,4HICW=,I8, 2X,4HIRW=,I8, &
        !               2X,5HTYPE=,I4, 2X,6HMTYPE=,I4, 2X,3HNV=,I4)')&
        !               NRKCO, ICW, IRW, IUNITF-2000, MTYPE, NV
        !write(IUNITF,'(4x,7hIVCOEF=,1000i11)')(J, J=1,NV)
        !write(IUNITF,'(4x,7hLabel1=,1000i11)')(LABEL(1,J), J=1,NV)
        !write(IUNITF,'(4x,7hLabel2=,1000i11)')(LABEL(2,J), J=1,NV)
        !write(IUNITF,'(4x,7hLabel3=,1000i11)')(LABEL(3,J), J=1,NV)
        !write(IUNITF,'(4x,7hLabel4=,1000i11)')(LABEL(4,J), J=1,NV)
        !write(IUNITF,'(4x,7hLabel5=,1000i11)')(LABEL(5,J), J=1,NV)
        !write(IUNITF,'(4x,7hCoeff =,1000(1pe11.3))')(Coeff(J), J=1,NV)
! Check
        IF (NRKCO /= NRKCOTMP) THEN
          WRITE(*,*)"ICW, IRW, NRKCO, NRKCOTMP =",            &
                    ICW, IRW, NRKCO,NRKCOTMP
          STOP "Error when reading V-coefficients: NRKCO /= NRKCOTMP!"
        ENDIF

! Record V-coefficients data into arrays, then return
200     CONTINUE
        NUMHP = NUMHP + 1
        IF (NUMHP > NDHGPV) THEN
          NDHGPV = NDHGPV + NDHGPV
          CALL RALLOC (NVHP,   NDHGPV, 'NVHP',   'READ_TV')
        ENDIF
        NV = NVCOEF
        NVHP(NUMHP) = NV

        IOFF = NUMHVK
        NUMHVK = NUMHVK + NV
        IF (NUMHVK > NDHVK) THEN
          NDHVK = NDHVK + NDHVK
          CALL RALLOC (LABVKH, 5, NDHVK, 'LABVKH',   'READ_TV')
          CALL RALLOC (VCOEFH, NDHVK, 'VCOEFH',   'READ_TV')
        ENDIF
        DO J = 1, NV
          !LABVKH(1:5, IOFF + J) = LABEL(1:5, J)
          LABVKH(1:5, IOFF + J) = LABCDK(1:5, J)
          VCOEFH(IOFF + J) = COEFF(J)
        ENDDO

        IF (NV > 0) CALL DALLOC(LABCDK, 'LABCDK', 'READ_TV')

        RETURN

       ELSE
! Read V-coefficients data from arrays, then return
        NUMHP = NUMHP + 1
        !IF (NUMHP > NDHGPV) STOP "Error in read_TV, NUMHP > NDHGPV ..."
      
        NVCOEF = NVHP(NUMHP)
        NV = NVCOEF
        !IF (NV > MAXNV) STOP "Error in read_TV, NV > MAXNV ..."

        IOFF = NUMHVK
        DO J = 1, NVCOEF
          LABEL(1:5, J) = LABVKH(1:5, IOFF + J)
          COEFF(J) = VCOEFH(IOFF + J)
        ENDDO
        NUMHVK = NUMHVK + NVCOEF
        !IF (NUMHVK > NDHVK) STOP "Error in read_TV, NUMHVK > NDHVK ..."
        RETURN
       ENDIF
!
!=======================================================================
      CASE DEFAULT
        STOP "Error IFLAG value in read_TV.f90 ..."
      END SELECT

      RETURN
      END SUBROUTINE READ_TV
