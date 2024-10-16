!**********************************************************************
!  Written by Chongyang Chen, Fudan university, Shanghai              *
!                                               Nov 2023              *
!**********************************************************************
      MODULE csfg_tv_C
      USE vast_kind_param, ONLY:  DOUBLE, LONG, SHORT, BYTE

! Logical variable LABTVFRST, which is used to record
! T- and V- coefficients, the orbitals of SLATER integrals involving X-
! and Y-potentials.
      LOGICAL :: LABTVFRST
! CYC: 2024/08/14
      LOGICAL :: LREADMCP 
! CYC: 2024/02/13
      LOGICAL  :: LCPOT, LCHM
!======================================================================
! Variables used for construction of potentials, and LAG.
!
! Potentials data involving one-electron integrals
!
! NTPT: Total Number of T-coefficients for Potential construction.
! NDTPT: Size of dynamical arrays TCOEFP, and NCRPST
      INTEGER       ::  NTPT, NDTPT, NATPT 
! Record NTTP for every Jp block
      INTEGER       ::  NBPT(0:50) = 0
! TCOEFP: T-coefficients involving potential caculations.
      REAL(DOUBLE), DIMENSION (:),   pointer :: TCOEFP        ! NDTPT 
! Number of the Column (IC) - Row (IR) pairs using the Same T value to
! calculate the Potentials, this is also the number of the interacting
! orbital-pair IA -- IB using same T coefficients. 
      INTEGER,      DIMENSION (:),   POINTER :: NCRPST        ! NDTPT
! Maximum value of array NCRPST(:)
      INTEGER       :: MAXNCRT
! NCRPT: Number of total Column (IC) - Row (IR) pairs involving Potentials 
! from one-body interaction with T-coefficient.
! NDCRPT: Size of dynamical arrays LABPT, and ICIRPT
      INTEGER       :: NCRPT, NDCRPT, NACRPT
! Record NCRPT for every Jp block
      INTEGER       :: NBCRPT(0:50) = 0
! Detailed IA - IB pairs stored by LABPT
      INTEGER,      DIMENSION (:), POINTER   :: LABPT         ! NDCRPT
! Total number of unique LAB = IA * KEY + IB
      INTEGER  :: NCNTGT
! Sorted arrays with uinque LAB = IA * KEY + IB
      INTEGER,      DIMENSION (:), POINTER   :: LABALLT       ! NCNTGT
! Match arrays LABPT to LABALLT                
      INTEGER,      DIMENSION (:), POINTER   :: MATLABT       ! NDCRPT
      REAL(DOUBLE), DIMENSION(:),  POINTER   :: TSUMG         ! NCOUNTT
! Detailed IC - IR pairs stored by ICIRPT
      INTEGER,      DIMENSION (:,:), POINTER :: ICIRPT        ! NDCRPT

!----------------------------------------------------------------------
! Potentials data involving two-electron integrals
!
! NTPV: Total Number of V-coefficients for Potential construction.
! NDTPV: Size of dynamical arrays VCOEFP, and NCRPSV
      INTEGER       ::  NTPV, NDTPV, NATPV
!
! Many arrasys are deleted, they are too memory-consuming, the needed
! data are reproduced by spinangularXX.f90
!
! Record NTPV for every Jp block
      INTEGER       ::  NBPV(0:50) = 0
! VCOEFP: V-coefficients involving potential caculations.
      REAL(DOUBLE), DIMENSION (:),   pointer :: VCOEFP        ! NDTPV 

      INTEGER :: KMAXV
      INTEGER(BYTE),  DIMENSION (:),  POINTER :: KVCOFP       ! NDTPV
! Number of the Column (IC) - Row (IR) pairs using the Same V value to
! calculate the Potentials, this is also the number of the interacting
! orbitals IA,IB,IC,ID using same V coefficient. 
      INTEGER,      DIMENSION (:),   POINTER :: NCRPSV        ! NDTPV
! Maximum value of array NCRPSV(:)
      INTEGER       :: MAXNCRV = 0

! Total number of V-K-Coefficients generated by one ICG-IRG pair, the
! demension is SUM(NPAIRS2,1:NBLOCK)), these VK terms may use the same
! drs, i.e. the general occupation number.
      INTEGER,      DIMENSION (:),   POINTER :: NVKPG

! NCRPV: Number of total Column (IC) - Row (IR) pairs involving Potentials 
! from TWO-body interaction with V-coefficient.
! NDCRPV: Size of dynamical arrays LABPKV, and ICIRPV
      INTEGER       :: NCRPV, NDCRPV
      INTEGER       :: NACRPV
! Record NCRPV for every Jp block
      INTEGER       :: NBCRPV(0:50) = 0
! Detailed IA,IB,IC,ID, and K pairs stored by LABPV
      INTEGER,      DIMENSION (:), POINTER :: LABPV           ! NDCRPV
! Total number of unique LAB = ..., see below
      INTEGER  :: NCNTGV
! Sorted arrays with uinque LAB = ((IA * KEY + IB) * KEY + IC )* KEY + ID
      INTEGER,      DIMENSION (:), POINTER   :: LABALLV       ! NCNTGV
      REAL(DOUBLE), DIMENSION(:),  POINTER   :: VSUMG         ! NCNTGV
! Match arrays LABPV to LABALLV
      INTEGER,      DIMENSION (:), POINTER   :: MATLABV       ! NDCRPV
      INTEGER,      DIMENSION (:), POINTER   :: MATLABK       ! KMAX
! Detailed IC - IR pairs stored by ICIRPV
      INTEGER,      DIMENSION (:,:), POINTER :: ICIRPV        ! NDCRPV

!======================================================================
!
!======================================================================
! Variables used for construction of H-matrix
! T coefficients needed to construct H matrix
      INTEGER  :: NUMHT, NDHGPT
! LABELs (IA, IB) 
      INTEGER, DIMENSION(:,:), pointer :: LABTH
! T coefficients stored in memory, read from IGMCP31 (csfg_mcpXXX.31)
      REAL(DOUBLE), DIMENSION(:), pointer :: TCOEFH

! V coefficients needed to construct H matrix
      INTEGER  :: NUMHP, NDHGPV, MAXNV = 0 ! ICC - IRR CSFG pairs
      INTEGER  :: NUMHVK, NDHVK           ! Label(1:5, NUMHV)  
! NVCOEF of every ICC - IRR CSFG pair
      INTEGER, DIMENSION(:), pointer :: NVHP
! LABELs (IA, IB, IC, ID, and K) 
      INTEGER(BYTE), DIMENSION(:,:), pointer :: LABVKH
! V coefficients stored in memory, read from IGMCP32 (csfg_mcpXXX.32)
      REAL(DOUBLE), DIMENSION(:), pointer :: VCOEFH
!======================================================================

! IGMCP30: csfg_mcpXXX.30
! IGMCP31: csfg_mcpXXX.31
      INTEGER :: IGMCP30=3030, IGMCP31=3031, IGMCP32=3032

! Number of calls ONESCALAR which generates Tab coefficients.
      INTEGER :: NONESCALAR
! Number of calls of RKCO_GG which generates Vabcd coefficients.
      INTEGER :: NRKCO
! Total number of Tab coefficients
      INTEGER :: LLISTT

! Total number of V(K) coefficients generated
      INTEGER, DIMENSION(:), pointer :: LLISTV

! Number of non-zero CSFG JA (IC) - JB (IR) pairs
      INTEGER :: NPAIRS(1:2, 0:50) = 0
! NODENCOLS ==>:
      INTEGER :: NCOLBLK(0:50) = 0

! Nonzero elements, NELMNTCC ==>
      INTEGER :: NZBLK(0:50) = 0

! Column indexes obtained from csfg_mcpXXX.30
      INTEGER, DIMENSION(:), POINTER :: NODECOLSTV
! IENDC
      INTEGER, DIMENSION(:), POINTER :: IENDCTV
! IROW
      INTEGER, DIMENSION(:), POINTER :: IROWTV

! End index of every interacting CSFG pairs (Column)
! ICGENDC(JA) : End position for COLUMN JA-CSFG, used in
! rangular_mpi_csfg and rmcdhf_mpi_csfg, interacting pairs.
      INTEGER, DIMENSION(:), POINTER :: ICGPEND

! Detailed JB CSFG, used in rangular_mpi_csfg and rmcdhf_mpi_csfg,
! interacting pairs.
      INTEGER, DIMENSION(:), POINTER :: ICGPDET
      INTEGER, DIMENSION(:), POINTER :: IRGPDET
! ncsfDF1, ncsfDF1_core
      INTEGER :: MCPDF1(50), MCPDF1C(50)

      INTEGER :: NCFGPAST, NCSFPAST, NZPAST, NCOLPAST, NROWGPAST
      INTEGER :: NPAIRPAST(2) 

! Column index recorded in setham_csfg
      INTEGER, DIMENSION(:), POINTER :: ICOLIDX


      END MODULE csfg_tv_C