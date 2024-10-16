      MODULE csfg_scf_C
      USE vast_kind_param,  ONLY: DOUBLE, LONG
      USE parameter_def,    ONLY: NNNW
!...Created by Pacific-Sierra Research 77to90  4.3E  06:38:40  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      REAL(DOUBLE), DIMENSION(NNNW) :: UCF
      INTEGER, DIMENSION(NNNW) :: METHOD
      REAL(DOUBLE), DIMENSION(NNNW) :: SCNSTY
      REAL(DOUBLE) :: EPSMIN, EPSMAX, EMIN, EMAX, ZINF
      INTEGER :: NDCOF, NXCOF, NYCOF, NDDIM, NXDIM, NYDIM
      REAL(DOUBLE), DIMENSION(:), pointer :: da, xa, ya
      INTEGER, DIMENSION(:), pointer :: nda, nxa, nya

! CYC: save da, xa, ya, and nda, nxa, nya for the .not. LFIX(J), then
! reload them within IMPROVmpi :: COFPOTmpi :: SETCOF_CYC
      INTEGER :: NDIMX
      INTEGER, DIMENSION (:), POINTER :: NDCOFOPT, NXCOFOPT, NYCOFOPT
      REAL(DOUBLE), DIMENSION(:,:), POINTER :: DAOPT, XAOPT, YAOPT
      INTEGER, DIMENSION(:,:), POINTER :: NDAOPT, NXAOPT, NYAOPT
      INTEGER :: JXIPOS(0:NNNW, NNNW) = 0  ! For fast search
      LOGICAL :: FIRSTCOF=.TRUE.
      INTEGER, DIMENSION(:,:,:), POINTER :: NXAIORB  ! For fast search

! CYC 2024, Aug, vectorize NXAOPT, NYAOPT, XAOPT, YAOPT, and the
! corresponding NXCOFOPT and NYCOFOPT
      INTEGER  :: NTOTXA, NTOTYA 
      INTEGER  :: NXCOFW(0:NNNW) = 0, NYCOFW(0:NNNW) = 0
      INTEGER,      DIMENSION (:), POINTER :: NXAWRK, NYAWRK
      REAL(DOUBLE), DIMENSION (:), POINTER :: XAWRK,  YAWRK
      REAL(DOUBLE), DIMENSION (:), POINTER :: TMPXYA
       


! CYC: prepare the data involving one-electron integrals
      INTEGER       ::  NCOUNTT, NDIMT
      INTEGER, DIMENSION (:), POINTER :: LABTALLA
      INTEGER, DIMENSION (:), POINTER :: LABTALLB
      REAL(DOUBLE), DIMENSION(:), POINTER :: TALLSUM

! CYC: prepare the data involving two-electron integrals
      INTEGER       :: NCOUNTV, NDIMV
      INTEGER, DIMENSION (:), POINTER :: LABVALL, KVALL 
      REAL(DOUBLE), DIMENSION(:), POINTER :: VALLSUM

! CYC: Skip the Potential calculation for the orbitals fixed and not
! involved any Lag-Multiply Factor
      LOGICAL, DIMENSION(NNNW) :: LSKIPPOT = .FALSE.

! CYC: 2023/01/10/, calculate and save dsubrr for the xa and ya
! unassociated with MCP coefficients.
      REAL(DOUBLE), DIMENSION(:), POINTER :: DSUBRR
      REAL(DOUBLE), DIMENSION(:,:), POINTER :: DDRS
      !LOGICAL, DIMENSION(:), POINTER :: LDSUBRR
      INTEGER, DIMENSION (:), POINTER :: KNUMIR, IRKPOS
      INTEGER, DIMENSION (:), POINTER :: IRDSUBRR
      REAL(DOUBLE), DIMENSION(:), POINTER :: FC0GC0
! CYC: 2023/02/21, H-matrix calculation
      INTEGER :: NHREC, NDIMHL, NDIMHT
      INTEGER, DIMENSION (:), POINTER :: NUMNH, NHACC, NHIR
      REAL(DOUBLE), DIMENSION(:), POINTER :: FGCH


      END MODULE csfg_scf_C
