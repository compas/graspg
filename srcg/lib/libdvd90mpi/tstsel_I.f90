      MODULE tstsel_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      LOGICAL FUNCTION tstsel (NLOOPS, KPASS, NUME, NEIG, ISELEC, SVEC, EIGVAL, ICV&
         , CRITE, CRITC, ROWLAST, IND, OLDVAL, NNCV, INCV)
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER, INTENT(IN) :: NLOOPS
      INTEGER, INTENT(IN) :: KPASS
      INTEGER, INTENT(IN) :: NUME
      INTEGER, INTENT(IN) :: NEIG
      INTEGER, DIMENSION(NEIG), INTENT(IN) :: ISELEC
      REAL(DOUBLE), DIMENSION(KPASS*NUME), INTENT(IN) :: SVEC
      REAL(DOUBLE), DIMENSION(NUME), INTENT(IN) :: EIGVAL
      INTEGER, DIMENSION(NUME), INTENT(OUT) :: ICV
      REAL(DOUBLE), INTENT(IN) :: CRITE
      REAL(DOUBLE), INTENT(IN) :: CRITC
      REAL(DOUBLE), DIMENSION(NEIG), INTENT(INOUT) :: ROWLAST
      INTEGER, DIMENSION(NEIG), INTENT(INOUT) :: IND
      REAL(DOUBLE), DIMENSION(NUME), INTENT(IN) :: OLDVAL
      INTEGER, INTENT(INOUT) :: NNCV
      INTEGER, DIMENSION(NEIG), INTENT(OUT) :: INCV
!VAST.../MPI/ MYID(IN)
!VAST...Calls: IDAMAX
!...This routine performs I/O.
      END FUNCTION
      END INTERFACE
      END MODULE
