      MODULE rkintc_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!   Modify  by G Gaigalas                         May 2021
      SUBROUTINE rkintc (IA0, IB0, IC0, ID0, K, TEGRAL)
      USE vast_kind_param,ONLY: DOUBLE
!      INTEGER KMAX
!      PARAMETER (KMAX = 20)
!CYC      INTEGER, INTENT(INOUT) :: IA0, IB0, IC0, ID0
      INTEGER, INTENT(IN) :: IA0, IB0, IC0, ID0
      INTEGER, INTENT(IN) :: K
      REAL(DOUBLE), INTENT(OUT) :: TEGRAL
      END SUBROUTINE
      END INTERFACE
      END MODULE
