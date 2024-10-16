      MODULE brint2_I
      INTERFACE
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!   Modify  by G Gaigalas                         May 2021
      SUBROUTINE BRINT2 (IA,IB,IC,ID,K,TEGRAL,icc,irr,ivv)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER, INTENT(IN) :: ia, ib, ic, id, k
      INTEGER, INTENT(IN) :: ICC, IRR, IVV
      REAL(DOUBLE), INTENT(out) :: tegral
      END SUBROUTINE
      END INTERFACE
      END MODULE
