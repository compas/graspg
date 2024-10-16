      MODULE brint1_I
      INTERFACE
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!   Modify  by G Gaigalas                         May 2021
      SUBROUTINE BRINT1 (IA,IB,IC,ID,K,TEGRAL,ICC,IRR,IVV)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER, INTENT(IN) :: ia, ib, ic, id, k
      INTEGER, INTENT(IN) :: ICC, IRR, IVV
      REAL(DOUBLE), INTENT(out) :: tegral
      END SUBROUTINE
      END INTERFACE
      END MODULE
