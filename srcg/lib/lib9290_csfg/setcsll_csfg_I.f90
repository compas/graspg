      MODULE setcsll_csfg_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  10:50:34   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      SUBROUTINE setcsll_csfg (NUNIT, NAME, NBLKIN, NBLOCK, NCFBLK, NCFTOT, IDBLK)
      INTEGER, INTENT(IN) :: NUNIT
      CHARACTER (LEN = *), INTENT(INOUT) :: NAME
      INTEGER, INTENT(IN) :: NBLKIN
      INTEGER, INTENT(OUT) :: NBLOCK
      INTEGER, DIMENSION(*), INTENT(INOUT) :: NCFBLK
      INTEGER, INTENT(OUT) :: NCFTOT
      CHARACTER (LEN = 8), DIMENSION(*), INTENT(OUT) :: IDBLK
!VAST...Calls: OPENFL
!...This routine performs I/O.
      END SUBROUTINE
      END INTERFACE
      END MODULE
