      MODULE blandb_I
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
      SUBROUTINE blandb (ORG, NMAX, VARMAX, LOCK, FIL, LOW, LIM, POSN, POSL&
         , MINJ, MAXJ)
      integer, DIMENSION(25,0:10), INTENT(INOUT) :: ORG
      integer, INTENT(IN) :: NMAX
      integer :: VARMAX
      logical, DIMENSION(25,0:10) :: LOCK
      integer, INTENT(IN) :: FIL
      integer, DIMENSION(25,0:10) :: LOW
      integer, DIMENSION(25), INTENT(IN) :: LIM
      integer, DIMENSION(220) :: POSN
!VAST...Dummy argument POSN is not referenced in this routine.
      integer, DIMENSION(220) :: POSL
!VAST...Dummy argument POSL is not referenced in this routine.
      integer :: MINJ
!VAST...Dummy argument MINJ is not referenced in this routine.
      integer :: MAXJ
!VAST...Dummy argument MAXJ is not referenced in this routine.
      END SUBROUTINE
      END INTERFACE
      END MODULE