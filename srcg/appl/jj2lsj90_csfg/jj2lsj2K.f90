!
!***********************************************************************
!                                                                      *
      PROGRAM jj2lsj2K
!                                                                      *
!     This MAIN program controls the transformation of atomic states,  *
!     which are given in a jj-coupled CSF basis, into an LS-coupled    *
!     basis. The program requires a jj-coupled basis in standard order *
!     where, if both subshells | n j = l-1/2> and | n j = l+1/2>       *
!     of a given shell (nl) occurs, they always follow successively    *
!     in this order. The LS-coupled basis, moreover, is given in       *
!     the same sequence of shells.                                     *
!                                                                      *
!     All LS-jj transformation coefficients are precalculated and      *
!     'stored' in the modules rabs_lsj_data_1, rabs_lsj_data_2 and     *
!     rabs_lsj_data_3.                                                 *
!                                                                      *
!     Calls:  FACTT, SETISO, JJ2LSJ, starttime, stoptime.              *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                                  May 2011   *
!     VILNIUS                                               May 2017   *
!                                                                      *
!     Modified by G. Gaigalas and C.Y. Chen                     2021   *
!     Modified by G. Gaigalas                                   2022   *
!     Last modification for CSFG,                                      *
!          and fixes to some minor bugs by CYC,             Dec 2023   *
!     Modified by G. Gaigalas                               Jan 2024   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE jj2lsj_code
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer           :: ncount1
!-----------------------------------------------
      call starttime (ncount1, 'jj2lsj_csfg')
      print *, " "
      print *,                                                         &
      "jj2lsj_csfg: Transformation of ASFs from a jj-coupled CSF basis"
      print *,                                                         &
      "             into an LS-coupled CSF basis  (Fortran 95 version)"
      print *,                                                         &
      "             (C) Copyright by   G. Gaigalas and Ch. F. Fischer,"
      print *, "             (2024)."
      print *,                                                         &
      "             Modifications using CSFGs      by C.Y. Chen (2023)"
      print *, "Input files: isodata, name.c, name.(c)m, name.g, name.l"
      print *, "          (optional)  name.lsj.T"
      print *, "         Ouput files: name.lsj.lbl,"
      print *, "          (optional)  name.lsj.c, name.lsj.j,"
      print *,                                                         &
      "                   name.uni.lsj.lbl, name.uni.lsj.sum,"
      print *, "          name.lsj.T"
      print *,                                                         &
      "             Optional file name.lsj.T is not available for"
      print *,                                                         &
      "             jj2lsj_csfg"
      print *, " "
!
!  Set up the table of logarithms of factorials
      call factt
!
      CALL setiso('isodata')
      CALL jj2lsj
!
      call stoptime (ncount1, 'jj2lsj_csfg')
      stop
      end program jj2lsj2K
