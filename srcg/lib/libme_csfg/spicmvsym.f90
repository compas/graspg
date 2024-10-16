!***********************************************************************
!                                                                      *
      SUBROUTINE SPICMVsym (N,M,B,C)
!                                                                      *
!   This routine now works for both rscfmpivu and rcimpivu. By removing*
!   the include statement and the call to gdsummpi and setting myid=0, *
!   nprocs=1, it also works for the corresponding serial program.      *
!   Matrix is stored in the mode of upper-triangle-by columns, or      *
!   you can say lower-triangle-by-rows. 98-08-06                       *
!                                                                      *
!   Matrix-matrix product: C = AB.  A  sparse  representation of the   *
!   lower triangle of the  (NxN)  matrix  A  is assumed available in   *
!   COMMON/HMAT/.                                                      *
!                                                                      *
!   This is an adaptation of  Andreas Stathopulos   routine  SPSBMV,   *
!   and is specific to GRASP2 derivatives.                             *
!                                                                      *
!   Call(s) to: [AUXBLAS]: DINIT/SINIT                                 *
!                                                                      *
!   F A Parpia and A Stathopoulos         Last revision: 19 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 29 Jul 1998   *
!                                                                      *
!   Modified from SPICMVmpi for symmetry-ordered-CSFs (CSFG)           *
!   Chongyang Chen, Fudan university, 202008                           *
!   Last modification by CYC     Dec  2023                             *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE, LONG
      use symmatrix_mod
      USE hmat_C
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: M
      REAL(DOUBLE), DIMENSION(N,M), INTENT(IN) :: B
      REAL(DOUBLE), DIMENSION(N,M), INTENT(INOUT) :: C
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: DIAG, DL
      INTEGER      :: IEND, IBEG
      INTEGER      :: I, ICOL, IV, NELC
!-----------------------------------------------
!
!   Initialise the result matrix; note that this is specific to the
!   data structure of DVDSON --- no overdimensioning
!
      CALL DINIT (N*M, 0.D0, C, 1)

      ibeg = 1
      !DO ICOL = myid + 1, N, nprocs
       DO I = 1, NODENCOLS
         ICOL = NODECOLS(I) 
            !IBEG = IENDC(ICOL-1)+1
            !IEND = IENDC(ICOL)
         IEND = IENDC(ICOL)
         NELC = IEND - IBEG + 1
         DO IV = 1, M
            DIAG =  C(ICOL,IV) + EMT(IEND)*B(ICOL,IV)
            CALL DMERGE (NELC-1,B(1,IV),C(1,IV),                     &
                         IROW(IBEG),EMT(IBEG),B(ICOL,IV),DL)
            C(ICOL,IV) = DIAG + DL
         ENDDO
         ibeg = iend + 1
      ENDDO

      CALL gdsummpi (C, N*M)

      RETURN
      END SUBROUTINE SPICMVsym 

!C
!C fsplited from spicmv.f
!C
!      SUBROUTINE dmerge (n, db, dc, idy, da, dconst, dl)
!C 
!C  this merge version has the advantage of loading da(i)
!C  and idy(i) only once.
!C
!      IMPLICIT REAL*8           (A-H,O-Z)
!      DIMENSION da(n), db(*), dc(*), idy(n)
!
!      dsum = 0.0
!      DO i = 1, n
!         dsum = dsum + da(i) * db(idy(i))
!         dc(idy(i)) = dc(idy(i)) + dconst * da(i)
!      ENDDO
!      dl = dsum
!
!      RETURN
!      END
