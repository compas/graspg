      SUBROUTINE DVDSON(IRC, IREV, N, LIM, NOC, ILOW, IHIGH, ISELEC, NIV, &
         MBLOCK, CRITE, CRITC, CRITR, MAXITER, WORK, IWRSZ, IWORK, IIWSZ, HIEND&
         , NLOOPS, IERR)
!=======================================================================
!
!     Author: Andreas Stathopoulos, Charlotte F. Fischer
!
!     Computer Science Department
!     Vanderbilt University
!     Nashville, TN 37212
!     andreas@vuse.vanderbilt.edu
!     cff@vuse.vanderbilt.edu
!
!     June 1995
!
!     Copyright (c) by Andreas Stathopoulos and Charlotte F. Fischer
!
!     References: A Stathopoulos and C F Fischer,
!          Comp. Phys. Commun. 79 (1994) 268-290.
!
!                   A Stathopoulos, Y. Saad and C F Fischer,
!                   J. Computational and applied Mathematics (in print).
! 
!                                                                      *
!     Modify to use MPI by Chongyang Chen, Fudan unversity, 2020       *
!     Also forcing dvdson to run at least TWO NLOOPS                   *
!                                                                      *
!     Last modification by Chongyang Chen               Aug 2024       *
!                                                                      *
!CYC-2024                                                              *
! Forcing at least TWO iterations performed within DVDSON              *
!                                                                      *
! Reason: WORK array is not properly re-initialized, or the            *
! calculations converge coincidentally by one iteration in the         *
! following case:                                                      *
!                                                                      *
! More eigen-pairs are searched from a smaller CSF expansion of Jp     *
! block indexed by IBLOCK, but less eigen-pairs are serached from a    *
! larger CSF expansion for another Jp block indexed by IBLOCK+1.       *
!                                                                      *
! For example, RCI calculation performed for all the 3l3l' states in Fe*
! XV, 5 eigen-pairs of 2- symmetry are searched from a expansion of    *
! 761882 CSFs, but only 2 eigen-pairs of 3+ symmetry are searched from *
! a expansion of 1010283 CSFs. The expansions are obtained by SD       *
! excitation from all 3l3l' confirgurations, considering all possible  *
! SD substitutions from 2s, 2p, 3l eletrons (VV, CV, CC).              *
!                                                                      *
! This issue could also exist in rmcdhf calculation. It dersves to be  *
! fiexed in GRASP2018.                                                 *
!                                                                      *
!CYC-2024 Correction End                                               *
!
!     DVDSON is a Fortran77 program that finds a few selected
!     eigenvalues and their eigenvectors at either end of spectrum of
!     a large, symmetric (and usually sparse) matrix, denoted as A.
!     The matrix A is only referenced externally through the user
!     supplied routine which implements a block matrix-vector
!     operation(see below) and the user supplied preconditioning
!       routine (also perfomed externally). Either the range of the
!       eigenvalues wanted or an array of the indices of selected ones
!       can be specified.
!     DVDSON is a front-end routine for setting up arrays, and initial
!     guess. It also performs detailed error checking.
!     DVDRVR is the driver routine that implements a version of the
!     Davidson algorithm. The characteristics of this version are:
!      o  Use of REVERSE COMMUNICATION to perform the preconditioning
!         and the matrix vector multiply.
!      o  All arrays used by the program are stored in MEMORY.
!      o  BLOCK method (many vectors may be targeted per iteration.)
!      o  Eigenvectors are targeted in an optimum way without
!         the need to compute all unconverged residuals,
!      o  Uses two Modified Gram-Schmidt passes for orthogonality loss.
!      o  Finds HIGHEST eigenpairs by using the negative of the A.
!      o  Finds SELECTED eigenpairs specified by the user.
!      o  If not enough initial vectors are available (NIV<NUME)
!           it builds the first NUME-NIV Lanczos from those given.
!      o  The user can provide STOPPING CRITERIA for eigenvalues,
!         and residuals and coefficients and he  can CONTROL block size.
!      o  On exit INFORMATION is given about the convergence status
!         of eigenpairs.
!
!     The program consists of the following routines:
!     DVDSON, INITDVD, DVDRVR, ADDS, TSTSEL, MULTBC, OVFLOW,
!       NEWVEC, MGS_NRM
!
!     It also calls some basic BLAS routines:
!     DCOPY, DSCAL, DDOT, DAXPY, IDAMAX, DGEMV, DINIT
!
!     For solving the small eigenproblem, the routine DSPEVX from
!     LAPACK is used. DSPEVX is obtainable from NETLIB, together
!     with a series of subroutines that it calls.
!
!     All the routines have IMPLICIT REAL*8          (A-H,O-Z)
!
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!   Editted by C. F. Fischer                                       5/10/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE initdvd_I
      USE dvdrvr_I
      !USE dscal_I
      !USE dcopy_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IRC
      INTEGER, INTENT(IN)    :: N
      INTEGER, INTENT(IN)    :: LIM
      INTEGER  :: NOC
      INTEGER, INTENT(INOUT) :: ILOW
      INTEGER, INTENT(INOUT) :: IHIGH
      INTEGER  :: NIV
      INTEGER  :: MBLOCK
      INTEGER  :: MAXITER
      INTEGER, INTENT(IN) :: IWRSZ
      INTEGER, INTENT(IN) :: IIWSZ
      INTEGER  :: NLOOPS
!      INTEGER, INTENT(INOUT) :: IERR
      INTEGER  :: IERR
      REAL(DOUBLE)  :: CRITE
      REAL(DOUBLE)  :: CRITC
      REAL(DOUBLE)  :: CRITR
      LOGICAL  :: HIEND
      INTEGER  :: IREV(*)
      INTEGER  :: ISELEC(LIM)
      INTEGER  :: IWORK(IIWSZ)
      REAL(DOUBLE)  :: WORK(IWRSZ)
      real(kind(0.0d0)) :: ddot
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NEIG, I, NUME, IORTHO, IBASIS, IEIGVAL, IAB, IS, ITEMPS, ISVEC&
         , ISCRA1, IOLDVAL, ISCRA2, ISCRA3, IICV
      SAVE
!      SAVE NEIG, I, NUME, IORTHO, IBASIS, IEIGVAL, IAB, IS, ITEMPS, ISVEC, &
!         ISCRA1, IOLDVAL, ISCRA2, ISCRA3, IICV
!-----------------------------------------------
!
!       (this is for the reverse communication)

!-----------------------------------------------------------------------
!  (Important to the following is the concept of NUME, the distance of
!   the index of the eigenpair wanted which is farthest from the
!   extremes,i.e.,
!      if  lowest  eigepairs i1<i2<...<ik are wanted, NUME=ik
!      if highest eigenpairs i1<i2<...<ik are wanted, NUME=N-i1+1
!   where i1,...,ik are the indices of the wanted eigenpairs.
!   Obviously, NUME.GE.(No. of EiGenpairs wanted). )

!   on entry
!   -------
!   IRC         REVERSE COMMUNICATION integer.
!               IRC=0   --->>  Exit
!               IRC=1   --->>  Preconditioning
!               IRC=2   --->>  Matrix vector for initial setups.
!               IRC=3   --->>  Matrix vector in DVDRVR
!        NOTE Usually in reverse communication the routine called
!        ---- gives out two arrays w1 and w2. The first is used
!      as input for calculation Aw1->w2, or solve Mw2=w1.
!      Then these are passed back to the routine.
!      TO avoid extra work space, AND to be able to
!      perform block operations (i.e., where w1 is m blocks)
!      we use the same WORK array, where BASIS and AB are
!      stored in. An integer array with 6 integers is passed out,
!      except IRC: IREV.
!      IREV(1)=NB Simply NB operations must be performed
!                 from either the preconditioner or the matvec.
!          Also the number of targeted eigenpairs.
!      IREV(2)=iw1
!      IREV(3)=iw2
!              A*WORK(iw1:iw1+NB*N) --> WORK(iw2:iw1+NB*N)
!              Solve NB rhsides:
!              M*WORK(iw2:iw1+NB*N) = WORK(iw1:iw1+NB*N)
!              Obviously: iw1=ibasis+("indexB"-1)*N
!                  iw2=iab   +("indexD"-1)*N
!      for the matvec, while for the preconditioner it is
!      in reverse order.
!      IREV(4)=ieigval  (To denote where the eigenvalues are)
!      IREV(5)=iscra3   (where the eigenvalue index starts)
!              Simply the targeted eigenvalues are:
!           do i=1,NB
!          indx = IWORK(iscra3+i-1)
!          L_i  = WORK(ieigval+indx-1)
!                          solve (M-L_i) w2 = w1
!       enddo
!               The following is additional information that advanced
!      users may need.
!      IREV(6)=iscra1   (where the residuals of the corresponding
!             eigenpairs are. According to the above index)
!      IREV(7)=iscra1+LIM  this will hold the \delta eps shifts
!      from robust preconditioning
!
!   N           Order of the matrix.
!   LIM  The upper limit on the dimension of the expanding basis.
!      NUME.LT.LIM.LE.N must hold. The case LIM=NUME is allowed
!      only for LIM=NUME=N. The choice of LIM depends on the
!      available workspace (see below). If the space is
!      available it is preferable to have a large LIM, but not
!      larger than NUME$+$40.
!   NOC  Number of Orthogonalization Constraints (vectors).
!             These NOC vectors are given on input in the first
!      N*NOC elements of WORK.
!   ILOW The index of the lowest eigepair to be computed. If
!      (ILOW.LE.0).or.(ILOW.GT.N), the selected eigenpairs
!      to be computed should be contained in array ISELEC.
!      (Modified on exit).
!   IHIGH  The index of the highest eigenpair to be computed.
!      Considered ONLY when ILOW is in the range
!      (0.LT.ILOW.LE.N). (Modified on exit).
!   ISELEC Array of size LIM holding the user specified indices
!      for the eigenpairs to be computed. Considered only when
!      (ILOW.LE.0).or.(ILOW.GT.N). The indices are read from
!      the first position until a non positive integer is met.
!         Example: if N=500, ILOW=0, and ISELEC(1)=495,
!         ISELEC(2)=497, ISELEC(3)=-1, the program will find
!         2 of the highest eigenpairs, pairs 495 and 497.
!      Any order of indices is acceptable (Modified on exit).
!   NIV  Number of Initial Vector estimates provided by the user.
!      If LIM > NIV > 0 then the NIV columns of size N of WORK
!      coming after the NOC constraint vectors,  should contain
!               the estimates (see below).
!               If NIV = 0, the vector (1,1,1,...,1)T is chosen.
!   MBLOCK Number of vectors to be targeted in each iteration.
!      1.LE.MBLOCK.LE.(No. EiGenpairs wanted) should hold.
!      Large block size reduces the number of iterations
!      (matrix acceses) but increases the matrix-vector
!      multiplies. It should be used when the matrix accese
!      is expensive (disc, recomputed or distributed).
!   CRITE Convergence threshold for eigenvalues.
!      If ABS(EIGVAL-VALOLD) is less than CRITE for all wanted
!      eigenvalues, convergence is signaled.
!   CRITC       Convergence threshold for the coefficients of the last
!               added basis vector(s). If all of those corresponding
!               to unconverged eigenpairs are less than CRITC convergence
!               is signaled.
!   CRITR Convergence threshold for residual vector norms. If
!      all the residual norms ||Ax_i-l_ix_i|| of the targeted
!      x_i are less than CRITR convergence is signaled.
!      If ANY of the criteria are satisfied the algorithm stops
!   MAXITER Upper bound on the number of iterations of the
!      algorithm. When MAXITER is exceeded the algorithm stops.
!      A typical MAXITER can be MAX(200,NUME*40), but it can
!      be increased as needed.
!   WORK Real array of size IWRSZ. Used for both input and output
!      If NIV is in (1.LE.(NIV).LE.(LIM)), on input, WORK
!      must have the NIV initial estimates. These NIV N-element
!      vectors start from WORK(1) and continue one after the
!      other. They must form an orthonormal basis.
!   IWRSZ The size of the real workspace. It must be at least as
!      large as:
!
!      Not this any more:
!           (2*LIM+NOC)*N + LIM*LIM + (NUME+10)*LIM + NUME
!      New workspace:
!     --> this       (2*LIM+NOC)*N + LIM*LIM*2 +      11*LIM + NUME
!
!   IWORK  Integer work array of size IIWSZ. Used as scrath array
!      for indices and for use in the LAPACK routines.
!   IIWSZ The size of the integer workspace. It must be at least
!      as large as:
!                    6*LIM + NUME
!
!      If LIM or NUME needs to be increased, the space should
!               also be increased accordingly. For given IWRSZ and
!      IIWSZ one can calculate how big a problem one can
!      solve (LIM,NUME).
!
!   on exit
!   -------
!   WORK(NOC*N+1)
!      The first NUME*N locations contain the approximations to
!      the NUME extreme eigenvectors. If the lowest eigenpairs
!      are required, (HIEND=false), eigenvectors appear in
!      ascending order, otherwise (HIEND=false), they appear in
!      descending order. If only some are requested, the order
!      is the above one for all the NUME extreme eigenvectors,
!      but convergence has been reached only for the selected
!      ones. The rest are the current approximations to the
!      non-selected eigenvectors.
!   WORK(NOC*N+NUME*N+1)
!      The next NUME locations contain the approximations to
!      the NUME extreme eigenvalues, corresponding to the above
!      NUME eigenvectors. The same ordering and convergence
!      status applies here as well.
!   WORK(NOC*N+NUME*N+NUME+1)
!      The next NUME locations contain the corresponding values
!      of ABS(EIGVAL-VALOLD) of the NUME above eigenvalues, of
!      the last step of the algorithm.
!   WORK(NOC*N+NUME*N+NUME+NUME+1)
!      The next NUME locations contain the corresponding
!      residual norms of the NUME above eigenvectors, of the
!      last step.
!   NLOOPS      Number of iterations (accesses to matrix)
!
!   HIEND Logical. If .true. on exit the highest eigenpairs are
!      found in descending order. Otherwise, the lowest
!      eigenpairs are arranged in ascending order.
!   IERR An integer denoting the completions status:
!      IERR = 0  denotes normal completion.
!      IERR = -k  denotes error in DSPEVX (k eigenpairs
!         not converged)
!      0<IERR<=2048 denotes some inconsistency as follows:
!      If (INT( MOD(IERR,  2)/1  ) N < LIM
!      If (INT( MOD(IERR,  4)/2  ) LIM < 1
!      If (INT( MOD(IERR,  8)/4  ) ISELEC(1)<1, and no range specified
!      If (INT( MOD(IERR, 16)/8  ) IHIGH > N (in range or ISELEC)
!      If (INT( MOD(IERR, 32)/16 ) IHIGH < ILOW (Invalid range)
!      If (INT( MOD(IERR, 64)/32 ) NEIG >= LIM (Too many wanted)
!      If (INT( MOD(IERR,128)/64 ) Probable duplication in ISELEC
!      If (INT( MOD(IERR,256)/128) NUME >= LIM (max eigen very far)
!      If (INT( MOD(IERR,512)/256) MBLOCK is out of bounds
!      If (INT( MOD(IERR,1024)/512) IWRSZ or IIWSZ is not enough
!      If (INT( MOD(IERR,2048)/1024) Orthogonalization Failed
!      If (INT( MOD(IERR,4096)/2048) NLOOPS > MAXITER
!      If (INT( MOD(IERR,8192)/4096) Not proper number of NIV
!
!      The program will also print an informative message to
!      the standard output when NIV is not proper
!-----------------------------------------------------------------------
!     :    ILOW,IHIGH,ISELEC,NIV,MBLOCK,
!     :    CRITE,CRITC,CRITR,MAXITER,
!     :    WORK,IWRSZ,IWORK,IIWSZ,
!     :    HIEND,NLOOPS,IERR)
!-----------------------------------------------------------------------
!       Reverse Communication
!-----------------------------------------------------------------------
      IF (IRC==1 .OR. IRC==3) THEN
!          ..(=1)preconditioning. Go into DVDRVR
!          ..(=3)or Matrix vector in DVDRVR
         GO TO 200
      ELSE IF (IRC == 2) THEN
!          ..matrix vector for setting up.
         GO TO 100
      ENDIF
!-----------------------------------------------------------------------
!
! Checking user input errors, and setting up the problem to solve.
!
      IERR = 0
      IF (LIM > N) IERR = IERR + 1
      IF (LIM <= 0) IERR = IERR + 2

      HIEND = .FALSE.

      IF (ILOW<=0 .OR. ILOW>N) THEN
!          ..Look for user choice of eigenpairs in ISELEC
         IF (ISELEC(1) <= 0) THEN
!             ..Nothing is given in ISELEC
            IERR = IERR + 4
         ELSE
!             ..Find number of eigenpairs wanted, and their
!           ..min/max indices
            NEIG = 1
            ILOW = ISELEC(1)
            IHIGH = ISELEC(1)
            DO I = 2, LIM
               IF (ISELEC(I) <= 0) EXIT
               ILOW = MIN(ILOW,ISELEC(I))
               IHIGH = MAX(IHIGH,ISELEC(I))
               NEIG = NEIG + 1
            END DO
!           ..Check if a very large index is asked for
            IF (IHIGH > N) IERR = IERR + 8
         ENDIF
      ELSE
!          ..Look for a range between ILOW and IHIGH
!          ..Invalid range. IHIGH>N
         IF (IHIGH > N) IERR = IERR + 8
         NEIG = IHIGH - ILOW + 1
!          ..Invalid range. IHIGH<ILOW
         IF (NEIG <= 0) IERR = IERR + 16
         IF (NEIG > LIM) THEN
!           ..Not enough Basis space. Increase LIM or decrease NEIG
            IERR = IERR + 32
         ELSE
!           ..Fill in the ISELEC with the required indices
            DO I = 1, NEIG
               ISELEC(I) = ILOW + I - 1
            END DO
         ENDIF
      ENDIF

      IF (IERR /= 0) RETURN
!CYC: Restart mode
      !NUME = IHIGH
!set NUME=NIV. In step-by-step diagonalization, IHIGH is changed in
!different calling dvdson
        NUME = NIV

!       ..Identify if few of the highest eigenpairs are wanted.
      IF (ILOW + IHIGH - 1 > N) THEN
         HIEND = .TRUE.
         NUME = N - ILOW + 1
!          ..Change the problem to a minimum eipenpairs one
!          ..by picking the corresponding eigenpairs on the
!          ..opposite side of the spectrum.
         I = 1
         IF (NEIG > 0) THEN
            ISELEC(:NEIG) = N - ISELEC(:NEIG) + 1
            I = NEIG + 1
         ENDIF
      ENDIF
!       ..duplications in ISELEC
      IF (NEIG > NUME) IERR = IERR + 64
!       ..Not enough Basis space. Increase LIM or decrease NUME
      IF (NUME>LIM .OR. NUME==LIM .AND. NUME/=N) IERR = IERR + 128
!     ..Size of Block out of bounds
      IF (MBLOCK<1 .OR. MBLOCK>NEIG) IERR = IERR + 256

!       ..Check for enough workspace for Dvdson
      IF (IWRSZ<LIM*(2*N + LIM + (NUME + 10)) + NUME + NOC*N .OR. IIWSZ<6*LIM+&
         NUME) IERR = IERR + 512

      IF (NIV > LIM) THEN
!        ..Check number of initial estimates NIV is lower than LIM.
         WRITE (6, *) 'WARNING: Too many initial estimates.?'
         WRITE (6, *) 'The routine needs at most:', LIM
         IERR = IERR + 4096
      ELSE IF (NIV < NUME) THEN
!        ..check if enough initial estimates.
         WRITE (6, *) 'WARNING: Not enough initial estimates'
         WRITE (6, *) NUME - NIV, ' Lanczos vectors will be added'
      ENDIF

      IF (IERR /= 0) RETURN
!
! Assigning space for the real work arrays
!
      IORTHO  = 1
      IBASIS  = IORTHO  + N*NOC
      IEIGVAL = IBASIS  + N*LIM
      IAB     = IEIGVAL + LIM
      IS      = IAB     + N*LIM
      ITEMPS  = IS      + LIM*(LIM + 1)/2
      ISVEC   = ITEMPS  + LIM*(LIM + 1)/2
!CC   iscra1  = iSvec   + LIM*(NUME+1)
      ISCRA1  = ISVEC   + LIM*LIM
      IOLDVAL = ISCRA1  + 8*LIM
!
! Assigning space for the integer work arrays
!
      ISCRA2 = 1
      ISCRA3 = ISCRA2 + 5*LIM
      IICV = ISCRA3 + LIM
!
! Initialize tha basis, the AB, the S.
!

  100 CONTINUE
      CALL INITDVD (IRC, IREV, N, NOC, NIV, NUME + 1, LIM, HIEND, WORK(ISCRA1)&
         , WORK(IORTHO), WORK(IBASIS), WORK(IAB), WORK(IS))
!       ----------------------------------------------------------------
!       ..Reverse Communication for possible matrix vector
      IF (IRC == 2) THEN
         IREV(2) = IBASIS + (IREV(2)-1)*N
         IREV(3) = IAB + (IREV(3)-1)*N
         RETURN
      ENDIF
!       ----------------------------------------------------------------
!
! Call main driver routine.
!
      NLOOPS = 1

  200 CONTINUE
      CALL DVDRVR (IRC, IREV, N, HIEND, LIM, MBLOCK, NOC, NUME, NIV, NEIG, &
         ISELEC, CRITE, CRITC, CRITR, MAXITER, WORK(IEIGVAL), WORK(IBASIS), &
         WORK(IORTHO), WORK(IAB), WORK(IS), WORK(ITEMPS), WORK(ISVEC), WORK(&
         ISCRA1), IWORK(ISCRA2), IWORK(ISCRA3), IWORK(IICV), WORK(IOLDVAL), &
         NLOOPS, IERR)
!       ----------------------------------------------------------------
!       some Reverse Communication
      IF (IRC == 1) THEN
!        ..Preconditioning
         IREV(2) = IAB + (IREV(2)-1)*N
         IREV(3) = IBASIS + (IREV(3)-1)*N
         IREV(4) = IEIGVAL
         IREV(5) = ISCRA3
         IREV(6) = ISCRA1
         IREV(7) = ISCRA1 + LIM
         RETURN
      ELSE IF (IRC == 3) THEN
!        ..Matrix-vector
         IREV(2) = IBASIS + (IREV(2)-1)*N
         IREV(3) = IAB + (IREV(3)-1)*N
         RETURN
      ENDIF
!       ----------------------------------------------------------------

      IF (HIEND) CALL DSCAL (NUME, -1.D0, WORK(IEIGVAL), 1)
!
! -Copy the eigenvalues after the eigenvectors
! -Next, copy the difference of eigenvalues between the last two steps
! -Next, copy the residuals for the first NUME estimates
!
      CALL DCOPY (NUME, WORK(IEIGVAL), 1, WORK(IBASIS+N*NUME), 1)
      CALL DCOPY (NUME, WORK(IOLDVAL), 1, WORK(IBASIS+(N+1)*NUME), 1)
      CALL DCOPY (NUME, WORK(ISCRA1), 1, WORK(IBASIS+(N+2)*NUME), 1)
!
! Set IRC=0 for normal exit with no reverse communication
!
      IRC = 0
      RETURN
      END SUBROUTINE DVDSON


!=======================================================================
      SUBROUTINE ADDS(N, LIM, KPASS, NNCV, BASIS, AB, S)
!=======================================================================
!     Called by: DVDSON
!
!     Calculates the new column in the S matrix. S has a
!     new row and column, but being symmetric only the new column is
!     stored. S(i,kpass+1)=B(i)^T D(kpass+1) for all i.
!
!     subroutines called:
!     DDOT, DSCAL
!
! Modified for MPI version, Chongyang Chen, 2020/01/01
!
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE ddot_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)    :: N
      INTEGER , INTENT(IN) :: LIM
      INTEGER , INTENT(IN) :: KPASS
      INTEGER , INTENT(IN) :: NNCV
      REAL(DOUBLE)  :: BASIS(N*LIM)
      REAL(DOUBLE)  :: AB(N*LIM)
      REAL(DOUBLE) , INTENT(OUT) :: S(LIM*(LIM + 1)/2)
      real(kind(0.0d0)) :: ddot
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IDSTART, ISSTART, IV, IBSTART, IBV
      REAL(DOUBLE) :: SS
!-----------------------------------------------
!-----------------------------------------------------------------------
!   on entry
!   -------
!   N  The order of the matrix A
!   kpass The current dimension of the expanding sub-basis
!   NNCV Number of new basis vectors.
!   Basis the basis vectors, including the new NNCV ones.
!   AB  The matrix D=AB. (with new NNCV columns)
!   on exit
!   -------
!   S           The small matrix with NNCV new columns at the last part
!-----------------------------------------------------------------------
!
! The new S is calculated by adding the new last columns
! S(new)=B^T D(new).
!
      IDSTART = KPASS*N + 1
      ISSTART = KPASS*(KPASS + 1)/2
      DO IV = 1, NNCV
         IBSTART = 1
         DO IBV = 1, KPASS + IV
            !SS = DDOT(N,BASIS(IBSTART),1,AB(IDSTART),1)
            call DDOTMPI(N,BASIS(IBSTART),AB(IDSTART),SS)
            S(ISSTART+IBV) = SS
            IBSTART = IBSTART + N
         END DO
         ISSTART = ISSTART + KPASS + IV
         IDSTART = IDSTART + N
      END DO

      RETURN
      END SUBROUTINE ADDS


!=======================================================================
      SUBROUTINE DVDRVR(IRC, IREV, N, HIEND, LIM, MBLOCK, NOC, NUME, NIV, NEIG&
         , ISELEC, CRITE, CRITC, CRITR, MAXITER, EIGVAL, BASIS, ORTHOBASIS, AB&
         , S, TEMPS, SVEC, SCRA1, ISCRA2, INCV, ICV, OLDVAL, NLOOPS, IERR)
!=======================================================================
!     called by DVDSON
!
!     Driver routine implementing Davidson's main loop. On entry it
!     is given the Basis, the work matrix D=AB and the small symmetric
!     matrix to be solved, S=B^TAB (as found by ADDS). In each step
!     the small problem is solved by calling DSPEVX.
!     TSTSEL tests for eigenvalue convergence and selects the next
!     pairs to be considered for targeting (as a block).
!       NEWVEC computes the new vectors (block) to be added in the
!     expanding basis, and tests for residual convergence.
!     After Preconditioning ((M-l_i)d=res), the basis is orthogonalized.
!       Then the matrix-vector multiplication finds the new vectors of D
!       (Dnew=ABnew), and the new small problem S, is calculated.
!       The algorithm is repeated.
!     In case of a large expanding basis (KPASS=LIM) the Basis, AB,
!     SVEC and S are collapsed.
!     At the end the current eigenvector estimates are computed as
!     well as the residuals and eigenvalue differences.
!
!     Subroutines called:
!     DSPEVX, MULTBC, TSTSEL, OVFLOW, NEWVEC, ADDS,
!     DCOPY, DDOT, DAXPY
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
!     USE MPI_C,           ONLY: MYID
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE tstsel_I
      USE newvec_I
      USE mgs_nrm_I
      USE adds_I
      USE multbc_I

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(INOUT) :: IRC
      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(IN)  :: LIM
      INTEGER  :: MBLOCK
      INTEGER, INTENT(IN) :: NOC
      INTEGER, INTENT(IN) :: NUME
      INTEGER, INTENT(IN) :: NIV
      INTEGER, INTENT(IN) :: NEIG
      INTEGER, INTENT(IN) :: MAXITER
      INTEGER, INTENT(INOUT) :: NLOOPS
      INTEGER, INTENT(INOUT) :: IERR
      REAL(DOUBLE)  :: CRITE
      REAL(DOUBLE)  :: CRITC
      REAL(DOUBLE)  :: CRITR
      LOGICAL, INTENT(IN) :: HIEND
      INTEGER, INTENT(OUT) :: IREV(*)
      INTEGER  :: ISELEC(NEIG)
      INTEGER  :: ISCRA2(5*LIM)
      INTEGER, INTENT(INOUT)  :: INCV(LIM)
      INTEGER, INTENT(OUT)  :: ICV(NUME + 1)
      REAL(DOUBLE), INTENT(INOUT)  :: EIGVAL(LIM)
      REAL(DOUBLE)  :: BASIS(N*LIM)
      REAL(DOUBLE)  :: ORTHOBASIS(N*LIM + NOC*N)
      REAL(DOUBLE)  :: AB(N*LIM)
      REAL(DOUBLE)  :: S(LIM*(LIM + 1)/2)
      REAL(DOUBLE)  :: TEMPS(LIM*(LIM + 1)/2)
      REAL(DOUBLE)  :: SVEC(LIM*LIM)
      REAL(DOUBLE), INTENT(INOUT)  :: SCRA1(8*LIM)
      REAL(DOUBLE), INTENT(INOUT)  :: OLDVAL(NUME + 1)
      real(kind(0.0d0)) :: ddot
      real::slamch
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: REST, I, KPASS, NNCV, IFIND, NFOUND, INFO, NEWSTART
      REAL(DOUBLE) :: TOL, DSUM
      LOGICAL :: FIRST, DONE

      SAVE FIRST, DONE, REST, I, KPASS, NNCV, IFIND, TOL, NFOUND, INFO, &
         NEWSTART
!-----------------------------------------------
!CC   DIMENSION SVEC(LIM*(NUME+1)),EIGVAL(LIM)
!
!       !include 'mpif.h'
!       (this is for the reverse communication)

!-----------------------------------------------------------------------
!
!   on entry
!   -------
!
!   N           The order of the matrix A
!   HIEND Logical. True only if the highest eigenpairs are needed.
!   LIM  The limit on the size of the expanding Basis
!   MBLOCK Number of vectors to be targeted in each iteration.
!   NOC  Number of orthogonalization constraints (vectors)
!   NUME        The largest index of the eigenvalues wanted.
!   NIV  Starting dimension of expanding basis.
!   NEIG  Number of eigenvalues wanted.
!   ISELEC Array containg the indices of those NEIG eigenpairs.
!   CRITE       Convergence thresholds for eigenvalues, coefficients
!   CRITC,CRITR and residuals.
!   BASIS  Array with the basis vectors.
!   AB        Array with the vectors D=AB
!   S  Array keeping the symmetric matrix of the small problem.
!   TEMPS scratch array
!   SVEC Array for holding the eigenvectors of S
!   SCRA1  Srcatch array used by DSPEVX. It also holds residuals
!      and eigenvalue shifts for preconditioning.
!   ISCRA2 Integer Srcatch array used by DSPEVX.
!   INCV Srcatch array used in DSPEVX. Also used in TSTSEL and
!      NEWVEC where it holds the Indices of uNConVerged pairs
!   ICV  It contains "1" to the locations of ConVerged eigenpairs
!   OLDVAL Array keeping the previous' step eigenvalue estimates.
!
!   on exit
!   -------
!
!   EIGVAL Array containing the NUME lowest eigenvalues of the
!      the matrix A (or -A if the highest are sought).
!   Basis  On exit Basis stores the NUME corresponding eigenvectors
!   OLDVAL  On exit it stores the final differences of eigenvalues.
!   SCRA1 On exit it stores the NUME corresponding residuals.
!   NLOOPS Number of loops taken by the algorithm
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Reverse Communication
      IF (IRC == 1) THEN
!          ..Came from preconditioning
         GO TO 100
      ELSE IF (IRC == 3) THEN
!          ..came from matrix vector multiply
         GO TO 200
      ENDIF
!-----------------------------------------------------------------------

      I = 1
      IF (NUME > 0) THEN
         EIGVAL(:NUME) = 1.D30
         ICV(:NUME) = 0
         I = NUME + 1
      ENDIF
      FIRST = .TRUE.
      KPASS = NIV
      NNCV = KPASS
!
! Decide HERE how many to restart with, not more than LIM-3.
!
      REST = MAX(LIM/2,2*NUME)
      REST = MIN(REST,LIM - 3)

   10 CONTINUE
      IFIND = MIN(KPASS,REST)
      TOL = SLAMCH('S')
      CALL DCOPY (IFIND, EIGVAL, 1, OLDVAL, 1)
      CALL DCOPY ((KPASS*(KPASS + 1))/2, S, 1, TEMPS, 1)
      CALL DSPEVX ('Vectors also', 'In a range', 'Upper triangular', KPASS, &
         TEMPS, -1., -1., 1, IFIND, TOL, NFOUND, EIGVAL, SVEC, KPASS, SCRA1, &
         ISCRA2, INCV, INFO)
      IERR = -ABS(INFO)
      IF (IERR /= 0) GO TO 60
!
! TeST for convergence on the absolute difference of eigenvalues between
! successive steps. Also SELect the unconverged eigenpairs and sort them
! by the largest magnitude in the last added NNCV rows of Svec.
!
      DONE = TSTSEL(NLOOPS,KPASS,NUME,NEIG,ISELEC,SVEC,EIGVAL,ICV,CRITE,CRITC,SCRA1,&
         ISCRA2,OLDVAL,NNCV,INCV)
      IF (DONE .OR. KPASS>=N) GO TO 30
!
! Maximum size for expanding basis. Truncate basis, D, and S, Svec
! Consider the basis vectors found in TSTSEL for the newvec. KPASS=NUME
!

!          Change suggested by Anreas, March 18, 1996
      IF (KPASS == LIM) THEN
!      PRINT*,'collapsing the basis: lim,nume,rest',lim,nume,rest
!        PRINT*,'myid = ', myid, '  nprocs = ', nprocs
!        PRINT*,'collapsing the basis: lim,nume,rest',lim,nume,rest,myid
!23456789012345678901234567890123456789012345678901234567890123456789012

         CALL OVFLOW (N, REST, KPASS, SCRA1, BASIS, AB, S, SVEC, EIGVAL)
!             CALL OVFLOW(N,NUME,KPASS,SCRA1,BASIS,AB,S,SVEC,EIGVAL)
      ENDIF
!
! Compute and add the new residuals. NNCV is set to the number of new
! vectors that have not converged. If none, DONE=true, exit.
!
      CALL NEWVEC (N, NLOOPS, NUME, LIM, MBLOCK, KPASS, CRITR, NNCV, INCV, SVEC, EIGVAL&
         , OLDVAL, AB, BASIS, ICV, SCRA1, SCRA1(LIM+1), DONE)

      IF (DONE) GO TO 30
!-----------------------------------------------------------------------
! Preconditioning the NNCV (residuals-deps x) stored in AB(kpass+1).
!          ..Robust Preconditioning (Eigenvalue shift).
!          ..Look for lowest, so Li-|eps_i|
      DO I = 1, NNCV
         EIGVAL(INCV(I)) = EIGVAL(INCV(I)) - SCRA1(LIM+INCV(I))
      END DO
!
!          ..Change sign of eigenvalues if HIEND.
      IF (HIEND) CALL DSCAL (NUME, -1.D0, EIGVAL, 1)

! Use of Reverse Communication.
!
      IREV(1) = NNCV
      IREV(2) = KPASS + 1
      IREV(3) = KPASS + 1
      IRC = 1
      IF (IRC == 1) RETURN
!          ..Continue from preconditioning
  100 CONTINUE
      IF (HIEND) CALL DSCAL (NUME, -1.D0, EIGVAL, 1)
!
!          ..Shift the eigenvalues back to what they were.
      DO I = 1, NNCV
         EIGVAL(INCV(I)) = EIGVAL(INCV(I)) + SCRA1(LIM+INCV(I))
      END DO
!-----------------------------------------------------------------------
!
! Orthonormalization of the previous vectors to the Basis and to any
! orthogonalization constraints. The not-yet-filled
! spaces of AB (from NEWSATRT onwards) are used for scratch.
!
      NEWSTART = KPASS*N + 1

      CALL MGS_NRM (N, NOC + KPASS, NNCV, SCRA1(LIM+1), ORTHOBASIS)
!-----------------------------------------------------------------------
! Use of Reverse Communication. Add new columns in D through matrix vector
! multiplication CALL OP(N,NNCV,BASIS(NEWSTART),AB(NEWSTART))
!
      IREV(1) = NNCV
      IREV(2) = KPASS + 1
      IREV(3) = KPASS + 1
      IRC = 3
      RETURN

!          ..Continue from matrix-vector multiply
  200 CONTINUE
      IF (HIEND) CALL DSCAL (N*NNCV, -1.D0, AB(NEWSTART), 1)
!
! Add new column in S, from the NNCV new vectors.
!
      CALL ADDS (N, LIM, KPASS, NNCV, BASIS, AB, S)

      KPASS = KPASS + NNCV
      NLOOPS = NLOOPS + 1

      IF (NLOOPS <= MAXITER) GO TO 10
      IERR = IERR + 2048
      NLOOPS = NLOOPS - 1
      KPASS = KPASS - NNCV
   30 CONTINUE
      DO I = 1, NUME
         OLDVAL(I) = ABS(OLDVAL(I)-EIGVAL(I))
      END DO

      CALL MULTBC (N, KPASS, NUME, SVEC, SCRA1, BASIS)
      CALL MULTBC (N, KPASS, NUME, SVEC, SCRA1, AB)
!
! i=1,NUME residual(i)= DCi-liBCi= newDi-linewBi
! temporarily stored in AB(NUME*N+1)
!
      DO I = 1, NUME
         CALL DCOPY (N, AB((I-1)*N+1), 1, AB(NUME*N+1), 1)
         CALL DAXPY (N, (-EIGVAL(I)),BASIS((I-1)*N+1), 1, AB(NUME*N+1), 1)
         !SCRA1(I) = DDOT(N,AB(NUME*N+1),1,AB(NUME*N+1),1)
         !SCRA1(I) = SQRT(SCRA1(I))
         call DDOTMPI(N,AB(NUME*N+1),AB(NUME*N+1),DSUM) 
         SCRA1(I) = SQRT(DSUM)
      END DO
!
! Set IRC=0 for normal exit with no reverse communication
!
!     IF (MYID == 0) WRITE (6, *) 'DVDSON::NLOOPS =', NLOOPS
   60 CONTINUE
      IRC = 0
      RETURN
      END SUBROUTINE DVDRVR


!=======================================================================
      SUBROUTINE INITDVD(IRC, IREV, N, NOC, NIV, NUME, LIM, HIEND, SCRA1, &
         ORTHOBASIS, BASIS, AB, S)
!=======================================================================
!     Initializes the basis and the auxiliary arrays AB and S.
!     If not enough initial estimates exist the basis will be
!     supplemented with Lanczos vectors of the current NIV vectors.
!     eg. (b1,b2) --> (b1,Ab1,AAb1, b2,Ab2,AAb2) for a NUME of 6.
!
!     OrthoBasis is the Basis with the NOC orthogonalization
!     constraint vectors in the begining. This equialence holds:
!           equivalence(Basis(1),OrthoBasis(NOC*N+1))
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dinit_I
      !USE dcopy_I
      USE mgs_nrm_I
      !USE dscal_I
      USE adds_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(INOUT) :: IRC
      INTEGER, INTENT(IN)    :: N
      INTEGER, INTENT(IN) :: NOC
      INTEGER, INTENT(INOUT)  :: NIV
      INTEGER, INTENT(IN) :: NUME
      INTEGER, INTENT(IN) :: LIM
      LOGICAL, INTENT(IN) :: HIEND
      INTEGER, INTENT(OUT) :: IREV(*)
      REAL(DOUBLE)  :: SCRA1(*)
      REAL(DOUBLE)  :: ORTHOBASIS(N*(NOC + LIM))
      REAL(DOUBLE)  :: BASIS(N*LIM)
      REAL(DOUBLE)  :: AB(N*LIM)
      REAL(DOUBLE)  :: S(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IST, IEND, NEW, KPASS

      SAVE IST, IEND, NEW, KPASS
!-----------------------------------------------
!
!-----------------------------------------------------------------------
!       ..Reverse Communication
      IF (IRC == 2) GO TO 100
!-----------------------------------------------------------------------
!
! If no initial estimates pick one
!
      IF (NIV == 0) THEN
         CALL DINIT (N, 1.D0/SQRT(DBLE(N)), BASIS, 1)
         NIV = 1
      ENDIF
!
! Compute AB. Also fill basis with orthonormalized ABs until enough NIVs.
!
      IST = 1
      IEND = NIV
      NEW = NIV
   10 CONTINUE
      IRC = 2
      IREV(1) = NEW
      IREV(2) = IST
      IREV(3) = IST
      RETURN
  100 CONTINUE
      IF (IEND < NUME) THEN
         NEW = MIN(NUME - IEND,IEND - IST + 1)
         CALL DCOPY (N*NEW, AB((IST-1)*N+1), 1, BASIS(1+IEND*N), 1)
!            ..orthonormalize OrthoBasis i.e.,
!            B-1..B-noc,B1,...,Biend,( Biend+1,...Biend+new )
         CALL MGS_NRM (N, NOC + IEND, NEW, SCRA1, ORTHOBASIS)
         IST = IEND + 1
         IEND = IEND + NEW
         GO TO 10
      ENDIF
      NIV = IEND
      !xhh print*, 'niv=',niv, 'nume=', nume
!
! Scale if HIEND for highest eigepairs
!
      IF (HIEND) CALL DSCAL (N*NIV, -1.D0, AB, 1)
!
! Also find the small matrix S = B^TAB.
!
      KPASS = 0
      CALL ADDS (N, LIM, KPASS, NIV, BASIS, AB, S)

      IRC = 0
      RETURN
      END SUBROUTINE INITDVD



!=======================================================================
      SUBROUTINE MULTBC(N, K, M, C, TEMP, B)
!=======================================================================
!     called by: DVDRVR
!
!     Multiplies B(N,K)*C(K,M) and stores it in B(N,M)
!     Used for collapsing the expanding basis to current estimates,
!     when basis becomes too large, or for returning the results back
!
!       Subroutines called
!       DINIT, DGEMV, DCOPY
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dgemv_I
      !USE dcopy_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(IN)  :: K
      INTEGER, INTENT(IN)  :: M
      REAL(DOUBLE)  :: C(K*M)
      REAL(DOUBLE)  :: TEMP(M)
      REAL(DOUBLE)  :: B(N*K)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IROW
!-----------------------------------------------
!-----------------------------------------------------------------------
      DO IROW = 1, N
         CALL DGEMV ('Transp', K, M, 1.D0, C, K, B(IROW), N, 0.D0, TEMP, 1)
         CALL DCOPY (M, TEMP, 1, B(IROW), N)
      END DO

      RETURN
      END SUBROUTINE MULTBC

!=======================================================================
      SUBROUTINE MULTBC_COL(N, K, M, C, TEMP, B)
!=======================================================================
! CYC: obtain the column of B(N,M), by using DGEMVMPI
!
!     called by: DVDRVR
!
!     Multiplies B(N,K)*C(K,M) and stores it in B(N,M)
!     Used for collapsing the expanding basis to current estimates,
!     when basis becomes too large, or for returning the results back
!
!       Subroutines called
!       DINIT, DGEMV, DCOPY
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dgemv_I
      !USE dcopy_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)    :: N
      INTEGER  :: K
      INTEGER  :: M
      REAL(DOUBLE)  :: C(K*M)
      REAL(DOUBLE)  :: TEMP(M)
      REAL(DOUBLE)  :: B(N*K)
!
! CYC: Temporary array
      REAL(DOUBLE)  :: TMPWORK(N*M)
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: JCOL
!-----------------------------------------------
!-----------------------------------------------------------------------
!      DO IROW = 1, N
!         CALL DGEMV ('Transp', K, M, 1.D0, C, K, B(IROW), N, 0.D0, TEMP, 1)
!         CALL DCOPY (M, TEMP, 1, B(IROW), N)
!      END DO

      DO JCOL = 1, M
         CALL DGEMVMPI ('N', N, K, 1.D0, B, N, C((JCOL-1)*K+1),1,       &
                           0.d0, TMPWORK((JCOL-1)*N+1), 1)
      END DO
      CALL DCOPY (N*M, TMPWORK, 1, B, 1)
      CALL DCOPY (M, B(N), N, TEMP, 1)
      RETURN
      END SUBROUTINE MULTBC_COL


!=======================================================================
      SUBROUTINE NEWVEC(N, NLOOPS, NUME, LIM, MBLOCK, KPASS, CRITR, NNCV, INCV, SVEC, &
         EIGVAL, OLDVAL, AB, BASIS, ICV, SCRA1, EPSIL, DONE)
!=======================================================================
!
!     Called by: DVDRVR
!
!     It calculates the new expansion vectors of the basis.
!     For each one of the vectors in INCV starting with the largest
!     megnitude one, calculate its residual Ri= DCi-liBCi and check
!     the ||Ri|| for convergence. If it is converged do not add it
!     but look for the immediate largest coefficient and its vector.
!     The above procedure continues until MBLOCK vectors have been
!     added to the basis, or the upper limit has been encountered.
!     Thus only  the required MBLOCK residuals are computed.
!     For robust Preconditioning, each residual is then modified as
!               \delta_eps x - res
!     where \delta_eps is an estimate of the eigenvalue error.
!     A different estimate is used for the eigenvalue shift.
!
!     Subroutines called:
!     DNRM2, DGEMV
!
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE MPI_C,           ONLY: MYID
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)    :: N
      INTEGER , INTENT(IN) :: NLOOPS
      INTEGER , INTENT(IN) :: NUME
      INTEGER , INTENT(IN) :: LIM
      INTEGER , INTENT(IN) :: MBLOCK
      INTEGER, INTENT(IN)  :: KPASS
      INTEGER , INTENT(INOUT) :: NNCV
      REAL(DOUBLE) , INTENT(IN) :: CRITR
      LOGICAL , INTENT(OUT) :: DONE
      INTEGER , INTENT(INOUT) :: INCV(NUME)
      INTEGER , INTENT(OUT) :: ICV(NUME)
      REAL(DOUBLE)  :: SVEC(LIM*NUME)
      REAL(DOUBLE) , INTENT(IN) :: EIGVAL(LIM)
      REAL(DOUBLE) , INTENT(IN) :: OLDVAL(NUME)
      REAL(DOUBLE)  :: AB(N*LIM)
      REAL(DOUBLE)  :: BASIS(N*LIM)
      REAL(DOUBLE) , INTENT(INOUT) :: SCRA1(LIM)
      REAL(DOUBLE) , INTENT(OUT) :: EPSIL(LIM)
       real(kind(0.0d0)) :: ddot
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NEWSTART, NADDED, ICVC, LIMADD, ICUR, I, INDX
      REAL(DOUBLE) :: SQRES, GAPLOW, GAP, GAPUP, DL, RES, RHSEPS
!-----------------------------------------------
!-----------------------------------------------------------------------
!   on entry
!   --------
!   N           The order of the matrix A
!   NUME        The largest index of the eigenvalues wanted.
!   LIM         The limit on the size of the expanding Basis
!   MBLOCK Maximum number of vectora to enter the basis
!   KPASS the current dimension of the expanding basis
!   CRITR Convergence threshold for residuals
!   NNCV Number of Non ConVerged pairs (MBLOCK will be targeted)
!   INCV        Index to the corresponding SVEC columns of these pairs.
!   SVEC,EIGVAL Arrays holding the eigenvectors and eigenvalues of S
!   AB  Array with the vectors D=AB
!   BASIS the expanding basis having kpass vectors
!   ICV  Index of converged eigenpairs (ICV(i)=1 <=>i converged)
!
!   on exit
!   -------
!   NNCV The number of vectors finally added to the basis.
!   BASIS  The eigenector estimated reside at the end.
!   AB          The right hand sides for preconditioning.
!   EPSIL       The eigenvalue shifts for preconditioning,
!   ICV  Index of converged eigenpairs (updated)
!   DONE  logical, if covergance has been reached.
!      the Basis should be collapsed to current approximations.
!-----------------------------------------------------------------------
!       !include 'mpif.h'
      DONE = .FALSE.
      NEWSTART = KPASS*N + 1
      NADDED = 0
      ICVC = 0
      LIMADD = MIN(LIM,MBLOCK + KPASS)
      ICUR = NEWSTART
!
! Compute RESIDUALS for the MBLOCK of the NNCV not converged vectors.
!
      DO I = 1, NNCV
         INDX = INCV(I)
!
!          ..Compute b = BASIS*Svec_indx
!        ..Compute d = AB*Svec_indx
!        ..Daxpy d'= d - eigval b  gives the residual

!         CALL DGEMV ('N', N, KPASS, 1.D0, BASIS, N, SVEC((INDX-1)*KPASS+1) &
!            , 1, 0.D0, BASIS(ICUR), 1)
!         CALL DGEMV ('N', N, KPASS, 1.D0, AB, N, SVEC((INDX-1)*KPASS+1), 1, &
!            0.D0, AB(ICUR), 1)
      CALL DGEMVMPI('N',N,KPASS,1.D0,BASIS,N,SVEC((INDX-1)*KPASS+1),1, &
                       0.d0,BASIS(ICUR),1)
      CALL DGEMVMPI('N',N,KPASS,1.D0,AB,N,SVEC((INDX-1)*KPASS+1),1,    &
                       0.d0,AB(ICUR),1)

         CALL DAXPY (N, (-EIGVAL(INDX)),BASIS(ICUR), 1, AB(ICUR), 1)
!
!        ..Compute the norm of the residual
!        ..and check for convergence
!
         !SQRES = DDOT(N,AB(ICUR),1,AB(ICUR),1)
         CALL DDOTMPI(N,AB(ICUR),AB(ICUR),SQRES)
         SCRA1(INDX) = SQRT(SQRES)

!        IF (MYID == 0) WRITE (6, '(A11,F22.16,I3,A10,F19.16)') &
!           ' EIGVAL(i) ', EIGVAL(INDX), INDX, ' Res.Norm ', SCRA1(INDX)

         IF (SCRA1(INDX) < CRITR .AND. NLOOPS > 2) THEN
!           ..Converged,do not add. Go for next non converged one
            ICVC = ICVC + 1
            ICV(INDX) = 1
            IF (ICVC < NNCV) CYCLE
!           ..All have converged.
            !IF (MYID == 0) print*, 'converged by critr',myid
            DONE = .TRUE.
            RETURN
         ELSE
!           ..Not converged. Consider it for preconditioning
!          ---------------  ROBUST MODIFICATION ---------------------
!           ..Daxpy d'= d'- \delta_eps b  gives the rhs for precond.
!           ..It is stored in AB.
!             ..Deps are also stored in EPSIL
!
            IF (INDX == 1) THEN
               GAPLOW = 1.0D+99
               GAP = ABS(EIGVAL(2)-EIGVAL(1))
!                 print*, 'h1,h2:',eigval(1),eigval(2)
            ELSE
               GAPLOW = ABS(EIGVAL(INDX)-EIGVAL(INDX-1))
               GAPUP = ABS(EIGVAL(INDX+1)-EIGVAL(INDX))
               GAP = MIN(GAPLOW,GAPUP)
            ENDIF
            DL = ABS(OLDVAL(INDX)-EIGVAL(INDX))
            RES = SCRA1(INDX)

            IF (GAP > RES) THEN
!                 EPSIL(indx) = min(dl,res,gaplow)
               RHSEPS = SQRT(DL*RES)
            ELSE
!                EPSIL(indx) = min( res, gaplow)
               RHSEPS = MIN(DL,SQRT(DL))
            ENDIF

            EPSIL(INDX) = 0.D0
            CALL DAXPY (N, (-RHSEPS), BASIS(ICUR), 1, AB(ICUR), 1)

!          -------------- END OF ROBUST MODIFICATIONS -----------------
!
            NADDED = NADDED + 1
            INCV(NADDED) = INDX
            IF (NADDED + KPASS == LIMADD) EXIT
!           ..More to be added in the block
            ICUR = ICUR + N
         ENDIF
      END DO

      NNCV = NADDED

      RETURN
      END SUBROUTINE NEWVEC


!=======================================================================
      SUBROUTINE OVFLOW(N, NUME, KPASS, SCRA1, BASIS, AB, S, SVEC, EIGVAL)
!=======================================================================
!     Called by: DVDRVR
!     Called when the upper limit (LIM) has been reached for the basis
!     expansion (by KPASS). The basis is truncated B'=BC and
!     similarly the array D'=DC. The new S is computed as
!     S'(i,j)=l(i)delta(i,j) where l(i) eigenvalues, and delta of
!     Kronecker, i,j=1,NUME. The new eigenvectors of the small matrix
!     are the unit vectors.
!
!     Subroutines called:
!     DCOPY, DINIT, MULTBC
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE multbc_I
      !USE dinit_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)    :: N
      INTEGER  :: NUME
      INTEGER  :: KPASS
!CFF   ... add dimensionto SCRA1
      REAL(DOUBLE)  :: SCRA1(NUME)
      REAL(DOUBLE)  :: BASIS(N*KPASS)
      REAL(DOUBLE)  :: AB(N*KPASS)
      REAL(DOUBLE)  :: S((KPASS*(KPASS + 1))/2)
      REAL(DOUBLE)  :: SVEC(KPASS*NUME)
      REAL(DOUBLE) , INTENT(IN) :: EIGVAL(KPASS)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IND, ICUR, I
!-----------------------------------------------
!   on entry
!   -------
!   NUME The largest index of eigenvalues wanted.
!   SVEC  the kpass eigenvectors of the smaller system solved
!   EIGVAL  the eigenvalues of this small system
!   on exit
!   -------
!   Basis
!   AB
!   S   The new small matrix to be solved.
!-----------------------------------------------------------------------
! Truncate  the basis and the AB array.
!
      !CALL MULTBC (N, KPASS, NUME, SVEC, SCRA1, BASIS)
      !CALL MULTBC (N, KPASS, NUME, SVEC, SCRA1, AB)
      CALL MULTBC_COL (N, KPASS, NUME, SVEC, SCRA1, BASIS)
      CALL MULTBC_COL (N, KPASS, NUME, SVEC, SCRA1, AB)
!
! calculation of the new upper S=diag(l1,...,l_NUME) and
! its matrix Svec of eigenvectors (e1,...,e_NUME)
!
      CALL DINIT ((NUME*(NUME + 1))/2, 0.D0, S, 1)
      CALL DINIT (NUME*NUME, 0.D0, SVEC, 1)
      IND = 0
      ICUR = 0
      DO I = 1, NUME
         S(IND+I) = EIGVAL(I)
         SVEC(ICUR+I) = 1
         ICUR = ICUR + NUME
         IND = IND + I
      END DO

      KPASS = NUME

      RETURN
      END SUBROUTINE OVFLOW


!=======================================================================
      LOGICAL FUNCTION TSTSEL (NLOOPS, KPASS, NUME, NEIG, ISELEC, SVEC, EIGVAL, ICV, &
         CRITE, CRITC, ROWLAST, IND, OLDVAL, NNCV, INCV)
!=======================================================================
!
!     Called by: DVDRVR
!
!     It first checks if the wanted eigenvalues have reached
!     convergence and updates OLDVAL. Second, for each wanted and non
!     converged eigenvector, it finds the largest absolute coefficient
!       of the NNCV last added vectors (from SVEC) and if not coverged,
!     places it in ROWLAST. IND has the corresponding indices.
!     Third, it sorts ROWLAST in decreasing order and places the
!     corresponding indices in the array INCV. The index array INCV
!       and the number of unconverged pairs NNCV, are passed to DVDRVR.
!      Later in NEWVEC only the first MBLOCK of NNCV pairs
!      will be targeted, since if ROWLAST(i)>ROWLAST(j)
!      then approximately RESIDUAL(i)>RESIDUAL(j)
!
!     Subroutines called
!     IDAMAX
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
     USE MPI_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE idamax_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NLOOPS
      INTEGER , INTENT(IN) :: KPASS
      INTEGER , INTENT(IN) :: NUME
      INTEGER , INTENT(IN) :: NEIG
      INTEGER , INTENT(INOUT) :: NNCV
      REAL(DOUBLE) , INTENT(IN) :: CRITE
      REAL(DOUBLE) , INTENT(IN) :: CRITC
      INTEGER , INTENT(IN) :: ISELEC(NEIG)
      INTEGER , INTENT(OUT) :: ICV(NUME)
      INTEGER , INTENT(INOUT) :: IND(NEIG)
      INTEGER , INTENT(OUT) :: INCV(NEIG)
      REAL(DOUBLE) , INTENT(IN) :: SVEC(KPASS*NUME)
      REAL(DOUBLE) , INTENT(IN) :: EIGVAL(NUME)
      REAL(DOUBLE), INTENT(INOUT)  :: ROWLAST(NEIG)
      REAL(DOUBLE) , INTENT(IN) :: OLDVAL(NUME)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NNCE, I, IVAL, ICNT, ICUR, L, INDX, ITEMP
      REAL(DOUBLE) :: TMAX, TEMP
      LOGICAL :: DONE, IFLAG_CONVERGE
      integer :: idamax
!-----------------------------------------------
!-----------------------------------------------------------------------
!
!   on entry
!   -------
!   KPASS       current dimension of the expanding Basis
!   NUME        Largest index of the wanted eigenvalues.
!   NEIG        number of wanted eigenvalues of original matrix
!   ISELEC index array of the wanted eigenvalues.
!   SVEC        the eigenvectors of the small system
!   EIGVAL The NUME lowest eigenvalues of the small problem
!   ICV  Index of converged eigenpairs.ICV(i)=1 iff eigenpair i
!      has converged, and ICV(i)=0 if eigenpair i has not.
!   CRITE,CRITC Convergence thresholds for eigenvalues and coefficients
!   ROWLAST    scratch array, keeping the largest absolute coefficient
!      of the NNCV last rows of Svec.
!   IND  scratch array, temporary keeping the indices of Rowlast
!   OLDVAL The previous iteration's eigenvalues.
!
!   on exit
!   -------
!   NNCV   Number of non converged eigenvectors (to be targeted)
!   INCV  Index to these columns in decreasing order of magnitude
!   TSTSEL     true if convergence has been reached
!
!-----------------------------------------------------------------------
!       !include 'mpif.h'

      DONE = .FALSE.
!
!CYC-2024 
! In some cases, energies converge coincidentally by few iterations,
! Remvoe convergence test by CRITE?
!
! Test all wanted eigenvalues for convergence under CRITE
!
!CYC-
!      IF (NLOOPS == 2 .AND. MYID == 0) THEN
!        WRITE(6, *)"NLOOPS = 2, the old and new eigvals: "
!        WRITE(6, *)'OLDVAL:', OLDVAL(1:NEIG)
!        WRITE(6, *)'EIGVAL:', EIGVAL(1:NEIG)
!      ENDIF
!CYC-End
      NNCE = 0
      DO I = 1, NEIG
         IVAL = ISELEC(I)
         IF (ABS(OLDVAL(IVAL)-EIGVAL(IVAL)) < CRITE .AND. &
             NLOOPS > 2) CYCLE
         NNCE = NNCE + 1
      END DO
      IF (NNCE == 0) THEN
         IF (NLOOPS > 2) TSTSEL = .TRUE.
         !IF (MYID == 0) WRITE (6, *) 'converged by crite'
!CYC-
         IF (NLOOPS == 3 .AND. MYID == 0) THEN
           WRITE(6, *)'Warning, calculations converged by 3 loops ...'
           WRITE(6, *)'OLDVAL:', OLDVAL(1:NEIG)
           WRITE(6, *)'EIGVAL:', EIGVAL(1:NEIG)
         ENDIF
!CYC-End
         RETURN
      ENDIF
!CYC-2024 End
!
! Find the maximum element of the last NNCV coefficients of unconverged
! eigenvectors. For those unconverged coefficients, put their indices
! to IND and find their number NNCV
!
      ICNT = 0
      DO I = 1, NEIG
!Rasa      IF (ICV(ISELEC(I)).EQ.0) THEN
!             ..Find coefficient and test for convergence
         ICUR = KPASS*ISELEC(I)
         TMAX = ABS(SVEC(ICUR))
         DO L = 1, NNCV - 1
            TMAX = MAX(TMAX,ABS(SVEC(ICUR-L)))
         END DO
         IF (TMAX < CRITC .AND. NLOOPS > 2) THEN
!                ..this  coefficient converged
            ICV(ISELEC(I)) = 1
         ELSE
!                ..Not converged. Add it to the list.
!Rasa -- start change
            ICV(ISELEC(I)) = 0
!Rasa -- end change
            ICNT = ICNT + 1
            IND(ICNT) = ISELEC(I)
            ROWLAST(ICNT) = TMAX
         ENDIF
!Rasa      ENDIF
      END DO

      NNCV = ICNT
      IF (NNCV == 0 .AND. NLOOPS > 2) THEN
         DONE = .TRUE.
         !IF (MYID == 0)     &
         !WRITE (6, *) 'converged by critc', KPASS
      ENDIF
!
! Sort the ROWLAST elements interchanging their indices as well
!
      DO I = 1, NNCV
         INDX = IDAMAX(NNCV - I + 1,ROWLAST(I),1)
         INCV(I) = IND(INDX+I-1)

         TEMP = ROWLAST(INDX+I-1)
         ROWLAST(INDX+I-1) = ROWLAST(I)
         ROWLAST(I) = TEMP
         ITEMP = IND(INDX+I-1)
         IND(INDX+I-1) = IND(I)
         IND(I) = ITEMP
      END DO

      TSTSEL = DONE
      RETURN
      END FUNCTION TSTSEL


!=======================================================================
      subroutine mgs_nrm(n, kp, new, scra, b)
!=======================================================================
!       Orthogonalizis the new vectors in B (from kp+1...kp+new)
!     to the previous kp B vectors and to themselves
!     using Modified Gram Schmidt. Then normalizes them.
!       The procedure is repeated twice.
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double
!     USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !use dgemv_I
      !use ddot_I
      !use dscal_I
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, INTENT(IN)  :: n
      integer, intent(in) :: kp
      integer, INTENT(IN)  :: new
      real(double), INTENT(IN)  :: scra(new)
      real(double), INTENT(INOUT)  :: b((kp + new)*n)
      real(kind(0.0d0)) :: ddot
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, newstart, kcur, k, jcur, j, mm
      real(double) :: dnm
!-----------------------------------------------

!
! MODIFIED GRAM-SCHMIDT  (twice)
!
        !if (myid.eq.0) then
      do i = 1, 2

         newstart = kp*n + 1
!
!  First record contribution from the kp vectors to the new ones.
!
         kcur = 1
         do k = 1, kp
            jcur = newstart
!            call dgemv ('T', n, new, 1.D0, b(jcur), n, b(kcur), 1, 0.D0, scra, &
!               1)
            CALL DGEMVMPI ('T', n, new, 1.D0, b(jcur), n, b(kcur), 1, &
                                0.D0, scra, 1)
            do j = 1, new
!           call daxpy(N,-scra(j),B(kcur),1,B(jcur),1)
               !b(jcur:jcur+n-1)=b(jcur:jcur+n-1)-scra(j)*b(jcur:jcur+n-1)
               do mm = 0, n - 1
                  b(jcur+mm) = b(jcur+mm) - scra(j)*b(kcur+mm)
               end do
               jcur = jcur + n
            end do
            kcur = kcur + n
         end do
!
!  Then orthogonalize the new ones among themselves.
!
         do k = 1, new
            jcur = kcur + n
!  The current vector should be normalized
!
            !dnm = ddot(n,b(kcur),1,b(kcur),1)
            CALL DDOTMPI(N,B(kcur),B(kcur),dnm)
            dnm = sqrt(dnm)
            call dscal (n, 1/dnm, b(kcur), 1)
!
            !call dgemv ('T', n, new - k, 1.D0, b(jcur), n, b(kcur), 1, 0.D0, &
            !   scra, 1)
            CALL DGEMVMPI ('T', n, new - k, 1.D0, b(jcur), n, b(kcur), &
                                1, 0.D0, scra, 1)
            do j = k + 1, new
!              call daxpy(N,-scra(j-k),B(kcur),1,B(jcur),1)
               do mm = 0, n - 1
                  b(jcur+mm) = b(jcur+mm) - scra(j-k)*b(kcur+mm)
               end do
               jcur = jcur + n
            end do
            kcur = kcur + n
         end do
      end do
        !end if

!        call MPI_BCAST(B,(kp+New)*n,MPI_DOUBLE_PRECISION,0,
!     :                 MPI_COMM_WORLD,ierr)

      return
      end subroutine mgs_nrm

!***********************************************************************
!cychen, 29/05/2020
      SUBROUTINE ddotmpi (n, x, y, ddot)
!
!-----------------------------------------------
!...Translated by Gediminas Gaigalas
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE MPI_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) :: x(n), y(n), ddot
      INTEGER :: n
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: dsum
      INTEGER :: i, nbar
!-----------------------------------------------

      ddot=0.0d0
      nbar = n/nprocs

      do i=myid*nbar+1, (myid+1)*nbar
         ddot=ddot+x(i)*y(i)
      enddo

      CALL MPI_Allreduce (ddot, dsum, 1, MPI_DOUBLE_PRECISION,        &
                                   MPI_SUM, MPI_COMM_WORLD, ierr)
      do i=nbar*nprocs+1, n
         dsum=dsum+x(i)*y(i)
      enddo
      ddot=dsum

      RETURN
      END SUBROUTINE ddotmpi

!***********************************************************************
!cychen, 29/05/2020
      subroutine dgemvmpi ( trans, m, n, alpha, a, lda, x, incx,     &
                         BETA, Y, INCY )
!
!cychen:
!DGEMV performs one of the matrix-vector operations by using mpi:
! y :=A*x, or y:=A'*x
!then, alpha and beta must be input as 1.0d0 and 0.0d0, respectively
!Also, for simplicit, incx and incy are modified to be both 1.
!
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!-----------------------------------------------
!...Translated by Gediminas Gaigalas
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
!     .. Scalar Arguments ..
      REAL(DOUBLE) :: ALPHA, BETA
      INTEGER      :: INCX, INCY, LDA, M, N
      CHARACTER*1  :: TRANS
!     .. Array Arguments ..
      REAL(DOUBLE) :: A( LDA, * ), X( * ), Y( * )
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: nbar
      REAL(DOUBLE) :: tm(m), tn
!-----------------------------------------------

!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.                         &
               .NOT.LSAME( TRANS, 'T' ).AND.                         &
               .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
!cychen      ELSE IF( INCX.EQ.0 )THEN
      ELSE IF( INCX.NE.1 )THEN
         INFO = 8
         if (myid.eq.0) print *, "INCX=",INCX
!cychen      ELSE IF( INCY.EQ.0 )THEN
      ELSE IF( INCY.NE.1 )THEN
         INFO = 11
         if (myid.eq.0) print *, "INCY=",INCY
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                               &
          ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )                &
         RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
!cychen      IF( BETA.NE.ONE )THEN
!cychen         IF( INCY.EQ.1 )THEN
!cychen            IF( BETA.EQ.ZERO )THEN
!cychen               DO 10, I = 1, LENY
!cychen                  Y( I ) = ZERO
!cychen   10          CONTINUE
!cychen            ELSE
!cychen               DO 20, I = 1, LENY
!cychen                  Y( I ) = BETA*Y( I )
!cychen   20          CONTINUE
!cychen            END IF
!cychen         ELSE
!cychen            IY = KY
!cychen            IF( BETA.EQ.ZERO )THEN
!cychen               DO 30, I = 1, LENY
!cychen                  Y( IY ) = ZERO
!cychen                  IY      = IY   + INCY
!cychen   30          CONTINUE
!cychen            ELSE
!cychen               DO 40, I = 1, LENY
!cychen                  Y( IY ) = BETA*Y( IY )
!cychen                  IY      = IY           + INCY
!cychen   40          CONTINUE
!cychen            END IF
!cychen         END IF
!cychen      END IF
!cychen      IF( ALPHA.EQ.ZERO )
!cychen     $   RETURN

      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
!         JX = KX
         IF( INCY.EQ.1 )THEN
            CALL DINIT(M, 0.d0, tm, 1)
            nbar = M / nprocs
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                do i = myid * nbar + 1, (myid + 1) * nbar
                   tm(i)=tm(i) + x(j) * A(i,j)
                enddo
               END IF
   60       CONTINUE
            CALL MPI_Allreduce (tm, y, m, MPI_DOUBLE_PRECISION,      &
                                   MPI_SUM, MPI_COMM_WORLD, ierr)
            do j = 1, N
               do i = nbar*nprocs + 1, M
                  y(i) = y(i) + x(j) * A(i,j)
               enddo
            enddo

         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
!
!        Form  y := alpha*A'*x + y.
!
         JY = KY
         IF( INCX.EQ.1 )THEN
            nbar = M / nprocs
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = myid * nbar + 1, (myid + 1) * nbar
                  TEMP = TEMP + A( I, J ) * X( I )
   90          CONTINUE
               CALL MPI_Allreduce (TEMP, tn, 1, MPI_DOUBLE_PRECISION,&
                                   MPI_SUM, MPI_COMM_WORLD, ierr)
               do i = nbar * nprocs + 1, M
                  tn = tn + A( I, J ) * X( I )
               enddo
               Y( J ) = tn
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMV .
!
      END subroutine dgemvmpi

!***********************************************************************
!cychen, 29/05/2020
!performs y := da*x + y
!For simplicit, modification only for incx=incy=1.

      subroutine daxpympi(n,da,dx,incx,dy,incy)
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
!-----------------------------------------------
!...Translated by Gediminas Gaigalas
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE mpi_C
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      REAL(DOUBLE) :: dx(1),dy(1),da
      INTEGER      :: n,incx,incy
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: i,ix,iy,m,mp1,nbar
      REAL(DOUBLE) :: tm(n)
!-----------------------------------------------

      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return

!     mpi code for both increments equal to 1
   20 CALL DINIT(n, 0.d0, tm, 1)
      nbar = n / nprocs
      do i = myid*nbar+1, (myid+1)*nbar
         tm(i) = dy(i) + da * dx(i)
      enddo
      CALL MPI_Allreduce (tm, dy, nprocs*nbar, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      do i = nprocs*nbar+1, n
         dy(i) = dy(i) + da*dx(i)
      enddo

!        m = mod(n,4)
!      if( m .eq. 0 ) go to 40
!      do 30 i = 1,m
!        dy(i) = dy(i) + da*dx(i)
!   30 continue
!      if( n .lt. 4 ) return
!   40 mp1 = m + 1
!      do 50 i = mp1,n,4
!        dy(i) = dy(i) + da*dx(i)
!        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
!        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
!        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
!   50 continue
      return
      end subroutine daxpympi
