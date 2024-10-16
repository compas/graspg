!**********************************************************************
!                                                                     *      
      SUBROUTINE TRANSFER_cyc(IC, IR)
!                                                                     *
!  This subroutine transfers the non-zero data in B(M2,N) along       *
!  with the corresponding ACCUMULATED row position to the arrays      *
!  A(M1,N) and IROW(M1,N).                                            *
!  IACOUNT(N) is a counter array that keeps track of                  *
!  the number of non-zero elements in A(M1,N) for each column         *
!                                                                     *
!  For a documentation see xxxx                                       *         
!                                                                     *
!  Per Jönsson, Malmö University,                April 2017         *
!                                                                     *
!     Modified by cychen                                              *
!     Record the non-zero matrixelements in EMTBLOCK into EMTSYM, also*
!     set IROWSYM values                                              *
!     Last modification by Chongyang Chen   Nov 2023                  *
!********************************************************************** 
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------        
      USE decide_C,        ONLY: LFORDR
      use symmatrix_mod,   ONLY: CUTOFF, IRTOT, LICCUT, NTYPE,         &
                                 EMTBLOCK, NELCSYM, IROWSYM, EMTSYM
      use csfg_decide_C,   ONLY: NDISCARDBR, LDISCARDBR,               &
                                 NDISBR, LDISBR, LSYMD
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER IC, IR
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER I, J, M, N, IR0
! CUTOFF is set within symmatrix_mod.f90
      !REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-20
!-----------------------------------------------
!

      M=NTYPE(2,IR)
      N=NTYPE(2,IC)
      IF ((.NOT.LFORDR) .OR. (.NOT.LICCUT) .OR. IC.NE.IR .OR. LSYMD)THEN
        ! Include all of the non-zero matrixelements
        ! LSYMD: Diagonal-Block interaction
        DO J = 1, N
           IR0 = IRTOT
           DO I = 1, M
              IR0 = IR0 + 1
              IF (DABS(EMTBLOCK(I,J)) .GT. CUTOFF) THEN
                 NELCSYM(J) = NELCSYM(J) + 1
                 EMTSYM(NELCSYM(J),J) = EMTBLOCK(I,J)
                 IROWSYM(NELCSYM(J),J)=IR0
              ENDIF
           ENDDO
        ENDDO
        RETURN
      ELSE
      ! Discard the CI within the Diagonal-Block
        DO I = 1, M
          J = I
          NELCSYM(J) = NELCSYM(J) + 1
          EMTSYM(NELCSYM(J),J) = EMTBLOCK(I,J)
          IROWSYM(NELCSYM(J),J)= IRTOT + I
        ENDDO
        RETURN
      ENDIF

      RETURN
      END SUBROUTINE TRANSFER_cyc
