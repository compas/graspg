!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXBLOCK13(ICSAV,IRSAV,NCOEC,INCOR,NCTEC,INC2,     &
                 NMCBP,NCORE,ELSTO)
!                                                                      *
!   This subroutine calls onescalar and computes one electron          *
!   matrix elements when IC and IR are of type 1 and 3                 *
!                                                                      *
!   Written by Per Jönsson & Kai Wang                        May 2017 *
!                                                                      *
!***********************************************************************
!   Modified by CHONG-YANG CHEN                              JUNE 2020 *
!   Last modification by C. Y. Chen                          Dec  2023 *
!***********************************************************************
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP
      use symmatrix_mod
      USE buffer_C,   ONLY: NVCOEF, LABEL, COEFF
      USE debug_C,    ONLY: IBUG1
      USE decide_C
      USE def_C
      USE orb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
      EXTERNAL BREID,CORD
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: ICSAV,IRSAV,NCOEC,INCOR,NCTEC,INC2,NMCBP,NCORE
      REAL(DOUBLE)        :: ELSTO
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!   Matrix elements smaller than CUTOFF are not accumulated
!
!      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-20
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: EMTTMP,ATWINV,ENONSYM,TCOEFF,VCOEFF
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER :: I,J,IC,IR,IV,N,NORBCOLL,NORBCOLU
      INTEGER :: LABVI,LABVJ,LABVK,LABVL
!-----------------------------------------------------------------------
      IF (LTRANSPOSE) THEN
       IC = IRSAV
       IR = ICSAV
      ELSE
       IC = ICSAV
       IR = IRSAV
      ENDIF

      ATWINV = 1.D0/EMN
      IBUG1 = 0

!   Set all matrix elements to zero
      !EMTBLOCK = 0.0d0
!
!   There is no onebody contribution for TYPE 1 - 3, 4, and 5
!
!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation; the latter
!   is computed first because the orbital indices may be
!   permuted by RKINTC
!
      NVCOEF = 0
      NORBCOLL = NTYPE(4,IC) - NTYPE(3,IC) + 1 
      NORBCOLU = NTYPE(6,IC) - NTYPE(5,IC) + 1 
!
      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
      IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
      ENONSYM = 0.0D0
      DO IV = 1, NVCOEF
        VCOEFF = COEFF(IV)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC) 

! Determine the position of the symmetry-ordered orbitals
          LABV=LABEL(1:6, IV)
          CALL ANALABV(IC, IR, IV)

          IF (NSYMCR .EQ. 0) THEN 
! There is no symmetry-ordered orbs for this interact pairs. It is common for
! those generated by the IC and IR symmetry-ordered CSF. It should not occur 
! for type 1 - 2, 3, 4, 5 matrix elements.
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')             &
                 'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***0***'

! Matrixelement between Type 1 - 3
! Loop for symmetry-ordered-orbitals
          ELSEIF (NSYMCR .EQ. 1) THEN 
            WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')             &
                 'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
            STOP 'Warning!!! NSYMCR ERROR ***1***'

          ELSEIF (NSYMCR .EQ. 2) THEN
            N = 0  ! Column index
            DO I = 1, NORBCOLL
              ! Replace the symmetry-ordered-orbital within LABV, the position
              ! is determied by IPSym(1) and IPSym(2)
              LABV(IPSym(1)) = NTYPE(3, IC) + I - 1
              DO J = 1, NORBCOLU
                N = N + 1 
                LABV(IPSym(2)) = NTYPE(5, IC) + J - 1
                CALL twobody_DC_SMS(ATWINV, VCOEFF, EMTTMP)
                EMTBLOCK(1,N) = EMTBLOCK(1,N) + EMTTMP
              ENDDO
            ENDDO

          ELSEIF (NSYMCR .EQ. 3) THEN
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***3***'

          ELSEIF (NSYMCR .EQ. 4) THEN
           WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')              &
                'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
           STOP 'Warning!!! NSYMCR ERROR ***4***'
          ENDIF 
        ENDIF
      ENDDO

      IBUG1 = 0
!
!   Accumulate the contribution from the two-electron
!   transverse interaction operator
      IF (LTRANS .AND. (INC2.EQ.1)) THEN
        NVCOEF = 0
        CALL RKCO_GG (IC, IR, BREID, 1, 2)
        IF (NVCOEF .GT. 0) LTRANSFER = .TRUE.
        DO 10 IV = 1, NVCOEF
          VCOEFF = COEFF(IV)
          IF (DABS (VCOEFF) .GT. CUTOFF) THEN
            NMCBP = NMCBP + NTYPE(2,IC) 
! Determine the position of the symmetry-ordered orbitals
            LABV = LABEL(1:6,IV)
            CALL ANALABV(IC, IR, IV)

            IF (NSYMCR .EQ. 0) THEN
! There is no symmetry-ordered orbs for this interact pairs. It is common for
! those generated by the IC and IR symmetry-ordered CSF. There is no such
! contributions for type 1 - 2, 3, 4, 5 matrix elements.
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')           &
                   'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
              STOP 'Warning!!! NSYMCR ERROR ***0***'

! Matrixelement between Type 1 - 2, 5, Loop for symmetry-ordered-orbitals
            ELSEIF (NSYMCR .EQ. 1) THEN
              WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')           &
                   'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
              STOP 'Warning!!! NSYMCR ERROR ***1***'

            ELSEIF (NSYMCR .EQ. 2) THEN
              IF (NSYMCOL.NE.1 .OR. NSYMROW.NE.1) THEN
                WRITE(*,*)'IC,IR,IV,LABV=',IC,IR,IV,LABV(1:6)
                STOP 'Unexpected NSYMCOL.NE.1 .OR. NSYMROW.NE.1 in B13 '
              ENDIF
              N = 0 ! Column index
              DO I = 1, NORBCOLL
                LABVI = NTYPE(3,IC) + I - 1
                DO J = 1, NORBCOLU
                  LABVJ = NTYPE(5,IC) + J - 1
                  N = N + 1
              ! Replace the symmetry-ordered-orbital within LABV, the position
              ! is determied by IPSym(1) and IPSym(2)
                  IF (LABEL(IPSym(1),IV).EQ.NTYPE(4,IC)) THEN
                    !LABV(IPSym(1)) = NTYPE(3,IC) + I - 1
                    !LABV(IPSym(2)) = NTYPE(5,IC) + J - 1
                    LABV(IPSym(1)) = LABVI
                    LABV(IPSym(2)) = LABVJ
                  ELSE
                    !LABV(IPSym(2)) = NTYPE(3,IC) + I - 1
                    !LABV(IPSym(1)) = NTYPE(5,IC) + J - 1
                    LABV(IPSym(1)) = LABVJ
                    LABV(IPSym(2)) = LABVI
                  ENDIF
                  CALL breit_d(IC,IR,IV,VCOEFF,NCORE,ELSTO,EMTTMP)
                  EMTBLOCK(1,N) = EMTBLOCK(1,N) + EMTTMP
                ENDDO
              ENDDO

            ELSEIF (NSYMCR .EQ. 3) THEN
             WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')            &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
             STOP 'Warning!!! NSYMCR ERROR ***3***'

            ELSEIF (NSYMCR .EQ. 4) THEN
             WRITE(*,'(a3, i7, 2x, a3, i7, 2x, a5, 5i3)')            &
                  'IC=', IC, 'IR=', IR, 'LABV=', LABV(1:5)
             STOP 'Warning!!! NSYMCR ERROR ***4***'
            ENDIF

          ENDIF
   10   CONTINUE
      ENDIF
!cyc  Transfer EMTBLOCK to EMTSYM 
      IF (LTRANSPOSE) EMTBLOCK=TRANSPOSE(EMTBLOCK)
      IF (LTRANSFER) Call transfer_cyc(ICSAV, IRSAV)

      RETURN
      END SUBROUTINE MATRIXBLOCK13