!=======================================================================
!                                                                      *
! Sorted Vk according to the original GRASP2018 orders                 *
!                                                                      *
! Written by Chongyang Chen,      Fudan University,          JAN 2024  *
!                                                                      *
!=======================================================================
      SUBROUTINE SORTOFFDV(NV, VCOEFF, LABCDK, LABP)

      USE vast_kind_param, ONLY: DOUBLE, BYTE
      USE buffer_C,        ONLY: NVCOEF, LABEL, COEFF

      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)    :: NV
      INTEGER(BYTE)          :: LABCDK(5,NV)
      INTEGER                :: LABP(NV)
      REAL(DOUBLE)           :: VCOEFF(NV)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER       :: IV, I, J
      INTEGER(BYTE) :: LABTMP(5)
      INTEGER       :: LPTMP
      REAL(DOUBLE)  :: VTMP

!-----------------------------------------------
      !!Sorted according to K
      DO I = 1, NV - 1
        DO J = I + 1, NV
          IF (LABCDK(5,I) .GT. LABCDK(5,J)) THEN
            LABTMP(1:5) = LABCDK(1:5,I)
            VTMP = VCOEFF(I)
            LPTMP = LABP(I)

            LABCDK(1:5,I) = LABCDK(1:5,J)
            VCOEFF(I) = VCOEFF(J)
            LABP(I) = LABP(J)

            LABCDK(1:5,J) = LABTMP(1:5)
            VCOEFF(J) = VTMP
            LABP(J) = LPTMP
          ENDIF
        ENDDO
      ENDDO

      !!Sorted according to LABP
      DO I = 1, NV - 1
        DO J = I + 1, NV
          IF (LABCDK(5,I) /= LABCDK(5,J)) EXIT
          IF (LABP(I) .GT. LABP(J)) THEN
            LABTMP(1:5) = LABCDK(1:5,I)
            VTMP = VCOEFF(I)
            LPTMP = LABP(I)

            LABCDK(1:5,I) = LABCDK(1:5,J)
            VCOEFF(I) = VCOEFF(J)
            LABP(I) = LABP(J)

            LABCDK(1:5,J) = LABTMP(1:5)
            VCOEFF(J) = VTMP
            LABP(J) = LPTMP
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE SORTOFFDV
