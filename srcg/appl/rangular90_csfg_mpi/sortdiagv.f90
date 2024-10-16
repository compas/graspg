!=======================================================================
! Sort the V-coefficients according to the origianl GRASP2018 orders   *
!                                                                      *
! Written by Chongyang Chen,      Fudan University,          Jan 2024  *
!                                                                      *
!=======================================================================
      SUBROUTINE SORTDIAGV(NV, VCOEFF, LABCDK, LABP)

      USE vast_kind_param, ONLY: DOUBLE, BYTE

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
      INTEGER       :: LPTMP, LSTMP
      REAL(DOUBLE)  :: VTMP

      INTEGER(BYTE) :: JARR(4), NDIFFJ
      INTEGER       :: ISORT(NV) 
!-----------------------------------------------
      DO IV = 1, NV
        JARR = -1
        JARR(1) = LABCDK(1,IV)
        NDIFFJ = 1
        DO I = 2, 4
          IF (COUNT(JARR == LABCDK(I,IV)) == 0) THEN
            NDIFFJ = NDIFFJ + 1
            JARR(NDIFFJ) = LABCDK(I,IV) 
          ENDIF
        ENDDO

        SELECT CASE (NDIFFJ)

        CASE (1)
          ISORT(IV) = 1

        CASE (2)
          IF (LABCDK(3,IV) < LABCDK(4,IV)) THEN
            ISORT(IV) = 2
          ELSE
            ISORT(IV) = 3
          ENDIF

        CASE (3)
          IF (LABCDK(3,IV) < LABCDK(4,IV)) THEN
            ISORT(IV) = 4
          ELSE
            ISORT(IV) = 5
          ENDIF 

        CASE (4)
          IF (LABCDK(3,IV) < LABCDK(4,IV)) THEN
            ISORT(IV) = 6
          ELSE
            ISORT(IV) = 7
          ENDIF 
        END SELECT 
      ENDDO

      !!Sorted according to K
      DO I = 1, NV - 1
        DO J = I + 1, NV
          IF (LABCDK(5,I) .GT. LABCDK(5,J)) THEN
            CALL EXCHANGE_DATA(NV, I, J, VCOEFF, LABCDK, LABP, ISORT)
          ENDIF
        ENDDO
      ENDDO

      !!Sorted according to ISORT
      DO I = 1, NV - 1
        DO J = I + 1, NV
          IF (LABCDK(5,I) /= LABCDK(5,J)) EXIT
          IF (ISORT(I) .GT. ISORT(J)) THEN
            CALL EXCHANGE_DATA(NV, I, J, VCOEFF, LABCDK, LABP, ISORT)
          ENDIF
        ENDDO
      ENDDO

      !!Sorted according to LABP
      DO I = 1, NV - 1
        DO J = I + 1, NV
          IF (ISORT(I) /= ISORT(J)) EXIT
          IF (LABCDK(5,I) /= LABCDK(5,J)) EXIT
          IF (LABP(I) .GT. LABP(J)) THEN
            CALL EXCHANGE_DATA(NV, I, J, VCOEFF, LABCDK, LABP, ISORT)
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE SORTDIAGV


!======================================================================= 
      SUBROUTINE EXCHANGE_DATA(NV, I, J, VCOEFF, LABCDK, LABP, ISORT)
      USE vast_kind_param, ONLY: DOUBLE, BYTE

      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)    :: NV, I, J
      REAL(DOUBLE)           :: VCOEFF(NV)
      INTEGER(BYTE)          :: LABCDK(5,NV)
      INTEGER                :: LABP(NV)
      INTEGER                :: ISORT(NV)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER(BYTE):: LABTMP(5)
      INTEGER      :: LPTMP, LSTMP
      REAL(DOUBLE) :: VTMP

      LABTMP(1:5) = LABCDK(1:5,I)
      VTMP = VCOEFF(I)
      LPTMP = LABP(I)
      LSTMP = ISORT(I)

      LABCDK(1:5,I) = LABCDK(1:5,J)
      VCOEFF(I) = VCOEFF(J)
      LABP(I) = LABP(J)
      ISORT(I) = ISORT(J)
 
      LABCDK(1:5,J) = LABTMP(1:5)
      VCOEFF(J) = VTMP
      LABP(J) = LPTMP
      ISORT(J) = LSTMP

      END SUBROUTINE EXCHANGE_DATA
