      MODULE Interact_CSF_I
      INTERFACE
!
      SUBROUTINE Interact_CSF(JA,JB,ICOLBREI,int_CSF)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER, INTENT(IN) :: JA,JB,ICOLBREI
      INTEGER, INTENT(OUT) :: int_CSF
      END SUBROUTINE
      END INTERFACE
      END MODULE