!=======================================================================
! Print buffer_C to figure out the relationship of the VCOEFF (two -   *
! electron interaction and BREIT interaction between different IC -    *
! IR pairs consisting of symbolic-orbital.                             *
! Written by Chongyang Chen,      Fudan University,          Jan 2022  *
!=======================================================================
      SUBROUTINE PRINTLABELV(IFLAG, ICW, IRW, NTXX, MTYPE, IA, IB, TSHELL1)
      USE symmatrix_mod, ONLY: MAP
      USE buffer_C,      ONLY: NVCOEF, LABEL, COEFF
      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER    :: IFLAG, ICW, IRW, NTXX, IA, IB, MTYPE
      REAL*8     :: TSHELL1

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER(LEN=256) :: STREAM
      CHARACTER(LEN=2)   :: CHAR2
      INTEGER            :: IW, J
!-----------------------------------------------
      IF (IFLAG == 0) THEN
        CHAR2 = 'BE'  ! Begin IC0 - IC, IR0 - IR block
      ELSEIF (IFLAG == 1) THEN
        CHAR2 = 'OB'  ! ONE - BODY
      ELSEIF (IFLAG == 2) THEN
        CHAR2 = 'TB'  ! TWO - BODY ELECTRON INTERACTION
      ELSEIF (IFLAG == 3) THEN
        CHAR2 = 'BR'  ! BREIT INTERACTION
      ENDIF
      WRITE(STREAM,'(A2,2x,A6,I3,2X,A4,I8,2X,A8,I8,2X,A4,I8,   &
                     2X,A8,I8,2X,A14,I3)')   &
            CHAR2,                                    &
            'MTYPE=', MTYPE,                          &
            'ICW=', ICW, 'ICW*3+5=', ICW*3+5,         &
            'IRW=', IRW, 'IRW*3+5=', IRW*3+5,         &
            'MatrixblockXX=', NTXX                    
      ! Write the data into different files for different NTXX
      IW = 1000+NTXX 

      ! Write the data into one file 
      !IW = 2543  

      IF (IFLAG == 0) WRITE(IW,*)
      WRITE(IW, '(A)')TRIM(STREAM)
      IF (IFLAG == 0) WRITE(IW,'(4X,A4,I8,2X,A8,I8, &
                                 2X,A4,I8,2X,A8,I8)')       &
        'IR0=',MAP(IRW-1)+1, 'IR0*3+5=',(MAP(IRW-1)+1)*3+5, &
        'IRE=',MAP(IRW  )  , 'IRE*3+5=',(MAP(IRW  )  )*3+5   

      IF (IFLAG == 0) WRITE(IW,'(4X,A4,I8,2X,A8,I8, &
                                 2X,A4,I8,2X,A8,I8)')       &
        'IC0=',MAP(ICW-1)+1, 'IC0*3+5=',(MAP(ICW-1)+1)*3+5, &
        'ICE=',MAP(ICW  )  , 'ICE*3+5=',(MAP(ICW  )  )*3+5   
          
      IF (IFLAG == 1) THEN
       WRITE(IW,'(4X, 3HIA=, I4, 2X, 3HIB=, I4, 2X, 8HTSHELL1=,1PE12.3)') &
                IA, IB, TSHELL1 
      ELSEIF (IFLAG == 2) THEN
       write(IW,'(4x,7hIVCOEF=,1000i11)')(J, J=1,NVCOEF)
       write(IW,'(4x,7hLabel1=,1000i11)')(Label(1,J), J=1,NVCOEF)
       write(IW,'(4x,7hLabel2=,1000i11)')(Label(2,J), J=1,NVCOEF)
       write(IW,'(4x,7hLabel3=,1000i11)')(Label(3,J), J=1,NVCOEF)
       write(IW,'(4x,7hLabel4=,1000i11)')(Label(4,J), J=1,NVCOEF)
       write(IW,'(4x,7hLabel5=,1000i11)')(Label(5,J), J=1,NVCOEF)
       write(IW,'(4x,7hCoeff =,1000(1pe11.3))')(Coeff(J), J=1,NVCOEF)
      ELSEIF (IFLAG == 3) THEN
       write(IW,'(4x,7hIVCOEF=,1000i11)')(J, J=1,NVCOEF)
       write(IW,'(4x,7hLabel1=,1000i11)')(Label(1,J), J=1,NVCOEF)
       write(IW,'(4x,7hLabel2=,1000i11)')(Label(2,J), J=1,NVCOEF)
       write(IW,'(4x,7hLabel3=,1000i11)')(Label(3,J), J=1,NVCOEF)
       write(IW,'(4x,7hLabel4=,1000i11)')(Label(4,J), J=1,NVCOEF)
       write(IW,'(4x,7hLabel5=,1000i11)')(Label(5,J), J=1,NVCOEF)
       write(IW,'(4x,7hLabel6=,1000i11)')(Label(6,J), J=1,NVCOEF)
       write(IW,'(4x,7hCoeff =,1000(1pe11.3))')(Coeff(J), J=1,NVCOEF)
      ENDIF
      END SUBROUTINE PRINTLABELV
