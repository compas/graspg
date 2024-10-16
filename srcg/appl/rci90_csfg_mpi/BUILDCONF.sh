EXE=rci_mpi_csfg
LIBRARIES="mpiu90 9290 9290_csfg mod rang90_csfg mod_csfg me_csfg dvd90mpi"
LAPACK=true
ISMPI=true
FILES="
      setham_gg.f90 
      setham_gg_I.f90 
      genmat.f90 
      genmat_I.f90 
      getcid.f90 
      getcid_I.f90 
      lodres.f90 
      lodres_I.f90 
      lodstate_sym.f90 
      lodstate_sym_I.f90 
      maneigmpi.f90 
      maneigmpi_I.f90 
      matrix.f90 
      matrix_I.f90 
      setres.f90 
      setres_I.f90 

      rci90mpi.f90 
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
