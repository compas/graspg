LIB=mod_csfg
LAPACK=true
LIBRARIES="mod 9290"
FILES="
    csfg_decide_C.f90  
    symexpand_mod.f90  
    symmatrix_mod.f90  
    symmatrix_restart_C.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
