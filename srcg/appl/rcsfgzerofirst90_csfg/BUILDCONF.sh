EXE=rcsfzerofirst_omp
LIBRARIES="9290 mod"
FILES="
        symzf_mod.f90
        removeblank_I.f90 
        lodcsl_Part_I.f90 
        lodcsl_Zero_I.f90 
        lodcsl_Full_I.f90 
        set_CSF_ZFlist_I.f90 
        set_CSF_number_I.f90 
        lodcsl_Part.f90 
        lodcsl_Zero.f90 
        lodcsl_Full.f90 
        set_CSF_ZFlist.f90 
        set_CSF_number.f90 
        removeblank.f90 
        set_CSF_numberF.f90 
        RCSFzerofirst.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
