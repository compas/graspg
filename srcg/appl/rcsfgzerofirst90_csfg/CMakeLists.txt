set(CMAKE_Fortran_FLAGS "-fcray-pointer -fallow-argument-mismatch ${OpenMP_Fortran_FLAGS}")
add_executable(rcsfzerofirst_omp
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
)
target_link_libraries_Fortran(rcsfzerofirst_omp PUBLIC 9290 mod)
install(TARGETS rcsfzerofirst_omp DESTINATION bin/)
