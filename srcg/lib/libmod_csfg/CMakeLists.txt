add_library(mod_csfg STATIC
    csfg_decide_C.f90 
    symexpand_mod.f90  
    symmatrix_mod.f90  
    symmatrix_restart_C.f90
)
setup_fortran_modules(mod_csfg)
target_link_libraries_Fortran(mod_csfg PRIVATE mod 9290)
install(TARGETS mod_csfg DESTINATION lib/)
