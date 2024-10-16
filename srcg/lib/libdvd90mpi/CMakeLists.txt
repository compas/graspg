add_library(dvd90mpi STATIC
    adds_I.f90
    dvdrvr_I.f90
    dvdson_I.f90
    initdvd_I.f90
    mgs_nrm_I.f90
    multbc_I.f90
    newvec_I.f90
    ovflow_I.f90
    tstsel_I.f90
    dvdson.f90
    gdvd.f90
    gdvd_I.f90
    iniest.f90
    iniest_I.f90
)
setup_fortran_modules(dvd90mpi)
target_link_libraries_Fortran(dvd90mpi PRIVATE mod mpiu90)
install(TARGETS dvd90mpi DESTINATION lib/)
