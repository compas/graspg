EXE=rmixaccumulate_csfg_basic
LIBRARIES=
FILES="
   rmixaccumulate_csfg.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
