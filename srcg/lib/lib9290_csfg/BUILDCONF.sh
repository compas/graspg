LIB=9290_csfg
LIBRARIES="mod 9290 mpiu90 mod_csfg"
FILES="
   analabv.f90
   analabv_I.f90
   csfg_ispar.f90
   csfg_ispar_I.f90
   csfg_itjpo.f90
   csfg_itjpo_I.f90
   csfg_jcup.f90
   csfg_jcup_I.f90
   csfg_setqna.f90
   csfg_setqna_I.f90
   fictious_csf.f90
   findtype.f90
   indexsym.f90
   indexsym_I.f90
   printfictcsf.f90
   printlabelv.f90
   rcsfsymexpand.f90
   rcsfsymexpand_I.f90
   setlaborb.f90
   setlaborb_I.f90
"
generate-makefile > ${MAKEFILE}
generate-cmakelists > ${CMAKELISTSTXT}
