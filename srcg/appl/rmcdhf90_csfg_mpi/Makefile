.SUFFIXES: .f90 .mod
# This code deal with very high-n orbitals, the weighted average energy varies slightly
EXE = ${GRASP}/bin/rmcdhf_csfg_mpi
#BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
#BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ${GRASP}/src/lib
MODDIR = ${SRCLIBDIR}/libmod
MODL92 = ${SRCLIBDIR}/lib9290
MODLMPIU90 = ${SRCLIBDIR}/mpi90
MODLMCP90 = ${SRCLIBDIR}/libmcp90
MODRAD   = ${SRCLIBDIR}/librang90

MODDVD = ${GRASP}/srcg/lib/libdvd90mpi
#MODDVD = ${GRASP}/srcg/lib/libdvd90
MODCSFG9290 = ${GRASP}/srcg/lib/lib9290_csfg
MODCSFGME   = ${GRASP}/srcg/lib/libme_csfg
MODCSFGMOD   = ${GRASP}/srcg/lib/libmod_csfg
GRASPLIBS =-lmod -l9290 -lmcp90 -lmpiu90 -ldvdmpi -lme_csfg -l9290_csfg -lmod_csfg -lrang90
LAPACK90_LIB = -llapack -lblas

APP_LIBS = -L${GRASPLIB} ${GRASPLIBS}

#   Define data types
VASTO = ${MODDIR}/vast_kind_param_M.o

APP_OBJ= \
     report_mem.o \
     csfg_iqa.o  \
     dsubrs.o  \
     dsubrs_I.o  \
     calddrs.o \
     calddrs_I.o \
     outbnd.o  \
     outbnd_I.o  \
     newe.o  \
     newe_I.o  \
     prwf.o  \
     prwf_I.o  \
     setxuv.o  \
     setxuv_I.o  \
     setxv.o  \
     setxv_I.o  \
     in.o  \
     in_I.o  \
     estim.o  \
     estim_I.o  \
     eigen.o  \
     eigen_I.o  \
     out.o  \
     out_I.o  \
     setxz.o  \
     setxz_I.o  \
     solve.o  \
     solve_I.o  \
     prtrsl.o  \
     prtrsl_I.o  \
     setcof_cyc.o  \
     setcof_cyc_I.o  \
     ypot.o  \
     ypot_I.o  \
     xpot.o  \
     xpot_I.o  \
     lagcon.o  \
     lagcon_I.o  \
     dacon.o  \
     dacon_I.o  \
     engoutgg.o  \
     engoutgg_I.o  \
     getoldwt.o  \
     getoldwt_I.o  \
     getoldmpi.o  \
     getoldmpi_I.o  \
     getaldwt.o  \
     getaldwt_I.o  \
     getaldmpi.o  \
     getaldmpi_I.o  \
     \
     cofpotmpi.o  \
     cofpotmpi_I.o  \
     consis.o  \
     consis_I.o  \
     csfg_expand_findtype.o  \
     csfwgt.o  \
     csfwgt_I.o  \
     dampck.o  \
     dampck_I.o  \
     dampor.o  \
     dampor_I.o  \
     defcor.o  \
     defcor_I.o  \
     endsum.o  \
     endsum_I.o  \
     getscdmpi.o  \
     getscdmpi_I.o  \
     hmoutmpi.o  \
     hmoutmpi_I.o  \
     \
     setlagmpi.o  \
     setlagmpi_I.o  \
     orthy.o  \
     orthy_I.o  \
     \
     improvmpi.o  \
     improvmpi_I.o  \
     ispar.o  \
     ispar_I.o  \
     itjpo.o  \
     itjpo_I.o  \
     locateall.o  \
     locateall_I.o  \
     locate.o  \
     locate_I.o  \
     lodcsh2GG.o  \
     lodcsh2GG_I.o  \
     lodcslmpiGG.o  \
     lodcslmpiGG_I.o  \
     maneig_csfg.o  \
     maneig_csfg_I.o  \
     maxarr.o  \
     maxarr_I.o  \
     newcompi.o  \
     newcompi_I.o  \
     nxanyampi.o  \
     nxanyampi_I.o  \
     onebody_DC.o  \
     orbout.o  \
     orbout_I.o  \
     orthor.o  \
     orthor_I.o  \
     printdxy.o  \
     read_TV.o  \
     setcof_nda_csfg.o  \
     setcof_xy_csfg.o  \
     setcof_xy_csfg_I.o \
     setcof_nxy_csfg.o  \
     setarrs_csfg.o  \
     setcslmpi.o  \
     setcslmpi_I.o  \
     setdbg.o  \
     setdbg_I.o  \
     setdbgmpi.o  \
     setdbgmpi_I.o  \
     setham_pot.o  \
     setham_pot_I.o  \
     setham_cyc.o  \
     setmcp_csfg.o  \
     setmcp_csfg_I.o  \
     setmix.o  \
     setmix_I.o  \
     setsum.o  \
     setsum_I.o  \
     matrix_csfg.o  \
     matrix_csfg_I.o  \
     transfer_csfg.o  \
     spinangular11.o  \
     spinangular12.o  \
     spinangular13.o  \
     spinangular14.o  \
     spinangular15.o  \
     spinangular1.o  \
     spinangular22.o  \
     spinangular23.o  \
     spinangular24.o  \
     spinangular25.o  \
     spinangular2.o  \
     spinangular33.o  \
     spinangular34.o  \
     spinangular35.o  \
     spinangular3.o  \
     spinangular44.o  \
     spinangular45.o  \
     spinangular4.o  \
     spinangular55.o  \
     spinangular5.o  \
     strsum.o  \
     strsum_I.o  \
     twobody_DC.o  \
     vector_xya.o \
     genintrk_csfg.o  \
     genintrk_csfg_I.o  \
     genintrkwrap_csfg.o  \
     genintrkwrap_csfg_I.o  \
     scfmpi_I.o  \
     scfmpi.o    \
     rscfmpivu.o  


$(EXE): $(APP_OBJ)
	$(FC_MPI) -o $(EXE) $(FC_MPILD) $(APP_OBJ) $(APP_LIBS) \
		$(APP_LIBS) ${LAPACK_LIBS} 


.f90.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -I $(MODDIR) -I ${MODL92} -I $(MODLMPIU90) -I ${MODDVD} \
               -I ${MODCSFG9290} -I ${MODCSFGME} -I ${MODCSFGMOD} -I ${MODRAD} -o $@

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
        

