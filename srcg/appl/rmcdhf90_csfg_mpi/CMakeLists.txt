add_executable(rmcdhf_mpi_csfg
     rmcdhf_mem_C.f90 
     consis_I.f90 csfwgt_I.f90 dacon_I.f90 dampck_I.f90 dampor_I.f90 defcor_I.f90 
     dsubrs_I.f90 eigen_I.f90 endsum_I.f90 engoutgg_I.f90 estim_I.f90 fco_I.f90 gco_I.f90 getaldwt_I.f90 
     getoldwt_I.f90 in_I.f90 lagcon_I.f90 maxarr_I.f90 
     newe_I.f90 orbout_I.f90 
     orthy_I.f90 nxanyampi_I.f90 
     out_I.f90 outbnd_I.f90 prtrsl_I.f90 prwf_I.f90 setallcof_I.f90 settvcof_I.f90 setcof_cyc_I.f90
     setham_I.f90 setmcp_I.f90 setmcp_csfg_I.f90 setmix_I.f90 setsum_I.f90 setxuv_I.f90 setxv_I.f90 setxz_I.f90 solve_I.f90 
     strsum_I.f90 xpot_I.f90 ypot_I.f90 
     cofpotmpi_I.f90 getaldmpi_I.f90 getoldmpi_I.f90 getscdmpi_I.f90 
     maneigmpi_I.f90 matrixmpi_I.f90 newcompi_I.f90 
     improvmpi_I.f90  nxanyampi.f90 
     scfmpi_I.f90  itjpo_I.f90 ispar_I.f90
     setcslmpi_I.f90 setdbgmpi_I.f90 setlagmpi_I.f90 lodcsh2GG_I.f90 lodcslmpiGG_I.f90 
     systemmem_I.f90 systemfreemem_I.f90 setmcp_mem_I.f90 read1_mem_I.f90 read2_mem_I.f90 read3_mem_I.f90 
     read02_mem_I.f90 read03_mem_I.f90 

     read_TV.f90 twobody_DC.f90 savep_tcoeff.f90  savep_tlabcr.f90  savep_vk.f90  savep_vlabcr.f90 
     printlabt.f90 printlabv.f90 setarrs_csfg.f90 matrix_csfg.f90 maneig_csfg.f90 maneig_csfg_I.f90
     setham_csfg.f90 setham_csfg_I.f90 csfg_expand_findtype.f90 setmcp_csfg.f90
     locate.f90 locate_I.f90 match_labt.f90 match_labv.f90 locateall_I.f90 locateall.f90
     settvcof_csfg.f90 setallcof_csfg.f90 printdxy.f90
     spinangular1.f90
     spinangular11.f90
     spinangular12.f90
     spinangular13.f90
     spinangular14.f90
     spinangular15.f90
     spinangular22.f90
     spinangular23.f90
     spinangular24.f90
     spinangular25.f90
     spinangular2.f90
     spinangular33.f90
     spinangular34.f90
     spinangular35.f90
     spinangular3.f90
     spinangular44.f90
     spinangular45.f90
     spinangular4.f90
     spinangular55.f90
     spinangular5.f90

     consis.f90 csfwgt.f90 dacon.f90 dampck.f90 dampor.f90 defcor.f90 
     dsubrs.f90 eigen.f90 endsum.f90 engoutgg.f90 estim.f90 fco.f90 gco.f90 getaldwt.f90 
     getoldwt.f90 in.f90 lagcon.f90 maxarr.f90 newe.f90 orbout.f90 
     orthor.f90 orthy.f90 out.f90 outbnd.f90 prtrsl.f90 prwf.f90 setallcof.f90 settvcof.f90 setcof_cyc.f90 setdbg.f90 
     setham.f90 setham_cyc.f90 setmcp.f90 setmix.f90 setsum.f90 setxuv.f90 setxv.f90 setxz.f90 solve.f90 
     strsum.f90 xpot.f90 ypot.f90 
     cofpotmpi.f90 getaldmpi.f90 getoldmpi.f90 getscdmpi.f90 hmoutmpi.f90  
     improvmpi.f90 maneigmpi.f90 matrixmpi.f90 newcompi.f90  rscfmpivu.f90 scfmpi.f90 
     setcslmpi.f90 setdbgmpi.f90 setlagmpi.f90 lodcsh2GG.f90 lodcslmpiGG.f90 itjpo.f90 ispar.f90 
     systemmem.f90 systemfreemem.f90 setmcp_mem.f90 read1_mem.f90 read2_mem.f90 read3_mem.f90 
     read02_mem.f90 read03_mem.f90
)
target_link_libraries_Fortran(rmcdhf_mpi_csfg PUBLIC mod mcp90 9290 mpiu90 mod_csfg rang90_csfg 9290_csfg me_csfg dvd90mpi)
install(TARGETS rmcdhf_mpi_csfg DESTINATION bin/)
