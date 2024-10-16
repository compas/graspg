#!/bin/sh
#=======================================================================
#
# Beryllium-like C from MR-SD VV+CV+CC expansions
#
# Scripts are provied by the Fudan Atomic Spectroscopy Theory Group 
#                             ****(FASTG)****                     
#
#=======================================================================
# rmcdhf, rci, rtransition Calculations for Be-like C III
#
# Directory csfg / csf use the CSFG / CSF lists, respectively
#
# csfg: rcsfggenerate_csfg, rcsfgsplit_csfg,  rwfnestimate_csfg, 
#       rangular_csfg_mpi,  rmcdhf_csfg_mpi,  rsave_csfg,
#       rci_csfg_mpi,       jj2lsj_csfg,      rlevels_csfg
#
# csf : rcsfgenerate,       rcsfsplit,        rwfnestimate,      
#       rangular_mpi,       rmcdhf_mpi,       rsave,
#       rci_mpi,            jj2lsj_2024,      rlevels
#
#csf_CSFGExp : 
#                                             rwfnestimate,
#       rangular_mpi,       rmcdhf_mpi,       rsave,
#       rci_mpi,            jj2lsj_csfg
#       rcsfgexpand_csfg
#
# Both: rnucleus, # rtransition, ...
#
# !!!!! All the above programs should be found in the environmental variable ${PATH} !!!
# 
# Differences between the scripts
#
# vimdiff csfg/sh_files_c  csf/sh_files_c
# vimdiff csfg/sh_initial  csf/sh_initial
# vimdiff csfg/sh_scf      csf/sh_scf
#
#=======================================================================

cd csfg
./clean
./sh_csfg > ../out_csfg 2>&1
cd ../

cd csf
./clean
./sh_csf  > ../out_csf 2>&1
cd ../

cd csf_CSFGExp
./clean
./sh_csf  > ../out_csf_CSFGExp 2>&1
cd ../

#
# A quick look for GRASPG's speedup over GRASP2018 can be seen from 
# the information reported by 'time':
# csfg:
# real    30m6.073s
# user    454m55.288s
# sys     17m35.590s

# csf:
# real    601m49.449s
# user    9592m29.408s
# sys     11m8.433s

# csf_CSFGExp:
# real    450m45.355s
# user    7175m37.235s
# sys     10m31.594s

# Size and real (wall) time of rangular calculation: 
#   see the output log 'out_sh_scf', search string 'mcp000.30'
# Real (wall) times of rmcdhf and rci calculation:
#   see the correspoding output logs
#!!!
# vimdiff csfg/out_sh_scf csf/out_sh_scf
#!!!

# Wall-Time differences could be seen from the corresponding log-flies
# For example, please search string 'Wall time:' in the following comparisons:
#!!!
# vimdiff csfg/out_rmcdhf_n9 csf/out_rmcdhf_n9
# vimdiff csfg/out_rci_n9    csf/out_rci_n9
#!!!

# Final level energies:
#!!!
# The csfg energy for 1s(2).2p.3s_1P differ from the csf one by 2 cm^{-1}.
# vimdiff csfg/n9.lev csf/n9.lev 
# vimdiff csfg/n9CI.clev csf/n9CI.clev
#
# All the csfg energies differ from the csf_CSFGExp ones by less than 0.1 cm^{-1}
# vimdiff csfg/n9.lev csf_CSFGExp/n9.lev 
# vimdiff csfg/n9CI.clev csf_CSFGExp/n9CI.clev
#!!!

# There might be negligible differences between the GRASPG and GRASP2018 results. 
# The main reasons are the effects of inadequate numerical accuracies on:
#
#  a) exchanging summation orders of the contributions from different SLATER integral 
#     to construct the H matrixelement and Potentials

#  b) most spin-angular V-coefficients are taken from the diagonal-CSFG-pair. 
#     Ideally, they should be exactly same as those calculated from the 
#     off-diagonal-CSF-pair, if the two pairs differ from each other only 
#     by principle quantum numbers. In many cases, the calculated (T- and) V-ceofficients
#     are different by a factor of about 1.0d-14 / 1.0d-15. These differences may accumulate.

#  c) The most possible reason is that the different order of CSFs (basis)
#     may cause the differences, if the CC corelation is taken into account.
#     see the comparisons of script2, and script3/csfg and script3/csf_CSFGExp.
#     Also, see the comparisons of script4. 
