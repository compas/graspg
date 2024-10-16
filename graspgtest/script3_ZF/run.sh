#!/bin/sh
#=======================================================================
#
# Beryllium-like C from MR-SD VV+CV+CC expansions
#
# Scripts are provied by the Fudan Atomic Spectroscopy Theory Group 
#                             ****(FASTG)****                     
#
# test condensation by rcsfginteract_csfg
# test zero-first approximation by rcsfgzerofirst_csfg
#
#=======================================================================
# rmcdhf, rci, rtransition Calculations for Be-like C III
#
# Directory csfg / csf use the CSFG / CSF lists, respectively
#
# csfg: rcsfggenerate_csfg, rcsfgsplit_csfg,  rwfnestimate_csfg, 
#       rangular_csfg_mpi,  rmcdhf_csfg_mpi,  rsave_csfg,
#       rci_csfg_mpi,       jj2lsj_csfg,      rlevels_csfg
#       rcsfginteract_csfg, rcsfgzerofirst_csfg
#
# csf : rcsfgenerate,       rcsfsplit,        rwfnestimate,      
#       rangular_mpi,       rmcdhf_mpi,       rsave,
#       rci_mpi,            jj2lsj_2024,      rlevels
#       rcsfinteract,       rcsfzerofirst
#
#csf_CSFGExp : 
#       rcsfggenerate_csfg, rcsfgsplit_csfg,  rwfnestimate,
#       rangular_mpi,       rmcdhf_mpi,       rsave,
#       rci_mpi,            jj2lsj_csfg
#       rcsfgexpand_csfg
#       rcsfinteract,       rcsfzerofirst
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
# real    15m51.340s
# user    230m41.725s
# sys     18m25.111s

# csf:
# real    211m39.617s
# user    3363m42.710s
# sys     11m23.404s

# csf_CSFGExp:
# real    212m40.297s
# user    3383m12.844s
# sys     11m43.971s

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
# vimdiff csfg/n9.lev csf/n9.lev 

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

