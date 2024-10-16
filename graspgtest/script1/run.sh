#=======================================================================
#
# Boron-like Mo from MR-SD VV+CV+CC expansions
#
# Scripts are provied by the Fudan Atomic Spectroscopy Theory Group 
#                             ****(FASTG)****                     
#
#=======================================================================
# rmcdhf, rci, rtransition Calculations for B-like Mo.
#
# Scripts are modified from the original ones of GRASP2018 package.
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

# A quick look for GRASPG's speedup over GRASP2018 can be seen from
# the information reported by 'time' in out_csfg and out_csf:
# real    4m30.842s
# user    16m59.065s
# sys     0m15.521s

# real    37m46.294s
# user    149m16.009s
# sys     0m42.522s

# Size and real (wall) time of rangular calculation: 
#   see the output log 'out_sh_scf' in directories csf/ and csfg, search
#   string 'mcp000.30'
# Real (wall) times of rmcdhf and rci calculation:
#   see the correspoding output logs
#!!!
# vimdiff csfg/out_sh_scf csf/out_sh_scf
#!!!

# Wall-Time differences could be seen from the corresponding log-flies
# For example, please search string 'Wall time:' in the following comparisons:
#!!!
# vimdiff csfg/out_rmcdhf_odd6       csf/out_rmcdhf_odd6
# vimdiff csfg/out_rmcdhf_even6      csf/out_rmcdhf_even6
# vimdiff csfg/out_rci_odd           csf/out_rci_odd 
# vimdiff csfg/out_rci_even          csf/out_rci_even
#!!!

# Final level energies:
#!!!
# vimdiff csfg/n6.lev    csf/n6.lev 
# vimdiff csfg/n6CI.clev csf/n6CI.clev
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

