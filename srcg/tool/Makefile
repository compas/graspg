########################################################################
# Use the original GRASP libraries
LIBS=
FC_MODULES=
LAPACK_LIBS=
########################################################################
FC=gfortran

install: all

all:  ${GRASP}/bin/rsave_csfg ${GRASP}/bin/rcsfgexpand_csfg ${GRASP}/bin/rcsfgsplit_csfg \
      ${GRASP}/bin/rdistHmatrix_csfg ${GRASP}/bin/rcsfgextend_csfg \
      ${GRASP}/bin/rcsfgblocksplit_csfg ${GRASP}/bin/rlevels_csfg

${GRASP}/bin/rsave_csfg: rsave_csfg
	cp $^ $@
	chmod u+x $@

${GRASP}/bin/rcsfgexpand_csfg: rcsfgexpand_csfg.o
	$(FC) -o $@ $?

${GRASP}/bin/rcsfgextend_csfg: rcsfgextend_csfg.o
	$(FC) -o $@ $?

${GRASP}/bin/rcsfgsplit_csfg: rcsfgsplit_csfg.o
	$(FC) -o $@ $?

${GRASP}/bin/rdistHmatrix_csfg: rdistHmatrix_csfg.o
	$(FC) -o $@ $?

${GRASP}/bin/rcsfgblocksplit_csfg: rcsfgblocksplit_csfg.o
	$(FC) -o $@ $?

${GRASP}/bin/rlevels_csfg: rlevels_csfg.o
	$(FC) -o $@ $?

%.o: %.f90
	$(FC) -c $(FC_FLAGS) -o $@ $<

clean:
	-rm -f  ${GRASP}/bin/rcsfgexpand_csfg ${GRASP}/bin/rcsfgsplit_csfg ${GRASP}/bin/rdistHmatrix_csfg ${GRASP}/bin/rcsfgextend_csfg ${GRASP}/bin/rsave_csfg ${GRASP}/bin/rcsfgblocksplit_csfg ${GRASP}/bin/rlevels_csfg
	-rm -f *.o *.mod
