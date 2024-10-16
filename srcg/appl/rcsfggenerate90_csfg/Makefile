EXE=${GRASP}/bin/rcsfggenerate_csfg
LIBS=
FC_MODULES=

OBJS= \
	starttime.o \
	starttime_I.o \
	stoptime.o \
	stoptime_I.o \
	lockad.o \
	lockad_I.o \
	reada.o \
	reada_I.o \
	lasa1.o \
	lasa1_I.o \
	lasa2.o \
	lasa2_I.o \
	adder.o \
	adder_I.o \
	slug.o \
	slug_I.o \
	sluggo.o \
	sluggo_I.o \
	kopp1.o \
	kopp1_I.o \
	kopp2.o \
	kopp2_I.o \
	genb.o \
	gen_I.o \
	copy7t9.o \
	copy7t9_I.o \
	fivefirst.o \
	fivefirst_I.o \
	fivelines.o \
	fivelines_I.o \
	lika.o \
	lika_I.o \
	matain.o \
	matain_I.o \
	matbin.o \
	matbin_I.o \
	matcin.o \
	matcin_I.o \
	test.o \
	test_I.o \
	merge.o \
	merge_I.o \
	mergeb.o \
	mergeb_I.o \
	open79.o \
	open79_I.o \
	reffa.o \
	reffa_I.o \
	blanda.o \
	blanda_I.o \
	blandb.o \
	blandb_I.o \
	blandc.o \
	blandc_I.o \
	jjgen15b.o \
	rcsfblock.o \
	rcsfexcitation.o \
	wrapper.o \
	rcsfsym_cyc.o \
	rcsfsym.o \
	removeblank.o

$(EXE): $(OBJS)
	$(FC) -o $@ $?

%.o: %.f90
	$(FC) -c $(FC_FLAGS) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
