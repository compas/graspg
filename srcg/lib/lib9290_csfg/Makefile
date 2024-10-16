LIBA=${GRASP}/lib/lib9290_csfg.a
MODULES_INSTALL=${GRASP}/lib/9290_csfg
FC_MODULES= -I ${GRASP}/src/lib/mpi90 -I ${GRASP}/srcg/lib/libmod_csfg -I ${GRASP}/src/lib/libmod -I ${GRASP}/src/lib/lib9290

OBJS= \
        indexsym.o  \
        indexsym_I.o   \
        analabv.o   \
        analabv_I.o  \
        csfg_ispar.o  \
        csfg_ispar_I.o  \
        labpack.o  \
        labpack_I.o \
        csfg_itjpo.o  \
        csfg_itjpo_I.o  \
        csfg_jcup.o  \
        csfg_jcup_I.o  \
        csfg_setqna.o  \
        csfg_setqna_I.o  \
        fictious_csf.o  \
        findtype.o  \
        printfictcsf.o \
        printlabelv.o  \
        setlaborb.o  \
        setlaborb_I.o \
        testfiles.o  \
        testfiles_I.o \
        lodstate_csfg.o  \
        lodstate_csfg_I.o \
        lodcsh_csfg.o  \
        lodcsh_csfg_I.o \
        lodcsh2_csfg.o  \
        lodcsh2_csfg_I.o \
        lodcslmpi_csfg.o  \
        lodcslmpi_csfg_I.o \
        setcsll_csfg.o  \
        setcsll_csfg_I.o \
        setcsh_csfg.o  \
        setcsh_csfg_I.o \
        cslhmpi_csfg.o   \
        cslhmpi_csfg_I.o   \
        rcsfsymexpand.o  \
        rcsfsymexpand_I.o  \

PHONY: install
install: $(LIBA)
	mkdir -p $(MODULES_INSTALL)
	cp -v *.mod $(MODULES_INSTALL)

$(LIBA): $(OBJS)
	@echo "Installing $@"
	ar -curs $@ $?

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(LIBA)
	-@rm *.o *.mod
	-@rm -R $(MODULES_INSTALL)
