LIBA=${GRASP}/lib/libmod_csfg.a
MODULES_INSTALL=${GRASP}/lib/mod_csfg
FC_MODULES= -I ${GRASP}/src/lib/libmod -I ${GRASP}/src/lib/lib9290

OBJS= \
	csfg_decide_C.o symexpand_mod.o symmatrix_restart_C.o symmatrix_mod.o csfg_memory_man.o csfg_tv_C.o csfg_scf_C.o

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
