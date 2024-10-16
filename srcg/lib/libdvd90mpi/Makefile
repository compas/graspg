LIBA=${GRASP}/lib/libdvdmpi.a
MODULES_INSTALL=${GRASP}/lib/dvd90mpi
FC_MODULES= -I ${GRASP}/src/lib/libmod -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/mpi90

OBJS= \
	adds_I.o \
	dvdrvr_I.o \
	dvdson_I.o \
	initdvd_I.o \
	mgs_nrm_I.o \
	multbc_I.o \
	newvec_I.o \
	ovflow_I.o \
	tstsel_I.o \
	dvdson.o \
	gdvd.o \
	gdvd_I.o \
	iniest.o \
	iniest_I.o

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
