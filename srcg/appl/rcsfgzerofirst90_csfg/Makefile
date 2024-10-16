EXE=${GRASP}/bin/rcsfgzerofirst_csfg
LIBS=-L ${GRASP}/lib/ -l9290 -lmod -l9290_csfg
FC_MODULES= -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod -I ${GRASP}/srcg/lib/lib9290_csfg

OBJS= \
	symzf_mod.o    \
	removeblank_I.o \
	lodcsl_Part_I.o \
	lodcsl_Zero_I.o \
	lodcsl_Full_I.o \
	set_CSF_ZFlist_I.o \
	set_CSF_number_I.o \
	lodcsl_Part.o \
	lodcsl_Zero.o \
	lodcsl_Full.o \
	set_CSF_ZFlist.o \
	set_CSF_number.o \
	removeblank.o \
	set_CSF_numberF.o \
	RCSFzerofirst.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
