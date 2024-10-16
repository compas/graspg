EXE=${GRASP}/bin/rmixaccumulate_csfg
#LIBS=
#FC_MODULES=
#LIBS=-L ${GRASP}/lib/ -lrang90 -lmcp90 -lmpiu90 -l9290 -lmod
#FC_MODULES= -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/libmcp90 -I ${GRASP}/src/lib/mpi90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod
LIBS=-L ${GRASP}/lib/ -l9290_csfg -lmod_csfg -L ${GRASP}/lib/ -lmod
FC_MODULES= -I ${GRASP}/src/lib/libmod -I ${GRASP}/srcg/lib/lib9290_csfg -I ${GRASP}/srcg/lib/libmod_csfg

OBJS= \
	rmixaccumulate_csfg.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? ${LIBS}

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
