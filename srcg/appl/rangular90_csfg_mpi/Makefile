EXE=${GRASP}/bin/rangular_csfg_mpi
LIBS=-L ${GRASP}/lib/ -lrang90_csfg -lmod_csfg -l9290_csfg -lmcp90 -lmpiu90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/srcg/lib/librang90_csfg -I ${GRASP}/src/lib/libmcp90 -I ${GRASP}/src/lib/mpi90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod -I ${GRASP}/srcg/lib/libmod_csfg -I ${GRASP}/srcg/lib/lib9290_csfg

OBJS= \
     sortoffdv.o \
     sortdiagv.o \
     fndbeg.o \
     fndbeg_I.o \
     getinf.o \
     getinf_I.o \
     transfer_csfg.o \
     outsdampi.o \
     outsdampi_I.o \
     setdbg.o \
     setdbg_I.o \
     setdbgmpi.o \
     setdbgmpi_I.o \
     setmcp_csfg.o \
     setmcp_csfg_I.o \
     setmcpmpi.o \
     setmcpmpi_I.o \
     setsda_csfg.o \
     setsda_csfg_I.o \
     setsum.o \
     setsum_I.o \
     settmp.o \
     settmp_I.o \
     sort.o \
     sort_I.o \
     strsum.o \
     strsum_I.o \
     savetv.o   \
     spinangular1.o \
     spinangular11.o \
     spinangular12.o \
     spinangular13.o \
     spinangular14.o \
     spinangular15.o \
     spinangular22.o \
     spinangular23.o \
     spinangular24.o \
     spinangular25.o \
     spinangular2.o \
     spinangular33.o \
     spinangular34.o \
     spinangular35.o \
     spinangular3.o \
     spinangular44.o \
     spinangular45.o \
     spinangular4.o \
     spinangular55.o \
     spinangular5.o \
     mcpmpi_csfg.o \
     mcpmpi_csfg_I.o \
     genmcpmpi.o  

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
