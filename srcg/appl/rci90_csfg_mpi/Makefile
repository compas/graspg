EXE=${GRASP}/bin/rci_csfg_mpi
LIBS=-L ${GRASP}/lib -lmpiu90 -l9290 -lmod
LIBSG=-L ${GRASP}/lib -lmod_csfg -lme_csfg -lrang90_csfg -l9290_csfg -ldvdmpi
FC_MODULES= -I ${GRASP}/src/lib/mpi90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod -I ${GRASP}/srcg/lib/libmod_csfg -I ${GRASP}/srcg/lib/libme_csfg -I ${GRASP}/srcg/lib/librang90_csfg -I ${GRASP}/srcg/lib/lib9290_csfg -I ${GRASP}/srcg/lib/libdvd90mpi

#	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBSG) ${LIBS} $(LAPACK_LIBS)
OBJS= \
      setham_gg.o  \
      setham_gg_I.o  \
      genmat.o  \
      genmat_I.o  \
      getcid.o  \
      getcid_I.o  \
      lodres.o  \
      lodres_I.o  \
      maneigmpi.o  \
      maneigmpi_I.o  \
      matrix.o  \
      matrix_I.o  \
      setres.o  \
      setres_I.o  \
      rci90mpi.o  

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(LIBSG) ${LIBS} $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
