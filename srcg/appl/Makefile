SUBDIR = jj2lsj90_csfg  \
         rci90_csfg_mpi \
         rci90_block_csfg_mpi \
         rcsfggenerate90_csfg \
         rmixaccumulate_csfg \
         rangular90_csfg_mpi \
         rmcdhf90_csfg_mpi \
	 rwfnestimate90_csfg \
         rcsfginteract90_csfg \
         rcsfgzerofirst90_csfg

TARGETS	 = install

$(TARGETS):
	@for i in $(SUBDIR) ;\
	do \
		echo "Current directory: $$i" ; \
	                cd $$i ; \
			make; \
echo "....Leaving: $$i" ; echo; echo;\
		cd .. ; \
	done

clean :
	@for i in $(SUBDIR) ;\
	do \
	echo "Current directory: $$i" ; \
           cd $$i ; \
           make clean; \
           cd .. ; \
        done

