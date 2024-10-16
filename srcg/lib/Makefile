

SUBDIR = libdvd90mpi libmod_csfg lib9290_csfg librang90_csfg libme_csfg 
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
clean:
	@for i in $(SUBDIR) ;\
        do \
                echo "Current directory: $$i" ; \
                        cd $$i ; \
                        $(MAKE) $@ ; \
                cd .. ; \
        done


