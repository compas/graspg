EXE=${GRASP}/bin/rcsfginteract_csfg
#LIBS=-L ${GRASP}/lib/ -lrang90 -l9290 -lmod
# LIBS=-L ${GRASP}/lib/ -lrang90 -l9290
#LIBS=-L ${GRASP}/lib/ -l9290
#FC_MODULES= -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod

FC_FLAGS=-fcray-pointer #-fallow-argument-mismatch

OBJS= \
        vast_kind_param_M.o  parameter_def_M.o \
        memory_man.o \
        blk_C.o cons_C.o debug_C.o default_C.o def_C.o dumx_C.o facts_C.o hblock_C.o iccu_C.o iounit_C.o m_C.o mcp_C.o mtjj2_C.o mtjj_C.o orb_C.o rang_Int_C.o stat_C.o terms_C.o trk_C.o \
        ittk.o  \
        ittk_I.o \
        ixjtik.o     \
        ixjtik_I.o     \
        snrc.o  \
        snrc_I.o \
        dracah.o \
        dracah_I.o \
        sixj1.o  \
        sixj1_I.o \
        sixj2.o  \
        sixj2_I.o \
        sixj3.o  \
        sixj3_I.o \
        sixj4.o  \
        sixj4_I.o \
        sixj35.o  \
        sixj35_I.o \
        sixj5.o  \
        sixj5_I.o \
        sixj.o  \
        sixj_I.o \
        nine0.o \
        nine0_I.o \
        nine.o \
        nine_I.o \
        diaga1.o     \
        diaga1_I.o   \
        diaga2.o \
        diaga2_I.o \
        diaga3.o \
        diaga3_I.o \
        diaga4.o \
        diaga4_I.o \
        convrt2.o     \
        convrt2_I.o     \
        convrt_double.o     \
        convrt_double_I.o     \
        convrt.o     \
        convrt_I.o     \
        clrx.o  \
        clrx_I.o  \
        cre.o   \
        cre_I.o \
        coulom.o     \
        coulom_I.o     \
        eile.o     \
        eile_I.o     \
        factt.o     \
        factt_I.o     \
        iq.o     \
        iq_I.o     \
        ichkq1.o     \
        ichkq1_I.o     \
        ichkq2.o     \
        ichkq2_I.o     \
        ichop.o     \
        ichop_I.o     \
        ispar.o     \
        ispar_I.o     \
        itjpo.o     \
        itjpo_I.o     \
        itrexg.o     \
        itrexg_I.o     \
        itrig.o     \
        itrig_I.o     \
        jcup.o     \
        jcup_I.o     \
        jqs.o     \
        jqs_I.o     \
        nmtejj.o    \
        nmtejj_I.o  \
        openfl.o     \
        openfl_I.o     \
        pack.o     \
        pack_I.o     \
        parsjl.o     \
        parsjl_I.o     \
        perko1.o     \
        perko1_I.o     \
        perko2.o     \
        perko2_I.o     \
        prsrcn.o     \
        prsrcn_I.o     \
        prsrsl.o     \
        prsrsl_I.o     \
        reco3.o     \
        reco3_I.o     \
        rec3.o     \
        rec3_I.o     \
        reco2.o     \
        reco2_I.o     \
        reco4.o     \
        reco4_I.o     \
        setqna.o \
        setqna_I.o \
        reco.o \
        reco_I.o \
        recoonescalar.o \
        recoonescalar_I.o \
	onescalar1INT.o \
	onescalar1INT_I.o \
	onescalar2INT.o \
	onescalar2INT_I.o \
	onescalarINT.o \
	onescalarINT_I.o \
	el1INT.o \
	el1INT_I.o \
	el2INT.o \
	el2INT_I.o \
	el31INT.o \
	el31INT_I.o \
	el32INT.o \
	el32INT_I.o \
	el33INT.o \
	el33INT_I.o \
	el3INT.o \
	el3INT_I.o \
	el41INT.o \
	el41INT_I.o \
	el4INT.o \
	el4INT_I.o \
	el51INT.o \
	el51INT_I.o \
	el52INT.o \
	el52INT_I.o \
	el53INT.o \
	el53INT_I.o \
	el5INT.o \
	el5INT_I.o \
        starttime.o \
        starttime_I.o \
        stoptime.o \
        stoptime_I.o \
	Interact_csf.o \
	Interact_csf_I.o \
        getyn.o \
        getyn_I.o \
	getinf.o \
	getinf_I.o \
	lodcsl_CSF.o \
	lodcsl_CSF_I.o \
	lodcsl_MR.o \
	lodcsl_MR_I.o \
	set_CSF_list.o \
	set_CSF_list_I.o \
	set_CSF_number.o \
	set_CSF_number_I.o \
	set_CSF_number_full.o \
	set_CSF_number_full_I.o \
	testfiles.o \
	testfiles_I.o \
	RCSFinteract.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? 

%.o: %.f90
	$(FC) -c $(FC_FLAGS) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
