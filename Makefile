#######################
#    OCAML Programs   #
#######################

OCAMLC = ocamlc
OCAMLOPT = ocamlopt
OCAMLDEP = ocamldep



#######################
#     OCAML Flags     #
#######################

INCLUDES =					#all relevant -I options here
OCAMLFLAGS = $(INCLUDES)	#add other options for ocamlc here
OCAMLOPTFLAGS = $(INCLUDES)	#add other options for ocamlopt here



#######################
#    Sources files    #
#######################

CMAFILES = bigarray.cma unix.cma
MLFILES = alphaadn.ml alphaprot.ml common.ml config.ml seq.ml traiteSignal.ml preserveorf.ml collection.ml donnees_base.ml normalisation_mrna.ml donnees.ml classifgene.ml db.ml traitementHSP.ml pretraitement1.ml pretraitement2.ml orfsearch.ml matrice.ml graph.ml graph_aux.ml traiteCollection.ml sigsearch.ml printing.ml premain.ml main.ml
MLIFILES = cast.mli
CMXAFILES = $(CMAFILES:%.cma=%.cmxa) 
CMOFILES = $(MLFILES:%.ml=%.cmo) 
CMXFILES = $(MLFILES:%.ml=%.cmx)
CMIFILES = $(MLFILES:%.ml=%.cmi) $(MLIFILES:%.mli=%.cmi)
OBJFILES = $(CMIFILES:%.cmi=%.o)
BINFILE = exogean
DEPENDFILE = .depend



#######################
#        Rules        #
#######################

exogean: $(CMXFILES)
	@echo "LNKOPT $(BINFILE)"
	@$(OCAMLOPT) $(CMXAFILES) $(CMXFILES) $(OCAMLOPTFLAGS) -o $(BINFILE)

exogeantmp: $(CMOFILES)
	@echo "LNK $(BINFILE)"
	@$(OCAMLC) $(CMAFILES) $(CMOFILES) $(OCAMLFLAGS) -o $(BINFILE)


# Common rules
.SUFFIXES: .ml .mli .cmo .cmi .cmx

.ml.cmo:
	@echo "OCAMLC $<"
	@$(OCAMLC) $(OCAMLFLAGS) -c $<

.mli.cmi:
	@echo "OCAMLC $<"
	@$(OCAMLC) $(OCAMLFLAGS) -c $<

.ml.cmx:
	@echo "OCAMLOPT $<"
	@$(OCAMLOPT) $(OCAMLOPTFLAGS) -c $<

# Clean up
clean:
	@echo "Cleaning .cmo .cmx .o .cmi and binary"
	@rm -f $(CMOFILES) $(CMXFILES) $(OBJFILES) $(CMIFILES)

# Dependencies
depend:
	@echo "Calculating dependencies"
	@$(OCAMLDEP) $(INCLUDES) $(MLFILES) $(MLIFILES) > $(DEPENDFILE)

include $(DEPENDFILE)


