(**********************************************************************************************************************************************************************)
(* Ce fichier presente toutes les fonctions necessaires pour recuperer le contexte a partir de ce que l'utilisateur a rentre (fichier de configuration, arguments ... *)
(* This file includes all the functions to retrieve the context from the user input (configuration file, arguments ...)                                               *)
(**********************************************************************************************************************************************************************)



(***********************************************************************************************************************************************************************)
(* Type indiquant le mode de sélection des signaux d'epissage : en fonction du nombre d'acides aminés manquant -> Aamiss, ou en fonction du score par matrice de poids *)
(* Type for the way of selection of splicing signals in case several are possible: either depending on the number of amino acids missing in the protein -> Aamiss,     *)
(* or on the weight matrix score of the signals                                                                                                                        *)
(***********************************************************************************************************************************************************************)
type wsigsearch = Aamiss | Wghtmat;; (* modif au 27/02/06 pas d'entier apres Aamiss car deja dans naa 
				        modified on 2006/02/27: no integer parameter for Aamiss since already in naa 
				     *)



(*********************************************************)
(* La structure qui contient les parametres du programme *)
(* Structure for Exogean parameters                      *)
(*********************************************************)
type 'a context_t =
	{ 
	  mutable gseqfile:string;			(* nom du fichier de la séquence génomique à annoter, au format fasta *)
                                                        (* name of the fasta file including the genomic sequence to annotate  *)
	  mutable bqpfile:string;			(* nom du fichier des séquences protéiques *) 
                                                        (* name of the fasta file including the sequences of the proteins against which the genomic sequence has been compared *)
	  mutable bqafile:string;			(* nom du fichier des séquences arnm matchés par gseq *) 
                                                        (* name of the fasta file including the sequences of the mrna against which the genomic sequence has been compared *)
	  mutable hsppfile:string;			(* nom du fichier des alignements protéiques *) 
	                                                (* name of the file including the protein alignments *)
	  mutable hspafile:string;			(* nom du fichier des alignements arnm *) 
	                                                (* name of the file including the mrna alignments *)
	  mutable outfile:string;			(* nom de base des fichiers de sortie d'Exogean, on rajoute l'extension .gene, .prot .... *)
	                                                (* base name for Exogean output files *)
	  mutable magfile:string;			(* nom du fichier contenant les sequences fasta centrees en ag *)
                                                        (* name of the optional fasta file including the 10nt sequences for acceptor signals *)
	  mutable mgtfile:string;			(* nom du fichier contenant les sequences fasta centrees en gt *)
	                                                (* name of the optional fasta file including the 10nt sequences for donor signals *)
	  mutable fusionthresh:int;			(* TODO: mettre a jour. Seuil en nb nts pour la séparation entre deux hsps à fusionner *)
	                                                (* TODO: update. Theshold in number of nucleotides for the separation between two consecutive HSPs to merge *)
	  mutable waysigsearch:wsigsearch;	        (* façon dont on veut rechercher les signaux *)
	                                                (* way to look for splicing signal, either by the number of amino acids that are missing in the protein, or by the splicing signal weight matrix *)
	  mutable alignement:int;			(* 1 si on veut un alignement supplémentaire % tblastn , 0 sinon *)
	                                                (* 1 if we want to perform an alignment more than the tblastn, 0 otherwise *)
	  mutable matrice:string;			(* nom du fichier contenant la matrice protéique pour l'alignement *)
                                                        (* name of the file including the alignment proteic matrix *)
	  mutable bplusa:int;				(* seuil inférieur pour la différence entre ce qu'il manque dans le nucléique et ce qu'il manque dans la protéine entre deux hsps pour lancer un nouvel alignement *)
	                                                (* *)
	  mutable diffaa:int;				(* seuil inférieur pour ce qu'il manque dans la protéine entre deux hsps pour lancer un nouvel alignement *)
	  mutable scut:int;				(* seuil inférieur pour le score des hsps à retenir lors du nouvel alignement *)
	  mutable beta:int;				(* paramètre beta, mesure le laxisme *)
	  mutable tmaxintron : int;			(* taille maximum d'un intron, perspective : blusp monoprotéique trop étendu en génomique*)
	  mutable tmaxintron2 : int;			(* taille maximum d'un intron no 2 pour le cas de l'hsp isolé = vieux pseudogène, utilisé dans blusping monoprot *)
	  mutable tminintron: int;			(* taille maximum d'un intron *)
	  mutable nbaamanq: int;			(* nombre d'acides amines manquants par defaut, utilisé pour la recherche de signaux  *)
	  mutable nninhsp : int;		        (* nombre de nucléotides duquel on rentre dans l'hsp protéique pour rechercher le signal *)
	  mutable nnsuppl : int;		        (* nombre d'acides aminés duquel on additionne le nombre d'aa manquants dans la protéine lors de la recherche du signal *)
	  mutable szisolh : int;                        (* hsp isolé dans bluspmonoprot, représentant un vieux pseudogène, sans retour en arriere avec le reste du blusp, c'est la taille min d'un hsp pour etre toléré *)
	  mutable szhsp : int;                          (* taille min d'hsp prot conservé *)
	  mutable szcdsnomet : int;			(* taille minimale en nt d'un CDS sans methionine *)
	  mutable szcdswithmet : int;			(* taille minimale en nt d'un CDS avec methionine*)
	  mutable propcdsmonoex  : int;		        (* pourcentage (entre 0 et 100) d'un CDS monoex par rapport à la taille cumulée des exons du transcrit *)
	  mutable pnntfus : int;
	  mutable pnaamq : int;
	  mutable propmol : int;			(* pourcentage (entre 0 et 100) du nb aa de l'hsp pale d'un blusp mono ou dicluspique par rapport à la taille de la proteine principale *)
	  mutable propmol2 : int;                       (* idem mais pour raffiner regle elimination monoprotart *)
	  mutable sizecdsmonoex : int;		        (* taille minimale d'un cds monoexonique *)
	  mutable sizecdsdiex : int;			(* taille minimale d'un cds diexonique *) 
	  mutable protmonohsp : bool;			(* est-ce qu'on garde les protéines mono hsps *)
	  mutable pipe : string;			(* Indique quel fichier de sortie doit etre redirige sur la sortie standard *)
	  mutable offset: int;
	  mutable ajuste : bool;                        (* pour rogner les hsps *)
    } 



(********************************************************)
(* Designe le contexte de l'application, des parametres *)
(********************************************************)
let context = 
	{	
		gseqfile = "";
		(*"/export/home1/ldog/djebali/SARAH/GENCODE/Sequence/"^reg^"_Build35.fa";*)
		(*"/export/home1/ldog/djebali/SARAH/EXOGEAN/OCAML/Khaled_srcprot/Build34_chr22.fa";*)
		(* Attention : ici on peut avoir des gdes et petites lettres, EXOGEAN les traite de la meme façon *)
 		bqpfile = "";
		(*"/export/home1/ldog/djebali/SARAH/IPIMm/Mar2005/ipi.MOUSE_nom_court.fasta";*)
		(*"/export/home1/ldog/djebali/SARAH/IPIMm/Jul2004/ipi.MOUSE_nom_court.fasta"; *)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/Tetraodon_Zebrafish/pep_Tn_ipi_Zebra_sans_stop.fa";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/Poulet_Tetraodon/ensfeb_GG_pep_Tn_sans_stop.fa";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/Tetraodon_Zebrafish/pep_Tn_ipi_Zebra_sans_stop.fa";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/Souris_Poulet/ipi_Mouse_ensfeb_GG_sans_stop.fa";*) 
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/All5species/all_ipiMM_ipiRat_ensGG_ipiZebra_pepTn_sans_stop.fa";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/RAT/ipi.RAT.nom_court.fa";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/POULET/Gallus_gallus.WASHUC1.feb.pep.nom_court.fa";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/ZEBRAFISH/ipi.ZEBRA.nom_court.fasta";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/TETRAODON/tetraodon_peptides_v1.nom_court.fa";*)
		(*"/export/home1/ldog/djebali/SARAH/IPIMm/Jul2004/ipi.MOUSE_nom_court.fasta"; *)
 		bqafile = "";
		hsppfile = "";
		(* TODO: ne pas oublier de rajouter le prefixe 'list:', 'psl:' ... *)
		(*"/export/home1/ldog/djebali/SARAH/GENCODE/Proteines/IPImouse_Mar05_vs_"^reg^"_ENCODE_Hg17_C_min25.list";*)
		(*"/export/home1/ldog/djebali/SARAH/EXOGEAN/OCAML/srccomb/IPImouse_vs_K22_C_min25.list"; *)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/Incremental_Mm/IPImouse_sub90.list";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/Tetraodon_Zebrafish/Tn_Zebra_hsps.list";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/Poulet_Tetraodon/GG_TN_hsps.list";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/Tetraodon_Zebrafish/Tn_Zebra_hsps.list";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/Souris_Poulet/MM_GG_hsps.list";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/All5species/all_MM_Rat_GG_Zebra_Tn_vs_K22_Bld34.list"; *)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/RAT/RATipi_vs_K22NCBI34.list";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/POULET/POULETens_vs_K22NCBI34.list";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/ZEBRAFISH/ZEBRAipi_vs_K22NCBI34.list";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Proteines/TETRAODON/TETRA_vs_K22NCBI34.list";*)
		(*"/export/home1/ldog/djebali/SARAH/EXOGEAN/OCAML/srccomb/IPImouse_vs_K22_C_min25.list";*)
		hspafile = "";
		(* TODO: ne pas oublier de rajouter le prefixe 'list:', 'psl:' ... *)
		(*"/export/home1/ldog/djebali/SARAH/GENCODE/Transcrits/all_mrna_mar05_vs_ENCODE_"^reg^"_Hg17_4_best.list";*)
		(*"/export/home1/ldog/djebali/SARAH/EXOGEAN/OCAML/srccomb/chr22_mrna_hrc.exog";*)
		(*"/export/home1/ldog/djebali/SARAH/K22_NCBI_34/Transcrits/Incremental/K22mRNA_sub90.list";*)
		(*"/export/home1/ldog/djebali/SARAH/EXOGEAN/OCAML/srccomb/chr22_mrna_hrc.exog";*)
		outfile = "";
		magfile = "/export/home1/ldog/djebali/SARAH/EXOGEAN/OCAML/Khaled_srcprot/sig_AG.fa";
		mgtfile = "/export/home1/ldog/djebali/SARAH/EXOGEAN/OCAML/Khaled_srcprot/sig_GT.fa";
		fusionthresh = 30;
		waysigsearch = Wghtmat;
		alignement = 0;
		matrice = "BLOSUM62";
		bplusa = 50;
		diffaa = 25;
		scut = 50;
		beta = 1;
		tmaxintron = 75000;
		tmaxintron2 = 10000;
		tminintron = 60;
		nbaamanq = 10;
		nninhsp = 8;
		nnsuppl = 15;
		szisolh = 45;
		szhsp = 25;
		szcdsnomet = 300;		  (* au départ 150 nts, paramètre de la mise en pseudogènes de transcrits en aval *)
		szcdswithmet = 300;		(* au départ 240 nts, paramètre de la mise en pseudogènesde transcrits en aval *)
		propcdsmonoex	= 20;    (* si on met 50 semble pas changer gd chose, paramètre de la mise en pseudogènesde transcrits en aval *)
		pnntfus = 10;
		pnaamq = 5;
		propmol = 33; 
                propmol2 = 20;       (* paramètres de l'elimination de blusps protéiques en amont *)
		sizecdsmonoex = 420;
		sizecdsdiex = 300;		 (* en fait pour l'instant on ne s'en sert pas car trop radical *)
		protmonohsp = false;
		pipe = "";
		offset = 0;
		ajuste=true;
	};;


(*************************************************************************************)
(* Message affiche lorque l'utilisateur s'est trompe dans les arguments du programme *)
(*************************************************************************************)

let usage_us =
"

      EXOGEAN: EXpert On GEne ANnotation ; version v2.0 (January 2010)
           
         (Note: this version of Exogean is a slightly improved version 
            of the one used for the EGASP project and published in
                    Djebali et al. Genome Biology 2006)


Usage : "^(Filename.basename (Sys.argv.(0)))^" -gseq DNA_file.fa  [options]


As an alternative to options being specified on the command line,  all the arguments 
given to the program can be recapitulated in a file called exogean.ini that must be 
in the directory where exogean is executed. In this case, simply execute "^(Filename.basename (Sys.argv.(0)))^" 
with no arguments. If arguments are nevertheless provided on the command line, they 
override the corresponding ones in the exogean.ini file. 

Exogean annotates protein-coding genes in a genomic DNA sequence based on information 
provided by alignments of protein and/or mRNA sequences to this DNA. The only mandatory 
argument is the DNA sequence. Obviously if nothing else is provided, Exogean will not 
annotate any genes and will not produce any files. Alignments are provided to exogean 
in separate files for mRNAs (-hspa) and proteins (-hspp). If a protein alignment file
has been provided, Exogean also requires the file of protein sequences that were used 
to compute the alignments, with the option -bqp.

Exogean uses rules extracted from human expertise to process the alignments and combine
this information to annotate protein coding genes in DNA. Exogean is able to identify 
several transcripts per genes if mRNA alignments are provided, but not if only protein
alignments were given. Depending on the evidence given, Exogean may predict no genes, 
genes and/or pseudogenes, genes with overlapping boundaries, genes nested in introns of 
other genes on the same strand, etc. 

Briefly, Exogean will process the information in 4 stages: pre-processing, DACM algorithm, 
CDS identification and post-processing. Options below allow the user to modulate the pre- and 
post-processing stages, while rules from the DACM algorithm and for the CDS identification
are currently hard coded. A new language is being developed to let the user modify all the rules 
and create new ones. 

The output of Exogean consists in six files:

1. The gene annotation file, in either gtf or bed format. The extension is 
   *.gene.gtf or *.gene.bed. For a description of these formats please see 
   http://genes.cs.wustl.edu/GTF2.html or 
   http://genome.ucsc.edu/goldenPath/help/customTrack.html#BED 
2. The cDNA sequence of each annotated transcript, in fasta format. The extension 
   is *.cdna.fa. 
3. The translated CDS sequences of each annotated transcript, in fasta format. 
   The extension is *.cds.fa
4. The gene annotations in html format. The extension is *.gene.html. This file 
   can be uploaded in a browser to allow for an easier navigation in the gene 
   annotations and their respective transcripts, and contains quality information.
5. Putative pseudogenes and/or artifacts in the same format as the gene annotation 
   file. The extension is either *.pseudo.bed or *.pseudo.gtf.
6. A list of mRNA alignments that passed the initial filters, before entering the 
   DACM algorithm. The file is empty if no mRNA alignments were provided. This 
   file can generally be ignored, but can be useful to understand why some genes 
   where poorly annotated even if mRNA sequences for these genes were provided.


List of file arguments: 
----------------------

-gseq string          string=DNA sequence fasta file (mandatory). Currently a single 
                      DNA sequence per file is supported.  

-hspp format:file     The file containing the protein alignments (optional). The 
                      alignments can be provided in one of three formats that 
                      can be specified before the semicolons: psl, gff or exf. If 
                      only the file name is indicated, the format is assumed to be PSL. 
                      The gff format required here is very specific. It must contain 9 
                      fields separated by tabs: (1) DNA sequence name (not used), 
                      (2) name of the program that has generated the alignment (not used), 
                      (3) feature name (not used), (4) start of the alignment on the 
                      genomic sequence (1-based), (5) end of the alignment on the genomic 
                      sequence (1-based), (6) score of the alignment, (7) strand of the 
                      genomic sequence where the alignment lies, (8) frame of the genomic 
                      sequence to which the alignment corresponds (not used), (9) field 
                      containing at least 4 subfields separated by spaces: (1) type of 
                      protein (not used), (2) name of protein, (3) start of the alignment 
                      on the protein (1-based), (4) end of the alignment on the protein 
                      (1-based). Be careful: the gff file provided must only contain 
                      alignment features and should not contain any comment line.
                      The exf format is a simple format designed for Exogean and 
                      consists in 8 fields separated by spaces, tabs or a combination 
                      of spaces and tabs: (1) DNA sequence name, (2) protein/mRNA name, 
                      (3) a score for the alignment, (4) start in genomic DNA, (5) end 
                      in genomic DNA, (6) start in molecule, (7) end in molecule, 
                      (8) a plus or minus sign for the strand. When the strand is 
                      negative, then the start and end of the alignment are given 
                      relative to the plus strand. The PSL format is typically produced 
                      by BLAT alignments and is described here: 
                      http://genome.ucsc.edu/goldenPath/help/customTrack.html#PSL. 
                      A description of the Generic Feature Format (GFF) can be found 
                      here: http://www.sanger.ac.uk/Software/formats/GFF/. If a protein 
                      alignment file is provided, then a protein fasta file must also 
                      be provided with the -bqp argument (see below).

-bqp string	      string = protein sequence fasta file (mandatory if the -hspp option 
                      is used). Typically this is the same file that was used to compute 
                      the alignments with the genomic DNA, although here the protein 
                      sequences are preferably unmasked. The protein names in the 
                      alignment file and in the fasta file must be identical. The fasta 
                      file must contain all the proteins listed in the alignment file, 
                      although additional proteins that did not align may also be 
                      included and will be ignored. Protein names *must not* contain
                      pipe signs like this  -> | <- .   

-hspa format:string   The file containing mRNA alignments (optional). The filename can 
                      be preceded by a semicolon and the format: psl, gff or exf. If 
                      only the file name is indicated, the format is assumed to be PSL. 
                      See the -hspp option for information on file formats. 

-o format:string      The base for naming output files (optional). The name can be 
                      preceded by semicolon and the desired format: gtf or bed. 
                      See the -hspp option for information on file formats. If this 
                      option is not used, output files will be named after the name of 
                      the DNA sequence and the date and time of day, 
                      (sequencename_month_day_hour_minute) and the default format will 
                      be GTF. 
                      See http://genome.ucsc.edu/goldenPath/help/customTrack.html#GTF
                      for a description of the GTF format. 

-offset integer       nucleotide offset to add to the genomic coordinates of Exogean 
                      annotations (optional; default = 0). If annotations are computed 
                      on a subsequence of a chromosome, this option can be used to 
                      replace the positions on chromosome coordinate.  

-pipe string          string = prot or cdna or gene or pseudo (optional; default= no 
                      redirection). Enables one of the result files to be redirected 
                      to standard output. The other five files will still be written. 

-v            	      verbose mode (optional)



List of arguments to customize the behaviour of Exogean for pre-processing alignments.
-------------------------------------------------------------------------------------
Most are needed when building transcripts from structures based on protein 
alignments. In the following, the term HSP (inherited from BLAST) is used to describe 
a local alignment of a protein sequence on a genomic sequence. A single protein sequence 
may produce several HSPs, (e.g. corresponding to different exons). Options are listed 
below in the order in which they are used by exogean.

-sh integer           integer = threshold in nucleotide (default=25).
                      Proteic HSPs less than this threshold are totally eliminated from
                      Exogean processing.

-pmono                If used, allows Exogean to use proteins that produce a single HSP 
                      on the genomic DNA sequence. Otherwise, Exogean only uses proteins 
                      that produce more than one HSP on the genomic sequence (optional).

-f integer            integer = threshold in nucleotides (default=30). If two neighbouring 
                      HSPs are separated by a distance Dg in the genomic DNA and a 
                      distance Dp in the protein sequence, the two alignments are fused 
                      if Dg=(Dp*3) +/- threshold. 

-I integer   	      integer = threshold in nucleotides (default=75000). If two protein 
                      alignments are contiguous in both the genomic DNA and the aligned 
                      molecule but separated by more than the threshold on the genomic 
                      DNA, they will be considered as being part of different transcripts.

-I2 integer           integer = threshold in nucleotide (default=10000).
-sih integer          integer = threshold in nucleotide (default=45). 
                      If two protein alignments are contiguous in both the genomic DNA 
                      and the aligned molecule, and separated by a distance less than 
                      the threshold given by the I parameter but more than the threshold 
                      given by the I2 paramater on the genomic DNA, and one of them is less 
                      than the threshold given by the sih paramter then they will be considered 
                      as being part of different transcripts.

-pm integer           integer = percentage between 1 and 100 (default=33). 
-pm2 integer          integer = percentage between 1 and 100 (default=20). 
-pnamq integer        integer = percentage between 1 and 100 (default=5). 
-pnfus integer        integer = percentage between 1 and 100 (default=10).
                      Once the protein alignments have been processed and just before entering 
                      the DACM algorithm, a last verification is made on the validity of the HSPs. 
                      If a preliminary gene structure satisfies either one of the two following 
                      conditions, it will be considered as an artifact and eliminated :
                      - it is based on one or two HSPs and the ratio between the sum of their size 
                        and the size of the original protein is less than the threshold given by 
                        the pm option;
                      - it is based on more than three HSPs and either one of the three following
                        conditions is satisfied :
                           * the ratio between the number of missing amino acids and the number of
                             matching amino acids is more than the threshold given by the pnamq option,
                           * the ratio between the sum of the HSP sizes and the size of the original protein 
                             is less than the threshold given by the pm2 option,
                           * the ratio between the number of nucleotides merged between all the HSPs and the 
                             total number of nucleotides formed by the HSPs is more than the threshold given by 
                             the pnfus option.
                           

-naa integer	      integer = distance in codons (default = 10). This option is used only 
                      for exons that have protein alignment support, but no mRNA alignments.
                      HSPs are often shorter than the real exons and in such cases exogean 
                      must look outwards from the HSP for suitable splice donor and acceptor 
                      sites. Exogean will use the distance given to limit the region where 
                      putative splice sites will be searched for, unless the aligned protein 
                      has poorly conserved residues that did not align between the two HSPs, 
                      in which case the number of such amino-acids is used as the distance. 

-ni integer	      integer = distance in nucleotides (default=8). This option is used 
                      only for exons that have protein alignment support, but no mRNA 
                      alignments. Sometimes protein alignments extend into introns. To allow 
                      for such cases, Exogean will also search putative splice sites within 
                      this distance inside of HSPs. 

-ns integer           integer = distance in nucleotides (default=15; must be a multiple 
                      of 3). This option is used only for exons that have protein alignment 
                      support, but no mRNA alignments. 
                      When looking for splice signals outside or inside of HSPs, exogean will 
                      either be guided by the aligned protein if adjacent HSPs are from the 
                      same protein, or use the values given by the -naa option to limit the 
                      search. In either cases, the value given with the -ns option is added
                      as a tolerance threshold to extend further outside of HSPs, or to 
                      extend the limits provided by the aligned protein. 

-i integer            integer = threshold in nucleotides (default=60).  This option is used 
                      only for exons that have protein alignment support, but no mRNA 
                      alignments.
                      When searching for putative splice sites between two adjacent HSPs, 
                      exogean will eliminate all combinations of donor and acceptor sites 
                      that will create an intron of size less than this threshold. An
                      important consequence is that when only protein alignments are provided,
                      final genes structures will always have introns larger than this 
                      threshold

-mgt string    	      string = name of fasta file (optional). The multi fasta file must 
                      contain a list of 10 bp sequences that are used to compute a 
                      position weight matrix for splice donor sites. Position 1-4 are in 
                      the exon and 5-10 are in the intron. Exogean comes with a sig_GT.fa 
                      file that can be used here for human genes (and also works well for 
                      other mammalian species) that was extracted from the HS3D database 
                      (http://www.sci.unisannio.it/docenti/rampone/).

-mag string	      string = name of fasta file (optional). The multi fasta file must 
                      contain a list of 10 bp sequences that are used to compute a 
                      position weight matrix for splice acceptor sites. Position 1-6 are 
                      in the intron and 7-10 are in the exon. Exogean comes with a 
                      sig_AG.fa file that can be used here for human genes (and also works 
                      well for other mammalian species) that was extracted from the HS3D 
                      database (http://www.sci.unisannio.it/docenti/rampone/).

-wmat | -waa          makes a choice between two alternative methods to select the best 
                      splice donor and acceptor sites when information flanking a putative 
                      intron is only based on proteins (default: wmat). The -wmat option 
                      uses a position weight matrix score (built from the files given with 
                      the -mgt and -mag arguments) to select the most favourable intron 
                      between two HSPs. The -waa option instead uses the protein alignment 
                      as guide (if available), and will fill in a number of codons at the 
                      ends of putative exons, to approach the number of amino acids that 
                      are in the protein at this position but that could not be aligned. 


                    

List of arguments that control the last validation of Exogean predictions (post-processing)
-------------------------------------------------------------------------------------------

These filters are applied on Exogean predictions to classify them as \"correct\" or 
artifacts/pseudogenes. Ideally, Exogean predicts genes with a CDS starting with a 
methionine and ending at a stop codon. When these signals cannot be found or the CDS is 
small, Exogean may filter out some predictions based on the arguments given below. 

-scnmint              integer = size in nucleotides. Must be multiple of 3 
                      (Default=300). If a CDS is found but without a starting methionine, 
                      it must be at least this size. 

-scmint               integer = size in nucleotides. Must be multiple of 3 
                      (Default=210). If a CDS is found with a starting methionine, it 
                      must nevertheless be larger than this size to be considered valid. 

-pcdsint              integer = percentage between 1 and 100 (Default=20). The size of the 
                      CDS of a mono-exonic transcript must represent at least this fraction 
                      of the total size of the transcript (UTRs included) to be considered 
                      valid.

-scdsmint             integer = size in nucleotides (Default=420). If a CDS spans a 
                      unique exon, it must be at least this size to be considered valid, 
                      regardless of the presence of a starting methionine.

-scdsdint             integer = size in nucleotides (Default=210). If a CDS spans only 
                      two exons, it must be at least this size to be considered valid, 
                      regardless of a starting methionine.



------------------------------------------------------------------------------
To cite Exogean, please use

Djebali, S., Delaplace, F., Roest Crollius, H., Exogean: a framework for annotating 
protein-coding genes in eukaryotic genomic DNA. Genome Biology, 2006, 7 Suppl 1:S7.1-10. 

For bug reporting, questions or comments, please email exogean@biologie.ens.fr.

Exogean is developped in the Dyogen Group at the Ecole Normale Supérieure
(CNRS UMR8541) in Paris in collaboration with the IBISC (CNRS FRE2873) in Evry. 
-------------------------------------------------------------------------------

";;

let usage_fr=usage_us;;


(*
Pour élimination blusp monoprot artefactuels 
---------------------------------------------
- proportion d'acides aminés manquants par rapport au nombre d'acides aminés matchés dans la protéine en dessous de laquelle on tolère un blusp monoprot ayant au moins trois HSPs :    * pnaamq = ... dans exogean.ini
   * - pnamq integer en ligne de commande
   * 5 par défaut

- proportion de nucléotides fusionnés par rapport au nombre de nucléotides constitués par les HSPs :
   * pnntfus = ... dans exogean.ini
   * -pnfus integer en ligne de commande
   * 10 par défaut

Autre
-----
- Taille minimum des hsp protéiques utilisés :     *  szhsp= ... dans exogean.ini
  * - sh integer en ligne de commande
  * 25 par défaut

Pour découpe blusp monoprot en vue artefact/pseudogenes
-------------------------------------------------------
-  Taille maximum des hsps isolés dans blusp monoprot (ancien pseudogène)
  *  szisolh = ... dans exogean.ini
  * -sih integer en ligne de commande
  * 45 par défaut

- Taille maximale séparant deux hsps contigus dans bluspmonoprot dont l'un représente un hsp isolé (ancien pseudogène)
   * tmaxintron2 = ... dans exogean.ini
  * -I2 integer en ligne de commande
  * 10000 par défaut

*)

(*
[The following options are not currently used but proposed for a future usage]
-a bool		: -a 1 si on veut un alignement supplémentaire, -a 0 sinon (default: 0)
-m string	: nom du fichier contenant la matrice protéique utilisée pour l'alignement (default: BLOSUM62)
-ba int		: seuil inférieur pour la différence entre ce qu'il manque dans le nucléique et ce qu'il manque dans la protéine entre deux hsps pour lancer un nouvel alignement (default: 50)
-aa int		: seuil inférieur pour ce qu'il manque dans la protéine entre deux hsps pour lancer un nouvel alignement (default: 25)
-s int		: seuil inférieur pour le score des hsps à retenir lors du nouvel alignement (default: 50)
*)




(******************************************************************)
(* Lit le fichier de configuration et met a jour la configuration *)
(******************************************************************)
let read_configfile () =

	(* Les noms de parametres qu'on peut lire dans le fichier *)
	let params_fichier = [| "tmaxintron2"; 
				"propmol2";
				"gseqfile"; 
				"bqpfile"; 
				"bqafile"; 
				"hsppfile"; 
				"hspafile"; 
				"fusionthresh"; 
				"alignement"; 
				"matrice"; 
				"bplusa"; 
				"diffaa";
				"scut"; 
				"tmaxintron";
				"nninhsp"; 
				"nnsuppl"; 
				"szcdsnomet"; 
				"szcdswithmet"; 
				"propcdsmonoex"; 
				"propmol"; 
				"sizecdsmonoex";
				"outfile"; 
				"sigsearch"; 
				"tminintron"; 
				"nbaamanq"; 
				"sizecdsdiex"; 
				"magfile"; 
				"mgtfile"; 
				"verbose"; 
				"protmonohsp"; 
				"pipe"; 
				"offset" ; 
				"ajuste"; 
				"pnntfus" ; 
				"pnaamq"; 
				"szhsp"; 
				"szisolh" |]

	(* Dit si chaine commence par la sous-chaine debut: permet de voir chaque argument dans le fichier de configuration *)
	in let rec commence_par debut chaine =
		try 
		  let c = chaine.[String.length debut] in
		    ((String.compare (String.sub chaine 0 (String.length debut)) debut) == 0) && ((c == '=') || (c == ' ') || (c == '\t'))
		with
		    Invalid_argument _ -> false
		      

	(* Extrait la valeur du parametre: on recherche la forme NOM_PARAM[ ]*=[ ]*VAL_PARAM *)
	and extrait_valeur debut chaine = 
	    let i = ref (String.length debut) in
	      begin
		while (String.get chaine !i) = ' ' do
		  i := !i + 1;
		done;
		if (String.get chaine !i) != '=' then
		  raise (Invalid_argument "Syntax error");
		i := !i + 1;
		if !i = String.length chaine then
		  ""
		else (
		  while (String.get chaine !i) = ' ' do
		    i := !i + 1;
		  done;
		  String.sub chaine !i ((String.length chaine) - !i)
		)
	      end
		

	(* En fonction du numero du parametre, met a jour le bon champ de la structure contexte *)
	and ecrit_valeur i valeur = 
	    if i = 0 then
	      context.tmaxintron2 <- (int_of_string valeur)
	    else
	      if i = 1 then
		context.propmol2 <- (int_of_string valeur)
	      else
		if i = 2 then
		  context.gseqfile <- valeur
		else 
		  if i = 3 then
		    context.bqpfile <- valeur
		  else 
		    if i = 4 then
		      context.bqafile <- valeur
		    else 
		      if i = 5 then
			context.hsppfile <- valeur
		      else 
			if i = 6 then
			  context.hspafile <- valeur
			else 
			  if i = 7 then
			    context.fusionthresh <- int_of_string valeur
			  else 
			    if i = 8 then
			      context.alignement <- int_of_string valeur
			    else 
			      if i = 9 then
				context.matrice <- valeur
			      else 
				if i = 10 then
				  context.bplusa <- int_of_string valeur
				else 
				  if i = 11 then
				    context.diffaa <- int_of_string valeur
				  else 
				    if i = 12 then
				      context.scut <- int_of_string valeur
				    else 
				      if i = 13 then
					context.tmaxintron <- int_of_string valeur
				      else 
					if i = 14 then
					  context.nninhsp <- int_of_string valeur
					else 
					  if i = 15 then
					    context.nnsuppl <- int_of_string valeur
					  else 
					    if i = 16 then
					      context.szcdsnomet <- int_of_string valeur
					    else 
					      if i = 17 then
						context.szcdswithmet <- int_of_string valeur
					      else 
						if i = 18 then
						  context.propcdsmonoex <- int_of_string valeur
						else 
						  if i = 19 then
						    context.propmol <- int_of_string valeur
						  else 
						    if i = 20 then
						      context.sizecdsmonoex <- int_of_string valeur
						    else 
						      if i = 21 then
							context.outfile <- valeur
						      else 
							if i = 22 then
							  begin
							    if commence_par "aa " valeur then
							      context.waysigsearch <- Aamiss 
								(*
								  (int_of_string (String.sub valeur 3 (String.length valeur - 3))) 
								  modif au 27/02/06 
								*)
							    else 
							      if valeur = "mat" then
								context.waysigsearch <- Wghtmat
							      else
								raise (Invalid_argument "waysigsearch")
							  end
							else 
							  if i = 23 then
							    context.tminintron <- int_of_string valeur
							  else 
							    if i = 24 then
							      context.nbaamanq <- int_of_string valeur
							    else 
							      if i = 25 then
								context.sizecdsdiex <- int_of_string valeur
							      else 
								if i = 26 then
								  context.magfile <- valeur
								else 
								  if i = 27 then
								    context.mgtfile <- valeur
								  else 
								    if i = 28 then
								      Common.verbose := (bool_of_string valeur)
								    else 
								      if i = 29 then
									context.protmonohsp <- (bool_of_string valeur)
								      else 
									if i = 30 then
									  context.pipe <- valeur
									else 
									  if i = 31 then
									    context.offset <- (int_of_string valeur)
									  else 
									    if i = 32 then
									      context.ajuste <- (bool_of_string valeur) (* ne marche pas *)
									    else 
									      if i = 33 then
										context.pnntfus <- (int_of_string valeur)
									      else  
										if i = 34 then
										  context.pnaamq <- (int_of_string valeur)
										else 
										  if i = 35 then
										    context.szhsp <- (int_of_string valeur)
										  else 
										    context.szisolh <- (int_of_string valeur)
										      
										      
										      
										      
	(* Le corps de la fonction *)
	in let ligne = ref 0 in try
	    (* On ouvre le fichier *)
	    let file = open_in "exogean.ini" and ok = ref true in
	      Common.print_log "Config file is starting to be read\n";
	      while !ok do

		(* On parse chaque ligne *)
		let s = input_line file in
		  incr ligne;
		  if not (commence_par "#" s) && s != "" then 
		    begin
		      (* On essaie chaque parametre *)
		      for i = 0 to (Array.length params_fichier)-1 do
			let deb = params_fichier.(i) in
			  if commence_par deb s then 
			    begin
			      (* On a trouve le bon, on extrait la valeur *)
			      let v = extrait_valeur deb s in 
				begin
				  ok := true;
				  if (String.length v) != 0 then 
				    ecrit_valeur i v;
				end
			    end
		      done
		    end
	      done;
	      raise (Invalid_argument "Unknown parameter");
	  with
	      (* Puisque un fichier a obligatoirement une fin *)
	    | End_of_file -> Common.print_log "Config file is read\n";
		(* Erreur de syntaxe *)
	    | Failure s
	    | Invalid_argument s -> Common.print_error ("Syntax error ("^s^") in configuration file, at line "^(string_of_int (!ligne))^"\n");
		(* Ici, on n'a pas reussi a ouvrir ou a lire le fichier *)
	    | Sys_error _ -> Common.print_error "Error while reading the configuration file\nAn exogean.ini file has to be present where exogean is run\n";;






(***********************************************************************)
(* Lit les arguments de la ligne de commande et met a jour le contexte *)
(***********************************************************************)
let read_commandline () =

	(* Le numero de l'argument qu'on est en train de lire *)
	let argnum = ref 0 and ok = ref true in
	

	(* Cette fonction renvoie le prochain argument, mustbe dit si l'argument doit exister *)
	let getarg mustbe =
		incr argnum; 
		if !argnum < (Array.length Sys.argv) then 
			Sys.argv.(!argnum)
		else if mustbe then 
			raise Not_found 
		else
			""

	(* Lecture de chacun des arguments *)
	in try while (!ok) do
		match (getarg false) with
			| "-hspp"	-> context.hsppfile <- getarg true
			| "-hspa"	-> context.hspafile <- getarg true
			| "-gseq"	-> context.gseqfile <- getarg true
			| "-bqa"	-> context.bqafile <- getarg true
			| "-bqp"	-> context.bqpfile <- getarg true
			| "-mag"	-> context.magfile <- getarg true
			| "-mgt"	-> context.mgtfile <- getarg true
			| "-o"	-> context.outfile <- getarg true
			| "-offset"	-> context.offset <- int_of_string (getarg true)
			| "-f"	-> context.fusionthresh <- int_of_string (getarg true)
			| "-a"	-> context.alignement <- int_of_string (getarg true)
			| "-m"	-> context.matrice <- getarg true
			| "-ba"	-> context.bplusa <- int_of_string (getarg true)
			| "-daa"	-> context.diffaa <- int_of_string (getarg true)
			| "-s"	-> context.scut <- int_of_string (getarg true)
			| "-I"	-> context.tmaxintron <- int_of_string (getarg true)
			| "-I2" -> context.tmaxintron2 <- int_of_string (getarg true)
			| "-i"	-> context.tminintron <- int_of_string (getarg true)
			| "-ni"	-> context.nninhsp <- int_of_string (getarg true)
			| "-ns"	-> context.nnsuppl <- int_of_string (getarg true)
			| "-naa" -> context.nbaamanq <- int_of_string (getarg true)
			| "-sih" -> context.szisolh <- int_of_string (getarg true)
			| "-sh" -> context.szhsp <- int_of_string (getarg true)
			| "-scnm"	-> context.szcdsnomet <- int_of_string (getarg true)
			| "-scm"	-> context.szcdswithmet <- int_of_string (getarg true)
			| "-pcds"	-> context.propcdsmonoex <- int_of_string (getarg true)
			| "-pm"	-> context.propmol <- int_of_string (getarg true)
                        | "-pm2" -> context.propmol2 <- int_of_string (getarg true)
			| "-scdsm"	-> context.sizecdsmonoex <- int_of_string (getarg true)
			| "-scdsd"	-> context.sizecdsdiex <- int_of_string (getarg true)
			| "-waa"	-> context.waysigsearch <- Aamiss (*
									    (int_of_string (getarg true))
									    modif au 27/02/05 
									  *)
			| "-wmat"	-> context.waysigsearch <- Wghtmat
			| "-v"	-> Common.verbose := true
			| "-pmono"	-> context.protmonohsp <- true
			| "-protlax"    -> context.ajuste <- false
			| "-pipe"	-> context.pipe <- getarg true
			| "-pnfus" -> context.pnntfus <- int_of_string (getarg true)
			| "-pnamq" -> context.pnaamq <- int_of_string (getarg true)
			| "-h"
			| "-he"
			| "-help"	-> Printf.printf "%s\n" usage_us; exit 1; (* eprintf et --help avant *)
			| "-hf"
			| "-help-fr"	-> Printf.printf "%s\n" usage_fr; exit 1; (* eprintf et --help-fr avant *)
			| ""	-> ok := false
			| s	-> failwith s
	done;
	Common.print_log "Command line is read\n";
	with
		(* En cas d'erreur, on quitte *)
		| Not_found -> Common.print_error ("Missing parameter at the end of the command line\n"^usage_us^"\n");
		| Failure "int_of_string" -> Common.print_error ("Syntax error (incorrect int value)\n"^usage_us^"\n");
		| Failure s -> Common.print_error ("Syntax error ("^s^")\n"^usage_us^"\n");;


