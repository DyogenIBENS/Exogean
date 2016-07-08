(* printing.ml : contient toutes le fonction d'ecriture et les fonctions auxiliaires à ces dernieres
   cad en general produisant des string  *)


open Alphaprot
open Donnees_base
open Donnees
open Collection
open TraiteCollection



let brin_to_string br =
  match br with
    |Forward -> "+";
    |Reverse -> "-"



let intlist2string l offset =
  let s = ref "" in
    List.iter (fun pos -> s:=!s^(string_of_int (pos+offset+1))^",") l;
    !s


let intarray2string t offset =
	let s = ref "" in
	Array.iter (fun pos -> s:=!s^(string_of_int (pos+offset+1))^",") t;
	String.sub !s 0 (String.length !s - 1);;


let molec2string = function
	| P s -> if s="" then "NoneP" else s
	| A s -> if s="" then "NoneA" else s
	| _ ->  "None";;



let met2string b =
  if b then 
    "ATG"
  else
    "NoATG"


let stop2string b =
  if b then 
    "STOP"
  else
    "NoSTOP"


let upprotext2string b =
  if b then
    "UpProtExt"
  else
    "NoUpProtExt";;


let downprotext2string b =
  if b then
    "DownProtExt"
  else
    "NoDownProtExt";;


let pbfusion2string b = 
  if b then
    "FusionPb"
  else
    "NoFusionPb";;


let fndsigu2string b =
  if b then
    "UpSigFound"
  else
    "UpSigNotFound";;


let fndsigd2string b =
  if b then
    "DownSigFound"
  else
    "DownSigNotFound";;



(* Cette fonction ecrit un transcrit ou de l'arn au format BED *)
(* Si c'est de l'arn, on n'ecrit que les exons (donc ni les CDS, ni les codons start/stop, ni la proteine principale) *)
let print_bed lg o is_mrna nom_seq offset num_gene num_transcrit g =

	(* Le brin sous forme de chaine *)
	let brin = brin_to_string (Gene.bring g) 

	(* Coordonnee pour les brin - *)
	and compcoord c = lg-c-1

	(* Le nom de l'arn principale et de la proteine principale *)
	and nom_arn = molec2string (Gene.arnpal g) and nom_prot = molec2string (Gene.protpale g)

	(* Les noms de chaque objet*)
	in let gene_id = nom_seq^"."^(string_of_int (num_gene+1))

	and nom_transcrit = nom_seq^"."^(string_of_int (num_gene+1))^"."^(string_of_int (num_transcrit+1))
	
	and transcript_id = nom_seq^"."^(string_of_int (num_gene+1))^"."^(string_of_int (num_transcrit+1))^"__"^nom_arn^"__"^nom_prot


	(* Ecrit une ligne dans le fichier BED, les indices commencent a 1 *)
	in let print_line_bed deb_gen fin_gen deb_cds fin_cds nbex deb_ex fin_ex =  
		if is_mrna then
			Printf.fprintf o "%s\t%i\t%i\t%s\t%i\t%s\t%i\t%i\t%i\t%i\t%s\t%s\n" nom_seq (deb_gen+offset) (fin_gen+offset+1) gene_id 1000 brin (deb_cds+offset) (fin_cds+offset) 0 nbex deb_ex fin_ex
		else
			Printf.fprintf o "%s\t%i\t%i\t%s\t%i\t%s\t%i\t%i\t%i\t%i\t%s\t%s\n" nom_seq (deb_gen+offset) (fin_gen+offset+1) transcript_id 1000 brin (deb_cds+offset) (fin_cds+offset) 0 nbex deb_ex fin_ex

	(* Affiche la petite introduction en commentaire *)
	and intro deb_gen fin_gen deb_cds fin_cds =
		if not is_mrna then begin
			Printf.fprintf o "## On sequence %s : gene %i, transcript %i called \"%s\";\n" nom_seq (num_gene+1+offset) (num_transcrit+1+offset) transcript_id;
			Printf.fprintf o "##  Located on strand %s from %i to %i , CDS from %i to %i, main RNA:\"%s\" , main protein:\"%s\"\n" brin (deb_gen+1+offset) (fin_gen+1+offset) (deb_cds+1+offset) (fin_cds+1+offset) nom_arn nom_prot;
		end else
			Printf.fprintf o "## Alignment of %s on %s : strand %s from %i to %i\n" nom_arn nom_seq brin (deb_gen+1+offset) (fin_gen+1+offset);

	and liste_taille tab_deb tab_fin = Array.mapi (fun i b -> tab_fin.(i)-b+1) tab_deb

	

	(* Affiche le gene en BED *)
	in match (Gene.bring g) with
    	| Forward -> begin
						intro (Gene.gbeg g) (Gene.gend g) (Gene.gbcds g) (Gene.gecds g);
						print_line_bed (Gene.gbeg g) (Gene.gend g) (Gene.gbcds g) (Gene.gecds g) (Gene.nbex g) (intarray2string (liste_taille (Gene.gbegex g) (Gene.gendex g)) 0) (intarray2string (Gene.gbegex g) (-(Gene.gbeg g)));						
						nom_transcrit;
					end;
		| Reverse -> begin
						intro (compcoord (Gene.gend g)) (compcoord (Gene.gbeg g)) (compcoord (Gene.gecds g)) (compcoord (Gene.gbcds g));
						print_line_bed (compcoord (Gene.gend g)) (compcoord (Gene.gbeg g)) (compcoord (Gene.gecds g)) (compcoord (Gene.gbcds g)) (Gene.nbex g) (intlist2string (List.rev (Array.to_list (liste_taille (Gene.gendex g) (Gene.gbegex g)))) 0) (intarray2string (Array.map compcoord (Gene.gendex g)) (compcoord (Gene.gend g)));						
						nom_transcrit;
					end;;


(* Cette fonction ecrit un transcrit ou de l'arn au format GTF *)
(* Si c'est de l'arn, on n'ecrit que les exons (donc ni les CDS, ni les codons start/stop, ni la proteine principale) *)
let print_gtf lg o is_mrna nom_seq offset num_gene num_transcrit g =

	(* Le brin sous forme de chaine *)
	let brin = brin_to_string (Gene.bring g) 

	(* Coordonnee pour les brin - *)
	and compcoord c = lg-c-1

	(* Le nom de l'arn principale et de la proteine principale *)
	and nom_arn = molec2string (Gene.arnpal g)
	and nom_prot = molec2string (Gene.protpale g)

	(* Les noms de chaque objet*)
	in let gene_id = nom_seq^"."^(string_of_int (num_gene+1))
	and nom_transcrit = nom_seq^"."^(string_of_int (num_gene+1))^"."^(string_of_int (num_transcrit+1))
	and transcript_id = nom_seq^"."^(string_of_int (num_gene+1))^"."^(string_of_int (num_transcrit+1))^"__"^nom_arn^"__"^nom_prot
	
	(* Ecrit une ligne dans le fichier GTF, les indices commencent a 1 *)
	in let ecrit typ deb fin =
		if is_mrna then
			Printf.fprintf o "%s\tExogean\t%s\t%i\t%i\t.\t%s\t.\n" nom_seq typ (deb+1+offset) (fin+1+offset) brin
		else
			Printf.fprintf o "%s\tExogean\t%s\t%i\t%i\t.\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\";\n" nom_seq typ (deb+1+offset) (fin+1+offset) brin gene_id transcript_id

	(* Affiche la petite introduction en commentaire, ainsi que les CDS + codons start/stop si ce n'est pas de l'arn *)
	in let intro deb_gen fin_gen deb_cds fin_cds =
		if not is_mrna then begin
			Printf.fprintf o "## On sequence %s : gene %i, transcript %i called \"%s\";\n" nom_seq (num_gene+1) (num_transcrit+1) transcript_id;
			Printf.fprintf o "##  Located on strand %s from %i to %i , CDS from %i to %i, main RNA:\"%s\" , main protein:\"%s\"\n" brin (deb_gen+1+offset) (fin_gen+1+offset) (deb_cds+1+offset) (fin_cds+1+offset) nom_arn nom_prot;
			(*match (Gene.bring g) with
				| Forward -> ecrit "start_codon" deb_cds (deb_cds+2);
							ecrit "stop_codon" (fin_cds+1) (fin_cds+3);
				| Reverse -> ecrit "start_codon" (fin_cds-2) fin_cds;
							ecrit "stop_codon" (deb_cds-3) (deb_cds-1);*)
		end else
			Printf.fprintf o "## Alignment of %s on %s : strand %s from %i to %i\n" nom_arn nom_seq brin (deb_gen+1+offset) (fin_gen+1+offset);

	(* Ecrit l'exon et le cds si necessaire *)
	and print_exon_cds deb fin deb_cds fin_cds =
		ecrit "exon" deb fin;
		if deb_cds <= fin && fin_cds >= deb && not is_mrna then ecrit "CDS" (max deb deb_cds) (min fin fin_cds);

	(* Parcours des exons et appel de print_exon_cds *)
	in match (Gene.bring g) with
		| Forward -> begin
				intro (Gene.gbeg g) (Gene.gend g) (Gene.gbcds g) (Gene.gecds g);
				for j = 0 to (Gene.nbex g - 1) do begin
					print_exon_cds (Gene.gbegex g).(j) (Gene.gendex g).(j) (Gene.gbcds g) (Gene.gecds g);
				end; done;
				nom_transcrit;
				end;
		| Reverse -> begin
				intro (compcoord (Gene.gend g)) (compcoord (Gene.gbeg g)) (compcoord (Gene.gecds g)) (compcoord (Gene.gbcds g));
				for j = (Gene.nbex g - 1) downto 0 do
					print_exon_cds (compcoord (Gene.gendex g).(j)) (compcoord (Gene.gbegex g).(j)) (compcoord (Gene.gecds g)) (compcoord (Gene.gbcds g))
				done;
				nom_transcrit;
				end;;



(* Extrait le type des fichiers de sortie: bed ou gtf et renvoie la bonne fonction a appeler pour ecrire le fichier de sortie *)
(* La syntaxe du nom de fichier de sortie est "bed:nom_du_fichier" ou "gtf:nom_du_fichier" *)
(* Si le format n'est pas valide, on utilise le gtf *)
(* Si le nom n'est pas valide, on utilise "(nom_sequence)_(mois)_(annee)_(heure)_(minutes)" *)
let get_type_sortie nom_seq nom_out = 
	let (f, ext, s) = try 
		match String.sub nom_out 0 (String.index nom_out ':') with
			| "bed" -> (print_bed, "bed", Common.suffix nom_out 4)
			| "gtf" -> (print_gtf, "gtf", Common.suffix nom_out 4)
			| _ -> raise Not_found
	with
		| Invalid_argument _
		| Not_found -> Common.print_log "Error in output file specification. Falling back to defaults (GTF)\n"; (print_gtf, "gtf", "")
	in if s="" then 
		let time = Unix.localtime (Unix.time ()) in 
			(f, ext, nom_seq^"_"^(string_of_int (1+time.Unix.tm_mon))^"_"^(string_of_int time.Unix.tm_mday)^"_"^(string_of_int time.Unix.tm_hour)^"_"^(string_of_int time.Unix.tm_min))
	else 
		(f, ext, s);;




(* Pour afficher tous les genes, et tous les transcrits associes a chaque gene *)
let print_transcrits print_format_sortie nom_seq lg o offset tabgene =
	(* Pour afficher tous les transcrits d'un meme gene *)
	let rec print_transcrits_meme_gene i = function
		| [] -> []
		| t::r when Array.length t = 0 -> print_transcrits_meme_gene i r
		| t::r -> let tab = Array.to_list (Array.mapi (print_format_sortie lg o false nom_seq offset i) t) in tab::(print_transcrits_meme_gene (i+1) r)
	in print_transcrits_meme_gene 0 (SegSeq.elements (cluspinggene tabgene));;


(* Pour afficher de l'arn *)
let print_arn print_format_sortie nom_seq lg offset o =
	fun i g -> Pervasives.ignore(print_format_sortie lg o true nom_seq offset i 0 g);;


(* Pour afficher des pseudogenes *)
let print_pseudo print_format_sortie nom_seq lg offset o =
	fun i g -> Pervasives.ignore(print_format_sortie lg o false nom_seq offset i 0 g);;




(* Ecrit une proteine avec le nom du gène correspondant et le nombre d'exons codants *)
let print_prot o gene_id (nbexincds,laa) =
  let rec str = ref "" and n = ref 0 and taille_ligne = 60 in
  str := prot2string laa;
  n := String.length !str;
  Printf.fprintf o ">%s nb_exons=%i nb_aa=%i\n" gene_id nbexincds !n;
  while !n > 0 do
    if (!n <= taille_ligne) then (
    	Printf.fprintf o "%s\n\n" !str;
	n := 0;
    ) else (
        Printf.fprintf o "%s\n" (String.sub !str 0 taille_ligne);
	n := !n - taille_ligne;
	str := String.sub !str taille_ligne !n;
    )
  done;;


open Alphaadn
(* Convertit un transcrit = liste de nucleotides 
   en chaine de caracteres en appelant successivement to_stringp *)
let trans2string l =
	let str = String.create (List.length l) and i = ref 0 in
	List.iter (fun p -> str.[!i] <- to_chara p; incr i) l;
	str;;

(* Ecrit un transcrit avec le nom du gène (=transcrit) correspondant et le nombre d'exons *)
let print_tr o gene_id (nbex,tr) =
  let rec str = ref "" and n = ref 0 and taille_ligne = 60 in
  str := trans2string tr;
  n := String.length !str;
  Printf.fprintf o ">%s nb_exons=%i nb_nt=%i\n" gene_id nbex !n;
  while !n > 0 do
    if (!n <= taille_ligne) then (
    	Printf.fprintf o "%s\n\n" !str;
	n := 0;
    ) else (
        Printf.fprintf o "%s\n" (String.sub !str 0 taille_ligne);
	n := !n - taille_ligne;
	str := String.sub !str taille_ligne !n;
    )
  done;;


let print_new_line o =
  Printf.fprintf o "\n";;


(* pour imprimer une chaine suivie d'un retour chariot dans le canal o *)
let print_line o s =
  Printf.fprintf o "%s\n" s;;


(* pour imprimer i tabulations dans le canal o *)
let rec print_indent o i =
  if(i<=0) then 
    begin
      Printf.fprintf o ""
    end
  else
    begin
      Printf.fprintf o "\t"; 
      print_indent o (i-1)
    end;;
    

(* pour imprimer l'entête de la sortie html d'exogean dans le canal o *)
let print_header_html o =
  print_line o "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\" \"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\">";
  print_line o "<HTML>";
  print_line o "<HEAD>";
  print_indent o 1;
  print_line o "<TITLE>Exogean output</TITLE>";
  print_line o "</HEAD>";
  print_line o "<BODY BGCOLOR=\"FFF9D2\" LINK=\"0000CC\" VLINK=\"#330066\" ALINK=\"#330066\">";
  print_new_line o;
  print_new_line o;
  print_new_line o;
  print_line o "<center><H1>Exogean results</H1></center>";
  print_new_line o;
  print_new_line o;
  print_new_line o;
  print_line o "<A NAME=\"TOP\"></A>";;


(* pour imprimer le bas de page de la sortie html d'exogean dans le canal o *)
let print_footer_html o =
  print_new_line o;
  print_new_line o;
  print_new_line o;
  print_line o "</BODY></HTML>";;



(* Fonction qui imprime les données prises en entrée d'exogean en html *)
let print_input_exog_html gseqfile hspafile hsppfile bqpfile nomgseq lgseq o offset tug =
  print_line o "<H1>Input</H1>";
  print_line o ("<H4><U>DNA sequence Name</U>: "^(nomgseq)^"</H4>");
  print_line o ("<H4><U>DNA sequence Size</U>: "^(string_of_int lgseq)^" bases</H4>");
  print_line o "<H4><U>Source files:</U></H4>";
  print_line o "<UL>";
  print_line o "<LI>DNA:<BR>"; 
  print_line o (gseqfile^("</LI>"));
  print_line o "<LI>mRNA alignement:<BR>";
  print_line o (hspafile^("</LI>"));
  print_line o "<LI>protein alignement:<BR>"; 
  print_line o (hsppfile^("</LI>"));
  print_line o "<LI>proteins:<BR>"; 
  print_line o (bqpfile^("</LI>"));
  print_line o "</UL>";
  print_new_line o;
  print_line o "<BR>";
  print_new_line o;;


let rec somme_l = function
  |[] -> 0
  |t::q -> t+(somme_l q);;




(* Fonction utile à la fois à print_bloc_ug_in_bloc_ug et à print_ug_with_cgn,
   retourne une chaine de caractères indiquant le nom du gène avec ses caractéristiques *) 
let string_ug offset i ug =
  let st = UniteGenique.brin ug and deb = (UniteGenique.gbeg ug)+offset+1 and fin = (UniteGenique.gend ug)+offset+1 and nbtr = UniteGenique.nbgn ug and gext = UniteGenique.gext ug in
  (*  output_string stdout ((string_of_int fin)^"\n"); *)
    "GENE"^(string_of_int (i+1))^"\t"^(brin_to_string st)^"\t"^(string_of_int deb)^"\t"^(string_of_int fin)^"\t"^(string_of_int nbtr)^"\t"^(string_of_int gext);;



(* Fonction auxiliaire à print_bloc_ug : imprime une seule unité génique sans l'indication de ses transcrits *)
let print_ug_in_bloc_ug offset o i ug =
  print_line o ("<A href=\"#GENE"^(string_of_int (i+1))^"\">"^(string_ug offset i ug)^"</A><BR>");;


(* Indication des unités géniques (=gènes) les unes à la suite des autres, sans l'indication des transcrits *)
let print_bloc_ug gseqfile hspafile hsppfile bqpfile nomgseq lgseq o offset tug =
  Array.iteri (print_ug_in_bloc_ug offset o) tug;;


(* Fonction utile à print_cgn : retourne une chaine de caractères contenant le nom du transcrit
   suivi de toutes ses caractéristiques *)
let string_cg offset j cg =
  let g = CompleteGene.gn cg and nbexcds = CompleteGene.nbexincds cg and nbextot = CompleteGene.nbextot cg and lgexcds = CompleteGene.lgexcds cg and lgextot = CompleteGene.lgextot cg and metg = CompleteGene.met cg and stopg = CompleteGene.stop cg in
  let gbeg = (Gene.gbeg g)+offset+1 and gend = (Gene.gend g)+offset+1 and gbcds = (Gene.gbcds g)+offset+1 and gecds = (Gene.gecds g)+offset+1  and br = Gene.bring g and nbex = Gene.nbex g and ppale = Gene.protpale g and apal = Gene.arnpal g and upprotext = Gene.upprotext g and downprotext = Gene.downprotext g in 

    "Transcript"^(string_of_int (j+1))^"\t"^(brin_to_string br)^"\t"^(string_of_int gbeg)^"\t"^(string_of_int gend)^"\t"^(string_of_int nbextot)^"\t"^(string_of_int lgextot)^"\t"^(string_of_int gbcds)^"\t"^(string_of_int gecds)^"\t"^(string_of_int nbexcds)^"\t"^(string_of_int lgexcds)^"\t"^(molec2string apal)^"\t"^(molec2string ppale)^"\t"^(met2string metg)^"\t"^(stop2string stopg)^"\t"^(upprotext2string upprotext)^"\t"^(downprotext2string downprotext);;



let rec default_string n =
    if(n<=0) then ""
    else
      ("-"^(default_string (n-1)));;


 
(*(default_string 18^"\t"^(default_string 18))*)
let string_fnd_sig_up_down fndsigu fndsigd lab =
  match lab with
    | Unique -> ("*\t*")
    | Initial -> ("*\t"^(fndsigd2string fndsigd))
    | Internal -> ((fndsigu2string fndsigu)^"\t"^(fndsigd2string fndsigd))
    | Terminal -> ((fndsigu2string fndsigu)^"\t*");;
 


(* fonction qui donne les caractéristiques propres à un exon *)
let string_exon offset k ex =
  let gb = (Exon.gbeg ex)+offset+1 and ge = (Exon.gend ex)+offset+1 and be = Exon.brine ex and lab = Exon.lbl ex and fndsigu = Exon.fndsigu ex and fndsigd = Exon.fndsigd ex and comme = Exon.comme ex and pbfus= Exon.pbfusion ex in
    "Exon"^(string_of_int (k+1))^"\t"^(brin_to_string be)^"\t"^(string_of_int gb)^"\t"^(string_of_int ge)^"\t"^(pbfusion2string pbfus)^"\t"^(string_fnd_sig_up_down fndsigu fndsigd lab);;


(* Fonction qui imprime les caractéristiques propres à un exon *)
let print_exon offset o k exon =
  print_line o ("<LI>"^(string_exon offset k exon)^"</LI>");;





(* Fonction auxiliaire à print_ug_with_cgn : imprime un seul completegene *) 
let print_cgn offset o j cg =
  let texon = Gene.texon (CompleteGene.gn cg) in
    print_indent o 1;
    print_line o ("<LI><H6>"^(string_cg offset j cg)^"</H6></LI>");
    print_indent o 2;
    print_line o "<UL>";
    Array.iteri (print_exon offset o) texon;
    print_indent o 2;
    print_line o "</UL>";
    print_new_line o;;



(* Fonction auxiliaire à print_bloc_ug_with_cgn : imprime une seule ug avec ses completegene *) 
let print_ug_with_cgn offset o i ug =
  let tcg = UniteGenique.tcgn ug in
    print_line o ("<H4><A name=GENE"^(string_of_int (i+1))^"></A>"^(string_ug offset i ug)^"</H4>");
    print_indent o 1;
    print_line o "<UL>";
    Array.iteri (print_cgn offset o) tcg;
    print_indent o 1;
    print_line o "</UL>";
    print_new_line o;;



(* Indication des unités géniques (=gènes) de façon plus détailée :
   - une unité génique avec ses caractéristiques suivie de chacun de
   ses <<completeGenes>> (=transcrits) avec leurs caractéristiques *)
let print_bloc_ug_with_cgn gseqfile hspafile hsppfile bqpfile nomgseq lgseq o offset tug =
  Array.iteri (print_ug_with_cgn offset o) tug;;



(* Fonction qui imprime la sortie d'exogean en html *)
let print_output_exog_html gseqfile hspafile hsppfile bqpfile nomgseq lgseq o offset tug =
  let nug = Array.length tug in
  let tngn = Array.map UniteGenique.nbgn tug in
  let nbgn = somme_l (Array.to_list tngn) in
  print_line o "<H1>Output</H1>";
  print_line o ("<H4><U>#Genes</U>: "^(string_of_int nug)^", <U>#Transcripts</U>: "^(string_of_int nbgn)^", <U>#Transcript/Gene</U>: "^(string_of_float ((float_of_int nbgn)/.(float_of_int nug)))^"</H4>");
  print_bloc_ug gseqfile hspafile hsppfile bqpfile nomgseq lgseq o offset tug;
  print_new_line o;
  print_bloc_ug_with_cgn gseqfile hspafile hsppfile bqpfile nomgseq lgseq o offset tug
;;



(* Fonction qui imprime les résultats d'exogean au format html *)
let print_sortie_html gseqfile hspafile hsppfile bqpfile nomgseq lgseq o offset tug = 
  print_header_html o;
  print_input_exog_html gseqfile hspafile hsppfile bqpfile nomgseq lgseq o offset tug; 
  print_output_exog_html gseqfile hspafile hsppfile bqpfile nomgseq lgseq o offset tug; 
  print_footer_html o;;



(* fonction auxiliaire à print_cpl_ArnBmtlist_ProtBmtlist affichant un bluspmt *)
let print_Bluspmt o bmt =
  Printf.fprintf o "Bluspmt :\n";
  List.iter (fun bmm -> Printf.fprintf o "%s\n" (molec2string (Bluspmm.namem bmm))) (Bluspmt.blist bmt)
  

(* pour afficher la sortie de l'association couple arn - couple proteine du premain
   (cad la sortie de la formation des transcrits terminaux)
   avant de tout applatir sous forme d'hsps et de faire le clusping 
   Cette fonction affiche toutes les molécules composant un transcrit *)
let print_cpl_ArnBmtlist_ProtBmtlist o cplalpl =
  let (al,pl) = cplalpl in
    if((al!=[])||(pl=[])) then
      begin
	Printf.fprintf o "Composition in mRNA/EST\n*************************\n";
	List.iter (fun a -> print_Bluspmt o a) al;
	Printf.fprintf o "Composition in Proteins\n***************************\n";
	List.iter (fun p -> print_Bluspmt o p) pl;
	Printf.fprintf o "\n-------------------------------------------------------------------------\n\n"
      end


(*
let rec sstar n = match n with
  |0 -> "" ;
  |n -> "*"^(sstar (n-1))


let rec sminus n = match n with
  |0 -> "" ;
  |n -> "-"^(sminus (n-1))


let pbstrand2string b =
  match b with
    |true -> "pbbrin";
    |false -> "okbrin"


let upsig2string b =
  match b with
    |true -> "okupsig";
    |false -> "pbupsig"


let downsig2string b =
  match b with
    |true -> "okdownsig";
    |false -> "pbdownsig"
*)

(* printgene prend en entrée le nom de la séquence génomique dans acedbn la racine des noms des gènes
   à mettre dans acedb, le nom de la méthode acedb des annotations à visualiser, ainsi qu'une liste de 
   couples correspondant aux positions des exons et un nom de gène, et écrit dans le fichier de canal o 
   les lignes acedb correspondant à cette annotation *)
(*let printgene o strand i (namegene,lposex) =  
  Printf.fprintf o "%s.%i %s %i %i %i %s %s\n" namegene (i+1) (brin_to_string strand) ((fst (List.hd lposex))+1) ((snd (last lposex))+1) (List.length lposex) (intlist2string (fst (List.split lposex))) (intlist2string (snd (List.split lposex)))

  *)
