(* db.ml :

   lecteur des hsps (read_base) et de fichiers fasta contenant des séquences 
   - génomique (readg)
   - et protéique (readp)

   les arguments pour ces deux dernières formes de lecture sont :

   beg_seq : détermine le  début d'une séquence à partir d'une chaine de caractères
   ATTENTION La ligne contenant le début de séquence est consommée 
   
   end_seq  : fonction booleenne déterminant la fin de séquences
   ATTENTION La ligne contenant la fin  de séquence est consommée.
   beg_seq  teste la même ligne que end_seq.
   
   il est à noter que pour faciliter la définition de ces fonctions
   l'exception Invalid_argument levée si on accède à des éléments hors de 
   la séquence, est capturée et retourne false.
   
   of_char : qui convertit un caractère en un type interne. Si le caractère ne correspond pas
   une exception Not_found est levée. Ce caractère sera ignoré mais la lecture continue
   
   init : valeur quelconque du type convertit 

   nb est le nombre maximum de séquences à lire 
   
   chan est le canal ouvert en lecture. Le fichier est un fichier ouvert en mode texte
 
 le résulat est une liste de séquences mais organisée dans l'ordre inverse
   du fichier.    
  
*)



open Donnees_base
open Seq
open Common

let sup l v = (String.length l > v)
let iscomment s = sup s 2 && s.[0] = '/' && s.[1]='/'
						    
let fasta = ((fun s-> (sup s 0) && s.[0]='>' ),(fun s->(sup s 0) &&  s.[0]='>' ))
	      
 (* format embl *)
let  embl = ((fun s->(sup s 1) && s.[0]='S' &&  s.[1]='Q'  ),(fun s->(sup s 1) && s.[0]='/' &&  s.[1]='/'  ))

 (*format gcg *)
let   gcg =   let f = 
		   (fun s->let l= String.length s in 
		      try 
			for  i=l-1 downto 1 do 
			  if  s.[i]='.' &&  s.[i-1]='.' then raise Exit 
			done;
			false;
		      with Exit -> true ) in (f,f)

 (* format genbank *)
let genbank = (
  (fun s-> (sup s 5) && (String.sub s 0 6 = "ORIGIN")),
  (fun s->(sup s 1) && s.[0]='/' &&  s.[1]='/'  ))


 (* format personnel *)
let  sfd =  ((fun s->(sup s 0) &&  s.[0]='{' ), (fun s->(sup s 0) &&  s.[0] = '}') )

let brinfr_of_string s =
	match s with
		| "++"
		| "+" -> Forward
		| "+-"
		| "-" -> Reverse
		|  _  -> raise Not_found;;


let brin_to_string = function
	| Forward -> "+"
	| Reverse -> "-";;

(* Cette fonction lit les lignes du fichier une par une jusqu'à la fin et les met dans une liste a l'envers *)
(* Attention, les listes de chaines de caracteres sont a l'envers car on utilise rev_split *)
let faitListe c inchan =
	let rec lit acc = 
		try
			lit ((rev_split c (input_line inchan))::acc)
		with End_of_file -> acc
	in lit []


(* record_of_line prend en entrée une ligne (lline) de la base, c'est-à-dire ici une liste 
   de chaines de caracteres, et renvoie un enregistrement de type hsp, dont les (sept) champs 
   ont été remplis par les valeurs de la ligne lline, et surtout sont bien de type voulu, 
   c'est-à-dire entiers pour la plupart (sauf brin).
   Attention : on soustrait 1 à toutes les positions : en nucléique comme en protéique, afin d'etre 
   coherent avec les calculs qui se feront en prenant comme référence 0 et non 1 comme tblastn

   tmol est Prot  ou ARNm et lline est la liste de chaines de caractères représentant la ligne psl
   renvoit une liste d'HSPs (contrairement à idem pour gff et list)
*)
let record_of_line_psl tmol lline =
	let tab = Array.of_list lline in
	if Array.length tab != 21 then Common.print_error "Erreur de syntaxe dans le fichier d'hsp : il n'y a pas 21 colonnes\n";
	
	let score = float_of_string tab.(20) (* au 13/03/06 au lieu de int_of_string *)
	and brin = brinfr_of_string tab.(12)
	and nom = Common.extrait_nom tab.(11) 0
	and sizegseq = int_of_string tab.(6) 
	and lst_blocs_taille = rev_split ',' tab.(2)
				 (* liste de couples (debutsurprot, debutsurgenomique) *)
	and lst_blocs_deb = List.combine (rev_split ',' tab.(1)) (rev_split ',' tab.(0))
	
	(* _t est la taille du bloc et le couple est le début sur la molécule et sur l'ADN *)
	in let cree t_ (deb_, debg_) = 
		let t = int_of_string t_ and deb = int_of_string deb_ and debg = int_of_string debg_ in
			
			(*Printf.printf "HSP avec %s => brin: %s, debg: %d, fing: %d, debp: %d, finp: %d\n" nom (brin_to_string brin) (debg) (debg+(if tmol = Prot then 3 else 1)*t-1) (deb) (deb+t-1);*)

			Hsp.create 
			  tmol
			  (Molecule.create (if tmol = Prot then P nom else A nom))
			  score
			  (if (tmol=Prot && brin=Reverse) then (sizegseq-(debg+t*3)) else debg)
			  (if (tmol=Prot && brin=Reverse) then (sizegseq-debg-1) else (debg+(if tmol = Prot then 3 else 1)*t-1))
			  (deb)
			  (deb+t-1)
			  brin
			  0
			  false
			  false
			  "";
	in List.rev_map2 cree lst_blocs_taille lst_blocs_deb;;



(* Cette fonction lit les hsp du fichier inchan et les traite, renvoie une liste d'Hsp
   tmol est Prot ou ARNm
   ici pour le psl on a plusieurs HSPs par ligne, d'ou le List.flatten
   renvoit une liste d'HSPs, comme les autres fonctions lit_base... *)
let lit_base_hsps_psl tmol inchan = 
	List.flatten (List.rev_map (record_of_line_psl tmol) (faitListe '\t' inchan));;





(* Idem mais pour le format gff 
   For the gff format to work for protein or mrna alignment we require:
   1) No comments (=line beginning with #
   2) tab separated from 1 to 9 and then space separated
   3) have at least 4 subfields in the 9 field
   4) expect name of protein or mRNA in 2 subfield of 9th field and begining and end of alignment
      in molecule resp in 3rd and 4th subfields of 9th field

   As gff files are usually not made like this, we remove the possibility of using gff format for alignments
*)
let record_of_line_gff tmol lline =
	let tab = Array.of_list lline in
	let tab2 = Array.of_list (rev_split ' ' (Common.clean_end_string (Common.clean_beg_string tab.(0)))) in				  
	  begin
	    if Array.length tab < 9 || Array.length tab2 < 4 then 
	      Common.print_error "Syntax error in the hsp file: you do not have the expected format\n";
	    let nom_tmp = tab2.(2) in 
	    let nom = Common.extrait_nom (String.sub nom_tmp 1 ((String.length nom_tmp) - 2)) 0 in 
	      Hsp.create 
		tmol
		(* nom replaced by nom_tmp on 02/09/2010 *)
		(Molecule.create (if tmol = Prot then P nom_tmp else A nom_tmp))
		(float_of_string tab.(3))(* au 13/03/06 au lieu de int_of_string *)
		(try ((int_of_string tab.(5))-1) with Failure _ -> Common.print_log (("genomic start ")^(tab.(5))^(" is not a valid integer position\n")); 0)
		(try ((int_of_string tab.(4))-1) with Failure _ -> Common.print_log (("genomic start ")^(tab.(5))^(" is not a valid integer position\n")); 0)
		(try ((int_of_string tab2.(1))-1) with Failure _ -> Common.print_log (("molecule start ")^(tab2.(1))^(" is not a valid integer position\n")); 0)
		(try ((int_of_string tab2.(0))-1) with Failure _ -> Common.print_log (("molecule end ")^(tab2.(0))^(" is not a valid integer position\n")); 0)
		((brinfr_of_string tab.(2)))
		0
		false
		false
		"";
	  end;;


(* Cette fonction lit les hsp du fichier inchan et les traite, renvoie une liste d'Hsp *)
let lit_base_hsps_gff tmol inchan = 
	 List.rev_map (record_of_line_gff tmol) (faitListe '\t' inchan);;


(* Idem pour le format special *)
let record_of_line_list tmol lline =
	let tab = Array.of_list lline in
	
	let nom = Common.extrait_nom tab.(6) 0 in
	
	(*Printf.printf "HSP avec %s => brin: %s, debg: %d, fing: %d, debp: %d, finp: %d\n" nom (brin_to_string (brinfr_of_string tab.(0))) ((int_of_string tab.(4))-1) ((int_of_string tab.(3))-1) ((int_of_string tab.(2))-1)   ((int_of_string tab.(1))-1);*)
	
	  Hsp.create 
	    tmol
	    (Molecule.create (if tmol = Prot then P nom else A nom))
	    (float_of_string tab.(5)) (* au 13/03/06 au lieu de int_of_string *)
	    ((int_of_string tab.(4))-1)
	    ((int_of_string tab.(3))-1)
	    ((int_of_string tab.(2))-1)
	    ((int_of_string tab.(1))-1) 
	    (brinfr_of_string tab.(0))
	    0
	    false
	    false
	    "";;

(* Cette fonction lit les hsp du fichier inchan et les traite, renvoie une liste d'Hsp *)
let lit_base_hsps_list tmol inchan = 
	 List.rev_map (record_of_line_list tmol) (faitListe ' ' inchan);;

(* Cette fonction renvoie la bonne fonction a utiliser pour lire un fichier 
   si aucun format n'est spécifié alors le format par défaut du fichier d'alignement est le psl *)
(* format:nom_du_fichier  *)
let get_type_entree s =
	try
		match String.sub s 0 (String.index s ':') with
			| "psl" -> (lit_base_hsps_psl, suffix s 4)
			| "gff" -> (lit_base_hsps_gff, suffix s 4)
			| "exf" -> (lit_base_hsps_list, suffix s 4)
			| _ -> raise (Invalid_argument "Unknown alignment format")
	with
		| Invalid_argument _
		| Not_found -> (lit_base_hsps_psl, s);;
	

(* get_nb_seq lit une premiere fois le fichier, et calcule les tailles des sequences, pour preparer readg et readp *)
let get_seq_sizes (beg_seq, end_seq) chan  = 
 
  let l = ref [] and s = ref ""  and size = ref 0 in 
    begin
		seek_in chan 0;
		s := input_line chan;
		try 
			while true do begin
				size := 0;
				while not (beg_seq !s) do   
					s := input_line chan ;
				done;
				s := input_line chan; 
	    		while not (end_seq !s) do begin
					size := !size + (String.length !s);
					s := input_line chan
				end done;
				l := !size::!l;
			end done;
		with End_of_file -> l := !size::!l;
	end;
	List.rev !l;;



(* readg : lit de la façon g (génomique) nb séquences biologiques délimitées par le format (beg_seq, end_seq)
   dans le canal d'entree chan, et dont la fonction de conversion des elements
   à partir du type char est of_char et l'element nul init;
   façon 1 = renvoit une liste de triplets (nomseq, longueur, bigseq) représentant les nb
   séquences lues dans chan;
   readg utilise en fait la liste des tailles renvoyee par get_seq_sizes pour etre plus efficace
*)
let readg (beg_seq, end_seq) of_char init chan  = 
  let l = ref [] and nom = ref "" and vseq = ref (Seq.create 0 init) and s = ref "" and i = ref 0 and lsize = ref [] in 
    begin
      lsize := get_seq_sizes (beg_seq, end_seq) chan;
      Common.print_log "I have computed the genomic sequence size\n";
      seek_in chan 0;
      s := input_line chan;
      try 
	while true do begin
	  i := 0;
	  vseq := Seq.create (List.hd !lsize) init;
	  while not (beg_seq !s) do   
	    s := input_line chan ;
	  done;
	  nom := Common.extrait_nom !s 1;
	  s := input_line chan; 
	  while not (end_seq !s) do begin
	    String.iter (fun c -> (!vseq).{!i} <- of_char c ;incr i) !s;
	    s := input_line chan
	  end done;
	  l := (!nom, Seq.length !vseq, !vseq)::!l;
	  lsize := List.tl !lsize;
	end done;
      with 
	| Failure "hd" -> (); 
	| End_of_file -> l := (!nom, Seq.length !vseq, !vseq)::!l;
	end;
    List.rev !l;;



(* readp : lit de la façon p (protéique) nbtot séquences biologiques délimitées par le format (beg_seq, end_seq)
   dans le canal d'entree chan, et dont la fonction de conversion des elements à 
   partir du type char est of_char et l'element nul init,
   et ne met dans une table de hachage que celles dont le nom est contenu dans lpwant
   façon p = renvoit une table de hachage dont la clé est le nom de la séquence
   et dont le contenu est la Seq de la séquence elle-meme
   représentant les nb séquences lues dans chan;
   readp utilise la liste des tailles fournie par get_seq_sizes
*)
let readp (beg_seq, end_seq) of_char init chan lpwant = 
 
	let tabpwant = Array.of_list lpwant in
	let nom = ref "" and vseq = ref (Seq.create 0 init) and s = ref "" and i = ref 0 and lsize = ref [] and nbpwant = Array.length tabpwant in 
	let h = Hashtbl.create nbpwant and tab_lg_prot = Array.create nbpwant 0 in
	let commit () = 
		Hashtbl.add h (Molecule.create (P !nom)) !vseq;
		tab_lg_prot.(Common.trouve_index tabpwant (P !nom)) <- Seq.length !vseq;
	in begin
		lsize := get_seq_sizes (beg_seq, end_seq) chan;
		seek_in chan 0;
		s := input_line chan;
		try 
			while true do begin
				i := 0;
				vseq := Seq.create (List.hd !lsize) init;
				while not (beg_seq !s) do   
					s := input_line chan ;
				done;
				nom := Common.extrait_nom !s 1;
				s := input_line chan; 
	    		while not (end_seq !s) do begin
					String.iter (fun c -> try ((!vseq).{!i} <- of_char c ; incr i) with Not_found -> ()) !s;
					s := input_line chan
				end done;
				if (List.mem (P !nom) lpwant) then commit ();
				lsize := List.tl !lsize;
			end done;
		with 
			| Failure "hd" -> (); 
			| End_of_file -> if (List.mem (P !nom) lpwant) then commit ();
	end;
	(h, tab_lg_prot);;


