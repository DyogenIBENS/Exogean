(*******************************************************************************)
(* Ce fichier contient des petites fonctions qui sont utilisees un peu partout *)
(*******************************************************************************)



(*************************)
(* Fonctions d'affichage *)
(*************************)

(* Indique si on est en mode verbeux 
   Attention : on ne peut avoir ce mode modifié qu'après avoir lu le fichier de config
   exogean.ini 
*)
let verbose = ref false;;

(* Permet d'afficher des logs seulement si le flag verbose est positionne *)
let print_log s =
	if !verbose then output_string stderr (("# ")^(s));
	flush stderr;;

(* Permet d'afficher un message d'erreur et de quitter *)
let print_error s =
	output_string stderr (("# ")^(s));
	flush stderr;
	exit 1;;




(*******************************************)
(* Fonctions sur les chaines de caracteres *)
(*******************************************)

(* Extrait le nom de la proteine, de la sequence ... en enlevant les commentaires
   ie: ce qu'il y a apres un espace ou un | pipe *)
let extrait_nom s deb =
	let t1 = try
		String.index s ' '
	with
		Not_found -> (String.length s) in
	let t2 = try
		min t1 (String.index s '|')
	with
		Not_found -> t1 in
	String.sub s deb (t2-deb);;


(* suffix s i retourne le suffixe d'une chaine s a partir de la position i *)
let suffix s i = 
	try
		String.sub s i ((String.length s)-i)
	with
		Invalid_argument "String.sub" -> "";;



(* split c s découpe la chaine de caractères s selon le caractere c *)
let split c s = 
	let rec split_from n = 
		try
			let p = String.index_from s n c
			in (String.sub s n (p-n)) :: (split_from (p+1))
		with Not_found -> [ suffix s n ] 
	in if s="" then [] else split_from 0 ;;

(* rev_split fait comme split mais renvoie la liste a l'envers -> recursivite terminale *)
(* rev_split enleve le dernier bout de la chaine s'il est egal a "" *)
let rev_split c s =
	let rec rev_split_from n acc =
		try
			let p = String.index_from s n c
			in rev_split_from (p+1) ((String.sub s n (p-n))::acc)
		with Not_found -> match suffix s n with 
							| "" -> acc
							| p -> p::acc
	in if s="" then [] else rev_split_from 0 [];;

(* clean_end_string eliminates all space characters present at the end of a string.
   This is useful to clean the last field of a gff file, since it can contain
   such spaces at the end. Note that it cannot contain any \t character since
   these characters are separating the 8 first fields only. 
*)
let clean_end_string s =
  let n = String.length s in
  let i= ref (n-1) in
  let unclean=ref(s.[!i]=' ') in
    while ((!unclean) && (!i>=1)) do
      decr i;
      unclean:=(s.[!i]=' ');
    done;
    if((!i=0)&&(!unclean)) then
      ""
    else
      String.sub s 0 (!i+1);;
    
(* clean_beg_string eliminates all space characters present at the begining of a string.
   This is useful to clean the last field of a gff file, since it can contain
   such spaces. Note that it cannot contain any \t character since
   these characters are separating the 8 first fields only. 
*)
let clean_beg_string s =
  let n = String.length s in
  let i= ref 0 in
  let unclean=ref(s.[!i]=' ') in
    while ((!unclean) && (!i<=n-2)) do
      incr i;
      unclean:=(s.[!i]=' ');
    done;
    if((!i=(n-1))&&(!unclean)) then
      ""
    else
      String.sub s (!i) (n-(!i));;
    



(****************************)
(* Fonctions sur les listes *)
(****************************)


(* Cette fontion renvoie le dernier element d'une liste *)
let rec last = function
	| [] -> failwith "Liste vide"
	| [x] -> x
	| _::l -> last l;;

(* on veut tri_fusionner deux listes, selon un ordre sur leurs elements donné par comp *)
(*let rec tri_fusion comp l1 l2 =
  match l1,l2 with
    |(l,[]) -> l;
    |([],l) -> l;
    |(t1::q1,t2::q2) -> if(comp t1 t2 <= 0) then 
	t1::(tri_fusion comp q1 (t2::q2)) 
      else 
	t2::(tri_fusion comp (t1::q1) q2);;*)
let tri_fusion = List.merge;;



(* ordonne prend une liste de listes d'hsps et renvoie une liste d'hsps contenant toutes les hsps de toutes 
   les sous-listes de la liste donnee en paramètre, mais triées selon leur debut dans le génomique croissant *)
(*let rec flat_and_order comp = function
	| [] -> [];
	| [l1] -> l1;
	| l1::l2::ql -> flat_and_order comp ((tri_fusion comp l1 l2)::ql);;*)
let flat_and_order comp l = List.fold_left (fun a b -> tri_fusion comp a b) [] l;;


(* Enlève la redondance selon comp dans une liste supposée triée selon comp *)
let rec remove_redond comp = function
	| [] -> []
	| [t] -> [t]
	| t1::t2::q -> if (comp t1 t2) = 0 then remove_redond comp (t2::q) else t1::(remove_redond comp (t2::q));;


(* insert_right_place_without_redond est uen fonction qui, comme son nom l'indique, place un objet o
   dans une liste d'objets lo selon l'ordre donné par la fonction comp et sans redondance *)
let rec insert_right_place_without_redond o lo comp = match lo with
	| [] -> [o];
	| t::q -> if (comp o t) < 0 then o::lo
	          else if (comp o t) > 0 then t::(insert_right_place_without_redond o q comp)
	          else t::q;; (* si les deux objets o et t sont egaux selon comp, ils sont considérés comme redondants et donc o n'est pas ajouté à la liste lo *)


(* retire la première occurence d'un élément dans une liste *)
let rec removelt h = function
	| [] -> []
	| t::q -> if t = h then q else t::(removelt h q);;

 
 
(* idem mais avec une liste d'elements à pb à retirer de l *)
(*let rec remove_all lpb l = 
	match lpb with
		| [] -> l
		| tpb::qpb -> if ((removelt tpb l)!=l) then 
						remove_all qpb (removelt tpb l)
					else
						match l with
						| [] -> []
						| t::q -> t::(remove_all qpb q);;*)

(* Supprime tous les elements a probleme d'une liste *)
let rec remove_all lpb l = 
	match lpb with
		| [] -> l
		| tpb::qpb -> remove_all qpb (removelt tpb l);;


let rec listmissnb ltrie = function
	| [] -> []
	| t::q -> if not (List.mem t ltrie) then t::(listmissnb ltrie q) else listmissnb ltrie q;;



let rec spread elt = function
	| [] -> []
	| t::q -> (elt::t)::(spread elt q);;





(******************************)
(* Fonctions sur les tableaux *)
(******************************)


(* trouve_index permet de trouver l'indice (ref : 1) d'un objet dans un tableau *)
let trouve_index tp p =
	let i = ref 0 and trouve = ref false and nbp = Array.length tp in
	while (!i<nbp && (not !trouve)) do
		trouve := (p = tp.(!i));
		incr i;
	done;
	!i-1;;




let inv_comp a b = -(Pervasives.compare a b);;


