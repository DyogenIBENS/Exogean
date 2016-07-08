(* pretraitement1.ml : réalise l'ajustement des hsps tblastn, initiales ou supplémentaires, independamment
   de la protéine de laquelle elles proviennent *)

open Alphaprot
open Common
open Donnees_base
open Donnees
open Seq



(* comparepls prend en entrée deux tableaux pl1 et pl2 supposées de meme taille et renvoit le tableau (de meme 
   taille) issue de la comparaison deux à deux des elts de pl1 et pl2 avec un 1 aux positions où le contenu 
   est identique, et avec un 0 là ou c'est différent.
   Remarque : cette fonction renverra Invalid_argument "List.map2" dans le cas où les deux listes pl1 et pl2 
   ne sont pas de taille identique *)
let comparepls pl1 pl2 =
  try
    Array.of_list (List.map2 (fun aa1 aa2 -> if(aa1=aa2) then 1 else 0) pl1 pl2)
  with
    |Invalid_argument s -> raise (Invalid_argument s); Array.create 1 0;;




(* adequation_mots est une fonction qui prend en entreé une hsp ainsi que le génomique et l'ensemble 
   de protéines de la comparaison tblastn initiale et renvoit un tableau de 0 et de 1 en fonction de 
   l'adéquation du mot protéique de h avec la traduction en acides aminés du mot génomique de h.
   Attention : gseq doit etre la seq complementaire de la seq à analyser pour une hsp h de brin - 
   car à ce niveau (prétraitement) les hsps sont traitées brin séparés et le coordonnées des hsps 
   de brin - est en référence au début du brin - *)
let adequation_mots h gseq bqp =
  let lgword = (Hsp.fing h)-(Hsp.debg h)+1 and lpword = (Hsp.finp h)-(Hsp.debp h)+1 in
    (* à ce niveau-la on est sensé avoir lgword = 3*lpword 
       faut-il le tester et faire qqchose en cas de probleme, ou le voir avant 
       lorsque le tblastn sera inclus dans EXOGEAN? *)
  try
    begin
      (* 
	 output_string stderr ("Coordonnées sur prot :\ndebp = "^(string_of_int (Hsp.debp h))^"\tfinp = "^(string_of_int (Hsp.finp h))^"\n"); 
	 output_string stderr ("Le mot génomique de l'hsp traduit est :\n"^(prot2string (translate (Seq.to_list (Seq.sub gseq (Hsp.debg h) lgword))))^"\n");
	 output_string stderr ("Le mot protéique de l'hsp est :\n"^(prot2string (Seq.to_list (Seq.sub (Hashtbl.find bqp (Hsp.nomp h)) (Hsp.debp h) lpword)))^"\n\n");
      *)
      comparepls (translate (Seq.to_list (Seq.sub gseq (Hsp.debg h) lgword)))
	(Seq.to_list (Seq.sub (Hashtbl.find bqp (Hsp.nomp h)) (Hsp.debp h) lpword))
    end
  with
    |Invalid_argument s -> raise (Invalid_argument s); [||]
	(* Attention beaucoup de Not_found : est fait pour proteger le Hashtbl.find bqp 
	   mais pourrait tres bien venir du translate ou du comparepls -> voir à la main
	   sans doute à cause de l'* à la fin de certaines proteines des banques -> à rajouter 
	   dans notre alphabet?*)
    |Not_found -> Common.print_log ("La protéine "^(molec2string (Hsp.nomp h))^" n'est pas présente dans bqp\n"); [||];;
    


(* correct_window est l'une des possibles propriétés p de search_window, qui consiste à reagrder ici un tableau de
   5 elements qui peuvent etre des 0 ou des 1 et de dire s'ils verifient bien la pte de débuter (suivant dir) par 
   un 1 et de comporter au maximum 2 "0" parmi les 4 autres élements *)
let correct_window t5 dir =
  match dir with
    |Down -> t5.(0)=1 && (t5.(1)+t5.(2)+t5.(3)+t5.(4)>=2);
    |Up -> t5.(4)=1 && (t5.(3)+t5.(2)+t5.(1)+t5.(0)>=2);;



(* search_window prend en entrée une propriéte p (agissant sur un tableau de 5 "0" et "1"), 
   un tableau t de "0" et de "1" représentant l'adéquation entre les mots génomique et protéique 
   d'une hsp, et un sens Dir et renvoit le nombre minimum de positions à enlever de t pour dir=Down 
   (resp. Up) à partir du début (resp. de la fin) pour obtenir en première (resp. dernière) position 
   de t ainsi tronqué une fenetre de 5 "0" et "1" vérifiant la propriété p *)
let search_window p t dir =
  let lgt = Array.length t and trouve=ref false in
    if (lgt=0) then
      0
    else
      begin
	match dir with
	  |Down -> let i=ref 0 in 
	      while((!i<=lgt-5) && (not !trouve)) do
		trouve:= p (Array.sub t !i 5) Down; 
		incr i;
	      done;
	      !i-1;
	  |Up -> let i=ref (lgt-1) in 
	      while((!i>=4) && (not !trouve)) do
		trouve:= p (Array.sub t (!i-4) 5) Up; 
		decr i;
	      done; 
	      lgt-1-(!i+1)
      end


(* adjust prend en entrée une propriété p (de type correct_window), une hsp ainsi que les séquences génomique 
   et protéiques dont on est parti pour l'obtenir, et renvoit la liste de une ou 0 hsp constituée par l'hsp 
   h éventuellement tronquée de qq aa selon que la propriété p n'est pas vérifiée en début et/ou en fin de h 
   Attention : si h est de brin -, il faudra donner comme gseq le complémentaire de la 
   séquence genomique initiale *)
let adjust p gseq bqp h =
  (* let u = Common.print_log "Before search_window left\n" in *)
  let nbaa_to_substract_left = search_window p (adequation_mots h gseq bqp) Down in
  (* let u = Common.print_log "After search_window left\n" in *)
    (* Nombre d'acides aminés à soustraire de h sur ses deux mots génomique et protéique en partant de la gauche :
       correspond à un nombre d'acides aminés à ajouter aux positions de début de h *)
  let nbaa_to_substract_right = search_window p (adequation_mots h gseq bqp) Up in
    (* Nombre d'acides aminés à soustraire de h sur ses deux mots génomique et protéique en partant de la droite :
       correspond à un nombre d'acides aminés à enlever aux positions de fin de h *)
  (* let u = Common.print_log "After search_window right\n" in *)
      

    if ( ((Hsp.debg h) + (3*nbaa_to_substract_left)) < ((Hsp.fing h) - (3*nbaa_to_substract_right)) ) then
      Hsp.create Prot (Hsp.nomp h) (Hsp.score h) ((Hsp.debg h) + (3*nbaa_to_substract_left)) ((Hsp.fing h) - (3*nbaa_to_substract_right)) ((Hsp.debp h) + nbaa_to_substract_left) ((Hsp.finp h) - nbaa_to_substract_right) (Hsp.brinh h) 0 false (Hsp.aligncomp h) ((Hsp.commh h)^(if ((nbaa_to_substract_left<>0) || (nbaa_to_substract_right<>0)) then "***Adjustement***\n" else ""))
    else
      Hsp.null
	
(*(Hsp.score h) - 12*(nbaa_to_substract_left+nbaa_to_substract_right)*)

let rec filtre_null lh =
  match lh with
    | [] -> []
    | t::q -> if t==Hsp.null then filtre_null q else t::(filtre_null q)



(* idem pour une liste d'hsps (List.map) *)
let adjust_lhsps p gseq bqp lhsps =
  filtre_null (List.map (adjust p gseq bqp) lhsps)




