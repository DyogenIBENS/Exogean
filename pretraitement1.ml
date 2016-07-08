(* pretraitement1.ml : r�alise l'ajustement des hsps tblastn, initiales ou suppl�mentaires, independamment
   de la prot�ine de laquelle elles proviennent *)

open Alphaprot
open Common
open Donnees_base
open Donnees
open Seq



(* comparepls prend en entr�e deux tableaux pl1 et pl2 suppos�es de meme taille et renvoit le tableau (de meme 
   taille) issue de la comparaison deux � deux des elts de pl1 et pl2 avec un 1 aux positions o� le contenu 
   est identique, et avec un 0 l� ou c'est diff�rent.
   Remarque : cette fonction renverra Invalid_argument "List.map2" dans le cas o� les deux listes pl1 et pl2 
   ne sont pas de taille identique *)
let comparepls pl1 pl2 =
  try
    Array.of_list (List.map2 (fun aa1 aa2 -> if(aa1=aa2) then 1 else 0) pl1 pl2)
  with
    |Invalid_argument s -> raise (Invalid_argument s); Array.create 1 0;;




(* adequation_mots est une fonction qui prend en entre� une hsp ainsi que le g�nomique et l'ensemble 
   de prot�ines de la comparaison tblastn initiale et renvoit un tableau de 0 et de 1 en fonction de 
   l'ad�quation du mot prot�ique de h avec la traduction en acides amin�s du mot g�nomique de h.
   Attention : gseq doit etre la seq complementaire de la seq � analyser pour une hsp h de brin - 
   car � ce niveau (pr�traitement) les hsps sont trait�es brin s�par�s et le coordonn�es des hsps 
   de brin - est en r�f�rence au d�but du brin - *)
let adequation_mots h gseq bqp =
  let lgword = (Hsp.fing h)-(Hsp.debg h)+1 and lpword = (Hsp.finp h)-(Hsp.debp h)+1 in
    (* � ce niveau-la on est sens� avoir lgword = 3*lpword 
       faut-il le tester et faire qqchose en cas de probleme, ou le voir avant 
       lorsque le tblastn sera inclus dans EXOGEAN? *)
  try
    begin
      (* 
	 output_string stderr ("Coordonn�es sur prot :\ndebp = "^(string_of_int (Hsp.debp h))^"\tfinp = "^(string_of_int (Hsp.finp h))^"\n"); 
	 output_string stderr ("Le mot g�nomique de l'hsp traduit est :\n"^(prot2string (translate (Seq.to_list (Seq.sub gseq (Hsp.debg h) lgword))))^"\n");
	 output_string stderr ("Le mot prot�ique de l'hsp est :\n"^(prot2string (Seq.to_list (Seq.sub (Hashtbl.find bqp (Hsp.nomp h)) (Hsp.debp h) lpword)))^"\n\n");
      *)
      comparepls (translate (Seq.to_list (Seq.sub gseq (Hsp.debg h) lgword)))
	(Seq.to_list (Seq.sub (Hashtbl.find bqp (Hsp.nomp h)) (Hsp.debp h) lpword))
    end
  with
    |Invalid_argument s -> raise (Invalid_argument s); [||]
	(* Attention beaucoup de Not_found : est fait pour proteger le Hashtbl.find bqp 
	   mais pourrait tres bien venir du translate ou du comparepls -> voir � la main
	   sans doute � cause de l'* � la fin de certaines proteines des banques -> � rajouter 
	   dans notre alphabet?*)
    |Not_found -> Common.print_log ("La prot�ine "^(molec2string (Hsp.nomp h))^" n'est pas pr�sente dans bqp\n"); [||];;
    


(* correct_window est l'une des possibles propri�t�s p de search_window, qui consiste � reagrder ici un tableau de
   5 elements qui peuvent etre des 0 ou des 1 et de dire s'ils verifient bien la pte de d�buter (suivant dir) par 
   un 1 et de comporter au maximum 2 "0" parmi les 4 autres �lements *)
let correct_window t5 dir =
  match dir with
    |Down -> t5.(0)=1 && (t5.(1)+t5.(2)+t5.(3)+t5.(4)>=2);
    |Up -> t5.(4)=1 && (t5.(3)+t5.(2)+t5.(1)+t5.(0)>=2);;



(* search_window prend en entr�e une propri�te p (agissant sur un tableau de 5 "0" et "1"), 
   un tableau t de "0" et de "1" repr�sentant l'ad�quation entre les mots g�nomique et prot�ique 
   d'une hsp, et un sens Dir et renvoit le nombre minimum de positions � enlever de t pour dir=Down 
   (resp. Up) � partir du d�but (resp. de la fin) pour obtenir en premi�re (resp. derni�re) position 
   de t ainsi tronqu� une fenetre de 5 "0" et "1" v�rifiant la propri�t� p *)
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


(* adjust prend en entr�e une propri�t� p (de type correct_window), une hsp ainsi que les s�quences g�nomique 
   et prot�iques dont on est parti pour l'obtenir, et renvoit la liste de une ou 0 hsp constitu�e par l'hsp 
   h �ventuellement tronqu�e de qq aa selon que la propri�t� p n'est pas v�rifi�e en d�but et/ou en fin de h 
   Attention : si h est de brin -, il faudra donner comme gseq le compl�mentaire de la 
   s�quence genomique initiale *)
let adjust p gseq bqp h =
  (* let u = Common.print_log "Before search_window left\n" in *)
  let nbaa_to_substract_left = search_window p (adequation_mots h gseq bqp) Down in
  (* let u = Common.print_log "After search_window left\n" in *)
    (* Nombre d'acides amin�s � soustraire de h sur ses deux mots g�nomique et prot�ique en partant de la gauche :
       correspond � un nombre d'acides amin�s � ajouter aux positions de d�but de h *)
  let nbaa_to_substract_right = search_window p (adequation_mots h gseq bqp) Up in
    (* Nombre d'acides amin�s � soustraire de h sur ses deux mots g�nomique et prot�ique en partant de la droite :
       correspond � un nombre d'acides amin�s � enlever aux positions de fin de h *)
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




