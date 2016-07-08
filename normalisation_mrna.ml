(* normalisation_mrna.ml *)


open Seq
open Collection
open Donnees_base
open Common


(* bad_splice est une fonction qui prend en entrée une séquence génomique, une hsp h et son label et 
   renvoit true s'il n'a pas les bons siganux d'epissage (= bornes internes) 
   Attention : pour ce qui est de tolérer des bornes d'epissage non canoniques, il faudra améliorer
   le List.hd récurrent *)
let bad_splice gseq labl h = match labl with
  |Initial -> not (List.mem (Seq.sub gseq ((Hsp.fing h)+1) 2) donnorplus);
  |Internal -> (not (List.mem (Seq.sub gseq ((Hsp.debg h)-2) 2) acceptplus)) || (not (List.mem (Seq.sub gseq ((Hsp.fing h)+1) 2) donnorplus));
  |Terminal -> not (List.mem (Seq.sub gseq ((Hsp.debg h)-2) 2) acceptplus);
  |Unique -> false
  



(* gseq est la séquence génomique sur laquelle se sont alignés les ARNms, lhsp_one_mrna est la liste 
   d'hsps correspondant à un seul alignement d'ARNm et la réponse est la liste des numéros des hsps 
   de lhsp_one_mrna qui ont au moins une borne d'épissage incorrecte, ie ayant repondu true à la 
   fonction bad_splice *)
let show_no_hsp_bad_splice gseq lhsp_one_mrna = 
  let n = List.length lhsp_one_mrna and lbad_no_hsp = ref [] in
    Array.iteri (fun i h -> if (bad_splice gseq (give_label i n) h)  then
		   begin 
		     lbad_no_hsp := insert_right_place_without_redond i !lbad_no_hsp Pervasives.compare;
		   end
		)
      (Array.of_list lhsp_one_mrna);
      !lbad_no_hsp
  



(* show_no_hsp_imbric est une fonction qui prend en entrée un tableau d'hsps d'arnms et
   est auxiliaire au prétraitement 2 des hsps d'arn qui consiste en supprimer les exons
   qui ont une distance < seuil avec leur hsp suivante *)
let show_no_hsp_imbric seuil thsp_one_mrna = 
  let n = Array.length thsp_one_mrna and i = ref 0 and lbad_no_hsp = ref [] in
    while (!i<n-1) do
      if (Hsp.debg thsp_one_mrna.(!i+1) <= Hsp.fing thsp_one_mrna.(!i) +seuil) then
	lbad_no_hsp := insert_right_place_without_redond !i !lbad_no_hsp Pervasives.compare;
      incr i;
    done;
    if (List.mem (n-2) !lbad_no_hsp) then
      lbad_no_hsp := insert_right_place_without_redond (n-1) !lbad_no_hsp Pervasives.compare;
    !lbad_no_hsp



let show_no_hsp_imbric_or_hsp_bad_splice seuil gseq lhsp_one_mrna =
  tri_fusion Pervasives.compare (show_no_hsp_bad_splice gseq lhsp_one_mrna) (show_no_hsp_imbric seuil (Array.of_list lhsp_one_mrna))



(* Fonction qui retient les nos dans lhb des éléments à pb, ie contenu dans lhpb 
     A faire tourner avec i = 0 au début *)
  let rec show_no_hpb_in_b lhb lhpb i =
    match lhb with
	[] -> []
      |hb::qhb -> 
	 if (List.mem hb lhpb) then
	   i::(show_no_hpb_in_b qhb lhpb (i+1))
	 else
	   show_no_hpb_in_b qhb lhpb (i+1)
  
  

(* cut_one_to_n_list prend en entreé une liste lbad_no des nos des elts de l qui posent pb et une liste l 
   (en l'occurence ici lbad_no concerne les hsps de bornes d'epissage incorrect  ou les hsps séparés 
   par un intron trop court), et renvoit une liste de liste d'hsps correspondant aux
   sous-listes maximales de l dont les elements consécutifs ne posent pas de pb. 
   Remarque : par la suite on sera peut-etre amené à modifier egalement le nom des arnms des hsps de ces
   sous-liste (meme nom pour une meme sous-liste) *)
let cut_one_to_n_list lbad_no l =
  let i = ref 0 and ianc = ref 0 and n = List.length l and ll = ref [] and ltmp = ref [] in
    while(!i<n) do
      begin
	while (!i< n && not (List.mem !i lbad_no)) do
	  incr i;
	done;

	if (!i<n) then (* c'est que l'on est à un endroit ou on veut éliminer l'element *)
	  begin
	    ltmp := Array.to_list (Array.sub (Array.of_list l) !ianc (!i-(!ianc)));
	    (* on extrait tout ce qui précède l'élément depuis la dernière coupure, cad de ianc à i,
	       et on le met dans une liste individuelle *)
	    ll:= (!ltmp)::(!ll);
	    (* on concatene cette liste à la liste de listes existante*)
	    incr i;
	    ianc := !i;
	    (* on retient que i est désormais la derniere position de cassure *)
	  end

	else (* c'est qu'on est simplement à la fin de la liste initiale : on extrait la sous-liste 
		obtenue depuis la dernière cassure et on la concatene à ll *)
	  begin
	    ltmp := Array.to_list (Array.sub (Array.of_list l) !ianc (!i-(!ianc)));
	    ll:= (!ltmp)::(!ll);
	  end
      end
    done;
    (* apres avoir éliminé les blusps vides on retourne la liste *)
    List.rev (List.filter (fun l -> l!=[]) !ll)




  
(* fonction qui permet de transformer un arn découpé (llh = liste de listes d'hsps, chacun correspondant
   à un sous arn coupé) en un seul arn mais où on a retenu où a eu lieu la coupure du pretraitement.
   Un arn virtuel est un couple (liste d'hsps, liste d'indices de la liste d'hsps qui precede un intron virtuel) *)
let arncut_to_arn_virt llh = 
  let lflat = ref [] and lnhbefcut = ref [] and lcour = ref [] and i = ref 0 and ncumul = ref (-1) and n= List.length llh in 
    while(!i < n) do 
      lcour := List.nth llh !i; 
      ncumul := !ncumul + (List.length !lcour);
      lflat := List.append !lflat !lcour;
      lnhbefcut := (!ncumul)::(!lnhbefcut);
      incr i;
    done;
    (!lflat, (List.rev (!lnhbefcut)))




(* fonction qui permet de transformer un arn virtuel (liste d'hsps, liste d'indices de la liste 
   d'hsps qui precede un intron virtuel) en un arn découpé cad en une liste de listes d'hsps, chacun 
   correspondant à un sous arn coupé 
   Renvoit une liste de listes d'hsps (qui repertorie en fait les pbs du pretraitement) *)
let arn_virt_to_arncut av = 
  let (lh,lhofivirt) = av in
    List.map (Array.to_list) (SegSeq.elements (SegSeq.make3 (Array.of_list lh) lhofivirt))




(* supandcut prend un elt à pb et une liste et renvoit une liste de deux listes ou 
   1) l'elt à pb a été supprimé
   2) la coupure s'est fait au niveau de cet elt 
   Ex : supandcut 4 [1;2;3;4;5;6] renvoit [[1; 2; 3]; [5; 6]]
   Rq : cette fonction est générale meme si ici on l'appliquera aux blusps d'hsps arn 
   Remarque : on suppose que lh contient hpb une seule fois 
   hpb représente une hsp correspodnant en fait à un intron non épissé *)
let supandcut hpb lh = 
  if (not(List.mem hpb lh)) then 
    [lh] 
  else
    begin
      let i = ref 0 and j = ref 0 and n = List.length lh and llh = ref [] and th = Array.of_list lh in
	while (!i < n) do
	  if (th.(!i)=hpb) then
	    j:=!i;
	  incr i;
	done;
	  (Array.to_list (Array.sub th 0 !j))::[Array.to_list (Array.sub th (!j+1) (n-(!j+1)))]
    end



(* idem mais ou on élimine l'hsp à pb partout dans la liste de liste (en théorie elle n'apparait
   que dans une seule sous-liste d'hsp mais pour la rigueur on regarde dans toutes les lh de llh 
   Produit une liste de liste d'hsps sans hpb et ou la coupure se fait aux anciennes coupures et en hpb 
   Ex : supandcutinllh 4 [[1; 2; 3];[4;5;6] ; [8;4;9]] produit [[1; 2; 3]; [5; 6]; [8]; [9]] *)
let supandcutinllh hpb llh =
  List.filter (fun l -> l!=[]) (List.flatten (List.map (supandcut hpb) llh))


(* un cran au dessus on supprime dans llh (en redécoupant comme il faut) toutes les hsps à pb d'intron 
   non épissé contenus dans lhpb. On obtient une liste de liste d'hsps 
   Ex : supandcutlhpbinllh [4;6] [[1; 2; 3];[4;5;6] ; [8;4;9]] produit [[1; 2; 3]; [5]; [8]; [9]] *)
let rec supandcutlhpbinllh lhpb llh =
  match lhpb with
      [] -> llh
    |hpb::qhpb -> supandcutlhpbinllh qhpb (supandcutinllh hpb llh)
