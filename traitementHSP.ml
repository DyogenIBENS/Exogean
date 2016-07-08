
(*###############################################################################################

Programme traitementHSP.ml permet de faire une série de traitements sur un fichier d'HSPs 

###############################################################################################*)


(* Découpage de chacune des lignes en différents champs : utilisation des fonctions "String.sub, String.length, 
String.index_from et split_from" *)

open Common
open Donnees_base


(* ################################################################################################ *)
(* 
   map_on_other_strand est une fonction qui prend en entrée une hsp h
   et la longueur de la séquence à analyser et renvoit la meme hsp mais avec :
   1. comme reference pour la position sur le génomique le début du brin opposé
   2. le debut sur le génomique plus petit que la fin 
   Rq : la reference est le 0
*)
let map_on_other_strand n h =
  Hsp.create (Hsp.tmol h) (Hsp.nomp h) (Hsp.score h) (n-(Hsp.fing h)) (n-(Hsp.debg h)) (Hsp.debp h) (Hsp.finp h) (Hsp.brinh h) (Hsp.nbntfus h) (Hsp.cutdwn h) (Hsp.aligncomp h) (Hsp.commh h);;


(* C'est la fonction Hsp.compare *)
let triDebFinGene h1 h2 = if(Pervasives.compare (Hsp.debg h1) (Hsp.debg h2))=0 then
  Pervasives.compare (Hsp.fing h1) (Hsp.fing h2)
else Pervasives.compare (Hsp.debg h1) (Hsp.debg h2);;


let rec insert_right_place o lo comp =
  match lo with
    |[] -> [o];
    |t::q -> if (comp o t)<=0 then o::t::q
      else
	t::(insert_right_place o q comp);;



(* partage_selon_brin prend une liste d'hsps lhsp et la longueur de la séquence génomique matchée
   et renvoit un couple de listes d'hsps, la premiere sur le brin +, la deuxieme sur le brin - *)
let rec partage_selon_brin lggseq lhsp =
  match lhsp with
    |[] -> ([],[]); 
    |h::q -> let (l1,l2) = partage_selon_brin lggseq q in
	match (Hsp.brinh h) with
	    |Forward -> (insert_right_place h l1 triDebFinGene,l2);
	    |Reverse -> (l1, insert_right_place (map_on_other_strand (lggseq-1) h) l2 triDebFinGene);;



(* Cette fonction réalise les taches suivantes :
   1 - elle départage une liste d'HSPs en deux listes :
   - La première est la liste des HSPs dont la protéine ne matche qu'une fois (HSPs artefactuelles)
   - La deuxième donne une liste de listes d'HSPs (une par protéine) dont la protéine matche plus 
   d'une fois (HSPs non artefactuelles) 
   
   2 - permet d'obtenir la liste des protéines artefactuelles (celles se trouvant dans la liste des HSPs artefactuelles)
   et la liste des protéines non artefactuelles (celles se trouvant dans la liste des HSPs non artefactuelles) qui passera en paramètre de la fonction de lecture des protéines (readp) 

   - filterlh porte sur une liste d'hsp (pex taille)
   - filterh porte sur un hsp (pex taille)
   
   Attention : 
   ----------
   - cette fonction suppose les hsps des thsp triées par ordre de nom de molécule 
   - cette fonction modifie les coordonnées des hsps du brin - de sorte que la référence
   soit le début du brin -
  
   Ccl :
   -----
   à utiliser avec beaucoup de précautions
   (préférer un blusping avec segseq si possible )
*) 


let partageListeHsp filterh filterlh tHsp lgs =
  let premProt = ref (Hsp.nomp (tHsp.(0))) in
  let lhspnonartefact = ref [] in
  let lprotnonartefact = ref [] in
  let n = Array.length tHsp in    
  let nbprotall = ref 0 in
  let nbprotnonart = ref 0 in
  let ltemp = ref [] in
  let i = ref 0 in 
  let lp = ref [] and lm = ref [] in
    
    while(!i<=(n-1)) do
      while (!i<=(n-1) && (!premProt  = (Hsp.nomp (tHsp.(!i))))) do

	(* afin de filtrer les hsps à pb, pex trop petits *)
	if (filterh tHsp.(!i)) then
	  ltemp:=tHsp.(!i)::!ltemp;
	incr(i)
      done;
      
      (* pour filtrer les listes d'hsps de meme prot trop courtes *)
      if(filterlh!ltemp) then
	begin
	  lhspnonartefact:=!ltemp::!lhspnonartefact;
	  lprotnonartefact:=!premProt::(!lprotnonartefact);
	  incr nbprotnonart;
	  lp:=(fst (partage_selon_brin lgs !ltemp))::!lp;
	  lm:=(snd (partage_selon_brin lgs !ltemp))::!lm
	end;
      
      if(!i<=(n-1)) then
	begin
	  ltemp := [];
	  premProt:= Hsp.nomp (tHsp.(!i));
	end;

      incr(nbprotall);
    done;
    ((!lprotnonartefact,!nbprotall),(!lp,!lm));;	 



