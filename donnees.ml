(* 
   Fichier donnees.ml :définit les modules correspondant aux objets necessaires à la construction de modeles
   de gènes EXOGEAN, i.e. dans l'ordre les clusps, les eclusps, les exons, les gènes et enfin les completegenes.
   Remarque importante : les coordonnées des elements codants du brin - sont donnés en fonction du début
   du brin - dans clusp, eclusp, exon, gene, completegene (les hsps sont définies dans donnees_base.ml)
   Ce n'est qu'au niveua de l'affichage que s'opère la mise des gènes du brin - sur le brin +
   (cad refce : debut du brin +)
*)

open Alphaprot
open Common
open Seq
open Preserveorf
open Collection
open Donnees_base
open Normalisation_mrna




module Clusp = 
struct 
  type t = {
    gbeg : int ;       (* début dans le génomique *) 
    gend : int ;       (* fin dans le génomique *)
    lhsp : Hsp.t list; (* Liste des hsps constituant le clusp *)
    mset : MolSet.t;   (* ensemble de molécules, MolSet.t est le type de l'ensemble *)
    gpos : (Molecule.t, (int*int)) Hashtbl.t  ;  (* ensemble de couples (debut,fin) du clusp 
						pour chaque protéine sur le génomique *)
    ppos : (Molecule.t, (int*int)) Hashtbl.t  ;  (* ensemble de couples (debut,fin) du clusp 
						pour chaque protéine sur cette protéine *)
    brinc : brinfr;  (* brin du clusp *)
    cutdwn : bool;   (* di si on au - un hsp prot de ce clusp a le champ cutdwn à true *)
    nbhsp : int;     (* nombre d'hsps d'un clusp, autrement dit nombre de protéines (normalement) *)
    lhspp : Hsp.t list; (* Liste d'hsps principale du clusp, cad les plus longues *)
    lhspu : Hsp.t list; (* Liste d'hsps up du clusp, cad les plus up *)
    lhspd : Hsp.t list; (* Liste d'hsps down du clusp, cad les plus down *)
    commc : string array;  (* tableau de commentaires sur chacune des hsps du clusp *)
    pbbrin : bool;         (* dit si le clusp possède des hsps de brin contradictoire, auquel 
		       cas on fait un traitement à part *)
}
  let create gb ge lh ps gp pp br cd nbh lhp lhu lhd cm pb = try ({gbeg = gb; gend = ge; lhsp = lh; mset = ps; gpos = gp; ppos = pp; brinc = br; cutdwn=cd; nbhsp = nbh; lhspp = lhp; lhspu = lhu; lhspd = lhd; commc = cm; pbbrin = pb}) with
      (Invalid_argument s) -> (output_string stderr "Pb dans Clusp.create\n"); {gbeg = 0; gend = 0; lhsp = [Hsp.null]; mset = MolSet.empty; gpos = Hashtbl.create 1; ppos = Hashtbl.create 1; brinc = Forward; cutdwn=false; nbhsp = 1; lhspp = [Hsp.null]; lhspu = [Hsp.null]; lhspd = [Hsp.null]; commc = Array.create 1 ""; pbbrin = false};;

  let nullc = create 0 0 [Hsp.null] MolSet.empty (Hashtbl.create 1) (Hashtbl.create 1) Forward false 1 [Hsp.null] [Hsp.null] [Hsp.null] (Array.create 1 "") false 

  
  let compare c1 c2 = Pervasives.compare c1.gbeg c2.gbeg  (* fonction qui dit dans quel ordre sur le génomique (tous brins confondus) sont deux clusps donnés *) 
  let gbeg c = c.gbeg
  let gend c = c.gend
  let lhsp c = c.lhsp
  let mset c = c.mset
  let gpos c = c.gpos
  let ppos c = c.ppos
  let debgen c pk = fst (Hashtbl.find c.gpos pk)
  let fingen c pk = snd (Hashtbl.find c.gpos pk)
  let debprot c pk = fst (Hashtbl.find c.ppos pk)
  let finprot c pk = snd (Hashtbl.find c.ppos pk)	      
  let brinc c = c.brinc
  let cutdwn c = c.cutdwn
  let nbhsp c = c.nbhsp
  let lhspp c = c.lhspp
  let lhspu c = c.lhspu
  let lhspd c = c.lhspd
  let commc c = c.commc
  let pbbrin c = c.pbbrin

  let size c = (gend c) - (gbeg c) + 1
end 





let print_clusp o c =
  Printf.fprintf o "Contenu du clusp\n****************************************\n";
  List.iter (print_hsp o) (Clusp.lhsp c) 

 
let print_clusp2 o c =
  Printf.fprintf o "gbeg = %i; gend = %i; brinc = %s; pbcutdwn=%s\n" (Clusp.gbeg c) (Clusp.gend c) (brin_to_string (Clusp.brinc c)) (if (Clusp.cutdwn c) then "vrai" else "faux")



(* à terme il faudra peut-etre créer un module intron, avec dedans sa fonction compare, au meme titre que 
   le module Hsp *)

type intron = Real of (int * int) | Virt of (int * int) | Nullint


(* pour comparer deux intervalles i1 et i2 qui sont des couples d'entiers ordonnés (base d'intron) *)
let intervalle_comp i1 i2 =
  let (p1,p2) = i1 and (p3,p4) = i2 in
    if(p1<p3) then
      -1 
    else
      begin
	if(p1>p3) then
	  1
	else
	  (* ici forcément p1=p3*)
	  begin
	    if(p2<p4) then
	      -1
	    else
	      begin
	      if(p2>p4) then
		1
	      else
		(* ici forcément i1 = i2 *)
		begin
		  0
		end
	      end
	  end
      end;;



(* Fonction qui dit si un intron est réel (par opposition à virtuel *)
let realint i = 
  match i with
      Real(_,_) -> true
    |_ -> false;;



(* Fonction qui dit si l'hsp h inclut totalement l'intron de coordonnnées (p1,p2) *)
let intinhsp (p1,p2) h = ((Hsp.debg h) <= p1) && (p2<= (Hsp.fing h));;
  

(* Fonction qui dit si l'hsp h inclut l'un au moins des introns de la liste lintron
   de positions (p1,p2) *)
let rec ishapbintron lintron h  =
  match lintron with
      [] -> false
    |t::q -> (intinhsp t h) || (ishapbintron q h);;




(* passage d'une liste d'hsps protéique à une liste d'introns *)
let rec lhspp2lint lhp = 
  match lhp with
      [] -> []
    |[h] -> []
    |h1::h2::q -> (Real ((Hsp.fing h1),(Hsp.debg h2)))::(lhspp2lint (h2::q))




(* passage d'une liste d'hsps arns (avec introns virtuels) à une liste d'introns 
     rq : on pourrait rester en liste si on avait modifié arncut_to_arnvirt pour qu'au
     lieu de donner un couple de listes elle donne une liste de couples avec en chaque hsp
     attachée le fait qu'elle précède un intron virtuel ou pas (+ de symétrie) *)
  let lhspa2lint lhavirt = 
    let (lha,lindhaavantint) = lhavirt in
    let tha = Array.of_list lha and i = ref 0 in
    let n = Array.length tha and lint = ref [] in
      while(!i<n-1) do
	if (List.mem !i lindhaavantint) then
	  lint := (Virt(Hsp.fing tha.(!i), Hsp.debg tha.(!i+1)))::!lint
	else
	  lint := (Real(Hsp.fing tha.(!i), Hsp.debg tha.(!i+1)))::!lint;
	incr i;
      done;
      List.rev !lint;;
  
	  
      


  

(************************ Bluspmm ***********************************************)


(* Bluspmm = Blusp monomolécule, utilisé dans les étiquettes du premier graphe hétérogène
   cad celui qui relie les molécules de meme type (prot ou arn) *)
module Bluspmm = 
struct 
  type t = {
    typemol : mol;
    namem : Molecule.t; (* nom de la molécule, arn ou prot *)
    gbh : int; (* début génomique de la 1ere hsp du blusp *)
    geh : int; (* fin génomique de la 1ere hsp du blusp *)
    gbi : int; (* début génomique du 1er intron du blusp *)
    gei : int; (* fin génomique du 1er intron du blusp *)
    hl : Hsp.t list; (* liste d'hsps de meme molécule et correctement ordonné gp + avec introns pas trop longs *)
    il : intron list; (* liste d'introns associés *)
    st : brinfr;  (* brin du blusp monomolécule *)     
}
  let create tm n gbhsp gehsp gbint geint h i s = { typemol = tm; namem = n; gbh = gbhsp; geh = gehsp; gbi = gbint ; gei = geint; hl = h; il = i; st = s}

  (* Elément bluspmm null *)						    
  let nullb = { typemol = Prot; namem =  Molecule.null; gbh = 0; geh = 0; gbi = 0; gei = 0; hl = []; il = []; st = Forward};;

  (* Fonction qui dit dans quel ordre sur le génomique (tous brins confondus) sont deux bluspmm donnés *) 
  let compare b1 b2 = 
    if ((Pervasives.compare b1.gbh b2.gbh)<>0) then 
      Pervasives.compare b1.gbh b2.gbh
    else
      Pervasives.compare b1.geh b2.geh;;
   
 
  (* b2 est supposé plus grand que b1, est-ce que b2 chevauche b1? *)
  let overlap b1 b2 =
    (b2.gbh >= b1.gbh) && (b2.gbh <= b1.geh)
  
  (* b2 est supposé plus grand que b1, est-ce que b2 étend b1? *)
  let etend b1 b2 =
    (overlap b1 b2) && (b2.geh >= b1.geh) 

  (* b2 est supposé plus grand que b1, est-ce que b2 est inclus dans b1? *)
  let inclus b1 b2 =
    (overlap b1 b2) && (b2.geh < b1.geh) 

  let typemol b = b.typemol
  let namem b = b.namem
  let gbh b = b.gbh
  let geh b = b.geh
  let gbi b = b.gbi
  let gei b = b.gei
  let hl b = b.hl
  let il b = b.il
  let st b = b.st								    
  
end 


(* Pour passer d'une liste d'hsps protéiques à un bluspmm protéique *)
let lh2bluspmm lh =
    if (lh = []) then Common.print_error "Liste vide dans lh2bluspmm";
    let hprem = List.hd lh and hder = last lh in
    let li = lhspp2lint lh in
      match li with
	  [] -> Bluspmm.create Prot (Hsp.nomp hprem) (Hsp.debg hprem) (Hsp.fing hder) 0 0 lh li (Hsp.brinh hprem)
	|li -> let iprem = List.hd li and ider = last li in
	    match (iprem,ider) with
		(Real(gbi,_),Real(_,gei)) -> Bluspmm.create Prot (Hsp.nomp hprem) (Hsp.debg hprem) (Hsp.fing hder) gbi gei lh li (Hsp.brinh hprem)






(************************ Bluspmt ***********************************************)

(* Module Bluspmt : Blusp monotype = ensemble de blusps monomolécules mais de meme type
   utilisé pour le graphe quotienté des feuilles arn (avec real intron list en +)
   et des racines protéiques (avec clusp list en plus) afin de pouvoir
   tracer les liens de type ArnProtLink entre feuille arn et racine protéique *)
module Bluspmt = 
struct 
  type t = {
    typemol : mol;
    blist : Bluspmm.t list;
    gbeg : int;
    gend : int;
    strand : brinfr
  }

  let create tm bl gb ge s = { typemol = tm; blist = bl; gbeg = gb; gend = ge; strand = s}

  (* Elément bluspmt null *)						    
  let nullb = { typemol = Prot; blist = []; gbeg = 0; gend = 0; strand = Forward};;

  (* Fonction qui dit dans quel ordre sur le génomique (tous brins confondus) sont deux bluspm donnés *) 
  let compare b1 b2 = 
    if ((Pervasives.compare b1.gbeg b2.gbeg)<>0) then 
      Pervasives.compare b1.gbeg b2.gbeg
    else
      Pervasives.compare b1.gend b2.gend;;
   
 
  (* b2 est supposé plus grand que b1, est-ce que b2 chevauche b1? *)
  let overlap b1 b2 =
    (b2.gbeg >= b1.gbeg) && (b2.gbeg <= b1.gend)
  
  (* b2 est supposé plus grand que b1, est-ce que b2 étend b1? *)
  let etend b1 b2 =
    (overlap b1 b2) && (b2.gend >= b1.gend) 

  (* b2 est supposé plus grand que b1, est-ce que b2 est inclus dans b1? *)
  let inclus b1 b2 =
    (overlap b1 b2) && (b2.gend < b1.gend) 

  let typemol b = b.typemol
  let blist b = b.blist
  let gbeg b = b.gbeg
  let gend b = b.gend
  let strand b = b.strand
end




(* Pour former un bluspmt protéique à partir d'un tableau de clusps protéiques triés par ordre croissant
   Utilisé dans premain pour passer des bluspmm aux bluspmt en sachant que l'on a retenu dans les hsps
   et dans les clusps les endroits de cassure 
   (si on voulait le faire pour d'autres types de molécules il faudrait sortir tmol en argument)
   laissé au 14/03/06
*)
let tclusp2bluspmt tclusp =
  let n = Array.length tclusp and cprem = tclusp.(0) in
  let gb = Clusp.gbeg cprem and ge = Clusp.gend tclusp.(n-1) in
  let st = Clusp.brinc cprem in 
  let lhsp = List.sort Hsp.compnameprot (List.flatten (Array.to_list (Array.map Clusp.lhsp tclusp))) in
  let llhoneprot = partage_selon_prot (fun l -> true) (Array.of_list lhsp) in
    Bluspmt.create Prot (List.map lh2bluspmm llhoneprot) gb ge st;;
 

(* let llhoneprot = List.map (fun c -> ) (Array.to_list tclusp) in
    lbmmprot = List.map lh2bluspmm llhoneprot
*)





(******************************** eltassoc *****************************************)

type eltassoc = ArnLeaf of (Bluspmt.t * (intron list)) | ProtRoot of (Bluspmt.t * (Clusp.t list)) | Nulle
(* pour une feuille arn on retient en plus du blusp monotype les introns qui y ont mené 
   et pour une racine prot en plus du blusp monotype les clusps qui descendent vers la feuille *)



let eltassoccomp eas1 eas2 = match (eas1,eas2) with
    (ProtRoot (bmt1,_),ProtRoot (bmt2,_)) -> Bluspmt.compare bmt1 bmt2
  |(ArnLeaf (bmt1,_),ArnLeaf (bmt2,_)) -> Bluspmt.compare bmt1 bmt2 
  |(ProtRoot (bmt1,_),ArnLeaf (bmt2,_)) -> Bluspmt.compare bmt1 bmt2
  |(ArnLeaf (bmt1,_),ProtRoot (bmt2,_)) -> Bluspmt.compare bmt1 bmt2
  |(Nulle,Nulle) -> 0
  |(Nulle,_) -> -1
  |(_,Nulle) -> 1

let debeltassoc eas = match eas with
    ProtRoot (bmt,_) -> Bluspmt.gbeg bmt
  |ArnLeaf (bmt,_) -> Bluspmt.gbeg bmt
  |Nulle -> 0

let fineltassoc eas = match eas with
    ProtRoot (bmt,_) -> Bluspmt.gend bmt
  |ArnLeaf (bmt,_) -> Bluspmt.gend bmt
  |Nulle -> 0

let isleafarn eas = match eas with
    ArnLeaf _ -> true
  |_ -> false


let easoverlap e1 e2 = 
  if(eltassoccomp e1 e2 <= 0) then
    ((debeltassoc e2) >= (debeltassoc e1)) && ((debeltassoc e2)<= (fineltassoc e1))
  else
    ((debeltassoc e1) >= (debeltassoc e2)) && ((debeltassoc e1)<= (fineltassoc e2))
    

let compcpleas (e11,e12) (e21,e22) =
  if((e11=e21) && (e12=e22)) then
    0 
  else
    1


(* Fonctions relatives au pb de clusp de protrac à l'interieur d'un intron d'arnleaf
   en fait c'est pas exactement ce que l'on veut car on veut qu'un exon prot chevauche
   un et un seul exon arn

   Pb de ces fonctions : pour les tester c'est difficile car elles s'appliquent
   à des objets d'un type complexe
*)
  let clprotinintarn ((Real(bi,ei)),clp) = 
    let bc = Clusp.gbeg clp and ec= Clusp.gend clp in
    (bi<=bc) && (ec <= ei)
    

  let existclpininta aeas peas =
    match (aeas,peas) with
	(ArnLeaf (_,il),ProtRoot (_,cl)) -> 
	  let allicpairs = List.flatten (List.map (fun i -> List.map (fun c -> (i,c)) cl) il) in
	    begin
	      (*List.iter (fun c -> 
			   Common.print_log ("cldeb = "^(string_of_int (Clusp.gbeg c))^" - clfin = "^(string_of_int (Clusp.gend c))^"\n")
			)
		cl;
	      List.iter (fun i -> match i with
			     Real(d,f) -> Common.print_log ("Realideb = "^(string_of_int d)^" - Realifin = "^(string_of_int f)^"\n");
			   |Virt(d,f) -> Common.print_log ("Virtideb = "^(string_of_int d)^" - Virtifin = "^(string_of_int f)^"\n");
			   |Nullint -> Common.print_log ("Nullint\n");
			)
		il;*)
	      List.exists clprotinintarn allicpairs 
	    end
      |_ -> true



module BmtSet = Set.Make(Bluspmt)





(* passage d'un arn virtuel à un bluspmm *)
let avirt2bluspmm lhavirt = 
  let (lha,lindhaavantint) = lhavirt in
    try
      (let haprem = List.hd lha and hader = last lha and li = lhspa2lint lhavirt in
	 match li with
	     [] -> Bluspmm.create Rna (Hsp.nomp haprem) (Hsp.debg haprem) (Hsp.fing hader) 0 0 lha [] (Hsp.brinh haprem)
	   | li -> let iprem = List.hd li and ider = last li in
	       match (iprem,ider) with
		 | (Real (gbi,_), Real (_,gei)) -> Bluspmm.create Rna (Hsp.nomp haprem) (Hsp.debg haprem) (Hsp.fing hader) gbi gei lha li (Hsp.brinh haprem)
		 | (Real (gbi,_), Virt (_,gei)) -> Bluspmm.create Rna (Hsp.nomp haprem) (Hsp.debg haprem) (Hsp.fing hader) gbi gei lha li (Hsp.brinh haprem)
		 | (Virt (gbi,_), Real (_,gei)) -> Bluspmm.create Rna (Hsp.nomp haprem) (Hsp.debg haprem) (Hsp.fing hader) gbi gei lha li (Hsp.brinh haprem)
		 | (Virt (gbi,_), Virt (_,gei)) -> Bluspmm.create Rna (Hsp.nomp haprem) (Hsp.debg haprem) (Hsp.fing hader) gbi gei lha li (Hsp.brinh haprem))
    with
	Failure s -> Bluspmm.nullb
		

(* passage d'un bluspmm à un arn virtuel à un bluspmm *)
let bluspmm2avirt b = 
  let lh = Bluspmm.hl b and li = Bluspmm.il b and lhofivirt = ref [] in
  let n = List.length lh in
    Array.iteri 
      (fun i elt -> 
	 match elt with 
	     Virt (_,_) -> lhofivirt:=[i]::(!lhofivirt)
	   |_ -> ()) 
      (Array.of_list li);
    lhofivirt:=[n-1]::(!lhofivirt);
    (lh,List.rev (List.flatten !lhofivirt))
    


(* passage d'un bluspmm à un arn virtuel dans la convention SegSeq 
   En réalité c'est celle-ci que l'on utilise afin de pouvoir se servir de elements de SegSeq 
   par la suite *)
let bluspmm2avirts b = 
  let lh = Bluspmm.hl b and li = Bluspmm.il b and lhofivirt = ref [[0]] in
  let n = List.length lh in
    Array.iteri 
      (fun i elt -> 
	 match elt with 
	     Virt (_,_) -> lhofivirt:=[(i+1)]::(!lhofivirt)
	   |_ -> ()) 
      (Array.of_list li);
    (lh,List.rev (List.flatten !lhofivirt))
    
		

(* passe d'un blusp b à un blusp bnew dans lequel les hsps qu'il contient et qui sont 
   aussi dans lhpb ont été enlevées (pb de l'intron non épissé)
   pb avant : ici on avait complétement oublie les pbs dus aux mauvaises bornes d'épissage 
   puisque l'on se basait seulement sur les hsps et pas sur les introns dont
   certains pouvaient etre virtuels à ce stade (pb du pretraitement) 
   Remarque : un arncut = une llh 
   Remarque : p est la pté de filtre sur une lh résultant de la coupure, 
   et peut etre par exemple : (fun lh -> (List.length lh)>=2)
   Remarque : la fonction supandcutlhpbinllh est dans normalisation_mrna.ml *)
  let b2bnew p lhpb b =
    if(lhpb=[]) then
      b
    else
      begin
	let llhinit = arn_virt_to_arncut (bluspmm2avirts b) in
	let llhafterpbintron = List.filter p (supandcutlhpbinllh lhpb llhinit) in
	let lhavirt = arncut_to_arn_virt llhafterpbintron in
	  avirt2bluspmm lhavirt
      end

  




(* nouveau type, enregistrement, pour distinguer la recherche d'acides aminés manquants dans le cas 
   ou une molécule (mol) nous guide (champ mol à true) et dans le cas ou aucune molécule ne nous guide 
   (champ mol à false) --> il faudra en ce cas choisir un nombre d'acides aminés par défaut (ex : 10) 
   et en ce cas il faudra minimiser le rajout de nucléotides plutot que la difference entre ce que l'on
   ajoute et ce qui manque. hsp est l'hsp du clusp sur laquelle on s'appuie pour débuter la recherche
   de signaux alors que prot indique si on a une protéine de chaque coté. Si prot est à false
   alors il sera inutile de prevenir les stop et d'ajouter un nb de nt multiple de 3 car en ce cas
   on ne sais pas dans quelle phase on est : prot est specifiquement là pour cela *)

type nbaa2search = {nbaa : int; mol : bool; hsp : Hsp.t; prot : bool}



(* Quand on forme un eclusp c'est que l'on a deja la découpe en blusps, la différence avec un exon 
   c'est le fait que plusieurs signaux up et down lui sont assignés *)
module Eclusp = 
struct
  type t = {
    gbeg : int;   (* début sur le génomique du Eclusp *)
    gend : int;   (* fin sur le génomique du Eclusp *)
    lhsp : Hsp.t list; (* Liste des hsps constituant le eclusp *)
    brinec : brinfr;     (* brin du Eclusp *)
    lbl : label;  (* catégorie du Eclusp : unique, initial, interne ou terminal *)
    mset : MolSet.t;  (* ensemble de molécules arn et protéines, MolSet.t est le type de l'ensemble *)
    ppos : (Molecule.t, (int*int)) Hashtbl.t  ;  (* ensemble de couples (debut,fin) du eclusp sur chaque protéine *)
    protpale : Molecule.t; (* la protéine principale du blusp contenant cet eclusp, cad celle qui possede le plus de match au niveau du blusp *)
    arnpal : Molecule.t; (* l'arn principal du blusp contenant cet eclusp, cad celle qui possede le plus de match au niveau du blusp *)
    hspp : Hsp.t; (* Hsp principale du Eclusp, cad parmi les plus longues du clusp, celle de protéine protp, sinon la première *)
    hspu : Hsp.t; (* Hsp up du Eclusp, cad la plus up (idem ppale) *)
    hspd : Hsp.t; (* Hsp down du Eclusp, cad la plus down (idem ppale) *)

    lsigposu : nbaa2search * (int list); (* couple (nbre d'acide aminés à atteindre, liste des positions sur S des signaux que l'on doit chercher en amont du Clusp, la référence 5'->3' étant le brin Forward de S, afin de former le futur exon) : lsigposu contiendra donc pour un clusp de brin forward interne l'ensemble des signaux accepteurs AG recherchés en amont par rapport au brin Forward de S donc par rapport au brin du Clusp, alors que pour un clusp de brin reverse interne il contiendra les signaux donneurs CA recherchés en amont par rapport au brin Forward de S, donc en aval par rapport au brin du Clusp *)

    lsigposd : nbaa2search * (int list); (* couple (nbre d'acide aminés à atteindre, liste des positions sur S des signaux aval du Clusp) (inverse de lposigu, voir plus haut) *)

    fndsigu: bool; (* dit si on a trouvé au moins un signal up *)
    fndsigd: bool; (* dit si on a trouvé au moins un signal down *)

    commec : string array;  (* commentaires sur le Eclusp (voir clusp), chaque case du tableau commec, étant, tout comme pour les clusps, le commentaire correspondant à la ieme hsp du clusp *)
    pbbrinec : bool; 
  }
          
  (* pour comparer l'ordre de deux Eclusps sur le génomique *)  
  let compare ec1 ec2 = 
    if((Pervasives.compare ec1.gbeg ec2.gbeg)=0) then
      Pervasives.compare ec1.gend ec2.gend
    else
      Pervasives.compare ec1.gbeg ec2.gbeg

			  
  let create gb ge lh br lb ms pps pp ap hp hu hd lu ld fsu fsd cm pb =    (* pour introduire un nouvel Eclusp *)
    {
      gbeg = gb;
      gend = ge;
      lhsp = lh;
      brinec = br;
      lbl = lb;
      mset = ms;
      ppos = pps;
      protpale = pp;
      arnpal = ap;
      hspp = hp;
      hspu = hu;
      hspd = hd;
      lsigposu = lu;
      lsigposd = ld;
      fndsigu = fsu;
      fndsigd = fsd;
      commec = cm;
      pbbrinec = pb;
    }               
  let make = create

  let gbeg ec = ec.gbeg
  let gend ec = ec.gend
  let lhsp ec = ec.lhsp
  let brinec ec = ec.brinec
  let lbl ec = ec.lbl
  let mset ec = ec.mset
  let ppos ec = ec.ppos
  let debprot ec pk = fst (Hashtbl.find ec.ppos pk)
  let finprot ec pk = snd (Hashtbl.find ec.ppos pk)
  let protpale ec = ec.protpale
  let arnpal ec = ec.arnpal
  let hspp ec = ec.hspp
  let hspu ec = ec.hspu
  let hspd ec = ec.hspd
  let lsigposu ec = ec.lsigposu
  let lsigposd ec = ec.lsigposd
  let fndsigu ec = ec.fndsigu 
  let fndsigd ec = ec.fndsigd		     
  let commec ec = ec.commec
  let pbbrinec ec = ec.pbbrinec      
end






type sens = Up | Down


module Exon =
struct
  type t = {
    gbeg : int;   (* début sur le génomique de l'Exon *)
    gend : int;   (* fin sur le génomique de l'Exon *)
    lhsp : Hsp.t list; (* Liste des hsps constituant l'exon *)
    brine : brinfr;     (* brin de l'Exon *)
    lbl : label;  (* catégorie de l'Exon : unique, initial, interne ou terminal *)
    mset : MolSet.t; (* ensemble de molécules arn et protéine *) 
    protpale : Molecule.t; (* la protéine principale du blusp contenant cet exon, cad celle qui possede le plus de match au niveau du blusp *)
    arnpal : Molecule.t; (* l'arn principal du blusp contenant cet exon, cad celle qui possede le plus de match au niveau du blusp *)

    fndsigu : bool; (* dit si on a trouvé un signal up *)
    fndsigd : bool; (* dit si on a trouvé un signal down *)

    comme : string array;  (* commentaires sur l'Exon (voir Clusp), chaque case du tableau comme, étant, tout comme pour les clusps, le commentaire correspondant à la ieme hsp du clusp *)
    pbbrine : bool; 
  }
          
  (* pour comparer l'ordre de deux exons sur le génomique *)  
  let compare e1 e2 = 
    if((Pervasives.compare e1.gbeg e2.gbeg)=0) then
      Pervasives.compare e1.gend e2.gend
    else
      Pervasives.compare e1.gbeg e2.gbeg
			  
  let create gb ge lh br lb ms pp ap fsu fsd cm pb =    (* pour introduire un nouvel Exon *) 
    {
      gbeg = gb;
      gend = ge;
      lhsp = lh;
      brine = br;
      lbl = lb;
      mset = ms;
      protpale = pp;
      arnpal = ap;
      fndsigu = fsu;
      fndsigd = fsd;
      comme = cm;
      pbbrine = pb;
    }             
  
  let make = create

  let gbeg e = e.gbeg
  let gend e = e.gend
  let lhsp e = e.lhsp
  let brine e = e.brine
  let lbl e = e.lbl
  let mset e = e.mset
  let protpale e = e.protpale
  let arnpal e = e.arnpal
  let fndsigu e = e.fndsigu
  let fndsigd e = e.fndsigd
  let comme e = e.comme
  let pbbrine e = e.pbbrine  

  let existprot e = List.exists (fun h -> (Hsp.tmol h)=Prot) (lhsp e)
  let existrna e = List.exists (fun h -> (Hsp.tmol h)=Rna) (lhsp e)	
  let pbfusion e = List.exists (fun s -> (s <> "")) (Array.to_list (comme e))	    
end



(* fonction legerement plus simple que hsp2exon qui suit
   car elle donne toujours la catégorie Internal à un exon
   et non pas sa vraie catégorie *)
let hsp2exon2 h =
  Exon.create
    (Hsp.debg h) 
    (Hsp.fing h) 
    [h] 
    (Hsp.brinh h) 
    Internal
    (MolSet.singleton (Hsp.nomp h)) 
    (Hsp.nomp h) 
    (Hsp.nomp h) 
    true 
    true
    (Array.create 1 "") 
    false



type typeCDS = MetStop | BegStop | MetEnd | BegEnd;;

let typeCDS2string tcds =
  match tcds with
      MetStop -> "MetStop"
    | BegStop -> "BegStop"
    | MetEnd -> "MetEnd"
    | BegEnd -> "BegEnd" ;;


module CDS =
struct
  type t = {
    ctype : typeCDS;  
    gbeg : int;
    gend : int;
  }

  let comparebeg c1 c2 = 
    let selondeb = Pervasives.compare c1.gbeg c2.gbeg in
      if(selondeb!=0) then
	selondeb
      else
	Pervasives.compare c1.gend c2.gend;;


  let size c = (c.gend - c.gbeg)+1;;

  (* renvoit -1 si la taille de c2 est plus petite que celle de c1
     pour ranger les CDS par taille décroissante *)
  let comparesize c1 c2 =
    let size1 = size c1 and size2 = size c2 in
      Pervasives.compare size2 size1

  let create ct gb ge =    (* pour introduire un nouveau Gène = Transcrit *) 
    {
      ctype = ct;
      gbeg = gb;
      gend = ge;
    }     

  let ctype c = c.ctype
  let gbeg c = c.gbeg
  let gend c = c.gend

end


      
(* module Gene, en fait biologiquement Transcrit *)
module Gene = 
struct
  type t = {
    gbeg : int;   (* début sur le génomique du Gène = Transcrit *)
    gend : int;   (* fin sur le génomique du Gène = Transcrit *)
    gbcds : int;   (* début sur le génomique du CDS ou ORF *)
    gecds : int;   (* fin sur le génomique du CDS ou ORF *)
    texon : Exon.t array; (* tableau des exons constituant le gène *)
    bring : brinfr;     (* brin du Gène *)
    nbex : int; (* nombre d'exons du Gène *)
    gbegex : int array; (* débuts des nbex exons *)
    gendex : int array; (* fins des nbex exons *)

    protpale : Molecule.t; (* la protéine principale du gène *)
    arnpal : Molecule.t; (* l'arn principal du gène *)

    fndsigex : (bool*bool) array; (* répertorie pour chaque exon si on a trouvé les bornes up et down *)

    commg : (string array) array;
    pbbring : bool; 
  }
          
  (* pour comparer l'ordre de deux gènes sur le génomique *)  
  let compare g1 g2 = Pervasives.compare g1.gbeg g2.gbeg  
			  
  let create gb ge bcds ecds tex br nb gbe gee pp ap fse cm pb =    (* pour introduire un nouveau Gène = Transcrit *) 
    {
      gbeg = gb;
      gend = ge;
      gbcds = bcds;
      gecds = ecds;
      texon = tex;
      bring = br;
      nbex = nb;
      gbegex = gbe;
      gendex = gee;
      protpale = pp;
      arnpal = ap;
      fndsigex=fse;
      commg = cm;
      pbbring = pb;
    }             
  
  let make = create
  let null = create 0 0 0 0 [||] Forward 0 [||] [||] Donnees_base.N Donnees_base.N [||] [||] false
  
  let setbegcds begc g = {g with gbeg=begc}
  let setendcds endc g = {g with gend=endc}

  let gbeg g = g.gbeg
  let gend g = g.gend
  let gbcds g = g.gbcds
  let gecds g = g.gecds	
  let texon g = g.texon
  let bring g = g.bring
  let nbex g = g.nbex
  let gbegex g = g.gbegex
  let gendex g = g.gendex
  let protpale g = g.protpale
  let arnpal g = g.arnpal
  let fndsigex g = g.fndsigex 		   
  let commg g = g.commg
  let pbbring g = g.pbbring   

  let isnull g = (g=null)
  let existprot g = List.fold_left (fun b ex -> (b || (Exon.existprot ex))) false (Array.to_list (texon g))
  let existrna g = List.fold_left (fun b ex -> (b || (Exon.existrna ex))) false (Array.to_list (texon g))

  let upprotext g = try
    (let e = (texon g).(0) in
       (not (Exon.existrna e)) && (Exon.existprot e)) 
  with
      Invalid_argument s -> false;;
     
  let downprotext g = 
    let tex = texon g in
    let n = Array.length tex in
      try
	(let e = tex.(n-1) in
	   (not (Exon.existrna e) && (Exon.existprot e)))
      with
	  Invalid_argument s -> false;;
end







(* hsp2exon est une fonction qui transforme un hsp en exon. Cela est en particulier utile
   pour les alignements d'ARNm (prétraités du moins) dont les hsps sont deja des exons 
   (du moins les internes).
   i est la position de l'hsp h dans la collection de n hsps 
*)
let hsp2exon n i h =
  (* attention ici il faudrait voir si Hsp.nomp désigne une molécule arn ou protéique 
     avant de savoir à quel champ on assigne cette valeur (ppale ou arnpal) *)
  Exon.create (Hsp.debg h) (Hsp.fing h) [h] (Hsp.brinh h) (give_label i n) (MolSet.singleton (Hsp.nomp h)) (Hsp.nomp h) (Molecule.null) true true [|Hsp.commh h|] false




(* texon2onegene crée un seul gène à partir du tableau d'exon texon et d'un couple cpl qui indique le no 
   dans texon du premier et du dernier exon du gène + 1. Attention, si i=0, j doit etre le nbre d'exons *)
let texon2onegene texon cpl =
  let i = fst cpl and j = snd cpl in
    if(i=j) then
      begin
	Gene.create (Exon.gbeg texon.(i)) (Exon.gend texon.(i)) (Exon.gbeg texon.(i)) (Exon.gend texon.(i)) (Array.sub texon i 1) (Exon.brine texon.(i)) 1 (Array.create 1 (Exon.gbeg (texon.(i)))) (Array.create 1 (Exon.gend (texon.(i)))) (Exon.protpale texon.(i)) (Exon.arnpal texon.(i)) ([|(Exon.fndsigu texon.(i)), (Exon.fndsigd texon.(i))|]) ([|Exon.comme texon.(i)|]) (Exon.pbbrine texon.(i))
      end
    else
      begin
	let nb = j-i in
	let tgbex = Array.create nb 0 and tgeex = Array.create nb 0 and tfndsigex = Array.create nb (false,false) and tcomex = Array.create nb [||] and pb = ref false in
	  for k=0 to nb-1 do
	    tgbex.(k) <- Exon.gbeg texon.(i+k);
	    tgeex.(k) <- Exon.gend texon.(i+k);
	    tfndsigex.(k) <- ((Exon.fndsigu texon.(i+k)), (Exon.fndsigd texon.(i+k)));
	    tcomex.(k) <- Exon.comme texon.(i+k);
	    pb := !pb && (Exon.pbbrine texon.(i+k));
	  done;
	  Gene.create (Exon.gbeg (texon.(i))) (Exon.gend texon.(j-1)) (Exon.gbeg (texon.(i))) (Exon.gend texon.(j-1)) (Array.sub texon i (j-i)) (Exon.brine texon.(i)) nb tgbex tgeex (Exon.protpale texon.(i)) (Exon.arnpal texon.(i)) tfndsigex tcomex !pb
      end



(* tExon_to_tGene transforme un tableau d'exons en un tableau de gènes *)

let tExon_to_tGene texon =
  let ldeb = ref [] in
  Array.iteri (fun i e -> if ((Exon.lbl e)==Unique || (Exon.lbl e)==Initial) then ldeb:= i::!ldeb) texon;
  
  Array.map (texon2onegene texon) (Array.of_list (SegSeq.intervalles (SegSeq.setsegment (SegSeq.make2 texon) (List.rev (!ldeb)))))
  


(* transforme un gene (=transcrit) en une sequence arn (en gardant le nucleotide T), 
   cad la liste des nucleotides constitués par le raboutage de ses exons
   A terme peut etre utilisé par gene2prot *)
let gene2rna str seq cgseq cplldebendex = 
  let lg = Seq.length seq and (ldebex,lfinex) = cplldebendex in
    try
      (match str with
	 |Forward -> List.flatten (List.map2 (fun b e -> Seq.to_list (Seq.sub seq b (e-b+1))) ldebex lfinex);
	 |Reverse -> List.flatten (List.map2 (fun b e -> Seq.to_list (Seq.sub cgseq b (e-b+1))) ldebex lfinex)
      )
    with
      |Invalid_argument s -> []
    



(* prend le brin du gène, la séquence, son complementaire, le gène et le couple de listes
   correspondant aux debuts et aux fins des exons codants, et renvoit, lorsque cela est possible, 
   la proteine correspondant à la traduction du CDS sous forme de liste d'acides aminés. 
   Il faut vérifier que ça marche bien pour le brin - aussi et surtout que ce n'est pas plus simple de considérer 
   que les gènes du brin - sont tout simplement sur une autre sequence! *)
let gene2prot str seq cgseq cplldebendcdsex = 
  let lg = Seq.length seq and (ldebcdsex,lfincdsex) = cplldebendcdsex in
    try
      (match str with
	 |Forward -> translate (List.flatten (List.map2 (fun b e -> Seq.to_list (Seq.sub seq b (e-b+1))) ldebcdsex lfincdsex));
	 |Reverse -> translate (List.flatten (List.map2 (fun b e -> Seq.to_list (Seq.sub cgseq b (e-b+1))) ldebcdsex lfincdsex))
      )
    with
      |Invalid_argument s -> []
    


(* Les deux fonctions suivantes servent à afficher les bluspmm ou SMC pretraités (arns en particulier) *)
let lhsp2gene lhsp = 
  let texon = Array.of_list (List.map hsp2exon2 lhsp) in
  texon2onegene texon (0,(Array.length texon))


(* attention : ici pour un seul bluspmm pretraité et à cause des introns virtuel
   on peut avoir plusieurs "gènes" associés *)
let bpmpret2genes bpm =
  let avs = bluspmm2avirts bpm in
  let arncut = arn_virt_to_arncut avs in
  List.map lhsp2gene arncut


(* fonction qui à partir d'un gène = transcrit renvoit le couple 
   (nbexincds, couple des listes des débuts et fins des exons codants *)
let gene2cds g =
  let lbcdex = ref [] and lecdex = ref [] and i = ref 0 and ib = ref 0 and n = Gene.nbex g in
  let ie = ref (n-1) and (tbex,teex)= ((Gene.gbegex g),(Gene.gendex g)) in
  let (bcds,ecds) = ((Gene.gbcds g),(Gene.gecds g)) in
    while ((!ib <= n-1) && (teex.(!ib) < bcds)) do
      incr ib;
    done;
    while ((!ie >= 0) && (tbex.(!ie) > ecds)) do
      decr ie;
    done;
    lbcdex := bcds::(!lbcdex);
    for i=(!ib+1) to (!ie) do
      lbcdex := (tbex.(i))::(!lbcdex); 
    done;
    for i=(!ib) to (!ie-1) do
      lecdex := (teex.(i))::(!lecdex); 
    done;
    lecdex := ecds::(!lecdex); 
    (!ie-(!ib)+1,((List.rev !lbcdex),(List.rev !lecdex)))





open Alphaadn

module CompleteGene = 
struct
  type t = 
      {
	gn : Gene.t;
	prot : Alphaprot.prot_t list;  (* protéine issue du CDS *)
	trans : Alphaadn.dna_t list;  (* ARN issu du transcrit (en nucleotides) *)
	nbexincds : int;
	nbextot : int;
	lgexcds : int; (* en nucléotides *)
	lgextot : int; (* en nucléotides *)
	met : bool;
	stop : bool;
	debexcds : int array;
	finexcds : int array;
      } 

    		  
  let create g p tr necds netot lgecds lgetot m s decds fecds =    (* pour introduire un nouveau CompleteGène  *) 
    {
      gn = g;
      prot = p;
      trans=tr;
      nbexincds = necds;
      nbextot  = netot;
      lgexcds  = lgecds;
      lgextot  = lgetot;
      met = m;
      stop = s;
      debexcds=decds;
      finexcds=fecds;
    }             
  
  let make = create

  let gn cg = cg.gn
  let prot cg = cg.prot
  let trans cg = cg.trans
  let nbexincds cg = cg.nbexincds
  let nbextot cg = cg.nbextot 	
  let lgexcds cg = cg.lgexcds
  let lgextot cg = cg.lgextot 
  let met cg = cg.met   
  let stop cg = cg.stop
  let debexcds cg = cg.debexcds
  let finexcds cg = cg.finexcds

  let gbeg cg = Gene.gbeg (gn cg)
  let gend cg = Gene.gend (gn cg)
  let brincg cg = Gene.bring (gn cg)
  let compare cg1 cg2 = Gene.compare (gn cg1) (gn cg2)
end


  
  
(* transforme un gène en completegene (cf donnees.ml)
   à mettre à terme dans le module Gene précédent *)
let gene2completegene seq cgseq g = 
  let (nbexcds,(ldebexcds,lfinexcds)) = gene2cds g in
  let (ldebex,lfinex) = (Array.to_list (Gene.gbegex g), Array.to_list (Gene.gendex g)) in
  let laa = gene2prot (Gene.bring g) seq cgseq (ldebexcds,lfinexcds) in
  let tr = gene2rna (Gene.bring g) seq cgseq (ldebex,lfinex) in
  let seqendcds = Seq.sub seq ((Gene.gecds g)+1) 3 in
  let stopfnd = List.exists (fun unstop -> seqendcds=unstop) stop in
  let tdebexcds = Array.of_list ldebexcds and tfinexcds = Array.of_list lfinexcds in
  let s = ref 0 in
    CompleteGene.create 
      g 
      laa 
      tr
      nbexcds
      (Gene.nbex g) 
      (3*(List.length laa)) 
      (List.iter2 (fun b e -> s:=(!s)+e-b+1) (Array.to_list (Gene.gbegex g)) (Array.to_list (Gene.gendex g)); !s)
      (try ((List.hd laa)=Alphaprot.M) with Failure s -> false)
      stopfnd
      tdebexcds
      tfinexcds





(* Module correspondant à l'objet biologique gène, cad ensemble maximal de transcrits (=ici gènes) 
   chevauchants sur le meme brin *)
module UniteGenique =
struct
  type t = 
      {
	tcgn : CompleteGene.t array;
	gbeg : int;
	gend : int;
	brin : brinfr;
	nbgn : int; 
	gext : int; 
      } 

  let create tcg gb ge br nbg gex =
    {
      tcgn = tcg;
      gbeg = gb;
      gend = ge;
      brin = br;
      nbgn = nbg;
      gext = gex;
    }
  let make = create
  let null = create [||] 0 0 Forward 0 0

  let tcgn ug = ug.tcgn
  let gbeg ug = ug.gbeg
  let gend ug = ug.gend
  let brin ug = ug.brin	
  let nbgn ug = ug.nbgn
  let gext ug = ug.gext


  (* Transforme un tableau de complete genes (=transcrits) chevauchants 
     (ensemble maximal cluspique) en une unité génique (=gène) *)
  let tcgchev2ug tcgchev =
    try
      (let n= Array.length tcgchev in
       let lgbeg = List.sort Pervasives.compare (Array.to_list (Array.map (CompleteGene.gbeg) tcgchev)) and lgend = List.sort inv_comp (Array.to_list (Array.map (CompleteGene.gend) tcgchev)) in
       let gb = List.hd lgbeg and ge = List.hd lgend in
       let brin = CompleteGene.brincg tcgchev.(0) and gex = ge - gb + 1 in
	 create tcgchev gb ge brin n gex
      )
    with
	(* pour le cas d'un tableau avec 0 cases (n-1 = -1) *)
	Invalid_argument s -> null 
end
