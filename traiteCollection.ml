
(* traiteCollection.ml : forme les clusps, puis les blusps, les 
   eclusps (clusp connaissant son blusp avec 
   plusieurs signaux up et down) puis les exons (clusp 
   connaissant son blusp avec un seul signal up et un 
   seul signal down) et enfin les gènes *)

open Seq
open Donnees_base
open Donnees
open Collection
open Graph
open Common
open Config


(*************************************
            Fonctions pour le clusping
*************************************)


(* clusping prend en entrée une séquence segmentée 
   et rend en sortie une séquence segmentée *)
(* la fonction lclusp auxiliaure prend en entrée 
   un tableau d'hsps fusionnées ordonnées thspfusord, 
   et renvoie une liste d'entiers (le premier etant 
   forcement 0) qui représente la partition de thspfusord 
   en sous-tableaux d'hsps fusionnées ordonnées 
   Attention : clusping fait des clusps à partir d'hsps 
   tous brins confondus, ce n'est qu'ensuite, que par un 
   vote à la majorité des brins de toutes les hsps d'un 
   clusp que l'on calcule le brin d'un clusp, et que l'on 
   repertorie un probleme de brin dans un champ du clusp en question *)

let clusping shsp =
  let lclusp thsppretrait =  
    let i = ref 0 and lpos = ref [0] and curend = ref 0 and lgtab = Array.length thsppretrait in
      while (!i<lgtab) do
	curend := Hsp.fing thsppretrait.(!i);

	while (!i<lgtab && ((Hsp.debg thsppretrait.(!i)) <= !curend)) do
	  if ((Hsp.fing thsppretrait.(!i)) >= !curend) then
	    curend:= Hsp.fing thsppretrait.(!i);
	  incr i;
	done;
	if(!i<lgtab && (Hsp.debg thsppretrait.(!i)) > !curend) then
	  lpos := !i::(!lpos);

      done;
      List.rev (!lpos) in
    SegSeq.setsegment shsp (lclusp (SegSeq.tank shsp))



(* clusping sur les bluspmm *)
let cluspingbpm sbpm =
  let lcluspbpm tbpmord =
    let i = ref 0 and lpos = ref [0] and curend = ref 0 and lgtab = Array.length tbpmord in
      while (!i<lgtab) do
	curend := Bluspmm.geh tbpmord.(!i);
	
	while (!i<lgtab && ((Bluspmm.gbh tbpmord.(!i)) <= !curend)) do
	  if ((Bluspmm.geh tbpmord.(!i)) >= !curend) then
	    curend:= Bluspmm.geh tbpmord.(!i);
	  incr i;
	done;
	if(!i<lgtab && (Bluspmm.gbh tbpmord.(!i)) > !curend) then
	  lpos := !i::(!lpos);

      done;
      List.rev (!lpos) in 
    SegSeq.setsegment sbpm (lcluspbpm (SegSeq.tank sbpm))



(* clusping sur les Graph.vertlabel *)
let cluspingvlab sbpm =
  let lcluspbppm tbpmord =
    let i = ref 0 and lpos = ref [0] and curend = ref 0 and lgtab = Array.length tbpmord in
      while (!i<lgtab) do
	curend := Bluspmm.geh tbpmord.(!i);
	
	while (!i<lgtab && ((Bluspmm.gbh tbpmord.(!i)) <= !curend)) do
	  if ((Bluspmm.geh tbpmord.(!i)) >= !curend) then
	    curend:= Bluspmm.geh tbpmord.(!i);
	  incr i;
	done;
	if(!i<lgtab && (Bluspmm.gbh tbpmord.(!i)) > !curend) then
	  lpos := !i::(!lpos);

      done;
      List.rev (!lpos) in 
    SegSeq.setsegment sbpm (lcluspbppm (Array.map (fun vlabelt -> match vlabelt with
						       Protn b -> b;
						     |Arn b -> b)
					  (SegSeq.tank sbpm)))



(* clusping pour les eltassoc *)
let cluspingeas seas =
  let lclusp teas =  
    let i = ref 0 and lpos = ref [0] and curend = ref 0 and lgtab = Array.length teas in
      while (!i<lgtab) do
	curend := fineltassoc teas.(!i);

	while (!i<lgtab && ((debeltassoc teas.(!i)) <= !curend)) do
	  if ((fineltassoc teas.(!i)) >= !curend) then
	    curend:= fineltassoc teas.(!i);
	  incr i;
	done;
	if(!i<lgtab && (debeltassoc teas.(!i)) > !curend) then
	  lpos := !i::(!lpos);

      done;
      List.rev (!lpos) in
    SegSeq.setsegment seas (lclusp (SegSeq.tank seas))


(* clusping pour les genes, necessaire pour la sortie dans le format GTF *)
(* C'est comme les autres cluspings sauf qu'on doit traiter les brins + et - separement
   attention : normalement présuppose que les gènes sont dans le bon order génomique *)
let cluspinggene tabgene =
	let i = ref 0 and lpos = ref [0] and curend = ref 0 and lgtab = Array.length tabgene and curbr = ref Forward in
	while (!i < lgtab) do

		curend := Gene.gend tabgene.(!i);
		curbr := Gene.bring tabgene.(!i);

		while (!i<lgtab) && ((Gene.bring tabgene.(!i)) == !curbr) && ((Gene.gbeg tabgene.(!i)) <= !curend) do
			if (Gene.gend tabgene.(!i)) >= !curend then
				curend:= Gene.gend tabgene.(!i);
			incr i;
		done;
		if !i<lgtab then
			lpos := !i::(!lpos);

	done;
	SegSeq.make3 tabgene (List.rev (!lpos));;


(* clusping pour les genes, necessaire pour la sortie dans le format GTF *)
(* C'est comme les autres cluspings sauf qu'on doit traiter les brins + et - separement
   attention : normalement présuppose que les gènes sont dans le bon order génomique *)
let cluspingcgene tabcgene =
	let i = ref 0 and lpos = ref [0] and curend = ref 0 and lgtab = Array.length tabcgene and curbr = ref Forward in
	while (!i < lgtab) do

		curend := CompleteGene.gend tabcgene.(!i);
		curbr := CompleteGene.brincg tabcgene.(!i);

		while (!i<lgtab) && ((CompleteGene.brincg tabcgene.(!i)) == !curbr) && ((CompleteGene.gbeg tabcgene.(!i)) <= !curend) do
			if (CompleteGene.gend tabcgene.(!i)) >= !curend then
				curend:= CompleteGene.gend tabcgene.(!i);
			incr i;
		done;
		if !i<lgtab then
			lpos := !i::(!lpos);

	done;
	SegSeq.make3 tabcgene (List.rev (!lpos));;




(*********************************************
  Fonctions pour convertir les tableaux d'hsps 
  pré-traitées chevauchantes en type Clusp proprement dit 
******************************************************)


(* Détermine la liste d'hsps principales d'un clusp, 
   cad les plus longues hsps du clusp *)

let lhsppale thspchev = 
  let i = ref 0 and n = (Array.length thspchev) in
    Array.sort (fun h1 h2 -> Pervasives.compare ((Hsp.fing h2) -(Hsp.debg h2)) ((Hsp.fing h1) -(Hsp.debg h1))) thspchev;
    while (!i<n && ((Hsp.fing thspchev.(!i))-(Hsp.debg thspchev.(!i)))=((Hsp.fing thspchev.(0))-(Hsp.debg thspchev.(0)))) do
      incr i;
    done;
    Array.to_list (Array.sub thspchev 0 !i)
    


(*Détermine les hsps qui ont donné leur borne gauche au clusp *)
let lhspleft thspchev = 
  let i = ref 0 and n = Array.length thspchev in
    Array.sort (fun h1 h2 -> Pervasives.compare (Hsp.debg h1) (Hsp.debg h2)) thspchev;
    while (!i<n && ((Hsp.debg thspchev.(!i))=(Hsp.debg thspchev.(0)))) do
      incr i;
    done;
      Array.to_list (Array.sub thspchev 0 !i)
   
 


(*Détermine les hsps qui ont donné leur borne droite au clusp *)
let lhspright thspchev = 
  let i = ref 0 and n = Array.length thspchev in
    Array.sort (fun h1 h2 -> Pervasives.compare (Hsp.fing h2) (Hsp.fing h1)) thspchev;
    while (!i<n && ((Hsp.fing thspchev.(!i))=(Hsp.fing thspchev.(0)))) do
      incr i;
    done;
      Array.to_list (Array.sub thspchev 0 !i)

    


(* renvoit l'ensemble des protéines d'un clusp *)
let unionprot thspchev =
  let protset = ref (MolSet.empty) in
    Array.iter (fun p -> (protset:= MolSet.add p !protset)) (Array.map (fun hsp -> Hsp.nomp hsp) thspchev);
    !protset


(* renvoit l'ensemble de début,fin, pour chaque 
   protéine, de son hsp sur le génomique, sous forme de table de hachage.
   Attention : on suppose qu'il n'y a pas plus 
   de 20 protéines pour un clusp donné *)

let debfingen thspchev =
  let n = Array.length thspchev in  (* nombre d'hsps du clusp, normalement égal au nombre de protéines du clusp *)
  let h = Hashtbl.create n in
    Array.iter (fun hsp -> Hashtbl.add h (Hsp.nomp hsp) (Hsp.debg hsp,Hsp.fing hsp)) thspchev;
    h




(* renvoit l'ensemble de début,fin, pour chaque protéine
   de son hsp sur cette protéine, sous forme de table de hachage.
   Attention : on suppose qu'il n'y a pas plus 
   de 20 protéines pour un clusp donné *)

let debfinprot thspchev =
  let n = Array.length thspchev in  (* nombre d'hsps du clusp, normalement égal au nombre de protéines du clusp *)
  let h = Hashtbl.create n in
    Array.iter (fun hsp -> Hashtbl.add h (Hsp.nomp hsp) (Hsp.debp hsp,Hsp.finp hsp)) thspchev;
    h




(* renvoie vrai si tous les elements d'une liste sont identiques, 
   peut servir pour détecter des brins differents dans 
   les différtentes hsps d'un clusp *)

let all_id l = 
  if l=[] then 
    true
  else 
    List.fold_left (fun bool x -> bool && (x=(List.hd l))) true l



(* Détermine le brin d'un clusp par un vote à la majorité 
   de l'ensemble des hsps qui le constituent *)
let brinclusp thspchev = 
  let part_hsp_selon_brin = List.partition (fun h -> (Hsp.brinh h) = Forward) (Array.to_list thspchev) in 
    if (List.length (fst part_hsp_selon_brin) >= List.length (snd part_hsp_selon_brin)) then 
      Forward 
    else 
      Reverse
  


(* Détermine le tableau des commentaires sur un clusp, 
   on saura ainsi de quelle hsp ils proviennent *)
let commclusp thspchev = 
  Array.map Hsp.commh thspchev



(* Détermine le booleen qui dit si le clusp possède des hsps de brin different
   attention : ici il y avait un bug dans v1 car fun h -> Hsp.brinh au lieu de Hsp.brinh *)
let pbbrinclusp thspchev =
  not (all_id (Array.to_list (Array.map Hsp.brinh thspchev)))




(* crée un clusp à partir d'un tableau d'hsps chevauchantes 
   maximal, Attention : ne fait PAS en sorte que le mot 
   génomique du clusp soit mod 3, meme si dans l'ideal 
   cela est important pour la recherche de signaux : en 
   effet il faut tenir compte des possible frameshifts 
   sur le génomique que l'on annote
   (et détectés lors de la fusion du pretraitement2) 
   qui entraineront obligatoirement la formation de 
   modèles de gènes non mod3 au final. Le mieux que  
   l'on puisse faire est de documenter ces frameshifts 
   dans un deuxieme fichier de sortie (contenant ce 
   genre de commentaires)*)
let thspchev_to_clusp thspchev = 
  try 
    (if (thspchev=[||]) then
       begin
	 (*Common.print_log ("je fais un clusp nul dans thspchev_to_clusp\n");*)
	 Clusp.nullc
       end
     else
       begin
	 let lhp = lhsppale thspchev and lhl = lhspleft thspchev and lhr = lhspright thspchev 
								 and nbhsp = Array.length thspchev in
	 let cgbegi = Hsp.debg (List.hd lhl) and cgendi = Hsp.fing (List.hd lhr) in
	   Clusp.create cgbegi cgendi (Array.to_list thspchev) (unionprot thspchev) (debfingen thspchev) (debfinprot thspchev) (brinclusp thspchev) (List.exists (fun h -> Hsp.cutdwn h) (Array.to_list thspchev)) nbhsp lhp lhl lhr (commclusp thspchev) (pbbrinclusp thspchev)
       end)
  with Invalid_argument s -> (output_string stderr "Pb dans thspchev_to_clusp\n"); Clusp.nullc;;



(***********************************
      Fonctions pour le blusping
************************************)


(* pour transformer un tableau d'ensembles de protéines en ensemble de protéines *)
let set_of_array aprotset = 
  let s = ref (MolSet.empty) in
    Array.iter (fun aps -> s:= MolSet.union !s aps) aprotset;
    !s

(* ensprot prend un tableau de clusps (tclusp), une 
   position de debut et une longueur dans ce tableau 
   (comme Array.sub),  et renvoit l'ensemble des protéines 
   de tous les clusps compris entre deb et 
   long+deb-1 dans le tableau tclusp *)

let ensprot tclusp deb long = set_of_array (Array.map (fun c -> Clusp.mset c) (Array.sub tclusp deb long))  


let ensprotparclusp tclusp = Array.map (fun c -> MolSet.elements (Clusp.mset c)) tclusp


(* la fonction blusping1 permet de constituer des blusps 
   de catégorie 1 à partir des clusps : c'est une partition 
   de l'ensemble des clusps en objets nommés blusps ayant comme propriétés :
   - tous les clusps d'un blusp sont adjacents sur le génomique
   - au sein d'un beta-blusp d'un blusp (cad tout sous-ensemble 
   {ci,ci+1,..., ci+beta-1} de beta clusps d'un blusp adjacents 
   sur le genomique), le premier clusp (ci) doit avoir une protéine 
   commune avec au moins l'un des beta-1 clusps suivants (ci+1,...,cbeta+i-1)
  
   Cette fonction prend en entree un tableau de clusps et renvoie 
   une liste d'entiers qui sont les positions de début de chaque 
   blusp de catégorie 1 dans ce tableau. Elle prend egalement une 
   propriété permettant d'avancer dans le tableau de clusps 
   (qui renvoie un booleen) *)
(* Remarque : on considère toujours le clusp i (et donc sa coupure) 
   à la lumière des beta clusps precedents, au début on prendra 
   beta=2 cad que l'on ne tolerera qu'un seul clusp à pb entre 
   deux clusps <<normaux>> -> pose pb dasn la recherche de signaux -> beta= 1*)

let blusping1 sclusp beta = 
  let lblusp tclusp beta =
    let j=ref 0 and i=ref 0 and nc = Array.length tclusp and stop = ref false and pprec = ref MolSet.empty and l = ref [0] in
      
      while(!j<=nc-beta) do
	
	pprec := ensprot tclusp !j beta;
	stop:=false;
	i:=!j+beta;
      
	while(!i<=nc-beta && not (!stop)) do
	  stop:= (MolSet.is_empty (MolSet.inter !pprec (Clusp.mset tclusp.(!i))));
	  
	  if (not (!stop)) then
	    begin
	      pprec:=MolSet.union (MolSet.diff !pprec (Clusp.mset tclusp.(!i-beta))) (Clusp.mset tclusp.(!i));
	      incr i;
	    end;
	  
	done;

	if (!stop) then 
	  l:=(!i)::(!l);

	j:=!i;
	
      done;
      List.rev !l in
    SegSeq.setsegment sclusp (lblusp (SegSeq.tank sclusp) beta)





(* idbrin permet de savoir si en une position k d'un blusp, on a le meme 
   brin qu'en le premier clusp i0 du blusp courant. Aide au blusping no 2, 
   de type blusping2 (tintmax sert seulement à ce que toutes les ptes p
   utilisées par blusping2 aient meme signature (à savoir idbrin, gporder et tintronok) *)

let idbrin k i0 tintmax tclusp beta =
  Clusp.brinc tclusp.(k) = Clusp.brinc tclusp.(i0)




(* La fonction finale gporder prend une position k, une position i et un entier tintmax
   (ne signifiant rien, seulement la pour que gporder, idbrin et tintronok
   aient la meme signature), un tableau de clusps et un entier 
   beta, et regarde si les clusps qui correspondent aux 
   positions k et k-beta, soit n'ont pas de protéine commune, 
   soit, pour chaque protéine qu'ils ont en commun, ont les hsps 
   correspondant aux proteines qu'ils ont en commun dans un ordre 
   coherent dans le génomique et la protéine (en sachant que sur 
   le brin -  le premier clusp devient le dernier et reciproquement.
   Attention : k doit etre plus grand que beta. Aide au blusping no 3, 
   de type blusping2 

   Attention : 
   - i est la position du premier clusp, ne sert que pour idbrin,
   - tinitmax est la taille d'intron max, ne sert que pour tintronok
*)

let gporderf k i tintmax tclusp beta =
  let cprem = tclusp.(k-beta) and cder = tclusp.(k) in
  let setpcommune = MolSet.inter (Clusp.mset cprem) (Clusp.mset cder) in
    MolSet.fold (fun p bool -> bool && (Clusp.debprot cprem p <= Clusp.debprot cder p)) setpcommune true


let gporderr k i tintmax tclusp beta = 
  let cprem = tclusp.(k-beta) and cder = tclusp.(k) in
  let setpcommune = MolSet.inter (Clusp.mset cprem) (Clusp.mset cder) in
    MolSet.fold (fun p bool -> bool && (Clusp.debprot cprem p >= Clusp.debprot cder p)) setpcommune true


let gporder k i tintmax tclusp beta =
  try
    (match (Clusp.brinc tclusp.(k)) with
      |Forward -> gporderf k i tintmax tclusp beta ;
      |Reverse -> gporderr k i tintmax tclusp beta)
  with
    |Invalid_argument _ -> true





(* pour la combinaison de sources et le blusping monomolécule protéique on est sur un seul brin à la fois
   d'ou l'utilisation de gporderf (et non gporder) *)
let gporderf2 k i tclusp beta =
  let cprem = tclusp.(k-beta) and cder = tclusp.(k) in
  let setpcommune = MolSet.inter (Clusp.mset cprem) (Clusp.mset cder) in
    MolSet.fold (fun p bool -> bool && ((Clusp.debprot cprem p < Clusp.debprot cder p) && (Clusp.finprot cprem p < 10+(Clusp.debprot cder p)))) setpcommune true


(* Dans le cas ou le blusping protéique a lieu avant le clusping il faut prévenir la mise
   dans le meme blusp d'hsps qui meme si dans le meme ordre gp ont un intron génomique trop grand 
   (seuil fixé dans le contexte) -> p de blusping2monomol sera gporderf && tintronok *)
let tintronok k i tintmax tclusp beta =
  let cprem = tclusp.(k-beta) and cder = tclusp.(k) in
  let setpcommune = MolSet.inter (Clusp.mset cprem) (Clusp.mset cder) in
    MolSet.fold (fun p bool -> bool && (((Clusp.debgen cder p)-(Clusp.fingen cprem p))<tintmax)) setpcommune true 


(* pb de l'hsp d'un blusp monoprot sans retour en arrière avec son précédant
   ou son suivant, mais petit et a une distance eloignée de son précédant ou 
   son suivant mais pas aussi grande que tmaxintron

   Ajouté en mars 2006
   Remarque : cette fonction et les deux précédantes pourraient s'appliquer sur un ensemble
   d'hsps plutot que de clusp car on est sur les blusp monoprot!!

   De plus beta est inutile, à éliminer partout
*)
let isolatednonok k i szisolh tintmax2 tclusp beta =
  let cprem = tclusp.(k-beta) and cder = tclusp.(k) in
  let boolcprem = ((Clusp.size cprem)<szisolh) and boolcder = ((Clusp.size cder)<szisolh) in
    (* hsp séparé du suivant de plus de seuil (ex : 10kb) 
       et celui-ci ou le suivant de taille < seuil (ex : 45 nt) 
    *)
    (not (tintronok k i tintmax2 tclusp beta)) && (boolcprem || boolcder)
    


(* Utilisé dans premain.ml pour faire les bluspmm protéiques *)
let gporder_tintok k i tintmax tclusp beta =
  (gporderf2 k i tclusp beta) && (tintronok k i tintmax tclusp beta)


(*  Ajouté en mars 2006 pour remplace la précédante fonction dans le blusping
    monoprot *)
let gporder_tintok_hspisolok k i szisolh tintmax tintmax2 tclusp beta =
  (gporder_tintok k i tintmax tclusp beta) && (not (isolatednonok k i szisolh tintmax2 tclusp beta))


(* Utilisé dans premain.ml pour faire les bluspmt protéiques *)
let noncutdwn k i tintmax tclusp beta =
  let cprem = tclusp.(k-beta) and cder = tclusp.(k) in
    not (Clusp.cutdwn cprem)

(* 
   lcoupe prend en entrée le tableau de clusps complet qui a deja été 
   passé à blusping1 (c'est pourquoi on part de la position k = i+beta), 
   une position de début et de fin de recherche (cette dernière non 
   incluse, rq : ce sont les positions des blusps de catégorie 1) 
   et le critere de laxisme beta. Elle renvoit une sous-liste 
   de positions sur tclusp, entre i et j, 
   eventuellement vide, qui represente l'ensemble des clusps k au niveau desquels il a été necessaire 
   de former un nouveau blusp du fait de la non-validation de la propriété p (sur tclusp, k, beta). 
   Pour p=gporder, cela équivaut à couper en la position k au sein du blusp i-j, ssi le clusp k et 
   le clusp k-beta ont une protéine commune, mais ont leurs hsps correspondant à cette protéine 
   ordonnées de façon opposée entre cette protéine et le génomique 
   Rq : à appeler avec p = gporder OU AUTRE fonction prenant obligatoirement en argument l'indice 
   i de début du blusp courant, le tableau de l'ensemble des clusps, des positions i et j de début 
   et de fin du blusp de catégorie au moins 1 (cad que tous les clusps d'un meme blusp ont au moins 
   UNE protéine commune) à analyser pour une éventuelle découpe, et un entier beta mesurant la 
   distance de découpe minimale (laxisme) 
*)
let lcoupe p tclusp tintmax i j beta = 
  if(j=i) then  (* cas où le dernier blusp est monocluspique *)
    []
  else
    begin
      let k = ref (i+beta) and k0 = ref (i+beta) and stop = ref false and l = ref [] in
	while(!k<>j) do
	  stop:= not (p !k !k0 tintmax tclusp beta);
	  if (!stop) then
	    begin
	      l:=(!k)::(!l); k0:=!k;
	    end;
	  incr k;
	done;
    List.rev !l;
    end



(* idem que lcoupe mais utilisé avec blusping2monomol2 *)
let lcoupe2 p tclusp szisolh tintmax tintmax2 i j beta = 
  if(j=i) then  (* cas où le dernier blusp est monocluspique *)
    []
  else
    begin
      let k = ref (i+beta) and k0 = ref (i+beta) and stop = ref false and l = ref [] in
	while(!k<>j) do
	  stop:= not (p !k !k0 szisolh tintmax tintmax2 tclusp beta);
	  if (!stop) then
	    begin
	      l:=(!k)::(!l); k0:=!k;
	    end;
	  incr k;
	done;
    List.rev !l;
    end



(* blusping2 prend le résultat de blusping1 (une segseq dont tous les clusps d'un blusp 
   ont au moins une protéine commune avec un laxisme de beta clusps) et renvoit une nouvelle 
   segseq, dont la segmentation a changé, cad qu'une propriete p a ete verifiee. 
   Si p=gporder alors la propriété est la suivante : "si ci et cbeta+i-1 ont une protéine commune, 
   alors le début sur cette protéine de l'hsp du clusp cbeta+i-1 matchant cette protéine, 
   doit etre plus grand que la fin sur cette protéine de l'hsp du clusp ci matchant cette protéine". 
   blusping2 forme des blusps de catégorie au moins 2, cad raffinés par rapport aux blusps de 
   catégorie 1 (si seul blusping1 a été appliqué sur s) ou meme par rapport à des blusps de 
   catégorie 2 (si blusping2(blusping1) p a déjà été appliqué sur s avec un p différent du p 
   actuel). par exemple on réalisera un deuxième blusping2 pour vérifier qu'au sein d'un 
   meme blusp, tous les clusps ont bien le meme brin 
   
   Attention : dans cette version et bcp des précédantes (des que ARNm en + de protéines)
   blusping2 n'est plus utilisée mais est remplacée par blusping2monomol
*)
let blusping2 p s tintmax beta =   
  let li = SegSeq.intervalles s in
    SegSeq.setsegment s (tri_fusion Pervasives.compare 
			   (SegSeq.ord s) 
			   (List.flatten 
			      (List.map 
				 (fun inter -> lcoupe p (SegSeq.tank s) tintmax (fst inter) (snd inter) beta) li)))
  



(* appelé dans premain.ml avec pour gporder_tintok 
   prend en entrée une segseq s et renvoit la meme segseq s mais dont
   le champ ord a été révu en fonction de la découpe que l'on souhaite faire
   si la propriété p n'est pas vérifiée entre deux éléments de s
   p prend cinq arguments et renvoit un booléen :
   - entier k 
   - entier k0 (inutile le plus souvent, sauf pour idbrin)
   - entier tintmax (utile que pour tintronok)
   - tableau de clusps tclusp (tjs utile)
   - entier beta 
*)
let blusping2monomol p s tintmax beta = 
  try
    (let li = [(0,Array.length (SegSeq.tank s))] in
       SegSeq.setsegment 
	 s  (* la segseq *) 
	 (tri_fusion 
	    Pervasives.compare 
	    (SegSeq.ord s) 
	    (List.flatten 
	       (List.map 
		  (fun inter -> lcoupe p (SegSeq.tank s)  tintmax (fst inter) (snd inter) beta)
		  li
	       )
	    )
	 )  (* la liste des positions des éléments ou on segmente la segseq s *)
    )
  with Invalid_argument _ -> Common.print_log "Pb dans blusping2monomol\n"; s
(* à terme on regroupera ces deux blusping2 en un en mettant un argument suppl, genre monomol à true ou false,
   et en fonction de cela li change (car c'est la seule chose à modifier).
   On s'est rendu compte que pour un blusp monomolecule il fallait pouvoir considérer le dernier 
   à moins qu'il faille changer Collection.SegSeq.intervalles en faisant aller jusqu'a length tank au
   lieu de length tank -1 ??*)




(* appelé dans premain.ml avec pour gporder_tintok_nonhspisol
   prend en entrée une segseq s et renvoit la meme segseq s mais dont
   le champ ord a été révu en fonction de la découpe que l'on souhaite faire
   si la propriété p n'est pas vérifiée entre deux éléments de s
   p prend sept arguments  :
   - entier k 
   - entier k0 (inutile le plus souvent, sauf pour idbrin)
   - entier szisolh
   - entier tintmax (utile que pour tintronok)
   - entier tintmax2 
   - tableau de clusps tclusp (tjs utile)
   - entier beta 
   
   et renvoit un booléen
*)
let blusping2monomol2 p s szisolh tintmax tintmax2 beta = 
  try
    (let li = [(0,Array.length (SegSeq.tank s))] in
       SegSeq.setsegment 
	 s  (* la segseq *) 
	 (tri_fusion 
	    Pervasives.compare 
	    (SegSeq.ord s) 
	    (List.flatten 
	       (List.map 
		  (fun inter -> lcoupe2 p (SegSeq.tank s) szisolh  tintmax tintmax2 (fst inter) (snd inter) beta)
		  li
	       )
	    )
	 )  (* la liste des positions des éléments ou on segmente la segseq s *)
    )
  with Invalid_argument _ -> Common.print_log "Pb dans blusping2monomol\n"; s


(* b est un blusp, donc un tableau de clusps *)
let bluspbrin b = Clusp.brinc b.(0)


let bluspdeb b = 
  Array.sort (fun c1 c2 -> Pervasives.compare (Clusp.gbeg c1) (Clusp.gbeg c2)) b;
  Clusp.gbeg b.(0)


let bluspfin b =
 Array.sort (fun c1 c2 -> Pervasives.compare (Clusp.gend c2) (Clusp.gend c1)) b;
 Clusp.gend b.(0)


(************************************************************************************************************
                   Fonctions de recherche de la protéine principale d'un blusp 
****************************************************************************************************************)



(* select_p prend en entrée une table de hachage h de clé une protéine et de contenu le nombre
   total d'acides aminés de cette protéine qui ont matché sur le génomique au sein du blusp courant.
   Cette fonction renvoit la protéine qui a le plus grans total d'acides aminés. *)
let select_p h =
  let pmax = ref (Molecule.null) and aamax = ref 0 in
    Hashtbl.iter (fun p aa -> if (aa > (!aamax)) then 
		    begin
		      aamax:=aa;
		      pmax:=p;
		    end) h;
    (!pmax,!aamax)


(* debfin prend en entree une segseq et un entier i et recherche parmi tous les intervalles
   (j,k) de sc ceux tels que j<=i<k *)
let debfin sc i =
  let li = SegSeq.intervalles sc and nc = Array.length (SegSeq.tank sc) in
  let nli = List.length li in
    if (i<(nc-1)) then
      List.hd (List.filter (fun (j,k) -> (j<=i && i<k)) li)
    else
      List.nth li (nli-1)



(* prend un couple de positions deb, fin sur une protéine et renvoit la différence en acides aminés *)
let diff_aa couple =
  (snd couple) - (fst couple) 



(* prend en entree un clusp et renvoit une table de hachage qui à chaque protéine du clusp associe
   le nombre d'acides aminés de cette protéine qui matchent sur le génomique S à annoter *)
let assoc_prot_nbaa c = 
  let h = Hashtbl.create (Clusp.nbhsp c) in
    Hashtbl.iter (fun prot y -> Hashtbl.add h prot (diff_aa (Hashtbl.find (Clusp.ppos c) prot))) (Clusp.ppos c);
    h


(* calcule la somme des elements entiers d'une liste *)
let rec lsum l =
  match l with 
    |[] -> 0;
    |t::q -> t+(lsum q)


let tsum t = lsum (Array.to_list t)


(* la protéine principale d'un clusp une fois connue la segmentation de la séquence en blusps et que i est le no du clusp c (en fait c = (SegSeq.tank sc).(i)), est celle qui a le plus d'acides aminés (qui matchent) au sein du Blusp contenant le clusp c de position i dans sc *)
let prot_pal sc i =
  try
    (let (ideb,ifin) = debfin sc i in  (* position de début et de fin du blusp courant (où est i) *)
       if(ideb=ifin) then
	 begin
	   let tclusp_ds_blusp = Array.sub (SegSeq.tank sc) ideb 1 in (* tableau des clusps du blusp *)
	   let ensprot_ds_blusp = ensprot tclusp_ds_blusp 0 1 in (* ensemble des protéines du blusp *)
	   let thassoc = Array.map assoc_prot_nbaa tclusp_ds_blusp in   (* tableau de tables de hachage, une par clusp
								  du blusp *)
	   let nbhsp_ds_blusp = tsum (Array.map (fun c -> Clusp.nbhsp c) tclusp_ds_blusp) in
	   let h_diff_aa = Hashtbl.create nbhsp_ds_blusp in
	     MolSet.iter (fun p -> Hashtbl.add h_diff_aa p (tsum (Array.map (fun h -> if (Hashtbl.mem h p) then Hashtbl.find h p else 0) thassoc))) ensprot_ds_blusp;
	     select_p h_diff_aa;
	 end
       else
	 begin
	   let tclusp_ds_blusp = Array.sub (SegSeq.tank sc) ideb (ifin-ideb) in (* tableau des clusps du blusp *)
	   let ensprot_ds_blusp = ensprot tclusp_ds_blusp 0 (ifin-ideb) in (* ensemble des protéines du blusp *)
	   let thassoc = Array.map assoc_prot_nbaa tclusp_ds_blusp in   (* tableau de tables de hachage, une par clusp
									   du blusp *)
	   let nbhsp_ds_blusp = tsum (Array.map (fun c -> Clusp.nbhsp c) tclusp_ds_blusp) in
	   let h_diff_aa = Hashtbl.create nbhsp_ds_blusp in
	     MolSet.iter (fun m -> Hashtbl.add h_diff_aa m (tsum (Array.map (fun h -> if (Hashtbl.mem h m) then Hashtbl.find h m else 0) thassoc))) ensprot_ds_blusp;
	     select_p h_diff_aa;
	 end)
  with
      (* si on ne trouve pas de protéine principale on renvoit la molécule nulle N *)
      Not_found -> (Donnees_base.N,0)
   



(* Fonction de recherche de la protéine ppale d'un blusp proprement dite.
   La meme que la précédente mais en prenant un tableau de clusps au lieu d'une segseq de clusps. 
   Attention : renvoit un couple constitué de la fameuse protéine, mais aussi de son nb aa sur le blusp 
   Sert dans l'élimination des blusps artéfactuels monocluspiques et dicluspiques pour les
   protéines seulement pour le moment *)
let prot_pal2 tcl =
  try
    (let n = Array.length tcl in
     let ensprot_ds_blusp = ensprot tcl 0 n in 
       (* ensemble des protéines du blusp *)
     let thassoc = Array.map assoc_prot_nbaa tcl in   
       (* tableau de tables de hachage, une par clusp du blusp, qui repertorie le nb d'aa matché 
	  par chaque hsp de chaque protéine du clusp *)
     let nbhsp_ds_blusp = tsum (Array.map (fun c -> Clusp.nbhsp c) tcl) in
       (* nombre d'hsps dans le blusp *)
     let h_diff_aa = Hashtbl.create (List.length (MolSet.elements ensprot_ds_blusp)) in
       (* table de hachage associant à chaque protéine du blusp le nb d'aa matchés par les différentes hsps
	  constituant les différents clusps du blusp *)
       MolSet.iter (fun p -> Hashtbl.add h_diff_aa p (tsum (Array.map (fun h -> if (Hashtbl.mem h p) then Hashtbl.find h p else 0) thassoc))) ensprot_ds_blusp;
       (* remplissage de la table de hachage puis selection de la protéine principale *)
       select_p h_diff_aa)
  with Not_found -> (Donnees_base.N,0)



(* 
   longp renvoit la longueur de la protéine prot dont la séquence est stockée dans la table de 
   hachage bq (de clé le nom de la protéine et de contenu la séquence protéique en question 
   Attention : renvoit l'exception Not_found si la protéine n'est pas dans la banque 
*)	
let longp prot bq = 
  match prot with
    |P s ->
       (try
	 Seq.length (Hashtbl.find bq prot)
       with
	   Not_found -> Common.print_log ("Cette protéine n'est pas dans la banque et c'est "^s^"\n"); 0)
    |_ -> 0



(* sélectionne l'hsp de la liste lhsp qui a pour protéine prot 
   Si aucune hsp de lhsp n'a comme protéine prot alors selectionne la premier hsp de lhsp *)
let select_hsp prot lhsp =
  let f = fun h -> (Hsp.nomp h = prot) in
  let lhprot = List.filter f lhsp in
    if (lhprot <> []) then
      List.hd lhprot
    else 
      List.hd lhsp





(*****************************************************************************************************************
                       Fonctions d'élimination des blusps artéfactuels selon une propriété p 
******************************************************************************************************************)





(* ici on fixe le seuil d'acceptation d'un blusp monocluspique : c'est en fonction de la quantité de la proteine de son hsp ppale qui est matché (>0 == pas de seuil, > lgp/10 ==> plus que 1 dixieme de la prot ...*)
let nonmonoeclfp propmol bqp tidxp tlgp ecl = 
  let labecl = Eclusp.lbl ecl and hsppecl = Eclusp.hspp ecl in
  let lgp = tlgp.(trouve_index tidxp (Hsp.nomp hsppecl)) in
    ((labecl <> Unique) || ((labecl = Unique) && (((Hsp.finp hsppecl)-(Hsp.debp hsppecl)+1) > lgp/3))) (* lgp/20 pour 5% lgp/10 pour 10%, lgp/2 pour 50%, lgp pour 100%*)
  



let pbcutdwnintcl tcl =
  List.exists (fun c -> (Clusp.cutdwn c)) (Array.to_list tcl)


(* nonmonoclart est une fonction qui vérifie si un blusp n'est pas monocluspique artéfactuel,
   cad s'il n'est pas monocluspique, ou s'il l'est que son hsp ppale matche plus d'un certain 
   pourcentage de sa protéine ppale. bqp est la banque de protéine sous forme de table de hachage, 
   tidxp est le tableau d'index des protéines, tlgp est le tableau de longueur des protéines et 
   tcl est un tableau de clusps, autrement dit le blusp à analyser. 
*(propmol/100)*)
let nonmonoclart propmol bqp tidxp tlgp tcl =
  try
  (
    let n = Array.length tcl and (pp,aa) = prot_pal2 tcl in   
      (* ajout en mars 2006 : s'il montre un pb de coupure on ne veut pas l'éliminer
	 car sinon on perd l'information de coupure
      *)
      if (n>=2 || (pbcutdwnintcl tcl)) then
	true
      else
	begin
	  let hsppcl = (select_hsp pp (Clusp.lhspp tcl.(0))) in
	  let lgp = tlgp.(trouve_index tidxp (Hsp.nomp hsppcl)) in
	    (* lgp/20 pour 5% lgp/10 pour 10%, lgp/2 pour 50%, lgp pour 100%
	       rajouter un ou si le champ ya til une hsp de meme proteine à proximité est à 0 *)
	    (((Hsp.finp hsppcl)-(Hsp.debp hsppcl)+1) > ((lgp*propmol)/100)) && (((Clusp.gend tcl.(0)) - (Clusp.gbeg tcl.(0)) +1) mod 3 =0);
	    (* on peut rajouter la clause de ne pas contenir de stop mais il faut alors regarder le brin du clusp
	       et en plus avoir à disposition la séquence génomique ici 
	       Remarque importante : pour bien faire il faudrait qd meme annoter les gènes monoex > seuil mais
	       non mod 3 car ce sont en gal des pseudogènes et notre annotation sera ainsi plus complete *)
	end 
  )
  with
      Invalid_argument _ -> Common.print_log "Pb dans nonmonoclart\n"; true



 
(* nondiclart est une fonction qui vérifie si un blusp n'est pas dicluspique artéfactuel,
   cad s'il n'est pas dicluspique, ou s'il l'est que son hsp ppale matche plus d'un certain 
   pourcentage de sa protéine ppale. bqp est la banque de protéine sous forme de table de hachage, 
   tidxp est le tableau d'index des protéines, tlgp est le tableau de longueur des protéines et 
   tcl est un tableau de clusps, autrement dit le blusp à analyser. 
*)
let nondiclart propmol bqp tidxp tlgp tcl = 
  try
    (
      let n = Array.length tcl and (pp,aa) = prot_pal2 tcl in
	(* ajout en mars 2006 : s'il montre un pb de coupure on ne veut pas l'éliminer
	   car sinon on perd l'information de coupure
	*)
	((n!=2) || (pbcutdwnintcl tcl) || (aa > ((tlgp.(trouve_index tidxp pp))*propmol)/100))
    )
  with
      Invalid_argument _ -> Common.print_log "Pb dans nondiclart\n"; true





(* calcule le nombre d'acides aminés manquants dans la protéine d'un blusp monoprot
   représenté par le tableau de clusp monohsp tcl (en absolu mais théoriquemement 
   ici ce nombre doit être positif à cause du découpage précédant)
   Attention : pour tcl de taille 1 on aura 0
   Attention : on suppose que tcl est un blusp monoprotéique = 1 seule protéine

   Ajouté en mars 2006 pour pallier au pb des blusp monoprot bien formés mais
   en dessous du seuil propmol -> on souhaite pour eux baisser ce seuil à 20%
   (seuil que l'on sortira plus tard en paramètre) *)
let nbaamanquant bqp tidxp tlgp tcl =
  let n = Array.length tcl in
    if (n<=1) then
      0
    else
      begin
	let prot= MolSet.choose (Clusp.mset tcl.(0)) and i = ref 0 and nbaa = ref 0 in
	  while(!i<n-1) do
	    nbaa:=(!nbaa)+ abs ((fst (Hashtbl.find (Clusp.ppos tcl.(!i+1)) prot)) - (snd (Hashtbl.find (Clusp.ppos tcl.(!i)) prot)) -1);
	    incr i;
	  done;
	  !nbaa;
      end
    


(* calcule le nombre de nucléotides présents et le nombre de nucléotides fusionnés
   dans un ensemble d'hsps contenu dans un thsp
*)
let nbnt_ntfus_inthsp thsp =
  let i = ref 0 and n = Array.length thsp and nbnt = ref 0 and nbntfus = ref 0 in
    while(!i<n) do
      nbnt := (!nbnt) + (Hsp.fing thsp.(!i)) - (Hsp.debg thsp.(!i)) + 1 ;
      nbntfus := (!nbntfus) + (Hsp.nbntfus thsp.(!i));
      incr i;
    done;
    (!nbnt,!nbntfus);;


(* fonction qui dit si un blusp mono protéique n'est pas artéfactuel
   pour cela il suffit :
   - soit qu'il ait répertorié un pb de cut down
   (car on ne veut pas perdre cette information qui nous permettra de casser
   un blusp multiprot en plusieurs plus loin),
   -soit que son nombre d'acides aminés fasse - de propmol% de la protéine 

   Attention : fonction qui a remplacé nonmonoclart et nondicluart en Mars 2006 (pb Cluster)
   Attention : propmol2 est de 20% et sert pour les blusps monoprot très bien formés (= peu d'acides 
   aminés manquants dans la protéine) 
*)
let nonmonoprotart pnntfus pnaamq propmol propmol2 bqp tidxp tlgp tcl = 
  try
    (
      (* pb si le clusp contient ici un nombre d'hsps différent de 1 
	 0 -> crash de List.hd (Clusp.lhsp c)
	 >=2 -> la fusion n'a pas bien fait son travail? *)
      let n = Array.length tcl and (pp,aa) = prot_pal2 tcl and nbaamq = nbaamanquant bqp tidxp tlgp tcl and thsp = Array.map (fun c -> List.hd (Clusp.lhsp c)) tcl in

      let (nbnt,nbntfus) = nbnt_ntfus_inthsp thsp in

	(* ajout en mars 2006 : s'il montre un pb de coupure on ne veut pas l'éliminer
	   et ce quel que soit son propmol car sinon on perd l'information de coupure 
	   utile au moment du DACM2 protéique et qui permet à des protéines de s'entre-aider
	   au niveau de l'information de coupure (si au moins un hsp prot a une coupure
	   alors on coupe à ce niveau le blusp multiprot)
	   
	   De plus on veut a la fois :
	   - tolérer les blusp monoprot de 3 hsps au moins auxquels il manque peu d'aa 
	   % la protéine comparée (on demande pour eux un seuil propmol2 < propmol)
	   - eliminer les blusp monoprot qui passent le seuil propmol à cause de bcp de 
	   fusions dans leurs hsps.

	   Plus formellement pour que ce soit un pseudogène il faut que :
	   - pas de pb de cut
                     ET
           - [(mono ou diex) ET (propaa < seuil1)]
                     OU
             [(au - triex) ET [
                               (bcp aa manquants) 
                                   OU 
                               (propaa < seuil2)
                                   OU
                                (bcp nt fusionnés)
                               ]
             ]

	*)
	(
	  (pbcutdwnintcl tcl) 	

	  || 
	  
	  ( ((n>=3)  || (aa > ((tlgp.(trouve_index tidxp pp))*propmol)/100)   )
          
	  &&

	    ((n<2) ||  ((((nbaamq*100)/aa)<=pnaamq) 
                        && 
                        (aa > ((tlgp.(trouve_index tidxp pp))*propmol2)/100) 
                        &&
                        (((nbntfus*100)/nbnt)<=pnntfus)))  
          )
	)

    )
	 
  with
      Invalid_argument _ -> Common.print_log "Pb dans nonmonoprotnonart\n"; true


(*
qui marchait bien :
	(
	  (pbcutdwnintcl tcl) 	

	  || ((n>=3) 
	      && (((nbaamq*100)/aa)<=pnaamq) 
	      && (aa > ((tlgp.(trouve_index tidxp pp))*propmol2)/100) 
	      && (((nbntfus*100)/nbnt)<=pnntfus))

          || ((aa > ((tlgp.(trouve_index tidxp pp))*propmol)/100)) 
	      (* && (((nbntfus*100)/nbnt)<=pnntfus)) *)
	)
*)
