
(* sigsearch.ml : comme son nom l'indique, regroupe les fonctions de recherche de signaux à partir des blusps,
   et regroupe donc les phases de passage de clusp aux eclusps, de eclups aux exons et des exons aux gènes. *)


open Common
open Seq
open TraiteSignal
open Preserveorf
open Donnees_base
open Donnees
open Collection
open Matrice
open TraiteCollection
open Config 



(* clabel prend en entree une segseq s de clusps (entree et sortie du blusping) et un entier i, et attribue au ieme clusp (en partant de 0) son label deduit du blusping, cad clusp unique, initial, interne ou terminal au sein du blusp *)
let clabel sc i =
  if List.mem i (SegSeq.ord sc) then (* donc initial *)
    begin
      if ((List.mem (i+1) (SegSeq.ord sc))||(i= ((Array.length (SegSeq.tank sc)) -1))) then
	Unique (* car initial ET terminal *)
      else
	Initial
    end
  else
    if ((List.mem (i+1) (SegSeq.ord sc)) || (i= ((Array.length (SegSeq.tank sc)) -1))) then 
      Terminal
    else
      Internal





(* prot2hsp prend une protéine p et un clusp c et renvoit un couple dont le premier
   element dit si la protéine p est représentée dans une hsp du clusp et dont le
   deuxieme element donne, dans le cas ou c'est la cas l'hsp de c de protéine p 
   Ressemble à select_hsp qui prend une proteine et une liste d'hsps mais pas exactement
   le meme type de résultat (select_hsp moins propre, il faudra fair eun merge avec celle-ci
   mais attention : repercution à bcp d'endroits) *)
let prot2hsp p c = 
  let lh = Clusp.lhsp c in
  let lhp = List.filter (fun h -> ((Hsp.nomp h)=p)) lh in
    if(lhp=[]) then
      (false,Hsp.null)
    else
      (true,List.hd lhp)

   


(* polytypec2protc transforme un clusp d'hsps multi type (prot et arn) en clusp avec 
   que des hsps protéiques. Le but est de pouvoir calculer la protéine principale
   sans etre perturbé par les arn *)
let polytypec2protc c =
  let lh = Clusp.lhsp c in
  let lhprot = List.filter (fun h -> ((Hsp.tmol h) = Prot)) lh in
    thspchev_to_clusp (Array.of_list lhprot)



(* polytypec2arnc transforme un clusp d'hsps multi type (prot et arn) en clusp avec 
   que des hsps arn. Juste là par symetrie, au cas ou!  *)
let polytypec2arnc c =
  let lh = Clusp.lhsp c in
  let lharn = List.filter (fun h -> ((Hsp.tmol h) = Rna)) lh in
    thspchev_to_clusp (Array.of_list lharn)



let filterprot lh = List.filter (fun h -> ((Hsp.tmol h) = Prot)) lh

let filterarn lh = List.filter (fun h -> ((Hsp.tmol h) = Rna)) lh


let existhspprot lh =
 (filterprot lh) != []

let existhsparn lh = 
   (filterarn lh) !=[]



let c2hspofmol c mol =
  let lh = Clusp.lhsp c in
    List.hd (List.filter (fun h -> ((Hsp.nomp h) = mol)) lh)



(* regarde s'il existe un arn commun aux clusps c et cnext et renvoit un couple bool*Hsp.t
   ou le premier indique si cet arn existe et le deuxieme renvoit l'hsp de l'un de ces arns communs *)
let existcommonarn c cnext =
  let molsetc = MolSet.filter (fun m -> match m with (A _) -> true | _ -> false) (Clusp.mset c) and molsetcnext = MolSet.filter (fun m -> match m with (A _) -> true | _ -> false) (Clusp.mset cnext) in
  let interset = (MolSet.inter molsetc molsetcnext) in
    if (MolSet.is_empty interset) then
      (false,Hsp.null)
    else
      (true,c2hspofmol c (MolSet.choose interset))
      




(* 
   nbaau renvoit le nombre d'acides aminés qu'il faut chercher en amont (up) du clusp c no i de 
   brin forward de la segseq sc de clusps. Rq : cette information est intrinsèque à c (ses protéines) 
   lorsque c est Unique ou Initial, mais est apportée par le clusp précédent si c est Internal ou Terminal) 
   Ajout au 17/03/06 :
   si hsp la + 5' est proteique et qu'il ne manque pas d'aa dans cette prot et que l'atg est là
   alors on le prend, implication : gseq en param de nbaau et pour homogénéiser en param de nbaad
*)
let nbaau c sc i p gseq bqp nbdef =
  let lhuc = (Clusp.lhspu c) and hspp = select_hsp p (Clusp.lhspp c) and (bp,hp) = prot2hsp p c in
  let hspu = select_hsp p lhuc and lhc = Clusp.lhsp c in
    match (clabel sc i) with
      | (Unique | Initial) -> 
	  if (existhsparn lhuc) then
	    {nbaa = 0; mol = true; hsp = List.hd (filterarn lhuc); prot = false}
	  else
	    begin
	      if (((Hsp.debp hspu)=0) && ((Seq.sub gseq (Hsp.debg hspu) 3))=(List.hd met)) then
		{nbaa = 0; mol= true; hsp = hspu; prot = true}
	      else 
		{nbaa = (Hsp.debp hspu)+40; mol= true; hsp = hspu; prot = true}
		  (* +40 ajouté par rapport à version MMuffato, à terme à sortir *)
	    end
      | (Internal | Terminal) -> 
	  let cprec = (SegSeq.tank sc).(i-1) in
	  let lhcprec = Clusp.lhsp cprec in
	  let hspprec = select_hsp p (Clusp.lhspp cprec) in
	  let (barncom,hsparn)= existcommonarn c cprec in
	    (* si c et cprec ont un arn en commun on se fonde sur lui *)
	    if barncom then
	      {nbaa = 0; mol = true; hsp = hsparn; prot = false}
	    else
	      begin
		let protinccprec = (existhspprot lhc) && (existhspprot lhcprec) in
		  (* si c et cprec ont une proteine de chaque coté -> comme pour exogprot *)
		  if protinccprec then
		    begin
		      let (bpprec,hpprec) = prot2hsp p cprec in
			if((bp && bpprec) && ((Hsp.debp hp)>=(Hsp.finp hpprec))) then
			  {nbaa = (Hsp.debp hp) - (Hsp.finp hpprec) -1; mol = true; hsp = hp; prot = true}
			else
			  begin
			    if (((Hsp.nomp hspp) = (Hsp.nomp hspprec))&&((Hsp.debp hspp)>=(Hsp.finp hspprec))) then
			      {nbaa = (Hsp.debp hspp) - (Hsp.finp hspprec) -1; mol = true; hsp = hspp; prot=true}
			    else 
			      {nbaa = nbdef; mol = false; hsp = hspp; prot=true}
			  end 
		    end
		  else
		    (* si c et cprec n'ont pas de proteine de chaque coté -> comme pour exogprot
		       mais sans regarder la phase et le stop *)
		    {nbaa = nbdef; mol = false; hsp = hspp; prot=false}
	      end


		  

	


(* 
   nbaad renvoit le nombre d'acides aminés qu'il faut chercher en aval (down) du clusp c no i 
   de brin forward de la segseq sc de clusps (information apportée par c et par la longueur de 
   la protéine de l'hsp droite pour Unique et Terminal, et par le clusp suivant pour Initial et Internal) 
   Attention hp peut etre l'hsp nulle (Hsp.null), dans le cas ou p n'apparait pas dans c
*)
let nbaad c sc i p gseq bqp nbdef =
  let lhdc = (Clusp.lhspd c) and hspp = select_hsp p (Clusp.lhspp c) and (bp,hp) = prot2hsp p c in
  let hspd = select_hsp p lhdc and lhc = Clusp.lhsp c in
    match (clabel sc i) with
      | (Unique | Terminal) ->
	  if (existhsparn lhdc) then
	    {nbaa = 0; mol = true; hsp = List.hd (filterarn lhdc); prot = false}
	  else
	    (* Ici il n'existe que des protéines dans les hsps les plus à droite
	       du clusp courant qui est unique ou terminal
	       Attention : ici longp peut rendre Not_found (si Hsp.nomp hspd est P "") *)
	    begin
	      (* let u = match (p,(Hsp.nomp hspd)) with
		|(P sp,P sphspd) -> 
		   begin
		     if ((sp="")||(sphspd="")) then
		       begin
			 Common.print_log ("p est "^sp^" et (Hsp.nomp hspd) est "^sphspd^"\n"); 
			 Donnees.print_clusp stderr c
		       end
		   end
		|_ -> 
		   begin 
		     Common.print_log ("Attention p ou (Hsp.nomp hspd) n'est même pas considérée comme une protéine\n"); 
		     Donnees.print_clusp stderr c 
		   end 
	      in *)
		{nbaa = (longp (Hsp.nomp hspd) bqp) - (Hsp.finp hspd)+40; mol= true; hsp = hspd; prot = true}
		  (* +40 ajouté par rapport à version Mmuffato, à terme à sortir *)
	    end
      | (Initial | Internal) -> 
	  let csuiv = (SegSeq.tank sc).(i+1) in
	  let lhcsuiv = Clusp.lhsp csuiv in
	  let hsppsuiv = select_hsp p (Clusp.lhspp csuiv) in
	  let (barncom,hsparn)= existcommonarn c csuiv in
	    (* si c et csuiv ont un arn en commun alors on se fonde sur lui *)
	    if barncom then
	      {nbaa = 0; mol = true; hsp = hsparn; prot = false}
	    else
	      begin
		let protinccsuiv = (existhspprot lhc) && (existhspprot lhcsuiv) in
		  (* si c et csuiv ont une proteine de chaque coté -> comme pour exogprot *)
		  if protinccsuiv then
		    begin
		      let (bpsuiv,hpsuiv) = prot2hsp p csuiv in
			if((bp && bpsuiv) && ((Hsp.debp hpsuiv)>=(Hsp.finp hp))) then
			  {nbaa = (Hsp.debp hpsuiv) - (Hsp.finp hp) -1; mol = true; hsp = hp; prot = true}
			else
			  begin
			    if (((Hsp.nomp hspp) = (Hsp.nomp hsppsuiv))&&((Hsp.debp hsppsuiv)>=(Hsp.finp hspp))) then
			      {nbaa = (Hsp.debp hsppsuiv) - (Hsp.finp hspp) -1; mol = true; hsp = hspp; prot=true}
			    else 
			      {nbaa = nbdef; mol = false; hsp = hspp; prot=true}
			  end 
		    end
		  else
		    (* si c et csuiv n'ont pas de proteine de chaque coté -> comme pour exogprot
		       mais sans regarder la phase et l'absence de stop *)
		    {nbaa = nbdef; mol = false; hsp = hspp; prot=false}
	      end





      
open Alphaadn





(****************        UP     *********************)

(* lposigstopu retourne la liste des STOPs situés entre le début du clusp c et 
   et clusp-(nbaa*3+7) (en amont). Ces STOPs sont en phase avec le début du clusp c 
   avant : nninhsp = 8 , nnsuppl = 15 , nnsuppl-nninhsp = 7 *)
let lposigstopu nninhsp nnsuppl c sc i p gbseq bqp nbdef =
  let aa_goal = nbaau c sc i p gbseq bqp nbdef in
    match (clabel sc i) with
      | Unique -> if((Clusp.gbeg c)-1<0) then [] else (sig_search_up gbseq ((Clusp.gbeg c)-1) (aa_goal.nbaa*3+nnsuppl-nninhsp) 3 stop) (* ajout 15/03/06, cas du clusp nul *)
      | Initial -> if((Clusp.gbeg c)-1<0) then [] else (sig_search_up gbseq ((Clusp.gbeg c)-1) (aa_goal.nbaa*3+nnsuppl-nninhsp) 3 stop) (* ajout 15/03/06, cas du clusp nul *)
      | Internal -> sig_search_up gbseq ((Hsp.debg aa_goal.hsp) + nninhsp) (aa_goal.nbaa*3 + nnsuppl) 3 stop
      | Terminal -> sig_search_up gbseq ((Hsp.debg aa_goal.hsp) + nninhsp) (aa_goal.nbaa*3 + nnsuppl) 3 stop



(* lposigtotalu retourne la liste totale des signaux situés entre le début du clusp c
   et "clusp - (nbaa*3+7)". Ces signaux peuvent être ATG ou AG (selon le label du clusp c) *)
let lposigtotalu nninhsp nnsuppl c sc i p gbseq bqp nbdef =
  let {nbaa = naa; mol = bm; hsp = h2follow; prot = bp} = nbaau c sc i p gbseq bqp nbdef in
    if(not bp && bm) then
      begin
	match (clabel sc i) with
	  | (Unique | Initial) -> [Hsp.debg h2follow]
	  | (Internal | Terminal) -> [(Hsp.debg h2follow)-2] (* car là on veut la position du AG *)
      end
    else
      begin
	match (clabel sc i) with
	  | Unique -> if((Clusp.gbeg c)-1<0) then [] else (sig_search_up gbseq ((Clusp.gbeg c)-1) (naa*3+nnsuppl-nninhsp) 3 met) (*+5,+13*) (* ajout 15/03/06, cas du clusp nul *)
	  | Initial -> if((Clusp.gbeg c)-1<0) then [] else (sig_search_up gbseq ((Clusp.gbeg c)-1) (naa*3+nnsuppl-nninhsp) 3 met) (*+5,+13*)(* ajout 15/03/06, cas du clusp nul *)
	  | Internal -> sig_search_up gbseq ((Hsp.debg h2follow) + nninhsp) (naa*3 + nnsuppl) 1 accept
	  | Terminal -> sig_search_up gbseq ((Hsp.debg h2follow) + nninhsp) (naa*3 + nnsuppl) 1 accept
      end


(* 
   lposigu prend en argument un clusp de brin forward, la segseq de clusps sc associée, 
   le no i du clusp dans ce tableau (tank de segseq), la séquence à annoter gseq, 
   ainsi que la banque de protéine bqp qui sert à l'analyse et renvoit l'ensemble 
   des positions de son signal gauche (up) qui ne sont pas précédées d'un STOP en phase 
   avec le début du clusp.  
*)
let lposigu nninhsp nnsuppl c sc i p gbseq bqp nbdef = 
  try
    (let aa_goal = nbaau c sc i p gbseq bqp nbdef and listeTotSigu = lposigtotalu nninhsp nnsuppl c sc i p gbseq bqp nbdef in
       (* si on se fonde sur une protéine alors on regarde phase + stop  *)
       if (aa_goal.prot) then
	 begin
	   let l = ref [] in
	   let listeStopu = lposigstopu nninhsp nnsuppl c sc i p gbseq bqp nbdef in
	     match listeStopu with
	       |[] -> (aa_goal, listeTotSigu);
	       |premstop::q -> (aa_goal,List.filter (fun pos -> (pos>=premstop)) listeTotSigu)
	 end
	   (* sinon il est inutile de regarder l'ORF car inconnue encore *)
       else
	 (aa_goal, listeTotSigu))
  with 
      Not_found ->  Common.print_log "pb dans lposigu\n"; ({nbaa =0; mol = false; hsp = Hsp.null; prot = false},[])
	    
	 




(****************        DOWN  *********************)

let lposigstopd nninhsp nnsuppl c sc i p gbseq bqp nbdef =
  let aa_goal = nbaad c sc i p gbseq bqp nbdef in
    match (clabel sc i) with
      | Unique   -> sig_search_down gbseq ((Clusp.gend c)+1) (aa_goal.nbaa*3+nnsuppl-nninhsp) 3 stop;
      | Initial  -> sig_search_down gbseq ((Hsp.fing aa_goal.hsp) - nninhsp) (aa_goal.nbaa*3 + nnsuppl) 3 stop;
      | Internal -> sig_search_down gbseq ((Hsp.fing aa_goal.hsp) - nninhsp) (aa_goal.nbaa*3 + nnsuppl) 3 stop;
      | Terminal -> sig_search_down gbseq ((Clusp.gend c)+1) (aa_goal.nbaa*3+nnsuppl-nninhsp) 3 stop



let lposigtotald nninhsp nnsuppl c sc i p gbseq bqp nbdef =
  let {nbaa = naa; mol = bm; hsp = h2follow; prot = bp} = nbaad c sc i p gbseq bqp nbdef in
     if(not bp && bm) then
      begin
	match (clabel sc i) with
	  | (Unique | Terminal) -> [(Hsp.fing h2follow)+1] (* on prend ici la fin du transcrit *)
	  | (Initial | Internal) -> [(Hsp.fing h2follow)+1] (* car là on veut la position du G du GT *)
      end
    else
      begin
	match (clabel sc i) with
	  | Unique   -> sig_search_down gbseq ((Clusp.gend c)+1) (naa*3+nnsuppl-nninhsp) 3 stop
	  | Initial  -> sig_search_down gbseq ((Hsp.fing h2follow) - nninhsp) (naa*3 + nnsuppl) 1 donnor
	  | Internal -> sig_search_down gbseq ((Hsp.fing h2follow) - nninhsp) (naa*3 + nnsuppl) 1 donnor
	  | Terminal -> sig_search_down gbseq ((Clusp.gend c)+1) (naa*3+nnsuppl-nninhsp) 3 stop
      end



(* Un Not_found proviendrait de longp si la protéine n'est pas dans la banque *)
let lposigd nninhsp nnsuppl c sc i p gbseq bqp nbdef = 
  let aa_goal = (try (nbaad c sc i p gbseq bqp nbdef ) with Not_found -> Common.print_log "pb dans lposigd-nbaad\n"; {nbaa =0; mol = false; hsp = Hsp.null; prot = false}) and listeTotSigd = (try (lposigtotald nninhsp nnsuppl c sc i p gbseq bqp nbdef ) with Not_found -> Common.print_log "pb dans lposigd-lposigtotald\n"; []) in
    (* si on se fonde sur une protéine alors on regarde phase + stop  *)
    if (aa_goal.prot) then
      begin
	let l = ref [] in
	let listeStopd = (try (lposigstopd nninhsp nnsuppl c sc i p gbseq bqp nbdef) with Not_found -> Common.print_log "pb dans lposigd-lposigtotald\n"; []) in
	  match listeStopd with
	    |[] -> (aa_goal, listeTotSigd);
	    |premstop::q -> (aa_goal,List.filter (fun pos -> (pos<=(premstop+2))) listeTotSigd)
      end
    else
      (aa_goal, listeTotSigd)
(*ici je prend les signaux GT dont la position est inferieure ou égale à celle du (premier Stop +2)
exemple: TAGTAG : il ya un GT mais ignoré à cause du premier Stop TAG *)	







(***************************************************************************************
  VERIFICATION QU'UN STOP N'EST PAS ENGENDRE DANS L'AJOUT DE NT QUI FORME L'INTRON
  REMARQUE : VALABLE QUE SI L'ON SUIT UNE PROTEINE DONNEE (PAS POUR LES ARNS) 
*****************************************************************************************)


(* complémentaire d'une seq génomique *)       
let complement seq =
  let len =  Seq.length seq in 
  let compseq = Seq.create len Alphaadn.nulla in
  let ub  = len -1 in 
  let ub2 = (len/2)-1  in
  let tmp = ref Alphaadn.nulla  in 

    for i = 0 to ub2  do
      begin 
	tmp := Alphaadn.complement seq.{i};
	compseq.{i} <- Alphaadn.complement seq.{ub-i}; 
	compseq.{ub-i} <- !tmp 
      end
    done;
    if (len  mod 2 = 1) then 
      compseq.{1+ub2} <- Alphaadn.complement seq.{1+ub2};
    compseq





(* Remarque importante : l'analyse se fait brin séparé du début à la fin d'exogean
   gseq prendra donc la valeur soit de la séquence génomique brin +, 
   soit de laséquence génomique brin -
   De ce fait cela simplifie la recherche de signaux car on recherche toujours GT,AG 
   (et non leur complémentaire, cf ec_to_ex)

   La fonction concat2Seg ne sert que si la recherche de signaux se fait grâce aux protéines
   en effet il faut alors vérifier que ce que l'on ajoute est bien modulo 3 et sans stop
  
   La fonction concat2Seg sert étant donnés deux eclusps voisins ec1 et ec2,
   leurs positions ainsi que les positions gp1 et gp2 des signaux GT, AG trouvés
   à renvoyer la séquence concaténée correspondant aux deux exons

   gseq = séquence à analyser
   gb1,ge1 = debut,fin du eclusp 1
   gb2,ge2 = debut,fin du eclusp 2
   gp1 = position du signal GT dans le eclusp 1 (niveau du G dans GT)
   gp2 = position du signal AG dans le eclusp 2 (niveau du A dans AG)
   *) 
let concat2Seg gb1 ge1 gp1 gp2 gb2 ge2 gseq =
  (* si les exons obtenus à gauche (1) et à droite (2) comportent au moins 1 nucléotide 
     cad  si ((gp1>=gb1+1) && (gp2<=ge2-(1+1))) *)
  if ((gp1>=gb1+1) && (gp2<=ge2-(1+1))) then 
    (true,Seq.concat (Seq.sub gseq gb1 (gp1-gb1)) (Seq.sub gseq (gp2+(1+1)) (ge2-gp2-1)))
  else 
    (false,(Seq.create 0 Alphaadn.A))





(*  tremplit prend en entree 2 listes l1 et l2 contenant le meme type d'objet et renvoit un tableau 
    contenant tous les couples possibles avec les elements de l1 et l2 *)
let tremplit l1 l2 = 
  let lg1 = List.length l1 and lg2 = List.length l2 in
  let lgt = lg1 * lg2 in 
  let t = Array.create lgt (0,0) and i = ref 0 and j = ref 0 and k = ref 0 in 
    while (!i < lg1) do
      while (!j < lg2) do
	t.(!k)<- ((List.nth l1 !i),(List.nth l2 !j));
	incr k;
	incr j;
      done;
      j:=0;
      incr i;
    done;
    t





(* filter_frame prend un tableau t de couples de positions bordant un intron (cad de signaux GT/AG) 
   et renvoit un tableau de couple de positions d'introns telles que le rajout de nts est mod3
   et qu'il n' ya de stop ni aux extensions ni à la jonction entre les deux exons
   filter_frame est donc un filtre qui préserve la phase donnée par les protéines 
   Remarque : cette fonction n'est utile que si l'on suit une proteine pour former l'intron et donc
   si le champ prot de nbaa2search est à true *)
let filter_frame gb1 ge1 gb2 ge2 gseq t = 
  Array.of_list (List.filter (fun (gp1,gp2) -> (pasStop (concat2Seg gb1 ge1 gp1 gp2 gb2 ge2 gseq))) (List.filter (fun (gp1,gp2) -> (((gp1-ge1-1 + gb2-gp2-2) mod 3)=0)) (Array.to_list t)))
    (* trie selon la valeur absolue du rajout de nucléotides croissant *)
    (* filtre parmi tous ces couples ceux dont le rajout est modulo 3 et qui ne contiennent pas de stop *)
    

let filter_size_intron seuil t = 
  Array.of_list (List.filter (fun (gp1,gp2) -> (gp2-gp1 > seuil)) (Array.to_list t))



(* selectsig_minaddlost prend en entree le nombre d'acides aminé que l'on souhaiterait atteindre entre
   les deux clusps, un brin, la fin (en prenant toujours comme référence le début du génomique d'entree) 
   ge1 du clusp 1, le début (en prenant toujours comme référence le début du génomique d'entree) gb2 du clusp 2, 
   la liste l1 des positions des signaux qui suivent (en prenant toujours comme référence le début du génomique 
   d'entree) ge1 et la liste l2 des positions des signaux qui précèdent (en prenant toujours comme référence 
   le début du génomique d'entree) gb2.
   Cette fonction renvoie la liste de tous les couples (pos1,pos2) provenant resp de l1 et de l2, 
   rangés par ordre croissant de rajout de nucléotides (en valeur absolue) par rapport à ge1 et gb2, 
   qui vérifie un nombre de nucléotides à rajouter multiple de 3.
*)
let selectsig_minaddlost nbaa_goal tcplsigfiltered gb1 ge1 gb2 ge2 gseq =
  (* trie selon la valeur absolue de la difference entre le nombre d'acides aminés à rajouter et 
     le nombre d'acides aminés manquants (car on sait que c'est mod3 cette fois) croissant *)
  Array.stable_sort (fun (gp1,gp2) (gp3,gp4) -> Pervasives.compare (abs ((gp1-ge1-1 + gb2-gp2-2)/3-nbaa_goal)) 
		       (abs ((gp3-ge1-1 + gb2-gp4-2)/3-nbaa_goal))) tcplsigfiltered;
  
  Array.to_list tcplsigfiltered






(* selectsig_minadding prend en entree un brin, la fin (en prenant toujours comme référence le début du 
   génomique d'entree) 
   ge1 du clusp 1, le début (en prenant toujours comme référence le début du génomique d'entree) gb2 du clusp 2, 
   la liste l1 des positions des signaux qui suivent (en prenant toujours comme référence le début du génomique 
   d'entree) ge1 et la liste l2 des positions des signaux qui précèdent (en prenant toujours comme référence 
   le début du génomique d'entree) gb2.
   Cette fonction renvoie la liste de tous les couples (pos1,pos2) provenant resp de l1 et de l2, 
   rangés par ordre croissant de rajout de nucléotides (en valeur absolue) par rapport à ge1 et gb2, 
   qui vérifie un nombre de nucléotides à rajouter multiple de 3.
*)
let selectsig_minadding tcplsigfiltered gb1 ge1 gb2 ge2 gseq =
  Array.stable_sort (fun (gp1,gp2) (gp3,gp4) -> 
		       Pervasives.compare (abs (gb2-gp2-2 + gp1-ge1-1))  
		       (abs (gb2-gp4-2 + gp3-ge1-1))) tcplsigfiltered;
  Array.to_list tcplsigfiltered
    	  





(* compare2cpls sert à comparer deux couples de floats, en cumulant leur somme *)
let compare2cpls (a,b) (c,d) =
  Pervasives.compare (a+.b) (c+.d)



(* selectsig_sc_w_matrix : fonction principale de selction de signaux par matrice de poids.
   Prend en entrée un tableau de couples de positions d'intron; les trie en fonction du score
   cumulé des signaux gauche et droit et le renvoit sous forme de liste *)
let selectsig_sc_w_matrix tcplsigfiltered prob gseq mAG mGT =
  try
    (Array.sort (fun (a,b) (c,d) -> compare2cpls ((score mGT prob (Seq.sub gseq (c-4) 10)),(score mAG prob (Seq.sub gseq (d-4) 10))) ((score mGT prob (Seq.sub gseq (a-4) 10)),(score mAG prob (Seq.sub gseq (b-4) 10)))) tcplsigfiltered;
     Array.to_list tcplsigfiltered)
  with
      Invalid_argument s -> Common.print_log ("Pb dans selectsig_sc_w_matrix et "^s^"\n"); []



let selectsig_arn_simple tcpl = Array.to_list tcpl




(*selectsig : fonction principale de selection de signaux, en fonction de la façon dont on veut le faire, renvoie
  le couple de signaux d'epissage bordant l'intron qui a été sélectionné (utilisé dans ec_to_ex par la suite)
  cpl1 et cpl2 sont des couples dont le premier element est le nombre d'acides aminés manquants dans la 
  prot ppale, 
  et le deuxieme est la liste des positions resp droite du clusp i et gauche du clusp i+1 *)
let selectsig n cpl1 cpl2 seuil_int mss gb1 ge1 gb2 ge2 prob gseq mAG mGT =
  (* quelle que soit la methode de selection des signaux (Aamiss ou Wghtmat), on commence toujours par appliquer 
     la fonction filter_frame, ceci afin de ne conserver que les rajouts de nts mod3 et de filtrer les stops *)
  let {nbaa = naa; mol = bm; hsp = h2follow; prot = bp} = fst cpl1 and tcpl = tremplit (snd cpl1) (snd cpl2) in
  let tcpl_int_ok = filter_size_intron seuil_int tcpl in
    (* si la protéine n'a pas d'importance *)
    if (not bp) then
      begin
	(* tout en ayant une molécule à suivre (donc un arn), cad que de part et d'autre de l'intron on a un meme arn *)
	if bm then
	  (try (selectsig_arn_simple tcpl_int_ok) with Invalid_argument s -> Common.print_log ("Pb dans le selectsig niveau selectsig_arn_simple et "^s^"\n");[])
	    (* si pas du tout de molécule à suivre alors on fait une séléction selon la matrice de poids *)
	else
	  (try (selectsig_sc_w_matrix tcpl_int_ok prob gseq mAG mGT) with Invalid_argument s -> Common.print_log ("Pb dans le selectsig niveau 1er selectsig_sc_w_matrix et "^s^"\n");[])
      end
    else
      (* si la protéine a une importance alors on fait comme dans exogprot cad que la recherche prend en compte 
	 la préservation de la phase donnée par la protéine (donc mod3 + stop) *)
      begin
	let t_mod3_without_stop = (try (filter_frame gb1 ge1 gb2 ge2 gseq tcpl_int_ok) with Invalid_argument s -> Common.print_log ("Pb dans le selectsig niveau filter_frame et "^s^"\n");[||]) in
	  match mss with
	    |Aamiss ->  (* modif au 27/02/06 -> le n est à chercher dans l'option nbaamanquant *) 
	       (* si on a une protéine a suivre et que le nbaa manquant est inférieur à un seuil *)
	       if (bm && (naa <= n)) then
		 (try (selectsig_minaddlost naa t_mod3_without_stop gb1 ge1 gb2 ge2 gseq) with Invalid_argument s -> Common.print_log ("Pb dans le selectsig niveau selectsig_min_add_lost et "^s^"\n");[])
		   (* si on n'a pas de protéine à suivre ou que le nombre d'acides aminés est superieur au seuil *)
	       else
		 (try (selectsig_minadding t_mod3_without_stop gb1 ge1 gb2 ge2 gseq) with Invalid_argument s -> Common.print_log ("Pb dans le selectsig niveau selectsig_minadding et "^s^"\n");[])
		   
	    |Wghtmat -> (try (selectsig_sc_w_matrix t_mod3_without_stop prob gseq mAG mGT) with Invalid_argument s -> Common.print_log ("Pb dans le selectsig niveau 2eme selectsig_sc_w_matrix et "^s^"\n");[])
      end

  (*with
      Invalid_argument s -> Common.print_log ("Pb dans le selectsig et "^s^"\n");[]*)
	
      

    
   





(* 
   sClusp_to_tEclusp prend en entree s qui est la segseq de clusps issue du blusping 
   (blusping1 seul ou de blusping2(blusping1), la bigseq représentant la séquence d'adn à analyser 
   et la banque de protéines sous forme de table de hachage, et rend un tableau de Eclusps, 
   cad un tableau de Clusps exonés ou encore étiquetés, et dont surtout la liste de signaux
   up et down (par rapport au gène) est exhaustive (cad tous les signaux satisfaisant les
   bonnes propriétées par exemple sur le nombre d'acides aminés séparant deux clusps dont 
   les hsps ppales proviennent de la meme protéine. 
   Attention : il faut calculer au vol la protéine principale du eclusp courant en additionnant
   les acides aminés de tous les matches d'une meme protéine sur le blusp courant.
   En fait bqp ne sert que pour la recherche d'ATG ou de Stop.
   Attention : désormais lposigu et lposigd sont des couples : (nbaarecherche, liste de signaux possibles)
   REMARQUE : dans la version mix arn prot on ne filtre pas les monoecluspiques pour l'instant 
*)   
let sClusp_to_tEclusp nninhsp nnsuppl sc gbseq bqp nbdef =
  Array.mapi 
    (fun i c -> 
       let ppale = fst (prot_pal (SegSeq.map polytypec2protc sc) i) and apal = fst (prot_pal (SegSeq.map polytypec2arnc sc) i) in
       Eclusp.create (Clusp.gbeg c) (Clusp.gend c) (Clusp.lhsp c) (Clusp.brinc c) 
       (clabel sc i) (Clusp.mset c) (Clusp.ppos c) 
       ppale
       apal
       (select_hsp ppale (Clusp.lhspp c)) 
       (select_hsp ppale (Clusp.lhspu c)) 
       (select_hsp ppale (Clusp.lhspd c)) 
       (lposigu nninhsp nnsuppl c sc i ppale gbseq bqp nbdef) 
       (lposigd nninhsp nnsuppl c sc i ppale gbseq bqp nbdef) 
       ((snd (lposigu nninhsp nnsuppl c sc i ppale gbseq bqp nbdef))<>[])
       ((snd (lposigd nninhsp nnsuppl c sc i ppale gbseq bqp nbdef))<>[])
       (Clusp.commc c) (Clusp.pbbrin c)) 
    (SegSeq.tank sc)








(* ec_to_ex prend un eclusp et retourne un exon, prob est un tableau de 5 floats repertoriant les probablilités
   génomiques de A, T, G, C, N dans gseq (calculé au préalable) 
   mss représente la façon dont on veut séléctionner les signaux *)
let ec_to_ex nbaamq tec seuil_int mss prob gseq mAG mGT i ec =
  match (Eclusp.lbl ec) with
    |Unique -> Exon.create 
       (* ATTENTION : ici et partout ou on recherche l'atg il faut prendre le dernier et non le premier
	  signal de la liste donnée par lsigposu car sinon on minimise au lieu de maximiser le nb d'aa 
	  manquant dans la protéine de l'hsp la plus 5'!!!!!!!!!!!!!!!!!! *)
       (try (last (snd (Eclusp.lsigposu ec))) with Failure _ -> (Eclusp.gbeg ec))
	  (* signal up : ATG *)
	  (try ((List.hd (snd (Eclusp.lsigposd ec)))-1) with Failure _ -> (Eclusp.gend ec)) 
	  (* signal down : Stop *)
	  (Eclusp.lhsp ec) (Eclusp.brinec ec) Unique (Eclusp.mset ec) (Eclusp.protpale ec) (Eclusp.arnpal ec) 
	  (snd (Eclusp.lsigposu ec)<>[]) (snd (Eclusp.lsigposd ec)<>[]) (Eclusp.commec ec) (Eclusp.pbbrinec ec);
	
    |Initial -> Exon.create 
       (* ATTENTION : ici et partout ou on recherche l'atg il faut prendre le dernier et non le premier
	  signal de la liste donnée par lsigposu car sinon on minimise au lieu de maximiser le nb d'aa 
	  manquant dans la protéine de l'hsp la plus 5'!!!!!!!!!!!!!!!!!! *)
       (try (last (snd (Eclusp.lsigposu ec))) with Failure _ -> (Eclusp.gbeg ec)) 
	  (* signal up : ATG *)
	  (try ((fst (List.hd (selectsig nbaamq (Eclusp.lsigposd ec) (Eclusp.lsigposu tec.(i+1)) seuil_int mss (Eclusp.gbeg ec) (Eclusp.gend ec) (Eclusp.gbeg tec.(i+1)) (Eclusp.gend tec.(i+1)) prob gseq  mAG mGT )))-1) with Failure _ -> (Eclusp.gend ec)) 
	  (* signal down : GT *)
	  (Eclusp.lhsp ec) (Eclusp.brinec ec) Initial (Eclusp.mset ec) (Eclusp.protpale ec) (Eclusp.arnpal ec) 
	  (snd (Eclusp.lsigposu ec)<>[]) ((selectsig nbaamq (Eclusp.lsigposd ec) (Eclusp.lsigposu tec.(i+1)) seuil_int mss (Eclusp.gbeg ec) (Eclusp.gend ec) (Eclusp.gbeg tec.(i+1)) (Eclusp.gend tec.(i+1)) prob gseq  mAG mGT )<>[]) (Eclusp.commec ec) (Eclusp.pbbrinec ec);
	
    |Internal -> Exon.create 
       (try ((snd (List.hd (selectsig nbaamq (Eclusp.lsigposd tec.(i-1)) (Eclusp.lsigposu ec) seuil_int mss 
			      (Eclusp.gbeg tec.(i-1)) (Eclusp.gend tec.(i-1)) (Eclusp.gbeg ec) (Eclusp.gend ec) prob gseq  mAG mGT )))+2) with Failure _ -> (Eclusp.gbeg ec))
	  (* signal up : AG *)
	  (try ((fst (List.hd (selectsig nbaamq (Eclusp.lsigposd ec) (Eclusp.lsigposu tec.(i+1)) seuil_int mss 
				 (Eclusp.gbeg ec)  (Eclusp.gend ec) (Eclusp.gbeg tec.(i+1)) (Eclusp.gend tec.(i+1)) prob gseq  mAG mGT )))-1) with Failure _ -> (Eclusp.gend ec))
	  (* signal down : GT *)
	  (Eclusp.lhsp ec) (Eclusp.brinec ec) Internal (Eclusp.mset ec) (Eclusp.protpale ec) (Eclusp.arnpal ec) 
	  ((selectsig nbaamq (Eclusp.lsigposd tec.(i-1)) (Eclusp.lsigposu ec) seuil_int  mss (Eclusp.gbeg tec.(i-1)) (Eclusp.gend tec.(i-1)) (Eclusp.gbeg ec) (Eclusp.gend ec) prob gseq  mAG mGT )<>[]) ((selectsig nbaamq (Eclusp.lsigposd ec) (Eclusp.lsigposu tec.(i+1)) seuil_int mss (Eclusp.gbeg ec)  (Eclusp.gend ec) (Eclusp.gbeg tec.(i+1)) (Eclusp.gend tec.(i+1)) prob gseq  mAG mGT )<>[])  (Eclusp.commec ec) (Eclusp.pbbrinec ec);
	
    |Terminal -> Exon.create 
       (try ((snd (List.hd (selectsig nbaamq (Eclusp.lsigposd tec.(i-1)) (Eclusp.lsigposu ec) seuil_int mss 
			      (Eclusp.gbeg tec.(i-1)) (Eclusp.gend tec.(i-1)) (Eclusp.gbeg ec) (Eclusp.gend ec) prob gseq  mAG mGT )))+2) with Failure _ -> (Eclusp.gbeg ec))
	  (* signal up : AG *)
	  
	  (try ((List.hd (snd (Eclusp.lsigposd ec)))-1) with Failure _ -> (Eclusp.gend ec)) 
	  (* signal down : Stop *)
	  (Eclusp.lhsp ec) (Eclusp.brinec ec) Terminal (Eclusp.mset ec) (Eclusp.protpale ec) (Eclusp.arnpal ec) 
	  ((selectsig nbaamq (Eclusp.lsigposd tec.(i-1)) (Eclusp.lsigposu ec) seuil_int mss (Eclusp.gbeg tec.(i-1)) (Eclusp.gend tec.(i-1)) (Eclusp.gbeg ec) (Eclusp.gend ec) prob gseq  mAG mGT )<>[]) ((snd (Eclusp.lsigposd ec))<>[])  (Eclusp.commec ec) (Eclusp.pbbrinec ec)
	




(* Comme son nom l'indique tEclusp_to_tExon transforme un tableau de eclusps en un tableau d'exons *)
let tEclusp_to_tExon nbaamq teclusp seuil_int mss prob gseq mAG mGT= 
  Array.mapi (ec_to_ex nbaamq teclusp seuil_int mss prob gseq mAG mGT) teclusp



    
