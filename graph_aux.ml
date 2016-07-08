(* graph_aux.ml : contient les fonctions auxiliaires à la construction du graphe acyclique coloré
   de blusps monoarn et monoproteine 
   Contient des fonctions relatives aux intervalles, aux conditions de liaison entre deux bmm arn

   INITIALEMENT ECRIT PAR MATTHIEU MANCENY au Lami, Evry, mmanceny@lami.univ-evry.fr
*)


open Donnees
open Graph
open Common



module Mg = Make(HashtblArray)


let vert2bmm v = match v with
      Arn b -> b
    |Protn b -> b
    |_ -> Bluspmm.nullb



(* QUELQUES FONCTIONS INTERMEDIAIRES A LA MISE DE LIENS ESTGENES *)

  (* i1 et i2 sont des intervalles, peu importe dans quel ordre et on veut savoir s'ils sont chevauchants.
     à mettre à terme dans donnees.ml au niveau des introns et de la comparaison d'intervalles*)
  let intervalle_chev i1 i2 =
    let (p1,p2) = i1 and (p3,p4) = i2 in
      if ((intervalle_comp i1 i2)<0) then
	(p3>=p1) && (p3<=p2)
      else
	(p1>=p3) && (p1<=p4) 



  (* cherche dans la liste de couples d'entiers (=intervalles introniques) rangés par ordre croissant
     le premier intervalle de li1 qui chevauche l'intervalle i2.
     Si ne trouve pas renvoit l'exception Not_found *)
  let rec find_prem_int_chev li1 i2 =
    match li1 with
	[] -> raise Not_found
      |i1::qi1 -> 
	 if (intervalle_chev i1 i2) then
	   i1
	 else
	   find_prem_int_chev qi1 i2



  (* parmi tous les introns réels de a2 (intervalles de lir2) cherche si ceux qui chevauchent
     un intron réel de a1 (liste lir1 d'intervalles) sont bien égaux à ces derniers  *)
  let rec tous_int_chev_egaux lir1 lir2 =
    match lir2 with
	[] -> true
	  (* si on aboutit à la fin de la liste lir1 c'est que l'on n'a pas trouvé de false -> c'est true *)
      |i2::qir2 -> 
	 try (let i1 = find_prem_int_chev lir1 i2 in
		if (i1=i2) then
		  tous_int_chev_egaux lir1 qir2
		else
		  false
		    (* si on trouve un seul intron de lir1 qui chevauche un intron de lir2 sasn y etre égal -> false *)
	     )
	 with
	     (* si on ne trouve pas d'introns dans lir1 qui chevauchent l'intron i2 de lir2 alors on
		poursuit la vérification dans le reste des introns réels de lir2, cad qir2 *)
	   Not_found -> tous_int_chev_egaux lir1 qir2 



  (* fonction qui vérifie que deux bluspmm b1 et b2 ne clashent pas, cad qu'aucun de leurs introns
     réels chevauchant ne clashent *)
  let do_not_clash b1 b2 =
    let lir1 = List.map (fun (Real x) -> x) (List.filter realint (Bluspmm.il b1)) and lir2 = List.map (fun (Real x) -> x) (List.filter realint (Bluspmm.il b2)) in
      tous_int_chev_egaux lir1 lir2
 

  (* Il est plus stratégique de voir la condition d'extension apres car la condtion d'introns egaux 
     sera plus vite fausse -> + rapide
     Attention : ici contrairement à donnees.ml c'est le premier element qui étend le deuxième
     ce qui est fait pour que la fonction for_all de Graph.pred puisse marcher
     donc b1 ancien = bisext
          b2 ancien = bext *)
  let cond_ext bext bisext =
    (do_not_clash bext bisext) && (Bluspmm.etend bisext bext)



  let cond_incl bincl bisincl =
    (do_not_clash bincl bisincl) && (Bluspmm.inclus bisincl bincl)
       
  

  (* tb est un tableau de Bluspmm arn, cad un cluster d'arn chevauchants 
     et graphe est le graphe dans lequel on souhaite mettre les liens 
     d'extension et d'inclusion entre les blsupmm arn (cf estgenes).
     Au départ il y a une seule feuille, c'est tb.(0) (on suppose qu'il 
     y a au moins un blusp b dans le cluster tb) puis on ajoute tb.(1)
     en regardant s'il étend ou est inclus dans tb.(0) ...etc.

     Il existe trois listes de Bluspmm feuilles : lf, lfext et lfincl,
     lf = liste des feuilles 
     lfext = listes des feuilles étendues jusque là (par tb.(i))
     lfincl = listes des feuilles qui incluent un element jusque là (tb.(i))
     
     Pour chaque nouvel element tb.(i) on va le comparer aux feuilles de lf:
     - quand il en etend une on ajoute un lien d'extension entre cette feuille
     et lui et on ajoute cette feuille à la liste lfext (des feuilles etendues)
     - quand il est inclus dans une on ajoute un lien d'inclusion entre cette feuille
     et lui et on ajoute cette feuille à la liste lfincl (des feuilles qui incluent)

     Quand on a fini de traiter l'élément courant tb.(i) :
     - on soustrait de lf (liste des feuilles) tous les élements de lfext
     (les feuilles étendues)
     - on n'ajoute tb.(i) à la liste des feuilles que si 
        * soit tb.(i) a étendu une feuille (nbext!=0)
        * soit tb.(i) n'a étendu personne et n'a été inclus dans personne (nfext=0 et nfincl = 0)
        * on réinitialise à 0 le nombre de feuilles étendues et incluses et on ajuste le nombre de feuilles
          totale (nf :=List.length (!lf)) et le nombre de feuille courante (j=0) car important pour la
          condition d'arret de la boucle interne (!j< (!nf))

*)
  let lie_arn graphe tb =
    try
      ( (* tableaux des labels des vertices et des bluspmm arn selon l'index *)
	let tvlab = Mg.Vertex.label graphe Nullv in
	let tblab = Array.map vert2bmm tvlab in
	  
	let i = ref 1 and n = Array.length tb and lf = ref [tb.(0)] and j = ref 0 and nf = ref 1 and lfext = ref [] and lfincl = ref [] and nfext = ref 0 and nfincl = ref 0 in
	  while (!i < n) do
	    while (!j< (!nf)) do
	      if ((cond_ext tb.(!i) (List.nth !lf !j)) 
		  && 
		  (Mg.Pred.forall (fun idx -> (do_not_clash tb.(!i) tblab.(idx))) graphe (Arn (List.nth !lf !j)))) then
		begin
		  Mg.Edge.set graphe (Arn (List.nth !lf !j)) (Arn tb.(!i)) Extension;
		  lfext := (List.nth !lf !j)::!lfext;
		  incr nfext;
		end
	      else
		begin
		  if ((cond_incl tb.(!i) (List.nth !lf !j)) 
		      && 
		      (Mg.Pred.forall (fun idx -> (do_not_clash tb.(!i) tblab.(idx))) graphe (Arn (List.nth !lf !j)))) then
		    begin
		      Mg.Edge.set graphe (Arn (List.nth !lf !j)) (Arn tb.(!i)) Inclusion;
		      lfincl := (List.nth !lf !j)::!lfincl;
		      incr nfincl;
		    end
		end;
	      incr j;
	    done;  (* fin du while (!j< (!nf)) *)
	    
	   (* Modification de la liste des feuilles courantes
	      ***********************************************
	      Retrait des feuilles d'ou partent des liens d'extension *)
	    lf := remove_all (!lfext) (!lf); 
	    
	    (* Ajout de tb.(i) = element courant parmi les nouvelles feuilles 
	     **************************************************************
	     * Se produit dans deux cas : 
	     * soit lorsque tb.(i) étend une feuille (!nfext!=0)
	     * soit lorsque tb.(i) n'a été lié par aucun lien aux autres (ie ni extension ni inclusion)
	     * En ce cas les autres restent feuilles et lui le devient *)
	    if((!nfext!=0) || ((!nfext=0) && (!nfincl=0))) then
	      lf := (tb.(!i))::(!lf);
	    
	    (* Réinitialisation pour chaque nouvel arn tb.(i) 
	     ************************************************  *)
	    nf := List.length (!lf); (* nombre total de feuilles *)
	    nfext := 0;  (* nombre de feuilles étendues *)
	    nfincl := 0;  (* nombre de feuilles incluantes *)
	    incr i; (* indice courant pour les blusps à traiter dans tb *)
	    j:=0; (* indice courant pour les feuilles, varie de 0 à nf strictement *)
	  done) (* fin du while (!i < n) *)
    with
	Invalid_argument s -> ()










 (* let cond_ext bext bisext =
    let lirisext = List.map (fun (Real x) -> x) (List.filter realint (Bluspmm.il bisext)) and lirext = List.map (fun (Real x) -> x) (List.filter realint (Bluspmm.il bext)) in
      (tous_int_chev_egaux lirisext lirext) && (Bluspmm.etend bisext bext)
 *)
 (* let cond_incl bincl bisincl =
    let lirisincl = List.map (fun (Real x) -> x) (List.filter realint (Bluspmm.il bisincl)) and lirincl = List.map (fun (Real x) -> x) (List.filter realint (Bluspmm.il bincl)) in
      (tous_int_chev_egaux lirisincl lirincl) && (Bluspmm.inclus bisincl bincl)
 *)
