(* premain.ml : fichier definissant l'interface utilisateur qu'est le contexte ainsi que la fonction essentielle
   find_modele qui trouve les modeles de genes sur la séquence génomique fournie dans le contexte c en fonction 
   des hsps fournies egalement dans c et donne egalement pour chaque modele sa traduction en protéine *)


open Alphaprot
open Alphaadn
open Collection
open Donnees_base
open Normalisation_mrna
open Donnees
open Classifgene
open Db
open TraitementHSP
open Pretraitement1
open Pretraitement2
open Orfsearch
open Matrice
open Graph
open Graph_aux
open TraiteCollection
open Sigsearch
open Common
open Config
open Seq
open Printing




(********************************************************************************************************************
             Fonction find_modeles auxiliaire à la fonction print_tous_genes principale de main.ml
*********************************************************************************************************************)

(* c est un contexte contenant 9 champs mutables vus plus haut (dont 3 obligatoires) *)
let find_modeles c =
  
  (* à terme il faut le mettre dans le contexte et vérifier que cela marche (car là non) *)
  let ajuste=true in

 (* CREATION DES MATRICES DE POIDS POUR LES SIGNAUX D'EPISSAGE *)
 let mAG = try (cree_matrix c.magfile 0) with
   | Sys_error _ -> mAGdefault
 and mGT = try (cree_matrix c.mgtfile 0) with
   | Sys_error _ -> mGTdefault in
    (* Rq : non nécessaire qd c.waysigsearch != WghtMat *)
   
 (* filtre des bouts d'arn apres découpe soit au pretraitement soit à l'intron non épissé 
    paramétrable car on souhaite etre - stringent avec les ests? *)
 let pfilter1 lh = (List.length lh) >=1 in
 let pfilter2 lh = (List.length lh) >=2 in


 (* LECTURE DE LA SEQUENCE GENOMIQUE, COMPLEMENTATION, PROB GENOMIQUES *)
    (* Lecture de la séquence génomique à annoter, 
       Attention : gseq est un tableau de nucléotides 
       Attention : Etape qui prend plus de temps que les autres *)
 let (nomgseq, lggseq, gseq) = List.hd (readg fasta of_chara nulla (open_in c.gseqfile)) in 

 let compgseq = complement gseq in
    

 (* CALCUL DES PROBABILITES GENOMIQUES *)
 let prob_seq = probATGC gseq in
   Common.print_log "I have read the genomic file\n";
   

   (* LECTURE DES HSPS PROT EN DISTINGUANT LES PROTEINES MONOHSP ET LE BRIN *)
   
   let filterlh = if (c.protmonohsp) then pfilter1 else pfilter2 in 
   let filterh h = ((Hsp.size h) >= c.szhsp) in
   
   let (fhspp,filehspp) = get_type_entree c.hsppfile in
   
    let ((protlist_nonart,nbprot_all),(llhpp,llhmp)) = 
     try (partageListeHsp filterh filterlh (Array.of_list (List.sort Hsp.compnameprot (fhspp Prot (open_in filehspp)))) lggseq) 
     with Invalid_argument _ | Sys_error _ -> (([],0),([],[])) in
     (* Liste de liste d'hsps *)
    
     Common.print_log ("I have read "^(string_of_int nbprot_all)^" proteins in the protein alignment file\n");
     Common.print_log ("I have kept "^(string_of_int (List.length (List.flatten llhpp)))^" proteic hsps on the + strand\n");
     Common.print_log ("I have kept "^(string_of_int (List.length (List.flatten llhmp)))^" proteic hsps on the - strand\n");


     (* LECTURE DES HSPS ARNM EN DISTINGUANT LES ARNM MONOHSP ET LE BRIN *) 
     let (fhspa,filehspa) = get_type_entree c.hspafile in
       (* quoiqu'il arrive on ne garde pas les arn monohsp d'ou le pfilter2
	  de plus on ne filtre pas sur la taille des hsps arn d'ou le fun h -> true 
       *)
     let ((arnlist_nonart,nbarn_all),(llhpa,llhma)) =  try (partageListeHsp (fun h -> true) pfilter2 (Array.of_list (List.sort Hsp.compnameprot (fhspa Rna (open_in filehspa)))) lggseq) with Invalid_argument _ | Sys_error _ -> (([],0),([],[])) in
       (* Liste de liste d'hsps *)
       

       Common.print_log ("I have read "^(string_of_int nbarn_all)^" lines in the mrna alignment file\n");
       
       let tidxp = Array.of_list protlist_nonart in
       let nbprot_nonart = List.length protlist_nonart and nbarn_nonart = List.length arnlist_nonart in
	 

  
       (* RQ : ICI S'ARRETE LA SIMILITUDE ENTRE LES TRAITEMENTS PROT ET ARNM *)
	 
       (* LECTURE DE LA BANQUE DE SEQUENCES PROTEIQUES *)
       let (bqp,tlgprot) = try (readp fasta of_charp nullp (open_in (c.bqpfile)) protlist_nonart) 
       with Sys_error _ -> ((Hashtbl.create 0),(Array.make 0 0)) in 
	 (* Lecture des protéines qui ont matché sur le génomique 
	    Paire : (Table de hachage des protéines utilisées pour l'analyse, tableau de longueur des protéines) 
	    Remarque : bqp est une table de hachage de clé un nom de protéine et de contenu le tableau
	    d'acides aminés constituant sa séquence
	    lpwant est la liste qui permet de passer d'un index protéique au nom de la protéine *)
	 Common.print_log ("I have read the proteins\n");
	  

	 (* PRETRAITEMENT DES HSPS PROTEIQUES *)
	 let thspajustp = if (ajuste) then 
	   (Array.map (adjust_lhsps correct_window gseq bqp) (Array.of_list llhpp))
	 else 
	   (Array.of_list llhpp)
	 
	 and thspajustm = if (ajuste) then 
	   (Array.map (adjust_lhsps correct_window compgseq bqp) (Array.of_list llhmp))
	 else
	   (Array.of_list llhmp)
	 in

	   Common.print_log ("There are "^(string_of_int (List.length (List.flatten (Array.to_list thspajustp))))^" protein hsps adjusted on the + strand\n");
	   Common.print_log ("There are "^(string_of_int (List.length (List.flatten (Array.to_list thspajustm))))^" protein hsps adjusted on the - strand\n");
	    
	  (* Array.iteri 
	     (fun i h -> 
		 Common.print_log 
		   (("Hsp ")^(string_of_int i)^(" : debg = ")^(string_of_int (Hsp.debg h))^(" ; fing = ")^(string_of_int (Hsp.fing h))^("  - debp = "^(string_of_int (Hsp.debp h))^" ; finp = "^(string_of_int (Hsp.finp h))^"\n")^(Hsp.commh h)^("\n"))
	     )
	     (Array.of_list thspajustm.(0));*)
	   
	   let thsppret1p = Array.map (fusremp Prot gseq regroupep c.bplusa c.diffaa c.fusionthresh bqp) thspajustp 
	   and thsppret1m = Array.map (fusremp Prot compgseq regroupep c.bplusa c.diffaa c.fusionthresh bqp) thspajustm in
	      (*
		Deux tableaux (un pour chaque brin) contenant des listes d'hsps protéiques prétraitées, 
		une pour chaque protéine non artéfactuelle, et triées selon le début génomique (et fin génomique si 
		deb identique)
	      *)
	      Common.print_log "I have pretreated the protein hsps\n";
	      
	      Common.print_log ("There are "^(string_of_int (List.length (List.flatten (Array.to_list thsppret1p))))^" pretreated protein hsps on the + strand\n");
	      Common.print_log ("There are "^(string_of_int (List.length (List.flatten (Array.to_list thsppret1m))))^" pretreated protein hsps on the - strand\n");
	 
	     (* Array.iteri 
		(fun i h ->  
		   Common.print_log 
		   (("Hsp ")^(string_of_int i)^(" : debg = ")^(string_of_int (Hsp.debg h))^(" ; fing = ")^(string_of_int (Hsp.fing h))^("  - debp = "^(string_of_int (Hsp.debp h))^" ; finp = "^(string_of_int (Hsp.finp h))^"\n")^(Hsp.commh h)^("\n"))
		)
		(Array.of_list (thsppret1m.(0)));
	     *)

	      (* BLUSPING MONOMOLECULE DES HSPS PROT PRETRAITEES = BLUSPING2 PARALLELISE *)
	      
	      (* D'abord clusping des hsps (artificiel ici mais pour pouvoir utiliser le blusping2) *)
	      let lshspp = List.filter (fun s -> not (SegSeq.isnull s)) (Array.to_list (Array.map SegSeq.single (Array.map (Array.of_list) thsppret1p))) 
	      and lshspm = List.filter (fun s -> not (SegSeq.isnull s)) (Array.to_list (Array.map SegSeq.single (Array.map (Array.of_list) thsppret1m))) in
		(* Deux listes de segseq d'hsps (un par brin), 
		   chaque segseq d'hsps représentant une protéine ayant matché *)
		
  
		
	      let lshsppseg = List.map clusping lshspp and lshspmseg = List.map clusping lshspm in 
		(* deux listes de segseq, une segseq par protéine *)
		
		Common.print_log "I have made one segmented sequence (segseq) for each protein\n";
		
		let llthspp = List.map SegSeq.elements lshsppseg and llthspm = List.map SegSeq.elements lshspmseg in
		  (* deux listes de listes de tableaux de 1 hsp chacun (artificiel) *) 
		
		  Common.print_log "I have done the elements on each segseq that represent a protein\n";
		  
		  let llcluspp = List.map (List.map thspchev_to_clusp) llthspp and llcluspm = List.map (List.map thspchev_to_clusp) llthspm in
		    (* deux listes de listes de clusps, une par protéine. Rq : la plupart du temps il y a une hsp par clusp sauf 
		     deux fois ici ou on a 2 hsps collées par clusp 
		       Rq : c'est apres la coupure du blusping que l'on doit faire l'elimination des artefacts *)
		    
		    Common.print_log "I have done the clusping of the protein hsps\n";
		  

		    (* Puis blusping monomolécule protéique, cad regroupement des hsps protéiques de meme protéine 
		       qui ont de plus les propriétés suivantes :
		       - taille intron génomique inférieure à un seuil
		       - ordre génomique/protéique des hsps adjacentes sur le génomique

		       Puis élimination des blusps mono et dicluspiques artéfactuels 
		      
		       Attention : au 15/03/06 : modif de monoclaret et diclart dans traiteCollection.ml
		       de façon à ce que si le mono ou dicluspique artefactuel porte l'indication
		       d'un problème de cutdown (du à intron trop grand ou ordre gx molécule)
		       on garde le mono ou di art et donc l'information
		       
		       Au 17/03/06 : 
		       (nonmonoclart c.propmol bqp tidxp tlgprot tcl) 
		       && 
		       (nondiclart c.propmol bqp tidxp tlgprot tcl)
		       
		       remplacé par 

		       (nonmonoprotart c.propmol bqp tidxp tlgprot tcl)
		    *)
		    let lltcp = List.map (fun lc -> 
					    List.filter 
					    (fun tcl -> 
					       (nonmonoprotart c.pnntfus c.pnaamq c.propmol c.propmol2 bqp tidxp tlgprot tcl)
					    )
					    (SegSeq.elements 
					       (blusping2monomol2 gporder_tintok_hspisolok
						  (SegSeq.make3 (Array.of_list lc) [0]) c.szisolh c.tmaxintron c.tmaxintron2 c.beta
					       )
					    )
					 )
				  llcluspp 
				  (* 
				     liste de listes de clusps de meme protéine ou liste de bluspmm protéique potentiel
				     au sein de chaque telle liste de clusps de meme protéine (un clusp = un hsp)
				     on regarde s'il faut encore couper entre les clusps dans les deux cas suivants :
				     - deux hsps contigus trop distantes sur le génomique
				     - deux hsps contigus avec un retour en arrière dans la protéine 
				     
				     Pour chaque liste de clusps représentant potentiellement un bluspmm protéique 
				     on obtient en sortie une liste de tableaux de clusp, et ce nombre de tableaux
				     correspond à 1 + le nombre de fois où il a fallu couper pour une des raisons 
				     précédantes
				  *)

		    and lltcm = List.map (fun lc -> 
					    List.filter 
					    (fun tcl -> 
					       (nonmonoprotart c.pnntfus c.pnaamq c.propmol c.propmol2 bqp tidxp tlgprot tcl) 
					    )
					    (SegSeq.elements 
					       (blusping2monomol2 gporder_tintok_hspisolok
						  (SegSeq.make3 (Array.of_list lc) [0]) c.szisolh c.tmaxintron c.tmaxintron2 c.beta
					       )
					    )
					 )
				  llcluspm in
		      
		      Common.print_log "I have done the single molecule clustering for the proteins\n";

		      (* deux listes de listes de tableaux de clusps, une pour chaque brin.
			 Chaque tableau de clusps représente un futur blusp monomolécule. 
			 On fait à ce stade une élimination des blusps protéiques monocluspiques
			 et discluspiques etranges!
			 Attention on utilise gporder_tintok qui lui-meme utilise gporderf2, 
			 pour les deux brins car les pos génomiques des données brin - sont en fonction 
			 du début du brin - *)

		      let llthp = List.map (fun ltc ->
					      List.map 
					      (Array.map 
					       (fun c -> List.hd (Clusp.lhspp c)
					       )
					      )
					      ltc
					   )
				  lltcp 
				    
		      and llthm = List.map (List.map (Array.map (fun c -> List.hd (Clusp.lhspp c)))) lltcm  in
			(* 
			   On repasse de clusp à hsp car c'est + simple (on loupe juste les deux hsps 
			   collées qui forment un seul clusp). 
			   Deux listes de listes de tableaux d'hsps, une par brin 
			   tous les blusps d'une meme protéine monobrin
			   sont donc représentés par une liste de tableaux d'hsps, 
			   un tableau par blusp de cette protéine 
			*)

		      let llthcutdwnp = List.map (fun lthp -> let longthsp = List.length lthp in
						    if(longthsp=1) then (* si un seul tableau d'hsps
									   signifie pas eu de coupure 
									 -> Rien à faire 
									*)
						      lthp
						    else  (* si au moins deux tableaux d'hsps pour une meme protéine *)
						    begin
						      Array.to_list 
							(Array.mapi 
							   (fun i thsp -> 
							      if (i!=(longthsp-1)) then  
								(* si longthsp tableaux d'hsps, signifie qu'il y a eu
								   longthsp-1 coupures, et il faut donc modifier 
								   les derniers hsps de chaque tableau d'hsp sauf 
								   le dernier hsp du dernier tableau
								   car ne représente pas une cassure
								*)
								begin
								  let longt = Array.length thsp in
								    Array.mapi 
								      (fun j h -> 
									 (* pour chaque tableau on ne modifie 
									    que son dernier hsp *)
									 if (j!=(longt-1)) then 
									   h
									 else
									   Hsp.setcutdwn h  (* j=longt-1, 
											       dernier hsp du tableau *)
								      )
								      thsp
								end
							      else  (* i= longthsp-1, dernier tableau, rien à y faire *)
								thsp
							   ) 
							   (Array.of_list lthp))  (* en fait ces deux cas 
										     (longthsp=1 ou !=1) peuvent 
										     etre regroupés en un seul *)
						    end
					       ) 
					  llthp in

		      let llthcutdwnm = List.map (fun lthp -> let longthsp = List.length lthp in
						    if(longthsp=1) then (* si un seul tableau d'hsps
									   signifie pas eu de coupure 
									 -> Rien à faire 
									*)
						      lthp
						    else  (* si au moins deux tableaux d'hsps pour une meme protéine *)
						    begin
						      Array.to_list 
							(Array.mapi 
							   (fun i thsp -> 
							      if (i!=(longthsp-1)) then  
								(* si longthsp tableaux d'hsps, signifie qu'il y a eu
								   longthsp-1 coupures, et il faut donc modifier 
								   les derniers hsps de chaque tableau d'hsp sauf 
								   le dernier hsp du dernier tableau
								   car ne représente pas une cassure
								*)
								begin
								  let longt = Array.length thsp in
								    Array.mapi 
								      (fun j h -> 
									 (* pour chaque tableau on ne modifie 
									    que son dernier hsp *)
									 if (j!=(longt-1)) then 
									   h
									 else
									   Hsp.setcutdwn h  (* j=longt-1, 
											       dernier hsp du tableau *)
								      )
								      thsp
								end
							      else  (* i= longthsp-1, dernier tableau, rien à y faire *)
								thsp
							   ) 
							   (Array.of_list lthp))  (* en fait ces deux cas 
										     (longthsp=1 ou !=1) peuvent 
										     etre regroupés en un seul *)
						    end
					       ) 
					  llthm in
			(*
			  Ici on souhaite retenir là où ont eu lieu les coupures au sein des bluspmm 
			  monoprot potentiels, dues à l'une ou l'autre des deux raisons suivantes :
			  - taille intron max 
			  - ordre gx-molécule

			  On va retenir ce pb au niveau de l'hsp qui a précédé une coupure
			  (champ cutdwn de l'hsp à true), et pour savoir quels hsps sont concernés, 
			  pour chaque liste de tableaux d'hsps il suffit de regarder si on a un tableau
			  ou plusieurs dans cette liste. En effet si on a un seul tableau d'hsps cela signifie
			  que le bluspmm protéique initil est resté intact (pas de coupures) 
			  alors que si on a n>=2 tableaux d'hsps dans une liste, cela signifie qu'il y a eu
			  n-1 coupures dues à l'une des deux raisons précédantes.
			  (en effet à l'étape suivante (une liste d'hsps par bluspmm après coupure) 
			  c'est trop tard car on a déjà applatit!)
			  
			*)
    
		      let llhp = List.flatten (List.map (fun lt -> List.map (Array.to_list) lt) llthcutdwnp) 
		      and llhm = List.flatten (List.map (fun lt -> List.map (Array.to_list) lt) llthcutdwnm) in
			(* On applatit chaque liste de tableau correspondant à une protéine 
			   car on souhaite former autant de blusps monoprotéiques qu'il y a de 
			   tableaux dans la liste représentant la protéine (du aux coupures)
		      *)
			
			(*	List.iter (print_hsp stderr) (List.flatten llhp); *)
				

			Common.print_log ("There are "^(string_of_int (List.length (List.flatten llhp)))^" protein hsps after the single molecule clustering and the elimination of the artifactual single molecule clusters on the + strand\n");
			Common.print_log ("There are "^(string_of_int (List.length (List.flatten llhm)))^" protein hsps after the single molecule clustering and the elimination of the artifactual single molecule clusters on the - strand\n");
			

			let lbp = List.map lh2bluspmm llhp and lbm = List.map lh2bluspmm llhm  in
			  (* On passe de listes d'hsps aux objets bluspmm *)
			  
			let lbpord = List.sort Bluspmm.compare lbp and lbmord = List.sort Bluspmm.compare lbm in
			  (* On ordonne les bluspmm proteiques suivant le génomique croissant, ie 5', puis 3', bien séparé *)
			  

			  (* Common.print_log "J ai fait les bluspmm monoprot\n"; *)
			  Common.print_log "I have done the single molecule clustering for proteins\n";
			  
			  


			  (*************** PRETRAITEMENT DES HSPS ARNM ET BLUSPING MONOMOLECULE ARN ***********************
			    le retour en arrière dans la molécule et l'alignement à des endroits tres différents
			    (à la fois en brin et en distance génomique) sont réglés dans le prétraitement et le blat
			    respectivement *)
			  
			  (* On élimine les exons de borne interne non correcte --> dans la liste de liste lp 
			     on cherchera GT/AG sur gseq, dans lm on cherchera aussi GT/AG mais sur comgseq 
			     remarque : ici on perd la notion de tri selon le transcrit car 1 meme alignement de transcrit
			     a pu donner lieu à plusieurs alignements de transcrits apres la verif du splice, et on ne fait pas 
			     cela parallement sur l'autre brin ==> du coup dans la liste lp_ok_splice il n'y a pas le meme
			     nombre de listes que dans la liste lm_ok_splice ==> on enlève les listes vides par un filtre 
			     Remarque : 1) l(p/m)avirt_ok_splice est une liste de couples (listes d'hsps arnm, no hsp précédant intron virtuel), un par transcrit non monoex + ok splice 
			     2) l(p/m)avirt_ok_imbric est une liste de couples (listes d'hsps arnm, no hsp précédant intron virtuel), un par transcrit non monoex + ok splice + ok imbric
			     3) Plus tard on pourra mettre dans une meme ligne les deux pretraitements (splice et intron) 
			  *)
			  
			let llhparninterm = List.map (fun lhsp_one_mrna -> let lhelimintpetmod3 = fusremp Rna gseq regroupea c.bplusa c.diffaa c.fusionthresh bqp lhsp_one_mrna in List.filter pfilter2 (cut_one_to_n_list (show_no_hsp_bad_splice gseq lhelimintpetmod3) lhelimintpetmod3)) llhpa 
					      
			and llhmarninterm = List.map (fun lhsp_one_mrna -> let lhelimintpetmod3 = fusremp Rna compgseq regroupea c.bplusa c.diffaa c.fusionthresh bqp lhsp_one_mrna in List.filter pfilter2 (cut_one_to_n_list (show_no_hsp_bad_splice compgseq lhelimintpetmod3) lhelimintpetmod3)) llhma in
			  (* liste de liste de liste d'hsps arn pretraitées (intron petit mod3 + bornes d'epissage)
			     auxquels on a éliminé les blusps avec trop peu d'hsps *)
			  
			 (* List.iter (fun lhpa -> List.iter (fun h -> ) lhpa
				    )
			    llhparninterm; *)
			let lhpavirt_ok_splice = List.filter (fun l -> pfilter2 (fst l)) (List.map arncut_to_arn_virt llhparninterm) 
			and lhmavirt_ok_splice = List.filter (fun l -> pfilter2 (fst l)) (List.map arncut_to_arn_virt llhmarninterm) in

			(* attention on a au milieu une liste de liste de listes d'hsps à cause du fait que l'on coupe 
			   une liste d'hsps de meme arn en plusieurs (cut_one_to_n_list) 
			   Rq : on avait prévu d'enlever aussi les hsps bordées par un intron < 60 pb mais en fait 
			   ok_splice les enleve deja en majorité, en plus il faudrait le faire avant de faire les 
			   introns virtuels, donc avoir une meme signature pour show_no_hsp_bad_splice et 
			   show_no_hsp_imbric, ce qui n'est pas le cas, ou bien le faire là,
			   apres lpavirt_ok_splice mais cela suppose que l'on manipule des couples dans lesquels 
			   on change le 2eme element en fonction de la découverte ou pas d'hsps trop rapprochées -> + tard 
			*)
			  

  
			let lbap = List.map avirt2bluspmm lhpavirt_ok_splice 
			and lbam = List.map avirt2bluspmm lhmavirt_ok_splice in
			  (* on passe des arn virtuels = couple de listes aux bluspmm, pour etre commun proteine *)
			  
			let lbapord = List.sort Bluspmm.compare lbap 
			and lbamord = List.sort Bluspmm.compare lbam in
			  (* on ordonne les blusps arn selon le génomique croissant, ie 5', puis 3', bien séparé 
			     2102 sur brin + et 1863 sur brin - *)



			(* PB DE L'INTRON NON EPISSE : HORS DU GRAPHE. 
			   Entre le pretraitement et le graphe dans tous les sens
			   pas un simple pretraitement comme le pb de l'intron trop court 
			   ou des bornes d'epissage incorrectes
			   car ici on doit comparer des arns chevauchants entre eux, 
			   mais comme la relation/comparaison n'est
			   pas orientée on le fait hors du graphe futur 
  
			   Problème d'intron non épissé (entre bluspmm arn uniquement) : en dehors du graphe car
			   1) changements voire élimination de bluspmm à sa suite
			   2) relation non orientée car il faut regarder l'inclusion de tous les introns du 1er 
			   avec tous les hsps du 2eme arn et l'inclusion de tous les introns du 2eme avec tous 
			   les hsps du 1er 
			*) 


			let llbp = List.map (Array.to_list) (SegSeq.elements (cluspingbpm (SegSeq.single (Array.of_list lbapord)))) 
			and llbm = List.map (Array.to_list) (SegSeq.elements (cluspingbpm (SegSeq.single (Array.of_list lbamord)))) in
			  (* liste des listes de bluspmm arn qui se chevauchent sur le génomique *)
  

			let lllhp = List.map (fun lb ->  List.map (Bluspmm.hl) lb) llbp 
			and lllip = List.map (List.map (Bluspmm.il)) llbp 

			and lllhm = List.map (fun lb ->  List.map (Bluspmm.hl) lb) llbm 
			and lllim = List.map (List.map (Bluspmm.il)) llbm in
			  (* Liste de liste de liste d'hsps, une par bluspmm, idem pour les introns *)


			let lllirealp = List.map (fun lli -> List.map (fun li -> List.filter realint li) lli) lllip 
			and lllirealm = List.map (fun lli -> List.map (fun li -> List.filter realint li) lli) lllim in

			let lllirealposp = List.map (fun lli -> List.map (fun li -> List.map (fun (Real(a,b)) -> (a,b)) li) lli) lllirealp 
			and lllirealposm = List.map (fun lli -> List.map (fun li -> List.map (fun (Real(a,b)) -> (a,b)) li) lli) lllirealm in


			let llhp = List.map (List.fold_left (tri_fusion Hsp.compare) []) lllhp 
			and llip = List.map (remove_redond Pervasives.compare) (List.map (List.fold_left (tri_fusion intervalle_comp) []) lllirealposp) 
			
			and llhm = List.map (List.fold_left (tri_fusion Hsp.compare) []) lllhm 
			and llim = List.map (remove_redond Pervasives.compare) (List.map (List.fold_left (tri_fusion intervalle_comp) []) lllirealposm) in
			  (* Pour chaque cluster de bluspmm arn on liste les hsps de tous les bluspmm arn dans llh
			     et les introns de tous les bluspmm arn dans lli
			     avant : List.map (remove_redond Hsp.compare) au début *)


			let llhpbp = List.map2 (fun lh li -> List.filter (ishapbintron li) lh) llhp llip 
			and llhpbm = List.map2 (fun lh li -> List.filter (ishapbintron li) lh) llhm llim in
			  (* Pour chaque bluspmm arn on a une liste des hsps correspondant à des introns non épissés
			     en tout cas des hsps qui incluent totalement un intron réel (voir + tard pour
			     les hsps terminales qui chevauchent partiellement un intron réel *)



			let llbnewp = List.map2 (fun lhpb lb -> List.filter (fun b -> b!= Bluspmm.nullb) (List.map (b2bnew pfilter2 lhpb) lb)) llhpbp llbp 
			and llbnewm = List.map2 (fun lhpb lb -> List.filter (fun b -> b!= Bluspmm.nullb) (List.map (b2bnew pfilter2 lhpb) lb)) llhpbm llbm in
			  (* une liste de bluspmm arn nouveaux par cluster selon génomique, 
			     en faisant attention de conserver les introns
			     virtuels precédemment crées par le pretraitement (en plus de ceux 
			     crées actuellement par l'intron non épissé *)


 
			let lbanewpord = List.flatten llbnewp and lbanewmord = List.flatten llbnewm in
			  (* la liste des nouveaux bluspsmm arn ordonnés sans pb d'intron non épissé
			     (rq : on a filtré les bluspmm nouveaux nullb -> 1 bluspmm en moins 
			     Rq : il est ici inutile d'ordonner à nouveau les blusps avant d'applatir
			     car dans les listes précedentes on avait des clusters de blusps chevauchants eux-memes
			     ordonnés sur le génomique 
			  *)

			 (* List.iter (fun ba -> List.iter (fun int -> match int with
							      Real(d,f) -> Common.print_log ("realideb = "^(string_of_int d)^" - realifin = "^(string_of_int f)^"\n");
							    | Virt(d,f) -> Common.print_log ("virtideb = "^(string_of_int d)^" - virtifin = "^(string_of_int f)^"\n");
							    | Nullint ->  Common.print_log ("Nullint\n");
							 ) 
				       (Bluspmm.il ba))
			    lbanewpord;
			 *)
			(* On constitue le tableau des arns pretraités à passer en sortie de find_modele et dans le main *) 
			let lparnpretgene = List.flatten (List.map bpmpret2genes lbanewpord) 
			and lmarnpretgene = List.flatten (List.map bpmpret2genes lbanewmord) in
			let tgenearnpret = Array.of_list (List.flatten (lparnpretgene::[lmarnpretgene])) in



			(* MISE DANS UN GRAPHE ARN PAR TABLE DE HACHAGE DE CES BLUSPS MONOMOLECULES (BLUSPMM) ARNM
			   EN TANT QUE SOMMETS (UN GRAPHE PAR BRIN) *)
  
			let nbpa = List.length lbanewpord and nbma = List.length lbanewmord in
			
			let gpa = Mg.make (nbpa+5) Null and gma = Mg.make (nbma+5) Null in
			  

			(* REMPLISSAGE DU GRAPHE ARN AVEC LES BLUSPMM ARN par ordre de génomique croissant 
			   => les index des sommets du graphe suivront aussi cet ordre *)
			let up = List.iter (fun b -> Pervasives.ignore (Mg.Vertex.set gpa (Arn b) Null)) lbanewpord 
			and um = List.iter (fun b -> Pervasives.ignore (Mg.Vertex.set gma (Arn b) Null)) lbanewmord in
			  
  
			(*****************************************************************************************
			  MISE DES LIENS ARN-ARN = LIENS EXTENSION OU INCLUSION 
			  FORMATION DES FEUILLES ARN = classe de transcrits de structure d'épissage compatible
			  ***************************************************************************************
			  = ALGO ESTGENES ADAPTE, SUR LES BLUSPMM ARN DE CHAQUE CLUSTER DE BLUSPMM HETEROGENES CHEVAUCHANTS
			  Principe : on liste tous les introns réels de tous les bluspmm arn de chaque cluster 
			  et on ne compare qu'eux : des que deux bluspmm arn se chevauchent et ont des introns 
			  chevauchants non identiques -> pas de lien "->" entre ces arn et la comparaison du nouvel 
			  arn a lieu avec eux deux ...
			  pb : pour les arn à un exon? et + également les clusters de bmm à un seul bmm?? -> il y en a 87!! 
			*)

  

			(* CLUSTERING SELON LE GENOMIQUE CROISSANT DES BLUSPMM ARN afin de pouvoir ne relier 
			   entre eux que des arns chevauchants (gain de temps car limitation des appels à lie_arn
			   (voir procédure plus loin) *)


			(* Clustering des bluspmm arns *)
			let sclbpp = cluspingbpm (SegSeq.single (Array.of_list lbanewpord)) 
			and sclbpm = cluspingbpm (SegSeq.single (Array.of_list lbanewmord)) in 
			  
  
			let lordpa = SegSeq.ord sclbpp and lordma = SegSeq.ord sclbpm in
			  (* liste des index des premiers bluspmm arn d'un cluster de bluspmm arn fondé sur le 
			     chevauchement génomique (meme rq que précédemment) *)

			let ltbap = SegSeq.elements sclbpp and ltbam = SegSeq.elements sclbpm in
			  (* liste des tableaux de bluspmm arn qui se chevauchent génomiquement dans le graphe gpa,
			     cad appartiennent au meme cluster *)

			let tbpa = Array.of_list (List.map (fun a -> Arn a) lbanewpord) 
			and tbma = Array.of_list (List.map (fun a -> Arn a) lbanewmord) in
			  (* pour avoir un tableau de type donné par Vertex.label *)
			  

			(* MISE DES LIENS ARN-ARN = D'EXTENSION ET D'INCLUSION, POUR CHAQUE CLUSTER inDEPENDAMMENT *)
			let up = List.iter (lie_arn gpa) ltbap and um = List.iter (lie_arn gma) ltbam in
			  
			(* On retire les noeuds arn inclus, inutiles au parcours du graphe -> voir avec Hugues si ok *)
			let up = Mg.Edge.iter (fun i1 i2 d12 -> match d12 with
						   Inclusion -> Mg.Vertex.remove gpa tbpa.(i2)
						 |_ -> ()) gpa 
			and um = Mg.Edge.iter (fun i1 i2 d12 -> match d12 with
						 Inclusion -> Mg.Vertex.remove gma tbma.(i2)
					       |_ -> ()) gma in 
			  
			let amemp = Mg.Vertex.selecti gpa (fun v -> Mg.Vertex.mem gpa v) 
			and amemm = Mg.Vertex.selecti gma (fun v -> Mg.Vertex.mem gma v) in
			
			  
			let lindokp = ref [] and lindokm = ref [] in
			let up = Array.iteri (fun i x -> if x then lindokp:=i::(!lindokp)) amemp 
			and um = Array.iteri (fun i x -> if x then lindokm:=i::(!lindokm)) amemm in

			let lindoktriep = List.sort (Pervasives.compare) !lindokp 
			and lindoktriem = List.sort (Pervasives.compare) !lindokm in
			let lvokp = List.map (fun i -> tbpa.(i)) lindoktriep 
			and lvokm = List.map (fun i -> tbma.(i)) lindoktriem in
  
			let nbokp = List.length lvokp and nbokm = List.length lvokm in
			  (* Dans lindokp/m sont les index des arns non inclus (donc soit etendent soit sont des racines)
			     idée : on refait un graphe avec seulement ceux-là *)

			let gpaok = Mg.make (nbokp+5) Null and gmaok = Mg.make (nbokm+5) Null in
			  

			(* REMPLISSAGE DU GRAPHE ARN AVEC LES BLUSPMM ARN par ordre de génomique croissant 
			   => les index des sommets du graphe suivront aussi cet ordre *)
			let up = List.iter (fun v -> Pervasives.ignore (Mg.Vertex.set gpaok v Null)) lvokp 
			and um = List.iter (fun v -> Pervasives.ignore (Mg.Vertex.set gmaok v Null)) lvokm in
			
			let up = Mg.Edge.iter (fun i1 i2 d12 -> let a = ref true in
						 if ((List.mem i1 lindoktriep) && (List.mem i2 lindoktriep)) then 
			     a:=Mg.Edge.set gpaok tbpa.(i1) tbpa.(i2) d12; ()) gpa 
			and um = Mg.Edge.iter (fun i1 i2 d12 -> let a = ref true in
						 if ((List.mem i1 lindoktriem) && (List.mem i2 lindoktriem)) then 
						   a:=Mg.Edge.set gmaok tbma.(i1) tbma.(i2) d12; ()) gma in
		

		      
			let nbvp = Mg.Vertex.num gpaok and nbvm = Mg.Vertex.num gmaok in
			let tlabp = Mg.Vertex.label gpaok Nullv and tlabm = Mg.Vertex.label gmaok Nullv in
			let ltp = List.rev (Mg.toposortkernel gpaok) and ltm = List.rev (Mg.toposortkernel gmaok) in
			  (* on part de l'ordre topologique inverse, cad des feuilles *)
			  
			let lttriep = List.sort (Pervasives.compare) ltp and lttriem = List.sort (Pervasives.compare) ltm in
 
			(*
			  let larci = Mg.Edge.elementsI gpaok in 
			  List.length larci in
			  let lvarc = ref [] in
			  List.iter (fun (a,b,c) -> (lvarc := a::b::(!lvarc)))  larci in
			  let lvarctrie = List.sort (Pervasives.compare) (!lvarc) in
			  let lvarcsansred = remove_redond (Pervasives.compare) lvarctrie in
			  List.length lvarcsansred in
			  580!! cad exactement le nombre de noeuds dans lttrie et meme :
			  lvarcsansred = lttrie in
			  - : bool = true
			  Conclusion : le tri topologique ne prend pas en compte les noeuds non connectés à un autre
			  (travaille sur les edges) -> il suffit donc d'ajouter au tri topo les index des sommets
			  présents mais non connecté car eux meme sont a la fois feuille et racine et constiueront
			  donc des chemins donc des formes alternatives 
			*)
			


  
			let ltotp = Array.to_list (Array.mapi (fun i elt -> i) (Array.create nbokp 0)) 
			and ltotm = Array.to_list (Array.mapi (fun i elt -> i) (Array.create nbokm 0)) in
			
			let liseulp = listmissnb lttriep ltotp and liseulm = listmissnb lttriem ltotm in
			let ltvoulup = List.flatten (ltp::[liseulp]) and ltvoulum = List.flatten (ltm::[liseulm]) in
			  (* 999 sur le brin + et 916 sur le brin - *)

  
			(* à terme mettre ici une fonction qui remplit tchem *)
			let tchemp = Array.create nbokp [] and tchemm = Array.create nbokm [] in
			  
			let up = List.iter (fun idx -> 
					      if ((Mg.Succ.numI gpaok idx)=0) then
						tchemp.(idx)<-[[tlabp.(idx)]]
					      else
						tchemp.(idx) <- spread tlabp.(idx) 
						  (List.flatten
						     (List.map 
							(fun idsucc -> tchemp.(idsucc))
							(Mg.Succ.getI gpaok idx)))) ltvoulup 
			and um = List.iter (fun idx -> 
					      if ((Mg.Succ.numI gmaok idx)=0) then
						tchemm.(idx)<-[[tlabm.(idx)]]
					      else
						tchemm.(idx) <- spread tlabm.(idx) 
						  (List.flatten
						     (List.map 
							(fun idsucc -> tchemm.(idsucc))
							(Mg.Succ.getI gmaok idx)))) ltvoulum in
			  
			let lrootp = List.filter (fun i -> ((Mg.Pred.numI gpaok i) = 0)) ltvoulup 
			and lrootm = List.filter (fun i -> ((Mg.Pred.numI gmaok i) = 0)) ltvoulum in
			  (* 605 racines pour le brin +, et 548 pour le brin - *)
			  
			let lchemp = List.flatten (List.fold_left
						     (fun ltmp idroot -> (tchemp.(idroot))::ltmp) 
						     []
						     lrootp) and
			  lchemm = List.flatten (List.fold_left
						   (fun ltmp idroot -> (tchemm.(idroot))::ltmp) 
						   []
						   lrootm) in
			  (* 605 chemins donc 605 formes alternatives 
			     (dont 250 avec au - 2 arn et 355 avec un seul arn) sur le brin +
			     
			     548 chemins donc 548 formes alternatives 
			     (dont 224 avec au - 2 arn et 324 avec un seul arn) sur le brin - *)

			  (* let lcheminterest = List.filter (fun l -> (List.length l)>=2) lchemp in
			     ou List.length (List.filter (fun l -> (List.length l)>=2) lchemp) in
			     ça a l'air de marcher!! ex no 30 de lcheminterest à 4 ARN est ok 
			     Remarque : à terme on caclulera les chemins dans une fonction du genre give_path 
			     mais il faut d'abord résoudre le pb du tri topologique qui ne prend pas les sommets
			     non connectés *)
			  
			  Common.print_log "I have found the paths in the graph for mrna data\n";
  

			  (*****************************************************************************************
			    CREATION DE TOUTES LES FORMES ALTERNATIVES (TRANSCRITS ALTERNATIFS) A PARTIR DES ARNS
                            = LES FEUILLES ARN 
			  **********************************************************************************************)

			  let llbafp = List.map (fun lv -> List.map vert2bmm lv) lchemp 
			  and llbafm = List.map (fun lv -> List.map vert2bmm lv) lchemm in

  
			  let larnleafp = List.map 
					    (fun lb -> ArnLeaf 
					       (
						 (Bluspmt.create Rna lb (Bluspmm.gbh (List.hd lb)) (Bluspmm.geh (last lb)) Forward), 
						 (List.map 
						    (fun (a,b) -> Real(a,b)) 
						    (remove_redond Pervasives.compare 
						       (List.sort Pervasives.compare (List.map (fun (Real(a,b)) -> (a,b)) (List.filter realint (List.flatten (List.fold_left (fun ltmp il -> il::ltmp) [] (List.map Bluspmm.il lb)))))))  
						 )
					       )
					    ) llbafp 

			  and larnleafm = List.map (fun lb -> ArnLeaf ((Bluspmt.create Rna lb (Bluspmm.gbh (List.hd lb)) (Bluspmm.geh (last lb)) Reverse), List.map (fun (a,b) -> Real(a,b)) (remove_redond Pervasives.compare (List.sort Pervasives.compare (List.map (fun (Real(a,b)) -> (a,b)) (List.filter realint (List.flatten (List.fold_left (fun ltmp il -> il::ltmp) [] (List.map Bluspmm.il lb)))))))  )) llbafm in
			    (* remarque : liste de bluspmt * intron list *)

					     (* List.iter (fun (ArnLeaf (a,il)) ->
							   List.iter (fun  (Real (d,f)) -> 
									Common.print_log ("ideb = "^(string_of_int d)^" - ifin = "^(string_of_int f)^"\n")) il
							)
						larnleafp;*)	

			  let larnleaftriep = List.sort (fun al1 al2 -> 
							   match (al1,al2) with
							       (ArnLeaf (bmt1,_),ArnLeaf (bmt2,_)) -> Bluspmt.compare bmt1 bmt2 
							     |_ -> 0
							) larnleafp 
						
			  and larnleaftriem = List.sort (fun al1 al2 -> 
							   match (al1,al2) with
							       (ArnLeaf (bmt1,_),ArnLeaf (bmt2,_)) -> Bluspmt.compare bmt1 bmt2 
							     |_ -> 0
							) larnleafm in
			    
			  (* arnleaf rangées par génomique croissant, brins séparés *)
			     






			  (*******************************************************************************************
			    
			    CREATION DES RACINES PROTEIQUES FONDEES SUR LE CHEVAUCHEMENT GENOMIQUE DES BLUSPS PROTEIQUES 
			    
			  ********************************************************************************************)


			  let nbp = List.length lbpord and nbm = List.length lbmord in
			    (* nombre de bluspmm protéiques sur chacun des brins *)
			    

			  (* Clustering des bluspmm protéique selon le chevauchement génomique *)
			  let sclbpp = cluspingbpm (SegSeq.single (Array.of_list lbpord)) 
			  and sclbpm = cluspingbpm (SegSeq.single (Array.of_list lbmord)) in
			    
			  let llb2p = List.map (Array.to_list) (SegSeq.elements sclbpp) 
			  and llb2m = List.map (Array.to_list) (SegSeq.elements sclbpm) in
			    (* clusters de bluspmm prot chevauchants sur le génomique
			       = liste de listes de bluspmm protéiques chevauchants sur le génomique 
			    *)
  

			  (* au 13/03/06 : le blusping multimolécule protéique est trop simpliste : il se fait uniquement
			     sur la base du chevauchement génomique, ce qui fait que les bluspmm d'une meme protéine
			     séparés au moment de la vérification du retour en arrière (dacm1 = blusping momolécule protéique),
			     peuvent etre de nouveau raccordés ici.
			     Au lieu de faire un protroot par liste de bluspmm chevauchants on va pouvoir en faire plusieurs
			     en fonction des pb de retour en arrière rencontrés dans les différents hsps des clusps de cet 
			     ensemble de bluspmm protéique
			  *)
			  let lprotrootp = List.sort eltassoccomp (List.flatten
					     (List.map 
						(fun lb -> 
						   let lhsp = List.flatten (List.map Bluspmm.hl lb) in
						   let lhsptrie = List.sort Hsp.compare lhsp in
						     
						   (* On fait les clusps des différents hsps des bluspmm protéiques 
						      Attention ici si un hsp prot du clusp a son champ cutdwn à true
						      cela est répercuté dans le champ cutdwn du clusp 
						   *)
						   let lclusp = List.map thspchev_to_clusp 
								  (SegSeq.elements (clusping (SegSeq.single (Array.of_list lhsptrie)))) in
						     
						   (*  List.iter (print_clusp2 stderr) lclusp; *)
						   (* liste de tableaux de clusps, afin de couper qd on avait repéré un intron avant *)
						   let ltcl = SegSeq.elements 
								(blusping2monomol
								   noncutdwn
								   (SegSeq.make3 (Array.of_list lclusp) [0]) c.tmaxintron c.beta) in
						     
						     List.map 
						       (
							 fun tcl -> ProtRoot ((tclusp2bluspmt tcl),(Array.to_list tcl))
						       )
						       ltcl) 
						
						llb2p)) 

			  and lprotrootm = List.sort eltassoccomp (List.flatten
					     (List.map 
						(fun lb -> 
						   let lhsp = List.flatten (List.map Bluspmm.hl lb) in
						   let lhsptrie = List.sort Hsp.compare lhsp in
						     
						   (* On fait les clusps des différents hsps des bluspmm protéiques 
						      Attention ici si un hsp prot du clusp a son champ cutdwn à true
						      cela est répercuté dans le champ cutdwn du clusp 
						   *)
						   let lclusp = List.map thspchev_to_clusp 
								  (SegSeq.elements (clusping (SegSeq.single (Array.of_list lhsptrie)))) in
						     
						   (*  List.iter (print_clusp2 stderr) lclusp; *)
						   (* liste de tableaux de clusps, afin de couper qd on avait repéré un intron avant *)
						   let ltcl = SegSeq.elements 
								(blusping2monomol
								   noncutdwn
								   (SegSeq.make3 (Array.of_list lclusp) [0]) c.tmaxintron c.beta) in
						     
						     List.map 
						       (
							 fun tcl -> ProtRoot ((tclusp2bluspmt tcl),(Array.to_list tcl))
						       )
						       ltcl) 
						
						llb2m)) in
  
			  (* liste de protroot rangés en fonction de leur début sur le génomique
			     = liste de blusp multi molécule protéiques 
			     Rq : chaque protroot contient à la fois le bluspmt (liste de bluspmm, deb, fin, brin...)
			     mais aussi la liste des clusps correspondants à ce bluspmt
			     
			     Attention : ici on retient les cassures qui ont eu lieu lors du blusping
			     monoprotéique, de sorte que l'on peut avoir plusieurs Bluspmt alors qu'avant
			     on n'en avait qu'un
			  *)
			    
   
			  (* REGROUPEMENT DES ARNLEAF AVEC LES PROTRAC dans un meme objet nommé eltassoc 
			     (attention : cet objet ne stocke pas les memes info pour un protroot ou arnleaf) *)
			  let leltassoctriep = List.merge eltassoccomp lprotrootp larnleaftriep 
			  and leltassoctriem = List.merge eltassoccomp lprotrootm larnleaftriem in
			    (* 1102 eltassoc rangées par génomique croissant sur le brin + 
			       ... sur le brin - *)
 


			  (******************************************************************************
			    POUR CHAQUE FEUILLE ARN TENTATIVE D'ASSOCIATION AVEC LES RACINES PROTEIQUES 
                            DE SON CLUSTER
			    Pour chaque cluster : Liste de couples 1 ARN - 1 PROT
			  ********************************************************************************)


			  let lleassocp = List.map Array.to_list (SegSeq.elements (cluspingeas (SegSeq.single (Array.of_list leltassoctriep)))) 
			  and lleassocm = List.map Array.to_list (SegSeq.elements (cluspingeas (SegSeq.single (Array.of_list leltassoctriem)))) in
			    (* clusters (ou listes de listes) d'elements assoc (cad ProtRoot ou ArnLeaf) 
			       brin séparé *)
  

			  (* Pour chaque lassoc, autrement dit cluster d'eltassoc (prot ou arn) recherche des associations
			     recevables entre arnleaf et protrac *)
			  let llcplapp = List.map 
					   (fun leas -> 
					      let (la,lp) = List.partition isleafarn leas in
						(* couple de listes la premiere d'arn et la deuxieme de prot *) 
						match (la,lp) with
						    (_,[]) -> List.map (fun a -> (a,Nulle)) la
						  |([],_) -> List.map (fun p -> (Nulle,p)) lp
						  |_ -> (* ici on est sur que la et lp sont toutes deux non vides *)
						     let cartprod = List.flatten (List.map (fun a -> List.map (fun p -> (a,p)) lp) la) in
						       (* produit cartésien de ces listes, cad liste de couples (arn,prot) *)
						       remove_redond 
							 compcpleas
							 (List.map (fun (a,p) -> 
								      if ((easoverlap a p) && (not (existclpininta a p))) then
									(a,p) 
								      else
									(a,Nulle)
								   ) cartprod)
							 (* filtrage des couple arn,prot qui ont le pb de clusp prot dans intron arn *)
					   )
					   lleassocp 

			  and llcplapm = List.map 
					   (fun leas -> 
					      let (la,lp) = List.partition isleafarn leas in
						(* couple de listes la premiere d'arn et la deuxieme de prot *) 
						match (la,lp) with
						    (_,[]) -> List.map (fun a -> (a,Nulle)) la
						  |([],_) -> List.map (fun p -> (Nulle,p)) lp
						  |_ -> (* ici on est sur que la et lp sont toutes deux non vides *)
						     let cartprod = List.flatten (List.map (fun a -> List.map (fun p -> (a,p)) lp) la) in
						       (* produit cartésien de ces listes, cad liste de couples (arn,prot) *)
						       remove_redond 
							 compcpleas
							 (List.map (fun (a,p) -> 
								      if ((easoverlap a p) && (not (existclpininta a p))) then
									(a,p) 
								      else
									(a,Nulle)
								   ) cartprod)
							 (* filtrage des couple arn,prot qui ont le pb de clusp prot dans intron arn *)
					   )
					   lleassocm 
			  in
			    
		

			    Common.print_log "I have done the list of pairs 1 mrna - 1 protein\n";
			    
			    (********************************************************************************
			      REGROUPEMENT DES PROTROOT RELIEES DE FACON VALIDE AU MEME ARNLEAF
			      Pour chaque cluster : Liste de couples 1 ARN - plusieures PROT
			    ***************************************************************************************)
			    

			    (* Pour chaque cluster on forme une table de hachage de clé un arnleaf et de valeur 
			       un protroot, avec find_all on récupère tous les protroot liées à un même arnleaf. 
			       llcplapset contient donc pour chaque cluster une liste de couples (a,pset)
			       ou le premier element est un arnleaf et le deuxieme un ensemble de bluspmt proteiques.
			       Remarque importante : on est passé aux BmtSet afin de pouvoir faire des intersections
			       de Bluspmt protéique afin de pouvoir relier deux arn non chevauchant ayant l'intersection
			       de leur BmtSet protéique non vide! *)
			    let llcplaplp = List.map 
					      (fun lcplap -> 
						 let n = List.length lcplap in
						 let h = Hashtbl.create n in
						 let lcplapl = ref [] in
						   (* pour chaque liste de couple (arn,prot) on met à la bonne place dans la table h *)
						   List.iter 
						     (fun (a,p) -> 
							match p with
							    ProtRoot(b,_) -> Hashtbl.add h a b
							  |Nulle -> Hashtbl.add h a (Bluspmt.nullb)
						     ) 
						     lcplap;
						   
						   (* pour chaque arn de la table h on remplit la liste de couples (arn,pset) lcplapset *)
						   Hashtbl.iter 
						     (fun a p -> 
							lcplapl := (a,(List.sort Bluspmt.compare (Hashtbl.find_all h a)))::(!lcplapl)
						     ) 
						     h;
						   (* on renvoit la liste de couples (arnleaf, protbmtset) triée par ordre de génomique
						      croissant sur les arnleaf et sans redondance *)
						   remove_redond (fun cpl1 cpl2 -> if(cpl1=cpl2) then 0 else 1) (List.sort (fun (a1,b1) (a2,b2) -> eltassoccomp a1 a2) (!lcplapl))
					      )
					      llcplapp 

			    and llcplaplm = List.map 
					      (fun lcplap -> 
						 let n = List.length lcplap in
						 let h = Hashtbl.create n in
						 let lcplapl = ref [] in
						   (* pour chaque liste de couple (arn,prot) on met à la bonne place dans la table h *)
						   List.iter 
						     (fun (a,p) -> 
							match p with
							    ProtRoot(b,_) -> Hashtbl.add h a b
							  |Nulle -> Hashtbl.add h a (Bluspmt.nullb)
						     ) 
						     lcplap;
						   
						   (* pour chaque arn de la table h on remplit la liste de couples (arn,pset) lcplapset *)
						   Hashtbl.iter 
						     (fun a p -> 
							lcplapl := (a,(List.sort Bluspmt.compare (Hashtbl.find_all h a)))::(!lcplapl)
						     ) 
						     h;
						   (* on renvoit la liste de couples (arnleaf, protbmtset) triée par ordre de génomique
						      croissant sur les arnleaf et sans redondance *)
						   remove_redond (fun cpl1 cpl2 -> if(cpl1=cpl2) then 0 else 1) (List.sort (fun (a1,b1) (a2,b2) -> eltassoccomp a1 a2) (!lcplapl))
					      )
					      llcplapm in
			      
			      Common.print_log "I have done the list of pairs 1 mrna - several proteins\n";

			      
			      (********************************************************************************
				REGROUPEMENT DES PROTROOT RELIEES DE FACON VALIDE AU MEME ARNLEAF
				Pour chaque cluster : Liste de couples plusieurs ARN - plusieures PROT
			      ***************************************************************************************)
			      

			      (* on sait que l'on a un seul arn dans chaque premiere partie du couple cplapset 
				 Remarque importante : on a simplifié le traitement dans le sens ou il faudrait 
				 rendre possible de rejoindre plusieurs arns non chevauchants par un intron virtuel
				 des lors qu'une protéine les relie (en fait un bluspmm protéique contenu dans 
				 le pset vu au-dessus 
				 Ici on passe artificiellement de l'arn à une liste d'arn (ici ne contenant que lui-meme)
				 afin de pouvoir faire le traitement précédant (regrouper deux arns par une prot) *)
			      let llcplalplp = List.map 
						 (fun lcplapl -> List.map 
						    (fun (xa,pl) -> match xa with
						       |ArnLeaf (a,_) -> ([a], pl)
						       |Nulle -> ([],pl))
						    lcplapl)
						 llcplaplp 
			      and llcplalplm = List.map 
						 (fun lcplapl -> List.map 
						    (fun (xa,pl) -> match xa with
						       |ArnLeaf (a,_) -> ([a], pl)
						       |Nulle -> ([],pl))
						    lcplapl)
						 llcplaplm in		



			      (* Ajout au 15/03/06, pour afficher quelles molécules entrent dans la composition 
				 des différents transcrits exogean
				 seul pb : nous n'avons pas ici les noms des transcrits 
			      *)
			      let (f, ext, s) = Printing.get_type_sortie "mysequence" context.outfile in
			      let ocomp = open_out ((s)^(".compmolec.txt")) in
				List.iter (fun lcpl -> List.iter (fun cpl -> Printing.print_cpl_ArnBmtlist_ProtBmtlist ocomp cpl) lcpl) llcplalplp;
				
				List.iter (fun lcpl -> List.iter (fun cpl -> Printing.print_cpl_ArnBmtlist_ProtBmtlist ocomp cpl) lcpl) llcplalplm;
				
				flush ocomp;
				close_out ocomp;				
				

			      let lllhap = List.map 
					     (fun lcplalpl -> 
						List.map
						(fun (al,pl) -> 
						   List.merge 
						   (Hsp.compare) 
						   (List.sort (Hsp.compare) (List.flatten (List.map (Bluspmm.hl) (List.flatten (List.map (Bluspmt.blist) al))))) 
						   (List.sort (Hsp.compare) (List.flatten (List.map (Bluspmm.hl) (List.flatten (List.map (Bluspmt.blist) pl)))))
						)
						lcplalpl 
					     )
					     llcplalplp 
  
			      and lllham = List.map 
					     (fun lcplalpl -> 
						List.map
						(fun (al,pl) -> 
						   List.merge 
						   (Hsp.compare) 
						   (List.sort (Hsp.compare) (List.flatten (List.map (Bluspmm.hl) (List.flatten (List.map (Bluspmt.blist) al))))) 
						   (List.sort (Hsp.compare) (List.flatten (List.map (Bluspmm.hl) (List.flatten (List.map (Bluspmt.blist) pl)))))
						)
						lcplalpl 
					     )
					     llcplalplm in
				
				
				Common.print_log "I have done the list of pairs several mrna - several proteins\n";
				
				(*****************************************************************************************************
				  FORMATION DES CLUSPS D'HSPS MIXTES PROT+ARN
				******************************************************************************************************)
				
				(* Liste de liste de segseq clusps  : pour chaque cluster et au sein d'un cluster 
				   pour chaque forme alternative on a une liste de segseq de clusps 
				*)
				let llscap = List.map 
					       (fun llhap -> 
						  List.map
						  (fun lhap -> 
						     SegSeq.make3 (Array.of_list (List.map thspchev_to_clusp (SegSeq.elements (clusping (SegSeq.single (Array.of_list lhap)))))) [0] )
						  (List.filter (fun lhap -> lhap!=[]) llhap)
					       )
					       lllhap and
				  
				  llscam = List.map 
					     (fun llhap -> 
						List.map
						(fun lhap -> 
						   SegSeq.make3 (Array.of_list (List.map thspchev_to_clusp (SegSeq.elements (clusping (SegSeq.single (Array.of_list lhap)))))) [0] )
						(List.filter (fun lhap -> lhap!=[]) llhap)
					     )
					     lllham in
						
				  Common.print_log "I have made the clusps\n";

				  (*
				    Array.iteri
				    (fun i cl -> 
				       Common.print_log 
				       (("Clusp ")^(string_of_int i)^(" : deb = ")^(string_of_int (Clusp.gbeg cl))^(" ; fin = ")^(string_of_int (Clusp.gend cl))^("\n"))
				    )
				    (SegSeq.tank (List.hd (List.hd llscam)));

				  *)	  
				  (**********************************************************************************************
                                    RECHERCHE DE SIGNAUX 
				  **********************************************************************************************)  
				  
				  (* Les eclusps : prend un peu de temps *)
				  let lltecp = List.map 
						 (fun lscap -> 
						    List.map 
						    (fun scap -> sClusp_to_tEclusp c.nninhsp c.nnsuppl scap gseq bqp c.nbaamanq)
						    lscap
						 ) 
						 llscap 
				    
				  and lltecm = List.map 
						 (fun lscap -> 
						    List.map 
						    (fun scap -> sClusp_to_tEclusp c.nninhsp c.nnsuppl scap compgseq bqp c.nbaamanq)
						    lscap
						 ) 
						 llscam in

				  
				    Common.print_log "I have made the eclusps\n";
				    
				    (* 
				       Array.iteri 
				       (fun i ec ->  
				       Common.print_log 
				       (("Eclusp ")^(string_of_int i)^(" : deb = ")^(string_of_int (Eclusp.gbeg ec))^(" ; fin = ")^(string_of_int (Eclusp.gend ec))^("\n"))
				       )
				       (List.hd (List.hd lltecm));
				   *)



				    (* Les exons : prend un peu de temps *)
				    let lltexp = List.map 
						   (fun ltec -> 
						      List.map 
						      (fun tec -> 
							 Array.mapi 
							 (ec_to_ex c.nbaamanq tec c.tminintron c.waysigsearch prob_seq gseq mAG mGT)
							 tec)
						      ltec
						   )
						   lltecp in

				      Common.print_log "I have made the + strand exons\n";
				      
				      let lltexm = List.map 
						     (fun ltec -> 
							List.map 
							(fun tec -> 
							   Array.mapi 
							   (ec_to_ex c.nbaamanq tec c.tminintron c.waysigsearch prob_seq compgseq mAG mGT)
							   tec)
							ltec
						     )
						     lltecm in
					
					Common.print_log "I have made the exons\n";
					
					
					(* Les gènes = transcrits sans l'indication de leur CDS
					   Chose étrange : il semble qu'il y ait ici des gènes =tr
					   avec arnpal=N et protpale =null sans etre nuls eux meme,
					   rq : vient quand pas de prot *)
					let lgenep = List.filter 
						       (fun g -> (((Gene.arnpal g)!= Donnees_base.N) || ((Gene.protpale g)!= Donnees_base.N))) 
						       (List.flatten 
							  (List.map 
							     (fun ltex -> 
								List.map 
								(fun tex -> 
								   texon2onegene tex (0,(Array.length tex)))
								ltex)
							     lltexp)
						       ) 
					  
					and lgenem = List.filter 
						       (fun g -> (((Gene.arnpal g)!= Donnees_base.N) || ((Gene.protpale g)!= Donnees_base.N))) 
						       (List.flatten 
							  (List.map 
							     (fun ltex -> 
								List.map 
								(fun tex -> 
								   texon2onegene tex (0,(Array.length tex)))
								ltex)
							     lltexm)
						       ) in
					  
					
					  (*Array.iteri 
					    (fun i ex -> 
					       Common.print_log 
					       (("Exon ")^(string_of_int i)^(" : deb = ")^(string_of_int (Exon.gbeg ex))^(" ; fin = ")^(string_of_int (Exon.gend ex))^("\n"))
					    )
					    (Gene.texon (List.hd lgenep));
					  *)

					(* les gènes (= transcrits en réalité), AVEC l'indication de leur CDS 
					   Rq : si on veut faire une recherche d'orf en fonction des protéines présentes on met le 
					   premier parametre de gene2genewithcds à true sinon on le laisse à false *)
					let pprot g = true and comparecds = CDS.comparesize in
					
					let lgenewithcdsp = List.map (gene2genewithCDS_mieux comparecds gseq) (List.filter (fun g -> not (Gene.isnull g)) lgenep) and lgenewithcdsm = List.map (gene2genewithCDS_mieux comparecds compgseq) (List.filter (fun g -> not (Gene.isnull g)) lgenem) in
					let lgenewithcds = List.flatten (lgenewithcdsp::[lgenewithcdsm]) in
					  
					let lcompletegene = List.map (gene2completegene gseq compgseq) lgenewithcds in
					  (* Le couple de listes (nonpseudogenes, pseudogènes) des completegene (cf donnees.ml) 
					     Attention : ici il faut mettre les gènes du brin - avec comme reference le début du brin + 
					     et non comme ici le début du brin - *)
					  (((List.partition (nonpseudogene c.szcdsnomet c.szcdswithmet c.propcdsmonoex c.sizecdsmonoex c.sizecdsdiex) lcompletegene),lggseq),tgenearnpret,nomgseq);;




(*
  Pour le pb de la recherche de signaux en s'aidant des protéines 
  un gène qui est juste sans les prot (ou avec les prot toutes de meme phase) et
  qui est faux en prenant seulement la phase de la premiere hsp proteique :
  En face Vega .492 :
  Arnppal = A "mRNA_GB:CR456366.5084"
  ppale =  P "IPI:IPI00133626.1"
  
  let g = List.hd (List.filter (fun g -> (Gene.arnpal g)= (A "mRNA_GB:CR456366.5084")) lgene) in
  
  
*)














