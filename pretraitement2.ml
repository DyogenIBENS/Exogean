(* pretraitement2.ml : réalise les operations de fusion et de remplacement de deux hsps contigues de meme proteine
   et de meme brin. Succede au pretraitement1 qui lui rognait eventuellement certaines hsps aux extremités *)

open Common
open Seq
open Preserveorf
open Donnees_base


let brin_to_string br =
  match br with
    |Forward -> "+";
    |Reverse -> "-";;




(********************************************************************************************************************
             Fonctions pour la lecture des hsps dans le fichier d'entree (resultats Exofish) (apres Db.read)
*********************************************************************************************************************)

let fstof3 tr =
match tr with (a,b,c) -> a;;

let sndof3 tr =
match tr with (a,b,c) -> b;;

let thdof3 tr =
match tr with (a,b,c) -> c;;






(********************************************************************************************************************
             Fonctions pour le prétraitement (essentiellement : Fusion, remplacement, ajustement et etiquetage alignement supplémentaire)
*********************************************************************************************************************)


(* ATTENTION : h1 et h2 sont ici supposés de meme brin car le pretraitement s'effectue brins séparés 
   k est le seuil en nucléotides pour la fusion des hsps de meme proteines proches de k nts sur le génomique 
   et assez proches sur la protéine *)

let regroupep cpa daa k h1 h2 =
  let fusremp = ref 0 in  (* fusremp est à 0 s'il ne faut ni fusionner ni remplacer h1 et h2, à 1 s'il faut fusionner h1 et h2 et à 2 s'il faut remplacer h1 et h2 par l'une des deux *) 
  let pquoi = ref "" in
  let nt = (Hsp.debg h2) - (Hsp.fing h1) -1 in  (* représente le nombre de nucléotides entre la fin du mot génomique de h1 et le debut du mot génomique de h2 *)
  let aa = (Hsp.debp h2) - (Hsp.finp h1) -1 in  (* représente le nombre d'acides aminés entre la fin du mot protéique de h1 et le debut du mot protéique de h2 *)
   (* Common.print_log ((string_of_int (Hsp.fing h1))^"\t"^(string_of_int (Hsp.debg h2))^" ; "^(string_of_int (Hsp.finp h1))^"\t"^(string_of_int (Hsp.debp h2))^"\nnt = "^(string_of_int nt)^" et aa = "^(string_of_int aa)^"\n"); *)
    if (nt<3*aa && nt<100) then (* rajouté au 16/02/06 pour Diatome mais il se peut que ne marche pas
				   bien pour eucaryotes supérieurs, veut-on nt>0 ici? hrc? *)
      begin
	fusremp := 1;
	pquoi := "***Fusion des deux hsps (nt<3*aa && nt<=100) due à une région de la protéine orthologue totalement absente dans la protéine annotée***\n"
      end
    else
      begin
	if(aa<=(-1)) then
	  begin
	    if(nt<=(-1)) then 
              begin
		fusremp := 1;
		pquoi := "***Fusion des deux hsps (aa<=-1 et nt<=-1) due à des répétitions dans la protéine que l'on annote et dans son orthologue***\n"
              end
	  end  (* fin de aa<=-1 *)
	else  
	  begin
	    if(aa=0) then
	      begin
		if(nt>=0 && nt <= 11) then 
		  begin
		    fusremp := 1;
		    pquoi := "***Fusion des deux hsps (aa=0 et 0<=nt<=11), due à insertion d'acides aminés dans la protéine que l'on annote"^(if ((nt mod 3) <> 0) then " + un frameshift***\n" else "***")
		  end
		else
		  begin
		    if(nt<=(-2)) then
		      begin
			fusremp := 2;
			pquoi := "***Remplacement de deux hsps (aa=0 et nt<=-2) par celle qui apparait la premiere dans leur protéine commune due à une répétition dans la séquence à annoter***\n"
		      end
		  end  (* fin de nt<=-2 (dans aa=-1) *)
	      end   (* fin de aa=-1 *)
	    else  (* aa>0 *)
	      begin 
		if(nt>0) then
		  begin 
		    if(nt=3*aa) then 
		      begin
			fusremp := 1;
			pquoi := "***Fusion des deux hsps (nt=3*aa), due à une région peu conservée entre ces deux hsps dans la protéine que l'on annote***\n"
		      end
		    else (* nt!= 3*aa *)
		      begin
			(* dans toutes les fusions qui ont lieu par la suite seul le commentaire 
			   change mais meme principe, sauf pour le dernier cas de fusion *)
			if(((nt-3*aa) mod 3 = 0) && ((nt-3*aa)>=(-3)) && ((nt-3*aa)<=k+9)) then (* avant k-3 
												   mais on souhaite
												   etre plus laxiste
												   si x est modulo 3 *)
			  begin
			    fusremp := 1;
			    pquoi := "***Fusion des deux hsps ((nt-3*aa) mod 3 == 0 && (nt-3*aa)>=-3 && (nt-3*aa)<=k-3), due à une région peu conservée dans un exon et qq acides aminés supplémentaires entre ces deux hsps dans la protéine que l'on annote***\n"
			  end 
			else 
			  begin
			    if((nt-3*aa) mod 3 <> 0 && (nt-3*aa)>=(-2)  && (nt-3*aa)<=5) then
			      begin
				fusremp := 1;
				pquoi := "***Fusion des deux hsps ((nt-3*aa) mod 3 <> 0 && (nt-3*aa)>=(-2)  && (nt-3*aa)<=5), due à un frameshift, avec une perte ou une insertion d'acide aminé dans la protéine que l'on annote***\n"
			      end  
			    else
			      begin
				if((nt-3*aa) mod 3 <> 0 && (nt-3*aa)>=6  && (nt-3*aa)<=k-1) then (* avant k=15 *)
				  begin
				    fusremp := 1;
				    pquoi := "***Fusion des deux hsps ((nt-3*aa) mod 3 <> 0 && (nt-3*aa)>=6  && (nt-3*aa)<=k-1), due à une région peu conservée au sein de l'exon courant, un frameshift et une insertion d'acides aminés dans la protéine que l'on annote***\n"
				  end
				else
				  begin
				    if((nt-3*aa)<=5) then (* avant -4 au lieu de 5 *)
				      begin
					(* Remarque : ici il faudra peut-etre mettre une borne inférieure pour nt-3*aa! *)
					fusremp := 1;
					pquoi := "***Fusion des deux hsps ((nt-3*aa) mod 3 <> 0 && (nt-3*aa)<=(-4)), due à une région peu conservée au sein de l'exon courant, "^(if ((nt-3*aa) mod 3 <> 0) then ("un frameshift ") else "")^"et une perte d'acides aminés dans la protéine que l'on annote***\n"
				      end
				    else
				      begin
					if (nt-3*aa >= cpa && aa>=daa) then
					  begin
					    fusremp := 0;
					    pquoi:= "***nécessité d'un alignement complémentaire avec son hsp suivante de meme protéine, aa>=0,nt>=0,a>=5***\n"
					  end
				      end
				  end
			      end
			  end
		      end
		  end (* fin de nt>0 dans aa>=0 *)
		else 
		  begin
		    if(nt=0) then
		      begin
			fusremp := 1;
			pquoi := "***Fusion des deux hsps due à une juxtaposition des deux hsps sur le génomique et une perte d'acides aminés dans la protéine orthologue (aa>0 et nt=0)***\n"
		      end
		    else
		      begin
			if(nt=(-1)) then
			  begin
			    fusremp := 1;
			    pquoi := "***Fusion des deux hsps (aa>0 et nt=-1), due à une délétion d'un nombre de nucléotides multiple de 3 dans la protéine que l'on annote***\n"
			  end
			else (* nt<=-2 *)
			  begin
			    fusremp := 2;
			    pquoi := "***Remplacement de deux hsps (aa>=0 et nt<=-2) par celle qui apparait la premiere dans leur protéine commune due à une répétition***\n"
			  end
		      end
		  end
	      end
	  end
      end;
    (*Common.print_log ("(Hsp.fing h1)="^(string_of_int(Hsp.fing h1))^" (Hsp.debg h2) = "^(string_of_int(Hsp.debg h2))^"\n");
    Common.print_log ("(Hsp.finp h1)="^(string_of_int(Hsp.finp h1))^" (Hsp.debp h2) = "^(string_of_int(Hsp.debp h2))^"\n");
    Common.print_log (!pquoi);*)
    (!fusremp,!pquoi) ;;



(* fonction qui détermine si h1 et h2 hsps arn sont à fusionner dans le cas d'un intron petit 
   de taille inferieure à 3 
   Remarque : cpa et daa ne servent à rien mais il fallait que regroupea ait la meme signature
   que regroupep
   || (((nt mod 3) = 0) && (nt <= szsmintr))  *)
let regroupea cpa daa szsmintr h1 h2 =
  let fusremp = ref 0 in  (* fusremp est à 0 s'il ne faut ni fusionner ni remplacer h1 et h2, 
			     à 1 s'il faut fusionner h1 et h2 
			     et à 2 s'il faut remplacer h1 et h2 par l'une des deux *) 
  let pquoi = ref "" in
  let nt = (Hsp.debg h2) - (Hsp.fing h1) -1 in  (* représente le nombre de nucléotides entre la fin 
						   du mot génomique de h1 et le debut du mot génomique de h2 *)
    if (nt<=3) then
      begin
	fusremp := 1;
	pquoi := "***Fusion des deux hsps arn due à un intron de taille tres petite (1,2, ou 3 nt)***"
      end;
    (!fusremp,!pquoi) ;;	





(* Fonction auxiliaire à regroupe_hsp destinée à fusionner deux hsps h1 et h2, mais attention : 
   seulement dans le cas où cett fusion n'engendre pas l'inclusion de stops pour les protéines 
   Cette fonction renvoit une liste d'hsps, idealement contenant une seule HSP*)
let fusion gseq h1 h2 s =
(*  Common.print_log "Je vais faire une fusion d'hsps\n"; *)
  Hsp.create (Hsp.tmol h1) (Hsp.nomp h1) ((Hsp.score h1) +. (Hsp.score h2)) (min (Hsp.debg h1) (Hsp.debg h2)) (max (Hsp.fing h1) (Hsp.fing h2)) (min (Hsp.debp h1) (Hsp.debp h2)) (max (Hsp.finp h1) (Hsp.finp h2)) (Hsp.brinh h1) ((Hsp.nbntfus h1)+abs ((Hsp.debg h2)-(Hsp.fing h1)-1)) false false ((Hsp.commh h1)^(Hsp.commh h2)^s);;

    

    

(* Fonction auxiliaire à regroupe_hsp destinée à remplacer h1,h2 par celle qui apparait la première dans la protéine. 
   Attention : dans cette fonction comme dans la précédente les hsps sont supposées provenir de la
   meme protéine et du meme brin génomique, on peut eventuellement le tester par une exception *)
let remplacement h1 h2 s =
  let h = if((Hsp.debp h1) <= (Hsp.debp h2)) then h1 else h2 in
    Hsp.create (Hsp.tmol h1) (Hsp.nomp h) (Hsp.score h) (Hsp.debg h) (Hsp.fing h) (Hsp.debp h) (Hsp.finp h) (Hsp.brinh h) (Hsp.nbntfus h1) false false ((Hsp.commh h1)^(Hsp.commh h2)^s);;


(* 
   Dans lhspfinale on trouvera la liste d'hsps apres une passe de fusion/remplacement sur la liste lhsps
   dans change on trouvera 1 ou 2 si un changement a été nécessaire. Remarque : les hsps de lhsps sont
   supposé provenir de la meme protéine. 
   regroupe est la fonction de regroupement d'hsps :
   - regroupep pour les hsps proteiques 
   - regroupea pour les hsps arns

   Attention lhspfinale est passée initialement vide et change égal à 0
*)
let rec regroupe_hsp tmol gseq regroupe cpa daa k lhsps lhspfinale change bqp =  
  match lhsps with
   |[] -> (!lhspfinale,0);
   |[hsp] -> if (not (List.mem (Hsp.debg hsp) (List.map Hsp.debg (!lhspfinale)))) then
       lhspfinale:= List.append !lhspfinale [hsp]; 
       (!lhspfinale,0);(*if(!change<=2) then !change else 0*)
   |h1::h2::q -> match (regroupe cpa daa k h1 h2) with
 
       |(1,s) ->  (* cas de fusion nécessaire *)
	  if (tmol=Rna || (pasStop (true,(Seq.sub gseq (Hsp.debg h1) ((Hsp.fing h2)-(Hsp.debg h1)+1))))) then
	    begin
	      let hfus=fusion gseq h1 h2 s in
		(* on ajoute la condition de si cet hsp n'est pas deja dans la liste 
		   mais on n'évite pas la redondance (h12 et h123), peut-etre vaudrait-il
		   mieux tester si le début de hfus est le meme qu'un début d'un h dans la liste? 
		*)
		if (not (List.mem (Hsp.debg hfus) (List.map Hsp.debg (!lhspfinale)))) then
		  lhspfinale:=  List.append (!lhspfinale) [hfus];
		change := 1; 
		regroupe_hsp tmol gseq regroupe cpa daa k (hfus::q) lhspfinale change bqp (* avant q *)
	    end
	  else
	    begin
	      if (not (List.mem (Hsp.debg h1) (List.map Hsp.debg (!lhspfinale)))) then
		lhspfinale:= List.append !lhspfinale [h1];
	      if (not (List.mem (Hsp.debg h2) (List.map Hsp.debg (!lhspfinale)))) then
		lhspfinale:= List.append !lhspfinale [h2];
	      change :=0;
	      regroupe_hsp tmol gseq regroupe cpa daa k q lhspfinale change bqp
	    end
	   (* Attention ici il faut faire attention de garder des hsps de mot génomique multiple de 3,
	      et ceci meme si on n'a plus g=3*p *)

       |(2,s) -> (* cas de remplacement nécessaire *)
	  let hremp = remplacement h1 h2 s in
            begin
	      if (not (List.mem (Hsp.debg hremp) (List.map Hsp.debg (!lhspfinale)))) then
                lhspfinale:= List.append !lhspfinale [hremp]; 
              change := 2;
	      regroupe_hsp tmol gseq regroupe cpa daa k (hremp::q) lhspfinale change bqp (* avant q *)
            end  
	   (* remplacement de h1,h2 par h1 ou h2 (celle qui apparait la première dans la protéine), 
	      iteration sur hremp::q *)

       |(0,s) -> (* ni fusion ni remplacement nécessaires *)
	  if (not (List.mem (Hsp.debg h1) (List.map Hsp.debg (!lhspfinale)))) then
	    lhspfinale:= List.append !lhspfinale [h1]; 
	   regroupe_hsp tmol gseq regroupe cpa daa k (h2::q) lhspfinale change bqp ;;  




(* 
   La fonction fusremp prend une liste d'hsps et remplit une liste d'hsps en fusionnant remplacant. 
   Tant que des fusions ou des remplacements sont possibles, on itere la fonction de
   regroupement regroupe_hsp sur la liste lhsps passée en entrée 
   Retourne une liste d'hsps.
*)
let fusremp tmol gseq regroupe cpa daa k bqp lhsps =
  let resregroupe = ref (regroupe_hsp tmol gseq regroupe cpa daa k lhsps (ref []) (ref 0)  bqp) in 
    while ((snd !resregroupe) <> 0) do
      resregroupe := regroupe_hsp tmol gseq regroupe cpa daa k (fst !resregroupe) (ref []) (ref 0)  bqp;
    done;
    fst !resregroupe;;





(********************************************************************************************************************
             Fonctions pour le rajout d'hsps par un tblastn relaxé (par rapport au tblastn initial)
*********************************************************************************************************************)



(* 
   print_into_mask permet d'ecrire dans le fichier masqué 3 lignes explicitant la requete de recherche d'hsps 
   dans la zone bien delimitée entre deux hsps de fin et debut sur génomique, protéique : e1g b2g e1p b2p
*)
let print_into_mask o idx e1g b2g e1p b2p =
  let print_une_phase out id eg bg ep bp phase = 
    Printf.fprintf out "1 %i %i %i %i 7 %i %i\n" phase (eg/3+1) (bg/3-1) id (ep+1) (bp-1) in
    
    (* +1 et -1 resp. pour dire que l'on veut exlure les bornes actuelles de la recherche *)
    List.iter (print_une_phase o idx e1g b2g e1p b2p) [0;1;2];;




(* Pour ajouter des hsps avant le début des matches d'une protéine *)
let remind_align_lhsp_bef o aa ip lh =
  try
    let h = List.hd lh in
    let nbaamq = Hsp.debp h in
      if (nbaamq >= aa) then 
	print_into_mask o ip ((Hsp.debg h)+1-2000) ((Hsp.debg h)+1) ((Hsp.debp h)+1-nbaamq) ((Hsp.debp h)+1);
      (* rq : les +1 pour revenir à la reference lassap blast *)
      (* rq : il faudra peut-etre mettre un nombre fonction de naamq au lieu de 2000 *)
  with
    |Failure _ -> ();;(* si on est dans une liste vide car pas d'hsps correspondant à cette protéine *)



(* Pour ajouter des hsps entre deux hsps contigues d'une meme protéine *)
let rec remind_align_lhsp_bet o ip lh =
  match lh with
    |h1::(h2::q as lq) when (Hsp.aligncomp h1) -> print_into_mask o ip ((Hsp.fing h1)+1) ((Hsp.debg h2)+1) 
       ((Hsp.finp h1)+1) ((Hsp.debp h2)+1); remind_align_lhsp_bet o ip lq; 
	(* rq : les +1 pour revenir à la reference lassap blast *)
    |_ -> ();;
	  


(* Pour ajouter des hsps après la fin des matches d'une protéine *)
let rec dernier l = 
  match l with
    |[] -> failwith "liste vide";
    |[x] -> x;
    |t::q -> dernier q;;

(* tlgp est le tableau des longueurs des protéines, ip est l'index de la protéine dans le fichier de protéines initial *)
let remind_align_lhsp_aft o aa tlgp ip lh =
  try
    let h = dernier lh and lgp = tlgp.(ip) in
      if (lgp > ((Hsp.finp h)+1)+aa) then
	print_into_mask o ip ((Hsp.fing h)+1) ((Hsp.fing h)+1+2000) ((Hsp.finp h)+1) lgp;
      (* rq : il faudra peut-etre mettre un nombre fonction de naamq au lieu de 2000 *)
  with
    |Failure _ -> ();;



let remind_align_lhsp o aa tlp ip lh =
  remind_align_lhsp_bef o aa ip lh;
  remind_align_lhsp_bet o ip lh;
  remind_align_lhsp_aft o aa tlp ip lh;;






