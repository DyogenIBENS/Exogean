(* orfsearch.ml *)

open Donnees_base
open Donnees
open Seq
open TraiteSignal
open Preserveorf
open Common
open Config


(* gene2gseq prend en entr�e un g�ne (ou du moins une succession d'exons) et la sequence g�nomique
   sur lequel il se trouve, et renvoit la sous-sequence g�nomique correspondant � la concat�nation 
   des exons de ce g�ne 
   Attention : il y a une difference fondamentale avec gene2prot de premain.ml, qui a l'air d'etre
   gene2seq suivi de translate, car en fait ici les exons ont leurs positions sur la s�quence
   du brin sur lequel ils sont, et non sur le brin Forward *)
let gene2gseq seq g =
  try
    List.flatten (List.map2 (fun b e -> Seq.to_list (Seq.sub seq b (e-b+1))) (Array.to_list (Gene.gbegex g)) (Array.to_list (Gene.gendex g)))
  with
    |Invalid_argument s -> [];;
   



(* dans le cas ou une s�quence a sa plus longue orf entre son dernier stop et sa fin, 
   comme la distance qui s�pare ce dernier stop (voire le d�but de la phase)
   et la fin de la s�quence n'est pas n�cessairement multiple de 3 on soustrait ce qui est en trop. *)
let adjust2mod3 n = n-(n mod 3);;


  
(* size_window_search prend une taille de fenetre de recherche de signal (tel stop ou atg)
   et la phase dans laquelle se fait la recherche, cad 0, 1 ou 2, et renvoit la taille maximale
   de la sous-fenetre modulo 3 dans laquelle cette recherche peut se faire (en vue de la passer
   � sig_search_down au sein de la fonction gene2genewithcds.
   Remarque : la taille mod 3 en question permet de savoir quelles sont les positions de debut 
   et de fin dans la sequence initiale : debut = orf, fin = (retour fonction)-orf+1
   Pb : la valeur de retour peut etre n�gative 

   Attention : par rapport aux autres versions on a interverti les param�tres
   et enlev� le mot stop du nom de cette fonction
*)
let size_window_search frame size_seq_search =
  if ((size_seq_search mod 3)==0) then
    begin
      if(frame==0) then 
	begin
	  size_seq_search
	end
      else
	begin
	  if (frame==1) then
	    begin
	      size_seq_search-3
	    end
	  else
	    begin
	      size_seq_search-3
	    end
	end
    end
  else
    if ((size_seq_search mod 3)==1) then
      begin
	if(frame==0) then 
	  begin
	    size_seq_search-1
	  end
	else
	  begin
	    if (frame==1) then
	      begin
		size_seq_search-1
	      end
	    else
	      begin
		size_seq_search-4
	      end
	  end
      end
    else
      if ((size_seq_search mod 3)==2) then
	begin
	  if(frame==0) then
	    begin
	      size_seq_search-2
	    end
	  else
	    begin
	      if (frame==1) then
		begin
		  size_seq_search-2
		end
	      else
		begin
		  size_seq_search-2
		end
	    end
	end
      else
	begin
	  raise (Invalid_argument "la phase doit etre 0, 1 ou 2\n"); 
	  0
	end;;
 


(* fonction qui prend une liste d'au moins 2 positions tri�s par ordre croissant
   et renvoit la liste des couples de positions adjacentes tri�e par ordre croissant
   Rq : ceci est une fonction tr�s g�n�rale, peut-etre � mettre dans un autre .ml

   Attention : on fait la supposition que lpos est une liste de positions tri�es
 *)
let listpos2listcplpos lpos =
  let i = ref 0 and n = List.length lpos in
  let tpos = Array.of_list lpos and lcplpos = ref [] in

    if(n<=1) then
      begin
	raise (Invalid_argument "listpos2listcplpos"); 
	[]
      end
    else (* au moins deux pos *)
      begin
	while (!i<n-1) do
	  lcplpos := (tpos.(!i),tpos.(!i+1))::(!lcplpos); 
	  incr i;
	done;
	List.rev (!lcplpos)
      end;;
      


(* fonction qui �tant donn� une liste de positions de stop et une liste de positions
   de met trouv�s dans la meme phase d'une miniseq, toutes deux suppos�es non vides,
   renvoit la liste des CDS de type MetStop initiaux, cad met entre le debut de la miniseq
   et le premier stop

   Attention : on suppose que les stop et met sont dans la MEME phase
   si cela n'est pas le cas alors on va trouver des CDS non mod 3 aberrants

   Rq : on a ce genre de CDS MetStop initial que lorsque la premi�re met de lposmet
   est avant le premier stop de lpostop, on aura autant de tels CDS que l'on a de posmet
   avant lpostop.(0)
*)
let listCDSMetStopInitial lpostop lposmet =
  if (lpostop=[] || lposmet=[]) then
    begin
      raise (Invalid_argument "listCDSMetStopInitial"); 
      []
    end
  else
    begin
      let posfstop = List.nth lpostop 0 in
      let lposmetok =List.filter (fun pmet -> pmet < posfstop) lposmet in
	List.map (fun pmet -> CDS.create MetStop pmet (posfstop-1)) lposmetok
    end;;
    


(* fonction qui �tant donn� une liste de positions de stop et une liste de positions
   de met trouv�s dans la meme phase d'une miniseq, toutes deux suppos�es non vides, 
   renvoit tous les CDS possibles de type MetStop internes (cad entre met entre deux stops)

   Attention : on suppose que les stop et met sont dans la MEME phase
   si cela n'est pas le cas alors on va trouver des CDS non mod 3 aberrants
*)
let listCDSMetStopInternal lpostop lposmet =
  if (lpostop=[] || lposmet=[]) then
    begin
      raise (Invalid_argument "listCDSMetStopInternal"); 
      []
    end
  else
    begin
      if ((List.length lpostop)=1) then
	[]
      else
	begin
	  let lcplpostop = listpos2listcplpos lpostop in
	    List.flatten 
	      (List.map
		 (fun (pos1,pos2) -> 
		    (* on ne garde que les met situ�es entre deux stops contigus *)
		    let lposmetok = List.filter (fun pmet -> (pos1 < pmet) && (pmet < pos2)) lposmet in
		      List.map (fun pmet -> CDS.create MetStop pmet (pos2-1)) lposmetok
		 )
		 lcplpostop
	      )
	end
    end;;
	  



(* fonction qui �tant donn� une liste de positions de stop et une liste de positions
   de met trouv�s dans la meme phase d'une miniseq, toutes deux suppos�es non vides,
   renvoit la liste des CDS de type MetEnd

   Attention : on suppose que les stop et met sont dans la MEME phase
   si cela n'est pas le cas alors on va trouver des CDS non mod 3 aberrants
*)
let listCDSMetEnd i posendsearch lpostop lposmet =
  if (lpostop=[] || lposmet=[]) then
    begin
      raise (Invalid_argument "listCDSMetEnd"); 
      []
    end
  else
    begin
      let poslastop = last lpostop in 
      let lposmetok = List.filter (fun posmet -> (poslastop < posmet)) lposmet in
	List.map (fun pmet -> CDS.create MetEnd pmet (posendsearch-i)) lposmetok
    end;;



(* fonction qui prend une miniseq (= liste de nt issus de la concat�nation des exons 
   d'un transcrit, ainsi qu'une phase i (0,1 ou 2) et renvoit la liste des objets 
   CDS qui se trouvent dans cette phase dans miniseq (de type MetStop, BegStop, MetEnd
   ou BegEnd)
*)
let listCDSinframei i miniseq =
  let szmseq = List.length miniseq in

  (* taille modulo 3 maximale dans miniseq qui d�bute en position i *)
    let sizesearch = size_window_search i szmseq in

    let lpostop = sig_search_down 
		  (Seq.of_list miniseq) 
		  i
		  sizesearch
		  3 
		  stop in

    let lposmet = sig_search_down 
		    (Seq.of_list miniseq) 
		    i
		    sizesearch
		    3 
		    met in
      
  let nstop = List.length lpostop and nmet = List.length lposmet in
    
    if (nstop=0) then
      begin
	if (nmet=0) then
	  begin
	    [CDS.create BegEnd i (sizesearch-i)]
	  end
	else (* nmet !=0 *)
	  begin
	    List.map (fun posmet -> CDS.create MetEnd posmet (sizesearch-i)) lposmet
	  end
      end

    else  (* nstop !=0 *)
      begin
	if (nmet=0) then
	  begin
	    [CDS.create BegStop i ((List.nth lpostop 0)-1)]
	  end
	else (* nmet !=0, cas plus compliqu� *)
	  begin
	    let lMetStopinit = listCDSMetStopInitial lpostop lposmet in
	    let lMetStopintern = listCDSMetStopInternal lpostop lposmet in
	    let lMetEnd = listCDSMetEnd i sizesearch lpostop lposmet in
	    let lMetStopintern_and_MetEnd = List.append lMetStopintern lMetEnd in
	      if (lMetStopinit = []) then (* pas de methionine avant le premier stop -> BegStop *)
		(CDS.create BegStop i ((List.nth lpostop 0)-1))::lMetStopintern_and_MetEnd
	      else
		List.append lMetStopinit lMetStopintern_and_MetEnd 
	  end
      end;;


(* tsizeex2tposexinminiseq prend en entr�e un tableau de tailles d'exons et 
   renvoit le tableau des couples (deb,fin) de ces exons dans la miniseq
   correspondant aux exons concat�n�s.
   Ex : tsizeex2tposexinminiseq [|5;8;2|] renvoit [|(0, 4); (5, 12); (13, 14)|] *)
let tsizeex2tposexinminiseq tsizeex =
  try
    (let i = ref 0 and n = Array.length tsizeex and scour = ref 0 in
     let tposexinminiseq = Array.create n (0,0) in
       while(!i<n) do
	 tposexinminiseq.(!i) <- (!scour,!scour + tsizeex.(!i) -1);
	 scour := !scour + tsizeex.(!i);
	 incr i;
       done;
       tposexinminiseq
    )
  with
      Invalid_argument _ ->  Common.print_log "Pb in the function tsizeex2tposexinminiseq \n"; [||]



(* n est une soustraction mod p, par cq s'il est st n�gatif, il faut lui ajouter
   un nombre suffisant de fois p pour qu'il redevienne positif 
   Ajout� en mars 2006 pour pallier au probl�me des hsps qui d�bordent des exons
   de la miniseq car ils avaient une phase n�gative 
*)
let rend_positif_modp n p =
  if(n>=0) then
    n 
  else
    begin
      let npos = ref n in
	while(!npos<0) do
	  npos:=!npos+p;
	done;
	!npos
    end;;



(* fonction qui prend un g�ne = transcrit et qui s'il y a des hsps prot�iques
   renvoit le couple bool�en * int = (sameframe,thisframe) correspondant � si les
   hsps prot sont dans la meme phase et cette phase (sinon phase premier hsp) 

   Extraite de l'ancienne fonction gene2genewithcds
*)
let protframe_in_gene genomseq g =
  let gseqofg = gene2gseq genomseq g in
  let sizegseqofg = List.length gseqofg in
    
  (* la liste des tailles des exons du g�ne (en fait transcrit) g *)
  let lsizeexofg = List.map2 
		     (fun b e -> e-b+1)
		     (Array.to_list (Gene.gbegex g)) 
		     (Array.to_list (Gene.gendex g)) in
    
  (* les tableau des positions et des d�buts des exons dans la miniseq = concat�nation de tous les exons de g *)
  let tposexinminiseq = try (tsizeex2tposexinminiseq (Array.of_list lsizeexofg)) with Invalid_argument _ ->  Common.print_log "in the function tsizeex2tposexinminiseq\n"; [||] in
  let tdebexinminiseq = Array.map fst tposexinminiseq in
   
  (* Liste des hsps prot�iques du g�ne accompagn�es des no des exons dont elles font partie *)
  let tex = Gene.texon g in
  let lhpnoex = try 
    (List.flatten 
       (Array.to_list 
	  (Array.mapi 
	     (fun i ex -> 
		List.map 
		(fun h -> (h,i))
		(List.sort 
		   Hsp.compare 
		   (List.filter Hsp.isprot (Exon.lhsp ex))))
	     tex)))
  with Invalid_argument _ ->  Common.print_log "Pb in lhpnoex\n"; [] in

  (* Recherche des nos des phases de toutes les hsp prot�iques : on s'attend a ce que ce soit le meme 
     -> sameorfallhp � true *)

  (* Liste des d�buts dans la miniseq de tous les hsps prot�iques du transcrit 
     attention : ici on peut avoir des d�buts n�gatifs 
  *)
  let lbeghpinminiseq = try (List.map 
			       (fun (hp,noex) -> tdebexinminiseq.(noex) + ((Hsp.debg hp)-(Exon.gbeg tex.(noex)))) 
			       lhpnoex
			    ) 
  with Invalid_argument _ -> Common.print_log "in the computation of lbeghpinminiseq\n"; [] in

  (* le dernier HSP prot�ique *)
  let lsthp = try (fst (last lhpnoex)) with Failure s -> Hsp.null in
    
  (* d�but du premier et fin du dernier HSP prot�ique dans la miniseq.
     Attention !!!!! : ces chiffres n'ont pas de sens quand les hsps prot�iques
     d�passent d'un exon de la miniseq *)
  let dp = (try (List.hd lbeghpinminiseq) with Failure s -> 0) 
  and fp = (try ((last lbeghpinminiseq) + ((Hsp.fing lsthp)-(Hsp.debg lsthp)+1)) with Failure s -> 0) in

  (* num�ros des phases des diff�rents hsps prot�iques *)
  let lnophasehp =   List.map (fun begh -> rend_positif_modp (begh mod 3) 3) lbeghpinminiseq in

  let nophasefsthp = try (List.hd lnophasehp) with Failure s -> 0 in
  let sameorfallhp = List.fold_left (fun b nophase -> b && (nophasefsthp=nophase)) true lnophasehp in
    (sameorfallhp,nophasefsthp);;


 
(* tposexinminiseq2noexwhereposis prend en entr�e une position dans une miniseq d'exons concat�n�s
   et le tableau des positions des exons dans cette miniseq (tableau de couples d'entiers)
   et renvoit le no de l'exon qui contient cette position.
   Ex : tposexinminiseq2noexwhereposis 7 [|(0, 4); (5, 12); (13, 14)|] renvoit 1, comme exon no 1 (ref 0) *)
let tposexinminiseq2noexwhereposis pos tposexinminiseq =
  try
    (let noexwhereposis = ref (-1) in
     let u = Array.iteri (fun i (gb,ge) -> 
			    if((gb<=pos) && (pos<=ge)) then noexwhereposis := i
			 )
	       tposexinminiseq in
       if(!noexwhereposis=(-1)) then
	 begin
	   (* Common.print_log "Pb in tposexinminiseq2noexwhereposis because of !noexwhereposis=(-1)\n"; *)
	   0
	 end
       else
	 !noexwhereposis
    )
  with
      Invalid_argument _ ->  Common.print_log "Invalid_arg in tposexinminiseq2noexwhereposis\n"; 0


let molec2string m =
  match m with
      P s -> if s="" then "NoneP" else s
    |A s -> if s="" then "NoneA" else s
    |_ ->  "None"



(* fonction qui a chaque gene = transcrit associe son CDS, ou plus pr�cis�ment son deb, fin de CDS *)
let gene2genewithCDS_mieux comparecds genomseq g =
  (* miniseq du gene = transcrit, cad concat�nation des exons *)
  let miniseq = gene2gseq genomseq g in
    
  (* la liste des tailles des exons du g�ne (en fait transcrit) g 
     cela sert pour passer en coordonn�es globales et non plus en coord
     dans la miniseq 
  *)
  let tgbegex = Gene.gbegex g 
  and tgendex = Gene.gendex g in
  
  let lsizeexofg = List.map2 
		     (fun b e -> e-b+1)
		     (Array.to_list tgbegex) 
		     (Array.to_list tgendex) in
  let tposexinminiseq = tsizeex2tposexinminiseq (Array.of_list lsizeexofg) in 
  let tdebexinminiseq = Array.map fst tposexinminiseq in
      

  let (phsp_in_same_frame,pframe) = protframe_in_gene genomseq g in
				      
    if (Gene.existprot g && phsp_in_same_frame) then
      begin
	if (not (Gene.existrna g)) then (* si que prot et prot ok alors pas de pred d'utr *)
	  g
	else (* si il y a prot et que ok niveau phase alors recherche des CDS uniquement dans CETTE phase *)
	  begin
	    let lcdsinpframe_sorted_bysize = List.sort comparecds (listCDSinframei pframe miniseq) in
	      if (lcdsinpframe_sorted_bysize = []) then
		Gene.null
	      else
		begin
		  let thecds = List.hd lcdsinpframe_sorted_bysize in
		  let begcds = CDS.gbeg thecds 
		  and endcds = CDS.gend thecds in
		  let noexbegcds = tposexinminiseq2noexwhereposis begcds tposexinminiseq 
		  and noexendcds = tposexinminiseq2noexwhereposis endcds tposexinminiseq in
		    (* Common.print_log ("Gene.arnpal g = "^(molec2string (Gene.arnpal g))^" ; Gene.protpale g = "^(molec2string (Gene.protpale g))^"\n");
		    Common.print_log ("thecds : type = "^(typeCDS2string (CDS.ctype thecds))^" ; deb = "^(string_of_int (CDS.gbeg thecds))^" ; fin = "^(string_of_int (CDS.gend thecds))^" ; noexbegcds = "^(string_of_int noexbegcds)^" ; noexendcds = "^(string_of_int noexendcds)^"\n");
		    *)
		    let globbegcds = tgbegex.(noexbegcds)+(begcds-tdebexinminiseq.(noexbegcds)) 	 
		    and globendcds = tgbegex.(noexendcds)+(endcds-tdebexinminiseq.(noexendcds)) in
		      (* Common.print_log ("thecds : globbegcds = "^(string_of_int globbegcds)^" ; globendcds = "^(string_of_int globendcds)^"\nGene.gbeg g = "^(string_of_int (Gene.gbeg g))^" ; Gene.gend g = "^(string_of_int (Gene.gend g))^"\n");
		      *)	
		      Gene.create (Gene.gbeg g) (Gene.gend g) globbegcds globendcds (Gene.texon g) (Gene.bring g) (Gene.nbex g) (Gene.gbegex g) (Gene.gendex g) (Gene.protpale g) (Gene.arnpal g) (Gene.fndsigex g) (Gene.commg g) (Gene.pbbring g)
		end
	  end
      end

    else  (* ici il n'y a pas de prot ou leurs hsps sont dans des phases diff�rentes 
	     On recherche donc les CDS dans les trois phases *)
      begin  
	let llcdsineachframe = List.map (fun i -> List.sort comparecds (listCDSinframei i miniseq)) [0;1;2] in
	let allcdsallframe_sorted_bysize = List.fold_left (List.merge comparecds) [] llcdsineachframe in
	  (*List.iter 
	    (fun cds -> Common.print_log ("cds : type = "^(typeCDS2string (CDS.ctype cds))^" ; deb = "^(string_of_int (CDS.gbeg cds))^" ; fin = "^(string_of_int (CDS.gend cds))^"\n")
	    ) 
	    allcdsallframe_sorted_bysize;*)	  

	  if (allcdsallframe_sorted_bysize = []) then
	    Gene.null
	  else
	    begin
	      let thecds = List.hd allcdsallframe_sorted_bysize in 
	      let begcds = CDS.gbeg thecds 
	      and endcds = CDS.gend thecds in
	      let noexbegcds = tposexinminiseq2noexwhereposis begcds tposexinminiseq 
	      and noexendcds = tposexinminiseq2noexwhereposis endcds tposexinminiseq in
		(* Common.print_log ("Gene.arnpal g = "^(molec2string (Gene.arnpal g))^" ; Gene.protpale g = "^(molec2string (Gene.protpale g))^"\n");
		Common.print_log ("thecds : type = "^(typeCDS2string (CDS.ctype thecds))^" ; deb = "^(string_of_int (CDS.gbeg thecds))^" ; fin = "^(string_of_int (CDS.gend thecds))^" ; noexbegcds = "^(string_of_int noexbegcds)^" ; noexendcds = "^(string_of_int noexendcds)^"\n");
		*)
		let globbegcds = tgbegex.(noexbegcds)+(begcds-tdebexinminiseq.(noexbegcds)) 	 
		and globendcds = tgbegex.(noexendcds)+(endcds-tdebexinminiseq.(noexendcds)) in
		  (* Common.print_log ("thecds : globbegcds = "^(string_of_int globbegcds)^" ; globendcds = "^(string_of_int globendcds)^"\nGene.gbeg g = "^(string_of_int (Gene.gbeg g))^" ; Gene.gend g = "^(string_of_int (Gene.gend g))^"\n");	
		  *)
		  Gene.create (Gene.gbeg g) (Gene.gend g) globbegcds globendcds (Gene.texon g) (Gene.bring g) (Gene.nbex g) (Gene.gbegex g) (Gene.gendex g) (Gene.protpale g) (Gene.arnpal g) (Gene.fndsigex g) (Gene.commg g) (Gene.pbbring g)
	    end
      end;;






























(**************************** Anciennes fonctions de recherche d'orf puis cds ************************************)
(**************************** Ne marchent pas bien ***************************************************************)







(* tpostop2longestorf prend en entr�e un tableau de positions de stops en phase frame (0,1 ou 2) 
   ainsi que la taille de la s�quence (concat d'exons) dans laquelle ces stops ont �t� trouv�s
   et renvoit le d�but et la fin de l'orf la plus longue. 
   C'est � l'int�rieur de celle-ci que se fera la recherche de la m�thionine start *)
let tpostop2longestorf size frame tpostop =
  try
    (let max = ref 0 and dfmaxorf = ref (0,0) in
       (* il n'y a pas de stop en phase frame durant size 
	  la plus longue orf est donc la s�quence entiere � partir de 
	  la position frame et sur une dist mod 3 = 0 *)
       if(tpostop=[||]) then
	 begin
	   max := adjust2mod3 (size-frame);
	   dfmaxorf := (frame,!max+frame-1);
	 end
	   (* il y a n stops en phase frame durant size 
	      en ce cas soit la plus longue orf se situe entre deux stops
	      soit (apres le while) entre le dernier stop (postop.(n-1)) et
	      size *)
       else
	 begin
	   let n = Array.length tpostop and i = ref 0 in
	     max := adjust2mod3 (tpostop.(0)-frame);
	     dfmaxorf := (frame,!max+frame-1);
	     while (!i<n-1) do
	       if ((tpostop.(!i+1)-(tpostop.(!i)+3))> (!max)) then
		 begin
		   max := (tpostop.(!i+1)-(tpostop.(!i)+3));
		   dfmaxorf := (tpostop.(!i)+3,tpostop.(!i+1)-1);
		 end;
	       incr i;
	     done;
	     if((adjust2mod3(size-(tpostop.(n-1)+3)))> (!max)) then
	       begin
		 max := adjust2mod3(size-(tpostop.(n-1)+3));
		 dfmaxorf := (tpostop.(n-1)+3,size-((size-(tpostop.(n-1)+3)) mod 3)-1)
	       end;
	 end;
       (!max,!dfmaxorf)
    )
  with
      Invalid_argument _ -> Common.print_log "Pb dans tpostop2longestorf\n"; (0,(0,0))



(* tpostop2longorfwithprot est identique � la pr�cedente sauf qu'elle prend en plus
   un d�but et une fin dans la miniseq de la premiere et de la derniere hsp prot�ique
   et qu'elle renvoit la longueur maximale et le debut et la fin dans la miniseq
   de la plus grande orf qui en plus contient dp-fp.
   S'il n'y a pas de telle orf alors on se rabat sur tpostop2longestorf dans la phase
   des hsps prot�iques (si la meme cf recherche de gene2genewithcds plus bas).
   C'est � l'int�rieur de l'orf trouv�e que se fera la recherche de la m�thionine start *)
let tpostop2longestorfwithprot dp fp size frame tpostop =
  try
    (let max = ref 0 and dfmaxorf = ref (0,0) in
       (* il n'y a pas de stop en phase frame durant size 
	  la plus longue orf est donc la s�quence entiere � partir de 
	  la position frame et sur une dist mod 3 = 0 *)
       if(tpostop=[||]) then
	 begin
	   (* ici pas de stop donc obligatoirement les positions dp-fp sont incluses dans l'orf *)
	   max := adjust2mod3 (size-frame);
	   dfmaxorf := (frame,!max+frame-1);
	 end
	   (* il y a n stops en phase frame durant size 
	      en ce cas soit la plus longue orf se situe entre deux stops
	      soit (apres le while) entre le dernier stop (postop.(n-1)) et
	      size *)
       else
	 begin
	   let n = Array.length tpostop and i = ref 0 in
	     if ((frame<=dp) && (fp<=(adjust2mod3 (tpostop.(0)-frame)+frame))) then
	       begin
		 max := adjust2mod3 (tpostop.(0)-frame);
		 dfmaxorf := (frame,!max+frame-1);
	       end;
	     while (!i<n-1) do
	       if (((tpostop.(!i+1)-(tpostop.(!i)+3))> (!max)) && (((tpostop.(!i)+3)<=dp) && (fp<=tpostop.(!i+1)))) then
		 begin
		   max := (tpostop.(!i+1)-(tpostop.(!i)+3));
		   dfmaxorf := (tpostop.(!i)+3,tpostop.(!i+1)-1);
		 end;
	       incr i;
	     done;
	     if (((adjust2mod3(size-(tpostop.(n-1)+3)))> (!max)) && ((tpostop.(n-1)+3<=dp) && (fp<=size-((size-(tpostop.(n-1)+3)) mod 3)))) then
	       begin
		 max := adjust2mod3(size-(tpostop.(n-1)+3));
		 dfmaxorf := (tpostop.(n-1)+3,size-((size-(tpostop.(n-1)+3)) mod 3)-1)
	       end;
	 end;
       (!max,!dfmaxorf)
    )
  with
      Invalid_argument _ -> Common.print_log "Pb dans tpostop2longestorfwithprot\n"; (0,(0,0))




(* gene2genewithcds prend en entr�e une s�quence g�nomique et un g�ne qui est sur son brin forward
   (dans le cas d'un g�ne sur le brin Reverse de genomseq on remplacera genomseq par son compl�mentaire)
   et renvoit le g�ne mais dans lequel les champs d�but et fin de cds ont �t� remplis correctement.
   La recherche basique (cad sans prot�ine pour orienter sur la phase) se fait en prenant la phase qui contient
   globalement la plus longue orf et dans cette phase on recherche le plus grand cds ou encore la methionine
   (start) la plus en amont. *)
let gene2genewithcds pprot genomseq g = 
  let gseqofg = gene2gseq genomseq g in
  let sizegseqofg = List.length gseqofg in

  (* la liste des tailles des exons du g�ne (en fait transcrit) g *)
  let lsizeexofg = List.map2 
		     (fun b e -> e-b+1)
		     (Array.to_list (Gene.gbegex g)) 
		     (Array.to_list (Gene.gendex g)) in
    
  (* les tableau des positions et des d�buts des exons dans la miniseq = concat�nation de tous les exons de g *)
  let tposexinminiseq = try (tsizeex2tposexinminiseq (Array.of_list lsizeexofg)) with Invalid_argument _ ->  Common.print_log "dans la fonction tsizeex2tposexinminiseq\n"; [||] in
  let tdebexinminiseq = Array.map fst tposexinminiseq in
   
  (* Liste des hsps prot�iques du g�ne accompagn�es des no des exons dont elles font partie *)
  let tex = Gene.texon g in
  let lhpnoex = try 
    (List.flatten 
       (Array.to_list 
	  (Array.mapi 
	     (fun i ex -> 
		List.map 
		(fun h -> (h,i))
		(List.sort 
		   Hsp.compare 
		   (List.filter Hsp.isprot (Exon.lhsp ex))))
	     tex)))
  with Invalid_argument _ ->  Common.print_log "Pb dans lhpnoex\n"; [] in



  (* Recherche des nos des phases de toutes les hsp prot�iques : on s'attend a ce que ce soit le meme 
     -> sameorfallhp � true *)

  (* Liste des d�buts dans la miniseq de tous les hsps prot�iques du transcrit 
     attention : ici on peut avoir des d�buts n�gatifs 
  *)
  let lbeghpinminiseq = try (List.map 
			       (fun (hp,noex) -> tdebexinminiseq.(noex) + ((Hsp.debg hp)-(Exon.gbeg tex.(noex)))) 
			       lhpnoex
			    ) 
  with Invalid_argument _ -> Common.print_log "dans le calcul de lbeghpinminiseq\n"; [] in

  (* le dernier HSP prot�ique *)
  let lsthp = try (fst (last lhpnoex)) with Failure s -> Hsp.null in
    
  (* d�but du premier et fin du dernier HSP prot�ique dans la miniseq.
     Attention !!!!! : ces chiffres n'ont pas de sens quand les hsps prot�iques
     d�passent d'un exon de la miniseq *)
   let dp = (try (List.hd lbeghpinminiseq) with Failure s -> 0) 
   and fp = (try ((last lbeghpinminiseq) + ((Hsp.fing lsthp)-(Hsp.debg lsthp)+1)) with Failure s -> 0) in

   (* num�ros des phases des diff�rents hsps prot�iques *)
   let lnophasehp =   List.map (fun begh -> rend_positif_modp (begh mod 3) 3) lbeghpinminiseq in

   let nophasefsthp = try (List.hd lnophasehp) with Failure s -> 0 in
   let sameorfallhp = List.fold_left (fun b nophase -> b && (nophasefsthp=nophase)) true lnophasehp in 


     (* Si on veut faire une recherche en fonction des prot�ines
	et qu'il exiqte une prot�ine dans le g�ne = transcrit, 
	et que toutes les hsps prot�iques sont dans la meme phase dans la miniseq 
	alors on n'inspecte que cette phase 
     *)
    if ((pprot g) && (Gene.existprot g) && (sameorfallhp)) then
      begin
	if (not (Gene.existrna g)) then  (* on n'a QUE des prot�ines et elles sont toutes dans la meme phase
					    alors je veux cr�er un g�ne (tr) SANS utr -> debcds=debtr, fincds=fintr
					 *)
	  begin
(*	    output_string stderr ("Je suis dans prot ok sans arn pour le transcrit d'arn ppal "^(molec2string (Gene.arnpal g))^" et de protpale "^(molec2string (Gene.protpale g))^"\n");*)
	    try (Gene.create (Gene.gbeg g) (Gene.gend g) (Gene.gbeg g) (Gene.gend g) (Gene.texon g) (Gene.bring g) (Gene.nbex g) (Gene.gbegex g) (Gene.gendex g) (Gene.protpale g) (Gene.arnpal g) (Gene.fndsigex g) (Gene.commg g) (Gene.pbbring g)) with Invalid_argument s -> Common.print_log "au niveau du premier Gene.create\n"; Gene.null  
	  end
	else
	  begin
	    (* Le tableau des positions des stops uniquement dans la phase donnee par la premiere hsp 
	       plus tard : donn�e par l'hsp qui a subi le moins de pb depuis le d�part!!! *)
	    let tstop = try
	      (if ((size_window_search nophasefsthp sizegseqofg )<0) then
		 [||]
	       else
		 (Array.of_list 
		    (sig_search_down 
		       (Seq.of_array (Array.of_list gseqofg)) 
		       nophasefsthp
		       (size_window_search nophasefsthp sizegseqofg ) (*(sizegseqofg-6+nophasefsthp) *)
		       3 
		       stop)
		 )
	      )
	    with Invalid_argument _ -> Common.print_log "Pb dans tstop de pprot\n"; [||] in
	      
	      
	    (* le d�but et la fin de l'orf la plus longue dans la phase donn�e par la premiere hsp prot�ique 
	       et si une telle orf n'existe pas alors on prend la plus grande orf dans cette phase 
	       Remarque : dp et fp sont le d�but et la fin resp de la 1ere et de la derniere hsp proteique *)
	    let (dorf,forf) = 
	      if ((fst (tpostop2longestorfwithprot dp fp sizegseqofg nophasefsthp tstop)) = 0) then 
		snd (tpostop2longestorf sizegseqofg nophasefsthp tstop)
	      else
		snd (tpostop2longestorfwithprot dp fp sizegseqofg nophasefsthp tstop)	  
	    in
	     
	    let lstart =  try
	      (if((forf-dorf)<0) then 
		 [] 
	       else
		 (sig_search_down 
		    (Seq.of_array (Array.of_list gseqofg)) 
		    dorf  
		    (forf-dorf) 
		    3
		    met)
	      )
	    with Invalid_argument _ -> Common.print_log "Pb dans lstart de pprot\n"; [] in
	      

	(* S'il y a un start alors on prend comme d�but de CDS ce start
	   et sinon on prend comme d�but de CDS le d�but de l'orf dans 
	   la bonne phase, cad dorf0. 
	   Comme fin de CDS on prend toujours la fin de l'orf ds la 
	   bonne phase, cad forf0 
	*)
	    let dcds = try (List.hd lstart) with Failure _ -> dorf and fcds = forf in
	    let noexdcds = tposexinminiseq2noexwhereposis dcds tposexinminiseq and noexfcds = tposexinminiseq2noexwhereposis fcds tposexinminiseq in
	      try (Gene.create (Gene.gbeg g) (Gene.gend g) ((Gene.gbegex g).(noexdcds)+(dcds-tdebexinminiseq.(noexdcds))) ((Gene.gbegex g).(noexfcds)+(fcds-tdebexinminiseq.(noexfcds))) (Gene.texon g) (Gene.bring g) (Gene.nbex g) (Gene.gbegex g) (Gene.gendex g) (Gene.protpale g) (Gene.arnpal g) (Gene.fndsigex g) (Gene.commg g) (Gene.pbbring g)) with Invalid_argument s -> Common.print_log "au niveau du premier Gene.create\n"; Gene.null
		
	  end 
      end (* fin de if ((pprot g) && (Gene.existprot g) && (sameorfallhp)) *) 

    else
      (* s'il n'exite pas de prot�ine dans le g�ne = transcrit, alors regarde toutes les phases
	 cad 0, 1 et 2 *)
      begin
	(* les tableaux des positions des stops dans chacune des trois phases *)
	let [t0;t1;t2] =  try
	  (List.map 
	     (fun i ->
		Array.of_list 
		(if((size_window_search i sizegseqofg)<0) then
		   []
		 else
		   (sig_search_down 
		      (Seq.of_array (Array.of_list gseqofg)) 
		      i 
		      (size_window_search i sizegseqofg) (*(sizegseqofg-6+i) *)
		      3 
		      stop)
		)
	     )
	     [0;1;2])
	with Invalid_argument _ -> Common.print_log "Pb dans [t0;t1;t2] de else pprot\n"; [] in
	  
	(* les couples (taille max orf,(debut,fin) de l'orf) dans chacune des trois phases *)
	let [(m0,(dorf0,forf0));(m1,(dorf1,forf1));(m2,(dorf2,forf2))] = 
	  List.map2 
	    (tpostop2longestorf sizegseqofg) 
	    [0;1;2] 
	    [t0;t1;t2] 
	in
	  
	(* dans la phase (0, 1 ou 2) dans laquelle l'orf la + grande a �t� trouv�e
	   on recherche un start et si on entrouve on prend le premier et sinon on prend
	   comme d�but de "CDS" le d�but de l'orf dans la phase s�l�ctionn�e *)
	let m = max (max m0 m1) m2 in
	  
	  if (m=m0) then
	    begin
	      let lstart = try
		(if((forf0-dorf0)<0) then
		   []
		 else
		   (sig_search_down 
		      (Seq.of_array (Array.of_list gseqofg)) 
		      dorf0  
		      (forf0-dorf0) 
		      3
		      met)
		)
	      with Invalid_argument _ -> Common.print_log "Pb dans lstart de m=m0\n"; [] in
		
		
	      (* S'il y a un start alors on prend comme d�but de CDS ce start
		 et sinon on prend comme d�but de CDS le d�but de l'orf dans 
		 la bonne phase, cad dorf0. 
		 Comme fin de CDS on prend toujours la fin de l'orf ds la 
		 bonne phase, cad forf0 *)
	      let dcds = try (List.hd lstart) with Failure _ -> dorf0 and fcds = forf0 in
	      let noexdcds = tposexinminiseq2noexwhereposis dcds tposexinminiseq and noexfcds = tposexinminiseq2noexwhereposis fcds tposexinminiseq in
		try (Gene.create (Gene.gbeg g) (Gene.gend g) ((Gene.gbegex g).(noexdcds)+(dcds-tdebexinminiseq.(noexdcds))) ((Gene.gbegex g).(noexfcds)+(fcds-tdebexinminiseq.(noexfcds))) (Gene.texon g) (Gene.bring g) (Gene.nbex g) (Gene.gbegex g) (Gene.gendex g) (Gene.protpale g) (Gene.arnpal g) (Gene.fndsigex g) (Gene.commg g) (Gene.pbbring g)) with Invalid_argument s -> Common.print_log "au niveau du troisi�me Gene.create\n"; Gene.null
	    end (* fin du if (m=m0) *)
		
	  else 
	    begin
	      if (m=m1) then
		begin
		  let lstart = try
		    (if ((forf1-dorf1)<0) then 
		       []
		     else
		      (sig_search_down 
			 (Seq.of_array (Array.of_list gseqofg)) 
			 dorf1 
			 (forf1-dorf1)  
			 3
			 met)
		    )
		  with
		      Invalid_argument _ -> Common.print_log "Pb dans lstart de m=m1\n"; [] in
		    
		    begin
		      let dcds = try (List.hd lstart) with Failure _ -> dorf1 and fcds = forf1 in
		      let noexdcds = tposexinminiseq2noexwhereposis dcds tposexinminiseq and noexfcds = tposexinminiseq2noexwhereposis fcds tposexinminiseq in
			try (Gene.create (Gene.gbeg g) (Gene.gend g) ((Gene.gbegex g).(noexdcds)+(dcds-tdebexinminiseq.(noexdcds))) ((Gene.gbegex g).(noexfcds)+(fcds-tdebexinminiseq.(noexfcds))) (Gene.texon g) (Gene.bring g) (Gene.nbex g) (Gene.gbegex g) (Gene.gendex g) (Gene.protpale g) (Gene.arnpal g) (Gene.fndsigex g) (Gene.commg g) (Gene.pbbring g)) with Invalid_argument s -> Common.print_log "au niveau du quatri�me Gene.create\n"; Gene.null	    
		    end
		end (* fin du if (m=m1) *)
		  
	      else
		let lstart = try
		  (if ((forf2-dorf2)<0) then 
		     []
		   else
		     (sig_search_down 
			(Seq.of_array (Array.of_list gseqofg)) 
			dorf2 
			(forf2-dorf2)  
			3
			met)
		  )
		with
		    Invalid_argument _ -> Common.print_log "Pb dans lstart de m=m2\n"; [] in  
		  
		  begin
		    let dcds = try (List.hd lstart) with Failure _ -> dorf2 and fcds = forf2 in
		    let noexdcds = tposexinminiseq2noexwhereposis dcds tposexinminiseq and noexfcds = tposexinminiseq2noexwhereposis fcds tposexinminiseq in
		      try (Gene.create (Gene.gbeg g) (Gene.gend g) ((Gene.gbegex g).(noexdcds)+(dcds-tdebexinminiseq.(noexdcds))) ((Gene.gbegex g).(noexfcds)+(fcds-tdebexinminiseq.(noexfcds))) (Gene.texon g) (Gene.bring g) (Gene.nbex g) (Gene.gbegex g) (Gene.gendex g) (Gene.protpale g) (Gene.arnpal g) (Gene.fndsigex g) (Gene.commg g) (Gene.pbbring g)) with Invalid_argument s -> Common.print_log "au niveau du cinq�me Gene.create\n"; Gene.null  			        
		  end
	    end
      end
	      
   
	









