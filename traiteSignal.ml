(* traiteSignal.ml : permet de rechercher un ensemble de signaux ou motifs (de meme longueur dans une meme recherche) au sein d'une séquence de nucléotides, ce en amont ou en aval d'une position donnée, et sur une distance donnée *)



open Seq
open Config


(* tsplit_down prend en argument une Seq seq et deux entiers que sont une taille de motif t_motif, 
   et un pas de coupure pas, et répertorie dans une table de hachage tous les motifs de taille t_motif 
   contenus dans la bigseq seq, par pas de step, en numérotant les positions d'apparition des motifs 
   en partant du début de seq, la position du début du signal étant rendue *)
let tsplit_down seq t_motif step =
  let n = Seq.length seq in
    try
      (let h = Hashtbl.create (n/step+1) and i = ref 0 in
	 while (!i+t_motif <= n) do
	   Hashtbl.add h !i (Seq.sub seq !i t_motif);
	   i:=!i+step;
	 done;
	 h)
    with 
      |Invalid_argument _ -> Common.print_log "c'est dans le tsplit_down qu'il y a un probleme\n"; Hashtbl.create 5;;





(* sig_search_down permet de chercher dans la séquence génomique seq (seq de nucléotides), 
   en aval de la position beg et pendant une distance de size, l'un des signaux génomiques inclus 
   dans la liste de signaux lseqsig (tous de meme taille, qui sont eux-memes des seq). 
   Elle permet d'avancer par pas de step, non necessairement égal à 1 (peut etre égal à 3 c
   omme pour la recherche d'un codon Stop). 
   Cette fonction renvoit la position dans gseq de la première occurence de l'un de ces motifs *)
let sig_search_down seq beg size step lseqsig =
  let t_motif = Seq.length (List.hd lseqsig) and lpos = ref [] in
    try
      (let hrechsig = tsplit_down (Seq.sub seq beg size) t_motif step in
      let i = ref 0 and iend = size-t_motif in
      let stop = ref false in


	while(!i<=iend) do
	
	  stop:= false;
	  
	  while(!i <= iend && not (!stop)) do
	    stop := List.fold_left (fun bool a -> bool || ((Hashtbl.find hrechsig (!i))=a)) false lseqsig;
	    i:=!i+step;
	  done;

	  if (!stop) then
	    begin
	      lpos := (beg+(!i-step))::(!lpos);
	    end;

	done; 
	List.rev !lpos)
  
    with
      |Invalid_argument s ->  Common.print_log ("pb dans sig_search_down avec beg = "^(string_of_int beg)^", size = "^(string_of_int size)^", step = "^(string_of_int step)^" et List.length lseqsig = "^(string_of_int (List.length lseqsig))^" et "^s^"\n"); List.rev !lpos;
      |Not_found -> List.rev !lpos;;
  




(* tsplit_up prend en argument une seq seq et deux entiers que sont une taille de motif t_motif, 
   et un pas de coupure pas, et répertorie dans une table de hachage tous les motifs de taille t_motif 
   contenus dans la Seq seq, par pas de step, en numérotant les positions d'apparition des motifs 
   en partant également du début de seq, la position du début du signal étant rendue *)
let tsplit_up seq t_motif step =
  let n = Seq.length seq in
    try
      (let h = Hashtbl.create (n/step+1) and i = ref (n-t_motif) in
	while (!i >= 0) do
	  Hashtbl.add h !i (Seq.sub seq !i t_motif);
	  i:=!i-step;
	done;
	h)
    with 
      |Invalid_argument _ -> Common.print_log "c'est dans le tsplit_up qu'il y a un probleme\n"; Hashtbl.create 5;;



(* Idem que sig_search_down mais en recherchant vers l'amont (up). Les positions rendues sont 
   celles de début de signal *)
let sig_search_up seq beg size step lseqsig =
  let t_motif = Seq.length (List.hd lseqsig) and lpos = ref [] in
    try
      (let hrechsig = tsplit_up (Seq.sub seq (beg-size+1) size) t_motif step in
       let i = ref (size-t_motif) and iend = (size-t_motif) mod step in
       let stop = ref false in

	 while(!i>=iend) do
      
	   stop:= false;
	
	   while((!i>=iend) && not (!stop)) do
	     stop := List.fold_left (fun bool a -> bool || ((Hashtbl.find hrechsig (!i))=a)) false lseqsig; 
	     i:=!i-step;
	   done;
	 
	   if (!stop) then
	     begin
	       lpos := (beg-size+1+(!i+step))::(!lpos);
	     end;
	   
	 done;
	 List.rev !lpos)
      
    with
      |Invalid_argument s ->  Common.print_log ("pb dans sig_search_up avec beg = "^(string_of_int beg)^", size = "^(string_of_int size)^", step = "^(string_of_int step)^" et List.length lseqsig = "^(string_of_int (List.length lseqsig))^" et "^s^"\n"); List.rev !lpos;
      |Not_found -> List.rev !lpos;;




