(* main.ml :fichier principal, interface utilisateur, permet d'afficher les modeles de gènes ainsi que leur traduction 
   en proteine dans deux fichiers resultat, qui pour l'instant ont les noms predeterminés RES_exogean1.gene et 
   RES_exogean1.prot respectivement
*)

open Config
open Collection
open Donnees
open TraiteCollection
open Premain
open Printing


(* Attention fonction principale du programme au 12/09/03 (car la paramétrisation de l'executable ne marche plus.
   hsppfile est la string correspondant au fichier des résultats tblastn sur le brin+, hspmfile est le string 
   correspondant au fichier des résultats tblastn sur le brin-, o est le channel de sortie du fichier 
   sur lequel on veut imprimer les résultats du blusping, ex : let o = open_out "nom_fich.txt", 
   ou simplement stdout 
   Attention : on suppose que l'on a seqfile.fasta, comp.fasta (son complementaire) et 
   la banque de proteines.fasta, en .ind, en .tab et en .seq. (resultat de lspextend)
   EXOGEAN doit connaitre la séquence entiere (sans les masquages) et doit pouvoir connaitre l'endroit 
   des masquages, soit par un fichier complementaire, soit de lui-meme en fonction -> soit des N ou n, 
   soit des majuscules/minuscules.
*)


let print_tous_genes () =
		read_configfile ();
		read_commandline ();
		if context.gseqfile = "" then Common.print_error "No genomic sequence specified !\nIf you need help, please type ./exogean --help\n";
		if (context.hsppfile = "") && (context.hspafile = "") then Common.print_error "No HSP specified !\nIf you need help, please type ./exogean --help\n";
		if not (context.hsppfile = "") && (context.bqpfile = "") then Common.print_error "No protein bank specified !\nIf you need help, please type ./exogean --help\n";
		(* lcgok est la liste des completegene ok, et lcgpseudo des pseudogenes *)
		let (((lcgok,lcgpseudo),lg),tgenearn,nomgseq) = find_modeles context in
		try
			Common.print_log "Exogean did its work ! Now writing output files\n";

			(* La fonction a appeler pour ecrire les fichiers de sortie suivant le type choisi *)
			let (f, ext, s) = get_type_sortie nomgseq context.outfile in 
			
			(* On prepare ce qui doit etre ecrit *)
			let tcgeneok = Array.of_list lcgok in
			let lcgok_sorted = List.sort CompleteGene.compare lcgok in
			let lcgpseudo_sorted = List.sort CompleteGene.compare lcgpseudo in

			let tgeneok = Array.of_list (List.map CompleteGene.gn lcgok) in
			let tcgeneok_sorted = Array.of_list lcgok_sorted in
			let tgeneok_sorted = Array.of_list (List.map CompleteGene.gn lcgok_sorted) in

			let lprotok = List.map (fun cg -> ((CompleteGene.nbexincds cg),(CompleteGene.prot cg))) lcgok in
			let lprotok_sorted = List.map (fun cg -> ((CompleteGene.nbexincds cg),(CompleteGene.prot cg))) lcgok_sorted in

			let ltransok = List.map (fun cg -> ((CompleteGene.nbextot cg),(CompleteGene.trans cg))) lcgok in
			let ltransok_sorted = List.map (fun cg -> ((CompleteGene.nbextot cg),(CompleteGene.trans cg))) lcgok_sorted in

			let tgenepseudo = Array.of_list (List.map CompleteGene.gn lcgpseudo) in
			let tgenepseudo_sorted = Array.of_list (List.map CompleteGene.gn lcgpseudo_sorted) in
		  
			  (* à terme il faudra mettre tcgeneok_sorted mais pour etre cohérent avec la sortie
			     gene.gtf on laisse tcgeneok *)
			let ttcgchev = Array.of_list (SegSeq.elements (cluspingcgene tcgeneok)) in
			let tugok = Array.map (UniteGenique.tcgchev2ug) ttcgchev in

			(* attention : à terme remplacer tgeneok, tprotok et tgenepseudo par idem_sorted
			   car en theorie le clusping présuppose que les objets soient triés 
			   par conséquent si on ne le fait pas on ne peut pas garantir que les gènes crées
			   soient bien crées sur le principe du chevauchement génomique
			*)

			(* on ecrit les modèles de gènes *)
			let file = if context.pipe = "gene" then stdout else open_out (s^".gene."^ext) in
			let noms_genes = List.flatten (print_transcrits f nomgseq lg file context.offset tgeneok) in
			flush file;
			close_out file;

			(* on ecrit les protéines en acide aminés dans un fichier fasta 
			   Attention : tprot est en réalité un tableau de couples (nbexincds,proteine) *)
			let file = if context.pipe = "prot" then stdout else open_out (s^".cds.fa") in
			List.iter2 (print_prot file) noms_genes lprotok;
			flush file;
			close_out file;

			(* on écrit les transcrits en nucleotides dans un fichier fasta *)
			let file = if context.pipe = "tr" then stdout else open_out (s^".cdna.fa") in
			List.iter2 (print_tr file) noms_genes ltransok;
			flush file;
			close_out file;

			(* on ecrit les pseudogenes *)
			let file = if context.pipe = "pseudo" then stdout else open_out (s^".pseudo."^ext) in
			Array.iteri (print_pseudo f nomgseq lg context.offset file) tgenepseudo;
			flush file;
			close_out file;

			(* on ecrit les arns pretraités *)
			let file = if context.pipe = "rna" then stdout else open_out (s^".rlgn."^ext) in
			Array.iteri (print_arn f nomgseq lg context.offset file) tgenearn; 
			flush file;
			close_out file;
			
			(* on écrit les gènes au format html (ajout de decembre 2005) *)
			let file = if context.pipe = "html" then stdout else open_out (s^".gene.html") in
			  print_sortie_html context.gseqfile context.hspafile context.hsppfile context.bqpfile nomgseq lg file context.offset tugok;  
			  flush file;
			  close_out file;
			
		with
			| Invalid_argument s -> Common.print_error ("Error while processing output files ("^s^")\n");
			    (* I do not understand this any more *)
			| Failure "hd" -> Common.print_error ("No annotation\n");
			| Sys_error s -> Common.print_error ("I/O Error while accessing an output file: \""^s^"\"\n");;


(* Application de la fonction principale du main *)
print_tous_genes ();; 



