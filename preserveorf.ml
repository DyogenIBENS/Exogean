(* preserveorf.ml : pour préserver l'orf en cas de fusion d'hsps ou de recherche de signal *)


open Alphaadn
open Seq
open TraiteSignal



let stop = List.map Seq.of_array [[|T;A;G|];[|T;G;A|];[|T;A;A|]]

(* retourne la liste des stops entre le début et la fin de la seq génomique mot  *)
let lstop_between_cl mot = 
  sig_search_down mot 0 ((Seq.length mot)) 3 stop




(* retourne true s'il n'y a pas de STOP et retourne false s'il y a STOP dans la seq génomique mot
   Attention : si le booleen b precedant le mot est à false, c'est qu'il y a eu un pb
   dans la taille de l'exon obtenue par choix des signaux gp1,gp2, et en ce cas on veut eliminer
   de tels signaux d'ou le false *)
let pasStop (b,mot) = 
  if b then
    begin
      let listeStop = lstop_between_cl mot in
	match listeStop with
	    [] -> true
	  |_ -> false
    end
  else
    false


