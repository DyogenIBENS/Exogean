(* ALPHABET de l'ADN (les nucl�otides ou bases) : PREMIERE INTERFACE D'ALPHABET AVEC LE RESTE DE L'APPLICATION *)


type dna_t = 
  |A
  |T
  |G
  |C
  |N


(* inutilis� par la suite *)
let astra =[|"A";"T";"G";"C";"N"|];;

let of_chara  = function
  |	('A' |'a') -> A
  |	('T' |'t') -> T
  |	('G' |'g') -> G
  |	('C' |'c') -> C
  |	('N' |'n') -> N
  | _   -> raise Not_found


(* si on se sert de cette fonction par la suite il faudra r�flechir � la distinction majuscule/minuscule *)	    
let to_chara = function 
  |	 A -> 'A'
  |	 T -> 'T'
  |	 G -> 'G'
  |	 C -> 'C'
  |	 N -> 'N'


let complement = function 
  | A -> T
  | T -> A
  | G -> C
  | C -> G 
  | N -> N 

let nulla = A

(* inutilis� par la suite *)
let carda = Array.length astra


