(* ALPHABET des PROTEINES (les acides aminés) : DEUXIEME INTERFACE D'ALPHABET AVEC LE RESTE DE L'APPLICATION *)


type prot_t = 
  |F 
  |L  
  |I 
  |M 
  |V 
  |S 
  |P 
  |T 
  |A 
  |Y 
  |H 
  |Q 
  |N 
  |K 
  |D 
  |E 
  |C 
  |W 
  |R 
  |G
  |Z  (* Z est là pour le codon stop *)
  |X  (* X est là pour la traduction en acide aminé de tout codon contenant un N *)
  

let astrp =[|"F";"L";"I";"M";"V";"S";"P";"T";"A";"Y";"H";"Q";"N";"K";"D";"E";"C";"W";"R";"G";"Z";"X";"B"|];;

let of_charp  = function
  |	'F' -> F
  |	'L' -> L
  |	'I' -> I
  |	'M' -> M
  |	'V' -> V
  |	'S' -> S
  |	'P' -> P
  |	'T' -> T
  |	'A' -> A
  |	'Y' -> Y
  |	'H' -> H
  |	'Q' -> Q
  |	'N' -> N
  |	'K' -> K
  |	'D' -> D
  |	'E' -> E
  |	'C' -> C
  |	'W' -> W
  |	'R' -> R
  |	'G' -> G
  |     'Z' -> Z 
  |     'X' -> X
  |     'x' -> X
  | _   -> raise Not_found


	    
let to_charp = function 
  |	F  -> 'F'
  |	L  -> 'L'
  |	I  -> 'I'
  |	M  -> 'M'
  |	V  -> 'V'
  |	S  -> 'S'
  |	P  -> 'P'
  |	T  -> 'T'
  |	A  -> 'A'
  |	Y  -> 'Y'
  |	H  -> 'H'
  |	Q  -> 'Q'
  |	N  -> 'N'
  |	K  -> 'K'
  |	D  -> 'D'
  |	E  -> 'E'
  |	C  -> 'C'
  |	W  -> 'W'
  |	R  -> 'R'
  |	G  -> 'G'
  |     Z  -> 'Z'
  |     X  -> 'X'
  

let to_stringp = function 
  |	F  -> "F"
  |	L  -> "L"
  |	I  -> "I"
  |	M  -> "M"
  |	V  -> "V"
  |	S  -> "S"
  |	P  -> "P"
  |	T  -> "T"
  |	A  -> "A"
  |	Y  -> "Y"
  |	H  -> "H"
  |	Q  -> "Q"
  |	N  -> "N"
  |	K  -> "K"
  |	D  -> "D"
  |	E  -> "E"
  |	C  -> "C"
  |	W  -> "W"
  |	R  -> "R"
  |	G  -> "G"
  |     Z  -> "Z"
  |     X  -> "X"

let nullp = F

let cardp = Array.length astrp


let prot2string l =
  let str = String.create (List.length l) and i = ref 0 in
    List.iter (fun p -> str.[!i] <- to_charp p; incr i) l;
    str;;



