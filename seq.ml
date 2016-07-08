(* seq.ml : includes the Seq module
   enables to put in an array sequences of size > Sys.max_array_length (= 4194303 
   for the 32 bits machines, 1814 billions for the 64 bits machines), thanks to
   the Bigarray module and more particularly to the submodule Array1 of Bigarray. 
   Be careful: in one of these bigseqs the type cannot be more than a simple variant type.
   (ie codable by integers) 
   Typically the sequences concerned here are the genomic sequences but for a question 
   of homogeneity we will use this module also for genomic sequences subsequences (eg splice signals),
   as well as for the proteic sequences
*)



(* Note: this module needs the bigarray.cma file to be compiled first *)
(* I am realizing on 06/26/2008 that the Seq module uses Alphaadn.nulla twice here
   whereas we wanted this module to be able to handle both genomic and proteic sequences
*)
module Seq = 
 struct
   external of_bigarray : (int, Bigarray.int8_unsigned_elt, Bigarray.c_layout) Bigarray.Array1.t -> 
     ('a, Bigarray.int8_unsigned_elt, Bigarray.c_layout) Bigarray.Array1.t = "%identity"


   let fill = Bigarray.Array1.fill
   
   let create n v  = let bseq = of_bigarray ( Bigarray.Array1.create Bigarray.int8_unsigned Bigarray.c_layout n) in (fill bseq v; bseq)
														   
   let length = Bigarray.Array1.dim
   
   let isnull s = ((length s)==0)
   
   let sub = Bigarray.Array1.sub
   
   
   let concat s1 s2 =
     let n1 = length s1 and n2 = length s2 in
       if((n1+n2)==0) then
	 begin
	   create 0 Alphaadn.nulla
	 end
       else
	 if(n1==0 && n2!=0) then
	 begin
	   let s= create n2 s2.{0} in
	     for i=0 to n2-1 do
	       s.{i} <- s2.{i};
	     done;
	     s;
	 end
       else
	 begin
	   let s = create (n1+n2) s1.{0} in
	     for i=1 to n1-1 do
	       s.{i} <- s1.{i};
	     done;
	     for i=0 to n2-1 do
	       s.{i+n1} <- s2.{i};
	     done;
	     s;
	 end
       
   
   
   let of_array arr = 
     try
       (let n = Array.length arr in
	let s = create n arr.(0) in 
	  for i=1 to n-1 do
	    s.{i} <- arr.(i);
	  done;
	  s
       )
     with
	 Invalid_argument _ -> create 0 Alphaadn.nulla


   let to_array bseq = 
     let n = length bseq in
     if(n==0) then 
       [||]
     else
       try
	 let arr = Array.create n bseq.{0} in 
	   for i=0 to n-1 do
	     arr.(i) <- bseq.{i};
	   done;
	   arr;
       with
	 |(Out_of_memory | Stack_overflow) -> Common.print_log "The value of the variable Sys.max_array_length is not big enough for what you need to do,\nyou need to change the architecture of your machine\n"; [||]
   
   let of_list lst = of_array (Array.of_list lst)

   (* lorsque cela le permet, cad si n est < 4194303, mais doit-on le faire ainsi? *)
   let to_list bseq = Array.to_list (to_array bseq)
 end
 


(*
module Seq :
  sig
    external of_bigarray :
      (int, Bigarray.int8_unsigned_elt, Bigarray.c_layout) Bigarray.Array1.t ->
      ('a, Bigarray.int8_unsigned_elt, Bigarray.c_layout) Bigarray.Array1.t
      = "%identity"
    val fill : ('a, 'b, 'c) Bigarray.Array1.t -> 'a -> unit
    val create :
      int ->
      'a ->
      ('a, Bigarray.int8_unsigned_elt, Bigarray.c_layout) Bigarray.Array1.t
    val length : ('a, 'b, 'c) Bigarray.Array1.t -> int
    val isnull : ('a, 'b, 'c) Bigarray.Array1.t -> bool
    val sub :
      ('a, 'b, 'c) Bigarray.Array1.t ->
      int -> int -> ('a, 'b, 'c) Bigarray.Array1.t
    val of_array :
      'a array ->
      ('a, Bigarray.int8_unsigned_elt, Bigarray.c_layout) Bigarray.Array1.t
    val to_array : ('a, 'b, 'c) Bigarray.Array1.t -> 'a array
    val to_list : ('a, 'b, 'c) Bigarray.Array1.t -> 'a list
  end

# let t = Seq.create 10 A
val t :
  (dna_t, Bigarray.int8_unsigned_elt, Bigarray.c_layout) Bigarray.Array1.t =
  <abstr>
# t.{7}
- : dna_t = A
#

*)


(* Deux fonctions du module bigseq "error prone" pour Franck
(* lorsque cela le permet, cad si n est < 4194303, mais doit-on le faire ainsi? *)
   let to_array bseq = let n = length bseq and v = bseq.{0} in
  try
  let arr = Array.create n v in 
  for i=0 to n-1 do
  arr.(i) <- bseq.{i};
  done;
  arr;
  with
  |(Out_of_memory | Stack_overflow) -> output_string stdout "la valeur de la variable Sys.max_array_length n'est pas assez importante\n"; [||]
  

   (* lorsque cela le permet, cad si n est < 4194303, mais doit-on le faire ainsi? *)
  let to_list bseq = Array.to_list (to_array bseq)
*)
