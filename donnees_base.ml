(* donnees_base.ml : donnees utiles aux listes et de base comme les hsps ou les molécules *)


open Seq
open Preserveorf
open Alphaprot
open Alphaadn

 

let code_gen t3nt =
  match t3nt with
    |[|A; A; A|] -> Alphaprot.K;
    |[|A; A; T|] -> Alphaprot.N;
    |[|A; A; G|] -> Alphaprot.K;
    |[|A; A; C|] -> Alphaprot.N;
    |[|A; T; A|] -> Alphaprot.I;
    |[|A; T; T|] -> Alphaprot.I;
    |[|A; T; G|] -> Alphaprot.M;
    |[|A; T; C|] -> Alphaprot.I;
    |[|A; G; A|] -> Alphaprot.R;
    |[|A; G; T|] -> Alphaprot.S;
    |[|A; G; G|] -> Alphaprot.R;
    |[|A; G; C|] -> Alphaprot.S;
    |[|A; C; A|] -> Alphaprot.T;
    |[|A; C; T|] -> Alphaprot.T;
    |[|A; C; G|] -> Alphaprot.T;
    |[|A; C; C|] -> Alphaprot.T;
    |[|T; A; A|] -> Alphaprot.Z;
    |[|T; A; T|] -> Alphaprot.Y;
    |[|T; A; G|] -> Alphaprot.Z;
    |[|T; A; C|] -> Alphaprot.Y;
    |[|T; T; A|] -> Alphaprot.L;
    |[|T; T; T|] -> Alphaprot.F;
    |[|T; T; G|] -> Alphaprot.L;
    |[|T; T; C|] -> Alphaprot.F;
    |[|T; G; A|] -> Alphaprot.Z;
    |[|T; G; T|] -> Alphaprot.C;
    |[|T; G; G|] -> Alphaprot.W;
    |[|T; G; C|] -> Alphaprot.C;
    |[|T; C; A|] -> Alphaprot.S;
    |[|T; C; T|] -> Alphaprot.S;
    |[|T; C; G|] -> Alphaprot.S;
    |[|T; C; C|] -> Alphaprot.S;
    |[|G; A; A|] -> Alphaprot.E;
    |[|G; A; T|] -> Alphaprot.D;
    |[|G; A; G|] -> Alphaprot.E;
    |[|G; A; C|] -> Alphaprot.D;
    |[|G; T; A|] -> Alphaprot.V;
    |[|G; T; T|] -> Alphaprot.V;
    |[|G; T; G|] -> Alphaprot.V;
    |[|G; T; C|] -> Alphaprot.V;
    |[|G; G; A|] -> Alphaprot.G;
    |[|G; G; T|] -> Alphaprot.G;
    |[|G; G; G|] -> Alphaprot.G;
    |[|G; G; C|] -> Alphaprot.G;
    |[|G; C; A|] -> Alphaprot.A;
    |[|G; C; T|] -> Alphaprot.A;
    |[|G; C; G|] -> Alphaprot.A;
    |[|G; C; C|] -> Alphaprot.A;
    |[|C; A; A|] -> Alphaprot.Q;
    |[|C; A; T|] -> Alphaprot.H;
    |[|C; A; G|] -> Alphaprot.Q;
    |[|C; A; C|] -> Alphaprot.H;
    |[|C; T; A|] -> Alphaprot.L;
    |[|C; T; T|] -> Alphaprot.L;
    |[|C; T; G|] -> Alphaprot.L;
    |[|C; T; C|] -> Alphaprot.L;
    |[|C; G; A|] -> Alphaprot.R;
    |[|C; G; T|] -> Alphaprot.R;
    |[|C; G; G|] -> Alphaprot.R;
    |[|C; G; C|] -> Alphaprot.R;
    |[|C; C; A|] -> Alphaprot.P;
    |[|C; C; T|] -> Alphaprot.P;
    |[|C; C; G|] -> Alphaprot.P;
    |[|C; C; C|] -> Alphaprot.P;
    |[|N; _; _|] -> Alphaprot.X; (* X représente un acide aminé hypothétique qui serait la traduction d'un codon comportant au moins un N *)
    |[|_; N; _|] -> Alphaprot.X;
    |[|_; _; N|] -> Alphaprot.X;
    | _ -> raise (Invalid_argument "code_gen");;


(* split_to_codon prend en entreé une liste de nucléotides lnt (mod 3) et renvoit une liste de tableaux
   de 3 nucléotides représentant les différents codons de cette liste *)
let rec split_to_codon lnt =
  match lnt with
    |n1::n2::n3::q -> [|n1;n2;n3|]::(split_to_codon q);
    |[] -> [];
    |_ -> raise (Invalid_argument "split_to_codon");;


(* translate traduit un mot génomique de taille modulo 3 en un mot protéique.
   Attention gword est une liste et on obtient en sortie une liste *)
let translate gsword =
  try 
    List.map code_gen (split_to_codon gsword)
  with
    |Invalid_argument s -> raise (Invalid_argument s); [];;



(* listes de segseq *)
let met = List.map Seq.of_array [[|A;T;G|]];;
let accept = List.map Seq.of_array [[|A;G|]];;
let donnor = List.map Seq.of_array [[|G;T|]];;
let acceptplus = List.map Seq.of_array [[|A;G|];[|A;T|]];;
let donnorplus = List.map Seq.of_array [[|G;T|];[|G;C|]];;


(* nouveau type, énuméré, pour étiqueter un eclusp selon sa position dans le blusp ou futur modèle de gène *)
type label = Unique | Initial | Internal | Terminal;;


type mol = 
    | Prot 
    | Rna 
    | Null

type molec = 
    | P of string 
    | A of string 
    | N


module Molecule = 
struct
  type t = molec                     (* a molecule identifier, its key, this is a string 
					(its name) preceeded by its type *)
  let compare = Pervasives.compare   (* for comparing the names of two molecules *)
  let create m  = m                  (* for creating a new molecule *)
  let make = create
  let null = N
end

module MolSet = Set.Make(Molecule)


type seqa = 
    {
      noma : string;         (* DNA sequence name *)
      lga : int;             (* DNA sequence length *)
      sa : dna_t array       (* array containing the nucleotides *)
    }     
type seqp = 
    {
      nomprot : Molecule.t;  (* molecule name *)
      lgp : int;             (* molecule length *)
      sp : prot_t array      (* array containing the amino acids *)
    }

type brinfr = | Forward  
              | Reverse
          


(* Hsp module: handles hsps.
   An hsp is a piece of alignment (for instance from the program blat) 
   between a genomic sequence that we want to annotate, and another molecule, 
   typically a protein or an mrna. *)
module Hsp = 
struct
  type t = 
      {
	tmol : mol;        (* molecule type (Prot or Rna) *)
	nomp : Molecule.t; (* molecule name *) 
	score : float;     (* hsp score, float replaces int on 03/13/06 *)   
	debg : int;        (* hsp genomic start *)
	fing : int;        (* hsp genomic end *)
	debp : int;        (* hsp start on the molecule *)
	finp : int;        (* hsp end on the molecule *)
	brinh : brinfr;    (* hsp strand *)
	nbntfus : int;     (* number of genomic nucleotides that were fused during the normalisation step 
			      this information is useful for eliminating monoproteic blusps that go over the
			      pm threshold due to too many fusions. Remark: if nt<0 we add -nt *)
	cutdwn : bool;     (* has monomolecule blusping seen a breaking necesary on the up direction of this hsp? *)
	aligncomp : bool;  (* says whether it is necesary to make a complentary alignment between this hsp
			      and the following one coming from the same molecule *)
	commh : string;    (* hsp comments *)
      }
                     
  let compare h1 h2 = 
    if (Pervasives.compare h1.debg h2.debg)!=0 then 
      Pervasives.compare h1.debg h2.debg   
    else
      (* here genomic starts are equal --> need to look at the genomic ends *)
      begin
	if (Pervasives.compare h1.fing h2.fing)!=0 then
	  (Pervasives.compare h1.fing h2.fing)
	else
	  (* genomic ends are also equal --> elements have same boundaries *)
	  0
      end

  let compnameprot  h1 h2 = Pervasives.compare h1.nomp h2.nomp
 
  let create tm np sc dg fg dp fp br nf cd tb cm  =  
    {
      tmol = tm;       (* molecule type (Prot or Rna) *)
      nomp = np;       (* molecule name *)
      score = sc;      (* hsp score *)   
      debg = dg;       (* hsp genomic start *)
      fing =fg;        (* hsp genomic end *)
      debp =dp;        (* hsp start on the molecule *)
      finp =fp;        (* hsp end on the molecule *)
      brinh =br;       (* hsp strand *)
      nbntfus =nf;     (* number of genomic nucleotides that were fused during the normalisation step *)
      cutdwn =cd;      (* has monomolecule blusping seen a breaking necesary on the up direction of this hsp? *)
      aligncomp = tb;  (* says whether it is necesary to make a complentary alignment between this hsp
			  and the following one coming from the same molecule *)
      commh =cm;       (* hsp comments *)
    } 
   
  (* In order to create a new hsp *)
  let make = create

  let null = create Null Molecule.null 0.0 0 0 0 0 Forward 0 false false ""
  let tmol h = h.tmol
  let nomp h = h.nomp
  let score h = h.score
  let debg h = h.debg
  let fing h = h.fing
  let debp h = h.debp
  let finp h = h.finp
  let brinh h = h.brinh
  let nbntfus h = h.nbntfus
  let cutdwn h = h.cutdwn
  let aligncomp h = h.aligncomp
  let commh h = h.commh

  let size h = (fing h) - (debg h) + 1

  let setcutdwn h = 
    {
      h with cutdwn=true
    }
  let isprot h = (tmol h = Prot)
  let isarn h = (tmol h = Rna)
	  
end



(* A partir d'un tableau d'hsps ordonnées selon la protéines fait une liste de liste d'hsps
   une par protéine (ne prend pas en compte le brin dans cette séparation)
   et ne conserve une liste d'hsp interne à cette liste que si la propriété p
   est vérifiée sur cette liste (par défaut p est vrai)
   utilisée par tclusp2bluspmt de donnees.ml 
   (pourrait être remplacée par un blusping)
*)

let partage_selon_prot p tHsp =
  let i = ref 0 and n = Array.length tHsp in
  let ltemp = ref [] and llh = ref [] in
  let premProt = ref (Hsp.nomp (tHsp.(0))) in

    while(!i<=(n-1)) do
      while (!i<=(n-1) && (!premProt  = (Hsp.nomp (tHsp.(!i))))) do
	ltemp:=tHsp.(!i)::!ltemp;
	incr(i)
      done;
      
      if(p !ltemp) then
	llh:=(!ltemp)::!llh;
	
      
      if(!i<=(n-1)) then
	begin
	  ltemp := [];
	  premProt:= Hsp.nomp (tHsp.(!i));
	end;

    done;
    !llh 


let brin_to_string br =
  match br with
    |Forward -> "+";
    |Reverse -> "-"



let intlist2string l offset =
  let s = ref "" in
    List.iter (fun pos -> s:=!s^(string_of_int (pos+offset))^",") l;
    !s


let intarray2string t offset =
	let s = ref "" in
	Array.iter (fun pos -> s:=!s^(string_of_int (pos+offset))^",") t;
	String.sub !s 0 (String.length !s - 1);;




let molec2string = function
	| P s -> if s="" then "NoneP" else s
	| A s -> if s="" then "NoneR" else s
	| _ ->  "None";;



let typemol2string = function
  |Prot -> "Protein"
  |Rna -> "mRNA"
  |Null -> "Null molecule"



let print_hsp o h =
  Printf.fprintf o "tmol = %s; nomp = %s; score = %f; debg = %i; fing = %i; debp = %i; finp = %i; brinh = %s; pbcutdwn=%s\n" (typemol2string (Hsp.tmol h)) (molec2string (Hsp.nomp h)) (Hsp.score h) (Hsp.debg h) (Hsp.fing h) (Hsp.debp h) (Hsp.finp h) (brin_to_string (Hsp.brinh h)) (if (Hsp.cutdwn h) then "true" else "false")




(* give_label est une fonction qui accorde un label à un objet de position i dans une collection 
   de n objets ordonnés (de no 0 à no n-1). Ces objets seront ici des hsps de meme alignement d'ARNm *)
let give_label i n = 
  if (i=0 && i=n-1) then Unique 
  (* normalement impossible à ce stade car on a éliminé les transcrits monohsps*)
  else
    if (i=0) then Initial
    else
      if (i=n-1) then Terminal
      else
	Internal;;

