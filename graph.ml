(* ========================================================================= 
   $Author: Matthieu MANCENY
   $Mail: mmanceny@lami.univ-evry.fr
   Modification : Sarah Djebali (pour adapter au pb de la combinaison de sources)
  ======================================================================== *)

(* Librairie permettant de gérer les graphes 

   ATTENTION les sommets sont obligatoirement 
labellés avec un type variant 
(cf http://caml.inria.fr/ocaml/htmlman/manual006.html#htoc35 
et http://pauillac.inria.fr/~remy/cours/appsem/ocaml051.html)

*)

(* ========================================================================= *)
(* =================================   TYPE   ============================== *)
(* ========================================================================= *)
open Donnees


(* Type utilisé pour donner un cout aux arcs et aux sommets *)
type cost = Extension  | Inclusion | Overlap | Null 
(* noir ext estgenes | noir incl estgenes | vert prot-prot | Null *)


(* Type utilisé pour l'étiquetage des sommets *)
type vertlab = Arn of Bluspmm.t | Protn of Bluspmm.t | Nullv


(* ========================================================================= *)
(* ============================   SIGNATURE   ============================== *)
(* ========================================================================= *)

(* SIGNATURE *)
module type EDGETANK =
  sig
    type 'a t
    val make : int -> int -> 'a -> 'a t
    val set : 'a t -> int -> int -> 'a -> unit
    val get : 'a t -> int -> int -> 'a
    val rem : 'a t -> int -> int -> unit
    val null : 'a t -> 'a
    module Pred :
      sig
        val get : 'a t -> int -> int list
	val iter : (int -> int -> 'a -> unit) -> 'a t -> int -> unit
	val fold : (int -> int -> 'a -> 'b -> 'b) -> 'a t -> int -> 'b -> 'b
      end
    module Succ :
      sig
        val get : 'a t -> int -> int list
	val iter : (int -> int -> 'a -> unit) -> 'a t -> int -> unit
	val fold : (int -> int -> 'a -> 'b -> 'b) -> 'a t -> int -> 'b -> 'b
      end
  end


(* ========================================================================= *)
(* ===============================   FONCTEUR   ============================ *)
(* ========================================================================= *)

(* FONCTEUR
   prend en argument un module dont la signature 
   est identique à celle de EDGETANK *)
module Make = 
  functor (Edges_tank:EDGETANK) ->
    struct

      (*'vertex est le format du sommet, `M1 of int par exemple *)
      type ('vertex) t = 
      {
	(* dernière valeur d'index attribuée *) 
	mutable index :int;

	(* table d'association entre le label d'un sommet et son index *)
	vertices : ('vertex,int)  Hashtbl.t;

	(* table d'association entre l'index d'un sommet et son cout *)
	vertexcost : (int, cost) Hashtbl.t;

	(* arcs *)
	mutable edges :  (cost) Edges_tank.t; 
      }

      (*ratio des tables de hachage*)
      let ratio_hash = 1
			 
      (* initialisation d'un graphe *)
      let make taille v  =
	{
	  index = 0 ;
	  (* les sommets *)
	  vertices = Hashtbl.create (max ratio_hash (taille/ratio_hash));
	  (* couts des sommets *)
	  vertexcost = Hashtbl.create (max ratio_hash (taille/ratio_hash));
	  (* les arcs *)
	  edges = Edges_tank.make taille taille v
	}
  
      let create = make 


      (* ********************   SOUS MODULE VERTEX   ******************* *)

      (* module spécialisé aux traitements sur les sommets *) 
      module Vertex = struct

	(* nombre de sommets dans le graphe g *)
	let num g = g.index

	(* vérifie si un sommet v, label, est mémorisé *)
	let mem graph v  = Hashtbl.mem graph.vertices v

	(* vérifie si un sommet v, index, est mémorisé *)
	let memI graph i  = Hashtbl.mem graph.vertexcost i

	(* retourne l'index du sommet v dans le graphe graph.
	   Si ce sommet n'existe pas 
	   alors l'exception Not_found est declenchee *)
	let index graph v = Hashtbl.find graph.vertices v

	(* retourne un vecteur, où la case i contient le label
	   du sommet d'index i *)
        let label graph null = 
	  let lab = Array.make graph.index null in
	  Hashtbl.iter 
	    (fun label index -> ( lab.(index) <- label) ) 
	    graph.vertices ;
	  lab

	(* cout du sommet v, v est le label du sommet*)
	let get graph v = Hashtbl.find graph.vertexcost (index graph v)

	(* cout du sommet i, i est un index de sommet *)
	let getI graph i = Hashtbl.find graph.vertexcost i

	(* liste [(label sommet, cout sommet)] *)
	let elements graph =
	  List.sort 
	    (Pervasives.compare) 
	    (Hashtbl.fold 
	       (fun label index l -> (label,getI graph index)::l) 
	       graph.vertices [])

	(* liste [(index sommet, cout sommet)] ordonnée *)
	let elementsI graph =
	  List.sort 
	    (Pervasives.compare) 
	    (Hashtbl.fold 
	       (fun index cost l -> (index,cost)::l) 
	       graph.vertexcost 
	       [])

	(* transforme une liste [(sommets indexés,cout)]
	   en liste [(sommets labellés, cout)] *)
	let elementsItoL graph l null =
	  let lab = label graph null in
	  List.map
	    (fun (index,cost) -> (lab.(index),cost)) 
	    l

	(* ajoute le sommet v au graphe graph si celui ci n'existe pas déjà 
	   renvoie un booleen pour dire si le sommet est créé*)
	let set graph v c=
	  try 
	    let i = index graph v in (false,i)  
	  with Not_found ->
	    let i =  graph.index in 
	      Hashtbl.add graph.vertices v i;
              Hashtbl.add graph.vertexcost i c;
	      graph.index <- graph.index+1 ; 
	      (true,i) 

	(* supprime le sommet v, label, du graphe g *)
	let remove graph v =
	  try
	    let i0 = index graph v in
	    Hashtbl.remove graph.vertices v ;
            Hashtbl.remove graph.vertexcost i0;
	    (*suppression des arcs associés au sommet v*)
	    Edges_tank.Pred.iter
	      (fun i1 i2 d12 -> Edges_tank.rem graph.edges i1 i2)
	      graph.edges
	      i0;
	    Edges_tank.Succ.iter
	      (fun i1 i2 d12 -> Edges_tank.rem graph.edges i1 i2)
	      graph.edges
	      i0;
	  with
	    Not_found -> ()

	(* iter f sur l'ensemble des sommets de graph 
	   f: int -> 'a -> unit *)
	let iter f graph = 
	  Hashtbl.iter 
	    (fun lab ind -> f ind (Hashtbl.find graph.vertexcost ind) ) 
	    graph.vertices
		  
	(* fold sur les sommets déclarés 
	   f: i -> ci -> init -> init *)
	let fold f graph init = 
	  Hashtbl.fold 
	    (fun lab ind l -> f ind (Hashtbl.find graph.vertexcost ind) l )
	    graph.vertices 
	    init

	(* Partie de traitement sur les sommets  - selon leur format interne *)
	let selecti g p = 
	  let a = Array.create g.index false in  
	  Hashtbl.iter (fun k i -> if p k then a.(i) <-true) g.vertices ; a

      end (*module Vertex*)


      (* ********************   SOUS MODULE EDGE   ********************* *)

      (* module spécialisé aux arcs *)
      module Edge = struct

	(* récupère la valeur d'un arc ; travail sur les labels
	   renvoie la valeur associée à null si l'arc n'existe pas *)
	let get graph v1 v2 = 
	  Edges_tank.get 
	    graph.edges 
	    (Vertex.index graph v1) (Vertex.index graph v2)

	(* récupère la valeur d'un arc ; travail sur les index
	   renvoie la valeur associée à null si l'arc n'existe pas *)
	let getI graph i j = 
	  Edges_tank.get graph.edges i j

	(* teste l'existence de l'arc v1->v2 ; travail sur les labels *)
	let exist graph v1 v2 = 
	  not( (get graph v1 v2) = Edges_tank.null graph.edges)

	(* teste l'existence de l'arc i-> j ; travail sur les index *)
	let existI graph i j = 
	  not( (getI graph i j) = Edges_tank.null graph.edges )

	(* fold f sur l'ensemble des arcs de graph 
	   f : (i1 i2 d12) (type de init) -> (type de init) *)
	let fold f graph init=
	  Hashtbl.fold
	    (fun lab ind l ->
	       Edges_tank.Pred.fold
	         f
	         graph.edges
	         ind
	         l)
	    graph.vertices 
	    init

	(* iter f sur l'ensemble des arcs de graph 
	   f : (i1 -> i2 -> d12) -> unit *)
	let iter f graph =
	  Hashtbl.iter
	    (fun lab ind ->
	      Edges_tank.Pred.iter
	        f
	        graph.edges
	        ind)
	    graph.vertices

	(* liste des arcs. La sortie est sous forme de liste [(i1,i2,d12)] *) 
	let elementsI graph = 
	  fold
	    (fun i1 i2 d12 l -> 
	       if (existI graph i1 i2) then (i1,i2,d12)::l else l) 
	    graph  
	    []

	(* transforme l, liste d'arcs indexés en liste d'arcs labellés *)
	let elementsItoL graph l null =
	  let lab = Vertex.label graph null in
	  List.map
	    ( fun (i,j,dij) -> (lab.(i),lab.(j),dij) ) 
	    l

	(* liste des arcs. La sortie est sous forme de liste [(v1,v2,d12)] *) 	
	let elements graph null = elementsItoL graph (elementsI graph) null;;

	(* ajoute l'arc v1->v2 si les deux sommets existent
	   et que l'arc n'existe pas déjà.
	   renvoie un booleen indiquant si l'arc a été créé *)
	let set graph v1 v2 d12 =
	  if ( Vertex.mem graph v1 && Vertex.mem graph v2 && not(exist graph v1 v2))
	  then 
	    (Edges_tank.set 
	       graph.edges 
	       (Vertex.index graph v1) (Vertex.index graph v2) d12; 
	     true)
          else (false)

	(* Ajoute l'arc i1->i2 de cout d12 ;
	   si un des index est trop grand, exception Invalid_argument*)
	let setI graph i1 i2 d12 =
	  if (i1<graph.index)&&(i2<graph.index)
	  then Edges_tank.set graph.edges i1 i2 d12
	  else invalid_arg "Index sommet trop grand, Graph.Edge.addI"

	(* supprime l'arc v1->v2 de g, travail sur les labels *)
	let remove graph v1 v2 = 
	  try
	    let i1 = Vertex.index graph v1 in
	    let i2 = Vertex.index graph v2 in
	    Edges_tank.rem graph.edges i1 i2
	  with
	    Not_found -> invalid_arg "Sommet inexistant, Graph.Edge.remove"

	(* supprime l'arc i->j du graphe, travail sur les index *)
	let removeI graph i j = 
	  let b=((i<graph.index)&&(j<graph.index)) in
	    if b 
	    then Edges_tank.rem graph.edges i j
	    else invalid_arg "Index sommet trop grand Graph.Edge.removeI"

      end(*module Edge*)


      (* ********************   SOUS MODULE PRED   ********************* *)

      (* graphe traité selon la relation prédécesseur *) 
      module Pred = struct

	(* liste des prédecesseurs *)
	let get graph v  = 
	  Edges_tank.Pred.get graph.edges (Hashtbl.find graph.vertices v)
	let getI graph j  = Edges_tank.Pred.get graph.edges j
			  
	(* nombre de prédecesseurs*)
	let num graph v =  List.length (get graph v)
	let numI graph j = List.length (getI graph j)

	(*p travaille sur les index*)
	let forall p graph v = List.for_all p (get graph v)
	let exists p graph v = List.exists p (get graph v)

	(* iter de f sur les prédécesseurs de j 
	   f: int -> 'a -> unit *)
	let iterverticesI f graph j = 
	  List.iter 
	    (fun i -> f i (Hashtbl.find graph.vertexcost i))
	    (getI graph j)

	(* fold f sur les prédécesseurs de j *)
	let foldverticesI f graph j init = 
	  List.fold_right 
	    (fun i -> f i (Hashtbl.find graph.vertexcost i) init) 
	    (getI graph j) 
	    init

	(* iter de f sur les arcs entrant de j, i->j *)
	let iteredgesI f graph j = 
	  Edges_tank.Pred.iter f graph.edges j

	(* fold f sur les arcs entrant de j, i->j *)
	let foldedgesI f graph j init = 
	  Edges_tank.Pred.fold f graph.edges j init

      end (*module Pred*)


      (* ********************   SOUS MODULE SUCC   ********************* *)

      (* graphe traité selon la relation successeur *) 
      module Succ = struct

	(* liste des successeurs *)
	let get graph v  = 
	  Edges_tank.Succ.get graph.edges (Hashtbl.find graph.vertices v)
	let getI graph i  = Edges_tank.Succ.get graph.edges i

        (* nombre de successeurs *)
	let num graph v =  List.length (get graph v)
	let numI graph i = List.length (getI graph i)

	let forall p graph v = List.for_all p (get graph v)
	let exists p graph v = List.exists p (get graph v)

	(* iter f sur les successeurs de i 
	   f: int -> 'a -> unit *)
	let iterverticesI f graph i = 
	  List.iter 
	    (fun j -> f j (Hashtbl.find graph.vertexcost j)) 
	    (getI graph i)

	(* fold f sur les prédécesseurs de j *)
	let foldverticesI f graph i init = 
	  List.fold_right 
	    (fun j -> f i j (Hashtbl.find graph.vertexcost j))
	    (getI graph i) 
	    init

	(* iter f sur les arcs i->j *)
	let iteredgesI f graph i = 
	  Edges_tank.Succ.iter f graph.edges i

	(* fold f sur les arcs i->j *)
	let foldedgesI f graph i init =
	  Edges_tank.Succ.fold f graph.edges i init

      end (*module Succ*)


      (* ************************   FONCTIONS   ************************ *)

      (* TRI TOPOLOGIQUE *)
		      
      (* liste des index des sommets triés par ordre topologique *)
      let toposortkernel g =
	begin
	
	  (* création de la liste des indices sommets *)
	  let order = ref [] in

	  (* tableau du nombre de successeurs *)
	  let succ = Array.make g.index (-1) in
	  Edge.iter
	    (fun i j dij -> 
	       if ( succ.(i) = (-1) )
	       then (succ.(i) <- 1)
	       else (succ.(i) <- (succ.(i)+1));
	       
	       if succ.(j) = (-1) then succ.(j)<-0
	    )
	    g;

          Array.iteri
	    (fun i succi -> 	    
	       if (succ.(i) = 0) 
	       then (order:=i::!order)
	    )
	    succ;

	  let rec subtopo s=
	    match s with
		[] -> []
	      | s -> begin
		  let new_s = ref [] in
		    List.iter
 		      (	fun j ->
 			  List.iter
 			  ( fun i ->
		 	      (succ.(i)<-succ.(i)-1;
			       if (succ.(i)=0)
			       then (order:=i::!order; new_s := i::!new_s)
			       else ()
 			      )
	 		  )
			  (Pred.getI g j)
		      )
 		      s;
 		    subtopo !new_s;
		end
	  in
	  subtopo !order;
	  !order;
	end (*toposortkernel*)


      (* liste des labels des sommets triés par ordre topologique *)
      let toposort g null =
	let lab = Vertex.label g null in
 	List.map 
 	  (fun i -> lab.(i)) 
	  (toposortkernel g)

  
      

    end (*module Graph*)


(* ========================================================================= *)
(* =======================   IMPLÉMENTATION DE EDGE_TANK   ================= *)
(* ====================   SOUS FORME DE MATRICE D'ADJACENCE   ============== *)
(* ========================================================================= *)

(* Structure implantant la matrice d'adjacence comportant les arcs *)
module Mat = struct

  type 'a t = { tank:'a array array; null : 'a} 

  let make i j v = {tank=Array.make_matrix i j v; null=v}

  (* création de l'arc i->j de distance v *)
  let set t i j v = t.tank.(i).(j) <- v

  (* distance de l'arc i->j *)
  let get t i j = t.tank.(i).(j)

  (* renvoie la valeur null *)
  let null t = t.null

  (* remise à 0 de la case i, j *)
  let rem t i j =
    t.tank.(i).(j) <- t.null

  module Pred = struct

    (* renvoie la liste des éléments, non nuls, de la ligne colonne j
       ie les prédecesseurs de j *)
    let get t j = 
      let len = Array.length t.tank  in 
      let rec sub i =
	if i=len
	then []
 	else (
	  if t.tank.(i).(j) <> t.null
	  then i::(sub (i+1))
	  else (sub (i+1)) )
      in sub 0

    (* iter de f sur les arcs entrant de j 
       f: int -> int -> 'a -> unit *)
    let iter f t j = 
      List.iter 
	(fun i -> f i j t.tank.(i).(j)) 
	(get t j)

    (* fold de f sur les arcs entrant de j
       f: int -> int -> 'a -> 'b -> 'b *)
    let fold f t j init = 
      List.fold_right 
	(fun i l -> f i j t.tank.(i).(j) l) 
	(get t j) 
	init

  end (* module Pred *)

  module Succ = struct
    (* renvoie la liste des éléments, non nuls, de la ligne i
       ie les successeurs de i *)
    let get t i =
      let len = Array.length t.tank.(i)  in
      let rec sub  j = 
	if j=len then 
	  []
	else
	  if (t.tank.(i).(j) <> t.null) then j::sub (j+1) 
          else sub (j+1)
      in sub  0

    (* iter de f sur les arcs sortant de i
       f: int -> int -> 'a -> unit *)
    let iter f t i = 
      let a = Array.mapi (fun ind elt -> (ind, elt)) t.tank.(i) in
      Array.iter 
	(fun (j,dij) -> if (dij = t.null) then () else (f i j dij)) 
	a

    (* fold de f sur les arcs sortant de i
       f: int -> int -> 'a -> 'b -> 'b *)
    let fold f t i init = 
      let a = Array.mapi (fun ind elt -> (ind, elt)) t.tank.(i) in
      Array.fold_right 
	(fun (j,dij) l -> if (dij = t.null) then (l) else (f i j dij l)) 
	a 
	init

  end (* module Succ *)


end (* module Mat *)


(* ========================================================================= *)
(* =======================   IMPLÉMENTATION DE EDGE_TANK   ================= *)
(* ==============   SOUS FORME DE TABLEAU DE TABLES DE HACHAGE   =========== *)
(* ========================================================================= *)

(* Structure implantant le vecteur de tables de hachage 
   la case i du vecteur contient une table de hachage 
   des prédécesseurs de i, sous la forme 
   ( pred(i) , distance de l'arc (pred(i) -> i) ) *)
module HashtblArray = struct

  type 'a t = { tank : (int, 'a ref) Hashtbl.t array ; null : 'a} 

  let make i j v = {tank=Array.init i (fun i -> Hashtbl.create 0); null=v}

  (* création de l'arc i->j de distance dij *)
  let set t i j dij = Hashtbl.add t.tank.(j) i (ref dij)

  (* distance de l'arc i->j *)
  let get t i j = 
    try !(Hashtbl.find t.tank.(j) i)
    with Not_found -> t.null

  (* renvoie la valeur null *)
  let null t = t.null

  (* suppression de l'arc i, j *)
  let rem t i j =
    Hashtbl.remove t.tank.(j) i


  module Pred = struct

    (*liste des index des prédécesseurs de j *)
    let get t j = 
      Hashtbl.fold 
	(fun k d l -> k::l)
	t.tank.(j)
	[]

    (* iter de f sur les arcs entrant de j
       f: int -> int -> 'a -> unit *)
    let iter f t j =
      Hashtbl.iter
	(fun i dij -> f i j (!dij)) 
	t.tank.(j)

    (* fold de f sur les arcs entrant de j
       f: int -> int -> 'a -> 'b -> 'b*)
    let fold f t j init = 
      Hashtbl.fold 
	(fun i dij l -> f i j (!dij) l) 
	t.tank.(j) 
	init

  end (*module Pred*)

  module Succ = struct 

    (* renvoie la liste des index des successeurs de i *)
    let get t i = 
      let l = ref [] in
	for j = 0 to (Array.length t.tank - 1) do
	  if (Hashtbl.mem t.tank.(j) i)
	  then ( l := j::!l )
	  else ()
	done;
	!l

    (* iter de f sur les arcs entrant de j
       f: int -> int -> 'a -> unit *)
    let iter f t i = 
      List.iter 
	(fun j -> f i j !(Hashtbl.find t.tank.(j) i)) 
	(get t i)

    (* fold de f sur les arcs sortant de i *)
    let fold f t i init = 
      List.fold_right 
	(fun j l -> f i j !(Hashtbl.find t.tank.(j) i) l) 
	(get t i) 
	init

  end (*module Succ*)

end (*module HashtblArray*)



(*
Un essai :
---------
#open Graph;;
# #load "graph.cmo";;                                                                                                # module Mg = Make(HashtblArray);;                                                                                    module Mg :                                                                                                                                                                 sig                                                                                                                                                                            type 'a t =                                                                                                                                                                    'a Graph.Make(Graph.HashtblArray).t = {                                                                                                                                      mutable index : int;                                                                                                                                                         vertices : ('a, int) Hashtbl.t;                                                                                                                                              vertexcost : (int, Graph.cost) Hashtbl.t;                                                                                                                               mutable edges : Graph.cost Graph.HashtblArray.t;                                                                                                                      }                                                                                                                                                                            val ratio_hash : int                                                                                                                                                         val make : int -> Graph.cost -> 'a t                                                                                                                                    val create : int -> Graph.cost -> 'a t                                                                                                                                  module Vertex :                                                                                                                                                                sig                                                                                                                                                                            val num : 'a t -> int                                                                                                                                                        val mem : 'a t -> 'a -> bool                                                                                                                                                 val memI : 'a t -> int -> bool                                                                                                                                               val index : 'a t -> 'a -> int                                                                                                                                                val label : 'a t -> 'a -> 'a array                                                                                                                                           val get : 'a t -> 'a -> Graph.cost                                                                                                                                      val getI : 'a t -> int -> Graph.cost                                                                                                                                    val elements : 'a t -> ('a * Graph.cost) list                                                                                                                           val elementsI : 'a t -> (int * Graph.cost) list                                                                                                                         val elementsItoL : 'a t -> (int * 'b) list -> 'a -> ('a * 'b) list                                                                                                           val set : 'a t -> 'a -> Graph.cost -> bool * int                                                                                                                        val remove : 'a t -> 'a -> unit                                                                                                                                              val iter : (int -> Graph.cost -> unit) -> 'a t -> unit                                                                                                                  val fold : (int -> Graph.cost -> 'a -> 'a) -> 'b t -> 'a -> 'a                                                                                                          val selecti : 'a t -> ('a -> bool) -> bool array                                                                                                                           end                                                                                                                                                                        module Edge :                                                                                                                                                                  sig                                                                                                                                                                            val get : 'a t -> 'a -> 'a -> Graph.cost                                                                                                                                val getI : 'a t -> int -> int -> Graph.cost                                                                                                                             val exist : 'a t -> 'a -> 'a -> bool                                                                                                                                         val existI : 'a t -> int -> int -> bool                                                                                                                                      val fold :                                                                                                                                                                     (int -> int -> Graph.cost -> 'a -> 'a) -> 'b t -> 'a -> 'a
        val iter : (int -> int -> Graph.cost -> unit) -> 'a t -> unit
        val elementsI : 'a t -> (int * int * Graph.cost) list
        val elementsItoL :
          'a t -> (int * int * 'b) list -> 'a -> ('a * 'a * 'b) list
        val elements : 'a t -> 'a -> ('a * 'a * Graph.cost) list
        val set : 'a t -> 'a -> 'a -> Graph.cost -> bool
        val setI : 'a t -> int -> int -> Graph.cost -> unit
        val remove : 'a t -> 'a -> 'a -> unit
        val removeI : 'a t -> int -> int -> unit
      end
    module Pred :
      sig
        val get : 'a t -> 'a -> int list
        val getI : 'a t -> int -> int list
        val num : 'a t -> 'a -> int
        val numI : 'a t -> int -> int
        val forall : (int -> bool) -> 'a t -> 'a -> bool
        val exists : (int -> bool) -> 'a t -> 'a -> bool
        val iterverticesI :
          (int -> Graph.cost -> unit) -> 'a t -> int -> unit
        val foldverticesI :
          (int -> Graph.cost -> 'a -> 'a -> 'a) ->
          'b t -> int -> 'a -> 'a
        val iteredgesI :
          (int -> int -> Graph.cost -> unit) -> 'a t -> int -> unit
        val foldedgesI :
          (int -> int -> Graph.cost -> 'a -> 'a) ->
          'b t -> int -> 'a -> 'a
      end
    module Succ :
      sig
        val get : 'a t -> 'a -> int list
        val getI : 'a t -> int -> int list
        val num : 'a t -> 'a -> int
        val numI : 'a t -> int -> int
        val forall : (int -> bool) -> 'a t -> 'a -> bool
        val exists : (int -> bool) -> 'a t -> 'a -> bool
        val iterverticesI :
          (int -> Graph.cost -> unit) -> 'a t -> int -> unit
        val foldverticesI :
          (int -> int -> Graph.cost -> 'a -> 'a) ->
          'b t -> int -> 'a -> 'a
        val iteredgesI :
          (int -> int -> Graph.cost -> unit) -> 'a t -> int -> unit
        val foldedgesI :
          (int -> int -> Graph.cost -> 'a -> 'a) ->
          'b t -> int -> 'a -> 'a
      end
    val toposortkernel : 'a t -> int list
    val toposort : 'a t -> 'a -> 'a list
  end

# let g = Mg.make 5 Extension;;
val g : '_a Mg.t =
  {Mg.index = 0; Mg.vertices = <abstr>; Mg.vertexcost = <abstr>;
   Mg.edges =
    {Graph.HashtblArray.tank =
      [|<abstr>; <abstr>; <abstr>; <abstr>; <abstr>|];
     Graph.HashtblArray.null = Extension}}
#let b1 = {hl = []; il= []; egh = (0,0); egi = (0,0); st = Forward};;                                                val b1 : Graph.blusp =                                                                                               {hl = []; il = []; egh = (0, 0); egi = (0, 0); st = Forward}
#let v1 = Arn b1;;
val v1 : Graph.vertices =
  Arn {hl = []; il = []; egh = (0, 0); egi = (0, 0); st = Forward}
# Mg.Vertex.set g v1 Inclusion;;
- : bool * int = (true, 0)
Mais on n'a pas besoin de cout pour les sommets!
# Mg.Vertex.index g v1;;
- : int = 0
#  Mg.Vertex.set g v1 Extension;;
- : bool * int = (false, 0)
ici il n'a pas créé de nouveau sommet car v1 existe deja dans le graphe
# Mg.Vertex.index g v1;;
- : int = 0
idem
# let b2 = {hl = []; il = []; egh =(1,0); egi = (0, 0); st = Forward};;
val b2 : Graph.blusp =
  {hl = []; il = []; egh = (1, 0); egi = (0, 0); st = Forward}
# let v2 = Prot b2;;
val v2 : Graph.vertices =
  Prot {hl = []; il = []; egh = (1, 0); egi = (0, 0); st = Forward}
# Mg.Vertex.set  g v2 Inclusion;;
- : bool * int = (true, 1)
ici il a créeé le nouveau sommet v2 car différent de v1
# Mg.Vertex.index g v2;;
- : int = 1
idem (son index est 1)
#  Mg.Vertex.num g;;
- : int = 2
il ya deux sommets dans le graphe
# Mg.Vertex.mem g v1;;
- : bool = true
# Mg.Vertex.mem g v2;;
- : bool = true
v1 et v2 sont bien dans le graphe (v1 et v2 sont des labels)
# Mg.Vertex.remove g v1;;
- : unit = ()
On supprime le sommet de label v1
# Mg.Vertex.index g v1;;
Exception: Not_found.
Il n'y est plus
# Mg.Vertex.set g v1 Extension;;
- : bool * int = (true, 2)
# Mg.Edge.set g v1 v2 Extension;;
- : bool = true
#
*)
