(* collection.ml : définit l'objet SegSeq (sous forme de module) ie les séquences segmentées que sont pex les
   ensembles de clusps, de blusps et d'exons. Une séquence segmentée est un ensemble d'objets departagés en 
   sous-ensembles (=segments) d'objets (un peu comme une partition). Pour representer cette segmentation
   en sous-ensembles d'objets, l'ensemble initial  est stocké dans un tableau tank, et on se sert ensuite 
   d'une liste d'entiers ord pour indiquer les indices i dans tank tels que l'objet tank.(i) est le délimiteur
   de début d'un segment/sous-ensemble d'objets *)



module SegSeq = 
struct
  type 'a t = { 
    tank : 'a array ; (*  tableau de données *)  
    ord : int list (*  Liste d'indicateur de début dans l'ordre inverse *)
  }

  let tank s = s.tank
  let ord s = s.ord
  let make size v = { tank= Array.make size v ; ord = []}
  let make2 t = {tank=t; ord = []}
  let make3 t o = {tank=t; ord = o}
  let create = make 

  let isnull s = (s.tank = [||])

  let hd s = 
    match s.ord with 
      | [] -> failwith "SegSeq.hd"
      | i::_ -> i
	  
  exception Seg of int 
    

  let setsegment s l   = {s with ord=l}
			     
  let partition p s = 
    let rec apart deb fin = 
      try 
	for i = deb to fin do 
	  if not (p s.tank i) then raise (Seg(i))  (* p prend en entree un tableau et la pos courante *)
	done;
	[]
      with Seg i -> i::(apart (i+1) fin)
 
    in
    let rec lpart = function
      |[] -> [] 
      |[i]-> [apart i ((Array.length s.tank) -1)]
      |i::(j::r as l) -> (apart i j)::(lpart l)
    in
	
      setsegment s (List.flatten (lpart s.ord))

 
  let intervalles s =
    let rec aux l n =  
      match l with
	|[] -> [];
	|[t] -> [(t,n)];
	|t1::t2::q -> (t1,t2)::(aux (t2::q) n) in
      aux s.ord ((Array.length s.tank)-1)


  let single a = {tank = a ; ord =[]}

 (*  let fold f s e =  *) 


(* convertit une segseq en liste de tableaux *) 
  let elements s = 
    let rec lambda = function 
      |[] -> []
      |[i]-> [Array.sub s.tank i ((Array.length s.tank) -i)]
      |i::(j::r as l)-> (Array.sub s.tank i (j-i))::(lambda l)
    in
      lambda (s.ord) 
 

  (* f est une fonction qui prend un élément de type ce qu'il y a dans tank *)
  let map f s =
    make3 (Array.map f s.tank) s.ord


  (* convertit une liste de tableaux en segseq *) 
  let tosegseq lt =
    let i = ref 0 and s = ref 0 and n = List.length lt in
    let tord = Array.create n 0 in
      
      while(!i<n) do
	tord.(!i) <- !s;
	s := !s + Array.length (List.nth lt !i);
	incr i;
      done;
      
      setsegment (make2 (Array.of_list (List.flatten (List.map Array.to_list lt)))) (Array.to_list tord)
	
end 
    

    

      
