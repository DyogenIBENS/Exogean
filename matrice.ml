(*###########################################################
  
  Ce programme permet de réaliser la matrice des poids pour A, 
  T, C et G à des distances de + ou - 4 bases soit du signal 
  AG (fichier IE.seq) soit du signal GT (fichier EI.seq).
  
  Le signal AG ou GT se trouve à la position 70 de chaque 
  ligne du fichier. Les bases prises en compte sont donc 
  celles se trouvant aux positions :
  66 67 68 69 - 70 71 - 72 73 74 75 .
  
  #############################################################*) 

(*  awk 'NR>=5{print ">\n"$NF}' file.seq > file.fa *)


(* Matrice des poids *)

open Alphaadn
open Seq
open Db
open Cast



let mAGdefault =
[|[|0.0594609613781606; 0.0580716865796054443; 0.207279799944429022;
      0.0383439844401222557; 0.999722145040289; 0.; 0.235065295915532102;
      0.208113364823562103; 0.226729647124201156; 0.22950819672131148|];
    [|0.411503195332036664; 0.472909141428174473; 0.207279799944429022;
      0.193942761878299536; 0.; 0.00027785495971103082; 0.105307029730480686;
      0.373159210891914395; 0.234231731036399; 0.237843845512642399|];
    [|0.0858571825507085246; 0.0819672131147540922; 0.264240066685190345;
      0.00333425951653237027; 0.00027785495971103082; 0.999444290080577891;
      0.490969713809391473; 0.205612670186162833; 0.284523478744095559;
      0.252570158377327048|];
    [|0.443178660739094177; 0.387051958877465963; 0.321200333425951667;
      0.764378994165045844; 0.; 0.00027785495971103082; 0.168657960544595725;
      0.213114754098360643; 0.254515143095304264; 0.280077799388719073|];
    [|0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25|]|];;

let mGTdefault =
[|[|0.267574326201722723; 0.310919699916643533; 0.627118644067796605;
      0.0819672131147540922; 0.000833564879133092567;
      0.000833564879133092567; 0.49736037788274523; 0.661572659071964386;
      0.0586273964990275051; 0.145873853848291185|];
    [|0.170880800222283968; 0.103639899972214511; 0.132536815782161699;
      0.0550152820227841066; 0.00027785495971103082; 0.995276465684912459;
      0.0166712975826618509; 0.0914142817449291462; 0.0600166712975826605;
      0.425951653237010286|];
    [|0.255070853014726318; 0.189497082522923022; 0.130591831064184483;
      0.823284245623784439; 0.998055015282022784; 0.; 0.45179216449013615;
      0.150875243123089753; 0.809947207557654925; 0.230619616560155588|];
    [|0.306474020561267; 0.395943317588218935; 0.109752709085857186;
      0.0397332592386774111; 0.000833564879133092567; 0.00388996943595443191;
      0.0341761600444567964; 0.0961378160600166731; 0.0714087246457349306;
      0.197554876354542941|];
    [|0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25|]|];;




let probATGC s =
  let n = Seq.length s and tprob=Array.create 5 0.0 and ntnonN = ref 0.0 in
    for i=0 to n-1 do
      match s.{i} with
	|Alphaadn.A -> (tprob.(0) <- tprob.(0) +. 1.0); (ntnonN := !ntnonN +. 1.0);
	|Alphaadn.T -> (tprob.(1) <- tprob.(1) +. 1.0); (ntnonN := !ntnonN +. 1.0);
	|Alphaadn.G -> (tprob.(2) <- tprob.(2) +. 1.0); (ntnonN := !ntnonN +. 1.0);
	|Alphaadn.C -> (tprob.(3) <- tprob.(3) +. 1.0); (ntnonN := !ntnonN +. 1.0);
	|Alphaadn.N -> ();   
    done;
    for j=0 to 3 do
      tprob.(j) <- tprob.(j) /. (!ntnonN); 
    done;
    tprob.(4) <- 0.25;
    tprob;;




let cree_matrix fseq posdeb = 
  let matrice = Array.make_matrix 5 10 0. in
  let lseq = List.map (fun s -> Seq.sub s posdeb 10) (List.map (fun (_,_,c) -> c) (readg fasta of_chara nulla (open_in fseq))) in
    
  let nseq = List.length lseq in
    List.iter (fun s -> Array.iteri (fun i b -> match b with |A -> matrice.(0).(i) <- matrice.(0).(i)+.1.; |T -> matrice.(1).(i) <- matrice.(1).(i)+.1.; |G -> matrice.(2).(i) <- matrice.(2).(i)+.1.; |C -> matrice.(3).(i) <- matrice.(3).(i)+.1.; |N -> ();) (Seq.to_array s)) lseq;
    
    (* cas de A, T, G, C *)
    for i = 0 to 3 do
      for j = 0 to 9 do
	matrice.(i).(j) <- matrice.(i).(j) /. float_of_int(nseq);
      done;
    done;

    (* cas du N *)
    for j = 0 to 9 do
      matrice.(4).(j) <- 0.25;
    done;
    matrice




(* Calcul du score d'un signal d'épissage de 10 nucléotides *)
let score matrice pnt seq =
  let score_seq = ref 0. in
    for i = 0 to 9 do
      score_seq:= !score_seq +. log((matrice.(int_of_symb (seq.{i})).(i)) /. (pnt.(int_of_symb (seq.{i}))));
    done;
    !score_seq


