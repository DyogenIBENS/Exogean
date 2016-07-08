(* classifgene.ml : fichier qui sert soit à eliminer les artefacts encore presents dans les completegene
   de sortie d'exogean, soit à classer certains complete gene en pseudogènes
   à terme nonmonoclart et nondiclart avec les fonctions prune_èartefact de traiteCollection.ml seront là *)

open Alphaprot
open Donnees_base
open Donnees




(* szcdsnomet est en gal à 240, szcdswithmet à 150 et propcdsmonoex à 20/100*)
let nonpseudogene szcdsnomet szcdswithmet propcdsmonoex sizecdsmonoex sizecdsdiex cg = 
  let g= CompleteGene.gn cg in
    ((CompleteGene.met cg) || (CompleteGene.lgexcds cg)>=szcdsnomet) 
    && 
    ((not (CompleteGene.met cg)) || (CompleteGene.lgexcds cg)>=szcdswithmet)
    &&
    (((CompleteGene.nbexincds cg)>1) || ((CompleteGene.lgexcds cg)>=((propcdsmonoex*(CompleteGene.lgextot cg))/100)))
    &&
    (((CompleteGene.nbexincds cg)!=1) || ((CompleteGene.lgexcds cg)>=sizecdsmonoex))
      
    (* en fait plus complexe, peut-etre taille intron?*)
    && 
    (((CompleteGene.nbexincds cg)!=2) || ((CompleteGene.lgexcds cg)>=sizecdsdiex))
    &&
    (* ajout au 17/03/06 : si que prot alors mod3 et pas de stop *)
    (((Gene.arnpal g)!=(Donnees_base.N)) || ((((List.length (CompleteGene.trans cg)) mod 3)=0) && (not (List.mem Alphaprot.Z (CompleteGene.prot cg)))))
 

