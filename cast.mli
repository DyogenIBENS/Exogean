(* conversion de type pour transformer des entiers en type somme ex : les nucleotides ou les acides amin�s *)

external  int_of_symb  : 'a -> int = "%identity"
external  symb_of_int  : int-> 'a  = "%identity"
