(* File generated from hgraph.idl *)

type internal

(*
 This file is part of the APRON Library, released under LGPL license.
 Please read the COPYING file packaged in the distribution.
*)

 
(** Shape abstract domain. *)
 


 
type t
(** Type of hgraphs.

Hgraphs are defined by conjunctions of equalities of the form
x ====> y  
Abstract values which are hgraphs have the type [t Apron.AbstractX.t].

Managers allocated for hgraphs have the type [t Apron.manager.t].
*)

 
(** Allocate a new manager to manipulate hgraphs. *)
external manager_alloc : unit -> t Apron.Manager.t
	= "camlidl_hgraph_hgraph_manager_alloc"



(**
{2 Compilation information}

{3 Bytecode compilation}
To compile to bytecode, you should first generate a custom interpreter with a
command which should look like:

[ocamlc -I $APRON_PREFIX/lib -make-runtime -o myrun bigarray.cma gmp.cma apron.cma hgraph.cma]

and then you compile and link your example [X.ml] with

[ocamlc -I $APRON_PREFIX/lib -c X.ml] and

[ocamlc -I $APRON_PREFIX/lib -use-runtime myrun -o X bigarray.cma gmp.cma apron.cma hgraph.cma X.cmo]

{b Comments:} The C libraries related to [gmp.cma] and [apron.cma] are
automatically looked for (thanks to the auto-linking feature provided by
[ocamlc]).  For [hgraph.cma], the library [libhgraph.a] is
selected by default. The [-noautolink] option should be used to select a
different version. See the C documentation of [Shape] library for details.

With the [-noautolink] option, the generation of the custom runtime executable
should be done with

[ocamlc -I $APRON_PREFIX/lib -noautolink -make-runtime -o myrun bigarray.cma gmp.cma apron.cma hgraph.cma -ccopt "-L$GMP_PREFIX/lib ..." -cclib "-lhgraph_caml -loct_caml -lapron_caml -lapron -lgmp_caml -lmpfr -lgmp -lbigarray -lcamlidl"]

{3 Native-code compilation}
You compile and link with

[ocamlopt -I $APRON_PREFIX/lib -c X.ml] and

[ocamlopt -I $APRON_PREFIX/lib -o X bigarray.cmxa gmp.cmxa apron.cmxa hgraph.cmxa X.cmx]

{b Comments:} Same as for bytecode compilation. With the
[-noautolink] option, the linking command becomes

[ocamlopt -I $APRON_PREFIX/lib -o X bigarray.cmxa gmp.cmxa apron.cmxa hgraph.cmxa -ccopt "-L$GMP_PREFIX/lib ..." -cclib "-lhgraph_caml -loct_caml -lapron_caml -lapron -lgmp_caml -lmpfr -lgmp -lbigarray -lcamlidl" X.cmx]
*)
