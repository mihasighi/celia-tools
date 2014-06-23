(* File generated from shape.idl *)

type internal

(*
 This file is part of the APRON Library, released under LGPL license.
 Please read the COPYING file packaged in the distribution.
*)

 
(** Shape abstract domain. *)
 


 
type t
(** Type of shapes.

Shapes are defined by conjunctions of equalities of the form
x ==l==> y  
and inequalites on length variables of the form
[+/-l_i +/- l_j >= 0].

Abstract values which are shapes have the type [t Apron.AbstractX.t].

Managers allocated for shapes have the type [t Apron.manager.t].
*)

 
(** Allocate a new manager to manipulate shapes. *)
external manager_alloc : unit -> t Apron.Manager.t
	= "camlidl_shape_shape_manager_alloc"



(**
{2 Compilation information}

{3 Bytecode compilation}
To compile to bytecode, you should first generate a custom interpreter with a
command which should look like:

[ocamlc -I $APRON_PREFIX/lib -make-runtime -o myrun bigarray.cma gmp.cma apron.cma shape.cma]

and then you compile and link your example [X.ml] with

[ocamlc -I $APRON_PREFIX/lib -c X.ml] and

[ocamlc -I $APRON_PREFIX/lib -use-runtime myrun -o X bigarray.cma gmp.cma apron.cma shape.cma X.cmo]

{b Comments:} The C libraries related to [gmp.cma] and [apron.cma] are
automatically looked for (thanks to the auto-linking feature provided by
[ocamlc]).  For [shape.cma], the library [libshape.a] is
selected by default. The [-noautolink] option should be used to select a
different version. See the C documentation of [Shape] library for details.

With the [-noautolink] option, the generation of the custom runtime executable
should be done with

[ocamlc -I $APRON_PREFIX/lib -noautolink -make-runtime -o myrun bigarray.cma gmp.cma apron.cma shape.cma -ccopt "-L$GMP_PREFIX/lib ..." -cclib "-lshape_caml -loct_caml -lapron_caml -lapron -lgmp_caml -lmpfr -lgmp -lbigarray -lcamlidl"]

{3 Native-code compilation}
You compile and link with

[ocamlopt -I $APRON_PREFIX/lib -c X.ml] and

[ocamlopt -I $APRON_PREFIX/lib -o X bigarray.cmxa gmp.cmxa apron.cmxa shape.cmxa X.cmx]

{b Comments:} Same as for bytecode compilation. With the
[-noautolink] option, the linking command becomes

[ocamlopt -I $APRON_PREFIX/lib -o X bigarray.cmxa gmp.cmxa apron.cmxa shape.cmxa -ccopt "-L$GMP_PREFIX/lib ..." -cclib "-lshape_caml -loct_caml -lapron_caml -lapron -lgmp_caml -lmpfr -lgmp -lbigarray -lcamlidl" X.cmx]
*)
