(* File generated from shape.idl *)

type internal

(*
 This file is part of the APRON Library, released under LGPL license.
 Please read the COPYING file packaged in the distribution.
*)

 
type t
(** Type of shapes.

Shapes are defined by conjunctions of equalities of the form
x ==l==> y  
and inequalites on length variables of the form
[+/-l_i +/- l_j >= 0].

Abstract values which are shapes have the type [t Apron.AbstractX.t].

Managers allocated for shapes have the type [t Apron.manager.t].
*)

 
external manager_alloc : unit -> t Apron.Manager.t
	= "camlidl_shape_shape_manager_alloc"

