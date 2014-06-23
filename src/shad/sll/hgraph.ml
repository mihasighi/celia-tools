(* File generated from hgraph.idl *)

type internal

(*
 This file is part of the APRON Library, released under LGPL license.
 Please read the COPYING file packaged in the distribution.
*)

 
type t
(** Type of hgraphs.

Hgraphs are defined by conjunctions of equalities of the form
x ====> y  
Abstract values which are hgraphs have the type [t Apron.AbstractX.t].

Managers allocated for hgraphs have the type [t Apron.manager.t].
*)

 
external manager_alloc : unit -> t Apron.Manager.t
	= "camlidl_hgraph_hgraph_manager_alloc"

