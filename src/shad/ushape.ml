(* File generated from ushape.idl *)

type internal

(*
 This file is part of the APRON Library, released under LGPL license.
 Please read the COPYING file packaged in the distribution.
*)

 
type t
(** Type of ushapes.

Hgraphs are defined by conjunctions of equalities of the form
x ====> y  
Abstract values which are ushapes have the type [t Apron.AbstractX.t].

Managers allocated for ushapes have the type [t Apron.manager.t].
*)

 
external manager_alloc : unit -> t Apron.Manager.t
	= "camlidl_ushape_ushape_manager_alloc"

