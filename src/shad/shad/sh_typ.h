/* SHAD - Library of shape abstract domains
 * Copyright (C) 2012-2013 LIAFA
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * Please read the COPYING file packaged in the distribution.
 */

/* 
 * sh_typ.h: types, fields and predicates defined in the program and
 *           passed in the global manager
 */

#ifndef _SH_TYP_H_
#define _SH_TYP_H_

#include "sh_loc.h"

#ifdef __cplusplus
extern "C"
  {
#endif

/* ********************************************************************** */
/* I. Types */
/* ********************************************************************** */

/* type for the unique identifiers of types, fields, preds, etc. */
#define size_t guint

/* ====================================================================== */
/* I.1 Kind of types */
/* ====================================================================== */

/* Integer types, taken from CIL */
typedef enum sh_typ_ikind_t
{
  SH_IKIND_BOOL = 0, /* [_Bool] */
  SH_IKIND_CHAR, /* [char] */
  SH_IKIND_SCHAR, /* [signed char] */
  SH_IKIND_UCHAR, /* [unsigned char] */
  SH_IKIND_INT, /* [int] */
  SH_IKIND_UINT, /* [unsigned int] */
  SH_IKIND_SHORT, /* [short] */
  SH_IKIND_USHORT, /* [unsigned short] */
  SH_IKIND_LONG, /* [long] */
  SH_IKIND_ULONG, /* [unsigned long] */
  SH_IKIND_LONGLONG, /* [long long] (or [_int64] on Microsoft Visual C) */
  SH_IKIND_ULONGLONG, /* [unsigned long long] */
  SH_IKIND_UNKNOWN
/* not to be used */
} sh_typ_ikind_t;

/* Float types, taken from CIL */
typedef enum sh_typ_fkind_t
{
  SH_FKIND_FLOAT = 0, /* [float] */
  SH_FKIND_DOUBLE, /* [double] */
  SH_FKIND_LONGDOUBLE, /* [long double] */
  SH_FKIND_UNKNOWN
/* not to be used */
} sh_typ_fkind_t;

/* General types, taken from CIL */
typedef enum sh_typ_discr_t
{
  SH_TYP_TVOID = 0, /* [void] */
  SH_TYP_TINT, /* integers */
  SH_TYP_TFLOAT, /* floating */
  SH_TYP_TPTR, /* pointer */
  SH_TYP_TARRAY, /* array */
  SH_TYP_TFUN, /* function */
  SH_TYP_TNAMED, /* typedef expression */
  SH_TYP_TRECORD, /* record type */
  SH_TYP_TUNION, /* union type */
  SH_TYP_TVALIST, /* variable size list */
  SH_TYP_TFIELD, /* field type, internal use */
  SH_TYP_TOBJECT, /* object type, internal use */
  SH_TYP_UNKNOWN
/* not to be used */
} sh_typ_discr_t;

/* ====================================================================== */
/* I.2 Types definition */
/* ====================================================================== */

/* Definition of pointer types */
typedef struct sh_typ_ptr_t
{
  size_t deref; /* number of * used to dereference the basic type >= 1 */
  size_t base_tid; /* the base type identifier */
} sh_typ_ptr_t;

/* Definition of array types */
typedef struct sh_typ_array_t
{
  size_t elem_tid; /* type of the elements */
  int size; /* number of elements, -1 if not given */
} sh_typ_array_t;

/* Definition of function types:
 * seen as pointer with deref = 1 to the result type */
typedef sh_typ_ptr_t sh_typ_fun_t;

/* Definition of the named types:
 * seens as a pointer with deref = 0 to the named type */
typedef sh_typ_ptr_t sh_typ_named_t;

/* Definition of struct or union types */
typedef struct sh_typ_comp_t
{
  GArray* fields; /* array of field identifiers */
  size_t rec; /* recursive degree, = 0 if not recursive,
   * >= 1 number of recursive fields */
} sh_typ_comp_t;

typedef struct sh_typ_t
{
  char* name; /* type name, as defined/used in the program */
  sh_loc_t* loc; /* location of the (first) definition */
  sh_typ_discr_t discr; /* discriminant for info */
  union
  {
    sh_typ_ikind_t ikind;
    sh_typ_fkind_t fkind;
    sh_typ_ptr_t ptr;
    sh_typ_array_t arr;
    sh_typ_fun_t fun;
    sh_typ_named_t named;
    sh_typ_comp_t comp;
  } info;
} sh_typ_t;

/* ====================================================================== */
/* I.3 Fields definition */
/* ====================================================================== */

typedef struct sh_fld_t
{
  char* name; /* field name */
  sh_loc_t* loc; /* location of its definition */
  size_t o_tid; /* type identifier where the field is defined */
  size_t tid; /* type identified of the type of the field */
  gboolean isrec; /* if recursive field */
} sh_fld_t;

/* ====================================================================== */
/* I.4 Predicate definition */
/* ====================================================================== */

typedef enum sh_pkind_t
{
  SH_PRED_LS = 0,
  SH_PRED_DLL,
  SH_PRED_SLL,
  SH_PRED_NLL,
  SH_PRED_SEP,
  SH_PRED_WSEP,
  SH_PRED_OTHER
} sh_pkind_t;

typedef struct sh_pred_t
{
  char* name; /* predicate name */
  sh_loc_t* loc; /* location of its definition */

  sh_pkind_t kind; /* particular case if identified,
   * kind of the separation operator otherwise */
  GArray* tids; /* types covered by the predicate,
   tids[0] is the specified type */
  size_t minlen; /* minimum length of the segment defined >= 0 */
  size_t arity; /* number of arguments,
   * the first two are recursive */
  size_t ex_vars; /* number of existential vars, the first one is recursive */
  /* The following array contains the definition of the predicate
   * in form of a graph with:
   * - vertices numbered in arity + ex_vars
   * - edges represented by the adjacency matrix def_mat
   *   def_mat[i] is the GPtrArray of sh_prededge_t
   *   of edges starting from i-th vertex
   */
  GPtrArray* def_mat;
} sh_pred_t;

typedef enum
{
  SH_PFORM_PTO, /* points-to */
  SH_PFORM_PRED, /* predicate call */
  SH_PFORM_UNKNOWN
} sh_prededge_e;

/* Edges used in predicate definitions */
typedef struct sh_prededge_s
{
  sh_prededge_e discr; /* pto or predicate */
  size_t src; /* source vertex */
  size_t dst; /* destination vertex */
  union
  {
    size_t fid; /* field for points-to */
    struct
    { /* predicate */
      size_t pid;
      GArray* args; /* Neighbor args, NULL if none */
    } pred;
  } info;
} sh_prededge_t;

/* ====================================================================== */
/* I.5 Environment for types and fields */
/* ====================================================================== */

typedef struct sh_typenv_t
{
  GPtrArray* tinfo; /* array of types; elements of type sh_typ_t* */
  GHashTable* tdict; /* dictionary for type names
   (type name -> int in tinfo) */
  GPtrArray* finfo; /* array of fields; elemets of type sh_fld_t* */
  GHashTable* fdict; /* dictionary for field names
   (type_name.field_name -> int in finfo) */
  GPtrArray* pinfo; /* array of predicates specifying types */
  GHashTable* pdict; /* dictionary for predicate name
   * (pred name -> int in pinfo) */
  int count; /* For memory management by reference counting */
  sh_manager_t* man; /* to retrieve informations about files, etc. */
} sh_typenv_t;

/* ********************************************************************** */
/* II. Functions */
/* ********************************************************************** */

/* ====================================================================== */
/* II.1 Memory management, constructors, destructors */
/* ====================================================================== */

sh_typenv_t*
sh_typenv_alloc_empty(sh_manager_t* man);
/* Build the empty type environment */
void
sh_typenv_free2(sh_typenv_t* tenv);
/* Free the environment
 (the structure itself and the memory pointed to by fields) */

static inline
void
sh_typenv_free(sh_typenv_t* tenv);
/* Decrement the reference counter and possibly deallocate */

void
sh_typ_free(sh_typ_t* data);
void
sh_fld_free(sh_fld_t* data);
void
sh_pred_free(sh_pred_t* data);
void
sh_prededge_free(sh_prededge_t* data);
/* ====================================================================== */
/* II.2 Constructors and destructors for components */
/* ====================================================================== */
/* Types */
size_t
sh_typenv_add_type_void(sh_typenv_t* tenv);
/* Add type void and returns its id even if it is already there */
size_t
sh_typenv_add_type_valist(sh_typenv_t* tenv);
/* Add type valist and returns its id even if it is already there */
size_t
sh_typenv_add_type_int(sh_typenv_t* tenv, sh_typ_ikind_t ikind);
size_t
sh_typenv_add_type_float(sh_typenv_t* tenv, sh_typ_fkind_t fkind);
/* Add a specific numeric type */
size_t
sh_typenv_add_type_ptr(sh_typenv_t* tenv, sh_loc_t* loc, char* tname,
    size_t deref);
/* Add the pointer type with a number of dereferencing to tname,
 * the name of the type appends deref * to tname
 */
size_t
sh_typenv_add_type_ptr0(sh_typenv_t* tenv, sh_loc_t* loc, size_t tid,
    size_t deref);
/* Same but the type identifier of the target type is given */
size_t
sh_typenv_add_type_array(sh_typenv_t* tenv, sh_loc_t* loc, char* tname,
    int size);
/* Add the array type with type of elements tname and size elements,
 * the name is tname[size]
 */
size_t
sh_typenv_add_type_array0(sh_typenv_t* tenv, sh_loc_t* loc, size_t tid,
    int size);
/* Same but the type identifier of the target type is given */
size_t
sh_typenv_add_type_fun(sh_typenv_t* tenv, sh_loc_t* loc, char* tname);
size_t
sh_typenv_add_type_fun0(sh_typenv_t* tenv, sh_loc_t* loc, size_t tid);
/* Add the function type whose result type is given,
 * the name is fun_tname_
 */
size_t
sh_typenv_add_type_named(sh_typenv_t* tenv, sh_loc_t* loc, char* tname,
    char* aname);
size_t
sh_typenv_add_type_named0(sh_typenv_t* tenv, sh_loc_t* loc, size_t tid,
    char* aname);
/* Add an alias type to aname of type tname/tid */
size_t
sh_typenv_add_type_record(sh_typenv_t* tenv, sh_loc_t* loc, char* tname);
/* Add a record type, fields come afterwards */
size_t
sh_typenv_add_type_union(sh_typenv_t* tenv, sh_loc_t* loc, char* tname);
/* Add an union type, fields come afterwards */
/* Fields */
size_t
sh_typenv_add_field(sh_typenv_t* tenv, sh_loc_t* loc, char* fname, char* tname,
    char* ftype);
/* Add a field of the record/union tname, of type ftype */
size_t
sh_typenv_add_field0(sh_typenv_t* tenv, sh_loc_t* loc, char* fname, size_t tid,
    size_t ftid);
/* The same but with identifiers for types */
/* Predicates */
size_t
sh_typenv_add_pred_ls(sh_typenv_t* tenv, sh_loc_t* loc, char* pname,
    size_t nextid);
/* Add a predicate definition pname of singly linked list
 with next field nextid */
size_t
sh_typenv_add_pred_dll(sh_typenv_t* tenv, sh_loc_t* loc, char* pname,
    size_t nextid, size_t previd);
/* Add a predicate definition pname of doubly linked list
 with next field nextid, and prev field previd */
size_t
sh_typenv_add_pred_sll(sh_typenv_t* tenv, sh_loc_t* loc, char* pname,
    size_t nextid, size_t previd);
/* Add a predicate definition pname of singly linked list
 * obtained from dll with all prev field at some node */
size_t
sh_typenv_add_pred_nll(sh_typenv_t* tenv, sh_loc_t* loc, char* pname,
    size_t nextid, size_t fid, size_t npred);
/* Add a predicate definition pname of nested linked list
 with next field nextid, change of level field fid,
 and nested predicate npred */

sh_prededge_t*
sh_prededge_pto(size_t from, size_t fid, size_t to);
/* Add the points-to constraint to the existing one, if possible. */

sh_prededge_t*
sh_prededge_pred(size_t pid, GArray* argv);
/* Build the predicate constraint */

/* ====================================================================== */
/* II.3 Printing */
/* ====================================================================== */

void
sh_typenv_fdump(FILE* stream, sh_typenv_t* tenv);
void
sh_typenv_fdump_types(FILE* stream, sh_typenv_t* tenv);
void
sh_typenv_fdump_fields(FILE* stream, sh_typenv_t* tenv);
void
sh_typenv_fdump_preds(FILE* stream, sh_typenv_t* tenv);
/* Print (in debug mode) a type environment and one or some components */

void
sh_typenv_fdump_tinfo(FILE* stream, sh_typenv_t* tenv, size_t tid);
/* Print type information for type index tid */
void
sh_typenv_fdump_finfo(FILE* stream, sh_typenv_t* tenv, size_t fid);
/* Print field information for field index fid */
void
sh_typenv_fdump_pinfo(FILE* stream, sh_typenv_t* tenv, size_t pid);
/* Print predicate information for predicate index pid */

void
sh_typ_fdump(FILE* stream, sh_typenv_t* tenv, sh_typ_t* t);
/* Print type information for type */
void
sh_fld_fdump(FILE* stream, sh_typenv_t* tenv, sh_fld_t* f);
/* Print field information for field */
void
sh_pred_fdump(FILE* stream, sh_typenv_t* tenv, sh_pred_t* p);
/* Print predicate information for predicate */

/* ====================================================================== */
/* II.4 Tests */
/* ====================================================================== */

/* ====================================================================== */
/* II.5 Getters and setters */
/* ====================================================================== */

/* Types */
sh_typ_t*
sh_typenv_get_type(sh_typenv_t* tenv, size_t tid);
/* Get the type information for the type identifier */

int
sh_typenv_lookup_type(sh_typenv_t* tenv, const char* tname);
/* Get the type identifier >= 0 associated to the name, -1 if absent. */

int
sh_typenv_lookup_type_record(sh_typenv_t* tenv, const char* tname);
/* Get the type identifier >= 0 associated to the name struct_tname
 , -1 if absent. */
size_t
sh_typenv_push_type(sh_typenv_t* tenv, sh_typ_t* t);
/* Add the type and return the identifier */

void
sh_typenv_set_record_recursive(sh_typenv_t* tenv, size_t tid, size_t degree);
/* Set the recursive attribute for record type */

gboolean
sh_typenv_is_type_field(sh_typenv_t* tenv, size_t tid, size_t fid);
/* Check if field fid belongs to type tid */

/* Fields */
sh_fld_t*
sh_typenv_get_field(sh_typenv_t* tenv, size_t fid);
/* Get the field information for the field identifier */

int
sh_typenv_lookup_field(sh_typenv_t* tenv, char* fname, char* tname);
/* Get the field identifier >= 0 associated to the name built
 as tname.fname, -1 if absent. */
size_t
sh_typenv_push_field(sh_typenv_t* tenv, sh_fld_t* f);
/* Add the field and return the identifier */

void
sh_typenv_set_field_recursive(sh_typenv_t* tenv, size_t fid, gboolean isrec);
/* Set the recursive attribute to field or record */

/* Predicates */
sh_pred_t*
sh_typenv_get_pred(sh_typenv_t* tenv, size_t pid);
/* Get the predicate information for the predicate identifier */

int
sh_typenv_lookup_pred(sh_typenv_t* tenv, char* pname);
/* Get the predicate identifier >= 0 associated to the name, -1 if absent. */
size_t
sh_typenv_push_pred(sh_typenv_t* tenv, sh_pred_t* p);
/* Add the predicate and return the identifier */

gboolean
sh_typenv_is_pred_field0(sh_typenv_t* tenv, size_t pid, size_t fid,
    gboolean dir);
/* Return true if field fid belongs to the 0 level of predicate pid
 * with direction forward (dir=TRUE) or backward (dir=FALSE)
 */

/* ====================================================================== */
/* III. Inline definitions */
/* ====================================================================== */

static inline
void
sh_typenv_free(sh_typenv_t* tenv)
{
  if (tenv->count <= 1)
    sh_typenv_free2(tenv);
  else
    tenv->count--;
}
static inline sh_typenv_t*
sh_typenv_copy(sh_typenv_t* tenv)
{
  tenv->count++;
  return tenv;
}

#ifdef __cplusplus
}
#endif

#endif /* _SH_TYP_H_ */
