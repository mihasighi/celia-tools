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

#include <assert.h>
#include <string.h>
#include "sh_typ.h"
#include "sh_manager.h"

/* ********************************************************************** */
/* I. Types */
/* ********************************************************************** */

/* ====================================================================== */
/* I.1 Name of predefined types */
/* ====================================================================== */

static const char* sh_typenv_name_int[] =
  { "_Bool", "char", "signed char", "unsigned char", "int", "unsigned int",
      "short", "unsigned short", "long", "unsigned long", "long long",
      "unsigned long long" };

static const char* sh_typenv_name_float[] =
  { "float", "double", "long double" };

static const char* sh_typenv_name_type[] =
  { "void", "int", "float", "*", "array of", "fun", "is a", "struct", "union",
      "valist", "field", "object" };

/* ********************************************************************** */
/* II. Functions */
/* ********************************************************************** */

/* ====================================================================== */
/* II.1 Memory management, constructors, destructors */
/* ====================================================================== */

sh_typenv_t*
sh_typenv_alloc_empty(sh_manager_t* man)
{
  sh_typenv_t* ret = g_new0(sh_typenv_t, 1);
  assert(NULL != ret);

  ret->tinfo = g_ptr_array_new_with_free_func((GDestroyNotify) sh_typ_free);
  assert(NULL != ret->tinfo);

  ret->tdict = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
  assert(NULL != ret->tdict);

  ret->finfo = g_ptr_array_new_with_free_func((GDestroyNotify) sh_fld_free);
  assert(NULL != ret->finfo);

  ret->fdict = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
  assert(NULL != ret->fdict);

  ret->pinfo = g_ptr_array_new_with_free_func((GDestroyNotify) sh_pred_free);
  assert(NULL != ret->pinfo);

  ret->pdict = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
  assert(NULL != ret->pdict);

  ret->count = 1;

  ret->man = man;
  return ret;

}

void
sh_typenv_free2(sh_typenv_t* tenv)
{
  assert(NULL != tenv);

  g_ptr_array_free(tenv->tinfo, TRUE); /* free the content using sh_typ_free */
  tenv->tinfo = NULL;
  g_hash_table_destroy(tenv->tdict); /* the names will be freed here */
  tenv->tdict = NULL;
  g_ptr_array_free(tenv->finfo, TRUE);
  tenv->finfo = NULL;
  g_hash_table_destroy(tenv->fdict);
  tenv->fdict = NULL;
  g_ptr_array_free(tenv->pinfo, TRUE);
  tenv->pinfo = NULL;
  g_hash_table_destroy(tenv->pdict);
  tenv->pdict = NULL;
  tenv->man = NULL;
  free(tenv);
}

void
sh_typ_free(sh_typ_t* data)
{
  assert (NULL != data);

  assert (NULL != data->name);
  free(data->name); // WARNING: the name is strdup in the dictionary
  data->name = NULL;

  if (data->discr == SH_TYP_TRECORD || data->discr == SH_TYP_TUNION)
    {
      assert (NULL != data->info.comp.fields);
      g_array_free(data->info.comp.fields, TRUE);
      data->info.comp.fields = NULL;
    }

  free(data);
}

void
sh_fld_free(sh_fld_t* data)
{
  assert (NULL != data);

  assert (NULL != data->name);
  free(data->name); // WARNING: the name is strdup in the dictionary
  data->name = NULL;

  free(data);
}

void
sh_pred_free(sh_pred_t* data)
{
  assert (NULL != data);

  assert (NULL != data->name);
  free(data->name); // WARNING: the name is strdup in the dictionary
  data->name = NULL;

  assert (NULL != data->tids);
  g_array_free(data->tids, TRUE);
  data->tids = NULL;

  assert (NULL != data->def_mat);
  g_ptr_array_free(data->def_mat, TRUE);
  data->def_mat = NULL;

  free(data);
}

void
sh_prededge_free(sh_prededge_t* data)
{
  assert (NULL != data);

  if (data->discr == SH_PFORM_PRED)
    {
      assert (NULL != data->info.pred.args);
      g_array_free(data->info.pred.args, TRUE);
      data->info.pred.args = NULL;
    }

  free(data);
}

/* ====================================================================== */
/* II.2 Constructors and destructors for components */
/* ====================================================================== */

size_t
sh_typenv_add_type_void(sh_typenv_t* tenv)
{
  assert (NULL != tenv);

  int tid;

  tid = sh_typenv_lookup_type(tenv, sh_typenv_name_type[SH_TYP_TVOID]);
  if (tid < 0)
    {
      sh_typ_t* tvoid = g_new(sh_typ_t, 1);
      assert (NULL != tvoid);

      tvoid->name = g_strdup(sh_typenv_name_type[SH_TYP_TVOID]);
      tvoid->loc = NULL;
      tvoid->discr = SH_TYP_TVOID;

      tid = sh_typenv_push_type(tenv, tvoid);
    }

  assert (NULL != g_ptr_array_index(tenv->tinfo,tid));

  return (size_t) tid;
}

size_t
sh_typenv_add_type_valist(sh_typenv_t* tenv)
{
  assert (NULL != tenv);

  int tid;

  tid = sh_typenv_lookup_type(tenv, sh_typenv_name_type[SH_TYP_TVALIST]);
  if (tid < 0)
    {
      sh_typ_t* tva = g_new(sh_typ_t, 1);
      assert (NULL != tva);

      tva->name = g_strdup(sh_typenv_name_type[SH_TYP_TVALIST]);
      tva->loc = NULL;
      tva->discr = SH_TYP_TVALIST;

      tid = sh_typenv_push_type(tenv, tva);
    }

  assert (NULL != g_ptr_array_index(tenv->tinfo,tid));

  return (size_t) tid;
}

size_t
sh_typenv_add_type_int(sh_typenv_t* tenv, sh_typ_ikind_t ikind)
{
  assert (NULL != tenv);

  int new_tid = sh_typenv_lookup_type(tenv, sh_typenv_name_int[ikind]);
  if (new_tid < 0)
    {
      sh_typ_t* tint = g_new(sh_typ_t, 1);
      assert (NULL != tint);

      tint->name = g_strdup(sh_typenv_name_int[ikind]);
      tint->loc = NULL;
      tint->discr = SH_TYP_TINT;
      tint->info.ikind = ikind;

      new_tid = sh_typenv_push_type(tenv, tint);
    }

  assert (NULL != g_ptr_array_index(tenv->tinfo,new_tid));

  return (size_t) new_tid;
}

size_t
sh_typenv_add_type_float(sh_typenv_t* tenv, sh_typ_fkind_t fkind)
{
  assert (NULL != tenv);

  int new_tid = sh_typenv_lookup_type(tenv, sh_typenv_name_float[fkind]);
  if (new_tid < 0)
    {
      sh_typ_t* tf = g_new(sh_typ_t, 1);
      assert (NULL != tf);

      tf->name = g_strdup(sh_typenv_name_int[fkind]);
      tf->loc = NULL;
      tf->discr = SH_TYP_TFLOAT;
      tf->info.fkind = fkind;

      new_tid = sh_typenv_push_type(tenv, tf);
    }

  assert (NULL != g_ptr_array_index(tenv->tinfo,new_tid));

  return (size_t) new_tid;
}

size_t
sh_typenv_add_type_ptr(sh_typenv_t* tenv, sh_loc_t* loc, char* tname,
    size_t deref)
{
  assert (NULL != tenv);
  assert (NULL != tname);
  assert (deref >= 1);

  /* build the name tname*...* (deref *) */
  int ptr_name_sz = strlen(tname) + deref + 1;
  char * ptr_name = g_new(char, ptr_name_sz);
  assert (NULL != ptr_name);
  snprintf (ptr_name, ptr_name_sz, "%s", tname);
  for (int i = deref; i >= 1; i--)
    ptr_name[ptr_name_sz - 1 - i] = '*';
  ptr_name[ptr_name_sz - 1] = '\0';

  int new_tid = sh_typenv_lookup_type(tenv, ptr_name);
  if (new_tid < 0)
    {
      sh_typ_t* t = g_new (sh_typ_t, 1);
      assert (NULL != t);

      t->name = ptr_name;
      t->loc = loc;
      t->discr = SH_TYP_TPTR;
      t->info.ptr.deref = deref;
      t->info.ptr.base_tid = sh_typenv_lookup_type(tenv, tname);

      new_tid = sh_typenv_push_type(tenv, t);
    }

  assert (NULL != g_ptr_array_index(tenv->tinfo,new_tid));

  return (size_t) new_tid;
}

size_t
sh_typenv_add_type_ptr0(sh_typenv_t* tenv, sh_loc_t* loc, size_t tid,
    size_t deref)
{
  assert (NULL != tenv);
  assert (deref >= 1);

  /* get the name of the tid */
  sh_typ_t* tbase = sh_typenv_get_type(tenv, tid);
  assert (NULL != tbase);
  char* tname = tbase->name;

  return sh_typenv_add_type_ptr(tenv, loc, tname, deref);
}

size_t
sh_typenv_add_type_array(sh_typenv_t* tenv, sh_loc_t* loc, char* tname,
    int size)
{
  assert (NULL != tenv);
  assert (NULL != tname);
  assert (size >= -1);

  /* build the name tname[] or tname[size] */
  int ptr_name_sz = strlen(tname) + ((size < 0) ? 0 : 4) + 3;
  char * ptr_name = g_new(char, ptr_name_sz);
  assert (NULL != ptr_name);
  if (size < 0)
    snprintf (ptr_name, ptr_name_sz, "%s[]", tname);
  else
    snprintf (ptr_name, ptr_name_sz, "%s[%4d]", tname, size);

  int new_tid = sh_typenv_lookup_type(tenv, ptr_name);
  if (new_tid < 0)
    {
      sh_typ_t* t = g_new (sh_typ_t, 1);
      assert (NULL != t);

      t->name = ptr_name;
      t->loc = loc;
      t->discr = SH_TYP_TARRAY;
      t->info.arr.elem_tid = sh_typenv_lookup_type(tenv, tname);
      t->info.arr.size = size;

      new_tid = sh_typenv_push_type(tenv, t);
    }

  assert (NULL != g_ptr_array_index(tenv->tinfo,new_tid));

  return (size_t) new_tid;
}

size_t
sh_typenv_add_type_array0(sh_typenv_t* tenv, sh_loc_t* loc, size_t tid,
    int size)
{
  assert (NULL != tenv);
  assert (size >= -1);

  /* get the name of the tid */
  sh_typ_t* tbase = sh_typenv_get_type(tenv, tid);
  assert (NULL != tbase);
  char* tname = tbase->name;

  return sh_typenv_add_type_array(tenv, loc, tname, size);
}

size_t
sh_typenv_add_type_fun(sh_typenv_t* tenv, sh_loc_t* loc, char* tname)
{
  assert (NULL != tenv);
  assert (NULL != tname);

  /* build the name fun_tname* */
  int ptr_name_sz = strlen(tname) + 5 + 1;
  char * ptr_name = g_new(char, ptr_name_sz);
  assert (NULL != ptr_name);
  snprintf (ptr_name, ptr_name_sz, "fun_%s*", tname);

  int new_tid = sh_typenv_lookup_type(tenv, ptr_name);
  if (new_tid < 0)
    {
      sh_typ_t* t = g_new (sh_typ_t, 1);
      assert (NULL != t);

      t->name = ptr_name;
      t->loc = loc;
      t->discr = SH_TYP_TFUN;
      t->info.ptr.deref = 1;
      t->info.ptr.base_tid = sh_typenv_lookup_type(tenv, tname);

      new_tid = sh_typenv_push_type(tenv, t);
    }

  assert (NULL != g_ptr_array_index(tenv->tinfo,new_tid));

  return (size_t) new_tid;
}

size_t
sh_typenv_add_type_fun0(sh_typenv_t* tenv, sh_loc_t* loc, size_t tid)
{
  assert (NULL != tenv);

  /* get the name of the tid */
  sh_typ_t* tbase = sh_typenv_get_type(tenv, tid);
  assert (NULL != tbase);
  char* tname = tbase->name;

  return sh_typenv_add_type_fun(tenv, loc, tname);
}

size_t
sh_typenv_add_type_named(sh_typenv_t* tenv, sh_loc_t* loc, char* tname,
    char* aname)
{
  assert (NULL != tenv);
  assert (NULL != tname);
  assert (NULL != aname);

  int new_tid = sh_typenv_lookup_type(tenv, aname);
  if (new_tid < 0)
    {
      sh_typ_t* t = g_new(sh_typ_t, 1);
      assert (NULL != t);

      t->name = aname;
      t->loc = loc;
      t->discr = SH_TYP_TNAMED;
      t->info.ptr.deref = 0;
      t->info.ptr.base_tid = sh_typenv_lookup_type(tenv, tname);

      new_tid = sh_typenv_push_type(tenv, t);
    }

  assert (NULL != g_ptr_array_index(tenv->tinfo,new_tid));

  return (size_t) new_tid;
}

size_t
sh_typenv_add_type_named0(sh_typenv_t* tenv, sh_loc_t* loc, size_t tid,
    char* aname)
{
  assert (NULL != tenv);
  assert (NULL != aname);

  /* get the name of the tid */
  sh_typ_t* tbase = sh_typenv_get_type(tenv, tid);
  assert (NULL != tbase);
  char* tname = tbase->name;

  return sh_typenv_add_type_named(tenv, loc, tname, aname);
}

size_t
sh_typenv_add_type_record(sh_typenv_t* tenv, sh_loc_t* loc, char* tname)
{
  assert (NULL != tenv);
  assert (NULL != tname);

  /* build the name struct_tname */
  int str_name_sz = strlen(tname) + 7 + 1;
  char * str_name = g_new(char, str_name_sz);
  assert (NULL != str_name);
  snprintf (str_name, str_name_sz, "struct_%s*", tname);

  int new_tid = sh_typenv_lookup_type(tenv, str_name);
  if (new_tid < 0)
    {
      sh_typ_t* t = g_new (sh_typ_t, 1);
      assert (NULL != t);

      t->name = str_name;
      t->loc = loc;
      t->discr = SH_TYP_TRECORD;
      // the array of fields is updated by sh_typenv_add_field*
      t->info.comp.fields = g_array_new(FALSE, FALSE, sizeof(size_t));
      t->info.comp.rec = 0; // to be updated by sh_typeenv_set_*_recursive

      new_tid = sh_typenv_push_type(tenv, t);
    }

  assert (NULL != g_ptr_array_index(tenv->tinfo,new_tid));

  return (size_t) new_tid;
}

size_t
sh_typenv_add_type_union(sh_typenv_t* tenv, sh_loc_t* loc, char* tname)
{
  assert (NULL != tenv);
  assert (NULL != tname);

  /* build the name union_tname */
  int union_name_sz = strlen(tname) + 7 + 1;
  char * union_name = g_new(char, union_name_sz);
  assert (NULL != union_name);
  snprintf (union_name, union_name_sz, "union_%s*", tname);

  int new_tid = sh_typenv_lookup_type(tenv, union_name);
  if (new_tid < 0)
    {
      sh_typ_t* t = g_new (sh_typ_t, 1);
      assert (NULL != t);

      t->name = union_name;
      t->loc = loc;
      t->discr = SH_TYP_TUNION;
      // the array of fields is updated by sh_typenv_add_field*
      t->info.comp.fields = g_array_new(FALSE, FALSE, sizeof(size_t));
      t->info.comp.rec = 0; // never recursive ?!

      new_tid = sh_typenv_push_type(tenv, t);
    }

  assert (NULL != g_ptr_array_index(tenv->tinfo,new_tid));

  return (size_t) new_tid;
}

/* Fields */
size_t
sh_typenv_add_field(sh_typenv_t* tenv, sh_loc_t* loc, char* fname, char* tname,
    char* ftype)
{
  assert (NULL != tenv);
  assert (NULL != loc);
  assert (NULL != fname);
  assert (NULL != tname); // name of the recode/union containing the field
  assert (NULL != ftype); // type given to this field in the record/union

  // test if the field already exists
  int id = sh_typenv_lookup_field(tenv, fname, tname);
  if (id >= 0)
    return (size_t) id;

  // type tname is a struct normaly
  int o_tid = sh_typenv_lookup_type_record(tenv, tname);
  assert (-1 != o_tid); // the owner shall be defined

  int f_tid = sh_typenv_lookup_type(tenv, ftype);
  assert (-1 != f_tid); // the type shall be defined

  sh_fld_t* fld = g_new(sh_fld_t, 1);
  assert (NULL != fld);
  fld->name = g_strdup(fname);
  fld->loc = loc;
  fld->o_tid = o_tid;
  fld->tid = f_tid;

  return (size_t) sh_typenv_push_field(tenv, fld);
}

size_t
sh_typenv_add_field0(sh_typenv_t* tenv, sh_loc_t* loc, char* fname, size_t tid,
    size_t ftid)
{
  assert (NULL != tenv);
  assert (NULL != fname);

  sh_fld_t* fld = g_new(sh_fld_t, 1);
  assert (NULL != fld);
  fld->name = g_strdup(fname);
  fld->loc = loc;
  fld->o_tid = tid;
  fld->tid = ftid;

  return (size_t) sh_typenv_push_field(tenv, fld);
}

/* Predicates */
size_t
sh_typenv_add_pred_ls(sh_typenv_t* tenv, sh_loc_t* loc, char* pname,
    size_t nextid)
{
  assert (NULL != tenv);
  assert (NULL != loc);
  assert (NULL != pname);

  sh_fld_t* fld_next = sh_typenv_get_field(tenv, nextid);
  assert (NULL != fld_next);
  size_t tid = fld_next->tid; // the type of the predicate

  sh_pred_t* p = g_new(sh_pred_t, 1);
  assert (NULL != p);
  p->name = g_strdup(pname);
  p->loc = loc;
  p->kind = SH_PRED_LS;
  p->tids = g_array_new(FALSE, FALSE, sizeof(size_t));
  g_array_append_val(p->tids, tid);
  p->minlen = 0;
  p->arity = 2;
  p->ex_vars = 1;
  p->def_mat
      = g_ptr_array_new_with_free_func((GDestroyNotify) g_ptr_array_free);
  /* edges from in: one */
  GPtrArray* e_in = g_ptr_array_new_with_free_func(
      (GDestroyNotify) sh_prededge_free);
  g_ptr_array_add(e_in, sh_prededge_pto(0, nextid, 2)); /* in --> (nextid, u) */
  g_ptr_array_add(p->def_mat, e_in);
  /* edges from out: none */;
  g_ptr_array_add(p->def_mat, NULL);
  /* edges from u: one */
  GPtrArray* e_u = g_ptr_array_new_with_free_func(
      (GDestroyNotify) sh_prededge_free);
  GArray* argv = g_array_new(FALSE, FALSE, sizeof(size_t));
  size_t arg = 2; /* u */
  g_array_append_val(argv, arg); /* need a left value for arg */
  arg = 1; /* out */
  g_array_append_val(argv, arg);
  int pos = tenv->pinfo->len;
  g_ptr_array_add(e_u, sh_prededge_pred(pos, argv)); /* p(u,out) */
  g_ptr_array_add(p->def_mat, e_u);

  sh_typenv_push_pred(tenv, p);

  return (size_t) pos;
}

size_t
sh_typenv_add_pred_dll(sh_typenv_t* tenv, sh_loc_t* loc, char* pname,
    size_t nextid, size_t previd)
{
  assert (NULL != tenv);
  assert (NULL != loc);
  assert (NULL != pname);

  sh_fld_t* fld_next = sh_typenv_get_field(tenv, nextid);
  size_t tid = fld_next->o_tid; // the type of the predicate

  sh_pred_t* p = g_new(sh_pred_t, 1);
  assert (NULL != p);
  p->name = g_strdup(pname);
  p->loc = loc;
  p->kind = SH_PRED_DLL;
  p->tids = g_array_new(FALSE, FALSE, sizeof(size_t));
  g_array_append_val(p->tids, tid);
  p->minlen = 0;
  p->arity = 2;
  p->ex_vars = 1;
  p->def_mat
      = g_ptr_array_new_with_free_func((GDestroyNotify) g_ptr_array_free);
  /* edges from in: one */
  GPtrArray* e_in = g_ptr_array_new_with_free_func(
      (GDestroyNotify) sh_prededge_free);
  g_ptr_array_add(e_in, sh_prededge_pto(0, nextid, 2)); /* in --> (nextid, u) */
  g_ptr_array_add(p->def_mat, e_in);
  /* edges from out: none */;
  g_ptr_array_add(p->def_mat, NULL);
  /* edges from u: two */
  GPtrArray* e_u = g_ptr_array_new_with_free_func(
      (GDestroyNotify) sh_prededge_free);
  g_ptr_array_add(e_in, sh_prededge_pto(2, previd, 0)); /* u --> (previd, in) */
  GArray* argv = g_array_new(FALSE, FALSE, sizeof(size_t));
  size_t arg = 2; /* u */
  g_array_append_val(argv, arg); /* need a left value for arg */
  arg = 1; /* out */
  g_array_append_val(argv, arg);
  int pos = tenv->pinfo->len;
  g_ptr_array_add(e_u, sh_prededge_pred(pos, argv)); /* p(u,out) */
  g_ptr_array_add(p->def_mat, e_u);

  sh_typenv_push_pred(tenv, p);

  return (size_t) pos;
}

size_t
sh_typenv_add_pred_sll(sh_typenv_t* tenv, sh_loc_t* loc, char* pname,
    size_t nextid, size_t previd)
{
  assert (NULL != tenv);
  assert (NULL != loc);
  assert (NULL != pname);

  sh_fld_t* fld_next = sh_typenv_get_field(tenv, nextid);
  size_t tid = fld_next->o_tid; // the type of the predicate

  sh_pred_t* p = g_new(sh_pred_t, 1);
  assert (NULL != p);
  p->name = g_strdup(pname);
  p->loc = loc;
  p->kind = SH_PRED_DLL;
  p->tids = g_array_new(FALSE, FALSE, sizeof(size_t));
  g_array_append_val(p->tids, tid);
  p->minlen = 0;
  p->arity = 3;
  p->ex_vars = 1;
  p->def_mat
      = g_ptr_array_new_with_free_func((GDestroyNotify) g_ptr_array_free);
  /* edges from in: two */
  GPtrArray* e_in = g_ptr_array_new_with_free_func(
      (GDestroyNotify) sh_prededge_free);
  g_ptr_array_add(e_in, sh_prededge_pto(0, nextid, 3)); /* in --> (nextid, u) */
  g_ptr_array_add(e_in, sh_prededge_pto(0, previd, 2)); /* in --> (previd, n) */
  g_ptr_array_add(p->def_mat, e_in);
  /* edges from out: none */;
  g_ptr_array_add(p->def_mat, NULL);
  /* edges from u: one */
  GPtrArray* e_u = g_ptr_array_new_with_free_func(
      (GDestroyNotify) sh_prededge_free);
  GArray* argv = g_array_new(FALSE, FALSE, sizeof(size_t));
  size_t arg = 2; /* u */
  g_array_append_val(argv, arg); /* need a left value for arg */
  arg = 1; /* out */
  g_array_append_val(argv, arg);
  arg = 2; /* neighbor */
  g_array_append_val(argv, arg);
  int pos = tenv->pinfo->len;
  g_ptr_array_add(e_u, sh_prededge_pred(pos, argv)); /* p(u,out,n) */
  g_ptr_array_add(p->def_mat, e_u);

  sh_typenv_push_pred(tenv, p);

  return (size_t) pos;
}

size_t
sh_typenv_add_pred_nll(sh_typenv_t* tenv, sh_loc_t* loc, char* pname,
    size_t nextid, size_t fid, size_t nestedpred)
{
  assert (NULL != tenv);
  assert (NULL != loc);
  assert (NULL != pname);

  sh_fld_t* fld_next = sh_typenv_get_field(tenv, nextid);
  size_t tid = fld_next->o_tid; // the type of the predicate

  sh_fld_t* fld_inside = sh_typenv_get_field(tenv, fid);
  size_t tid_inside = fld_inside->tid; // the type of the nested list

  sh_pred_t* p = g_new(sh_pred_t, 1);
  assert (NULL != p);
  p->name = g_strdup(pname);
  p->loc = loc;
  p->kind = SH_PRED_NLL;
  p->tids = g_array_new(FALSE, FALSE, sizeof(size_t));
  g_array_append_val(p->tids, tid);
  g_array_append_val(p->tids, tid_inside);
  p->minlen = 0;
  p->arity = 3;
  p->ex_vars = 2; /* u = 3, v = 4 */

  p->def_mat
      = g_ptr_array_new_with_free_func((GDestroyNotify) g_ptr_array_free);
  /* edges from in: two */
  GPtrArray* e_in = g_ptr_array_new_with_free_func(
      (GDestroyNotify) sh_prededge_free);
  g_ptr_array_add(e_in, sh_prededge_pto(0, nextid, 3)); /* in --> (nextid, u) */
  g_ptr_array_add(e_in, sh_prededge_pto(0, fid, 4)); /* in --> (fid, v) */
  g_ptr_array_add(p->def_mat, e_in);
  /* edges from out: none */;
  g_ptr_array_add(p->def_mat, NULL);
  /* edges from neighbor: none */;
  g_ptr_array_add(p->def_mat, NULL);
  /* edges from u: one pname(u,out,n) */
  GPtrArray* e_u = g_ptr_array_new_with_free_func(
      (GDestroyNotify) sh_prededge_free);
  GArray* argv = g_array_new(FALSE, FALSE, sizeof(size_t));
  size_t arg = 3; /* u */
  g_array_append_val(argv, arg); /* need a left value for arg */
  arg = 1; /* out */
  g_array_append_val(argv, arg);
  arg = 2; /* neighbor */
  g_array_append_val(argv, arg);
  int pos = tenv->pinfo->len;
  g_ptr_array_add(e_u, sh_prededge_pred(pos, argv)); /* p(u,out,n) */
  g_ptr_array_add(p->def_mat, e_u);
  /* edges from v: one nestedpred(v,n) */
  GPtrArray* e_v = g_ptr_array_new_with_free_func(
      (GDestroyNotify) sh_prededge_free);
  argv = g_array_new(FALSE, FALSE, sizeof(size_t));
  arg = 4; /* v */
  g_array_append_val(argv, arg);
  arg = 2; /* v */
  g_array_append_val(argv, arg);
  g_ptr_array_add(e_u, sh_prededge_pred(nestedpred, argv)); /* p(u,out,n) */
  g_ptr_array_add(p->def_mat, e_v);

  sh_typenv_push_pred(tenv, p);

  return (size_t) pos;
}

sh_prededge_t*
sh_prededge_pto(size_t from, size_t fid, size_t to)
{
  sh_prededge_t* r = g_new(sh_prededge_t, 1);
  r->discr = SH_PFORM_PTO;
  r->src = from;
  r->info.fid = fid;
  r->dst = to;
  return r;
}

sh_prededge_t*
sh_prededge_pred(size_t pid, GArray* argv)
{
  sh_prededge_t* r = g_new(sh_prededge_t, 1);
  r->discr = SH_PFORM_PRED;
  r->info.pred.pid = pid;
  r->info.pred.args = argv;
  r->src = g_array_index(argv,size_t,0);
  r->dst = g_array_index(argv,size_t,1);

  return r;
}

/* ====================================================================== */
/* II.3 Printing */
/* ====================================================================== */

void
sh_typenv_fdump(FILE* stream, sh_typenv_t* tenv)
{
  assert (NULL != stream);
  assert (NULL != tenv);

  fprintf(stream, "Type env = [[\n");

  sh_typenv_fdump_types(stream, tenv);
  sh_typenv_fdump_fields(stream, tenv);
  sh_typenv_fdump_preds(stream, tenv);

  fprintf(stream, "]]\n");
}

void
sh_typenv_fdump_types(FILE* stream, sh_typenv_t* tenv)
{
  assert (NULL != stream);
  assert (NULL != tenv);

  fprintf(stream, "Types => {\n");
  if (tenv->tinfo == NULL)
    fprintf(stream, "(null)");
  for (guint i = 0; i < tenv->tinfo->len; i++)
    {
      fprintf(stream, "[typ-%d] ", i);
      sh_typ_fdump(stream, tenv, g_ptr_array_index(tenv->tinfo, i));
    }
  fprintf(stream, "},\n");
}

void
sh_typenv_fdump_fields(FILE* stream, sh_typenv_t* tenv)
{
  assert (NULL != stream);
  assert (NULL != tenv);

  fprintf(stream, "Fields => {\n");
  if (tenv->finfo == NULL)
    fprintf(stream, "(null)");
  for (guint i = 0; i < tenv->finfo->len; i++)
    {
      fprintf(stream, "[fld-%d] ", i);
      sh_fld_fdump(stream, tenv, g_ptr_array_index(tenv->finfo, i));
    }
  fprintf(stream, "},\n");

}

void
sh_typenv_fdump_preds(FILE* stream, sh_typenv_t* tenv)
{
  assert (NULL != stream);
  assert (NULL != tenv);

  fprintf(stream, "Preds => {\n");
  if (tenv->pinfo == NULL)
    fprintf(stream, "(null)");
  for (guint i = 0; i < tenv->pinfo->len; i++)
    sh_pred_fdump(stream, tenv, g_ptr_array_index(tenv->pinfo, i));
  fprintf(stream, "}\n");

}

void
sh_typenv_fdump_tinfo(FILE* stream, sh_typenv_t* tenv, size_t tid)
{
  assert (NULL != tenv);

  if (tenv->tinfo == NULL || tenv->tinfo->len <= tid)
    fprintf(stream, "(bad type id)\n");

  sh_typ_t* t = g_ptr_array_index(tenv->tinfo, tid);
  sh_typ_fdump(stream, tenv, t);

}

void
sh_typenv_fdump_finfo(FILE* stream, sh_typenv_t* tenv, size_t fid)
{
  assert (NULL != tenv);

  if (tenv->finfo == NULL || tenv->finfo->len <= fid)
    fprintf(stream, "(bad field id)\n");

  sh_fld_t* f = g_ptr_array_index(tenv->finfo, fid);
  sh_fld_fdump(stream, tenv, f);

}

void
sh_typenv_fdump_pinfo(FILE* stream, sh_typenv_t* tenv, size_t pid)
{
  assert (NULL != tenv);

  if (tenv->pinfo == NULL || tenv->pinfo->len <= pid)
    fprintf(stream, "(bad pred id)\n");

  sh_pred_t* p = g_ptr_array_index(tenv->pinfo, pid);
  sh_pred_fdump(stream, tenv, p);

}

void
sh_typ_fdump(FILE* stream, sh_typenv_t* tenv, sh_typ_t* t)
{
  fprintf(stream, "%s ", t->name);
  sh_loc_fprint(stream, sh_manager_get_filenv(tenv->man), t->loc);
  fprintf(stream, " =  ");
  switch (t->discr)
    {
  case SH_TYP_TVOID:
    fprintf(stream, " (void) ");
    break;
  case SH_TYP_TINT:
    fprintf(stream, " (int) "); // TODO detailed print
    break;
  case SH_TYP_TFLOAT:
    fprintf(stream, " (float) "); // TODO detailed print
    break;
  case SH_TYP_TPTR:
    fprintf(stream, " (*) "); // TODO detailed print
    break;
  case SH_TYP_TARRAY:
    fprintf(stream, " array of "); // TODO detailed print
    break;
  case SH_TYP_TFUN:
    fprintf(stream, " fun "); // TODO detailed print
    break;
  case SH_TYP_TNAMED:
    fprintf(stream, " alias "); // TODO detailed print
    break;
  case SH_TYP_TRECORD:
    fprintf(stream, " record "); // TODO detailed print
    break;
  case SH_TYP_TUNION:
    fprintf(stream, " union "); // TODO detailed print
    break;
  case SH_TYP_TVALIST:
    fprintf(stream, " valist "); // TODO detailed print
    break;
  case SH_TYP_TFIELD:
  case SH_TYP_TOBJECT:
    fprintf(stream, " (internal) "); // TODO detailed print
    break;
  default:
    fprintf(stream, "(error)");
    }
  fprintf(stream, "\n");
}

void
sh_fld_fdump(FILE* stream, sh_typenv_t* tenv, sh_fld_t* f)
{
  assert (NULL != tenv);

  if (f == NULL)
    {
      fprintf(stream, "(null),\n");
      return;
    }
  fprintf(stream, " %s: typ-%d --> typ-%d ", f->name, f->o_tid, f->tid); // TODO: print type names
  sh_loc_fprint(stream, sh_manager_get_filenv(tenv->man), f->loc);
  fprintf(stream, "\n");

}

void
sh_pred_fdump(FILE* stream, sh_typenv_t* tenv, sh_pred_t* p)
{

  assert (NULL != tenv);

  if (p == NULL)
    {
      fprintf(stream, "(null),\n");
      return;
    }
  fprintf(stream, " %s(x[%d]): typ-%d ", p->name, p->arity, g_array_index(
      p->tids, size_t, 0));
  sh_loc_fprint(stream, sh_manager_get_filenv(tenv->man), p->loc);
  fprintf(stream, "\n");

}

/* ====================================================================== */
/* II.5 Getters and setters */
/* ====================================================================== */

sh_typ_t*
sh_typenv_get_type(sh_typenv_t* tenv, size_t tid)
{
  assert (NULL != tenv);

  if (tenv->tinfo == NULL || tenv->tinfo->len <= tid)
    return (sh_typ_t*) NULL;

  return g_ptr_array_index(tenv->tinfo, tid);
}

int
sh_typenv_lookup_type(sh_typenv_t* tenv, const char* tname)
{
  assert(NULL != tenv);
  assert(NULL != tname);
  gpointer orig, lookup;

  if (g_hash_table_lookup_extended(tenv->tdict, (gconstpointer) tname, &orig,
      &lookup))
    return GPOINTER_TO_INT(lookup);
  else
    return -1;
}

int
sh_typenv_lookup_type_record(sh_typenv_t* tenv, const char* tname)
{
  assert(NULL != tenv);
  assert(NULL != tname);

  /* build the name struct_tname */
  int str_name_sz = strlen(tname) + 7 + 1;
  char * str_name = g_new(char, str_name_sz);
  assert (NULL != str_name);
  snprintf (str_name, str_name_sz, "struct_%s*", tname);

  gpointer orig, lookup;
  int r;

  if (g_hash_table_lookup_extended(tenv->tdict, (gconstpointer) str_name,
      &orig, &lookup))
    r = GPOINTER_TO_INT(lookup);
  else
    r = -1;

  g_free(str_name);
  return r;
}

size_t
sh_typenv_push_type(sh_typenv_t* tenv, sh_typ_t* t)
{
  assert (NULL != tenv);
  assert (NULL != tenv->tinfo);
  assert (NULL != t);

  int pos = tenv->tinfo->len;
  g_ptr_array_add(tenv->tinfo, t);
  g_hash_table_insert(tenv->tdict, t->name, GINT_TO_POINTER(pos));

  return (size_t) pos;
}

void
sh_typenv_set_record_recursive(sh_typenv_t* tenv, size_t tid, size_t degree)
{
  sh_typ_t* t = sh_typenv_get_type(tenv, tid);
  assert (t != NULL);
  assert (t->discr == SH_TYP_TRECORD);

  t->info.comp.rec = degree;

  return;
}

gboolean
sh_typenv_is_type_field(sh_typenv_t* tenv, size_t tid, size_t fid)
{
  sh_fld_t* fld = sh_typenv_get_field(tenv, fid);
  size_t o_tid = fld->o_tid;
  size_t btid = tid;
  sh_typ_t* ty = sh_typenv_get_type(tenv, btid);
  while ((o_tid != tid) && (ty != NULL) && (ty->discr == SH_TYP_TPTR
      || ty->discr == SH_TYP_TNAMED))
    {
      btid = ty->info.ptr.base_tid;
      ty = sh_typenv_get_type(tenv, btid);
    }
  if (btid != o_tid)
    return FALSE;
  return TRUE;
}

sh_fld_t*
sh_typenv_get_field(sh_typenv_t* tenv, size_t fid)
{
  assert (NULL != tenv);

  if (tenv->finfo == NULL || tenv->finfo->len <= fid)
    return (sh_fld_t*) NULL;

  return g_ptr_array_index(tenv->finfo, fid);
}

int
sh_typenv_lookup_field(sh_typenv_t* tenv, char* fname, char* tname)
{
  assert(NULL != tenv);
  assert(NULL != fname);
  assert(NULL != tname);

  /* build the name in the dictionary tname.fname */
  int fullname_sz = strlen(fname) + strlen(tname) + 2; /* +1 ., +1 \0 */
  char* fullname = g_new(char, fullname_sz);
  snprintf(fullname, fullname_sz, "%s.%s", tname, fname);

  gpointer orig, lookup;
  int res;
  if (g_hash_table_lookup_extended(tenv->fdict, (gconstpointer) fullname,
      &orig, &lookup))
    res = GPOINTER_TO_INT(lookup);
  else
    res = -1;

  g_free(fullname);
  return res;
}

size_t
sh_typenv_push_field(sh_typenv_t* tenv, sh_fld_t* f)
{
  assert (NULL != tenv);
  assert (NULL != f);

  /* push record in the array */
  int pos = tenv->finfo->len;
  g_ptr_array_add(tenv->finfo, f);

  sh_typ_t* t = sh_typenv_get_type(tenv, f->o_tid);

  /* build the name in the dictionary tname.fname */
  int fullname_sz = strlen(f->name) + strlen(t->name) + 2; /* +1 ., +1 \0 */
  char* fullname = g_new(char, fullname_sz);
  snprintf(fullname, fullname_sz, "%s.%s", t->name, f->name);
  g_hash_table_insert(tenv->fdict, fullname, GINT_TO_POINTER(pos));

  return (size_t) pos;
}

void
sh_typenv_set_field_recursive(sh_typenv_t* tenv, size_t fid, gboolean isrec)
{
  sh_fld_t* f = sh_typenv_get_field(tenv, fid);
  assert (f != NULL);

  f->isrec = isrec;

  return;
}

sh_pred_t*
sh_typenv_get_pred(sh_typenv_t* tenv, size_t pid)
{
  assert (NULL != tenv);

  if (tenv->pinfo == NULL || tenv->pinfo->len <= pid)
    return (sh_pred_t*) NULL;

  return g_ptr_array_index(tenv->pinfo, pid);
}

gboolean
sh_typenv_is_pred_field0(sh_typenv_t* tenv, size_t pid, size_t fid,
    gboolean dir)
{
  // get predicate definition
  sh_pred_t* pred = g_ptr_array_index(tenv->pinfo,pid);
  if (dir == TRUE)
    {
      // forward field, search in edges from in
      GPtrArray* edges_in = g_ptr_array_index(pred->def_mat,0);
      for (guint i = 0; i < edges_in->len; i++)
        {
          sh_prededge_t* pei = g_ptr_array_index(edges_in,i);
          if (pei->discr == SH_PFORM_PTO && (pei->info.fid == fid))
            return TRUE;
        }
    }
  else
    {
      // backward field, search in edges from the first existential var to in
      GPtrArray* edges_u = g_ptr_array_index(pred->def_mat,pred->arity);
      for (guint i = 0; i < edges_u->len; i++)
        {
          sh_prededge_t* pei = g_ptr_array_index(edges_u,i);
          if (pei->discr == SH_PFORM_PTO && (pei->dst == 0) && (pei->info.fid
              == fid))
            return TRUE;
        }
    }
  return FALSE;
}

int
sh_typenv_lookup_pred(sh_typenv_t* tenv, char* pname)
{
  assert(NULL != tenv);
  assert(NULL != pname);

  /* search the name in the dictionary */
  gpointer orig, lookup;
  int res;
  if (g_hash_table_lookup_extended(tenv->pdict, (gconstpointer) pname, &orig,
      &lookup))
    res = GPOINTER_TO_INT(lookup);
  else
    res = -1;

  return res;
}

size_t
sh_typenv_push_pred(sh_typenv_t* tenv, sh_pred_t* p)
{
  assert (NULL != tenv);
  assert (NULL != p);

  /* push record in the array */
  int pos = tenv->pinfo->len;
  g_ptr_array_add(tenv->pinfo, p);

  g_hash_table_insert(tenv->pdict, p->name, GINT_TO_POINTER(pos));

  return (size_t) pos;
}
