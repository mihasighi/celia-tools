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
 * sh_utils.c: utility functions for testing
 */

#include <assert.h>
#include <glib.h>
#include <string.h>
#include "sh_manager.h"

/* ********************************************************************** */
/* I. Globals */
/* ********************************************************************** */

/* ********************************************************************** */
/* II. Functions */
/* ********************************************************************** */

/* ====================================================================== */
/* II.1 Files */
/* ====================================================================== */

void
sh_file_register(sh_manager_t* man, char* fname)
{
  assert (NULL != man);
  assert (NULL != fname);

  sh_filenv_push_fname(sh_manager_get_filenv(man), fname);
}

void
sh_type_register_ls(sh_manager_t* man, char* fname, guint line, char* tname,
    char* nfname, char* dfname)
{
  assert (NULL != man);
  assert (NULL != fname);
  assert (NULL != tname);
  assert (NULL != nfname);
  assert (NULL != dfname);

  sh_filenv_t* fenv = sh_manager_get_filenv(man);
  guint file = sh_filenv_push_fname(fenv, fname);

  sh_typenv_t* tenv = sh_manager_get_typenv(man);
  /* type of data field */
  size_t tid_int = sh_typenv_add_type_int(tenv, SH_IKIND_INT);
  /* type for record list = the struct list */
  sh_loc_t* loc_str_list = sh_loc_alloc(file, line, 0);
  size_t tid_str_list = sh_typenv_add_type_record(tenv, loc_str_list, tname);
  /* type for recursive field = struct list* */
  sh_loc_t* loc_ref_list = sh_loc_alloc(file, line + 2, 0);
  size_t tid_ref_list = sh_typenv_add_type_ptr0(tenv, loc_str_list,
      tid_str_list, 1);

  /* register field data */
  sh_loc_t* loc_fld_data = sh_loc_alloc(file, line + 1, 0);
  size_t fid_data = sh_typenv_add_field0(tenv, loc_fld_data, dfname,
      tid_str_list, tid_int);

  /* register field next */
  sh_loc_t* loc_fld_next = sh_loc_alloc(file, line + 2, 0);
  size_t fid_next = sh_typenv_add_field0(tenv, loc_fld_next, nfname,
      tid_str_list, tid_ref_list);
  sh_typenv_set_field_recursive(tenv, fid_next, TRUE);
  sh_typenv_set_record_recursive(tenv, tid_str_list, 1);

  /* register predicate for this type */
  sh_loc_t* loc_pred_ls = sh_loc_alloc(file, line + 3, 0);
  char* pname = g_new(char, strlen(tname) + 4);
  sprintf(pname, "ls_%s", tname);
  size_t pid_ls = sh_typenv_add_pred_ls(tenv, loc_pred_ls, pname, fid_next);
  g_free(pname); // strdup in the function above

}

void
sh_pred_register_ls(sh_manager_t* man, char* pname, char* tname, char* fname)
{
  assert (NULL != man);
  assert (NULL != pname);
  assert (NULL != tname);
  assert (NULL != fname);

  sh_typenv_t* tenv = sh_manager_get_typenv(man);

  /* search field and its type, needed by the predicate definition */
  int tid = sh_typenv_lookup_type(tenv, tname);
  assert (-1 != tid);
  sh_typ_t* t = sh_typenv_get_type(tenv, tid);
  /* t shall be a pointer to some struct ... */
  glong
      str_tid =
          (t->discr == SH_TYP_TPTR && t->info.ptr.deref == 1) ? ((glong) t->info.ptr.base_tid)
              : ((glong) -1);
  if (str_tid == -1)
    {
      fprintf(stderr, "Bad type (no *) for predicate sll! Quit.\n");
      return;
    }
  sh_typ_t* str_t = sh_typenv_get_type(tenv, (size_t) str_tid);
  while (str_t->discr == SH_TYP_TNAMED)
    {
      str_tid = str_t->info.named.base_tid;
      str_t = sh_typenv_get_type(tenv, (size_t) str_tid);
    }
  if (str_t->discr != SH_TYP_TRECORD)
    {
      fprintf(stderr, "Bad type (no struct) for predicate sll! Quit.\n");
      return;
    }
  int fid = sh_typenv_lookup_field(tenv, fname, str_t->name);
  assert (-1 != fid);
  sh_loc_t* loc_pred = sh_loc_copy(t->loc);
  assert (loc_pred != NULL);
  size_t pid = sh_typenv_add_pred_ls(tenv, loc_pred, pname, fid);
  return;
}

void
sh_type_register_dll(sh_manager_t* man, char* fname, guint line, char* tname,
    char* nfname, char* pfname, char* dfname)
{
  assert (NULL != man);
  assert (NULL != fname);
  assert (NULL != tname);
  assert (NULL != nfname);
  assert (NULL != dfname);

  sh_filenv_t* fenv = sh_manager_get_filenv(man);
  guint file = sh_filenv_push_fname(fenv, fname);

  sh_typenv_t* tenv = sh_manager_get_typenv(man);
  /* type of data field */
  size_t tid_int = sh_typenv_add_type_int(tenv, SH_IKIND_INT);
  /* type for record list = the struct list */
  sh_loc_t* loc_str_list = sh_loc_alloc(file, line, 0);
  size_t tid_str_list = sh_typenv_add_type_record(tenv, loc_str_list, tname);
  /* type for recursive field = struct list* */
  sh_loc_t* loc_ref_list = sh_loc_alloc(file, line + 2, 0);
  size_t tid_ref_list = sh_typenv_add_type_ptr0(tenv, loc_str_list,
      tid_str_list, 1);

  /* register field data */
  sh_loc_t* loc_fld_data = sh_loc_alloc(file, line + 1, 0);
  size_t fid_data = sh_typenv_add_field0(tenv, loc_fld_data, dfname,
      tid_str_list, tid_int);

  /* register field next */
  sh_loc_t* loc_fld_next = sh_loc_alloc(file, line + 2, 0);
  size_t fid_next = sh_typenv_add_field0(tenv, loc_fld_next, nfname,
      tid_str_list, tid_ref_list);
  sh_typenv_set_field_recursive(tenv, fid_next, TRUE);

  /* register field prev */
  sh_loc_t* loc_fld_prev = sh_loc_alloc(file, line + 3, 0);
  size_t fid_prev = sh_typenv_add_field0(tenv, loc_fld_prev, pfname,
      tid_str_list, tid_ref_list);
  sh_typenv_set_field_recursive(tenv, fid_prev, TRUE);

  sh_typenv_set_record_recursive(tenv, tid_str_list, 2);

  /* register predicate for this type */
  sh_loc_t* loc_pred_dll = sh_loc_alloc(file, line + 4, 0);
  char* pname = g_new(char, strlen(tname) + 5);
  sprintf(pname, "dll_%s", tname);
  size_t pid_dll = sh_typenv_add_pred_dll(tenv, loc_pred_dll, pname, fid_next,
      fid_prev);
  g_free(pname); // strdup in the function above

}

