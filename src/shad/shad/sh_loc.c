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
 * sh_loc.c: file names and locations types defined in the program and
 *           passed in the global manager
 */

#include <assert.h>
#include "sh_loc.h"

/* ********************************************************************** */
/* II. Functions */
/* ********************************************************************** */

/* ====================================================================== */
/* II.1 Memory management, constructors, destructors */
/* ====================================================================== */

sh_filenv_t*
sh_filenv_alloc_empty(sh_manager_t* man)
{
  sh_filenv_t* fenv;

  fenv = g_new (sh_filenv_t, 1);
  assert (fenv != NULL);
  fenv->filenames = g_ptr_array_new_with_free_func (free);
  assert (fenv->filenames != NULL);
  fenv->max_files = 0; 
  fenv->man = man;

  return fenv;
}

void 
sh_filenv_free(sh_filenv_t* fenv)
{
  assert (fenv != NULL);
  
  if (fenv->filenames != NULL)
    g_ptr_array_free(fenv->filenames, TRUE);
  free (fenv);
}

sh_loc_t* 
sh_loc_alloc(guint fid, gulong lin, gulong col)
{
  sh_loc_t* loc;

  loc = g_new(sh_loc_t,1);
  assert (loc != NULL);
  loc->fid = fid;
  loc->lin = lin;
  loc->col = col;
  return loc;  
}

sh_loc_t* 
sh_loc_copy (sh_loc_t* l)
{
  if (l == NULL) return NULL;

  sh_loc_t* new_l = sh_loc_alloc(l->fid, l->lin, l->col);
  return new_l;
}

void 
sh_loc_free(sh_loc_t* loc)
{
  assert (loc != NULL);
  free(loc);
}


/* ====================================================================== */
/* II.2 Printing */
/* ====================================================================== */

void 
sh_filenv_fdump(FILE* stream, sh_filenv_t* fenv)
{
  assert (stream != NULL);
  fprintf (stream, "File environment: (size = %d)\n", 
	   ((fenv != NULL) ? fenv->max_files : 0));
  if (fenv == NULL || fenv->filenames == NULL)
    fprintf (stream, "empty\n");
  else
    for (guint i = 0; i < fenv->max_files; i++)
      fprintf (stream, "(%d:) %s\n", i,
	       (char*) g_ptr_array_index(fenv->filenames, i));
}

void 
sh_loc_fprint(FILE* stream, sh_filenv_t* fenv, sh_loc_t* l)
{
  assert (stream != NULL);

  if (l == NULL)
    {
      fprintf (stream, "(null)");
      return;
    }

  char* fname = sh_filenv_get_fname(fenv, l->fid); /* DO NOT FREE */
  fprintf (stream, "(%s:%d,%d)", 
	   (fname == NULL) ? "null" : fname,
	   l->lin, l->col);
}

/* ====================================================================== */
/* II.3 Tests */
/* ====================================================================== */

gint 
sh_loc_cmp(sh_loc_t* a, sh_loc_t* b)
{
  assert (a!=NULL);
  assert (b!=NULL);

  if (a->fid < b->fid)
    return -1;
  if (a->fid > b->fid)
    return 1;

  if (a->lin < b->lin)
    return -1;
  if (a->lin > b->lin)
    return 1;

  if (a->col < b->col)
    return -1;
  if (a->col > b->col)
    return 1;
  
  return 0;
}
 
gboolean 
sh_loc_cmp_file(sh_loc_t* a, sh_loc_t* b)
{  
  assert (a!=NULL);
  assert (b!=NULL);

  return (a->fid != b->fid) ? FALSE : TRUE;
}

/* ====================================================================== */
/* II.4 Getters and setters */
/* ====================================================================== */

char* 
sh_filenv_get_fname (sh_filenv_t* fenv, guint fid)
{
  if (fenv == NULL || 
      fenv->filenames == NULL ||
      fenv->max_files <= fid)
    return NULL;

  char* name = (char*) g_ptr_array_index(fenv->filenames, fid);
  return name;
}

guint 
sh_filenv_push_fname (sh_filenv_t* fenv, char* fname)
{
  assert (fenv != NULL);
  assert (fname != NULL);

  for (guint i = 0; i < fenv->max_files; i++)
    {
      char* fname_i = g_ptr_array_index (fenv->filenames, i);
      if (g_strcmp0(fname, fname_i) == 0)
	return i;
    }
  // else, fname is not present
  char * lfname = g_strdup(fname);
  assert (NULL != lfname);
  guint fid = fenv->max_files;
  g_ptr_array_add (fenv->filenames, lfname);
  fenv->max_files++;
  return fid;
}
