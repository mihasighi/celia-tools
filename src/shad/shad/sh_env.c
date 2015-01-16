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
 * sh_env.c: management of the program stack, binding of vars to dimensions
 */

#include <assert.h>
#include "sh_env.h"
#include "sh_manager.h"

sh_var_t sh_var_null =
  { "NULL", SH_TYP_TVOID, 0 };

/* ********************************************************************** */
/* II. User Functions */
/* ********************************************************************** */

/* ====================================================================== */
/* II.1 Variables */
/* ====================================================================== */

sh_var_t*
sh_var_new(char* name, size_t frame, size_t type)
{
  assert (NULL != name);
  sh_var_t* new_v = g_new(sh_var_t,1);
  assert (NULL != new_v);
  new_v->name = g_strdup(name);
  new_v->type = type;
  new_v->frame = frame;
  return new_v;
}

void
sh_var_free(sh_var_t* a)
{
  assert (NULL != a);
  assert (NULL != a->name);

  g_free(a->name);
  g_free(a);
}

sh_var_t*
sh_var_copy(sh_var_t* a)
{
  if (a == NULL)
    return NULL;
  sh_var_t* new_v = g_new(sh_var_t,1);
  new_v->name = (a->name == NULL) ? NULL : g_strdup(a->name);
  new_v->type = a->type;
  new_v->frame = (a->frame <= 0) ? 0 : a->frame;
  return new_v;
}

int
sh_var_cmp(sh_var_t* a, sh_var_t* b)
{
  assert (NULL != a);
  assert (NULL != b);

  if (a->frame < b->frame)
    return -1;
  if (a->frame > b->frame)
    return 1;
  return g_strcmp0(a->name, b->name);
}

int
sh_var_comparefunc(sh_var_t** a, sh_var_t** b)
{
  return sh_var_cmp(*a,*b);
}

guint
sh_var_hash(sh_var_t* a)
{
  assert (NULL != a);
  return g_direct_hash(a);
}

void
sh_var_fdump(FILE* stream, sh_var_t* v, sh_typenv_t* tenv)
{
  if (v == NULL)
    {
      fprintf(stream, "(null)");
      return;
    }
  fprintf(stream, "(frame %d) %s : ", v->frame, v->name);
  if (tenv == NULL)
    fprintf(stream, "type-id %d", v->type);
  else
    sh_typenv_fdump_tinfo(stream, tenv, v->type);
  // fflush(stream);
}

/* ====================================================================== */
/* II.2 Frames */
/* ====================================================================== */

sh_frame_t*
sh_frame_new(sh_loc_t* start, sh_loc_t* end)
{
  assert (NULL != start);
  assert (NULL != end);

  sh_frame_t* new_frame = g_new(sh_frame_t,1);
  assert (NULL != new_frame);
  new_frame->start = start;
  new_frame->end = end;
  new_frame->vars
      = g_ptr_array_new_with_free_func((GDestroyNotify) sh_var_free);
  return new_frame;
}

void
sh_frame_free(sh_frame_t* f)
{
  assert (NULL != f);

  if (f->vars != NULL)
    g_ptr_array_free(f->vars, TRUE);
  f->vars = NULL;
  sh_loc_free(f->start);
  f->start = NULL;
  sh_loc_free(f->end);
  f->end = NULL;
  free(f);
}

void
sh_frame_push(sh_frame_t* f, sh_var_t* v)
{
  assert (NULL != f);
  assert (NULL != v);

  /* this check is done elsewhere, thus it is an internal error */
  assert (sh_frame_get_var(f,v->name) == NULL);

  g_ptr_array_add(f->vars, v);
  g_ptr_array_sort(f->vars, (GCompareFunc) sh_var_comparefunc);

  //#ifdef NDEBUG
  sh_frame_fdump(stdout, f, NULL, NULL);
  //#endif
}

sh_var_t*
sh_frame_get_var(sh_frame_t* f, char* vname)
{
  size_t vid = sh_frame_get_var_id(f, vname);
  if (vid == SH_VAR_NULL)
    return NULL;
  else
    return sh_frame_get_varinfo(f, vid);
}

size_t
sh_frame_get_var_id(sh_frame_t* f, char* vname)
{
  if (f == NULL || f->vars == NULL || f->vars->len == 0)
    return SH_VAR_NULL;
  /* binary search, f->vars shall be sorted */
  size_t lo = 0;
  size_t hi = f->vars->len - 1;
  while (lo <= hi)
    {
      size_t mid = lo + (hi - lo) / 2;
      sh_var_t* vmid = sh_frame_get_varinfo(f, mid);
      int cmp = g_strcmp0((char*) vname, (char*) vmid->name);
      if (cmp < 0)
        hi = mid - 1;
      else if (cmp > 0)
        lo = mid + 1;
      else
        return mid;
    }
  return SH_VAR_NULL;
}

char*
sh_frame_get_varname(sh_frame_t* f, size_t vid)
{
  sh_var_t* var = sh_frame_get_varinfo(f, vid);
  return var->name;
}

void
sh_frame_fdump(FILE* stream, sh_frame_t* f, sh_typenv_t* tenv,
    sh_filenv_t* fenv)
{
  if (f == NULL)
    {
      fprintf(stream, "(frame)null\n");
      return;
    }
  fprintf(stream, "Frame: ");
  /* Print block information */
  fprintf(stream, "[start=>");
  sh_loc_fprint(stream, fenv, f->start);
  fprintf(stream, ",end=>");
  sh_loc_fprint(stream, fenv, f->end);
  /* Print variables information */
  fprintf(stream, "] [[");
  if (f->vars != NULL)
    for (guint i = 0; i < f->vars->len; i++) {
      sh_var_fdump(stream, g_ptr_array_index(f->vars,i), tenv);
      fprintf(stream,", ");
    }
  fprintf(stream, "]]\n");
}

/* ====================================================================== */
/* II.3 Stacks */
/* ====================================================================== */

sh_stack_t*
sh_stack_alloc_empty(sh_manager_t* man)
{
  sh_stack_t* stack = g_new(sh_stack_t, 1);
  stack->frames
      = g_ptr_array_new_with_free_func((GDestroyNotify) sh_frame_free);
  stack->man = man;
  return stack;
}

void
sh_stack_free(sh_stack_t* s)
{
  assert (NULL != s);
  if (s->frames != NULL)
    g_ptr_array_free(s->frames, TRUE);
  s->man = NULL;
  g_free(s);
}

size_t
sh_stack_push_frame(sh_stack_t* s, sh_frame_t* f)
{
  assert (NULL != s);
  assert (NULL != f);
  size_t fid = (s->frames == NULL) ? 0 : s->frames->len;
  g_ptr_array_add(s->frames, f);
  return fid;
}

void
sh_stack_push_var(sh_stack_t* s, size_t f, sh_var_t* v)
{
  assert (NULL != s);
  assert (NULL != v);
  size_t max = (s->frames == NULL) ? 0 : s->frames->len;
  assert (f < max);
  v->frame = f;
  sh_frame_t* frame = g_ptr_array_index (s->frames, f);
  sh_frame_push(frame, v);
  return;
}

sh_var_t*
sh_stack_get_var(sh_stack_t* s, size_t fid, size_t vid)
{
  assert (NULL != s);
  size_t max = (s->frames == NULL) ? 0 : s->frames->len;
  if (fid >= max)
    return NULL;

  sh_frame_t* frame = g_ptr_array_index (s->frames, fid);
  if (frame == NULL || frame->vars == NULL || vid >= frame->vars->len)
    return NULL;
  sh_var_t* v = g_ptr_array_index (frame->vars, vid);
  return v;
}

void
sh_stack_fdump(FILE* stream, sh_stack_t* s)
{
  if (s == NULL)
    {
      fprintf(stream, "(stack)null\n");
      return;
    }
  fprintf(stream, "Stack: size=%d\n", (s->frames == NULL) ? 0 : s->frames->len);
  /* Print frame information */
  sh_typenv_t* tenv = sh_manager_get_typenv(s->man);
  sh_filenv_t* fenv = sh_manager_get_filenv(s->man);
  if (s->frames != NULL)
    for (guint i = 0; i < s->frames->len; i++)
      {
        fprintf(stream, "==== %d-", i);
        sh_frame_fdump(stream, g_ptr_array_index(s->frames,i), tenv, fenv);
      }
}
