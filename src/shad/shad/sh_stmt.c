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
 * sh_stmt.c: management of the program statements
 */

#include <assert.h>
#include <string.h>
#include "sh_stmt.h"

/* ********************************************************************** */
/* II. User Functions */
/* ********************************************************************** */

sh_stmt_t*
sh_stmt_alloc(sh_manager_t* man, sh_stmt_e kind, sh_loc_t* loc)
{
  assert (NULL != man);
  assert (NULL != loc);

  sh_stmt_t* new_s = g_new(sh_stmt_t, 1);
  assert (NULL != new_s);
  memset(new_s, 0, sizeof(sh_stmt_t));
  new_s->kind = kind;
  new_s->loc = loc;
  int fid = sh_manager_get_frame_of_loc(man, loc);
  if (fid >= 0)
    {
      new_s->frame = g_new(size_t, 1);
      *(new_s->frame) = fid;
    }
  else
    new_s->frame = NULL;
  return new_s;
}

void
sh_stmt_free(sh_stmt_t* a)
{
  assert (NULL != a);

  sh_loc_free(a->loc);

  if (a->kind == SH_STMT_PCALL)
    {
      g_free(a->info.call.pname);
      g_array_free(&a->info.call.args, FALSE);
    }
  g_free(a);
}

sh_stmt_t*
sh_stmt_new_skip(sh_manager_t* man, sh_loc_t* loc)
{
  assert (NULL != man);
  assert (NULL != loc);

  sh_stmt_t* r = sh_stmt_alloc(man, SH_STMT_SKIP, loc);
  assert (r != NULL);
  assert (r->frame != NULL);

  return r;
}

sh_stmt_t*
sh_stmt_new_assume(sh_manager_t* man, sh_loc_t* loc, sh_stmt_e kind, size_t x,
    size_t y)
{
  assert (NULL != man);
  assert (NULL != loc);
  assert (kind == SH_STMT_ASSUME_EQ || kind == SH_STMT_ASSUME_NE);

  sh_stmt_t* r = sh_stmt_alloc(man, kind, loc);
  assert (r != NULL);
  assert (r->frame != NULL);

  /* check x and y before continue */
  sh_frame_t* frame = g_ptr_array_index((sh_manager_get_varenv(man))->frames,
      *(r->frame));
  assert (sh_frame_is_var(frame,x) == TRUE);
  assert (sh_frame_is_var(frame,y) == TRUE);

  /* fill r */
  r->info.binary.left = x;
  r->info.binary.right = y;

  return r;
}

sh_stmt_t*
sh_stmt_new_new(sh_manager_t* man, sh_loc_t* loc, size_t x)
{
  assert (NULL != man);
  assert (NULL != loc);

  sh_stmt_t* r = sh_stmt_alloc(man, SH_STMT_NEW, loc);
  assert (r != NULL);
  assert (r->frame != NULL);

  /* check x */
  sh_frame_t* frame = g_ptr_array_index((sh_manager_get_varenv(man))->frames,
      *(r->frame));
  assert (sh_frame_is_var(frame,x) == TRUE);
  assert (x != SH_VAR_NULL);

  /* fill r */
  sh_var_t* v_x = sh_frame_get_varinfo(frame, x);
  r->tid = v_x->type;
  r->info.binary.left = x;

  return r;
}

sh_stmt_t*
sh_stmt_new_free(sh_manager_t* man, sh_loc_t* loc, size_t x)
{
  assert (NULL != man);
  assert (NULL != loc);

  sh_stmt_t* r = sh_stmt_alloc(man, SH_STMT_FREE, loc);
  assert (r != NULL);
  assert (r->frame != NULL);

  /* check x */
  sh_frame_t* frame = g_ptr_array_index((sh_manager_get_varenv(man))->frames,
      *(r->frame));
  assert (sh_frame_is_var(frame,x) == TRUE);
  assert (x != SH_VAR_NULL);

  /* fill r */
  sh_var_t* v_x = sh_frame_get_varinfo(frame, x);
  r->tid = v_x->type;
  r->info.binary.left = x;

  return r;
}

sh_stmt_t*
sh_stmt_new_assign_x_y(sh_manager_t* man, sh_loc_t* loc, size_t x, size_t y)
{
  assert (NULL != man);
  assert (NULL != loc);

  sh_stmt_t* r = sh_stmt_alloc(man, SH_STMT_ASSIGN, loc);
  assert (r != NULL);
  assert (r->frame != NULL);

  /* check x and y before continue */
  sh_frame_t* frame = g_ptr_array_index((sh_manager_get_varenv(man))->frames,
      *(r->frame));
  assert (sh_frame_is_var(frame,x) == TRUE);
  assert (x != SH_VAR_NULL);
  assert (sh_frame_is_var(frame,y) == TRUE);

  /* fill r */
  sh_var_t* v_x = sh_frame_get_varinfo(frame, x);
  r->tid = v_x->type;
  r->info.binary.left = x;
  r->info.binary.offset_l = NULL;
  r->info.binary.right = y;
  r->info.binary.offset_r = NULL;

  return r;
}

sh_stmt_t*
sh_stmt_new_assign_x_f_y(sh_manager_t* man, sh_loc_t* loc, size_t x, size_t f,
    size_t y)
{
  assert (NULL != man);
  assert (NULL != loc);

  sh_stmt_t* r = sh_stmt_alloc(man, SH_STMT_ASSIGN, loc);
  assert (r != NULL);
  assert (r->frame != NULL);

  /* check x and y before continue */
  sh_frame_t* frame = g_ptr_array_index((sh_manager_get_varenv(man))->frames,
      *(r->frame));
  assert (sh_frame_is_var(frame,x) == TRUE);
  assert (x != SH_VAR_NULL);
  assert (sh_frame_is_var(frame,y) == TRUE);

  /* fill r */
  sh_fld_t* fld = sh_typenv_get_field(sh_manager_get_typenv(man), f);
  r->tid = fld->tid;
  r->info.binary.left = x;
  r->info.binary.offset_l = g_array_new(FALSE, FALSE, sizeof(size_t));
  g_array_append_val(r->info.binary.offset_l,f);
  r->info.binary.right = y;

  return r;
}

sh_stmt_t*
sh_stmt_new_assign_x_y_f(sh_manager_t* man, sh_loc_t* loc, size_t x, size_t y,
    size_t f)
{
  assert (NULL != man);
  assert (NULL != loc);

  sh_stmt_t* r = sh_stmt_alloc(man, SH_STMT_ASSIGN, loc);
  assert (r != NULL);
  assert (r->frame != NULL);

  /* check x and y before continue */
  sh_frame_t* frame = g_ptr_array_index((sh_manager_get_varenv(man))->frames,
      *(r->frame));
  assert (sh_frame_is_var(frame,x) == TRUE);
  assert (x != SH_VAR_NULL);
  assert (sh_frame_is_var(frame,y) == TRUE);
  assert (y != SH_VAR_NULL);

  /* fill r */
  sh_var_t* v_x = sh_frame_get_varinfo(frame, x);
  r->tid = v_x->type;
  r->info.binary.left = x;
  r->info.binary.offset_l = NULL;
  r->info.binary.right = y;
  r->info.binary.offset_r = g_array_new(FALSE, FALSE, sizeof(size_t));
  g_array_append_val(r->info.binary.offset_r,f);

  return r;
}

void
sh_stmt_fdump(FILE* stream, sh_manager_t* man, sh_stmt_t* s)
{
  assert (NULL != stream);
  assert (NULL != man);

  if (s == NULL)
    {
      fprintf(stream, "(null): (null)\n");
      return;
    }

  sh_loc_fprint(stream, sh_manager_get_filenv(man), s->loc);
  fprintf(stream, ": (typ-%d): ", s->tid);
  switch (s->kind)
    {
  case SH_STMT_SKIP:
    fprintf(stream, "skip");
    break;
  case SH_STMT_ASSUME_EQ:
  case SH_STMT_ASSUME_NE:
    {
      sh_var_t* x = sh_stack_get_var(sh_manager_get_varenv(man), *s->frame,
          s->info.binary.left);
      sh_var_t* y = sh_stack_get_var(sh_manager_get_varenv(man), *s->frame,
          s->info.binary.right);
      fprintf(stream, "assume (%s %s= %s)", x->name, (s->kind
          == SH_STMT_ASSUME_EQ) ? "=" : "!", y->name);
      break;
    }
  case SH_STMT_NEW:
    {
      sh_var_t* x = sh_stack_get_var(sh_manager_get_varenv(man), *s->frame,
          s->info.binary.left);
      sh_typ_t* ty = sh_typenv_get_type(sh_manager_get_typenv(man), x->type);
      fprintf(stream, "%s = new %s()", x->name, ty->name);
      break;
    }
  case SH_STMT_FREE:
    {
      sh_var_t* x = sh_stack_get_var(sh_manager_get_varenv(man), *s->frame,
          s->info.binary.left);
      fprintf(stream, "free(%s)", x->name);
      break;
    }
  case SH_STMT_ASSIGN:
    {
      sh_var_t* x = sh_stack_get_var(sh_manager_get_varenv(man), *s->frame,
          s->info.binary.left);
      sh_var_t* y = sh_stack_get_var(sh_manager_get_varenv(man), *s->frame,
          s->info.binary.right);
      /* print left hand side */
      fprintf(stream, "%s", x->name);
      if (s->info.binary.offset_l != NULL)
        for (guint i = 0; i < s->info.binary.offset_l->len; i++)
          {
            size_t offi = g_array_index(s->info.binary.offset_l,size_t,i);
            sh_fld_t* fldi = sh_typenv_get_field(sh_manager_get_typenv(man),
                offi);
            fprintf(stream, ".%s", fldi->name);
          }

      fprintf(stream, " = ");
      /* print right hand side */
      if (y == NULL)
        fprintf(stream, "NULL");
      else
        fprintf(stream, "%s", y->name);
      if (s->info.binary.offset_r != NULL)
        for (guint i = 0; i < s->info.binary.offset_r->len; i++)
          {
            size_t offi = g_array_index(s->info.binary.offset_r,size_t,i);
            sh_fld_t* fldi = sh_typenv_get_field(sh_manager_get_typenv(man),
                offi);
            fprintf(stream, ".%s", fldi->name);
          }
      break;
    }
  case SH_STMT_PCALL:
    {
      if (s->info.call.ret != NULL)
        {
          sh_var_t* x = sh_stack_get_var(sh_manager_get_varenv(man), *s->frame,
              *s->info.call.ret);
          fprintf(stream, "%s = ", x->name);
        }
      fprintf(stream, "%s(", s->info.call.pname);
      for (guint i = 0; i < s->info.call.args.len; i++)
        {
          sh_var_t* x = sh_stack_get_var(sh_manager_get_varenv(man), *s->frame,
              g_array_index(&(s->info.call.args),size_t,i));
          gboolean isval = g_array_index(&(s->info.call.argk),gboolean,i);
          fprintf(stream, " %s%s,", (isval == TRUE) ? "" : "&", x->name);
        }
      fprintf(stream, ")");
      break;
    }
  default:
    fprintf(stream, "unknown");
    }
  fprintf(stream, "\n");
}
