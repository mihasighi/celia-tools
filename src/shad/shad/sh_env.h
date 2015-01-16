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
 * sh_env.h: management of the program stack, binding of vars to dimensions
 */

#ifndef _SH_ENV_H_
#define _SH_ENV_H_

#include "sh_typ.h"

#ifdef __cplusplus
extern "C"
  {
#endif

/* ********************************************************************** */
/* I. Types */
/* ********************************************************************** */

/* ====================================================================== */
/* I.1 Variables */
/* ====================================================================== */

/* Variables are declared with their name, type and its definition frame */
typedef struct sh_var_t
{
  char* name;
  size_t type;
  size_t frame;
} sh_var_t;

/* Special variable NULL */
//extern  sh_var_t sh_var_null;
#define SH_VAR_NULL ((size_t) -1)

/* ====================================================================== */
/* I.2 Frame */
/* ====================================================================== */

/* Frames are arrays or variables; variables are uniquely identified 
 by its frame and its entry in the frame.
 */
typedef struct sh_frame_t
{
  sh_loc_t* start; /* start location */
  sh_loc_t* end; /* end location */
  GPtrArray* vars; /* sorted array of variables */
} sh_frame_t;

/* ====================================================================== */
/* I.3 Stack */
/* ====================================================================== */

/* Stacks are collection of frames; a frame identifier is the frame
 * number in the stack.
 */
typedef struct sh_stack_t
{
  GPtrArray* frames; /* array of frames */
  sh_manager_t* man;
} sh_stack_t;

/* ********************************************************************** */
/* II. User Functions */
/* ********************************************************************** */

/* ====================================================================== */
/* II.1 Variables */
/* ====================================================================== */

sh_var_t*
sh_var_new(char* name, size_t frame, size_t type);
/* Constructor a new variable */
void
sh_var_free(sh_var_t* a);
/* Destructor */
sh_var_t*
sh_var_copy(sh_var_t* v);
/* Copy a variable */
int
sh_var_cmp(sh_var_t* a, sh_var_t* b);
/* Comparison function on (frame, var) */
guint
sh_var_hash(sh_var_t* a);
/* Hash functions on (frame, type, var) */

void
sh_var_fdump(FILE* stream, sh_var_t* v, sh_typenv_t* tenv);
/* Dumping information about a variable */

/* ====================================================================== */
/* II.2 Frames */
/* ====================================================================== */

sh_frame_t*
sh_frame_new(sh_loc_t* start, sh_loc_t* end);
/* Constructor for a new frame, empty set of variables. */
void
sh_frame_free(sh_frame_t* f);
/* Destructor, free also the variables allocated. */
void
sh_frame_push(sh_frame_t* f, sh_var_t* v);
/* Add a variable to the frame, error if it exists. */

sh_var_t*
sh_frame_get_var(sh_frame_t* f, char* vname);
/* Returns the information about a variable, given its name */

static inline sh_var_t*
sh_frame_get_varinfo(sh_frame_t* f, size_t vid);
/* Returns the information about variable, given its id */
size_t
sh_frame_get_var_id(sh_frame_t* f, char* vname);
/* Returns the id of variable, given its name */

char*
sh_frame_get_varname(sh_frame_t* f, size_t vid);
/* Returns the name of variable, given its id */

static inline gboolean
sh_frame_is_var(sh_frame_t* f, size_t v);
/* Check that v is in frame */

void
sh_frame_fdump(FILE* stream, sh_frame_t* f, sh_typenv_t* tenv,
    sh_filenv_t* fenv);
/* Dump information about the frame */

/* ====================================================================== */
/* II.3 Stacks */
/* ====================================================================== */

sh_stack_t*
sh_stack_alloc_empty(sh_manager_t* man);
/* Constructor for a new stack, no frames. */
void
sh_stack_free(sh_stack_t* s);
/* Destructor, free also the frames. */
size_t
sh_stack_push_frame(sh_stack_t* s, sh_frame_t* f);
/* Add a frame, return its identifier. */

void
sh_stack_push_var(sh_stack_t* s, size_t fid, sh_var_t* v);
/* Add a variable to the frame */

sh_var_t*
sh_stack_get_var(sh_stack_t* s, size_t fid, size_t vid);
/* Returns the information about a variable,
 given its frame and its identifier */

static inline sh_frame_t*
sh_stack_get_frame(sh_stack_t* s, size_t fid);
/* Returns the information about a frame,
 given its frame and its identifier */

void
sh_stack_fdump(FILE* stream, sh_stack_t* s);
/* Dump information about the stack */

/* ====================================================================== */
/* Inline function definitions */
/* ====================================================================== */

static inline sh_var_t*
sh_frame_get_varinfo(sh_frame_t* f, size_t vid)
{
  return g_ptr_array_index(f->vars, vid);
}

static inline gboolean
sh_frame_is_var(sh_frame_t* f, size_t x)
{
  return (f != NULL) && (f->vars != NULL) && (x < f->vars->len || x
      == SH_VAR_NULL);
}

static inline sh_frame_t*
sh_stack_get_frame(sh_stack_t* s, size_t fid)
{
  return g_ptr_array_index(s->frames,fid);
}

#ifdef __cplusplus
}
#endif

#endif
