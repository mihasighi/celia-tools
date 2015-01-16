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
 * sh_loc.h: file names and locations types defined in the program and
 *           passed in the global manager
 */

#ifndef _SH_LOC_H_
#define _SH_LOC_H_

#include <stdlib.h>
#include <stdio.h>
#include <glib.h>
#include "shad.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ********************************************************************** */
/* I. Types */
/* ********************************************************************** */

/* ====================================================================== */
/* I.1 Environment of files */
/* ====================================================================== */

  /* File environment is an array of file names.
   * File ids are indexes in this array.
   */
  typedef struct sh_filenv_t {
    GPtrArray* filenames;
    guint max_files;
    sh_manager_t* man;
  } sh_filenv_t;
  
/* ====================================================================== */
/* I.2 Locations */
/* ====================================================================== */
 
  /* Definition of pointer types */
  typedef struct sh_loc_t {
    guint  fid;     /* file identifier, 0 is NA */
    guint  lin;     /* line number */
    guint  col;     /* column number */
  } sh_loc_t;

/* ********************************************************************** */
/* II. Functions */
/* ********************************************************************** */

/* ====================================================================== */
/* II.1 Memory management, constructors, destructors */
/* ====================================================================== */

  sh_filenv_t* sh_filenv_alloc_empty(sh_manager_t* man);
  /* Build the empty file environment */
  void sh_filenv_free(sh_filenv_t* fenv);
  /* Free the environment
     (the structure itself and the memory pointed to by fields) */

  sh_loc_t* sh_loc_alloc(guint fid, gulong line, gulong col);
  /* Constructor for locations */
  sh_loc_t* sh_loc_copy (sh_loc_t* l);
  /* Copy constructor */
  void sh_loc_free(sh_loc_t* loc);
  /* Destructor for locations */


/* ====================================================================== */
/* II.2 Printing */
/* ====================================================================== */

  void sh_filenv_fdump(FILE* stream, sh_filenv_t* fenv);
  /* Print (in debug mode) a file environment */

  void sh_loc_fprint(FILE* stream, sh_filenv_t* fenv, sh_loc_t* l);
  /* Print predicate information for predicate index pid */

/* ====================================================================== */
/* II.3 Tests */
/* ====================================================================== */

/* Locations */
  gint sh_loc_cmp(sh_loc_t* a, sh_loc_t* b);
  /* Comparison between location; lexicographic fid,lin,col;
     -1, 0, 1 if a < b, a == b, a > b
  */
  gboolean sh_loc_cmp_file(sh_loc_t* a, sh_loc_t* b);
  /* Comparison for exactly the same file   */

/* ====================================================================== */
/* II.4 Getters and setters */
/* ====================================================================== */

  char* sh_filenv_get_fname (sh_filenv_t* fenv, guint fid);
  /* Lookup for fid and return the pointer to the name */

  guint sh_filenv_push_fname (sh_filenv_t* fenv, char* fname);
  /* Lookup for fname and add (a copy) it if not present; 
     return the identifier */
 

#ifdef __cplusplus
}
#endif

#endif /* _SH_LOC_H_ */
