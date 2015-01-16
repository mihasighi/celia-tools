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
 * sh_utils.h: utility functions for testing
 */

#ifndef _SH_UTILS_H_
#define _SH_UTILS_H_

#include <assert.h>
#include "sh_manager.h"
#include "sh_stmt.h"

#ifdef __cplusplus
extern "C"
  {
#endif

/* ********************************************************************** */
/* I. Types */
/* ********************************************************************** */

/* ICFG vertex */
typedef guint sh_icfg_vertex_t;

/* ICFG edge */
typedef struct sh_icfg_edge_s
{
  sh_icfg_vertex_t from_call; /* for return edges */
  sh_icfg_vertex_t to;
  gboolean isprj; /* project the local variables */
  sh_stmt_t* stmt;
} sh_icfg_edge_t;

/* ICFG vertex annotation */
typedef struct sh_icfg_info_s
{
  size_t frame; /* frame of this vertex */
  sh_loc_t* loc; /* of the vertex */
  void* val; /* computed invariant */
  GPtrArray* edges; /* array of sh_icfg_edge_t* */
} sh_icfg_info_t;


/* ********************************************************************** */
/* II. Functions */
/* ********************************************************************** */

void
sh_file_register(sh_manager_t* man, char* fname);
/* Register file with name fname */

void
sh_type_register_ls(sh_manager_t* man, char* fname, guint line, char* tname,
    char* nfname, char* dfname);
/* Register type for singly linked list and its predicate */

void
sh_type_register_dll(sh_manager_t* man, char* fname, guint line, char* tname,
    char* nfname, char* pfname, char* dfname);
/* Register type for doubly linked list and its predicate */

void
sh_pred_register_ls(sh_manager_t* man, char* pname, char* tname, char* fname);
/* Register predicate for the singly linked list type tname */

#ifdef __cplusplus
}
#endif

#endif /* _SH_UTILS_H_ */
