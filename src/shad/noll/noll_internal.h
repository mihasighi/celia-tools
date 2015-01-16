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
 * noll_internal.h: internal representation of noll graphs and
 *                  additional operations called by the API implementation
 */

#ifndef _NOLL_INT_H_
#define _NOLL_INT_H_

#include <stdlib.h>
#include <stdio.h>
#include <glib.h>

#include "sh_manager.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ********************************************************************** */
/* I. Types */
/* ********************************************************************** */

/* ====================================================================== */
/* I.1 Edges */
/* ====================================================================== */

typedef enum {
	NOLL_EDGE_DIFF = 0, NOLL_EDGE_PTO, NOLL_EDGE_LS, NOLL_EDGE_UNKNOWN
/* NOT TO BE USED */
} noll_edge_e;

typedef struct noll_edge_s {
	size_t eid; /* unique identifier */
	size_t vrt_src; /* source node */
	size_t vrt_dst; /* destination node */
	noll_edge_e kind; /* kind of edge */
	union {
		size_t fid; /* for NOLL_EDGE_PTO */
		struct {
			size_t pid; /* predicate id in manager */
			GArray* ngb; /* neighbors vertex-ids */
			gboolean isrev; /* reverse of another edge */
		} pred; /* for NOLL_EDGE_LS */
	} info;
} noll_edge_t;

#define NOLL_EDGE_ID_UNKNOWN ((size_t) -1)

/* ====================================================================== */
/* I.2 Set of location variables */
/* ====================================================================== */

typedef struct noll_svar_s {
	char* name; /* optional */
	GArray* tids; /* collection of type-ids of locations, sorted array */
	size_t bto_eid; /* bound to edge-id */
} noll_svar_t;

/* ====================================================================== */
/* I.3 Set of location constraints */
/* ====================================================================== */

typedef enum {
	NOLL_OP_IN, NOLL_OP_SUBSETEQ, NOLL_OP_EQ, NOLL_OP_UNKNOWN
} noll_op_e;

typedef struct noll_share_s {
	noll_op_e op;
	size_t tid; /* type of this constraint */
	/* Term 1 */
	/* components are NULL if they are empty */
	GArray* vars_t1; /* array of prg var-ids sorted */
	GArray* lvars_t1; /* array of local var-ids sorted */
	GArray* svars_t1; /* array of location var-ids sorted */
	/* Term t2 */
	/* components are NULL if they are empty */
	GArray* vars_t2; /* array of prg var-ids sorted */
	GArray* lvars_t2; /* array of local var-ids sorted */
	GArray* svars_t2; /* array of location var-ids sorted */
} noll_share_t;

/* ====================================================================== */
/* I.4 Graph and abstract elements */
/* ====================================================================== */

/* Noll graph, an element of the abstract domain */
typedef struct noll_graph_s {
	/* Variables */
	size_t frame; /* index of the frame of program variables, in manager */
	GArray* lvars; /* array of local vars, mapping to their type-id */
	GPtrArray* svars; /* array of set of loc variables of type noll_svar_t */
	/* Vertices */
	size_t vrt_size; /* number of vertices in the graph, including nil=0 */
	GArray* vrt_type; /* type-id of each vertex in the graph */
	GArray* var2vrt; /* mapping of vars in frame and lvars to vertices */
	/* notice that var of id vid is at position vid+1 and null is at position 0 */
	gboolean isprecise; /* TRUE if it does not contain true node */

	/* Edges */
	GPtrArray* e_mat; /* adjacency matrix edges[v] = array of edge ids from v */
	GPtrArray* e_rmat; /* reverse adjacency matrix edges[v] = array of edge ids to v */
	GPtrArray* edges; /* array indexed by edge ids contains pointers to edges */

	/* Vertex and edge strong separation constraints
	 * represented by a up to diagonal triangular matrix
	 * of size given by the number of edges */
	GPtrArray* ssep; /* array indexed by edges, regions of a sorted array of edge_id */

	/* Vertex and edge sharing constraints */
	GPtrArray* wsep; /* array of noll_share_t sorted by first term */

} noll_graph_t;

/* Abstract domain, a set of graphs */
struct noll_val_s {
	size_t frame;
	GPtrArray* set; /* array of pointers to graphs */
};

/* in noll_fun.h */
/* typedef struct noll_val_s noll_val_t; */

/* ====================================================================== */
/* I.5 Graph homomorphism */
/* ====================================================================== */

typedef struct noll_hom_s {
	GArray* vrt; /* indexed by V(G1), elements in V(G2) */
	GPtrArray* used; /* indexed by E(G1), elements are arrays of noll_edgflds_t* */
} noll_hom_t;

typedef struct noll_edgflds_s {
	size_t eid; /* edge of pair */
	GArray* flds; /* array of fld-id sorted */
} noll_edgflds_t;

/* ********************************************************************** */
/* II. Constructor, getters, tests and property extraction */
/* ********************************************************************** */

/* see noll_build.c */

/* Manipulation of glib arrays */
#define g_ptr_array_copy(dst,src,copyFun,freeFun) \
  if (src == NULL) dst = NULL; \
  else { \
  dst = g_ptr_array_new_with_free_func((GDestroyNotify) freeFun);  \
  for (guint i = 0; i < src->len; i++) \
  g_ptr_array_add(dst, copyFun(g_ptr_array_index(src,i))); \
  }

#define g_array_copy(dst,src,ty)	\
  if (src == NULL) dst = NULL; \
  else { \
    dst = g_array_new (FALSE,FALSE,sizeof(ty));	\
    for (guint i = 0; i < src->len; i++) {	      \
    ty v = g_array_index(src,ty,i); \
    g_array_append_val(dst, v); }   \
  }

#define g_array_cmp(r,a,b,ty)	\
  if (a == NULL && b == NULL) r = 0; \
  else if (a == NULL && b != NULL) r = -1; \
  else if (a != NULL && b == NULL) r = 1; \
  else { \
    guint i = 0;	\
    r = 0; \
    for (; i < a->len && i < b->len && r==0; i++) {	      \
       ty ai = g_array_index(a,ty,i); \
       ty bi = g_array_index(b,ty,i); \
       if (ai < bi) r = -1; \
       else if (ai > bi) r = 1; \
    } \
    if (r == 0) { \
     if (i >= a->len && i < b->len) r = -1; \
     else if (i < a->len && i >= b->len) r = 1; \
    } \
  }

/* ============================================================ */
/* II.1 Basic constructors */
/* ============================================================ */

/* see noll_build.c */

/* Edges */
noll_edge_t* noll_edge_new(size_t id, size_t src, size_t dst, noll_edge_e kind);
void noll_edge_free(noll_edge_t* a);
noll_edge_t* noll_edge_copy(noll_edge_t* a);
int noll_edge_cmp(noll_edge_t* a, noll_edge_t* b);

/* Set of location variables */
noll_svar_t* noll_svar_new(char* name, GArray* tids, size_t e);
void noll_svar_free(noll_svar_t* a);
noll_svar_t* noll_svar_copy(noll_svar_t* a);

/* Sharing constraints */
noll_share_t* noll_share_new(noll_op_e op);
noll_share_t* noll_share_copy(noll_share_t* a);
void noll_share_free(noll_share_t* a);
int noll_share_cmp(noll_share_t* a, noll_share_t* b);

/* Graph */
/*noll_graph_t* noll_graph_alloc_null(sh_manager_t* man);
 Build the empty graph for the program variables,
 * only nil node labeled by all vars */
noll_graph_t* noll_graph_emp(sh_manager_t* man, size_t fid);
/* Build he graph corresponding to the emp heap */
noll_graph_t* noll_graph_true(sh_manager_t* man, size_t fid);
/* Build he graph corresponding to the true heap */
noll_graph_t* noll_graph_copy(noll_graph_t* a);
/* copy function */
void noll_graph_free(noll_graph_t* a);
/* deallocation function */

/* Values */
void noll_free_set(noll_val_t* a);
/* Free only the set structure, the graphs are re-used */

/* Homomorphism */
noll_hom_t* noll_hom_alloc(size_t vrt, size_t edg);
void noll_hom_free(noll_hom_t* a);
/* Allocation with size of mappings, deallocation */
/* ============================================================ */
/* II.2 Getters */
/* ============================================================ */
/* in noll_query.c */
size_t noll_graph_get_edge_from(noll_graph_t* a, size_t vrt, noll_edge_e kind,
		size_t id);
/* Return the id of the edge from vrt of kind and
 * info (field id or pred id) */
size_t noll_graph_get_label(noll_graph_t* a, size_t var);
/* Return the vertex labeled by var */

/* ============================================================ */
/* II.3 Tests */
/* ============================================================ */

/* see noll_query.c */

gboolean noll_is_sat(sh_manager_t * man, noll_val_t * a);
/* does the abstract value is satisfiable? */

gboolean noll_graph_is_sat(sh_manager_t * man, noll_graph_t * a);
/* does the abstract value is satisfiable? */

gboolean noll_graph_is_emp(sh_manager_t * man, noll_graph_t * a);
/* top check */

gboolean noll_graph_is_true(sh_manager_t * man, noll_graph_t * a);
/* top check */

gboolean noll_graph_is_le(sh_manager_t * man, noll_graph_t * a1,
		noll_graph_t * a2);
/* inclusion check */

gboolean noll_graph_is_eq(sh_manager_t * man, noll_graph_t * a1,
		noll_graph_t * a2);
/* equality check */

noll_hom_t* noll_graph_is_homomorphic(sh_manager_t * man, noll_graph_t * a1,
		noll_graph_t * a2);
/* homomorphism search */

gboolean noll_hom_is_isom(noll_hom_t* h);
/* isomorphism test */

gboolean noll_hom_is_consistent_diff(noll_hom_t* h, noll_graph_t* a1,
		noll_graph_t* a2);
/* h respects difference edges? */

gboolean noll_hom_is_consistent_ssep(noll_hom_t* h, noll_graph_t* a1,
		noll_graph_t* a2);
/* h respects separation constraints? */

/* ********************************************************************** */
/* III. Operations */
/* ********************************************************************** */

/* ============================================================ */
/* III.1 Meet and join */
/* ============================================================ */

/* see noll_nary.c */

noll_val_t* noll_join_graph(sh_manager_t* man, gboolean destructive,
		noll_val_t* a, gboolean copy, noll_graph_t* g);
/* Join the graph to the set of graphs */

noll_graph_t* noll_graph_join(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a1, noll_graph_t* a2);
/* Join the two graphs to find upper bound */

/* Widening */

noll_graph_t *
noll_graph_fold(sh_manager_t * man,
		noll_graph_t * a, size_t annon, GArray *preds);

void apply_edge_renameing(noll_graph_t *g, GArray *a);

void replace_unfolding_with_pred(sh_manager_t *man, noll_graph_t *g,
		GPtrArray *redges, guint src, guint dst, guint pid);

GPtrArray * check_unfolding_of_pred(GArray* nodes,noll_graph_t* g,
		sh_manager_t* man, guint pid);

GPtrArray* get_mat(sh_manager_t* man, sh_pred_t* p);

gboolean gint_array_contains(GArray *a,gint val);

/* ============================================================ */
/* III.2 Assignment */
/* ============================================================ */

/* see noll_transfer.c */

noll_val_t* noll_graph_expand(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t x, GArray* offs, GArray* edges);
/* Expand graph a on dimension x with offsets given and
 * store the last edge for each materialization */

gboolean
noll_graph_expand_1(sh_manager_t* man, gboolean destructive, noll_graph_t* a,
		size_t vrt, size_t fld, noll_val_t* res, GArray* edges);
/* Expand graph at vertex vrt for field fld and
 * append resulting graphs to res and resulting edges at edges
 */

/* ============================================================ */
/* III.5 Graph operations */
/* ============================================================ */

/* see noll_build.c */

noll_graph_t* noll_graph_add_vertex(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t tid);
/* Add a new vertex of type tid */

noll_graph_t* noll_graph_set_label(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t var, size_t vrt);
/* Change label of var to vrt, junk detected */

noll_graph_t* noll_graph_set_edge_dst(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t eid, size_t vrt);
/* Change destination of edge eid to vrt, junk detected */

noll_graph_t* noll_graph_add_edge_diff(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t src, size_t dst);
/* Add a difference edge to a. */

noll_graph_t* noll_graph_add_edge_pto(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t src, size_t fid, size_t dst, size_t* eid);
/* Add a points to edge to a and return its identifier in eid */

noll_graph_t* noll_graph_remove_junk(sh_manager_t* man, noll_graph_t* a,
		size_t vrt);
/* Remove the junk from vertex vrt in a; return a if vrt is not junk.
 */

noll_graph_t* noll_graph_collapse_edge(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t eid, size_t* vrt_f, size_t* eid_f);
/* Collapse the (list segment) edge eid in a and update values
 * for followed vertex vrt_f and edge eid_f.
 */

noll_graph_t* noll_graph_unfold_start(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t eid);
/* Unfold the (list segment) edge eid in a from start */

noll_graph_t* noll_graph_unfold_end(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t eid);
/* Unfold the (list segment) edge eid in a from end */

noll_graph_t* noll_graph_unfold_middle(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t eid, size_t tid);
/* Unfold the (list segment) edge eid in a until the type tid */

GQueue* noll_hom_fill_pvars(sh_manager_t* man, noll_hom_t* h, noll_graph_t* a1,
		noll_graph_t* a2);
/* Fill h using the program vars labeling,
 * return the queue of vertices in a1 mapped by h */

gboolean noll_hom_fill_pto(noll_hom_t* h, GQueue* v1mapped, noll_graph_t* a1,
		noll_graph_t* a2);
/* Fill h using the PTO edges in a1,
 * return TRUE only if a consistent homomorphism is obtained */

GPtrArray*
noll_hom_fill_pred(noll_hom_t* h, noll_graph_t* a1, noll_graph_t* a2);
/* Fill h using the predicate edges in a1,
 * return the substitution of svars in a1 with terms in a2
 * or NULL, if not consistent homomorphism found */

/* ********************************************************************** */
/* IV. Printing */
/* ********************************************************************** */

void noll_graph_fdump(FILE * stream, sh_manager_t * man, noll_graph_t * a);
void noll_graph_fprint(FILE * stream, sh_manager_t * man, noll_graph_t * a);
void noll_share_fdump(FILE * stream, sh_manager_t * man, noll_share_t * a);
void noll_edge_fdump(FILE * stream, sh_manager_t * man, noll_edge_t * a);

#ifdef __cplusplus
}
#endif

#endif /* _NOLL_INT_H_ */
