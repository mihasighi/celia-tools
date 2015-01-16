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
 * noll_build.c: Implementation of functions related to the
 *               data values allocation/copy/free
 */

#include <assert.h>
#include "noll.h"
#include "noll_fun.h"
#include "noll_internal.h"

/* ============================================================ */
/* I.1 Memory */
/* ============================================================ */

/* -------------------------------------------------------------*/
/* I.1.1 Noll values */
/* -------------------------------------------------------------*/

noll_val_t *
noll_copy(sh_manager_t * man, noll_val_t * a) {
	if (a == NULL)
		return NULL;

	noll_val_t* r = g_new(noll_val_t,1);
	g_ptr_array_copy(r->set, a->set, noll_graph_copy, noll_graph_free);
	return r;
}

void noll_free(sh_manager_t * man, noll_val_t * a) {
	if (a == NULL)
		return;
	g_ptr_array_free(a->set, TRUE); // noll_graph_free applied also
	return;
}

void noll_free_set(noll_val_t* a) {
	if (a == NULL)
		return;
	g_ptr_array_free(a->set, FALSE); // noll_graph_free is not applied!
	return;
}

size_t noll_size(sh_manager_t * man, noll_val_t * a) {
	assert (NULL != man);

	if (a == NULL)
		return 0;
	return (a->set == NULL) ? 0 : a->set->len;
}

/* -------------------------------------------------------------*/
/* I.1.2 Noll graph values */
/* -------------------------------------------------------------*/

noll_graph_t* noll_graph_emp(sh_manager_t* man, size_t fid) {

	noll_graph_t* r = g_new(noll_graph_t,1);
	/* variables frame */
	r->frame = fid;
	sh_frame_t* frame = sh_stack_get_frame(sh_manager_get_varenv(man), fid);
	size_t frame_size = (frame->vars == NULL) ? 0 : frame->vars->len;
	r->lvars = NULL;
	r->svars = NULL;

	/* vertices */
	r->vrt_size = 1; /* only NULL */
	r->vrt_type = g_array_new(FALSE, FALSE, sizeof(size_t));
	size_t v = SH_VAR_NULL;
	g_array_append_val(r->vrt_type,v); // type of NULL
	/* fill labeling of vertices by program variables */
	r->var2vrt
			= g_array_sized_new(FALSE, FALSE, sizeof(size_t), frame_size + 1);
	v = 0;
	g_array_append_val(r->var2vrt,v);
	for (guint i = 0; i < frame_size; i++)
		g_array_append_val(r->var2vrt,v);

	r->isprecise = TRUE;

	/* edges */
	r->e_mat = g_ptr_array_new();
	g_ptr_array_add(r->e_mat, NULL); // entry for NULL
	r->e_rmat = g_ptr_array_new();
	g_ptr_array_add(r->e_rmat, NULL); // reverse entry for NULL
	r->edges = NULL;

	/* separation */
	r->ssep = NULL;
	r->wsep = NULL;

	return r;
}

noll_graph_t* noll_graph_true(sh_manager_t* man, size_t fid) {

	noll_graph_t* r = g_new(noll_graph_t,1);
	/* variables frame */
	r->frame = fid;
	r->lvars = NULL;
	r->svars = NULL;

	/* vertices */
	r->vrt_size = 0; /* nothing */
	r->vrt_type = NULL;
	r->var2vrt = NULL;
	r->isprecise = FALSE;

	/* edges */
	r->e_mat = NULL;
	r->e_rmat = NULL;
	r->edges = NULL;

	/* separation */
	r->ssep = NULL;
	r->wsep = NULL;

	return r;
}

noll_graph_t*
noll_graph_copy(noll_graph_t* a) {
	if (a == NULL)
		return NULL;

	noll_graph_t* r = g_new(noll_graph_t,1);

	r->frame = a->frame;
	g_array_copy(r->lvars,a->lvars,size_t);
	g_ptr_array_copy(r->svars,a->svars,noll_svar_copy,noll_svar_free);

	r->vrt_size = a->vrt_size;
	g_array_copy(r->vrt_type,a->vrt_type,size_t);
	g_array_copy(r->var2vrt,a->var2vrt,size_t);
	r->isprecise = a->isprecise;

	g_ptr_array_copy(r->edges,a->edges,noll_edge_copy,noll_edge_free);
	/* special case, fill with arrays of ids */
	r->e_mat = NULL;
	if (a->e_mat != NULL) {
		r->e_mat = g_ptr_array_new(); /* no free */
		for (guint i = 0; i < a->e_mat->len; i++) {
			GArray* a_i = g_ptr_array_index(a->e_mat,i);
			/* copy this array */
			GArray* r_i = NULL;
			g_array_copy(r_i,a_i,size_t);
			g_ptr_array_add(r->e_mat, r_i);
		}
	}
	/* special case, fill with arrays of ids */
	r->e_rmat = NULL;
	if (a->e_rmat != NULL) {
		r->e_rmat = g_ptr_array_new(); /* no free */
		for (guint i = 0; i < a->e_rmat->len; i++) {
			GArray* a_i = g_ptr_array_index(a->e_rmat,i);
			/* copy this array */
			GArray* r_i = NULL;
			g_array_copy(r_i,a_i,size_t);
			g_ptr_array_add(r->e_rmat, r_i);
		}
	}

	r->ssep = NULL;
	if (a->ssep != NULL) {
		r->ssep = g_ptr_array_new();
		for (guint i = 0; i < a->ssep->len; i++) {
			GArray* a_i = g_ptr_array_index(a->ssep,i);
			GArray* r_i = NULL;
			g_array_copy(r_i,a_i,size_t);
			g_ptr_array_add(r->ssep, r_i);
		}
	}

	g_ptr_array_copy(r->wsep, a->wsep, noll_share_copy, noll_share_free);

	return r;
}

void noll_graph_free(noll_graph_t* a) {
	if (a == NULL)
		return;
	g_array_free(a->lvars, TRUE);
	a->lvars = NULL;
	g_ptr_array_free(a->svars, TRUE);
	a->svars = NULL;
	g_array_free(a->vrt_type, TRUE);
	a->vrt_type = NULL;
	g_array_free(a->var2vrt, TRUE);
	a->var2vrt = NULL;

	for (guint i = 0; i < a->e_mat->len; i++) {
		GArray* a_i = g_ptr_array_index(a->e_mat,i);
		g_array_free(a_i, TRUE);
		g_ptr_array_index(a->e_mat,i) = NULL;
		a_i = g_ptr_array_index(a->e_rmat,i);
		g_array_free(a_i, TRUE);
		g_ptr_array_index(a->e_rmat,i) = NULL;
	}
	g_ptr_array_free(a->e_mat, TRUE);
	g_ptr_array_free(a->e_rmat, TRUE);
	g_ptr_array_free(a->edges, TRUE);
	a->edges = NULL;
}

/* -------------------------------------------------------------*/
/* I.1.3 Sharing constraints */
/* -------------------------------------------------------------*/

noll_share_t*
noll_share_new(noll_op_e op) {
	noll_share_t* r = g_new(noll_share_t,1);
	r->op = op;
	r->op = SH_VAR_NULL;
	r->vars_t1 = r->lvars_t1 = r->svars_t1 = NULL;
	r->vars_t2 = r->lvars_t2 = r->svars_t2 = NULL;
	return r;
}

noll_share_t*
noll_share_copy(noll_share_t* a) {
	if (a == NULL)
		return NULL;
	noll_share_t* r = noll_share_new(a->op);
	r->op = a->op;
	g_array_copy(r->vars_t1, a->vars_t1, size_t);
	g_array_copy(r->lvars_t1, a->lvars_t1, size_t);
	g_array_copy(r->svars_t1, a->svars_t1, size_t);
	g_array_copy(r->vars_t2, a->vars_t2, size_t);
	g_array_copy(r->lvars_t2, a->lvars_t2, size_t);
	g_array_copy(r->svars_t2, a->svars_t2, size_t);
	return r;
}

void noll_share_free(noll_share_t* a) {
	if (a == NULL)
		return;
	g_array_free(a->vars_t1, TRUE);
	g_array_free(a->lvars_t1, TRUE);
	g_array_free(a->svars_t1, TRUE);
	g_array_free(a->vars_t2, TRUE);
	g_array_free(a->lvars_t2, TRUE);
	g_free(a);
	return;
}

/*
 * Lexicographic comparison over (t1, t2, op).
 * Notice that GArray are all sorted!
 * Result a-b.
 */
int noll_share_cmp(noll_share_t* a, noll_share_t* b) {
	if (a == NULL && b != NULL)
		return -1;
	if (a == NULL && b == NULL)
		return 0;
	if (a != NULL && b == NULL)
		return 1;
	// a != NULL && b != NULL
	/* compare types */
	if (a->tid < b->tid)
		return -1;
	else if (a->tid > b->tid)
		return 1;
	/* Compare vars_t1 */
	int r;
	g_array_cmp(r, a->vars_t1, b->vars_t1, size_t);
	if (r < 0 || r > 0)
		return r;
	/* Compare lvars_t1 */
	g_array_cmp(r, a->lvars_t1, b->lvars_t1, size_t);
	if (r < 0 || r > 0)
		return r;
	/* Compare svars_t1 */
	g_array_cmp(r, a->svars_t1, b->svars_t1, size_t);
	if (r < 0 || r > 0)
		return r;
	/* Compare vars_t2 */
	g_array_cmp(r, a->vars_t2, b->vars_t2, size_t);
	if (r < 0 || r > 0)
		return r;
	/* Compare lvars_t2 */
	g_array_cmp(r,a->lvars_t2, b->lvars_t2, size_t);
	if (r < 0 || r > 0)
		return r;
	/* Compare svars_t1 */
	g_array_cmp(r,a->svars_t2, b->svars_t2, size_t);
	if (r < 0 || r > 0)
		return r;
	if (a->op < b->op)
		return -1;
	if (a->op > b->op)
		return 1;
	return 0;
}

/* -------------------------------------------------------------*/
/* I.1.4 Set of locations values */
/* -------------------------------------------------------------*/

noll_svar_t*
noll_svar_new(char* name, GArray* tids, size_t e) {
	noll_svar_t* r = g_new(noll_svar_t,1);
	r->name = (name == NULL) ? NULL : g_strdup(name);
	g_array_copy(r->tids, tids, size_t);
	r->bto_eid = e;
	return r;
}

void noll_svar_free(noll_svar_t* a) {
	if (a == NULL)
		return;
	if (a->name != NULL)
		g_free(a->name);
	if (a->tids != NULL)
		g_array_free(a->tids, TRUE);
	g_free(a);
}

noll_svar_t*
noll_svar_copy(noll_svar_t* a) {
	if (a == NULL)
		return NULL;
	return noll_svar_new(a->name, a->tids, a->bto_eid);
}

/* -------------------------------------------------------------*/
/* I.1.5 Edges */
/* -------------------------------------------------------------*/

noll_edge_t*
noll_edge_new(size_t id, size_t src, size_t dst, noll_edge_e kind) {
	noll_edge_t* r = g_new(noll_edge_t,1);
	r->eid = id;
	r->kind = kind;
	r->vrt_src = src;
	r->vrt_dst = dst;
	return r;
}

void noll_edge_free(noll_edge_t* a) {
	if (a == NULL)
		return;
	if (a->kind == NOLL_EDGE_LS)
		g_array_free(a->info.pred.ngb, TRUE);
	g_free(a);
}

noll_edge_t*
noll_edge_copy(noll_edge_t* a) {
	if (a == NULL)
		return NULL;
	noll_edge_t* r = g_new(noll_edge_t,1);
	r->eid = a->eid;
	r->kind = a->kind;
	r->vrt_src = a->vrt_src;
	r->vrt_dst = a->vrt_dst;
	if (a->kind == NOLL_EDGE_PTO)
		r->info.fid = a->info.fid;
	else if (a->kind == NOLL_EDGE_LS) {
		r->info.pred.isrev = a->info.pred.isrev;
		r->info.pred.pid = a->info.pred.pid;
		r->info.pred.ngb = g_array_new(FALSE,FALSE,sizeof(size_t));
		for (guint i = 0; i < a->info.pred.ngb->len; i++) {
			size_t v = g_array_index(a->info.pred.ngb,size_t,i);
			g_array_append_val(r->info.pred.ngb, v);
		}
	}
	return r;
}

/*
 * Lexicographical comparison in order (src,kind,dst,fid|pid)
 */
int noll_edge_cmp(noll_edge_t* a, noll_edge_t* b) {
	if (a == NULL && b != NULL)
		return -1;
	if (a == NULL && b == NULL)
		return 0;
	if (a != NULL && b == NULL)
		return 1;
	// a != NULL && b != NULL
	/* compare src */
	int r = a->vrt_src - b->vrt_src;
	if (r != 0)
		return (r < 0) ? -1 : 1;
	/* compare kind */
	r = a->kind - b->kind;
	if (r != 0)
		return (r < 0) ? -1 : 1;
	/* compare dst */
	r = a->vrt_dst - b->vrt_dst;
	if (r != 0)
		return (r < 0) ? -1 : 1;
	if (a->kind != NOLL_EDGE_PTO && a->kind != NOLL_EDGE_LS)
		return 0;
	/* compare for pto */
	if (a->kind == NOLL_EDGE_PTO)
		return (int) (a->info.fid - b->info.fid);
	/* compare for pred */
	r = a->info.pred.pid - b->info.pred.pid;
	if (r != 0)
		return (r < 0) ? -1 : 1;
	/* compare args */
	GArray* va = a->info.pred.ngb;
	GArray* vb = b->info.pred.ngb;
	g_array_cmp(r,va,vb,size_t);
	return r;
}

/* -------------------------------------------------------------*/
/* I.1.5 Homomorphism */
/* -------------------------------------------------------------*/

noll_hom_t*
noll_hom_alloc(size_t vrt, size_t edg) {
	if (vrt >= SH_VAR_NULL && edg >= NOLL_EDGE_ID_UNKNOWN)
		return NULL;
	noll_hom_t* h = g_new(noll_hom_t,1);
	h->vrt = g_array_sized_new(FALSE, FALSE, vrt, sizeof(size_t));
	size_t vrt2 = SH_VAR_NULL; // unknown, by convention
	for (guint i = 0; i < vrt; i++)
		g_array_append_val(h->vrt,vrt2);
	if (edg > 0) {
		h->used = g_ptr_array_sized_new(edg);
		for (guint i = 0; i < edg; i++)
			g_ptr_array_add(h->used, NULL);
	} else
		h->used = NULL;
	return h;
}

void noll_hom_free(noll_hom_t* h) {
	if (h == NULL)
		return;
	if (h->vrt != NULL)
		g_array_free(h->vrt, TRUE);
	h->vrt = NULL;
	if (h->used != NULL)
		g_ptr_array_free(h->used, TRUE);
	h->used = NULL;
	g_free(h);
}

GQueue*
noll_hom_fill_pvars(sh_manager_t* man, noll_hom_t* h, noll_graph_t* a1,
		noll_graph_t* a2) {
	if (a1 == NULL || a2 == NULL)
		return NULL;

	/* check is done on the same frame */
	assert (a1->frame == a2->frame);
	sh_frame_t* frame = sh_stack_get_frame(sh_manager_get_varenv(man),
			a1->frame);

	/* check that h is already allocated */
	assert (h->vrt != NULL);
	assert (h->vrt->len == a1->vrt_size);

	GQueue* vrt_added = g_queue_new();

	/* start by mapping vertices using program variables */
	g_array_index(h->vrt, size_t, 0) = 0; // NULL
	for (guint vi = 0; vi < frame->vars->len; vi++) {
		size_t n1_i = g_array_index(a1->var2vrt,size_t,vi);
		size_t n2_i = g_array_index(a2->var2vrt,size_t,vi);
		if (n1_i == 0 && n2_i != 0) {
			// null is mapped to a not null vertex
			g_queue_free(vrt_added);
			return NULL;
		}
		g_array_index(h->vrt, size_t, n1_i) = n2_i;
		g_queue_push_tail(vrt_added, GSIZE_TO_POINTER(n1_i));
	}
	return vrt_added;
}

gboolean noll_hom_fill_pto(noll_hom_t* h, GQueue* vrt_added, noll_graph_t* a1,
		noll_graph_t* a2) {

	GArray* vrt_visited = g_array_sized_new(FALSE, FALSE, sizeof(size_t),
			a1->vrt_size);
	/* set flags for already explored vertices */
	size_t flag = 0;
	for (guint i = 0; i < a1->vrt_size; i++)
		g_array_append_val(vrt_visited,flag); // 1 if vertex already explored

	gboolean error = FALSE;
	fill_pto_err: if (error == TRUE) {
		if (vrt_visited != NULL)
			g_array_free(vrt_visited, TRUE);
		return FALSE;
	}

	while (g_queue_is_empty(vrt_added) == FALSE) {
		size_t n1_i = GPOINTER_TO_UINT(g_queue_pop_head(vrt_added));
		if (g_array_index(vrt_visited,size_t,n1_i) == 1) continue;
		size_t n2_i = g_array_index(h->vrt,size_t,n1_i); // h(n1_i)
		// look at pto edges from n1_i
		GArray* edg_n1_i = g_ptr_array_index(a1->e_mat,n1_i);
		// mark n1_i as visited
		g_array_index(vrt_visited,size_t,n1_i) = 1;
		if (edg_n1_i == NULL)
			continue;
		for (guint ei = 0; ei < edg_n1_i->len; ei++) {
			size_t eid_1i = g_array_index(edg_n1_i,size_t,ei);
			noll_edge_t* e1i = g_ptr_array_index(a1->edges,eid_1i);
			size_t dst_1i = e1i->vrt_dst;
			if (e1i->kind == NOLL_EDGE_PTO) {
				size_t eid_2i = noll_graph_get_edge_from(a2, n2_i,
						NOLL_EDGE_PTO, e1i->info.fid);
				if (eid_2i == NOLL_EDGE_ID_UNKNOWN)
					goto fill_pto_err;

				// else, map destinations by h
				noll_edge_t* e2i = g_ptr_array_index(a2->edges,eid_2i);
				assert(e2i != NULL);
				size_t dst_2i = e2i->vrt_dst;
				// map destination node
				g_array_index(h->vrt,size_t,dst_1i) = dst_2i;
				// add destination to queue if not already visited
				if (g_array_index(vrt_visited,size_t,dst_1i) == 0)
					g_queue_push_tail(vrt_added, GUINT_TO_POINTER(dst_1i));
				// add the edges to the homomorphism
				GPtrArray* used_e1 = g_ptr_array_index(h->used,eid_1i);
				if (used_e1 == NULL) {
					used_e1 = g_ptr_array_sized_new(a2->edges->len);
					for (guint i = 0; i < a2->edges->len; i++)
						g_ptr_array_add(used_e1, NULL);
					g_ptr_array_index(h->used,eid_1i) = used_e1;
				}
				noll_edgflds_t* used_e1_e2 = g_ptr_array_index(used_e1,eid_2i);
				if (used_e1_e2 == NULL) {
					used_e1_e2 = g_new(noll_edgflds_t,1);
					used_e1_e2->eid = eid_2i;
					used_e1_e2->flds
							= g_array_new(FALSE, FALSE, sizeof(size_t));
				}
				g_array_append_val(used_e1_e2->flds, e1i->info.fid);
				// TODO: sort array of fields
				g_ptr_array_index(used_e1,eid_2i) = used_e1_e2;
			}
		}
	}

	/* see if this does not imply inconsistency in distinct edges */
	if (noll_hom_is_consistent_diff(h, a1, a2) == FALSE)
	{
		error = TRUE;
		goto fill_pto_err;
	}
	/* see if strong separation constraints are filled */
	if (noll_hom_is_consistent_ssep(h, a1, a2) == FALSE)
	{
		error = TRUE;
		goto fill_pto_err;
	}

	return TRUE;
}

GPtrArray*
noll_hom_fill_pred(noll_hom_t* h, noll_graph_t* a1, noll_graph_t* a2) {

	GPtrArray* sigma = NULL; // substitution of svars
	/* which is an array of noll_share_t* where only term t1 is filled */
	if (a1->svars != NULL && a1->svars->len > 0) {
		sigma = g_ptr_array_sized_new(a1->svars->len);
		for (guint i = 0; i < a1->svars->len; i++)
			g_ptr_array_add(sigma, NULL);
	}
	return sigma; // TODO
}

/* ============================================================ */
/* II.1 Basic constructors */
/* ============================================================ */

noll_val_t *
noll_bottom(sh_manager_t * man, size_t fid) {
	return NULL;
}

noll_val_t *
noll_top(sh_manager_t * man, size_t fid) {
	/* different from NULL, it contains the empty set */
	noll_val_t* r = g_new(noll_val_t,1);
	r->frame = fid;
	r->set = g_ptr_array_new();
	return r;
}

noll_val_t *
noll_empty(sh_manager_t * man, size_t fid) {
	/* different from NULL, it contains the empty graph */
	noll_val_t* r = g_new(noll_val_t,1);
	r->frame = fid;
	r->set = g_ptr_array_new_with_free_func((GDestroyNotify) noll_graph_free);
	g_ptr_array_add(r->set, noll_graph_emp(man, fid));
	return r;
}

noll_val_t *
noll_of_formula(sh_manager_t* man, size_t fid, void* f) {
	return NULL; // TODO
}

/* ============================================================ */
/* II.2 Changing graph */
/* ============================================================ */

noll_graph_t*
noll_graph_add_vertex(sh_manager_t* man, gboolean destructive, noll_graph_t* a,
		size_t tid) {
	noll_graph_t* r = (destructive == TRUE) ? a : noll_graph_copy(a);
	if (r != NULL) {
		size_t vrt = r->vrt_size++;
		assert (r->vrt_type->len == vrt);
		g_array_append_val (r->vrt_type,tid);
		// add also the entry in the adjacency matrix
		assert (r->e_mat->len == vrt);
		g_ptr_array_add(r->e_mat, NULL);
		assert (r->e_rmat->len == vrt);
		g_ptr_array_add(r->e_rmat, NULL);
	}
	return r;
}

noll_graph_t*
noll_graph_set_label(sh_manager_t* man, gboolean destructive, noll_graph_t* a,
		size_t var, size_t vrt) {
	noll_graph_t* r = (destructive == TRUE) ? a : noll_graph_copy(a);
	if (r != NULL) {
		size_t ivar = var + 1;
		size_t old_vrt = g_array_index(r->var2vrt,size_t,ivar);
		g_array_index(r->var2vrt,size_t,ivar) = vrt;
		noll_graph_t* r_clean = noll_graph_remove_junk(man, r, old_vrt);
		if (r_clean != r) {
			// TODO: better signal junk
			fprintf(stderr, "(frame-%d): Junk detected from old var-%d!\n",
					a->frame, var);
		}
	}
	return r;
}

noll_graph_t* noll_graph_set_edge_dst(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t eid, size_t vrt) {
	noll_graph_t* r = (destructive == TRUE) ? a : noll_graph_copy(a);
	if (r != NULL) {
		assert (eid < r->edges->len);
		noll_edge_t* e = g_ptr_array_index(r->edges,eid);
		size_t old_vrt = e->vrt_dst;
		e->vrt_dst = vrt; // in place change for edges

		if (old_vrt != 0) {
			// remove the edge from r->e_rmat[old_vrt]
			GArray* from_edges = g_ptr_array_index(r->e_rmat,old_vrt);
			assert (from_edges != NULL);
			guint i = 0;
			for (; i < from_edges->len; i++)
				if (g_array_index(from_edges,size_t,i) == eid)
					break;
			assert (i < from_edges->len);
			g_ptr_array_index(r->e_rmat,old_vrt) = g_array_remove_index(
					from_edges, i);
		}

		// add the edge to  r->e_rmat[vrt]
		GArray* from_edges = g_ptr_array_index(r->e_rmat,vrt);
		if (from_edges == NULL)
			from_edges = g_array_new(FALSE, FALSE, sizeof(size_t));
		g_array_append_val(from_edges,eid);
		g_ptr_array_index(r->e_rmat,vrt) = from_edges;

		if (old_vrt != 0) {
			// detect junk
			noll_graph_t* r_clean = noll_graph_remove_junk(man, r, old_vrt);
			if (r_clean != r) {
				// TODO: better signal junk
				fprintf(stderr, "(frame-%d): Junk detected from old vrt-%d!\n",
						a->frame, old_vrt);
			}
		}
	}
	return r;
}

noll_graph_t*
noll_graph_add_edge_diff(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t src, size_t dst) {
	noll_graph_t* r = (destructive == TRUE) ? a : noll_graph_copy(a);
	assert (r != NULL);
	assert (src < r->vrt_size);
	assert (dst < r->vrt_size);
	/* add the edge in the general table */
	size_t eid = 0;
	if (r->edges == NULL)
		r->edges = g_ptr_array_new();
	eid = r->edges->len;
	noll_edge_t* e = noll_edge_new(eid, src, dst, NOLL_EDGE_DIFF);
	g_ptr_array_add(r->edges, e);

	/* add the edge in the adjacency matrix */
	GArray* edges_src = g_ptr_array_index(r->e_mat,src);
	if (edges_src == NULL)
		edges_src = g_array_new(FALSE, FALSE, sizeof(size_t));
	g_array_append_val(edges_src,eid);
	g_ptr_array_index(r->e_mat,src) = edges_src;

	GArray* edges_dst = g_ptr_array_index(r->e_mat,dst);
	if (edges_dst == NULL)
		edges_dst = g_array_new(FALSE, FALSE, sizeof(size_t));
	g_array_append_val(edges_dst,eid);
	g_ptr_array_index(r->e_mat,dst) = edges_dst;
	// no need for the reverse matrix

	return r;
}

noll_graph_t* noll_graph_add_edge_pto(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t src, size_t fid, size_t dst, size_t* eid) {

	noll_graph_t* r = (destructive == TRUE) ? a : noll_graph_copy(a);
	assert (r != NULL);
	assert (src < r->vrt_size);
	assert (dst < r->vrt_size);
	/* add the edge in the general table */
	if (r->edges == NULL)
		r->edges = g_ptr_array_new();
	*eid = r->edges->len;
	noll_edge_t* e = noll_edge_new(*eid, src, dst, NOLL_EDGE_PTO);
	e->info.fid = fid;
	g_ptr_array_add(r->edges, e);

	/* add the edge in the adjacency matrix */
	GArray* edges_src = g_ptr_array_index(r->e_mat,src);
	if (edges_src == NULL)
		edges_src = g_array_new(FALSE, FALSE, sizeof(size_t));
	g_array_append_val(edges_src,*eid);
	g_ptr_array_index(r->e_mat,src) = edges_src;

	// add to the reverse matrix if dst != NULL
	if (dst != 0) {
		GArray* edges_dst = g_ptr_array_index(r->e_rmat,dst);
		if (edges_dst == NULL)
			edges_dst = g_array_new(FALSE, FALSE, sizeof(size_t));
		g_array_append_val(edges_dst, *eid);
		g_ptr_array_index(r->e_rmat,dst) = edges_dst;
	}

	return r;
}

noll_graph_t*
noll_graph_remove_junk(sh_manager_t* man, noll_graph_t* a, size_t vrt) {
	if (a == NULL)
		return a;
	/* check if there is some var labeling vrt */
	gboolean found = FALSE;
	for (guint i = 0; i < a->var2vrt->len && (found == FALSE); i++) {
		size_t vi = g_array_index(a->var2vrt,size_t,i);
		if (vi == vrt)
			found = TRUE;
	}
	if (found == TRUE)
		return a;
	/* check if there is some edge going to vrt */
	GArray* to_vrt = g_ptr_array_index(a->e_rmat,vrt);
	if (to_vrt != NULL)
		return a;
	/* TODO: vrt is garbage, remove all edges from it */
	return NULL;
}

noll_graph_t* noll_graph_collapse_edge(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t eid, size_t* vrt_f, size_t* eid_f) {
	if (destructive == TRUE)
		return a; // TODO
	else
		return noll_graph_copy(a); // TODO
}

noll_graph_t* noll_graph_unfold_start(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t eid) {
	if (destructive == TRUE)
		return a; // TODO
	else
		return noll_graph_copy(a); // TODO
}
noll_graph_t* noll_graph_unfold_end(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t eid) {
	if (destructive == TRUE)
		return a; // TODO
	else
		return noll_graph_copy(a); // TODO
}

noll_graph_t* noll_graph_unfold_middle(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t eid, size_t tid) {
	if (destructive == TRUE)
		return a; // TODO
	else
		return noll_graph_copy(a); // TODO
}
