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
 * noll_representation.c: Implementation of functions related to the
 *                        data structure representation
 */

#include <assert.h>
#include "noll.h"
#include "noll_fun.h"
#include "noll_internal.h"

/* ============================================================ */
/* II.2 Getters */
/* ============================================================ */

size_t noll_dimension(sh_manager_t * man, noll_val_t * a) {
	if (a == NULL)
		return 0;
	return a->frame;
}

size_t noll_graph_get_edge_from(noll_graph_t* a, size_t vrt, noll_edge_e kind,
		size_t id) {
	if ((a == NULL) || (a->vrt_size <= vrt) || (a->e_mat == NULL)
			|| (a->e_mat->len <= vrt))
		return NOLL_EDGE_ID_UNKNOWN;

	GArray* edg_vrt = g_ptr_array_index(a->e_mat, vrt);
	for (guint i = 0; i < edg_vrt->len; i++) {
		size_t eid = g_array_index(edg_vrt,size_t,i);
		noll_edge_t* e = g_ptr_array_index(a->edges,eid);
		if (e->kind != kind)
			continue;
		if ((kind == NOLL_EDGE_PTO && e->info.fid == id) || (kind
				== NOLL_EDGE_LS && e->info.pred.pid == id))
			return eid;
	}
	return NOLL_EDGE_ID_UNKNOWN;
}

size_t noll_graph_get_label(noll_graph_t* a, size_t var) {
	size_t var_entry = var + 1; // for NULL
	return g_array_index(a->var2vrt,size_t,var_entry);
}

/* ============================================================ */
/* II.3 Tests */
/* ============================================================ */

/* ------------------------------------------------------------ */
/* II.3.1 Values */
/* ------------------------------------------------------------ */

gboolean noll_is_bottom(sh_manager_t * man, noll_val_t * a) {
	if (a == NULL)
		return TRUE;
	/* check consistency */
	return (noll_is_sat(man, a) == TRUE) ? FALSE : TRUE;
}

/*
 * Over-approximate the exact test.
 */
gboolean noll_is_top(sh_manager_t * man, noll_val_t * a) {
	if (a == NULL)
		return FALSE;
	if (a->set == NULL || a->set->len == 0)
		return TRUE;
	for (guint i = 0; i < a->set->len; i++) {
		noll_graph_t* gi = g_ptr_array_index(a->set,i);
		if (noll_graph_is_true(man, gi) == FALSE)
			return FALSE;
	}
	return TRUE;
}

gboolean noll_is_le(sh_manager_t * man, noll_val_t * a1, noll_val_t * a2) {
	if (a1 == NULL)
		return TRUE;
	if (a2 == NULL)
		return FALSE;
	/* a1 != NULL, a2 != NULL */

	/* compare equal frames */
	assert(a1->frame == a2->frame);

	/* powerset comparison */
	for (guint i = 0; i < a1->set->len; i++) {
		noll_graph_t* g1 = g_ptr_array_index(a1->set,i);
		gboolean found = FALSE;
		for (guint j = 0; j < a2->set->len && (found == FALSE); j++) {
			noll_graph_t* g2 = g_ptr_array_index(a2->set,j);
			found = noll_graph_is_le(man, g1, g2);
		}
		if (found == FALSE)
			return FALSE;
	}
	return TRUE;
}

gboolean noll_is_eq(sh_manager_t * man, noll_val_t * a1, noll_val_t * a2) {
	if (a1 == NULL)
		return TRUE;
	if (a2 == NULL)
		return FALSE;
	/* a1 != NULL, a2 != NULL */

	/* compare equal frames */
	assert(a1->frame == a2->frame);

	/* powerset comparison */
	for (guint i = 0; i < a1->set->len; i++) {
		noll_graph_t* g1 = g_ptr_array_index(a1->set,i);
		gboolean found = FALSE;
		for (guint j = 0; j < a2->set->len && (found == FALSE); j++) {
			noll_graph_t* g2 = g_ptr_array_index(a2->set,j);
			found = noll_graph_is_eq(man, g1, g2);
		}
		if (found == FALSE)
			return FALSE;
	}
	return TRUE;
}

gboolean noll_sat_constraint(sh_manager_t * man, noll_val_t * a, size_t x,
		size_t y, gboolean iseq) {
	return FALSE; // TODO
}

gboolean noll_is_dimension_unconstrained(sh_manager_t * man, noll_val_t * a,
		size_t x) {
	return FALSE; // TODO
}

gboolean noll_is_dimension_shared(sh_manager_t * man, noll_val_t * a, size_t x) {
	return FALSE; // TODO
}

gboolean noll_is_dimension_nil(sh_manager_t * man, noll_val_t * a, size_t x) {
	return FALSE; // TODO
}

gboolean noll_is_sat(sh_manager_t * man, noll_val_t * a) {
	if (a == NULL)
		return FALSE;
	if (a->set == NULL || a->set->len == 0)
		return TRUE;
	for (guint i = 0; i < a->set->len; i++) {
		noll_graph_t* gi = g_ptr_array_index(a->set,i);
		if (noll_graph_is_sat(man, gi) == TRUE)
			return TRUE;
	}
	return FALSE;
}

/* ------------------------------------------------------------ */
/* II.3.2 Graphs */
/* ------------------------------------------------------------ */

gboolean noll_graph_is_sat(sh_manager_t * man, noll_graph_t * a) {
	if (a == NULL)
		return FALSE;

	/* simple test, not difference edges between nodes */
	for (guint i = 0; i < a->edges->len; i++) {
		noll_edge_t* ei = g_ptr_array_index(a->edges,i);
		if (ei->kind == NOLL_EDGE_DIFF && ei->vrt_src == ei->vrt_dst)
			return FALSE;
	}
	// TODO: complete the sat test
	return TRUE;
}

gboolean noll_graph_is_emp(sh_manager_t * man, noll_graph_t * a) {
	if (a == NULL)
		return FALSE;
	return (a->vrt_size == 1 && a->isprecise == TRUE);
}

gboolean noll_graph_is_true(sh_manager_t * man, noll_graph_t * a) {
	if (a == NULL)
		return FALSE;
	return (a->vrt_size == 0 && a->isprecise == FALSE);
}

gboolean noll_graph_is_top(sh_manager_t * man, noll_graph_t * a) {
	if (a == NULL)
		return FALSE;
	if (a->isprecise == TRUE || a->vrt_size >= 2)
		return FALSE;
	return TRUE; // TOD: approx test, do better?
}

gboolean noll_graph_is_le(sh_manager_t * man, noll_graph_t * a1,
		noll_graph_t * a2) {
	/* call procedure for homomorphism computation */
	noll_hom_t* h = noll_graph_is_homomorphic(man, a2, a1);
	if (h == NULL)
		return FALSE;
	noll_hom_free(h);
	return TRUE;
}

gboolean noll_graph_is_eq(sh_manager_t * man, noll_graph_t * a1,
		noll_graph_t * a2) {
	if (a1 == NULL && a2 == NULL)
		return TRUE;
	if (a1 == NULL || a2 == NULL)
		return FALSE;
	/* a1 != NULL, a2 != NULL */

	/* compare equal frames */
	assert(a1->frame == a2->frame);

	/* shall have same number of edges */
	if (a1->vrt_size != a2->vrt_size)
		return FALSE;

	/* call procedure for homomorphism computation */
	noll_hom_t* h = noll_graph_is_homomorphic(man, a1, a2);
	if (h == NULL)
		return FALSE;

	/* check that it is an isomorphism, i.e., injective */
	gboolean isom = noll_hom_is_isom(h);
	noll_hom_free(h);
	return isom;
}

/* ------------------------------------------------------------ */
/* II.3.3 Homomorphism */
/* ------------------------------------------------------------ */

noll_hom_t*
noll_graph_is_homomorphic(sh_manager_t * man, noll_graph_t * a1,
		noll_graph_t * a2) {
	if (a1 == NULL || a2 == NULL)
		return NULL;

	/* check is done on the same frame */
	assert (a1->frame == a2->frame);
	size_t vrt1_size = a1->vrt_size;
	size_t edg1_size = (a1->edges == NULL) ? 0 : a1->edges->len;
	noll_hom_t* h = noll_hom_alloc(vrt1_size, edg1_size);

	/*
	 * all declarations here are added to return in case of error
	 */
	gboolean error = FALSE;
	GQueue* vrt_added = NULL; // queue of vertices to be explored
	GPtrArray* sigma = NULL; // substitution of svars

	herr: if (error == TRUE) {
		if (h != NULL)
			noll_hom_free(h);
		if (vrt_added != NULL)
			g_queue_free(vrt_added);
		if (sigma != NULL)
			g_ptr_array_free(sigma, TRUE);
		return NULL;
	}

	/* start by mapping vertices using program variables */
	vrt_added = noll_hom_fill_pvars(man, h, a1, a2);
	if (vrt_added == NULL || g_queue_is_empty(vrt_added) == TRUE) {
		error = TRUE;
		goto herr;
	}

	/*
	 * fill then looking at PTO edges starting from labeled nodes
	 */
	if (noll_hom_fill_pto(h, vrt_added, a1, a2) == FALSE) {
		error = TRUE;
		goto herr;
	}
	// all mappings related with PTO edges have been explored

	/*
	 * fill homomorphism by looking at predicate edges
	 */
	sigma = noll_hom_fill_pred(h, a1, a2);
	if (sigma == NULL) {
		error = TRUE;
		goto herr;
	}
	/* sigma gives the substitution of svars in a1 by terms in a2 */

	/* TODO: check that all nodes in a1 are mapped */
	/* TODO: check that sharing constraints are implied */

	return h; // TODO
}

gboolean noll_hom_is_consistent_diff(noll_hom_t* h, noll_graph_t* a1,
		noll_graph_t* a2) {
	if (h == NULL || a1 == NULL || a2 == NULL)
		return FALSE;
	/* go through all edges in a1 */
	for (guint e1i = 0; e1i < a1->edges->len; e1i++) {
		noll_edge_t* e1 = g_ptr_array_index(a1->edges,e1i);
		if (e1->kind != NOLL_EDGE_DIFF)
			continue;
		size_t src_1 = e1->vrt_src;
		size_t dst_1 = e1->vrt_dst;
		size_t src_2 = g_array_index(h->vrt,size_t,e1->vrt_src);
		size_t dst_2 = g_array_index(h->vrt,size_t,e1->vrt_dst);
		if (src_2 == SH_VAR_NULL || dst_2 == SH_VAR_NULL)
			continue; /* edge not yet mapped */
		/* at this point, the edge is mapped */
		/* search difference edges from src_2 */
		GArray* edg_2 = g_ptr_array_index(a2->e_mat,src_2);
		noll_edge_t* h_e1 = NULL;
		for (guint j = 0; j < edg_2->len; j++) {
			size_t e2j = g_array_index(edg_2,size_t,j);
			noll_edge_t *e2 = g_ptr_array_index(a2->edges,e2j);
			if (e2->kind != NOLL_EDGE_DIFF)
				continue;
			if (e2->vrt_src == src_2 && e2->vrt_dst == dst_2)
				h_e1 = e2;
		}
		if (h_e1 == NULL)
			return FALSE;
	}
	return TRUE;
}

gboolean noll_hom_is_consistent_ssep(noll_hom_t* h, noll_graph_t* a1,
		noll_graph_t* a2) {

	if (h == NULL || a1 == NULL || a2 == NULL)
		return FALSE;

	if (a1->ssep == NULL)
		return TRUE;

	// VERY BAD ALGORITHM in O(|edges|^4)
	// go through separation constraints in a1
	for (guint ei = 0; ei < a1->ssep->len; ei++) {
		/* array of a1->edges - ei of flags */
		GArray* ssep_ei = g_ptr_array_index(a1->ssep,ei);
		GPtrArray* used_ei = g_ptr_array_index(h->used,ei);
		for (guint j = 0; j < ssep_ei->len; j++)
			if (g_array_index(ssep_ei,size_t,j) > 0) {
				// ei and ej are separated in a1
				size_t ej = ei + j + 1;
				GPtrArray* used_ej = g_ptr_array_index(h->used,ej);
				// check that all edges in a2 image of ei (e2k)
				// are separated from all edges in a2 image of ej (e2l)
				for (guint e2k = 0; e2k < a2->edges->len; e2k++) {
					noll_edgflds_t* used_ei_e2k =
							g_ptr_array_index(used_ei,e2k);
					if (used_ei_e2k == NULL)
						continue; /* e2k is not used */
					for (guint e2l = 0; e2l < a2->edges->len; e2l++)
						if (e2k != e2l) {
							noll_edgflds_t* used_ej_e2l =
									g_ptr_array_index(used_ej,e2l);
							if (used_ej_e2l == NULL)
								continue; /* e2l is not used */
							// here, e2k is used in ei and e2l in ej
							// check if they are separated in a2
							if (e2k == e2l)
								return FALSE;
							size_t min = (e2k < e2l) ? e2k : e2l;
							size_t max = (e2k < e2l) ? e2l : e2k;
							GArray* ssep_min = g_ptr_array_index(a2->ssep,min);
							if (g_array_index(ssep_min,size_t,max-min-1) == 0)
								return FALSE;
						}
				}
			}
	}
	return TRUE;
}

gboolean noll_hom_is_isom(noll_hom_t* h) {
	if (h == NULL || h->vrt == NULL)
		return FALSE;
	gboolean isom = TRUE;
	guint size = h->vrt->len;
	for (guint i = 0; i < size && (isom == TRUE); i++) {
		size_t v2_i = g_array_index(h->vrt,size_t,i);
		for (guint j = 0; j < size && (isom == TRUE); j++)
			if (j != i) {
				size_t v2_j = g_array_index(h->vrt,size_t,j);
				isom = (v2_i == v2_j) ? FALSE : TRUE;
			}
	}
	return isom;
}
