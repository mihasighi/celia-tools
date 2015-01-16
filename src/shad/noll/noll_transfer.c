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

#include <stdio.h>
#include <assert.h>

#include "noll.h"
#include "noll_fun.h"
#include "noll_internal.h"

/* ============================================================ */
/* III.2 Assignment */
/* ============================================================ */

noll_val_t*
noll_graph_assign_post_assume_eq(sh_manager_t * man, gboolean destructive,
		noll_graph_t * a, sh_stmt_t* stmt) {
	assert (NULL != stmt);
	assert (NULL != man);
	assert (stmt->kind == SH_STMT_ASSUME_EQ);
	/* remove edges between eq nodes, if possible */
	/* some materialization steps may be needed before */

	return NULL; // TODO
}

noll_val_t*
noll_graph_assign_post_assume_ne(sh_manager_t * man, gboolean destructive,
		noll_graph_t * a, sh_stmt_t* stmt) {
	assert (NULL != stmt);
	assert (NULL != man);
	assert (stmt->kind == SH_STMT_ASSUME_NE);

	/* introduce edge and then check local consistency */
	/* some materialization steps may be needed before */

	return NULL; // TODO
}

noll_val_t*
noll_graph_assign_post_new(sh_manager_t * man, gboolean destructive,
		noll_graph_t * a, sh_stmt_t* stmt) {
	assert (NULL != stmt);
	assert (NULL != man);
	assert (stmt->kind == SH_STMT_NEW);

	noll_val_t* res = NULL;
	/* add a vertex */
	noll_graph_t* gnew = noll_graph_add_vertex(man, destructive, a, stmt->tid);
	/* it is the last vertex */
	size_t vrt_new = gnew->vrt_size - 1;
	/* if no offsets, do vertex relabeling */
	if (stmt->info.binary.offset_l == NULL) {
		/* label the added vertex by stmt->info.binary.left */
		gnew = noll_graph_set_label(man, TRUE, gnew, stmt->info.binary.left,
				vrt_new);
		/* mark the new node different from any other node in the graph */
		for (guint i = 0; i < vrt_new; i++)
			gnew = noll_graph_add_edge_diff(man, TRUE, gnew, i, vrt_new);
		res = noll_join_graph(man, TRUE, (noll_val_t*) NULL, FALSE, gnew);
	} else {
		/* do materialization step for the left hand side */
		GArray* edges = g_array_new(FALSE, FALSE, sizeof(size_t));
		noll_val_t* r_mat = noll_graph_expand(man, destructive, gnew,
				stmt->info.binary.left, stmt->info.binary.offset_l, edges);
		/* execute new on each graph in the value and then join to the result */
		for (guint i = 0; i < r_mat->set->len; i++) {
			noll_graph_t* gi = g_ptr_array_index(r_mat->set,i);
			/* edge to be changed for this graph */
			size_t eid = g_array_index(edges, size_t, i);
			noll_edge_t* edge_i = g_ptr_array_index(gi->edges, eid);
			edge_i->vrt_dst = vrt_new; /* in place change! */
			/* mark the new node different from any other node in the graph */
			for (guint i = 0; i < vrt_new; i++)
				gnew = noll_graph_add_edge_diff(man, TRUE, gnew, i, vrt_new);
			res = noll_join_graph(man, TRUE, res, FALSE, gi);
			/* in place change of r with no copy of g */
		}
		/* all graphs in r_mat have been reused, reset entries to NULL */
		noll_free_set(r_mat);
	}

	return res;
}

noll_val_t*
noll_graph_assign_post_free(sh_manager_t * man, gboolean destructive,
		noll_graph_t * a, sh_stmt_t* stmt) {
	assert (NULL != stmt);
	assert (NULL != man);
	assert (stmt->kind == SH_STMT_FREE);

	return NULL; // TODO
}

noll_val_t*
noll_graph_assign_post_assign(sh_manager_t * man, gboolean destructive,
		noll_graph_t * a, sh_stmt_t* stmt) {

	assert (NULL != stmt);
	assert (NULL != man);
	assert (stmt->kind == SH_STMT_ASSIGN);

	noll_val_t* res = NULL;
	noll_graph_t* g = (destructive == TRUE) ? a : noll_graph_copy(a);
	size_t vrt_lhs_orig = noll_graph_get_label(g, stmt->info.binary.left);
	size_t vrt_rhs_orig = noll_graph_get_label(g, stmt->info.binary.right);

	/* edges obtained from materialization of lhs */
	GArray* edges_lhs = g_array_new(FALSE, FALSE, sizeof(size_t));
	/* graphs obtained from materialization of lhs */
	noll_val_t* r_mat_lhs = NULL;
	/* vertices obtained from materialization of rhs */
	GArray* vrt_rhs = g_array_new(FALSE, FALSE, sizeof(size_t));
	/* graphs obtained from materialization of rhs */
	noll_val_t* r_mat_rhs = NULL;

	/* lhs */
	if (stmt->info.binary.offset_l == NULL) {
		// no materialisation
		r_mat_lhs = noll_join_graph(man, TRUE, r_mat_lhs, FALSE, g);
	} else {
		/* do materialization step for the left hand side */
		/* edges to be changed by the assignment are in edges */
		r_mat_lhs = noll_graph_expand(man, destructive, g,
				stmt->info.binary.left, stmt->info.binary.offset_l, edges_lhs);
		if (r_mat_lhs != NULL)
			assert(r_mat_lhs->set->len == edges_lhs->len);
	}
	gboolean err = (r_mat_lhs == NULL) ? TRUE : FALSE;

	post_assign_err: if (err == TRUE) {
		if (res != NULL) {
			noll_free(man, res);
			res = NULL;
		}
		if (g != NULL) {
			noll_graph_free(g);
			g = NULL;
		}
		if (edges_lhs != NULL) {
			g_array_free(edges_lhs, TRUE);
			edges_lhs = NULL;
		}
		if (vrt_rhs != NULL) {
			g_array_free(vrt_rhs, TRUE);
			vrt_rhs = NULL;
		}
		if (r_mat_lhs != NULL) {
			noll_free(man, r_mat_lhs);
			r_mat_lhs = NULL;
		}
		if (r_mat_rhs != NULL) {
			noll_free(man, r_mat_rhs);
			r_mat_rhs = NULL;
		}
		return NULL;
	}

	/* rhs */
	if (stmt->info.binary.offset_r == NULL) {
		// no rhs materialization
		// add r_mat_lhs->set->len elements to vrt_rhs to continue
		for (guint i = 0; i < r_mat_lhs->set->len; i++)
			g_array_append_val(vrt_rhs, vrt_rhs_orig);
		r_mat_rhs = r_mat_lhs;
	} else {
		/* do materialization step for the right hand side */
		GArray* edges_lhs_new = g_array_new(FALSE, FALSE, sizeof(size_t));
		for (guint i = 0; i < r_mat_lhs->set->len; i++) {
			// the graph materialized for lhs
			noll_graph_t* gi = g_ptr_array_index(r_mat_lhs->set,i);
			// the corresponding edge materialized in gi
			size_t eid_lhs_i = g_array_index(edges_lhs,size_t,i);
			// prepare rhs materialisation
			GArray * edges_rhs_i = g_array_new(FALSE, FALSE, sizeof(size_t));
			noll_val_t* r_mat_rhs_i = noll_graph_expand(man, destructive, gi,
					stmt->info.binary.right, stmt->info.binary.offset_r,
					edges_rhs_i);
			if (r_mat_rhs_i == NULL) {
				err = TRUE;
				g_array_free(edges_rhs_i, TRUE);
				g_array_free(edges_lhs_new, TRUE);
				goto post_assign_err;
			}
			// go through the result and add the graph built and its edge/vertex
			for (guint j = 0; j < r_mat_rhs_i->set->len; j++) {
				noll_graph_t* gj = g_ptr_array_index(r_mat_rhs_i->set,j);
				size_t eid_rhs_j = g_array_index(edges_rhs_i,size_t,j);
				noll_edge_t* e_rhs_j = g_ptr_array_index(gj->edges,eid_rhs_j);
				assert (e_rhs_j != NULL);
				size_t vrt_rhs_j = e_rhs_j->vrt_dst;
				// add the graph
				r_mat_rhs = noll_join_graph(man, TRUE, r_mat_rhs, FALSE, gj);
				// add the lhs edge
				g_array_append_val(edges_lhs_new,eid_lhs_i);
				// add the rhs node
				g_array_append_val(vrt_rhs,vrt_rhs_j);
			}
			// free intermediate results
			g_array_free(edges_rhs_i, TRUE);
			noll_free_set(r_mat_rhs_i); // do not free the graphs
		}
		// end of rhs materialisation
		// results in r_mat_rhs, edges_lhs_new, vrt_rhs
		assert (r_mat_rhs->set->len == edges_lhs_new->len);
		assert (vrt_rhs->len == edges_lhs_new->len);
		g_array_free(edges_lhs, TRUE);
		edges_lhs = edges_lhs_new;
	}
	// do the update of dst for edges in edges_lhs to vrt in vrt_rhs
	// fill the result and check for garbage
	for (guint i = 0; i < edges_lhs->len; i++) {
		noll_graph_t* gi = g_ptr_array_index(r_mat_rhs->set,i);
		size_t ei = g_array_index(edges_lhs,size_t,i);
		size_t vrti = g_array_index(vrt_rhs,size_t,i);
		gi = noll_graph_set_edge_dst(man, TRUE, gi, ei, vrti);
		if (gi == NULL) {
			err = TRUE;
			goto post_assign_err;
		}
		res = noll_join_graph(man, TRUE, res, FALSE, gi);
	}
	// free the additional memory
	g_array_free(edges_lhs, TRUE);
	g_array_free(vrt_rhs, TRUE);
	noll_free_set(r_mat_rhs);
	return res;
}

noll_val_t*
noll_graph_assign_post_pcall(sh_manager_t * man, gboolean destructive,
		noll_graph_t * a, sh_stmt_t* stmt) {
	assert (NULL != stmt);
	assert (NULL != man);
	assert (stmt->kind == SH_STMT_PCALL);

	/* do the call, change of frame, etc */

	return NULL; // TODO
}

noll_val_t*
noll_graph_assign_post_preturn(sh_manager_t * man, gboolean destructive,
		noll_graph_t * a, sh_stmt_t* stmt) {
	assert (NULL != stmt);
	assert (NULL != man);
	assert (stmt->kind == SH_STMT_PRETURN);

	/* do the call, change of frame, etc */

	return NULL; // TODO
}

noll_val_t *
noll_graph_assign_post(sh_manager_t * man, gboolean destructive,
		noll_graph_t * a, sh_stmt_t* stmt) {
	if (a == NULL)
		return NULL;
	assert (NULL != stmt);
	assert (NULL != man);

#ifdef NDEBUG
	fprintf(stdout, "\n==== post(");
	sh_stmt_fdump(stdout, man, stmt);
	fprintf(stdout, ", ");
	noll_graph_fdump(stdout, man, a);
	fprintf(stdout, ") = ");
#endif

	noll_val_t* r = NULL;
	switch (stmt->kind) {
	case SH_STMT_ASSUME_EQ:
		r = noll_graph_assign_post_assume_eq(man, destructive, a, stmt);
		break;
	case SH_STMT_ASSUME_NE:
		r = noll_graph_assign_post_assume_ne(man, destructive, a, stmt);
		break;
	case SH_STMT_NEW:
		r = noll_graph_assign_post_new(man, destructive, a, stmt);
		break;
	case SH_STMT_FREE:
		r = noll_graph_assign_post_free(man, destructive, a, stmt);
		break;
	case SH_STMT_ASSIGN:
		r = noll_graph_assign_post_assign(man, destructive, a, stmt);
		break;
	case SH_STMT_PCALL:
		r = noll_graph_assign_post_pcall(man, destructive, a, stmt);
		break;
	case SH_STMT_PRETURN:
		r = noll_graph_assign_post_preturn(man, destructive, a, stmt);
		break;
	default:
		assert(TRUE);
	}
#ifdef NDEBUG
	noll_fdump(stdout, man, r);
	fprintf(stdout, "\n====\n");
#endif

	return r;
}

noll_val_t *
noll_assign_post(sh_manager_t * man, gboolean destructive, noll_val_t * a,
		sh_stmt_t* stmt, noll_val_t * dest) {
	if (a == NULL)
		return NULL;
	assert (NULL != stmt);
	assert (NULL != man);

	noll_val_t* r = NULL;
	/* applies post on each graph in the value and then join to the result */
	for (guint i = 0; i < a->set->len; i++) {
		noll_val_t* r1 = noll_graph_assign_post(man, destructive,
				g_ptr_array_index(a->set,i), stmt);
		r = noll_join(man, TRUE, r, r1);
	}
	if (dest != NULL && noll_is_top(man, dest) == TRUE)
		r = noll_meet(man, TRUE, r, dest); // TODO: see if needed
	return r;
}

noll_val_t *
noll_assign_pre(sh_manager_t * man, gboolean destructive, noll_val_t * a,
		sh_stmt_t* stmt, noll_val_t * dest) {
	return NULL; // TODO
}

noll_val_t *
noll_expand(sh_manager_t * man, gboolean destructive, noll_val_t * a, size_t x,
		size_t f) {
	return NULL; // TODO
}

noll_val_t*
noll_graph_expand(sh_manager_t* man, gboolean destructive, noll_graph_t* a,
		size_t x, GArray* offs, GArray* edges) {
	if (a == NULL)
		return NULL;
	assert (offs != NULL); // done elsewhere

	noll_val_t* res = noll_top(man, a->frame);
	size_t vrt_x = noll_graph_get_label(a, x);
	size_t fld = g_array_index(offs,size_t,0);
	/* Materialization of the edge starting from vrt_x with field fld.
	 * After the first materialization step, we obtain
	 * an array of graphs and
	 * an array of edges last materialized in each graph
	 */
	GArray* edges_res = g_array_new(FALSE, FALSE, sizeof(size_t));
	gboolean success = noll_graph_expand_1(man, FALSE, a, vrt_x, fld, res,
			edges_res);
	// success == FALSE implies that zero expansion can been done
	graph_expand_err: if (success == FALSE) {
		if (res != NULL) {
			noll_free(man, res);
			res = NULL;
		}
		if (edges_res != NULL) {
			g_array_free(edges_res, TRUE);
			edges_res = NULL;
		}
		return NULL;
	}

	/* next steps uses the result of the previous steps */
	for (guint i = 1; i < offs->len; i++) {
		// field of expansion
		size_t fldi = g_array_index(offs,size_t,i);
		// expansion results
		noll_val_t* res_tmp = noll_top(man, a->frame);
		GArray* edges_tmp = g_array_new(FALSE, FALSE, sizeof(size_t));
		for (guint j = 0; j < res->set->len; j++) {
			noll_graph_t* gj = g_ptr_array_index(res->set,j);
			size_t eidj = g_array_index(edges_res,size_t,j);
			noll_edge_t* ej = g_ptr_array_index(gj->edges, eidj);
			size_t vrtj = ej->vrt_dst;
			gboolean successj = noll_graph_expand_1(man, TRUE, gj, vrtj, fldi,
					res_tmp, edges_tmp);
			// successj == FALSE implies that zero expansion can be done for this graph
			// continue with the other graphs
		}
		assert (res_tmp->set->len == edges_tmp->len);
		if (res_tmp->set->len > 0) {
			// continue with this set of graphs
			noll_free(man, res); // graphs are not destroyed since reused in expand_1
			g_array_free(edges_res, TRUE);
			res = res_tmp;
			edges_res = edges_tmp;
		} else {
			// zero graphs generated at this step
			noll_free(man, res_tmp);
			g_array_free(edges_tmp, TRUE);
			goto graph_expand_err;
		}
	}
	// end of expansion, results in res and edges_res
	// fill edges with edges_res
	for (guint j = 0; j < res->set->len; j++) {
		size_t eidj = g_array_index(edges_res,size_t,j);
		g_array_append_val(edges,eidj);
	}
	g_array_free(edges_res, TRUE);

	return res;
}

gboolean noll_graph_expand_1(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a, size_t vrt, size_t fld, noll_val_t* res, GArray* edges) {
	assert (res->set->len == edges->len);
	// check that the vertex has the type of fld
	size_t vrt_typ = g_array_index(a->vrt_type,size_t,vrt);
	if (sh_typenv_is_type_field(sh_manager_get_typenv(man), vrt_typ, fld)
			== FALSE)
		return FALSE;

	gboolean success = FALSE; // a pto exists in a
	gboolean found = FALSE; // something exists in a
	// graphs are put directly in res
	/* 1. Search in vertices from vrt */
	GArray* edges_vrt = g_ptr_array_index(a->e_mat,vrt);
	if (edges_vrt != NULL) {
		/* 1.1: search for pto or predicate edges with field fld */
		for (guint i = 0; (i < edges_vrt->len) && (success == FALSE); i++) {
			size_t eid = g_array_index(edges_vrt,size_t,i);
			noll_edge_t* e = g_ptr_array_index(a->edges, eid);
			if (e->kind == NOLL_EDGE_PTO && e->info.fid == fld) {
				noll_graph_t* g = noll_graph_copy(a);
				g_ptr_array_add(res->set, g);
				g_array_append_val(edges,eid);
				success = TRUE;
				// since resulting graphs are deterministic,
				// mark with success that a precise unfolding is found
			}
		}

		/* 1.2: if a pto exist, the other pred edges starting from vrt
		 * shall be collapsed
		 */
		if (success == TRUE) {
			gboolean docollapse = FALSE;
			// works on the last added graph in place
			size_t lst = res->set->len - 1;
			noll_graph_t* g = g_ptr_array_index(res->set,lst);
			size_t eid_new = g_array_index(edges,size_t,lst);
			size_t vrt_new = vrt;
			do {
				GArray* edges_new = g_ptr_array_index(g->e_mat,vrt_new);
				for (guint i = 0; i < edges_new->len; i++) {
					size_t eid = g_array_index(edges_new,size_t,i);
					noll_edge_t* e = g_ptr_array_index(g->edges, eid);
					if (e->kind == NOLL_EDGE_LS) {
						size_t pid = e->info.pred.pid;
						if (sh_typenv_is_pred_field0(
								sh_manager_get_typenv(man), pid, fld, TRUE)
								== TRUE) {
							// do collapse and update vrt and edge value
							noll_graph_collapse_edge(man, TRUE, g, eid,
									&vrt_new, &eid_new);
							docollapse = TRUE;
						}

					}
				}
			} while (docollapse == TRUE);
			// update the last edge
			g_array_index(edges,size_t,lst) = eid_new;
		} else {
			/* 1.3: otherwise, try to unfold predicate edges */
			for (guint i = 0; i < edges_vrt->len; i++) {
				size_t eid = g_array_index(edges_vrt,size_t,i);
				noll_edge_t* e = g_ptr_array_index(a->edges, eid);
				if (e->kind == NOLL_EDGE_LS) {
					// Case 1: search to obtain the edge from this predicate
					size_t pid = e->info.pred.pid;
					if (sh_typenv_is_pred_field0(sh_manager_get_typenv(man),
							pid, fld, TRUE) == TRUE) {
						noll_graph_t* a_new = noll_graph_unfold_start(man,
								FALSE, a, eid);
						// search on the resulting graph the expansion
						found = found || noll_graph_expand_1(man, TRUE, a_new,
								vrt, fld, res, edges);
					}
					// Case 2: the edge is obtained from collapsing this edge
					size_t tid = g_array_index(a->vrt_type,size_t,e->vrt_dst);
					if (sh_typenv_is_type_field(sh_manager_get_typenv(man),
							tid, fld) == TRUE) {
						size_t vrt_new = vrt;
						noll_graph_t* a_new = noll_graph_collapse_edge(man,
								FALSE, a, eid, &vrt_new, NULL);
						found = found || noll_graph_expand_1(man, TRUE, a_new,
								vrt_new, fld, res, edges);
					}
				}
			}
		}
	}

	if (success == TRUE || found == TRUE)
		goto expand_1_true;

	/* 2. Search in vertices to vrt_x with predicates type dll */
	edges_vrt = g_ptr_array_index(a->e_rmat,vrt);
	if (edges_vrt != NULL) {
		// search predicate edges with reverse field fld
		for (guint i = 0; i < edges_vrt->len; i++) {
			size_t eid = g_array_index(edges_vrt,size_t,i);
			noll_edge_t* e = g_ptr_array_index(a->edges, eid);
			if (e->kind == NOLL_EDGE_LS) {
				size_t pid = e->info.pred.pid;
				if (sh_typenv_is_pred_field0(sh_manager_get_typenv(man), pid,
						fld, FALSE) == TRUE) {
					noll_graph_t* a_new = noll_graph_unfold_end(man, FALSE, a,
							eid);
					// search on the resulting graph the expansion
					found = found || noll_graph_expand_1(man, TRUE, a_new, vrt,
							fld, res, edges);
				}
			}
		}
	}

	if (found == TRUE)
		goto expand_1_true;
	/* TODO: 3. Search in the sharing constraints */
	/* 4. If not found, it means that the edge is not initialized
	 *    then add it in the graph pointing to NULL
	 */
	size_t eid = 0;
	noll_graph_t* g = noll_graph_add_edge_pto(man, FALSE, a, vrt, fld, 0, &eid);
	g_ptr_array_add(res->set, g);
	g_array_append_val(edges,eid);

	expand_1_true: if (destructive)
		noll_graph_free(a);
	return TRUE;
}
