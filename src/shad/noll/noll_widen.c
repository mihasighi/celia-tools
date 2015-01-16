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

#include <string.h>
#include "noll.h"
#include "noll_fun.h"
#include "noll_internal.h"
#include "glib.h"

/* ============================================================ */
/* III.4 Widening */
/* ============================================================ */

noll_val_t *
noll_widening(sh_manager_t * man, noll_val_t * a1, noll_val_t * a2) {
	return NULL; // TODO
}

noll_val_t *
noll_fold(sh_manager_t * man, gboolean destructive, noll_val_t * a,
		size_t annon, GArray *preds) {
	size_t i;
	noll_val_t *r = noll_top(man, a->frame);

	for (i = 0; i < a->set->len; i++) {
		noll_graph_t *gr = noll_graph_fold(man, g_ptr_array_index(a->set,i),
				annon, preds);
		if (gr != NULL && noll_graph_is_true(man, gr) == FALSE) {
			g_ptr_array_add(r->set, gr);
		}
	}
	//TODO: eliminate homomorphic graphs in r->set
	if (destructive)
		noll_free(man, a);
	return r;
}

GPtrArray *node2var(noll_graph_t* a) {
	GPtrArray * vtx2var = g_ptr_array_sized_new(a->vrt_size);
	guint vi;
	for (vi = 0; vi < a->var2vrt->len; vi++) {
		size_t ni = g_array_index(a->var2vrt,size_t,vi);
		//	printf("%zu ",ni);

		if (g_ptr_array_index(vtx2var,ni) == NULL) {
			GPtrArray *aux = g_ptr_array_new();
			g_ptr_array_add(aux, GUINT_TO_POINTER(vi));
			g_ptr_array_index(vtx2var,ni) = aux;
		} else {
			GPtrArray *aux = g_ptr_array_index(vtx2var,ni);
			g_ptr_array_add(aux, GUINT_TO_POINTER(vi));
		}

	}
	//	printf("\n "); fflush(stdout);
	//
	//	for(guint i = 0; i<a->vrt_size; i++){
	//		GPtrArray* aux = (GPtrArray*)g_ptr_array_index(vtx2var,i);
	//		for(guint j=0;j<aux->len;j++){
	//			size_t vi = GPOINTER_TO_UINT(g_ptr_array_index(aux,j));
	//			printf("%zu ",vi);
	//		}
	//		printf("\n ");
	//		fflush(stdout);
	//	}
	return vtx2var;

}

/*
 * takes a graph as parameter and a set of edges, redges.
 * It replaces the subgraph defined by the edges in redges
 * with a predicate edge src == P ==>dst,
 * where the predicate P is identified by its id, i.e., pid.
 */
void replace_unfolding_with_pred(sh_manager_t *man, noll_graph_t *g,
		GPtrArray *redges, guint src, guint dst, guint pid) {

	// TODO: eliminate also nodes if they are isolated

	size_t i;
	for (i = 0; i < redges->len; i++) {
		//Step1: remove from g->e_mat;
		noll_edge_t* e = g_ptr_array_index(redges,i);

		GArray* edg_v = g_ptr_array_index(g->e_mat,e->vrt_src);
		if (edg_v != NULL)
			for (guint j = 0; j < edg_v->len; j++) {
				size_t ce = g_array_index(edg_v,size_t,j);
				if (ce == e->eid) {
					edg_v = g_array_remove_index(edg_v, j);
#ifndef NDEBUG1
					printf("\n is %zu", ce);
					printf("\n equal with ");
					noll_edge_fdump(stdout,man,e);
					fflush(stdout);
#endif
				}
			}
		//Step2: remove from g->e_rmat
		edg_v = g_ptr_array_index(g->e_rmat,e->vrt_dst);
		if (edg_v != NULL)
			for (guint j = 0; j < edg_v->len; j++) {
				size_t ce = (size_t) g_array_index(edg_v,size_t,j);
				if (ce == e->eid) {
					edg_v = g_array_remove_index(edg_v, j);
#ifndef NDEBUG1
					printf("\n is %d ",ce);
					printf("\n equal with ");
					noll_edge_fdump(stdout,man,e);
					fflush(stdout);
#endif
				}
			}
		//Step3: remove from edges g->edges
		g_ptr_array_index(g->edges,e->eid) = NULL;
	}

	//Step: permutation on edges due due remove edge
	// Step 3: and remove from g->edges
	GArray *eids = g_array_new(FALSE, FALSE, sizeof(gint));
	gint av = -1;
	for (guint j = 0; j < g->edges->len; j++) {
		g_array_append_val(eids,av);
	}
#ifndef NDEBUG1
	for(guint i=0;i<eids->len; i++) {
		printf("%d ", g_array_index(eids,gint,i));
	}
	fflush(stdout);
#endif

	size_t gap = 0;
	for (guint i = 0; i < g->edges->len; i++) {
		if (g_ptr_array_index(g->edges,i) == NULL) {
			printf("%d", i);
			fflush(stdout);
			g_ptr_array_remove_index(g->edges, i);
			i = i - 1;
			gap++;
		} else
			g_array_index(eids,gint,i+gap) = i;
	}
	for (guint i = 0; i < eids->len; i++) {
		printf(" eids[%d]==%d ", i, g_array_index(eids,gint,i));
	}
	fflush(stdout);

#ifndef NDEBUG1
	printf("\n intermediary step of replace ");
	fflush(stdout);
	noll_graph_fdump(stdout, man, g);
#endif
	apply_edge_renameing(g, eids);

	//Step3: create a new predicate edge
	noll_edge_t* ne = g_new(noll_edge_t,1);
	ne->eid = g->edges->len;
	ne->kind = NOLL_EDGE_LS;
	ne->vrt_src = src;
	ne->vrt_dst = dst;

	ne->info.pred.isrev = TRUE;
	ne->info.pred.pid = pid;

	ne->info.pred.ngb = g_array_new(FALSE, FALSE, sizeof(size_t));
	g_array_append_val(ne->info.pred.ngb,src);
	g_array_append_val(ne->info.pred.ngb,dst);

	/* add to array of edges */
	g_ptr_array_add(g->edges, ne);

	/* add edge to adjacency matrix for src and dst,
	 * WARNING: each entry of e_mat is an GArray* */
	g_array_append_val(g_ptr_array_index(g->e_mat,src), ne->eid);
	g_array_append_val(g_ptr_array_index(g->e_rmat,dst), ne->eid);

#ifndef NDEBUG
	printf("\n replace unfolding returns \n ");
	fflush(stdout);
	noll_graph_fdump(stdout, man, g);
#endif

}

void apply_edge_renameing(noll_graph_t *g, GArray *a) {
	guint i, j;
	// rename g->e_mat
	for (i = 0; i < g->e_mat->len; i++) {
		GArray* edgs = g_ptr_array_index(g->e_mat,i);
		if (edgs != NULL) {
			for (j = 0; j < edgs->len; j++) {
				size_t eid = g_array_index(edgs,size_t,j);
				gint n_eid = g_array_index(a,gint,eid);

				if (g_array_index(a,gint,eid) != -1)
					g_array_index(edgs,size_t,j)
							= (size_t) g_array_index(a,gint,eid);
				else {
					printf("\n ERROR!!!\n");
					fflush(stdout);
				}
			}
		}
	}
	//rename g->e_rmat
	for (i = 0; i < g->e_rmat->len; i++) {
		GArray* edgs = g_ptr_array_index(g->e_rmat,i);
		if (edgs != NULL) {
			for (j = 0; j < edgs->len; j++) {
				size_t eid = g_array_index(edgs,size_t,j);
				gint n_eid = g_array_index(a,gint,eid);
				if (g_array_index(a,gint,eid) != -1)
					g_array_index(edgs,size_t,j)
							= (size_t) g_array_index(a,gint,eid);
				else {
					printf("\n ERROR!!!\n");
					fflush(stdout);
				}
			}
		}
	}

	//rename g->edges
	for (i = 0; i < g->edges->len; i++) {
		noll_edge_t* e = g_ptr_array_index(g->edges,i);
		if (g_array_index(a,gint,e->eid) != -1)
			e->eid = (size_t) g_array_index(a,gint,(guint)e->eid);
		else {
			printf("\n ERROR!!!\n");
			fflush(stdout);
		}
	}
}

noll_graph_t *
noll_graph_fold(sh_manager_t * man, noll_graph_t * a, size_t annon,
		GArray *preds) {
	guint vi;

	//noll_graph_t * a = noll_graph_copy(aa);
	GPtrArray * vtx2var = node2var(a);

	//Step0 : select predicate pid and nodes, a tuple of nodes
	//Step1: check that pid(nodes) holds in the graph a

	guint i;
	//Step1: TODO: selection of nodes/predicates
	guint pid = 0;

	sh_pred_t * p = g_ptr_array_index(man->typenv->pinfo, pid);

	/*******************  First instance of P *************************************************/
	GArray * nodes = g_array_sized_new(FALSE, FALSE, sizeof(gint), ((p->arity)
			+ (p->ex_vars)));
	gint val = 2;
	size_t ll;
	g_array_insert_val(nodes,0,val);
	val = 3;
	g_array_insert_val(nodes,1,val);

	val = -1;
	for (ll = p->arity; ll < p->arity + p->ex_vars; ll++)
		g_array_insert_val(nodes,ll,val);

	//Step2
	GPtrArray* new_redges = check_unfolding_of_pred(nodes, a, man, pid);

	if (new_redges != NULL) {
		for (size_t ll = 0; ll < new_redges->len; ll++) {
			noll_edge_t* ee = g_ptr_array_index(new_redges,ll);
			fprintf(stdout, "edge[%d]->id = %d", ll, ee->eid);
			fflush(stdout);
		}
	}

	replace_unfolding_with_pred(man, a, new_redges, 2, 3, 0);

	g_array_free(nodes, TRUE);
	// g_array_free(new_redges,TRUE);

	/*******************  Second instance of P  *************************************************/

	nodes = g_array_sized_new(FALSE, FALSE, sizeof(gint), ((p->arity)
			+ (p->ex_vars)));
	val = 3;
	g_array_insert_val(nodes,0,val);
	val = 2;
	g_array_insert_val(nodes,1,val);

	val = -1;
	for (ll = p->arity; ll < p->arity + p->ex_vars; ll++)
		g_array_insert_val(nodes,ll,val);

	//Step2
	new_redges = check_unfolding_of_pred(nodes, a, man, pid);

	if (new_redges != NULL) {
		for (size_t ll = 0; ll < new_redges->len; ll++) {
			noll_edge_t* ee = g_ptr_array_index(new_redges,ll);
			fprintf(stdout, "edge[%d]->id = %d", ll, ee->eid);
			fflush(stdout);
		}
	}

	replace_unfolding_with_pred(man, a, new_redges, 3, 2, 0);

	g_array_free(nodes, TRUE);
	//	g_array_free(new_redges,TRUE);

	/*******************  End of instances  *************************************************/

	printf("\n fold_graph returns: \n");
	noll_graph_fdump(stdout, man, a);
	fflush(stdout);

	return a;

	//	sh_frame_t* frame = sh_stack_get_frame(sh_manager_get_varenv(man),
	//			a->frame);


	//	size_t vrt_typ = g_array_index(a->vrt_type,size_t,vrt);
	//	if (sh_typenv_is_type_field(sh_manager_get_typenv(man), vrt_typ, fld)
	//			== FALSE)
	//		return FALSE;
	//
	//	gboolean success = FALSE; // a pto exists in a
	//	gboolean found = FALSE; // something exists in a
	//	// graphs are put directly in res
	//	/* 1. Search in vertices from vrt */
	//	GArray* edges_vrt = g_ptr_array_index(a->e_mat,vrt);
}

GPtrArray * check_unfolding_of_pred(GArray * nodes, noll_graph_t *g,
		sh_manager_t * man, guint pid) {

	sh_pred_t * p = g_ptr_array_index(man->typenv->pinfo, pid);

	GPtrArray * pdef_mat = get_mat(man, p);

	g_assert((p->arity) + (p->ex_vars) == nodes->len);

	//nodes[i] = x[i], where x[i] is the i-th argument of the predicate i
	//nodes->len == p->arity;

	//noll_graph_t *g = noll_graph_copy(gg);
	GPtrArray* redges = NULL; // edges to replace
	redges = g_ptr_array_new();

	/* check */
	size_t i, j;

	/* for each edge in the predicate look for a match in the graph */

	/* start with node[0];
	 * use a queue for the explored nodes
	 * once no edge is mapped backtrack
	 * until Q is empty;
	 * */

	gint xn;
	size_t pn;
	GQueue *q = g_queue_new();//queue of nodes in the predicate

	/* add predicate parameters to explore */
	for (i = 0; i < p->arity; i++) {
		g_queue_push_head(q, GSIZE_TO_POINTER(i));
	}
	//re-allocate nodes complete it with -1
	while (!g_queue_is_empty(q)) {
		pn = GPOINTER_TO_SIZE (g_queue_pop_tail(q));
		//if(pn<p->arity){
		xn = g_array_index(nodes,gint,pn);
		//tot ce e in coada e si in nodes
		GPtrArray* eds = g_ptr_array_index(pdef_mat,pn);
		size_t nb_of_matched_edges = 0;
		for (j = 0; j < eds->len; j++) {
			//for each edge starting in n find match and add the other end to the queue
			sh_prededge_t *e = g_ptr_array_index(eds,j);
			if (e->src != pn)
				break;
			gint yn = g_array_index(nodes,size_t,e->dst);

			if (yn > -1) {
				//check that (nodes[e->src], nodes[e->dst]) belongs to the graph; e->dst is already mapped to one of the nodes in the graph g
				if (xn == yn && e->discr == SH_PFORM_PRED) {
					nb_of_matched_edges++;
					//carefull how it is called !
					// if e->dst && e->src are both less then p->parity is a problem
				} else {
					GArray * aux = g_ptr_array_index(g->e_mat,xn);
					gboolean found = FALSE;
					for (guint k = 0; k < aux->len && !found; k++) {
						size_t eid = g_array_index(aux,size_t, k); //the edge index
						noll_edge_t* edge = g_ptr_array_index(g->edges,eid);//get the actual edge

						if (edge->vrt_dst == yn)
							if (e->discr == SH_PFORM_PTO && edge->kind
									== NOLL_EDGE_PTO) {
								if (e->info.fid == edge->info.fid) {
									g_ptr_array_add(redges, edge);
									k = aux->len + 1;
									found = TRUE;
									nb_of_matched_edges++;
								}
							} else {
								if (e->discr == SH_PFORM_PRED && edge->kind
										== NOLL_EDGE_LS) {
									if (e->info.pred.pid == edge->info.pred.pid) {
										g_ptr_array_add(redges, edge);
										k = aux->len;
										found = TRUE;
										nb_of_matched_edges++;
									}
								}

							}
					}//end for k
					if (!found)
						return NULL;
					//one edge not matched so pred not unfolded
				}
			} else {
				GArray * aux = g_ptr_array_index(g->e_mat,xn);
				gboolean found = FALSE;
				for (guint k = 0; k < aux->len && !found; k++) {
					size_t eid = g_array_index(aux,size_t, k); //the edge index
					noll_edge_t* edge = g_ptr_array_index(g->edges,eid);//get the actual edge

					if (e->discr == SH_PFORM_PTO && edge->kind == NOLL_EDGE_PTO) {
						if (e->info.fid == edge->info.fid) {
							g_ptr_array_add(redges, edge);
							k = aux->len + 1;
							found = TRUE;
							g_array_index(nodes,gint,e->dst) = edge->vrt_dst;
							//TODO: don't insert twice the same element in the queue! non-termination in case predicate has cycles !
							//if(!gint_array_contains(nodes,e->dst))
							g_queue_push_head(q, e->dst);
							nb_of_matched_edges++;
						}
					} else {
						if (e->discr == SH_PFORM_PRED && edge->kind
								== NOLL_EDGE_LS) {
							if (e->info.pred.pid == edge->info.pred.pid) {
								g_ptr_array_add(redges, edge);
								k = aux->len + 1;
								found = TRUE;
								g_array_index(nodes,gint,e->dst)
										= edge->vrt_dst;
								//TODO: don't insert twice the same element in the queue! non-termination in case predicate has cycles !
								// if e->dst not in nodes add to the queue
								//if(!gint_array_contains(nodes,e->dst))
								g_queue_push_head(q, e->dst);
								nb_of_matched_edges++;
							}
						}
					}
				}//end find a match for e->dst in the graph. Add it to the queue for exploration.
			}//end case yn=-1
		}//end exploration of predicate edges starting in xn
		if (nb_of_matched_edges != eds->len)
			return NULL;
		//not all edges where matched from xn so no folding of p with these parameters
		//}else{
		//pn>=p->arity; pn is an annonimous
	}
	return redges;
}

gboolean gint_array_contains(GArray *a, gint val) {
	for (guint l = 0; l < a->len; l++)
		if (g_array_index(a,gint,l) == val)
			return TRUE;
	return FALSE;
}

GPtrArray* get_mat(sh_manager_t *man, sh_pred_t *p) {
	GPtrArray* newmat = g_ptr_array_new();

	guint j;
	for (j = 0; j < (p->arity + p->ex_vars); j++) {
		GPtrArray* ej = g_ptr_array_new();

		for (guint ll = 0; ll < p->def_mat->len; ll++) {
			GPtrArray* edges = g_ptr_array_index(p->def_mat,ll);
			if (edges != NULL) {
				for (guint i = 0; i < edges->len; i++) {
					sh_prededge_t *e = g_ptr_array_index(edges,i);

					if (e->src == j)
						g_ptr_array_add(ej, e);
					//					else if(e->discr == SH_PFORM_PRED){
					//						if(g_array_index(e->info.pred.args,size_t,0)==j)
					//							g_ptr_array_add(ej,e);
					//					}
				}
			}
		}
		g_ptr_array_add(newmat, ej);
	}

#ifndef NDEBUG
	for(guint l=0;l<newmat->len; l++) {
		fprintf(stdout,"\n");
		GPtrArray *a = g_ptr_array_index(newmat,l);
		if(a!=NULL) {
			for(guint k=0; k<a->len; k++) {
				sh_prededge_t* ak = g_ptr_array_index(a,k);
				fprintf(stdout, "(%d,%d) \t",ak->src,ak->dst);
				fflush(stdout);
			}
		}
	}
#endif
	return newmat;
}

