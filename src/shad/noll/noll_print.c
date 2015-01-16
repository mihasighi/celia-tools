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
 * noll_print.c: Printing function
 */

#include "noll.h"
#include "noll_fun.h"
#include "noll_internal.h"

/* ============================================================ */
/* I.3 Printing */
/* ============================================================ */

void noll_edge_fdump(FILE * stream, sh_manager_t * man, noll_edge_t * a) {
	if (a == NULL) {
		fprintf(stream, "(edge:null)");
		return;
	}
	fprintf(stream, "(edge:%d) ", a->eid);
	switch (a->kind) {
	case NOLL_EDGE_DIFF:
		fprintf(stream, "n%d != n%d", a->vrt_src, a->vrt_dst);
		break;
	case NOLL_EDGE_PTO:
		fprintf(stream, "n%d --f%d--> n%d", a->vrt_src, a->info.fid, a->vrt_dst);
		break;
	case NOLL_EDGE_LS:
		fprintf(stream, "n%d ==", a->vrt_src);
		if (a->info.pred.isrev == TRUE)
			fprintf(stream, "!");
		fprintf(stream, "p%d", a->info.pred.pid);
		if (a->info.pred.ngb->len > 0) {
			fprintf(stream, "(");
			GArray* args = a->info.pred.ngb;
			for (guint i = 0; i < a->info.pred.ngb->len; i++){
				size_t x = g_array_index(args,size_t,i);
			    fprintf(stream, "n%d,", g_array_index(args,size_t,i));
			}
			fprintf(stream, ")");
		}
		fprintf(stream, "==> n%d", a->vrt_dst);
		break;
	default:
		fprintf(stream, "unknown(n%d,n%d)", a->vrt_src, a->vrt_dst);
		break;
	}
}

void noll_share_fdump(FILE * stream, sh_manager_t * man, noll_share_t * a) {
	if (a == NULL) {
		fprintf(stream, "(null)");
		return;
	}
	fprintf(stream, "(over t%d loc) ", a->tid);
	/* print term t1 */
	if (a->vars_t1 != NULL) {
		for (guint i = 0; i < a->vars_t1->len; i++)
			fprintf(stream, "{v%d} U ", g_array_index(a->vars_t1,size_t,i));
	}
	if (a->lvars_t1 != NULL) {
		for (guint i = 0; i < a->lvars_t1->len; i++)
			fprintf(stream, "{l%d} U ", g_array_index(a->lvars_t1,size_t,i));
	}
	if (a->svars_t1 != NULL) {
		for (guint i = 0; i < a->svars_t1->len; i++)
			fprintf(stream, "s%d U ", g_array_index(a->svars_t1,size_t,i));
	}
	fprintf(stream, "%c", 0xe2);
	/* print operation */
	switch (a->op) {
	case NOLL_OP_IN:
		fprintf(stream, " in ");
		break;
	case NOLL_OP_SUBSETEQ:
		fprintf(stream, " <= ");
		break;
	case NOLL_OP_EQ:
		fprintf(stream, " = ");
		break;
	default:
		fprintf(stream, " (unknown) ");
		break;
	}
	/* print term t2 */
	if (a->vars_t2 != NULL) {
		for (guint i = 0; i < a->vars_t2->len; i++)
			fprintf(stream, "{v%d} U ", g_array_index(a->vars_t2,size_t,i));
	}
	if (a->lvars_t2 != NULL) {
		for (guint i = 0; i < a->lvars_t2->len; i++)
			fprintf(stream, "{l%d} U ", g_array_index(a->lvars_t2,size_t,i));
	}
	if (a->svars_t2 != NULL) {
		for (guint i = 0; i < a->svars_t2->len; i++)
			fprintf(stream, "s%d U ", g_array_index(a->svars_t2,size_t,i));
	}
	fprintf(stream, "%c", 0xe2);
}

void noll_graph_fdump(FILE * stream, sh_manager_t * man, noll_graph_t * a) {
	if (a == NULL)
		fprintf(stream, "false");
	/* for printing variables */
	sh_frame_t* frame = g_ptr_array_index(man->varenv->frames,a->frame);
	fprintf(stream, "Graph {\n(frame = %d) vars = (", a->frame);
	sh_frame_fdump(stream, frame, sh_manager_get_typenv(man),
			sh_manager_get_filenv(man));
	fprintf(stream, ")\n");
	/* print number of nodes */
	fprintf(stream, "nodes = %d,\n", a->vrt_size);
	/* print labeling of vars */
	fprintf(stream, "labeling = [NULL(n0)");
	for (size_t v = 0; v < frame->vars->len; v++) {
		fprintf(stream, ", %s(n%d)", sh_frame_get_varname(frame, v),
				noll_graph_get_label(a, v));
	}
	if (a->lvars != NULL)
		for (size_t vp = 0; vp < a->lvars->len; vp++) {
			fprintf(stream, ", v%d'(n%d)", vp, noll_graph_get_label(a, vp));
		}
	fprintf(stream, "],\n");
	/* print preciseness */
	fprintf(stream, "precise ");
	if (a->isprecise == TRUE)
		fprintf(stream, "TRUE");
	else
		fprintf(stream, "FALSE");
	fprintf(stream, ",\n");
	/* print edges */
	fprintf(stream, "edges = (%d) [\n", (a->edges == NULL) ? 0 : a->edges->len);
	if (a->edges != NULL)
		for (guint e = 0; e < a->edges->len; e++) {
			noll_edge_fdump(stream, man, g_ptr_array_index(a->edges,e));
			fprintf(stream, ",\n");
		}
	fprintf(stream, "],\n");
#ifndef NDEBUG
	/* print adjacency matrices */
	fprintf(stream, "e_mat = (%d) [\n", (a->e_mat == NULL) ? 0 : a->e_mat->len);
	if (a->e_mat != NULL)
	for (guint v = 0; v < a->e_mat->len; v++) {
		fprintf(stream, "(n-%d): [", v);
		GArray* edg_v = g_ptr_array_index(a->e_mat,v);
		if (edg_v == NULL) fprintf(stream, "(null)");
		else
		for (guint e = 0; e < edg_v->len; e++)
		fprintf(stream, "(edge:%d),", g_array_index(edg_v,size_t,e));
		fprintf(stream, "],\n");
	}
	fprintf(stream, "],\n");
	fprintf(stream, "e_rmat = (%d) [\n", (a->e_rmat == NULL) ? 0 : a->e_rmat->len);
	if (a->e_rmat != NULL)
	for (guint v = 0; v < a->e_rmat->len; v++) {
		fprintf(stream, "(n-%d): [", v);
		GArray* edg_v = g_ptr_array_index(a->e_rmat,v);
		if (edg_v == NULL) fprintf(stream, "(null)");
		else
		for (guint e = 0; e < edg_v->len; e++)
		fprintf(stream, "(edge:%d),", g_array_index(edg_v,size_t,e));
		fprintf(stream, "],\n");
	}
	fprintf(stream, "],\n");
#endif
	/* print separation constraints */
	fprintf(stream, "strongly separated edges = (%d) [\n",
			((a->ssep == NULL) ? 0 : a->ssep->len));
	if (a->ssep != NULL) {
		for (guint r = 0; r < a->ssep->len; r++) {
			fprintf(stream, "{");
			GArray* rr = g_ptr_array_index(a->ssep,r);
			for (guint e = 0; e < rr->len; e++)
				fprintf(stream, "edge-%d,", e);
			fprintf(stream, "},\n");
		}
	}
	fprintf(stream, "],\n");
	/* print sharing constraints */
	fprintf(stream, "sharing = (%d) [\n",
			((a->wsep == NULL) ? 0 : a->wsep->len));
	if (a->wsep != NULL) {
		for (guint s = 0; s < a->wsep->len; s++) {
			noll_share_t* constr = g_ptr_array_index(a->wsep,s);
			noll_share_fdump(stream, man, constr);
			fprintf(stream, "  and \n");
		}
		fprintf(stream, "true\n");
	}
	fprintf(stream, "]\n}\n");
}

/*
 * Print a formula
 */
void noll_graph_fprint(FILE * stream, sh_manager_t * man, noll_graph_t * a) {
	if (a == NULL)
		fprintf(stream, "false");
	/* for printing variables */
	sh_frame_t* frame = g_ptr_array_index(man->varenv->frames,a->frame);
	/* print the equalities of variables */
	for (size_t n = 0; n < a->vrt_size; n++) {
		size_t nv = 0;
		size_t vinit = 0;
		for (size_t v = 0; v < frame->vars->len; v++) {
			if (g_array_index(a->var2vrt,size_t,v) == n) {
				if (nv == 0)
					vinit = v;
				else if (nv == 1)
					fprintf(stream, "%s == %s", sh_frame_get_varname(frame,
							vinit), sh_frame_get_varname(frame, v));
				else
					fprintf(stream, " == %s", sh_frame_get_varname(frame, v));
				nv++;
			}
		}
		if (nv > 1)
			fprintf(stream, " and \n");
	}
	fprintf(stream, "true,\n");

	/* TODO: print the spatial constraints */
	if (a->isprecise)
		fprintf(stream, "emp *\n");
	else
		fprintf(stream, "true *\n");

	/* TODO: print sharing constraints */

}

void noll_fprint(FILE * stream, sh_manager_t * man, noll_val_t * a) {
	if (a == NULL)
		fprintf(stream, "bottom");
	else if (a->set == 0 || a->set == NULL || a->set->len == 0)
		fprintf(stream, "top (frame=%d)", a->frame);
	else {
		fprintf(stream, "\{ (frame=%d, size=%d)\n", a->frame, a->set->len);
		noll_graph_fprint(stream, man, g_ptr_array_index(a->set,0));
		for (guint i = 1; i < a->set->len; i++) {
			fprintf(stream, "\nor\n");
			noll_graph_fprint(stream, man, g_ptr_array_index(a->set,i));
		}
	}
}

void noll_fdump(FILE * stream, sh_manager_t * man, noll_val_t * a) {
	if (a == NULL)
		fprintf(stream, "(null)");
	else if (a->set == 0 || a->set == NULL || a->set->len == 0)
		fprintf(stream, "(frame=%d,size=0)", a->frame);
	else {
		fprintf(stream, "\{ (frame=%d, size=%d)\n", a->frame, a->set->len);
		noll_graph_fdump(stream, man, g_ptr_array_index(a->set,0));
		for (guint i = 1; i < a->set->len; i++) {
			fprintf(stream, "\nor\n");
			noll_graph_fdump(stream, man, g_ptr_array_index(a->set,i));
		}
	}
}

