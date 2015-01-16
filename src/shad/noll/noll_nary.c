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

/* ********************************************************************** */
/* III. Operations */
/* ********************************************************************** */

/* ============================================================ */
/* III.1 Meet and Join */
/* ============================================================ */

noll_val_t *
noll_meet(sh_manager_t * man, gboolean destructive, noll_val_t * a1,
		noll_val_t * a2) {
	return NULL; // TODO
}

noll_val_t *
noll_join(sh_manager_t * man, gboolean destructive, noll_val_t * a1,
		noll_val_t * a2) {
	noll_val_t *r;
	size_t i, j;
	gboolean *isin1, *isin2;
	/* simple cases */
	if (a1 == NULL && a2 == NULL)
		return NULL;
	if (noll_is_bottom(man, a1) || noll_is_top(man, a2)) {
		if (destructive == TRUE) {
			noll_free(man, a1);
			return a2;
		} else
			return noll_copy(man, a2);
	}
	if (noll_is_bottom(man, a2) || noll_is_top(man, a1)) {
		if (destructive == TRUE) {
			noll_free(man, a2);
			return a1;
		} else
			return noll_copy(man, a1);
	}
#ifndef NDEBUG
	fprintf (stdout, "\n====shape_join:\n ");
#endif
	/* general case, i.e., each set has at least one element, do the union */
	r = noll_top(man, a1->frame);
	isin1 = (gboolean *) malloc(a1->set->len * sizeof(gboolean));
	memset(isin1, FALSE, a1->set->len * sizeof(gboolean));
	isin2 = (gboolean *) malloc(a2->set->len * sizeof(gboolean));
	memset(isin2, FALSE, a2->set->len * sizeof(gboolean));
	for (j = 0; j < a2->set->len; j++) {
		for (i = 0; i < a1->set->len; i++)
			if (isin1[i] == FALSE) {
				noll_graph_t *gr = noll_graph_join(man, FALSE,
						g_ptr_array_index(a1->set,i),
						g_ptr_array_index(a2->set,j));
				if (gr != NULL && noll_graph_is_true(man, gr) == FALSE) {
					isin1[i] = TRUE;
					isin2[j] = TRUE;
					g_ptr_array_add(r->set, gr);
				}
			}
		if (!isin2[j]) {
			// nothing added for a2->set[i]
			noll_graph_t *gr = noll_graph_copy(g_ptr_array_index(a2->set,j));
			g_ptr_array_add(r->set, gr);
		}
	}
	for (i = 0; i < a1->set->len; i++)
		if (isin1[i] == FALSE) {
			r = noll_join_graph(man, TRUE, r, TRUE,
					g_ptr_array_index(a1->set,i));
		}
	if (destructive == TRUE) {
		noll_free(man, a1);
		noll_free(man, a2);
	}
#ifndef NDEBUG1
	fprintf (stdout, "\n====shape_join returns: ");
	noll_fdump(stdout, man, r);
#endif
	return r;
}

noll_val_t *
noll_meet_cons_array(sh_manager_t * man, gboolean destructive, noll_val_t * a,
		GPtrArray* array) {
	return NULL; // TODO
}


noll_val_t *
noll_join_graph(sh_manager_t * man, gboolean destructive, noll_val_t * a1,
		gboolean copy, noll_graph_t * a2) {
	noll_val_t* res = (destructive == TRUE) ? a1 : noll_copy(man, a1);
	noll_graph_t* g2 = (copy == TRUE) ? noll_graph_copy(a2) : a2;
	if (res == NULL) {
		res = noll_top(man, a2->frame);
		g_ptr_array_add(res->set, g2);
		return res;
	}
	/* check that there is no graph in res greater that g2 */
	gboolean found = FALSE;
	for (guint i = 0; i < res->set->len && (found == FALSE); i++) {
		noll_graph_t* gi = g_ptr_array_index(res->set,i);
		if (noll_graph_is_le(man, g2, gi) == TRUE)
			found = TRUE;
	}
	if (found == TRUE)
		noll_graph_free(g2);
	else {
		g_ptr_array_add(res->set, g2);
	}
	return res;
}

noll_graph_t*
noll_graph_join(sh_manager_t* man, gboolean destructive,
		noll_graph_t* a1, noll_graph_t* a2) {

	if (noll_graph_is_le(man, a1, a2) == TRUE)
		return (destructive == TRUE) ? a2 : noll_graph_copy(a2);

	if (noll_graph_is_le(man, a2, a1) == TRUE)
		return (destructive == TRUE) ? a1 : noll_graph_copy(a1);

	return NULL;
}
