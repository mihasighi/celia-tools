/*
 * ucons_print.c
 *
 * Printing and serialization
 *
 * APRON Library / Shape Domain
 *
 * Copyright (C) LIAFA 2009
 *
 */

/*
 * This file is part of the APRON Library, released under LGPL license.
 * Please read the COPYING file packaged in the distribution.
 */

#include "ucons.h"
#include "ucons_fun.h"
#include "ucons_internal.h"
#include "apron2shape.h"

/* ============================================================ */
/* Printing */
/* ============================================================ */

void pattern_key_fprint(FILE * stream, ucons_internal_t *pr, pattern_key_t * a,
		char **name_of_dim) {
	size_t i, j, e_seg;
	if (!a)
		fprintf(stream, "(null)");
	else {
		// print using the form of the pattern
		switch (pr->PI[a->type].kind) {
		case pattern_1:
			fprintf(stream, "y1 in [%s], l[%s]>=2 ",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[0]]);
			break;
		case pattern_1_2:
			fprintf(stream, "y1, y2 in [%s] y2y1  l[%s]>=3 ",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[0]]);
			break;
		case pattern_2_1:
			fprintf(stream,
					"y1 in [%s], y2 in [%s], y1==y2 \n /\\ l[%s]>=2 /\\ l[%s]>=2 ",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[1]],
					name_of_dim[a->segments[0]], name_of_dim[a->segments[1]]);
			break;
		case pattern_3_1:
			fprintf(stream,
					"y1 in [%s], y2 in [%s], y3 in [%s], y1==y2==y3\n /\\ l[%s]>=2 /\\ l[%s]>=2 /\\ l[%s]>=2 ",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[1]],
					name_of_dim[a->segments[2]], name_of_dim[a->segments[0]],
					name_of_dim[a->segments[1]], name_of_dim[a->segments[2]]);
			break;
		case pattern_1_lx: {
			e_seg = pr->PI[a->type].e_seg;
			if (e_seg == 1) {
				fprintf(stream, "y1 in [%s] /\\ l[%s]>=2 /\\ y1=l[%s]",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]]);
			} else if (e_seg == 2)
				fprintf(stream,
						"y1 in [%s] /\\ l[%s]>=2 \n/\\ y1=l[%s] + l[%s]",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]]);
			else
				fprintf(stream,
						"y1 in [%s] /\\ l[%s]>=2 \n/\\ y1=l[%s] + l[%s] +l[%s]",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]],
						name_of_dim[a->segments[3]]);
			break;
		}
		case pattern_1_lx_1:
			fprintf(stream, "y1 in [%s] /\\ l[%s]>=2 /\\ y1=l[%s]-1",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[0]],
					name_of_dim[a->segments[1]]);
			break;
		case pattern_1_l1:
			fprintf(stream, "y1 in [%s] /\\ l[%s]>=2 /\\ y1=1",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[0]]);
			break;
		case pattern_2_1_lx: {
			e_seg = pr->PI[a->type].e_seg;
			if (e_seg == 1) {
				fprintf(stream,
						"y1 in [%s], y2 in [%s], \nl[%s]>=2 /\\ l[%s]>=2 /\\ y1 = l[%s] + y2",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]]);
			} else if (e_seg == 2) {
				fprintf(stream,
						"y1 in [%s], y2 in [%s],\n l[%s]>=2 /\\ l[%s]>=2 /\\ y1 = l[%s] + l[%s] + y2",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]],
						name_of_dim[a->segments[3]]);
			} else {
				fprintf(stream,
						"y1 in [%s], y2 in [%s],\n l[%s]>=2 /\\ l[%s]>=2 /\\ y1 = l[%s] + l[%s] + l[%s] + y2",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]],
						name_of_dim[a->segments[3]],
						name_of_dim[a->segments[4]]);
			}
			break;
		}
		case pattern_2_1_mlx: {
			e_seg = pr->PI[a->type].e_seg;
			if (e_seg == 1) {
				fprintf(stream,
						"y1 in [%s], y2 in [%s], \nl[%s]>=2 /\\ l[%s]>=2 /\\ y1 + l[%s] = y2",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]]);
			} else if (e_seg == 2) {
				fprintf(stream,
						"y1 in [%s], y2 in [%s],\n l[%s]>=2 /\\ l[%s]>=2 /\\ y1 + l[%s] + l[%s] = y2",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]],
						name_of_dim[a->segments[3]]);
			} else {
				fprintf(stream,
						"y1 in [%s], y2 in [%s],\n l[%s]>=2 /\\ l[%s]>=2 /\\ y1 + l[%s] + l[%s] + l[%s] = y2",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]],
						name_of_dim[a->segments[3]],
						name_of_dim[a->segments[4]]);
			}
			break;
		}
		default:
			fprintf(stream, "(pattern unknown)");
			break;
		}
	}
	fflush(stream);
}
void dot_pattern_key_fprint(FILE * stream, ucons_internal_t *pr,
		pattern_key_t * a, char **name_of_dim) {
	size_t i, j, e_seg;
	if (!a)
		fprintf(stream, "EMPTY");
	else {
		// print using the form of the pattern
		switch (pr->PI[a->type].kind) {
		case pattern_1:
			fprintf(stream, "forall y1. \n y1 in [%s], l[%s] >= 2 ",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[0]]);
			break;
		case pattern_1_2:
			fprintf(stream,
					"forall y1, y2. \n y1, y2 in [%s],\n y2 >= y1,\n l[%s] >= 2 ",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[0]]);
			break;
		case pattern_2_1:
			fprintf(stream,
					"forall y1,y2. \n y1 in [%s], y2 in [%s]\n y1=y2 \n  l[%s] >= 2 \n l[%s] >= 2  ",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[1]],
					name_of_dim[a->segments[0]], name_of_dim[a->segments[1]]);
			break;
		case pattern_3_1:
			fprintf(stream,
					"y1 in [%s], y2 in [%s], y3 in [%s]\n y1 = y2 y2 = y3\n "
							" l[%s] >= 2 \n l[%s] >= 2 \n l[%s] >= 2 \n",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[1]],
					name_of_dim[a->segments[2]], name_of_dim[a->segments[0]],
					name_of_dim[a->segments[1]], name_of_dim[a->segments[2]]);
			break;
		case pattern_1_lx: {
			e_seg = pr->PI[a->type].e_seg;
			if (e_seg == 1) {
				fprintf(stream,
						"forall y1. \n y1 in [%s]\n l[%s] >= 2 \n y1 = l[%s]",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]]);
			} else if (e_seg == 2)
				fprintf(stream,
						"forall y1. \n y1 in [%s] l[%s] >= 2 \n y1 = l[%s] + l[%s]",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]]);
			else
				fprintf(stream,
						"forall y1. \n y1 in [%s] l[%s] >= 2 \n y1 = l[%s] + l[%s] +l[%s]",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]],
						name_of_dim[a->segments[3]]);
			break;
		}
		case pattern_1_lx_1:
			fprintf(stream,
					"forall y1. \n y1 in [%s] \n l[%s] >= 2 \n y1 = l[%s]-1",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[0]],
					name_of_dim[a->segments[0]]);
			break;
		case pattern_1_l1:
			fprintf(stream, "forall y1. \n y1 in [%s] \n l[%s] >= 2 \n y1 = 1",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[0]]);
			break;
		case pattern_succ_1_2:
			fprintf(stream,
					"forall y1,y2. \n y1,y2 in [%s] \n y1 + 1 = y2 \n l[%s] >= 3",
					name_of_dim[a->segments[0]], name_of_dim[a->segments[0]]);
			break;
		case pattern_2_1_lx: {
			e_seg = pr->PI[a->type].e_seg;
			if (e_seg == 1) {
				fprintf(stream,
						"forall y1,y2. \n y1 in [%s], y2 in [%s], \nl[%s] >= 2\n  l[%s] >= 2 \n y1 = l[%s] + y2",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]]);
			} else if (e_seg == 2) {
				fprintf(stream,
						"forall y1,y2. \n y1 in [%s], y2 in [%s],\n l[%s] >= 2\n  l[%s] >= 2 \n y1 = l[%s] + l[%s] + y2",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]],
						name_of_dim[a->segments[3]]);
			} else {
				fprintf(stream,
						"forall y1,y2. \n y1 in [%s], y2 in [%s],\n l[%s] >= 2\n l[%s] >= 2 \n y1 = l[%s] + l[%s] + l[%s] + y2",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]],
						name_of_dim[a->segments[3]],
						name_of_dim[a->segments[4]]);
			}
			break;
		}
		case pattern_2_1_mlx: {
			e_seg = pr->PI[a->type].e_seg;
			if (e_seg == 1) {
				fprintf(stream,
						"forall y1,y2. \n y1 in [%s], y2 in [%s] \nl[%s] >= 2 \n l[%s] >= 2 \n y1 + l[%s] = y2",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]]);
			} else if (e_seg == 2) {
				fprintf(stream,
						"forall y1,y2. \n y1 in [%s], y2 in [%s]\n l[%s] >= 2 \n l[%s] >= 2 \n y1 + l[%s] + l[%s] = y2",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]],
						name_of_dim[a->segments[3]]);
			} else {
				fprintf(stream,
						"forall y1,y2. \n y1 in [%s], y2 in [%s]\n l[%s] >= 2 \n  l[%s] >= 2 \n y1 + l[%s] + l[%s] + l[%s] = y2",
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[0]],
						name_of_dim[a->segments[1]],
						name_of_dim[a->segments[2]],
						name_of_dim[a->segments[3]],
						name_of_dim[a->segments[4]]);
			}
			break;
		}
		default:
			fprintf(stream, "(pattern unknown)");
			break;
		}
	}
	fflush(stream);
}

void ucons_fprint(FILE * stream, ap_manager_t * man, ucons_t * a,
		char **name_of_dim) {
	ucons_internal_t *pr = ucons_init_from_manager(man, AP_FUNID_FPRINT, 0);
	size_t i;
	fprintf(stream,
			"\n\tsubgraph cluster_ucons_%zu {\n\tnode [shape=Mrecord] ;\n",
			ushape_number);
	if (!a)
		//fprintf (stream, "ucons EMPTY\n");
		fprintf(stream, "\tlabel = \"ucons %zu EMPTY\" ;\n }\n", ushape_number);
	else {
		if (!a->econs)
			fprintf(stream, "\tucons_dcons_%zu [label=\"bot\"] ;\n",
					ushape_number);
		else {
			char **dname, **name;
			size_t i, size;
			if (!name_of_dim)
				shape_init_name_of_dim(a->datadim, a->segmentdim);
			fprintf(stream,
					"\tlabel = \"ucons %zu of (datadim=%zu, segmdim=%zu)\" ;\n",
					ushape_number, a->datadim, a->segmentdim);
			//fprintf (stream, "ucons of [datadim=%zu, ptrdim=%zu] \n", a->datadim, a->segmentdim);
			// prepare names for the domain
			size = (a->datadim + 2 * a->segmentdim + 2 * 3); /* TODO: max 3 universally quantified vars */
			dname = (char **) malloc(size * sizeof(char *));
			for (i = 0; i < a->datadim; i++)
				dname[i] =
						(name_of_dim) ? name_of_dim[i] : shape_name_of_dim(i);

			name = (char **) malloc(size * sizeof(char *));
			for (i = 0; i < a->segmentdim; i++) {

				char *n =
						(name_of_dim) ?
								name_of_dim[a->datadim + i] :
								shape_name_of_dim(a->datadim + i);
				size_t lsize = (4 + strlen(n));
				name[a->datadim + i] = (char *) malloc(lsize * sizeof(char));
				snprintf(name[a->datadim + i], lsize, "%s", n);
				dname[a->datadim + i] = (char *) malloc(lsize * sizeof(char));
				snprintf(dname[a->datadim + i], lsize, "d(%s)", n);
				dname[a->datadim + a->segmentdim + i] = (char *) malloc(
						lsize * sizeof(char));
				snprintf(dname[a->datadim + a->segmentdim + i], lsize, "l[%s]",
						n);
			}
			for (i = 0; i < 6; i++) // TODO: max 6
				dname[a->datadim + 2 * a->segmentdim + i] = (char *) malloc(
						8 * sizeof(char));
			// print the existential constraints
			fprintf(stream,
					"\n\tsubgraph cluster_econs_%zu {\n\tnode [shape=Mrecord] ;\n",
					ushape_number);
			fprintf(stream, "\t ucons_econs_%zu [label=<<table><tr><td>econs: ",
					ushape_number);
			ap_abstract0_fprint(stream, pr->man_dcons, a->econs, dname);
			fprintf(stream, "</td></tr></table>> ] ;\n");
			fprintf(stream, "\t}\n");

			fprintf(stream,
					"\n\tsubgraph cluster_formulas_%zu{\n\tnode [shape=Mrecord] ;\n",
					ushape_number);
			// print the universal constraints
			pattern_t * r = a->udcons;
			size_t ii = 0;
			while (r != NULL) {
				//fprintf(stream, "\n\tsubgraph cluster_ucons %zu{\n\tnode [shape=Mrecord] ;\n", ii);
				ii += 1;

				fprintf(stream, "\tpattern_%zu%zu [label=<<table><tr><td> ", ii,
						ushape_number);
				dot_pattern_key_fprint(stream, pr, &r->key, &name[a->datadim]);
				//fprintf(stream, " ==> ");
				fprintf(stream, "</td></tr></table>> ] ;\n");
				fprintf(stream,
						"\tucons_%zu%zu [label=<<table><tr><td>ucons%zu%zu: ",
						ii, ushape_number, ii, ushape_number);
				if (r->dcons) {
					// prepare dname for this pattern
					for (i = 0; i < pr->PI[r->key.type].nr_y && i < 3; i++) // TODO: max 3
							{
						snprintf(dname[a->datadim + 2 * a->segmentdim + i], 8,
								"y%zu", (i + 1));
						snprintf(
								dname[a->datadim + 2 * a->segmentdim
										+ pr->PI[r->key.type].nr_y + i], 8,
								"d(y%zu)", (i + 1));
					}
					// print with the names
					ap_abstract0_fprint(stream, pr->man_dcons, r->dcons, dname);
				} else
					fprintf(stream, "true");
				r = r->hh.next;
				fprintf(stream, "</td></tr></table>> ] ;\n");

				fprintf(stream,
						"\tpattern_%zu%zu -> ucons_%zu%zu [label = implies ];\n",
						ii, ushape_number, ii, ushape_number);
				//fprintf(stream, "\t}\n");
			}
			fprintf(stream, "\t}\n");
			/*	
			 fprintf(stream, "\n\tsubgraph cluster_nodes_%zu {\n\tnode [shape=Mrecord] ;\n",ushape_number);
			 if(a->n2p!=NULL){
			 size_t i,j;
			 fprintf(stream, "\tucons_nodes_%zu [label=<<table><tr><td>node to patterns: \n ",
			 ushape_number);
			 //	printf("  \n\t node to pattern : \n");
			 for( i=0;i<a->segmentdim ; i++)
			 {
			 if (a->n2p[i].size>0) fprintf(stream, "\t n%zu:\n", i);
			 for( j=0;j<a->n2p[i].size ; j++){
			 dot_pattern_key_fprint(stream, pr, a->n2p[i].p[j], &name[a->datadim]);
			 fprintf(stream,"\n");
			 }
			 }
			 fprintf(stream, "</td></tr></table>> ] ;\n");

			 }
			 //fprintf(stream, "</td></tr></table>> ] ;\n");
			 for (i = 0; i < 6; i++) // TODO: max 6
			 free(dname[a->datadim + 2*a->segmentdim+i]);
			 free(dname);
			 fprintf(stream, "\t}\n");
			 */
		}
	} //end !a->econs
	fprintf(stream, "\t}\n");
}

void ucons_fprint_econs(FILE * stream, ap_manager_t * man, ucons_t * a,
		char **name_of_dim) {
	ucons_internal_t *pr = ucons_init_from_manager(man, AP_FUNID_FPRINT, 0);
	size_t i;
	fprintf(stream,
			"\n\tsubgraph cluster_ucons_%zu {\n\tnode [shape=Mrecord] ;\n",
			ushape_number);
	if (!a)
		//fprintf (stream, "ucons EMPTY\n");
		fprintf(stream, "\tlabel = \"ucons %zu EMPTY\" ;\n }\n", ushape_number);
	else {
		if (!a->econs)
			fprintf(stream, "\tucons_dcons_%zu [label=\"bot\"] ;\n",
					ushape_number);
		else {
			char **dname, **name;
			size_t i, size;
			if (!name_of_dim)
				shape_init_name_of_dim(a->datadim, a->segmentdim);
			fprintf(stream,
					"\tlabel = \"ucons %zu of (datadim=%zu, segmdim=%zu)\" ;\n",
					ushape_number, a->datadim, a->segmentdim);
			//fprintf (stream, "ucons of [datadim=%zu, ptrdim=%zu] \n", a->datadim, a->segmentdim);
			// prepare names for the domain
			size = (a->datadim + 2 * a->segmentdim + 2 * 3); /* TODO: max 3 universally quantified vars */
			dname = (char **) malloc(size * sizeof(char *));
			for (i = 0; i < a->datadim; i++)
				dname[i] =
						(name_of_dim) ? name_of_dim[i] : shape_name_of_dim(i);

			name = (char **) malloc(size * sizeof(char *));
			for (i = 0; i < a->segmentdim; i++) {

				char *n =
						(name_of_dim) ?
								name_of_dim[a->datadim + i] :
								shape_name_of_dim(a->datadim + i);
				size_t lsize = (4 + strlen(n));
				name[a->datadim + i] = (char *) malloc(lsize * sizeof(char));
				snprintf(name[a->datadim + i], lsize, "%s", n);
				dname[a->datadim + i] = (char *) malloc(lsize * sizeof(char));
				snprintf(dname[a->datadim + i], lsize, "d(%s)", n);
				dname[a->datadim + a->segmentdim + i] = (char *) malloc(
						lsize * sizeof(char));
				snprintf(dname[a->datadim + a->segmentdim + i], lsize, "l[%s]",
						n);
			}
			for (i = 0; i < 6; i++) // TODO: max 6
				dname[a->datadim + 2 * a->segmentdim + i] = (char *) malloc(
						8 * sizeof(char));
			// print the existential constraints
			fprintf(stream,
					"\n\tsubgraph cluster_econs_%zu {\n\tnode [shape=Mrecord] ;\n",
					ushape_number);
			fprintf(stream, "\t ucons_econs_%zu [label=<<table><tr><td>econs: ",
					ushape_number);
			ap_abstract0_fprint(stream, pr->man_dcons, a->econs, dname);
			fprintf(stream, "</td></tr></table>> ] ;\n");
			fprintf(stream, "\t}\n");

		}
	} //end !a->econs
	fprintf(stream, "\t}\n");
}

void ucons_fprint_dcons(FILE * stream, ap_manager_t * man, ucons_t * a,
		char **name_of_dim, pattern_key_t *key) {
	ucons_internal_t *pr = ucons_init_from_manager(man, AP_FUNID_FPRINT, 0);
	size_t i;
	fprintf(stream,
			"\n\tsubgraph cluster_ucons_%zu {\n\tnode [shape=Mrecord] ;\n",
			ushape_number);
	if (!a)
		//fprintf (stream, "ucons EMPTY\n");
		fprintf(stream, "\tlabel = \"ucons %zu EMPTY\" ;\n }\n", ushape_number);
	else {
		if (!a->econs)
			fprintf(stream, "\tucons_dcons_%zu [label=\"bot\"] ;\n",
					ushape_number);
		else {
			char **dname, **name;
			size_t i, size;
			if (!name_of_dim)
				shape_init_name_of_dim(a->datadim, a->segmentdim);
			fprintf(stream,
					"\tlabel = \"ucons %zu of (datadim=%zu, segmdim=%zu)\" ;\n",
					ushape_number, a->datadim, a->segmentdim);
			//fprintf (stream, "ucons of [datadim=%zu, ptrdim=%zu] \n", a->datadim, a->segmentdim);
			// prepare names for the domain
			size = (a->datadim + 2 * a->segmentdim + 2 * 3); /* TODO: max 3 universally quantified vars */
			dname = (char **) malloc(size * sizeof(char *));
			for (i = 0; i < a->datadim; i++)
				dname[i] =
						(name_of_dim) ? name_of_dim[i] : shape_name_of_dim(i);

			name = (char **) malloc(size * sizeof(char *));
			for (i = 0; i < a->segmentdim; i++) {

				char *n =
						(name_of_dim) ?
								name_of_dim[a->datadim + i] :
								shape_name_of_dim(a->datadim + i);
				size_t lsize = (4 + strlen(n));
				name[a->datadim + i] = (char *) malloc(lsize * sizeof(char));
				snprintf(name[a->datadim + i], lsize, "%s", n);
				dname[a->datadim + i] = (char *) malloc(lsize * sizeof(char));
				snprintf(dname[a->datadim + i], lsize, "d(%s)", n);
				dname[a->datadim + a->segmentdim + i] = (char *) malloc(
						lsize * sizeof(char));
				snprintf(dname[a->datadim + a->segmentdim + i], lsize, "l[%s]",
						n);
			}
			for (i = 0; i < 6; i++) // TODO: max 6
				dname[a->datadim + 2 * a->segmentdim + i] = (char *) malloc(
						8 * sizeof(char));

			fprintf(stream,
					"\n\tsubgraph cluster_formulas_%zu{\n\tnode [shape=Mrecord] ;\n",
					ushape_number);
			// print the universal constraints
			pattern_t * r = NULL;

			unsigned keylen = (pr->PI[key->type].u_seg) * sizeof(size_t)
					+ sizeof(pattern_key_t);
			HASH_FIND(hh, a->udcons, key, keylen, r);
			size_t ii = 0;
			if (r != NULL) {
				//fprintf(stream, "\n\tsubgraph cluster_ucons %zu{\n\tnode [shape=Mrecord] ;\n", ii);
				ii += 1;

				fprintf(stream, "\tpattern_%zu%zu [label=<<table><tr><td> ", ii,
						ushape_number);
				dot_pattern_key_fprint(stream, pr, &r->key, &name[a->datadim]);
				//fprintf(stream, " ==> ");
				fprintf(stream, "</td></tr></table>> ] ;\n");
				fprintf(stream,
						"\tucons_%zu%zu [label=<<table><tr><td>ucons%zu%zu: ",
						ii, ushape_number, ii, ushape_number);
				if (r->dcons) {
					// prepare dname for this pattern
					for (i = 0; i < pr->PI[r->key.type].nr_y && i < 3; i++) // TODO: max 3
							{
						snprintf(dname[a->datadim + 2 * a->segmentdim + i], 8,
								"y%zu", (i + 1));
						snprintf(
								dname[a->datadim + 2 * a->segmentdim
										+ pr->PI[r->key.type].nr_y + i], 8,
								"d(y%zu)", (i + 1));
					}
					// print with the names
					ap_abstract0_fprint(stream, pr->man_dcons, r->dcons, dname);
				} else
					fprintf(stream, "true");
				r = r->hh.next;
				fprintf(stream, "</td></tr></table>> ] ;\n");

				fprintf(stream,
						"\tpattern_%zu%zu -> ucons_%zu%zu [label = implies ];\n",
						ii, ushape_number, ii, ushape_number);
				//fprintf(stream, "\t}\n");
			}
			fprintf(stream, "\t}\n");

		}
	} //end !a->econs
	fprintf(stream, "\t}\n");
}

void ucons_fprintdiff(FILE * stream, ap_manager_t * man, ucons_t * a1,
		ucons_t * a2, char **name_of_dim) {
	ucons_internal_t *pr = ucons_init_from_manager(man, AP_FUNID_FPRINTDIFF, 0);
	fprintf(stream, "ucons1 (size %zu) and ucons2 (size %zu)\n", a1->size,
			a2->size);
	/* TODO: priority 1 */
}

void ucons_fdump(FILE * stream, ap_manager_t * man, ucons_t * a) {
	ucons_fprint(stream, man, a, NULL);
}

/* ============================================================ */
/* Serialization */
/* ============================================================ */

/* TODO: priority 0 */
/* NOT IMPLEMENTED: do nothing */
ap_membuf_t ucons_serialize_raw(ap_manager_t * man, ucons_t * a) {
	ucons_internal_t *pr = ucons_init_from_manager(man, AP_FUNID_SERIALIZE_RAW,
			0);
	ap_membuf_t buf;
	buf.size = 0;
	buf.ptr = NULL;
	ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return buf;
}

/* Builds a value from a code!
 * The ptr is an array of 3 size_t [datadim, segmdim, codeop]
 */
ucons_t *
ucons_deserialize_raw(ap_manager_t * man, void *ptr, size_t * size) {
	ucons_internal_t *pr = ucons_init_from_manager(man,
			AP_FUNID_DESERIALIZE_RAW, 0);

	size_t *ucons_raw = (size_t *) ptr;
	size_t datadim = ucons_raw[0];
	size_t segmdim = ucons_raw[1];
	size_t code = ucons_raw[2];
	ap_dim_t l = 0;
	size_t nodex = 1;
	size_t nodey = 2;
	size_t nodez = 3;
	ucons_t *r = ucons_top(man, datadim, segmdim);

	// for initial pointers set to null
	if(code == 5) return r;

	/* common constraints l[x]==_l and _l>=1 */
	if (code != 1) {
		assert(segmdim >= 2 && datadim >= 1);
		ap_lincons0_array_t arr = ap_lincons0_array_make(2);
		arr.p[0] = shape_lincons_x_y_v_cst(AP_CONS_EQ, OFFSET_LEN, 1,
				datadim + nodex, 0, 0, -1, l, 0, datadim, segmdim);
		arr.p[1] = shape_lincons_x_y_v_cst(AP_CONS_SUPEQ, OFFSET_LEN, 1,
				datadim + nodex, 0, 0, 0, 0, -1, datadim, segmdim);
		r = ucons_meet_lincons_array(man, true, r, &arr);
		ap_lincons0_array_clear(&arr);
	}
	/* common constraint l[x]==l[y] and l[y]>=1 */
	if (code == 2 || code == 4) {
		assert(segmdim >= 3 && datadim >= 1);
		ap_lincons0_array_t arr = ap_lincons0_array_make(2);
		arr.p[0] = shape_lincons_x_y_v_cst(AP_CONS_EQ, OFFSET_LEN, 1,
				datadim + nodex, -1, datadim + nodey, 0, 0, 0, datadim,
				segmdim);
		arr.p[1] = shape_lincons_x_y_v_cst(AP_CONS_SUPEQ, OFFSET_LEN, 1,
				datadim + nodey, 0, 0, 0, 0, -1, datadim, segmdim);
		//      arr.p[2] =
		//              shape_lincons_x_y_v_cst (AP_CONS_EQMOD, var_ptr_null, 1,
		//                                       nodey, -1, nodex, 0, 0, 0, datadim,
		//                                       segmdim);
		r = ucons_meet_lincons_array(man, true, r, &arr);
		ap_lincons0_array_clear(&arr);
	}
	/* specific constraints */
	switch (code) {
	case 1: /* l[x]+l[y]==_l and l[x]>=1 and l[y]>=1 */
	{
		assert(segmdim >= 2 && datadim >= 1);
		ap_lincons0_array_t arr = ap_lincons0_array_make(3);
		arr.p[0] = shape_lincons_x_y_v_cst(AP_CONS_EQ, OFFSET_LEN, 1,
				datadim + nodex, 1, datadim + nodey, -1, l, 0, datadim,
				segmdim);
		arr.p[1] = shape_lincons_x_y_v_cst(AP_CONS_SUPEQ, OFFSET_LEN, 1,
				datadim + nodex, 0, 0, 0, 0, -1, datadim, segmdim);
		arr.p[2] = shape_lincons_x_y_v_cst(AP_CONS_SUPEQ, OFFSET_LEN, 1,
				datadim + nodey, 0, 0, 0, 0, -1, datadim, segmdim);
		r = ucons_meet_lincons_array(man, true, r, &arr);
		ap_lincons0_array_clear(&arr);
		break;
	}
	case 3: /* non-equal length lists: l[y]+1<=_l and l[y]>=1 */
	{
		assert(segmdim >= 3 && datadim >= 1);
		ap_lincons0_array_t arr = ap_lincons0_array_make(2);
		arr.p[0] = shape_lincons_x_y_v_cst(AP_CONS_SUPEQ, OFFSET_LEN, -1,
				datadim + nodey, 0, 0, 1, l, -1, datadim, segmdim);
		arr.p[1] = shape_lincons_x_y_v_cst(AP_CONS_SUPEQ, OFFSET_LEN, 1,
				datadim + nodey, 0, 0, 0, 0, -1, datadim, segmdim);
		r = ucons_meet_lincons_array(man, true, r, &arr);
		ap_lincons0_array_clear(&arr);
		break;
	}
	case 4: /* equal length lists: l[z]=l[x] */
	{
		assert(segmdim >= 4);
		ap_lincons0_array_t arr = ap_lincons0_array_make(1);
		arr.p[0] = shape_lincons_x_y_v_cst(AP_CONS_EQ, OFFSET_LEN, -1,
				datadim + nodex, 1, datadim + nodez, 0, 0, 0, datadim, segmdim);
		r = ucons_meet_lincons_array(man, true, r, &arr);
		ap_lincons0_array_clear(&arr);
		break;
	}
	default:
		break;
	}
	return r;
}
