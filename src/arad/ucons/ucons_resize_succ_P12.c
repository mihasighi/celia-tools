/**************************************************************************/
/*                                                                        */
/*  CINV Library / Shape Domain                                           */
/*                                                                        */
/*  Copyright (C) 2009-2011                                               */
/*    LIAFA (University of Paris Diderot and CNRS)                        */
/*                                                                        */
/*                                                                        */
/*  you can redistribute it and/or modify it under the terms of the GNU   */
/*  Lesser General Public License as published by the Free Software       */
/*  Foundation, version 3.                                                */
/*                                                                        */
/*  It is distributed in the hope that it will be useful,                 */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/*  GNU Lesser General Public License for more details.                   */
/*                                                                        */
/*  See the GNU Lesser General Public License version 3.                  */
/*  for more details (enclosed in the file LICENSE).                      */
/*                                                                        */
/**************************************************************************/


/*
 * ucons_resize_succ_P12.c
 *
 * CINV Library / Shape Domain    
 *  Created on: Jan 20, 2011
 *      Author: cezara
 */
#include "uthash.h"
#include "ucons.h"
#include "ucons_fun.h"
#include "ucons_internal.h"
#include "shape_macros.h"
#include "apron2shape.h"


ucons_t * split_with_pattern_succ_P12(ucons_internal_t* pr, ucons_t *r,
		pattern_key_t *pattern_j, ap_dim_t n1, ap_dim_t n2){


#ifndef NDEBUG
	fprintf(stdout,"  \n\t split_with_pattern_succ_P12 on : \n");
	fprintf(stdout," \n n1 =%zu n2 = %zu \n ",n1 ,n2);
	fprintf(stdout,"with pattern->type = %zu \n", pattern_j->type);
	ucons_fprint(stdout,pr->man,r,NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

	//the master pattern handles the auxiliary ones
	//they will be eliminated by y1 + 1 = y2
	if(pr->PI[pattern_j->type].kind == pattern_1_l1 ||
			pr->PI[pattern_j->type].kind == pattern_1_lx_1)
		return r;

	pattern_key_t * pat_trasf;
	unsigned keylen;
	size_t u_seg, e_seg, nr_y;
	size_t y1, dy1, y2, dy2, dn1, dn2;

	dn2 = r->datadim + n2;
	dn1 = r->datadim + n1;
	y1 = r->datadim + 2*r->segmentdim;
	y2 = r->datadim + 2*r->segmentdim + 1;
	dy1 = r->datadim + 2*r->segmentdim + 2;
	dy2 = r->datadim + 2*r->segmentdim + 3;

	u_seg = pr->PI[pattern_j->type].u_seg;
	e_seg = pr->PI[pattern_j->type].e_seg;
	nr_y = pr->PI[pattern_j->type].nr_y;

	keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

	pattern_t * jdcons = NULL; // the value in the hash table associated with pattern_j
	HASH_FIND(hh, r->udcons, pattern_j, keylen, jdcons);
	if(jdcons){
		//		 there is a constraint attached to this node
		if(jdcons->dcons!=NULL){
			/* test pattern satisfiability */
			if(!test_l_leq_v(pr->man_dcons, r->econs,
					r->datadim, r->segmentdim, pattern_j->segments[0], 2)){
				//pattern has at least one instance


				/*	step 1.
						build y1 = 1 ==> ap_jdcons_y1
				 */
				ap_abstract0_t * ap_jdcons_y1 = ap_abstract0_copy(pr->man_dcons, jdcons->dcons);

				ap_lincons0_array_t arr;
				arr = ap_lincons0_array_make (3);
				/*
				 *  y1-1 == 0
				 */
				arr.p[0].constyp = AP_CONS_EQ;
				arr.p[0].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
				arr.p[0].scalar = NULL;
				ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, y1, 1);
				ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);

				/* dy1 - dn2 == 0
				 */
				arr.p[1].constyp = AP_CONS_EQ;
				arr.p[1].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
				arr.p[1].scalar = NULL;
				ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy1, -1);
				ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dn2, 1);
				ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);
				/*
				 *  y2-2 == 0
				 */
				arr.p[2].constyp = AP_CONS_EQ;
				arr.p[2].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
				arr.p[2].scalar = NULL;
				ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0, y2, 1);
				ap_linexpr0_set_cst_scalar_int (arr.p[2].linexpr0, -2);

				ap_jdcons_y1 = ap_abstract0_meet_lincons_array (pr->man_dcons,
						true, ap_jdcons_y1, &arr);

				ap_lincons0_array_clear (&arr);

				// eliminate y1, d(y1)
				ap_dimchange_t dimrm;
				dimrm.realdim = 0;
				dimrm.intdim = 2;
				dimrm.dim = (ap_dim_t *) malloc (dimrm.intdim * sizeof (ap_dim_t));
				dimrm.dim[0] = y1 ;
				dimrm.dim[1] = dy1 ;

				ap_jdcons_y1=
						ap_abstract0_remove_dimensions (pr->man_dcons, true, ap_jdcons_y1, &dimrm);

				ap_dimchange_clear (&dimrm);
				// substitute y1 with y1+1, so y1 becomes one
				/* in ap_jdcons_y1 y1 := y1 + 1*/

				ap_linexpr0_t *expr;
				expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
				ap_linexpr0_set_coeff_scalar_int (expr, y1, 1);
				ap_linexpr0_set_cst_scalar_int (expr, 1);

				ap_jdcons_y1 = ap_abstract0_substitute_linexpr( pr->man_dcons, true,
						ap_jdcons_y1, y1, expr, NULL );

				ap_linexpr0_free(expr);

				pattern_key_t *look = NULL;
				checked_malloc(look, pattern_key_t, 1,
						sizeof(pattern_key_t) + u_seg*sizeof(size_t), return NULL;);
				look->type = get_pattern_type(pr, 1, 0, 1, pattern_1_l1);
				//look for y1 =1
				look->segments[0] = n2;

				pattern_t * dcons_y1 = NULL;

				HASH_FIND(hh,r->udcons, look, keylen, dcons_y1);
				if(!dcons_y1){
					checked_malloc(dcons_y1,pattern_t,1,
							sizeof(pattern_t)+(u_seg)*sizeof(size_t),return NULL;);
					memset(dcons_y1, 0, sizeof(pattern_t)+(u_seg)*sizeof(size_t));
					dcons_y1->key.type = look->type;
					//for (size_t i=0 ; i<(u_seg); i++)
					dcons_y1->key.segments[0] = n2;
					dcons_y1->dcons = ap_jdcons_y1;
					HASH_ADD(hh,r->udcons,key,keylen,dcons_y1);
					r = add_pattern_n2p(pr, r, look);
				}
				else{
					if(dcons_y1->dcons != NULL)
						ap_abstract0_free(pr->man_dcons, dcons_y1->dcons);
					dcons_y1->dcons = ap_jdcons_y1;
				}
				free(look);
				look = NULL;

				//step 2 add \forall y \in n2 . y = l[n2] - 1 ==> U
				//copy U from \forall y \in n1. y = l[n1] - 1 ==> U
// and  y = l[n2] -1 meet y = 1 if len(n1) <= 3
				checked_malloc(look, pattern_key_t, 1,
						sizeof(pattern_key_t) + u_seg*sizeof(size_t), return NULL;);
				look->type = get_pattern_type(pr, 1, 0, 1, pattern_1_lx_1);
				//look for y = l[n2] - 1
				look->segments[0] = n1;

				pattern_t * n1_dcons_ylx = NULL;
				pattern_t * n2_dcons_ylx = NULL;

				HASH_FIND(hh,r->udcons, look, keylen, n1_dcons_ylx);
				if(n1_dcons_ylx && n1_dcons_ylx->dcons ){

					look->segments[0] = n2;
					HASH_FIND(hh,r->udcons, look, keylen, n2_dcons_ylx);
					if(!n2_dcons_ylx){
						checked_malloc(n2_dcons_ylx,pattern_t,1,
								sizeof(pattern_t)+(u_seg)*sizeof(size_t),return NULL;);
						memset(n2_dcons_ylx, 0, sizeof(pattern_t)+(u_seg)*sizeof(size_t));
						n2_dcons_ylx->key.type = look->type;
						//for (size_t i=0 ; i<(u_seg); i++)
						n2_dcons_ylx->key.segments[0] = n2;
						n2_dcons_ylx->dcons = ap_abstract0_copy(pr->man_dcons,
								n1_dcons_ylx->dcons);
						HASH_ADD(hh,r->udcons,key,keylen,n2_dcons_ylx);
						r = add_pattern_n2p(pr, r, look);
					}
					else{
						if(n2_dcons_ylx->dcons != NULL)
							ap_abstract0_free(pr->man_dcons, n2_dcons_ylx->dcons);
						n2_dcons_ylx->dcons = ap_abstract0_copy(pr->man_dcons,
								n1_dcons_ylx->dcons);
					}
				}
				free(look);
				look = NULL;

				// step 3. transport y1 + 1= y2 on n2 if len(n1) > 3
				// and  y = l[n2] -1 meet y = 1 if len(n1) <= 3
				if(test_l_leq_v(pr->man_dcons, r->econs,
						r->datadim, r->segmentdim, pattern_j->segments[0], 3)){
					//len is less or equal then 3
					// pattern has only one instance
					//case 1

					checked_malloc(look, pattern_key_t, 1,
							sizeof(pattern_key_t) + u_seg*sizeof(size_t), return NULL;);
					look->type = get_pattern_type(pr, 1, 0, 1, pattern_1_l1);
					//look for y1 =1 ==> a1 and y = l[n2] -1 ==> a2
					// meet a1 with a2
					look->segments[0] = n2;
					pattern_t * a1 = NULL; // y1 = 1
					pattern_t * a2 = NULL; // y1 = l[n2] - 1

					HASH_FIND(hh,r->udcons, look, keylen, a1);

					look->type = get_pattern_type(pr, 1, 0, 1, pattern_1_lx_1);
					HASH_FIND(hh,r->udcons, look, keylen, a2);

					if((a1!=NULL) && (a2!=NULL)){
						ap_abstract0_t * aux = ap_abstract0_meet(pr->man_dcons,
								true, a1->dcons, a2->dcons);
						a1->dcons = ap_abstract0_copy(pr->man_dcons, aux);
						a2->dcons = ap_abstract0_copy(pr->man_dcons, aux);
						ap_abstract0_free(pr->man_dcons, aux);
					}
					else{
						if(a1 != NULL){

							checked_malloc(a2,pattern_t,1,
									sizeof(pattern_t)+(u_seg)*sizeof(size_t),return NULL;);
							memset(a2, 0, sizeof(pattern_t)+(u_seg)*sizeof(size_t));
							a2->key.type = get_pattern_type(pr, 1, 0, 1, pattern_1_lx_1);
							a2->key.segments[0] = n2;
							a2->dcons = ap_abstract0_copy(pr->man_dcons,
									a1->dcons);
							HASH_ADD(hh,r->udcons,key,keylen,a2);
							r = add_pattern_n2p(pr, r, &a2->key);
						}
						if(a2 != NULL){
							checked_malloc(a1,pattern_t,1,
									sizeof(pattern_t)+(u_seg)*sizeof(size_t),return NULL;);
							memset(a1, 0, sizeof(pattern_t)+(u_seg)*sizeof(size_t));
							a1->key.type = get_pattern_type(pr, 1, 0, 1, pattern_1_l1);;
							a1->key.segments[0] = n2;
							a1->dcons = ap_abstract0_copy(pr->man_dcons,
									a2->dcons);
							HASH_ADD(hh, r->udcons, key, keylen, a1);
							r = add_pattern_n2p(pr, r, &a1->key);
						}
					}

				}//end len(n1) <= 3
				else{
					//len might be strictly greater then 3
					// pattern has more then one instance
					//case 2
					//step3  transfer on n2  pattern y1 + 1 = y2 == > U1

					pattern_key_t *look2 = NULL;
					pattern_t *a2 = NULL;
					checked_malloc(look2, pattern_key_t, 1,
							sizeof(pattern_key_t) + u_seg*sizeof(size_t), return NULL;);
					look2->type = pattern_j->type;
					look2->segments[0] = n2;
					keylen = u_seg*sizeof(size_t) + sizeof(pattern_key_t);

					HASH_FIND(hh,r->udcons, look2, keylen, a2);
					if(!a2){
						checked_malloc(a2,pattern_t,1,sizeof(pattern_t)+(u_seg)*sizeof(size_t),return NULL;);
						memset(a2, 0, sizeof(pattern_t)+(u_seg)*sizeof(size_t));
						a2->key.type = pattern_j->type;
						a2->key.segments[0] = n2;
						a2->dcons = NULL;
						HASH_ADD(hh,r->udcons,key,keylen,a2);
						r=add_pattern_n2p(pr,r,look2);
					}
					HASH_FIND(hh,r->udcons, look2, keylen, a2);

					ap_dim_t * ydim = NULL;
					ap_linexpr0_t **yexpr = NULL;

					checked_malloc(ydim, ap_dim_t, nr_y,sizeof(ap_dim_t), return NULL;);
					checked_malloc(yexpr, ap_linexpr0_t*, nr_y,sizeof(ap_linexpr0_t*), return NULL;);


					a2->dcons = ap_abstract0_copy (pr->man_dcons, jdcons->dcons);
					for(size_t i = 0 ; i < nr_y; i++){
						ydim[i] = r->datadim + 2 * r->segmentdim + i;

						yexpr[i] = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2 * nr_y );
						ap_linexpr0_set_coeff_scalar_int (yexpr[i], ydim[i], 1);
						ap_linexpr0_set_cst_scalar_int (yexpr[i], 1);
						a2->dcons= ap_abstract0_substitute_linexpr( pr->man_dcons,
								true, a2->dcons, ydim[i], yexpr[i], NULL );
					}

					for(size_t i = 0 ; i < nr_y; i++){
						ap_linexpr0_free(yexpr[i]);
					}
					free(ydim);
					free(yexpr);

					ap_lincons0_array_t arr;
					arr = ap_lincons0_array_make (2);
					/*
					 *  y1-1 >= 0
					 */
					arr.p[0].constyp = AP_CONS_SUPEQ;
					arr.p[0].linexpr0 =
							ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim +  2 * nr_y);
					arr.p[0].scalar = NULL;
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, y1, 1);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);
					/*
					 *  l[n1] - 3 - y1 >=0
					 */
					arr.p[1].constyp = AP_CONS_SUPEQ;
					arr.p[1].linexpr0 =
							ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim +  2 * nr_y);
					arr.p[1].scalar = NULL;
					ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, r->datadim + r->segmentdim + n1, 1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, y1, -1);
					ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, -3);

					a2->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons,
							true, a2->dcons, &arr);

					ap_lincons0_array_clear (&arr);


				}
			}//end pattern satisfiable
		}//end jdcons->dcons ! = NULL
		else{
			fprintf(stdout,"\n pattern with null numerical constraint found in hash table ");
			fflush(stdout);
			ap_manager_raise_exception (pr->man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
					"not implemented");
		}
		HASH_DEL(r->udcons,jdcons);
		remove_pattern_n2p(pr,r,&jdcons->key);
		free(jdcons);
	}

	pattern_key_t * loox = NULL;
	pattern_t * aux = NULL;

	checked_malloc(loox, pattern_key_t, 1,
			sizeof(pattern_key_t) + u_seg*sizeof(size_t), return NULL;);
	loox->type = get_pattern_type(pr, 1, 0, 1, pattern_1_lx_1);
	loox->segments[0] = n1;
	keylen = u_seg*sizeof(size_t) + sizeof(pattern_key_t);

	HASH_FIND(hh,r->udcons, loox, keylen, aux);

	if(aux){
		HASH_DEL(r->udcons,aux);
		remove_pattern_n2p(pr,r,&aux->key);
		free(aux);
	}

	loox->type = get_pattern_type(pr, 1, 0, 1, pattern_1_l1);
	aux  = NULL;
	HASH_FIND(hh,r->udcons, loox, keylen, aux);
	if(aux){
		HASH_DEL(r->udcons,aux);
		remove_pattern_n2p(pr,r,&aux->key);
		free(aux);
	}

#ifndef NDEBUG
	fprintf(stdout,"  \n\t split_with_pattern_succ_P12: \n");
	fprintf(stdout,"n1 =%zu n2 = %zu ",n1 ,n2);
	fprintf(stdout,"\n");
	ucons_fprint(stdout,pr->man,r,NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif
	return r;

}
ucons_t * fold_with_closure_succ_P12( ucons_internal_t *pr,
		ucons_t * a, ap_dim_t * tdim, size_t size,bool update_lenght){


	size_t i,size_node_tdim;


	ap_dim_t *node_tdim; // the nodes to concatenate from dimensions in tdim

#ifndef NDEBUG
	fprintf(stdout,"  \n\t fold_with_pattern_succ_P12 \n tdim :");
	for(size_t q=0 ; q<size; q++)
		fprintf(stdout," %zu ",tdim[q]);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

	size_node_tdim = 0;
	node_tdim = NULL;
	for(i=0; i<size; i++){
		if(tdim[i] >= a->datadim)
			size_node_tdim++;
		checked_realloc(node_tdim,ap_dim_t,sizeof(ap_dim_t),size_node_tdim,return NULL;);
		node_tdim[size_node_tdim-1] = tdim[i] - a->datadim;
	}


	ucons_t *r = ucons_copy_internal (pr, a);

#ifndef NDEBUG
	fprintf(stdout,"  \n\t fold_with_pattern_succ_P12 \n tdim :");
	for(size_t q= 0 ; q<size_node_tdim; q++)
		fprintf(stdout," %zu ",node_tdim[q]);
	fprintf(stdout,"\n on \n");
	ucons_fprint(stdout,pr->man,a,NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif


	ap_abstract0_t *dcons = NULL; //the new value with y1 + 1 = y2;
	ap_abstract0_t *dcons_y1 = NULL; //the new value with y1 =1
	ap_abstract0_t *dcons_ylx = NULL; // the new value with y2 = l[tdim[0]] - 1

	ap_dim_t y1, y2, dy1, dy2;
	y1 = r->datadim + 2*r->segmentdim;
	y2 = r->datadim + 2*r->segmentdim + 1;
	dy1 = r->datadim + 2*r->segmentdim + 2;
	dy2 = r->datadim + 2*r->segmentdim + 3;

	bool sgl = test_singleton(pr->man_dcons, r->econs, r->datadim,
			r->segmentdim, node_tdim[0]);

	pattern_key_t * look;
	checked_malloc(look, pattern_key_t, 1,
			sizeof(pattern_key_t) + 1*sizeof(size_t), return NULL;);
	look->type = get_pattern_type(pr, 1, 0, 2, pattern_succ_1_2);
	look->segments[0] = node_tdim[0];
	unsigned keylen = 1*sizeof(size_t) + sizeof(pattern_key_t);

	pattern_t* aux = NULL; //the old value with succ_P12in tdim[0]

	HASH_FIND(hh,r->udcons, look, keylen, aux);
	if(!sgl && (aux == NULL || aux->dcons == NULL)){
		if (aux!=NULL){
			HASH_DEL(r->udcons,aux);
			remove_pattern_n2p(pr,r,&aux->key);
			free(aux);
		}
		if(update_lenght) update_lenghts(pr, r, node_tdim, size_node_tdim );
#ifndef NDEBUG
		fprintf(stdout,"  \n\t fold_with_pattern_succ_P12: \n");
		fprintf(stdout," returns with pattern info insufficient \n");
		fflush(stdout);
#endif
		//TODO update the other patterns
		return r;
	}
	//define dcons_y1

	if(sgl){
		ap_dimchange_t dimadd;
		ap_dimchange_init (&dimadd, 2, 0);
		dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
		dimadd.dim[0] = r->datadim + 2*r->segmentdim;
		dimadd.dim[1] = r->datadim + 2*r->segmentdim;

		dcons_y1 = ap_abstract0_add_dimensions(pr->man_dcons,false,r->econs,&dimadd,false);

#ifndef NDEBUG
	fprintf(stdout,"  \n\t r->econs : \n");
	ap_abstract0_fprint(stdout,pr->man_dcons,r->econs,NULL);
	fprintf(stdout,"\n add 2 dims is \n ");
	ap_abstract0_fprint(stdout,pr->man_dcons,dcons_y1,NULL);
	fflush(stdout);
#endif

		ap_dimchange_clear (&dimadd);

		ap_lincons0_array_t arr;
		arr = ap_lincons0_array_make (2);
		/*
		 *  y1-1 == 0
		 */
		arr.p[0].constyp = AP_CONS_EQ;
		arr.p[0].linexpr0 =
				ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
		arr.p[0].scalar = NULL;
		ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, y1, 1);
		ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);

		/* dy1 - d[node_tdim[1]] == 0
		 */
		arr.p[1].constyp = AP_CONS_EQ;
		arr.p[1].linexpr0 =
				ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
		arr.p[1].scalar = NULL;
		dy1 = r->datadim + 2*r->segmentdim + 1;
		ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy1, 1);
		ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, r->datadim + node_tdim[1], -1);
		ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);

		dcons_y1 = ap_abstract0_meet_lincons_array (pr->man_dcons,
				true, dcons_y1, &arr);

		ap_lincons0_array_clear (&arr);

#ifndef NDEBUG
	fprintf(stdout,"\n dcons_y1 is from 1  \n ");
	ap_abstract0_fprint(stdout,pr->man_dcons,dcons_y1,NULL);
	fprintf(stdout," \n");
	fflush(stdout);
#endif

	}
	else{
		//search the old value with y1 = 1
		pattern_key_t *look1 = NULL;
		checked_malloc(look1, pattern_key_t, 1,
				sizeof(pattern_key_t) + 1*sizeof(size_t), return NULL;);
		look1->type = get_pattern_type(pr, 1, 0, 1, pattern_1_l1);
		look1->segments[0] = node_tdim[0];
		keylen = 1*sizeof(size_t) + sizeof(pattern_key_t);

		pattern_t *aux1;
		HASH_FIND(hh,r->udcons, look1, keylen, aux1);

		if((aux1 != NULL) && (aux1->dcons != NULL)){
			dcons_y1 = ap_abstract0_copy (pr->man_dcons, aux1->dcons);
		}

#ifndef NDEBUG
	fprintf(stdout,"\n dcons_y1 is \n ");
	if(dcons_y1) ap_abstract0_fprint(stdout,pr->man_dcons,dcons_y1,NULL);
	fprintf(stdout," \n");
	fflush(stdout);
#endif
	} //end definition of dcons_y1

	ap_abstract0_t *aux_dcons = NULL;

	for(size_t i = 1 ; i < size; i++){
		//concatenate tdim[i]
		bool sgli_1 = test_singleton(pr->man_dcons, r->econs,
				r->datadim, r->segmentdim, node_tdim[i-1]);
		bool sgli = test_singleton(pr->man_dcons, r->econs,
				r->datadim, r->segmentdim, node_tdim[i]);
		//Connect node_tdim[i] with his left neighbor
		if( (i>1) && sgli_1){
			ap_abstract0_t *aux_dcons = NULL;
			//(y1,y2) == (node_tdim[i-1], node_tdim[i])
			aux_dcons = create1_succ_P12(pr, r, node_tdim, size_node_tdim, i);
			if(dcons !=NULL )
				dcons=ap_abstract0_join(pr->man_dcons,true, dcons, aux_dcons);
			else {
				dcons = ap_abstract0_copy(pr->man_dcons, aux_dcons);
				ap_abstract0_free(pr->man_dcons, aux_dcons);
			}
#ifndef NDEBUG
	fprintf(stdout,"  \n\t dcons : \n");
	ap_abstract0_fprint(stdout,pr->man_dcons,dcons,NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

		}
		if(!sgli_1){
			ap_abstract0_t *aux_dcons = NULL;
			//(y1,y2) == (last(node_tdim[i-1]), node_tdim[i])
			aux_dcons = create2_succ_P12(pr, r, node_tdim, size_node_tdim, i);

			if(aux_dcons == NULL){
				//MS: produces core dump
				//dcons = ap_abstract0_top(pr->man_dcons, r->datadim +2* r->segmentdim + 4,0);
				aux_dcons = ap_abstract0_top(pr->man_dcons, r->datadim +2* r->segmentdim + 4,0);

				//TODO define what happens with the other 2 patterns in the closure

#ifndef NDEBUG
				fprintf(stdout,"  \n\t fold_with_pattern_succ_P12: \n");
				fprintf(stdout," builds top because pattern info insufficient \n");
				fprintf(stdout," Warning!! not all the possible incorrect patterns are eliminated \n");
				fflush(stdout);
#endif
				//TODO Update first and last and exit !!!
				//	if(update_lenght) update_lenghts(pr, r, node_tdim, size_node_tdim );
				//return r;
			}
			if(dcons !=NULL )
				dcons=ap_abstract0_join(pr->man_dcons,true, dcons, aux_dcons);
			else {
				dcons = ap_abstract0_copy(pr->man_dcons, aux_dcons);
				ap_abstract0_free(pr->man_dcons, aux_dcons);
			}
#ifndef NDEBUG
	fprintf(stdout,"  \n\t dcons : \n");
	ap_abstract0_fprint(stdout,pr->man_dcons,dcons,NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif
		}

		if(sgli && (i==size-1)){
			/* define last element on the new segment */
			ap_dimchange_t dimadd;
			ap_dimchange_init (&dimadd, 2, 0);
			dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
			dimadd.dim[0] = r->datadim + 2*r->segmentdim;
			dimadd.dim[1] = r->datadim + 2*r->segmentdim;

			dcons_ylx = ap_abstract0_add_dimensions(pr->man_dcons,false,
					r->econs, &dimadd, false);
			ap_dimchange_clear (&dimadd);

#ifndef NDEBUG
	fprintf(stdout,"  \n\t r->econs : \n");
	ap_abstract0_fprint(stdout,pr->man_dcons,r->econs,NULL);
	fprintf(stdout,"\n add 2 dims is \n ");
	ap_abstract0_fprint(stdout,pr->man_dcons,dcons_ylx, NULL);
	fflush(stdout);
#endif
			ap_lincons0_array_t arr = ap_lincons0_array_make (2);
			// l[y1] - l(node_tdim[0] - ... - l(node_tdim[i-1]) == 0
			arr.p[0].constyp = AP_CONS_EQ;
			arr.p[0].linexpr0 =
					ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2 );
			arr.p[0].scalar = NULL;
			ap_dim_t li;
			for(size_t j=0; j < i; j++){
				li = r->datadim + r->segmentdim + node_tdim[j] ;
				ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, li, -1);
			}
			ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, y1, 1);
			ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

			// d(y1) - d(node_tdim[i]) ==0
			arr.p[1].constyp = AP_CONS_EQ;
			arr.p[1].linexpr0 =
					ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
			arr.p[1].scalar = NULL;
			ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
					r->datadim + node_tdim[i], 1);
			dy1 = r->datadim + 2*r->segmentdim + 1;
			ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy1, -1);
			//ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);

			dcons_ylx = ap_abstract0_meet_lincons_array (pr->man_dcons,
					true, dcons_ylx, &arr);
#ifndef NDEBUG
	fprintf(stdout,"\n dcons_ylx is \n ");
	ap_abstract0_fprint(stdout,pr->man_dcons,dcons_ylx, NULL);
	fprintf(stdout," \n");
	fflush(stdout);
#endif
			ap_lincons0_array_clear (&arr);

		}

		if(!sgli){
			ap_abstract0_t *aux_dcons = NULL;
			/* (y1,y2) \in node_tdim[i] _substitute of succ_P12  */
			aux_dcons = create3_succ_P12(pr, r, node_tdim, size_node_tdim, i);

			if(aux_dcons==NULL)
				dcons = ap_abstract0_top(pr->man_dcons,
						r->datadim + 2* r->segmentdim + 4,0);
			else{
				if(dcons !=NULL )
					dcons=ap_abstract0_join(pr->man_dcons,true, dcons, aux_dcons);
				else {
					dcons = ap_abstract0_copy(pr->man_dcons, aux_dcons);
					ap_abstract0_free(pr->man_dcons, aux_dcons);
				}
			}
			/* (y1,y2) = (tdim[i],y=1 ) */
			aux_dcons = create4_succ_P12(pr, r, node_tdim, size_node_tdim, i);
			if(aux_dcons==NULL)
				dcons = ap_abstract0_top(pr->man_dcons,
						r->datadim + 2* r->segmentdim + 4,0);
			else{
				if(dcons !=NULL )
					dcons=ap_abstract0_join(pr->man_dcons,true, dcons, aux_dcons);
				else {
					dcons = ap_abstract0_copy(pr->man_dcons, aux_dcons);
					ap_abstract0_free(pr->man_dcons, aux_dcons);
				}
			}

#ifndef NDEBUG
	fprintf(stdout,"  \n\t dcons : \n");
	ap_abstract0_fprint(stdout,pr->man_dcons,dcons,NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif
			if(i==size-1){
				//TODO update last
				pattern_key_t * lookl = NULL;
				checked_malloc(lookl, pattern_key_t, 1,
						sizeof(pattern_key_t) + 1*sizeof(size_t), return NULL;);
				lookl->type = get_pattern_type(pr, 1, 0, 1, pattern_1_lx_1);
				lookl->segments[0] = node_tdim[i];
				keylen = 1*sizeof(size_t) + sizeof(pattern_key_t);

				pattern_t * last = NULL;
				HASH_FIND(hh,r->udcons, look, keylen, last);
				if ((last==NULL) || (last->dcons == NULL))
					dcons_ylx = ap_abstract0_top(pr->man_dcons,
							r->datadim + 2* r->segmentdim + 2,0);
				else{

					dcons_ylx = ap_abstract0_copy(pr->man_dcons, last->dcons);

					ap_linexpr0_t *yexpr = NULL;

					yexpr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);

					ap_linexpr0_set_coeff_scalar_int (yexpr, y1, 1);
					ap_dim_t li;
					for(size_t j=0; j < i; j++){
						li = r->datadim + r->segmentdim + node_tdim[j] ;
						ap_linexpr0_set_coeff_scalar_int (yexpr, li, -1);
					}
					ap_linexpr0_set_cst_scalar_int (yexpr, 0);
					dcons_ylx = ap_abstract0_substitute_linexpr( pr->man_dcons,
							true, dcons_ylx, y1, yexpr, NULL );

					ap_linexpr0_free(yexpr);

				}
#ifndef NDEBUG
	fprintf(stdout,"\n ax_last \n ");
	ucons_fprint(stdout, pr->man, r, NULL);
	fprintf(stdout," \n");
	fflush(stdout);
#endif



			}

		}


	}

	//add the computed values to the table

	if(dcons_y1!=NULL){
		pattern_key_t * look1 = NULL;
		checked_malloc(look1, pattern_key_t, 1,
				sizeof(pattern_key_t) + 1*sizeof(size_t), return NULL;);
		look1->type = get_pattern_type(pr, 1, 0, 1, pattern_1_l1);
		look1->segments[0] = node_tdim[0];
		keylen = 1*sizeof(size_t) + sizeof(pattern_key_t);

		pattern_t * first = NULL;
		HASH_FIND(hh,r->udcons, look1, keylen, first);
		if(first){
			if(first->dcons!= NULL) ap_abstract0_free(pr->man_dcons,first->dcons);
			first->dcons = dcons_y1;
		}
		else{
			checked_malloc(first,pattern_t,1,
					sizeof(pattern_t)+(1)*sizeof(size_t),return NULL;);
			memset(first, 0, sizeof(pattern_t)+(1)*sizeof(size_t));
			first->key.type = get_pattern_type(pr, 1, 0, 1, pattern_1_l1);
			first->key.segments[0] = node_tdim[0];
			first->dcons = dcons_y1;
			HASH_ADD(hh, r->udcons, key, keylen, first);
			r = add_pattern_n2p(pr, r, &first->key);
		}
		free(look1);
		look1 = NULL;

#ifndef NDEBUG
	fprintf(stdout,"\n ax \n ");
	ucons_fprint(stdout, pr->man, r, NULL);
	fprintf(stdout," \n");
	fflush(stdout);
#endif
	}

	if(dcons_ylx!=NULL){
		pattern_key_t * look1 = NULL;
		checked_malloc(look1, pattern_key_t, 1,
				sizeof(pattern_key_t) + 1*sizeof(size_t), return NULL;);
		look1->type = get_pattern_type(pr, 1, 0, 1, pattern_1_lx_1);
		look1->segments[0] = node_tdim[0];
		keylen = 1*sizeof(size_t) + sizeof(pattern_key_t);

		pattern_t * last = NULL;
		HASH_FIND(hh,r->udcons, look1, keylen, last);
		if(last){
			if(last->dcons!= NULL) ap_abstract0_free(pr->man_dcons,last->dcons);
			last->dcons = dcons_ylx;
		}
		else{
			checked_malloc(last,pattern_t,1,
					sizeof(pattern_t)+(1)*sizeof(size_t),return NULL;);
			memset(last, 0, sizeof(pattern_t)+(1)*sizeof(size_t));
			last->key.type = get_pattern_type(pr, 1, 0, 1, pattern_1_lx_1);
			last->key.segments[0] = node_tdim[0];
			last->dcons = dcons_ylx;
			HASH_ADD(hh, r->udcons, key, keylen, last);
			r = add_pattern_n2p(pr, r, &last->key);
		}
		free(look1);
		look1 = NULL;

#ifndef NDEBUG
	fprintf(stdout,"\n ax \n ");
	ucons_fprint(stdout, pr->man, r, NULL);
	fprintf(stdout," \n");
	fflush(stdout);
#endif

	}

	if(dcons!=NULL){
		pattern_key_t * look1 = NULL;
		checked_malloc(look1, pattern_key_t, 1,
				sizeof(pattern_key_t) + 1*sizeof(size_t), return NULL;);
		look1->type = get_pattern_type(pr, 1, 0, 2, pattern_succ_1_2);
		look1->segments[0] = node_tdim[0];
		keylen = 1*sizeof(size_t) + sizeof(pattern_key_t);

		pattern_t * mid = NULL;
		HASH_FIND(hh,r->udcons, look1, keylen, mid);
		if(mid && mid->dcons){
			ap_abstract0_t * ax = ap_abstract0_copy(pr->man_dcons,mid->dcons);

#ifndef NDEBUG
	fprintf(stdout,"\n ax \n ");
	//ap_abstract0_fprint(stdout,pr->man_dcons,ax, NULL);
	ucons_fprint_dcons(stdout, pr->man, r, NULL, look1);
	fprintf(stdout,"\n dcons \n ");
	ap_abstract0_fprint(stdout,pr->man_dcons,dcons, NULL);
	fprintf(stdout," \n");
	fflush(stdout);
#endif
			mid->dcons = ap_abstract0_join(pr->man_dcons, true, dcons, ax);
		}
		else{
			checked_malloc(mid,pattern_t,1,
					sizeof(pattern_t)+(1)*sizeof(size_t),return NULL;);
			memset(mid, 0, sizeof(pattern_t)+(1)*sizeof(size_t));
			mid->key.type = get_pattern_type(pr, 1, 0, 2, pattern_succ_1_2);
			mid->key.segments[0] = node_tdim[0];
			mid->dcons = ap_abstract0_copy(pr->man_dcons,dcons);
			ap_abstract0_free(pr->man_dcons,dcons);
			HASH_ADD(hh, r->udcons, key, keylen, mid);
			r = add_pattern_n2p(pr, r, &mid->key);
		}
		free(look1);
		look1 = NULL;
	}

	if(update_lenght) update_lenghts(pr, r, node_tdim, size_node_tdim );

#ifndef NDEBUG
	fprintf(stdout,"  \n\t fold_with_pattern_succ_P12 returns:\n ");
	ucons_fprint(stdout,pr->man,r,NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

	return r;
}

/* (y1,y2) = (tdim[i],y=1 ) */
ap_abstract0_t * create4_succ_P12( ucons_internal_t *pr,
		ucons_t * r, ap_dim_t* node_tdim, size_t node_tdim_size, size_t i)
		{
	ap_dim_t y1, y2, dy1, dy2;
	y1 = r->datadim + 2*r->segmentdim;
	y2 = r->datadim + 2*r->segmentdim + 1;
	dy1 = r->datadim + 2*r->segmentdim + 2;
	dy2 = r->datadim + 2*r->segmentdim + 3;

	ap_abstract0_t * aux_dcons = NULL;



	pattern_t *a1 = NULL; //the value with succ_P12

	pattern_key_t * look = NULL;
	checked_malloc(look, pattern_key_t, 1,
			sizeof(pattern_key_t) + 1*sizeof(size_t), return NULL;);
	look->type = get_pattern_type(pr, 1, 0, 1, pattern_1_l1);
	look->segments[0] = node_tdim[i];
	unsigned keylen = 1*sizeof(size_t) + sizeof(pattern_key_t);

	HASH_FIND(hh,r->udcons, look, keylen, a1);

	if((a1==NULL) || (a1->dcons == NULL))
		return NULL;

	ap_dimchange_t dimadd;
	ap_dimchange_init (&dimadd, 2, 0);
	dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
	dimadd.dim[0] = r->datadim + 2*r->segmentdim + 0;
	dimadd.dim[1] = r->datadim + 2*r->segmentdim + 1;

	ap_abstract0_t * right = ap_abstract0_add_dimensions(pr->man_dcons,false,
			a1->dcons, &dimadd, false);
	ap_dimchange_clear (&dimadd);


#ifndef NDEBUG
	fprintf(stdout,"  \n\t create 4 add dim y1  d(y1)  on : \n");
	ap_abstract0_fprint(stdout,pr->man_dcons,a1->dcons,NULL);
	fprintf(stdout,"\n result is \n ");
	ap_abstract0_fprint(stdout,pr->man_dcons,right,NULL);
	fflush(stdout);
#endif

	ap_linexpr0_t *yexpr = NULL;

	yexpr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);

	ap_linexpr0_set_coeff_scalar_int (yexpr, y2, 1);
	ap_dim_t li;
	for(size_t j=0; j < i; j++){
		li = r->datadim + r->segmentdim + node_tdim[j] ;
		ap_linexpr0_set_coeff_scalar_int (yexpr, li, -1);
	}
	ap_linexpr0_set_cst_scalar_int (yexpr, 0);
	right = ap_abstract0_substitute_linexpr( pr->man_dcons,
			true, right, y2, yexpr, NULL );

	ap_linexpr0_free(yexpr);


	ap_dimchange_init (&dimadd, 4, 0);
	dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
	dimadd.dim[0] = r->datadim + 2*r->segmentdim ;
	dimadd.dim[1] = r->datadim + 2*r->segmentdim ;
	dimadd.dim[2] = r->datadim + 2*r->segmentdim ;
	dimadd.dim[3] = r->datadim + 2*r->segmentdim ;

	ap_abstract0_t * left = ap_abstract0_add_dimensions(pr->man_dcons,false,
			r->econs, &dimadd, false);
	ap_dimchange_clear (&dimadd);

#ifndef NDEBUG
	fprintf(stdout,"  \n\t create 4 add dim y1  d(y1)  on : \n");
	ap_abstract0_fprint(stdout,pr->man_dcons,r->econs,NULL);
	fprintf(stdout,"\n result is \n ");
	ap_abstract0_fprint(stdout,pr->man_dcons,left,NULL);
	fflush(stdout);
#endif

	ap_lincons0_array_t arr = ap_lincons0_array_make (2);
	// l[y1] - l(node_tdim[0] - ... - l(node_tdim[i-1]) ==0
	arr.p[0].constyp = AP_CONS_EQ;
	arr.p[0].linexpr0 =
			ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4 );
	arr.p[0].scalar = NULL;

	for(size_t j=0; j < i; j++){
		li = r->datadim + r->segmentdim + node_tdim[j] ;
		ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, li, 1);
	}
	ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, y1, -1);
	//ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

	// d(y1) - d(node_tdim[i]) ==0
	arr.p[1].constyp = AP_CONS_EQ;
	arr.p[1].linexpr0 =
			ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
	arr.p[1].scalar = NULL;
	ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
			r->datadim + node_tdim[i], 1);
	ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy1, -1);
	//ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);

	left = ap_abstract0_meet_lincons_array (pr->man_dcons,
			true, left, &arr);

	ap_lincons0_array_clear (&arr);

	aux_dcons = ap_abstract0_meet(pr->man_dcons, true, left, right);

#ifndef NDEBUG
	fprintf(stdout,"  \n\t create 4 returns aux_dcons : \n");
	ap_abstract0_fprint(stdout,pr->man_dcons,aux_dcons,NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif
	return aux_dcons;

		}

/* (y1,y2) _substitute of succ_P12  */
ap_abstract0_t * create3_succ_P12( ucons_internal_t *pr,
		ucons_t * r, ap_dim_t* node_tdim, size_t node_tdim_size, size_t i)
		{

	ap_dim_t y1, y2, dy1, dy2;
	y1 = r->datadim + 2*r->segmentdim;
	y2 = r->datadim + 2*r->segmentdim + 1;
	dy1 = r->datadim + 2*r->segmentdim + 2;
	dy2 = r->datadim + 2*r->segmentdim + 3;

	ap_abstract0_t * aux_dcons = NULL;



	pattern_t *a1 = NULL; //the value with succ_P12

	pattern_key_t * look = NULL;
	checked_malloc(look, pattern_key_t, 1,
			sizeof(pattern_key_t) + 1*sizeof(size_t), return NULL;);
	look->type = get_pattern_type(pr, 1, 0, 2, pattern_succ_1_2);
	look->segments[0] = node_tdim[i];
	unsigned keylen = 1*sizeof(size_t) + sizeof(pattern_key_t);

	HASH_FIND(hh,r->udcons, look, keylen, a1);

	if((a1==NULL) || (a1->dcons ==NULL))
		return NULL;

	aux_dcons = ap_abstract0_copy(pr->man_dcons, a1->dcons);

#ifndef NDEBUG
	fprintf(stdout,"  \n\t aux_dcons : \n");
	ap_abstract0_fprint(stdout,pr->man_dcons,aux_dcons,NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

	ap_linexpr0_t *yexpr = NULL;

	/* y1 <-- y1 + l[0] + ... + l[i-1] */
	yexpr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);

	ap_linexpr0_set_coeff_scalar_int (yexpr, y1, 1);
	ap_dim_t li;
	for(size_t j=0; j < i; j++){
		li = r->datadim + r->segmentdim + node_tdim[j] ;
		ap_linexpr0_set_coeff_scalar_int (yexpr, li, -1);
	}
	ap_linexpr0_set_cst_scalar_int (yexpr, 0);
	aux_dcons = ap_abstract0_substitute_linexpr( pr->man_dcons,
			true, aux_dcons, y1, yexpr, NULL );

	ap_linexpr0_free(yexpr);


	ap_linexpr0_t *yexpr2 = NULL;

	yexpr2 = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);

	ap_linexpr0_set_coeff_scalar_int (yexpr2, y2, 1);
	for(size_t j=0; j < i; j++){
		li = r->datadim + r->segmentdim + node_tdim[j] ;
		ap_linexpr0_set_coeff_scalar_int (yexpr2, li, -1);
	}
	ap_linexpr0_set_cst_scalar_int (yexpr2, 0);
	aux_dcons = ap_abstract0_substitute_linexpr( pr->man_dcons,
			true, aux_dcons, y2, yexpr2, NULL );

	ap_linexpr0_free(yexpr2);


#ifndef NDEBUG
	fprintf(stdout,"  \n\t create 3 returns : \n");
	ap_abstract0_fprint(stdout,pr->man_dcons,aux_dcons,NULL);
	fprintf(stdout," \n ");
	ucons_fprint(stdout,pr->man,r , NULL);
	fprintf(stdout," \n ");

	fflush(stdout);
#endif

	return aux_dcons;

		}



/*
 * creates the constraint corresponding to (y1,y2) == (last(node_tdim[i-1]), node_tdim[i])
 * (y1,y2) == (last(node_tdim[i-1]), node_tdim[i])
 */
ap_abstract0_t * create2_succ_P12( ucons_internal_t *pr,
		ucons_t * r, ap_dim_t* node_tdim, size_t node_tdim_size, size_t i)
		{
	ap_dim_t y1, y2, dy1, dy2;
	y1 = r->datadim + 2*r->segmentdim;
	y2 = r->datadim + 2*r->segmentdim + 1;
	dy1 = r->datadim + 2*r->segmentdim + 2;
	dy2 = r->datadim + 2*r->segmentdim + 3;

	ap_abstract0_t * aux_dcons = NULL;

	//get last of node_tdim[i-1]

	pattern_t *last = NULL;

	pattern_key_t *look = NULL;
	checked_malloc(look, pattern_key_t, 1,
			sizeof(pattern_key_t) + 1*sizeof(size_t), return NULL;);
	look->type = get_pattern_type(pr, 1, 0, 1, pattern_1_lx_1);
	look->segments[0] = node_tdim[i-1];
	unsigned keylen = 1*sizeof(size_t) + sizeof(pattern_key_t);

	HASH_FIND(hh,r->udcons, look, keylen, last);

	if((last==NULL) || (last->dcons ==NULL)){
#ifndef NDEBUG
		fprintf(stdout,"  \n\t foldt_with_pattern_succ_P12: \n");
		fprintf(stdout," returns with pattern info insufficient \n");
		fflush(stdout);
#endif
		return NULL;
	}

	ap_dimchange_t dimadd;
	ap_dimchange_init (&dimadd, 2, 0);
	dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
	dimadd.dim[0] = r->datadim + 2*r->segmentdim + 1;
	dimadd.dim[1] = r->datadim + 2*r->segmentdim + 2;

	ap_abstract0_t * last_val = ap_abstract0_add_dimensions(pr->man_dcons,false,
			last->dcons, &dimadd, false);
	ap_dimchange_clear (&dimadd);

	//apply substitution on last val y <-- y - l[node_tdim[0]] - ... - l[node_tdim[i-2]]
	ap_linexpr0_t *yexpr = NULL;

	yexpr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);

	ap_linexpr0_set_coeff_scalar_int (yexpr, y1, 1);
	ap_dim_t li;
	for(size_t j=0; j < i-1; j++){
		li = r->datadim + r->segmentdim + node_tdim[j] ;
		ap_linexpr0_set_coeff_scalar_int (yexpr, li, -1);
	}
	ap_linexpr0_set_cst_scalar_int (yexpr, 0);
	last_val= ap_abstract0_substitute_linexpr( pr->man_dcons,
			true, last_val, y1, yexpr, NULL );

	ap_linexpr0_free(yexpr);

	// make the constraint on node_tdim[i]
	ap_dimchange_init (&dimadd, 4, 0);
	dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
	dimadd.dim[0] = r->datadim + 2*r->segmentdim;
	dimadd.dim[1] = r->datadim + 2*r->segmentdim;
	dimadd.dim[2] = r->datadim + 2*r->segmentdim;
	dimadd.dim[3] = r->datadim + 2*r->segmentdim;

	aux_dcons = ap_abstract0_add_dimensions(pr->man_dcons,false,
			r->econs,&dimadd,false);

	ap_dimchange_clear (&dimadd);

	ap_lincons0_array_t arr = ap_lincons0_array_make (3);
	// l[y2] - l(node_tdim[0] - ... - l(node_tdim[i-1]) ==0
	arr.p[0].constyp = AP_CONS_EQ;
	arr.p[0].linexpr0 =
			ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4 );
	arr.p[0].scalar = NULL;
	for(size_t j=0; j < i; j++){
		li = r->datadim + r->segmentdim + node_tdim[j] ;
		ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, li, 1);
	}
	ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, y2, -1);
	//ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

	// d(y2) - d(node_tdim[i]) ==0
	arr.p[1].constyp = AP_CONS_EQ;
	arr.p[1].linexpr0 =
			ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
	arr.p[1].scalar = NULL;
	ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
			r->datadim + node_tdim[i], 1);
	ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy2, -1);
	//ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);

	//y2 - y1 -1  ==0
	arr.p[2].constyp = AP_CONS_EQ;
	arr.p[2].linexpr0 =
			ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
	arr.p[2].scalar = NULL;
	ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0,
			y1, -1);
	ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0, y2, 1);
	ap_linexpr0_set_cst_scalar_int (arr.p[2].linexpr0, -1);

	aux_dcons = ap_abstract0_meet_lincons_array (pr->man_dcons,
			true, aux_dcons, &arr);

	ap_lincons0_array_clear (&arr);


	aux_dcons = ap_abstract0_meet(pr->man_dcons,true,aux_dcons, last_val);


#ifndef NDEBUG
	fprintf(stdout,"  \n\t create 2 returns : \n");
	ap_abstract0_fprint(stdout,pr->man_dcons,aux_dcons,NULL);
	fprintf(stdout," \n ");
	ucons_fprint(stdout,pr->man,r , NULL);
	fprintf(stdout," \n ");
	fflush(stdout);
#endif

	return aux_dcons;

		}
/*(y1,y2) == (node_tdim[i-1], node_tdim[i])*/
ap_abstract0_t * create1_succ_P12( ucons_internal_t *pr,
		ucons_t * r, ap_dim_t* node_tdim, size_t node_tdim_size, size_t i)
		{


	ap_dim_t y1, y2, dy1, dy2;
	y1 = r->datadim + 2*r->segmentdim;
	y2 = r->datadim + 2*r->segmentdim + 1;
	dy1 = r->datadim + 2*r->segmentdim + 2;
	dy2 = r->datadim + 2*r->segmentdim + 3;


	ap_abstract0_t * aux_dcons = NULL;
	ap_dimchange_t dimadd;
	ap_dimchange_init (&dimadd, 4, 0);
	dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
	dimadd.dim[0] = r->datadim + 2*r->segmentdim;
	dimadd.dim[1] = r->datadim + 2*r->segmentdim;
	dimadd.dim[2] = r->datadim + 2*r->segmentdim;
	dimadd.dim[3] = r->datadim + 2*r->segmentdim;

	aux_dcons = ap_abstract0_add_dimensions(pr->man_dcons,false,
			r->econs,&dimadd,false);

	ap_dimchange_clear (&dimadd);


	ap_lincons0_array_t arr = ap_lincons0_array_make (5);

	// l[y1] - l(node_tdim[0] - ... - l(node_tdim[i-2]) ==0
	arr.p[0].constyp = AP_CONS_EQ;
	arr.p[0].linexpr0 =
			ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4 );
	arr.p[0].scalar = NULL;
	ap_dim_t li;
	for(size_t j=0; j < i-1; j++){
		li = r->datadim + r->segmentdim + node_tdim[j] ;
		ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, li, -1);
	}
	ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, y1, 1);
	//ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

	// d(y1) - d(node_tdim[i-1]) ==0
	arr.p[1].constyp = AP_CONS_EQ;
	arr.p[1].linexpr0 =
			ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
	arr.p[1].scalar = NULL;
	ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
			r->datadim + node_tdim[i-1], 1);
	ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy1, -1);
	//ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);


	// l[y2] - l(node_tdim[0] - ... - l(node_tdim[i-1]) ==0
	arr.p[2].constyp = AP_CONS_EQ;
	arr.p[2].linexpr0 =
			ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4 );
	arr.p[2].scalar = NULL;

	for(size_t j=0; j < i; j++){
		li = r->datadim + r->segmentdim + node_tdim[j] ;
		ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0, li, 1);
	}
	ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0, y2, -1);
	//ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

	// d(y2) - d(node_tdim[i]) ==0
	arr.p[3].constyp = AP_CONS_EQ;
	arr.p[3].linexpr0 =
			ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
	arr.p[3].scalar = NULL;
	ap_linexpr0_set_coeff_scalar_int (arr.p[3].linexpr0,
			r->datadim + node_tdim[i], 1);
	ap_linexpr0_set_coeff_scalar_int (arr.p[3].linexpr0, dy2, -1);
	//ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);

	// y2 - y1 -1 ==0
	arr.p[4].constyp = AP_CONS_EQ;
	arr.p[4].linexpr0 =
			ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
	arr.p[4].scalar = NULL;
	ap_linexpr0_set_coeff_scalar_int (arr.p[4].linexpr0,
			y2 , 1);
	ap_linexpr0_set_coeff_scalar_int (arr.p[4].linexpr0, y1, -1);
	ap_linexpr0_set_cst_scalar_int (arr.p[4].linexpr0, -1 );

	aux_dcons = ap_abstract0_meet_lincons_array (pr->man_dcons,
			true, aux_dcons, &arr);

	ap_lincons0_array_clear (&arr);

	return aux_dcons;
		}
