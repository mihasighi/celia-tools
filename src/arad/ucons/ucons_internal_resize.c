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

#include "uthash.h"
#include "ucons.h"
#include "ucons_fun.h"
#include "ucons_internal.h"
#include "shape_macros.h"
#include "apron2shape.h"



ucons_t * fold_with_closure_P11_or_P12( ucons_internal_t *pr,
		ucons_t * a, ap_dim_t * tdim, size_t size,bool update_lenght)
{

	if(size==1)
		/*TODO: make patterns for multiple segments */
		return a;

	ucons_t *r = ucons_copy_internal (pr, a);

#ifndef NDEBUG1
	printf("  \n\t fold : \n");
	ucons_fprint(stdout,pr->man,r,NULL);
	printf("\n");
	fprintf(stdout,"\t tdim =  ");
	for(size_t o = 0;o<size;o++)
		fprintf(stdout,"\t %zu ",tdim[o]-r->datadim);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

	size_t i,j,k,start,size_cdim,size_gdim,size_node_tdim;
	ap_dim_t li,ni;
	ap_linexpr0_t *expr;
	ap_lincons0_t cons;


	ap_dim_t *node_tdim; // the nodes to concatenate from dimensions in tdim

	ap_dim_t *gdim; /* sub-vectors of dimensions (without universals) used to generate data properties */
	ap_dim_t *cdim; /* Dimensions with universals to concatenate */

	size_node_tdim = 0;
	node_tdim = NULL;
	for(i=0; i<size; i++){
		if(tdim[i]>=a->datadim)
			size_node_tdim++;
		checked_realloc(node_tdim,ap_dim_t,sizeof(ap_dim_t),size_node_tdim,return NULL;);
		node_tdim[size_node_tdim-1] = tdim[i] - r->datadim;
	}


	/*tdim splited into sub-vectors gdim such that gdim=ni_start ... ni with l(ni_start)=1 ... l(ni)=1 */
	size_cdim = 0 ;
	cdim = NULL;

	i=0;
	start=i;
	while(i<size_node_tdim){
		ni = node_tdim[i];
		li = r->datadim + r->segmentdim + ni;
		/* li-1=0 */
		expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
		ap_linexpr0_set_coeff_scalar_int (expr, li, 1);
		ap_linexpr0_set_cst_scalar_int (expr, -1);
		cons = ap_lincons0_make(AP_CONS_EQ, expr, NULL);
		if(ap_abstract0_sat_lincons(pr->man_dcons,r->econs,&cons)){
			i++;
		}
		else{

			if (i-start==0){
				checked_realloc(cdim, ap_dim_t, sizeof(ap_dim_t),(size_cdim+1),return NULL;);
				cdim[size_cdim]=node_tdim[start];
				size_cdim++;
				i++;
				start=i;
			}else if (i==start+1){
				checked_realloc(cdim, ap_dim_t, sizeof(ap_dim_t),(size_cdim+2),return NULL;);
				cdim[size_cdim]=node_tdim[start];
				cdim[size_cdim+1]=node_tdim[i];
				size_cdim+=2;
				i++;
				start=i;
			}else{
				size_gdim = i-start;
				checked_malloc(gdim, ap_dim_t, sizeof(ap_dim_t),size_gdim,return NULL;);
				for(j=0;j<size_gdim;j++)
					gdim[j]=node_tdim[j+start+1];

				//#if defined(UCONS_DCONS_OCT_P11) || defined(UCONS_DCONS_POLY_P11)
				if( pr->active_patterns[4]){
					r=merge_succesors_2(pr, r, gdim, size_gdim,update_lenght);
				}
				else{
					if(pr->active_patterns[1]){
						r=merge_succesors_1(pr, r, gdim, size_gdim,update_lenght);
					}
				}
				//#endif
				free(gdim);

				checked_realloc(cdim, ap_dim_t, sizeof(ap_dim_t),(size_cdim+2),return NULL;);
				cdim[size_cdim]=node_tdim[start];
				cdim[size_cdim+1]=node_tdim[i];
				size_cdim+=2;

				i++;
				start=i;
			}
		}

		ap_lincons0_clear(&cons);
		expr=NULL;
	}

	if(i!=start){
		if (i==start+1){
			checked_realloc(cdim, ap_dim_t, sizeof(ap_dim_t),(size_cdim+2),return NULL;);
			cdim[size_cdim]=node_tdim[start];
			size_cdim++;
		}else{
			size_gdim = i-start;
			checked_malloc(gdim, ap_dim_t, sizeof(ap_dim_t),size_gdim,return NULL;);
			for(j=0;j<size_gdim;j++)
				gdim[j]=node_tdim[j+start];

			//#if defined(UCONS_DCONS_OCT_P11) || defined(UCONS_DCONS_POLY_P11)




			if( pr->active_patterns[4]){
				//#elif defined(UCONS_DCONS_OCT_P12) || defined(UCONS_DCONS_POLY_P12)
				r=merge_succesors_2(pr, r, gdim, size_gdim,update_lenght);
			}else{
				if(pr->active_patterns[1]){
					r=merge_succesors_1(pr, r, gdim, size_gdim,update_lenght);
				}
			}
			//#endif

			free(gdim);

			checked_realloc(cdim, ap_dim_t, sizeof(ap_dim_t),size_cdim+1,return NULL;);
			cdim[size_cdim]=node_tdim[start];;
			size_cdim++;
		}
	}


	//#elif defined(UCONS_DCONS_OCT_P12) || defined(UCONS_DCONS_POLY_P12)
	if( pr->active_patterns[4]){
		r=concat_nodes_2(pr,r,cdim,size_cdim,update_lenght);
	}else{
		//#if defined(UCONS_DCONS_OCT_P11) || defined(UCONS_DCONS_POLY_P11)
		if(pr->active_patterns[1]){
			r=concat_nodes_1(pr,r,cdim,size_cdim,update_lenght);
		}
	}
	//#endif

	free(cdim);


#ifndef NDEBUG1
	fprintf(stdout,"  \n\t fold returns : \n");
	ucons_fprint(stdout,pr->man,r,NULL);
	fprintf(stdout,"\n");
	fprintf(stdout,"\t tdim =  ");
	for(size_t o = 0;o<size;o++)
		fprintf(stdout,"\t %zu ",tdim[o]-r->datadim);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif


	return r;
}


ucons_t * split_with_pattern_P12_P11(ucons_internal_t* pr, ucons_t *r,
		pattern_key_t *pattern_j, ap_dim_t n1, ap_dim_t n2){

	pattern_key_t * pat_trasf;
	unsigned keylen;
	size_t u_seg,e_seg,nr_y;
	size_t y_pos,dy_pos , dn2;

	dn2 = r->datadim + n2;


	ap_lincons0_array_t arr;
	pattern_t *aux_find2, *aux, *sub_pattern_data;
	size_t i,ii,j,k;
	ap_abstract0_t *n2_aux,*econs_aux;

	ap_dim_t * tdim;
	size_t size,contor;
	size_t pos_y=0; /* the dimension in r->udcons of y to be eliminated , data(y) = r->datadim + 2*r->segmentdim + nr_y + pos_y*/
	ap_dimperm_t perm;
	ap_dim_t pos_i;

	ap_dimchange_t* dimchange;

	bool pattern_sat;
	//size_t ii, nr_y_ii;
	ap_dim_t lii;
	ap_lincons0_t cons;
	bool found_pos_y;

	u_seg = pr->PI[pattern_j->type].u_seg;
	e_seg = pr->PI[pattern_j->type].e_seg;
	keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

	/* test pattern validity */
	nr_y = pr->PI[pattern_j->type].nr_y;

	HASH_FIND(hh,r->udcons,pattern_j,keylen,aux);
	if(aux){

		pattern_t *aux_find2, *sub_pattern_data;
		ap_dim_t *ydim;
		ap_linexpr0_t ** yexpr;
		pattern_key_t *pat_trasf, *sub_pattern;
		ap_abstract0_t *temp_data;

		/*
		 * pat_tranf = original_pattern : n1_pat_set_p[j]->type
		 */

		checked_malloc(pat_trasf, pattern_key_t, 1,
				sizeof(pattern_key_t) + u_seg*sizeof(size_t), return NULL;);

		pat_trasf->type = pattern_j->type;
		pat_trasf->segments[0] = n2;

		keylen = u_seg*sizeof(size_t) + sizeof(pattern_key_t);
		HASH_FIND(hh,r->udcons, pat_trasf, keylen, aux_find2);
		if(!aux_find2){
/*ifndef NDEBUG1
  				if(pattern_j->type==0){
				fprintf(stdout,"\n pattern 0 added splitting %d on %d ", n1, n2);
				ucons_fprint(stdout,pr->man,r,NULL);
				fflush(stdout);
			}
endif
*/

			checked_malloc(aux_find2,pattern_t,1,sizeof(pattern_t)+(u_seg)*sizeof(size_t),return NULL;);
			memset(aux_find2, 0, sizeof(pattern_t)+(u_seg)*sizeof(size_t));
			aux_find2->key.type = pattern_j->type;
			//for (size_t i=0 ; i<(u_seg); i++)
			aux_find2->key.segments[0] = n2;
			aux_find2->dcons = NULL;
			HASH_ADD(hh,r->udcons,key,keylen,aux_find2);
			r=add_pattern_n2p(pr,r,pat_trasf);
		}
		HASH_FIND(hh,r->udcons, pat_trasf, keylen, aux_find2);
		//		if(aux_find2){
		/* transfer property on n2 */
		/* y substitutes with  y+1  for all y universally quantified */

		checked_malloc(ydim, ap_dim_t, nr_y,sizeof(ap_dim_t), return NULL;);
		checked_malloc(yexpr, ap_linexpr0_t*, nr_y,sizeof(ap_linexpr0_t*), return NULL;);


		ap_abstract0_t *found_val=NULL;

		if(aux_find2->dcons!=NULL)
			found_val = ap_abstract0_copy (pr->man_dcons, aux_find2->dcons);
		aux_find2->dcons = ap_abstract0_copy (pr->man_dcons, aux->dcons);
		if(found_val!=NULL)
				aux_find2->dcons = ap_abstract0_meet (pr->man_dcons, true, aux_find2->dcons,found_val);

		for(i = 0 ; i < nr_y; i++){
			ydim[i] = r->datadim + 2 * r->segmentdim + i;

			yexpr[i] = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2 * nr_y );
			ap_linexpr0_set_coeff_scalar_int (yexpr[i], ydim[i], 1);
			ap_linexpr0_set_cst_scalar_int (yexpr[i], 1);
			aux_find2->dcons= ap_abstract0_substitute_linexpr( pr->man_dcons, true, aux_find2->dcons, ydim[i],
					yexpr[i], NULL );
		}

		for(i = 0 ; i < nr_y; i++){
			ap_linexpr0_free(yexpr[i]);
		}
		free(ydim);
		free(yexpr);

		bool flag = false;
		ap_linexpr0_t *expr;
		ap_lincons0_t cons;

		ap_dim_t dim_y1 = r->datadim + 2*r->segmentdim;
		ap_dim_t dim_y2 = r->datadim + 2*r->segmentdim + 1;

		expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2*nr_y);

		if(pattern_j->type == 0){
			/*
			 * y1 - 1 >= 0
			 */
			/*
			 * y1 < len(n2)
			 */
			ap_linexpr0_set_coeff_scalar_int (expr, dim_y1, 1);
			ap_linexpr0_set_cst_scalar_int (expr, -1);

			ap_linexpr0_t *exprrd0 =
					ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
			ap_linexpr0_set_coeff_scalar_int (exprrd0, dim_y1, -1);
			ap_linexpr0_set_coeff_scalar_int (exprrd0, n2 + r->datadim + r->segmentdim, 1);
			ap_linexpr0_set_cst_scalar_int (exprrd0, -1);

			arr = ap_lincons0_array_make (2);
			arr.p[0] = ap_lincons0_make (AP_CONS_SUPEQ, expr, NULL);
			arr.p[1] = ap_lincons0_make (AP_CONS_SUPEQ, exprrd0, NULL);

			aux_find2->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, aux_find2->dcons, &arr);
			ap_lincons0_array_clear (&arr);
		}
		else  {

			/* test !(y2 <= 1)  */
			ap_linexpr0_set_coeff_scalar_int (expr, dim_y2, -1);
			ap_linexpr0_set_cst_scalar_int (expr, 2);
			cons = ap_lincons0_make(AP_CONS_SUP, expr, NULL);


			flag = ap_abstract0_sat_lincons(pr->man_dcons,aux_find2->dcons,&cons);
			ap_lincons0_clear(&cons);
			/* if aux_find2->dcons implies y2 == 1 */
			if (flag) aux_find2->dcons=
					ap_abstract0_bottom(pr->man_dcons,r->datadim + 2*r->segmentdim + 2*nr_y,0);

			/* impose (y1 >= 1)  */
			expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2*nr_y);
			ap_linexpr0_set_coeff_scalar_int (expr, dim_y1, 1);
			ap_linexpr0_set_cst_scalar_int (expr, -1);

			/* impose (y2 <= ln2 - 1) */
			ap_linexpr0_t *exprrd =
					ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2*nr_y);
			ap_linexpr0_set_coeff_scalar_int (exprrd, dim_y2, -1);
			ap_linexpr0_set_coeff_scalar_int (exprrd, n2 + r->datadim + r->segmentdim, 1);
			ap_linexpr0_set_cst_scalar_int (exprrd, -1);
			arr = ap_lincons0_array_make (2);

			arr.p[0] = ap_lincons0_make (AP_CONS_SUPEQ, expr, NULL);
			arr.p[1] = ap_lincons0_make (AP_CONS_SUPEQ, exprrd, NULL);


			aux_find2->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, aux_find2->dcons, &arr);
			ap_lincons0_array_clear (&arr);

		}

		if (nr_y>1){
			/*
			 * y2 - y1 >= 0
			 */
			expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2*nr_y);
			ap_linexpr0_set_coeff_scalar_int (expr, dim_y1, -1);

			ap_dim_t dim_y2 = r->datadim + 2*r->segmentdim + 1;
			ap_linexpr0_set_coeff_scalar_int (expr, dim_y2, 1);

			arr = ap_lincons0_array_make (1);
			arr.p[0] = ap_lincons0_make (AP_CONS_SUPEQ, expr, NULL);

			aux_find2->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, aux_find2->dcons, &arr);
			ap_lincons0_array_clear (&arr);

		}
#ifndef NDEBUG1

		fprintf(stdout,"\n  after transfering ctr on n2 %d on %d ", n1, n2);
		ucons_fprint(stdout,pr->man,r,NULL);
		fflush(stdout);

		/*fprintf(stdout,"\t pattern: \n");
		//pattern_key_fprint (stdout, pr, pattern_j, NULL);
		fprintf(stdout,"\n type %d segment %d ", pattern_j->type, pattern_j->segments[0]);
		fflush(stdout);
		fprintf(stdout,"\t => aux : \n");
		ap_abstract0_fprint(stdout,pr->man_dcons,aux_find2->dcons, NULL);
		fprintf(stdout,"\n");*/
#endif



		//}

		free(pat_trasf);

		/** TODO: case two: test to see if sub-patterns needed */
		//#endif
		if(pr->active_patterns[4]){
			//#if defined (UCONS_DCONS_OCT_P12) || defined (UCONS_DCONS_POLY_P12)

			//			ap_dim_t segment = pattern_j->segment[0];
			//
			//			nr_y = pr->PI[pattern_j->type].nr_y;
			//
			//			expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2*nr_y);
			//			ap_linexpr0_set_coeff_scalar_int (expr, ln1, 1);
			//			ap_linexpr0_set_cst_scalar_int (expr, -1);
			//			arr = ap_lincons0_array_make (1);
			//			arr.p[0] = ap_lincons0_make (AP_CONS_EQ, expr, NULL);
			//
			//			bool flag_ln1 = ap_abstract0_sat_lincons()

			if ((pattern_j->type == 3) &&
					//	test_pattern_sat(pr,r,pattern_j,0)){
					!ap_abstract0_is_bottom(pr->man_dcons, aux->dcons)){
				/*
				 * pat_tranf != original_pattern (sub-pattern of the original one)
				 * for the pattern P12 y1<y2 \in n1
				 * modify          P11 y2\in n2  by considering d(y1) = d(n2) and y1 = 1
				 */
				size_t y1, y2, dy1, dy2;

				// if pattern is P12
				// TODO: add constraint and order pattern add general case for
				// 		 any number of order constraints

				checked_malloc(sub_pattern, pattern_key_t, 1,
						sizeof(pattern_key_t) + 1*sizeof(size_t), return NULL;);
				sub_pattern->type = 0; // P11
				sub_pattern->segments[0] = n2;

				/* temp_data size is r->datadim + 2*r->segmentdim + 4 */
				temp_data = ap_abstract0_copy(pr->man_dcons, aux->dcons);

				y1 = r->datadim + 2*r->segmentdim;
				y2 = r->datadim + 2*r->segmentdim + 1;
				dy1 = r->datadim + 2*r->segmentdim + 2;
				dy2 = r->datadim + 2*r->segmentdim + 3;

				arr = ap_lincons0_array_make (2);
				/*
				 *  y1-1 == 0
				 */
				arr.p[0].constyp = AP_CONS_EQ;
				arr.p[0].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
				arr.p[0].scalar = NULL;
				ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,y1, 1);
				ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);

				/* dy1 - dn2 == 0
				 */
				arr.p[1].constyp = AP_CONS_EQ;
				arr.p[1].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
				arr.p[1].scalar = NULL;
				ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, (dy1), -1);
				ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dn2, 1);
				ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);


				temp_data = ap_abstract0_meet_lincons_array (pr->man_dcons, true, temp_data, &arr);
				ap_lincons0_array_clear (&arr);


#ifndef NDEBUG1
				fprintf(stdout,"\n \t temp_data after y1 - 1 = 0 and dy1 - dn2 = 0 : \n");
				ap_abstract0_fprint(stdout,pr->man_dcons,temp_data, NULL);
				fprintf(stdout,"\n");
				fflush(stdout);
#endif


				/* eliminate y1 and dy1 to get constraints with the pattern P11(n2) */
				checked_malloc(tdim,ap_dim_t,2,sizeof(ap_dim_t),return NULL;);
				tdim[0] = y1;
				tdim[1] = dy1;
				temp_data = ap_abstract0_forget_array(pr->man_dcons, true, temp_data, tdim, 1, true);
				//  remove dimension y1, dy1
				dimchange=ap_dimchange_alloc(2,0);
				dimchange->dim[0] = y1;
				dimchange->dim[1] = dy1;
				temp_data = ap_abstract0_remove_dimensions (pr->man_dcons, true, temp_data, dimchange);
				ap_dimchange_free(dimchange);

#ifndef NDEBUG1
				fprintf(stdout,"\n \t temp_data after removed y1 and dy1 \n");
				ap_abstract0_fprint(stdout,pr->man_dcons,temp_data, NULL);
				fprintf(stdout,"\n");
				fflush(stdout);
#endif

				/* in temp_data y1 := y1 + 1*/

				ap_dim_t dim_y1 = r->datadim + 2*r->segmentdim;
				ap_linexpr0_t *expr;
				expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
				ap_linexpr0_set_coeff_scalar_int (expr, dim_y1, 1);
				ap_linexpr0_set_cst_scalar_int (expr, 1);

				temp_data = ap_abstract0_substitute_linexpr( pr->man_dcons, true, temp_data, dim_y1,
						expr, NULL );

				ap_linexpr0_free(expr);

				/* insure y1>=1 */
				ap_lincons0_array_t arr = ap_lincons0_array_make (1);
				// y1 - 1 >=0
				arr.p[0].constyp = AP_CONS_SUPEQ;
				arr.p[0].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
				arr.p[0].scalar = NULL;
				ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
						dim_y1, 1);
				ap_linexpr0_set_cst_scalar_int(arr.p[0].linexpr0,-1);

				temp_data = ap_abstract0_meet_lincons_array (pr->man_dcons, true, temp_data, &arr);
				ap_lincons0_array_clear (&arr);




#ifndef NDEBUG1
				fprintf(stdout,"\n \t temp_data after incremeted y1 \n");
				ap_abstract0_fprint(stdout,pr->man_dcons,temp_data, NULL);
				fprintf(stdout,"\n");
				fflush(stdout);
#endif
				/*
				 *  modify P12 strengthen it with P11
				 */
				if(aux_find2){
					/*y2 - y1 >= 0 */
					ap_lincons0_array_t arr = ap_lincons0_array_make (1);
					arr.p[0].constyp = AP_CONS_SUPEQ;
					arr.p[0].linexpr0 =
							ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
					arr.p[0].scalar = NULL;
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
												dim_y1, -1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
												dim_y2, 1);

					aux_find2->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, aux_find2->dcons, &arr);
					ap_lincons0_array_clear (&arr);


					ap_abstract0_t * temp_data_p12 = ap_abstract0_copy(pr->man_dcons,temp_data);
					/*add back y1 and dy1 normaly also y2 and dy2
					 */

					ap_dimchange_t dimchange;
					/* add y1, d(y1) to ucons1->dcons and y2 = y , d(y2) = d(y) */
					dimchange.intdim = 2;
					dimchange.realdim = 0;
					dimchange.dim = (ap_dim_t*) malloc (2*sizeof(ap_dim_t));
					dimchange.dim[0] = r->datadim + 2*r->segmentdim + 1;
					dimchange.dim[1] = r->datadim + 2*r->segmentdim + 2;

					temp_data_p12 = ap_abstract0_add_dimensions(pr->man_dcons, false, temp_data_p12, &dimchange, false);
					free(dimchange.dim);

					/* y2 - y1 >= 0*/
					arr = ap_lincons0_array_make (1);
					arr.p[0].constyp = AP_CONS_SUPEQ;
					arr.p[0].linexpr0 =
							ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
					arr.p[0].scalar = NULL;
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
												dim_y1, -1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
												dim_y2, 1);

					temp_data_p12 = ap_abstract0_meet_lincons_array (pr->man_dcons, true, temp_data_p12, &arr);
					ap_lincons0_array_clear (&arr);

					aux_find2->dcons = ap_abstract0_meet(pr->man_dcons,true,aux_find2->dcons,temp_data_p12);

#ifndef NDEBUG1
					fprintf(stdout,"\n pattern 0 added splitting %d on %d ", n1, n2);
					ucons_fprint(stdout,pr->man,r,NULL);
					fflush(stdout);

			/*
					fprintf(stdout,"\n \t temp_data2 for P12  \n");
ap_abstract0_fprint(stdout,pr->man_dcons,temp_data_p12, NULL);
fprintf(stdout,"\n");
fflush(stdout)*/;
#endif
				}


				keylen = 1*sizeof(size_t) + sizeof(pattern_key_t);
				HASH_FIND(hh,r->udcons, sub_pattern, keylen, sub_pattern_data);

				if(sub_pattern_data){
					if (sub_pattern_data->dcons!=NULL){
						sub_pattern_data->dcons =
								ap_abstract0_meet(pr->man_dcons,true,sub_pattern_data->dcons,temp_data);
/*
						fprintf(stdout,"\n sub_pattern_data!= null \n sub-pattern => ");
						ap_abstract0_fprint(stdout,pr->man_dcons,sub_pattern_data->dcons, NULL);
						fprintf(stdout,"\n meet with temp_data: ");
						ap_abstract0_fprint(stdout,pr->man_dcons,temp_data, NULL);
						fprintf(stdout,"\n");
						fflush(stdout);
*/
					}
					else{
						sub_pattern_data->dcons =
								ap_abstract0_copy(pr->man_dcons,temp_data);
/*						fprintf(stdout,"\n sub_pattern_data= null \n ");
						fprintf(stdout,"\n");
						fflush(stdout);
*/
						ap_abstract0_free(pr->man_dcons,temp_data);
					}
				}
				else{
/*
					fprintf(stdout,"\n no sub_pattern \n ");
					fprintf(stdout,"\n");
					fflush(stdout);

*/
					checked_malloc(sub_pattern_data,pattern_t,1,sizeof(pattern_t)+1*sizeof(size_t),return NULL;);
					memset(sub_pattern_data, 0, sizeof(pattern_t)+(1)*sizeof(size_t));
					sub_pattern_data->key.type = sub_pattern->type;
					for (size_t i=0 ; i<1; i++)
						sub_pattern_data->key.segments[i] = sub_pattern->segments[i];
					sub_pattern_data->dcons =ap_abstract0_copy(pr->man_dcons,temp_data);
					ap_abstract0_free(pr->man_dcons,temp_data);

					HASH_ADD(hh,r->udcons,key,keylen,sub_pattern_data);
					r=add_pattern_n2p(pr,r,sub_pattern);
/*

					fprintf(stdout,"  \n\t ucons_afer add sub-pattern: \n");
					fprintf(stdout,"\n");
					ucons_fprint(stdout,pr->man,r,NULL);
					printf("\n");
					fflush(stdout);
*/

				}

#ifndef NDEBUG1
				fprintf(stdout,"\t pattern: \n");
				//pattern_key_fprint (stdout, pr, pattern_j, NULL);
				fprintf(stdout,"\n type %d segment %d ", sub_pattern->type, sub_pattern->segments[0]);
				fflush(stdout);
				fprintf(stdout,"\t => aux : \n");
				ap_abstract0_fprint(stdout,pr->man_dcons,sub_pattern_data->dcons, NULL);
				fprintf(stdout,"\n");
#endif
				free(sub_pattern);
			}
			//#endif
		}//end if active_pattern P21


		//}else{

		//} old else
		//}//end if pattern correct


		/* set old constraint to bottom */
		ap_abstract0_free(pr->man_dcons,aux->dcons);
		aux->dcons = NULL;

		remove_pattern_n2p(pr,r,&aux->key);
		HASH_DEL(r->udcons,aux);
		free(aux);

		//			aux->dcons=ap_abstract0_bottom(pr->man_dcons, r->datadim + 2*r->segmentdim + 2*nr_y, 0);
		//			size_t total_dim = r->datadim + 2*r->segmentdim + 2*nr_y;
		//
		//			arr = ap_lincons0_array_make (1);
		//			// l[dim]-1 == 0
		//			arr.p[0].constyp = AP_CONS_EQ;
		//			arr.p[0].linexpr0 =
		//					ap_linexpr0_alloc (AP_LINEXPR_DENSE, total_dim);
		//			arr.p[0].scalar = NULL;
		//			ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
		//					r->datadim + r->segmentdim + n1, 1);
		//			ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);
		//
		//
		//			aux->dcons =
		//					ap_abstract0_meet_lincons_array (pr->man_dcons, true, aux->dcons, &arr);
		//
		//			ap_lincons0_array_clear (&arr);

	}//end checking one pattern

	return r;
}



ucons_t *  ucons_change_pattern_2_1(ucons_internal_t * pr, ucons_t *r, ucons_t * a, ap_dim_t * tdim,
		size_t size_tdim, pattern_key_t  * pattern_j){

	pattern_kind kind = pr->PI[pattern_j->type].kind;
	size_t nr_y = pr->PI[pattern_j->type].nr_y; /* 2 */
	size_t e_seg = pr->PI[pattern_j->type].e_seg; /*>=1*/
	size_t u_seg = pr->PI[pattern_j->type].u_seg; /* 2 */
	size_t hd = tdim[0];
	pattern_t *aux = NULL;

	if(e_seg >= size_tdim){
		bool found = false;
		size_t pos_hd = 0;
		for(size_t jj=u_seg; jj<e_seg+u_seg && !found; jj++){
			if (pattern_j->segments[jj] == hd) { found = true; pos_hd = jj; }
		}
		if(found){
			size_t jj = pos_hd;
			for(size_t kk=1; kk<size_tdim && found; kk++){
				if (!pattern_j->segments[jj] == tdim[kk]) found = false;
				jj++;
			}
			if(found){
				/* update pattern */
				pattern_key_t *look;
				checked_malloc(look, pattern_key_t, sizeof(pattern_key_t) + (u_seg+(e_seg-size_tdim+1))*sizeof(size_t), 1, return NULL;);
				memset(look, 0, sizeof(pattern_key_t) +(u_seg+(e_seg-size_tdim+1))*sizeof(size_t));
				look->type = get_pattern_type(pr, u_seg,(e_seg-size_tdim+1),nr_y,kind );
				unsigned keylen = (u_seg+(e_seg-size_tdim+1))*sizeof(size_t) + sizeof(pattern_key_t);
				look->segments[0] = pattern_j->segments[0];
				look->segments[1] = pattern_j->segments[1];
				for(size_t jj=u_seg;jj<(u_seg+e_seg-size_tdim+1); jj++){
					if(jj<=pos_hd) look->segments[jj] = pattern_j->segments[jj];
					else look->segments[jj] = pattern_j->segments[jj+size_tdim];
				}

#ifndef NDEBUG1
				printf("  \n\t debug chainge_pattern_2_1 \n look: \n");
				for(size_t ii = 0; ii<u_seg+e_seg-size_tdim+1 ;ii ++)
					printf("look[%zu]=%zu \n",ii,look->segments[ii]);
				printf("\n");
#endif
				pattern_t *res, *curr;
				HASH_FIND(hh,r->udcons,look,keylen,res);
				unsigned keylen_curr = (u_seg+e_seg)*sizeof(size_t) + sizeof(pattern_key_t);
				HASH_FIND(hh,a->udcons,pattern_j,keylen_curr,curr);
				if(curr){
					if(res)
						res->dcons = ap_abstract0_join(pr->man_dcons,false,curr->dcons,res->dcons);
					else{
						checked_malloc(res,pattern_t,1,
								(sizeof(pattern_t)+(u_seg+(e_seg-size_tdim+1))*sizeof(size_t)),return NULL;);
						memset(res, 0, (sizeof(pattern_t)+(u_seg+(e_seg-size_tdim+1))*sizeof(size_t)));

						res->dcons = ap_abstract0_copy(pr->man_dcons, curr->dcons);
						res->key.type= look->type;

						for(size_t ii = 0; ii <(u_seg+(e_seg-size_tdim+1)); ii++){
							res->key.segments[ii] = look->segments[ii];
						}
						HASH_ADD(hh,r->udcons,key,keylen,res);
						add_pattern_n2p(pr,r,look);
					}
				}
			}
		}
	}
	r=remove_pattern_n2p(pr, r, pattern_j);
	unsigned keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);
	HASH_FIND(hh,r->udcons,pattern_j,keylen,aux);
	if(aux) {
		HASH_DEL(r->udcons,aux);
		free(aux);
	}
#ifndef NDEBUG1
	printf("  \n\t ucons_fold_pattern_1_lx returns : \n");
	ucons_fprint(stdout,pr->man,r,NULL);
	printf("\n");
#endif


	return r;
}

ucons_t *  ucons_fold_pattern_2_1(ucons_internal_t * pr, ucons_t *r, ucons_t * a, ap_dim_t * tdim,
		size_t size_tdim, pattern_key_t  * pattern_j){

	ap_dim_t hd = tdim[0];
	ap_dim_t hn;
	size_t upos_hd;
	size_t i,ctr;

	pattern_kind kind = pr->PI[pattern_j->type].kind;
	size_t nr_y = pr->PI[pattern_j->type].nr_y; /* 2 */
	size_t e_seg = pr->PI[pattern_j->type].e_seg;
	size_t u_seg = pr->PI[pattern_j->type].u_seg; /* 2 */



	/* fold */
	bool ord = true; /* true if tdim[0] < hn */
	if (pattern_j->segments[0] == hd) {
		hn =  pattern_j->segments[1];
		ord = true;
	}
	else{
		hn =  pattern_j->segments[0];
		ord = false;
	}

	ap_dim_t * len_e_seg = NULL;
	size_t le_seg = e_seg;
	checked_malloc(len_e_seg, ap_dim_t, e_seg, sizeof(ap_dim_t), return NULL;);
	for(i=0; i<e_seg; i++)
		len_e_seg[i] = pattern_j->segments[u_seg+i];

	pattern_key_t *look1, *look2;
	pattern_t *ucons1, *ucons2, *ucons_res;

	ap_abstract0_t *dcons=NULL;
	ap_abstract0_t *dcons2, *dcons1;

	unsigned keylen,keylen1, keylen2;

	keylen = (u_seg+e_seg)*sizeof(size_t) + sizeof(pattern_key_t);
	HASH_FIND(hh, a->udcons, pattern_j, keylen, ucons_res);

	if(!test_singleton(pr->man_dcons,a->econs,a->datadim,a->segmentdim,tdim[0])){
		if(ucons_res==NULL){
			/**/
			HASH_FIND(hh, r->udcons, pattern_j, keylen, ucons_res);
			if(ucons_res){
				r=remove_pattern_n2p(pr, r, pattern_j);
				HASH_DEL(r->udcons,ucons_res);
				free(ucons_res);
				ucons_res=NULL;
			}
			return r;
		}
		else
			dcons=ap_abstract0_copy(pr->man_dcons,ucons_res->dcons);
	}


	for(ctr = 1; ctr<size_tdim; ctr++){

		/* add len[tdim[ctr-1]] to the constraints */
		checked_realloc(len_e_seg, ap_dim_t,sizeof(ap_dim_t) , (le_seg + 1), return NULL;);
		size_t ii = le_seg;
		while((ii>=1) && (len_e_seg[ii-1]>tdim[ctr-1])){
			len_e_seg[ii] = len_e_seg[ii-1];
			ii--;
		}
#ifndef NDEBUG1
		fprintf(stdout,"\n 731 le_seg =%zu ii= %zu ctr=%zu \n", le_seg,ii,ctr);
		fflush(stdout);
#endif
		len_e_seg[ii] = tdim[ctr-1];
		le_seg += 1;

		/* look1 : y_hn = l[tdim[ctr-1] + .... */
		checked_malloc(look1, pattern_key_t,
				sizeof(pattern_key_t) + (1+le_seg)*sizeof(size_t), 1, return NULL;);
		look1->type = get_pattern_type(pr,1,le_seg,1,pattern_1_lx);
		look1->segments[0] = hn;

		for(size_t kk = 0;kk < le_seg; kk++){
			look1->segments[kk+1] = len_e_seg[kk];
		}
		keylen1 = (1+le_seg)*sizeof(size_t) + sizeof(pattern_key_t);
		HASH_FIND(hh, a->udcons, look1, keylen1, ucons1);
		if(!ucons1){
			/* no constraint with pattern in closure returns top with the analyzed pattern */
			keylen = (u_seg+e_seg)*sizeof(size_t) + sizeof(pattern_key_t);
			HASH_FIND(hh, r->udcons, pattern_j, keylen, ucons_res);
			if(ucons_res){
				HASH_DEL(r->udcons,ucons_res);
				free(ucons_res);
				ucons_res=NULL;
			}
			return r;
		}
		/* check for second pattern only if len of tdim[ctr] not one */


		if(!test_singleton(pr->man_dcons,r->econs,r->datadim,r->segmentdim, tdim[ctr])){
			bool ordCtr=true; /*true if samme order as between tdim[0] and hn  */
			checked_malloc(look2, pattern_key_t,
					sizeof(pattern_key_t) + (u_seg+le_seg)*sizeof(size_t), 1, return NULL;);
			if(tdim[ctr] < hn){
				look2->segments[0] = tdim[ctr];
				look2->segments[1] = hn;
				look2->type = get_pattern_type(pr,u_seg,le_seg,2,pattern_2_1_lx);
				if (ord) ordCtr = true;
				else ordCtr = false ;
			}
			else{
				look2->segments[0] = hn;
				look2->segments[1] = tdim[ctr];
				look2->type = get_pattern_type(pr,u_seg,le_seg,2,pattern_2_1_mlx);
				if (ord) ordCtr = false;
				else ordCtr = true ;
			}

			for(size_t kk = 0;kk < le_seg; kk++){
				look2->segments[kk+2] = len_e_seg[kk];
			}
			keylen2 = (u_seg+le_seg)*sizeof(size_t) + sizeof(pattern_key_t);
			HASH_FIND(hh, a->udcons, look2, keylen2, ucons2);

			if(!ucons2){
				/* no constraint with pattern in closure returns top with the analyzed pattern */
				keylen = (u_seg+e_seg)*sizeof(size_t) + sizeof(pattern_key_t);
				HASH_FIND(hh, r->udcons, pattern_j, keylen, ucons_res);
				if(ucons_res){
					HASH_DEL(r->udcons,ucons_res);
					free(ucons_res);
					ucons_res=NULL;
				}
				return r;
			}
			if (ordCtr == false){
				/* permute y1,y2 and d(y1),d(y2)*/
				ap_dimperm_t consperm;
				consperm.size = r->datadim + 2 * r->segmentdim + 2 * nr_y;
				consperm.dim = (ap_dim_t *) malloc (consperm.size * sizeof (ap_dim_t));
				ap_dimperm_set_id (&consperm);

				for(i=0;i<u_seg;i++){
					consperm.dim[r->datadim + 2 * r->segmentdim + 2*i] = r->datadim + 2 * r->segmentdim + (2*i)+1;
					consperm.dim[r->datadim + 2 * r->segmentdim + 2*i + 1] = r->datadim + 2 * r->segmentdim + (2*i);
				}
				dcons2 = ap_abstract0_permute_dimensions(pr->man_dcons,false,ucons2->dcons,&consperm);
				free(consperm.dim);
			}
			else
				dcons2 = ap_abstract0_copy(pr->man_dcons,ucons2->dcons);
		}
		else{
			dcons2 = ap_abstract0_bottom(pr->man_dcons,r->datadim + 2*r->segmentdim + 4,0);
		}

#ifndef NDEBUG1
		printf("  \n\t dcons2 : \n");
		ap_abstract0_fprint(stdout,pr->man_dcons,dcons2,NULL);
		printf("\n");
#endif
		ap_lincons0_array_t arr = ap_lincons0_array_make (2);
		// l[y1] == l[fst_node_dim[0]] + ... l[fst_node_dim[i-1]]
		arr.p[0].constyp = AP_CONS_EQ;
		arr.p[0].linexpr0 =
				ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4 );
		arr.p[0].scalar = NULL;

		ap_dim_t lj, ly, dy;
		for(size_t j = 0; j < le_seg; j++){
			lj = r->datadim + r->segmentdim + len_e_seg[j] ;
			ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, lj, 1);
		}
		ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

		// d(y) - d(tdim[crt]) ==0
		ap_dim_t dctr =  tdim[ctr];

		arr.p[1].constyp = AP_CONS_EQ;
		arr.p[1].linexpr0 =
				ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
		arr.p[1].scalar = NULL;
		ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
				r->datadim + dctr, 1);

		ap_dimchange_t dimchange;
		if(ord){
			/* add y1, d(y1) to ucons1->dcons and y2 = y , d(y2) = d(y) */
			dimchange.intdim = 2;
			dimchange.realdim = 0;
			dimchange.dim = (ap_dim_t*) malloc (2*sizeof(ap_dim_t));
			dimchange.dim[0] = r->datadim + 2*r->segmentdim ;
			dimchange.dim[1] = r->datadim + 2*r->segmentdim + 1;

			dcons1 = ap_abstract0_add_dimensions(pr->man_dcons, false, ucons1->dcons, &dimchange, false);

			free(dimchange.dim);
			/* y1 = len(len_e_seg[0]) + len(len_e_seg[0])+ ... +len(len_e_seg[..]) */
			/* d(y1) = d(tdim[crt]) */
			ly = r->datadim + 2*r->segmentdim;
			dy = r->datadim + 2*r->segmentdim + 2;

		}
		else{
			/* add y2, d(y2) to ucons1->dcons and y1 = y , d(y1) = d(y) */
			dimchange.intdim = 2;
			dimchange.realdim = 0;
			dimchange.dim = (ap_dim_t*) malloc (2*sizeof(ap_dim_t));
			dimchange.dim[0] = r->datadim + 2*r->segmentdim +1;
			dimchange.dim[1] = r->datadim + 2*r->segmentdim + 2;

			dcons1 = ap_abstract0_add_dimensions(pr->man_dcons, false, ucons1->dcons, &dimchange, false);
			free(dimchange.dim);
			/* y2 = len(len_e_seg[0]) + len(len_e_seg[0])+ ... +len(len_e_seg[..]) */
			/* d(y2) = d(tdim[crt]) */
			ly = r->datadim + 2*r->segmentdim + 1;
			dy = r->datadim + 2*r->segmentdim + 3;
		}

		ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,ly, -1);
		ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy, -1);
		dcons1 = ap_abstract0_meet_lincons_array (pr->man_dcons,true, dcons1,&arr);
		ap_lincons0_array_clear (&arr);

#ifndef NDEBUG1
		printf("  \n\t dcons1 : \n");
		ap_abstract0_fprint(stdout,pr->man_dcons,dcons1,NULL);
		printf("\n");
#endif

		if(dcons){
			dcons = ap_abstract0_join(pr->man_dcons,false,dcons,dcons1);
			dcons = ap_abstract0_join(pr->man_dcons,false,dcons,dcons2);
		}
		else
			dcons = ap_abstract0_join(pr->man_dcons,false,dcons1,dcons2);
		ap_abstract0_free(pr->man_dcons,dcons1);
		ap_abstract0_free(pr->man_dcons,dcons2);

#ifndef NDEBUG1
		printf("  \n\t dcons : \n");
		ap_abstract0_fprint(stdout,pr->man_dcons,dcons,NULL);
		printf("\n");
#endif
	}

	/* joined all the nodes in tdimp[1..size_tdim] */
	/* set the new constraint on tdim[0] */

	keylen = (u_seg+e_seg)*sizeof(size_t) + sizeof(pattern_key_t);
	HASH_FIND(hh, r->udcons, pattern_j, keylen, ucons_res);

	if(ucons_res){
		ap_abstract0_free(pr->man_dcons, ucons_res->dcons);
		ucons_res->dcons = NULL;
		ucons_res->dcons = dcons;
	}
	else{
		checked_malloc(ucons_res,pattern_t,1,
				(sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t)),return NULL;);
		memset(ucons_res, 0, (sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t)));

		ucons_res->dcons = dcons;
		ucons_res->key.type= pattern_j->type;

		for(size_t ii = 0; ii <(u_seg+e_seg); ii++){
			ucons_res->key.segments[ii] = pattern_j->segments[ii];
		}
		HASH_ADD(hh,r->udcons,key,keylen,ucons_res);
		add_pattern_n2p(pr,r,pattern_j);

	}

#ifndef NDEBUG1
	printf("  \n\t generate_pattern_2_1 returns : \n");
	ucons_fprint(stdout,pr->man,r,NULL);
	printf("\n");
#endif

	return r;
}
//ucons_t *  ucons_fold_pattern_2_1_lx(ucons_internal_t * pr, ucons_t *r, ucons_t * a, ap_dim_t * tdim,
//		size_t size_tdim, pattern_key_t  * pattern_j){
//	return r;
//}
//ucons_t *  ucons_fold_pattern_2_1_mlx(ucons_internal_t * pr, ucons_t *r, ucons_t * a, ap_dim_t * tdim,
//		size_t size_tdim, pattern_key_t  * pattern_j){
//	return r;
//}

ucons_t *  ucons_fold_pattern_1_lx(ucons_internal_t * pr, ucons_t *r, ucons_t * a, ap_dim_t * tdim,
		size_t size_tdim, pattern_key_t  * pattern_j){

	ap_dim_t hd = tdim[0];
	pattern_t * aux;

	pattern_kind kind = pr->PI[pattern_j->type].kind;
	size_t nr_y = pr->PI[pattern_j->type].nr_y; //1//
	size_t e_seg = pr->PI[pattern_j->type].e_seg;
	size_t u_seg = pr->PI[pattern_j->type].u_seg; //1//
	if(pattern_j->segments[0] == hd){
		/* case 1: y\in hd. y = l[] + l[] + ...*/
		bool found = false;
		for(size_t jj=0; jj<e_seg && !found; jj++){
			for(size_t kk = 1; (kk < size_tdim && !found); kk++)
				found = (pattern_j->segments[jj+1] == tdim[kk])? true : false;
		}
		if(found){
			r=remove_pattern_n2p(pr, r, pattern_j);
			unsigned keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);
			HASH_FIND(hh,r->udcons,pattern_j,keylen,aux);
			if(aux) {
				HASH_DEL(r->udcons,aux);
				free(aux);
			}
		}
	}
	else{
		/* case 2 y\in hm. y = l[] + l[hd] + ...*/
		/* if y = l[.. ]+ l[tdim[0]] + l[tdim[1]] + ... + l[tdim[size_tdim-1]] + l[...]
		 * then transport the constraint on y = l[..] + l[hd] + l[..]
		 * */
		if(e_seg >= size_tdim){
			bool found = false;
			size_t pos_hd = 0;
			for(size_t jj=u_seg; jj<(e_seg+u_seg) && !found; jj++){
				if (pattern_j->segments[jj] == hd) { found = true; pos_hd = jj; }
			}
			if(found){
				size_t jj = pos_hd;
				for(size_t kk=1; kk<size_tdim && found; kk++){
					if (!pattern_j->segments[jj] == tdim[kk]) found = false;
					jj++;
				}
				if(found){
					/* update pattern */
					pattern_key_t *look;
					checked_malloc(look, pattern_key_t, sizeof(pattern_key_t) + (u_seg+(e_seg-size_tdim+1))*sizeof(size_t), 1, return NULL;);
					memset(look, 0, sizeof(pattern_key_t) +(u_seg+(e_seg-size_tdim+1))*sizeof(size_t));
					look->type = get_pattern_type(pr, u_seg,(e_seg-size_tdim+1),1,kind );
					unsigned keylen = (u_seg+(e_seg-size_tdim+1))*sizeof(size_t) + sizeof(pattern_key_t);
					look->segments[0] = pattern_j->segments[0];
					for(size_t jj=u_seg;jj<(u_seg+e_seg-size_tdim+1); jj++){
						if(jj<=pos_hd) look->segments[jj] = pattern_j->segments[jj];
						else look->segments[jj] = pattern_j->segments[jj+size_tdim];
					}

#ifndef NDEBUG1
					printf("  \n\t debug fold_pattern_1_lx \n look: \n");
					for(size_t ii = 0; ii<u_seg+e_seg-size_tdim+1 ;ii ++)
						printf("look[%zu]=%zu \n",ii,look->segments[ii]);
					printf("\n");
#endif
					pattern_t *res, *curr;
					HASH_FIND(hh,r->udcons,look,keylen,res);
					unsigned keylen_curr = (u_seg+e_seg)*sizeof(size_t) + sizeof(pattern_key_t);
					HASH_FIND(hh,a->udcons,pattern_j,keylen_curr,curr);
					if(curr){
						if(res)
							res->dcons = ap_abstract0_join(pr->man_dcons,false,curr->dcons,res->dcons);
						else{
							checked_malloc(res,pattern_t,1,
									(sizeof(pattern_t)+(u_seg+(e_seg-size_tdim+1))*sizeof(size_t)),return NULL;);
							memset(res, 0, (sizeof(pattern_t)+(u_seg+(e_seg-size_tdim+1))*sizeof(size_t)));

							res->dcons = ap_abstract0_copy(pr->man_dcons, curr->dcons);
							res->key.type= look->type;

							for(size_t ii = 0; ii <(u_seg+(e_seg-size_tdim+1)); ii++){
								res->key.segments[ii] = look->segments[ii];
							}
							HASH_ADD(hh,r->udcons,key,keylen,res);
							add_pattern_n2p(pr,r,look);
						}
					}
				}
			}
		}
		r=remove_pattern_n2p(pr, r, pattern_j);
		unsigned keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);
		HASH_FIND(hh,r->udcons,pattern_j,keylen,aux);
		if(aux) {
			HASH_DEL(r->udcons,aux);
			free(aux);
		}
	}
#ifndef NDEBUG1
	printf("  \n\t ucons_fold_pattern_1_lx returns : \n");
	ucons_fprint(stdout,pr->man,r,NULL);
	printf("\n");
#endif
	return r;
}


/* when fold_with_closoure_of_P21 is called from fold_without_closoure_of_P21 then
 * update_length == false
 * 	     and
 * minus_dim > 0
 * */
ucons_t * fold_with_closoure_of_P21(ucons_internal_t * pr, ucons_t * a, ap_dim_t * node_tdim,
		size_t size_node_tdim, ap_dim_t minus_dim, bool update_length){


#ifndef NDEBUG1
	fprintf(stdout,"  \n\t fold_with_closure_ P21 : \n");
	ucons_fprint(stdout,pr->man,a,NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

	arg_assert((update_length ||(!update_length&& ((minus_dim>0)|| (pr->active_patterns[1])))), return NULL;);

	if(size_node_tdim==1)
		/* no segments to fold*/
		return a;

	ucons_t *r = ucons_copy_internal (pr, a);
	size_t i,j,k,start;

	ap_dim_t li,ni;

	ap_dim_t *tdim; // the nodes to concatenate from dimensions in tdim
	size_t size_tdim=0;

	pattern_t *aux;

	size_tdim = 0;
	tdim = NULL;
	for(i=0; i<size_node_tdim; i++){
		//de ce f\u{a}ceam asa ??
		if(update_length == true || (!update_length&&pr->active_patterns[1])){ // PB
			// incoding of 0 2
			if(node_tdim[i]>=a->datadim){
				size_tdim++;
				checked_realloc(tdim,ap_dim_t,sizeof(ap_dim_t),size_node_tdim,return NULL;);
				tdim[size_tdim-1] = node_tdim[i] - r->datadim;
			}
		}
		else{
#ifndef NDEBUG1
			fprintf(stdout,"  \njust test3 \n");
			fflush(stdout);
#endif
			size_tdim++;
			checked_realloc(tdim,ap_dim_t,sizeof(ap_dim_t),size_node_tdim,return NULL;);
			tdim[size_tdim-1] = node_tdim[i];
		}
	}

	ap_dim_t hd = tdim[0];
	/* look for segment hm s.t. hm, hd are related by a pattern with two universals
	 * look also for segments hm st hd and hm are related by a pattern with one universal and hd is an e_seg
	 * 	todo for more precision, identify all segments ehead of equal length with hd, and
	 *       do the previous search for ehead instead of hd
	 *       particular saturation function. Deduces more patterns and constraints. Build it sepaatelly and
	 *       apply it before the folding procedure.
	 * 				 */

	/* r = saturate_lengths(pr, r, tdim, size_tdim); */

	bool hd_sgl;

	hd_sgl = test_singleton(pr->man_dcons,a->econs,a->datadim,a->segmentdim,hd);

	pattern_key_set_t patterns_on_hd = a->n2p[hd];

	bool ignore_pattern;
	for(j = 0; j < patterns_on_hd.size; j++){


		ignore_pattern = false;

		if (patterns_on_hd.p[j]!=NULL &&
				patterns_on_hd.p[j]->type == 1 &&
				((patterns_on_hd.p[j]->segments[0] == hd && patterns_on_hd.p[j]->segments[1] == minus_dim )
						||
						(patterns_on_hd.p[j]->segments[1] == hd && patterns_on_hd.p[j]->segments[0] == minus_dim )))
			ignore_pattern =true;

		/* if the fold_with_closure is called from fold_without_closure
		 * then we must ignore (not delete or remove the pattern
		 * analised by fold_without_closure namely
		 * \forall y1\in hd, y2\in minus_dim y1=y2
		 * */
		/*todo check !=NULL is needed */
		if(patterns_on_hd.p[j]!=NULL && ignore_pattern == false){


			pattern_key_t *	pattern_j = patterns_on_hd.p[j];
			size_t type = pattern_j->type;
			size_t u_seg = pr->PI[pattern_j->type].u_seg;
			size_t e_seg = pr->PI[pattern_j->type].e_seg;
#ifndef NDEBUG2
			fprintf(stdout,"  \n\t pattern analysed fold_with_closure_ P21 : \n");
			fprintf(stdout,"type= %zu seg[0] = %zu ",patterns_on_hd.p[j]->type,patterns_on_hd.p[j]->segments[0]);
			if (u_seg>1) fprintf(stdout,"seg[1] = %zu",patterns_on_hd.p[j]->segments[1]);
			fprintf(stdout,"\n ");
			fflush(stdout);
#endif
			/*
			 *  ((pr->PI[pattern_j->type].kind == pattern_2_1)||
					(pr->PI[pattern_j->type].kind == pattern_2_1_lx)||
							(pr->PI[pattern_j->type].kind == pattern_2_1_mlx)||
			 */

			if(pr->PI[pattern_j->type].kind == pattern_1_lx){
				r = ucons_fold_pattern_1_lx( pr, r,  a,  tdim, size_tdim,  pattern_j);
				if(hd_sgl){
					if(pr->PI[pattern_j->type].e_seg == 1){

						pattern_key_t *new_key;
						checked_malloc(new_key, pattern_key_t, sizeof(pattern_key_t) + (1+u_seg)*sizeof(size_t), 1, return NULL;);
						memset(new_key, 0, sizeof(pattern_key_t) + (1+u_seg)*sizeof(size_t));
						new_key->type = get_pattern_type(pr,2,0,2,pattern_2_1);
						if(pattern_j->segments[0]<pattern_j->segments[1]){
							new_key->segments[0] = pattern_j->segments[0];
							new_key->segments[1] = pattern_j->segments[1];
						}
						else{
							new_key->segments[0] = pattern_j->segments[1];
							new_key->segments[1] = pattern_j->segments[0];
						}
						/*
						 * don't do anything on minus_dim and hd
						 */

						if (!((new_key->segments[0]==hd &&new_key->segments[1]==minus_dim)||
								(new_key->segments[1]==hd &&new_key->segments[0]==minus_dim))){

#ifndef NDEBUG2
							fprintf(stdout,"  \n\t hd = sgl pattern analysed fold_with_closure_ P21 : \n");
							fprintf(stdout,"type= %zu new_key[0] = %zu ",new_key->type,new_key->segments[0]);
							fprintf(stdout,"new_key[1] = %zu",new_key->segments[1]);
							fprintf(stdout,"\n");
							fflush(stdout);
#endif
							r = ucons_fold_pattern_2_1(pr, r,  a, tdim,size_tdim, new_key);
						}
						/* do not free new_key, possibly added to n2p in generate */
					}
					else{
						bool ntdim=false;
						size_t pos_hd = 0;
						for(size_t ii = 1; (ii<1+pr->PI[pattern_j->type].e_seg)&& !ntdim; ii++){
							if(pattern_j->segments[ii]==tdim[0]) pos_hd = ii;
							for(size_t jj=1;((jj<size_tdim) && !ntdim); jj++)
								if(pattern_j->segments[ii]==tdim[jj])
									ntdim = true;
						}
						if(!ntdim){
							pattern_key_t *new_key;
							checked_malloc(new_key, pattern_key_t, sizeof(pattern_key_t) + u_seg*sizeof(size_t), 1, return NULL;);
							memset(new_key, 0, sizeof(pattern_key_t) + (u_seg+e_seg-1)*sizeof(size_t));
							if(pattern_j->segments[0]<pattern_j->segments[1]){
								new_key->segments[0] = pattern_j->segments[0];
								new_key->segments[1] = pattern_j->segments[1];
							}
							else{
								new_key->segments[0] = pattern_j->segments[1];
								new_key->segments[1] = pattern_j->segments[0];
							}
							size_t nk=2;
							for(size_t kk=0;kk<e_seg; kk++){
								if(pattern_j->segments[kk] == tdim[0]) kk++;
								else new_key->segments[nk++] = pattern_j->segments[kk];
							}
							new_key->type = get_pattern_type(pr,2,0,2,pattern_2_1_lx);
							r = ucons_fold_pattern_2_1(pr, r,  a, tdim,size_tdim, new_key);
							new_key->type = get_pattern_type(pr,2,0,2,pattern_2_1_mlx);
							r = ucons_fold_pattern_2_1(pr, r,  a, tdim,size_tdim, new_key);
						}
					}
				}
			}
			if(pr->PI[pattern_j->type].kind == pattern_2_1){
				if ((pattern_j->segments[0] == tdim[0] )|| (pattern_j->segments[1] == tdim[0] ))
					r = ucons_fold_pattern_2_1(pr, r,  a, tdim,size_tdim, pattern_j);
				else
					r = ucons_change_pattern_2_1(pr, r,  a, tdim,size_tdim, pattern_j);
			}
			if(pr->PI[pattern_j->type].kind == pattern_2_1_lx){
				if ((pattern_j->segments[0] == tdim[0] )|| (pattern_j->segments[1] == tdim[0] ))
					r = ucons_fold_pattern_2_1(pr, r,  a, tdim,size_tdim, pattern_j);
				else
					r = ucons_change_pattern_2_1(pr, r,  a, tdim,size_tdim, pattern_j);
			}
			if(pr->PI[pattern_j->type].kind == pattern_2_1_mlx){
				if ((pattern_j->segments[0] == tdim[0] )|| (pattern_j->segments[1] == tdim[0] ))
					r = ucons_fold_pattern_2_1(pr, r,  a, tdim,size_tdim, pattern_j);
				else
					r = ucons_change_pattern_2_1(pr, r,  a, tdim,size_tdim, pattern_j);
			}
		}
	}
	/* update length constraints in r->econs and r->udcons */
	/* recalculate length of tdim[0]  */
	if(update_length){
		ap_dim_t lhd;
		/* update length for the merged segments */
		/* l(tdim[0]) = l(tdim[0]) + l(tdim[1]) + ... + l(tdim[last])*/

		/* for econs */

		ap_linexpr0_t *expr =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
		ap_linexpr0_set_cst_scalar_int (expr, 0);
		for(i = 0;i < size_tdim; i++){
			li = r->datadim + r->segmentdim + tdim[i] ;
			ap_linexpr0_set_coeff_scalar_int (expr, li, 1);
		}
		/* for udcons */
		ap_linexpr0_t *expr_1y =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
		ap_linexpr0_set_cst_scalar_int (expr_1y, 0);
		for(i=0;i<size_tdim;i++){
			li = r->datadim + r->segmentdim + tdim[i] ;
			ap_linexpr0_set_coeff_scalar_int (expr_1y, li, 1);
		}
		ap_linexpr0_t *expr_2y =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
		ap_linexpr0_set_cst_scalar_int (expr_2y, 0);
		for(i=0;i<size_tdim;i++){
			li = r->datadim + r->segmentdim + tdim[i] ;
			ap_linexpr0_set_coeff_scalar_int (expr_2y, li, 1);
		}

		lhd = r->datadim + r->segmentdim + tdim[0] ;

		r->econs =  ap_abstract0_assign_linexpr (pr->man_dcons, true, r->econs,lhd , expr, NULL);
		ap_linexpr0_free(expr);

		pattern_t *s;
		for(s=r->udcons; s!=NULL; s=s->hh.next){
			if(pr->PI[s->key.type].nr_y == 1){
				s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons,lhd , expr_1y, NULL);
			}
			else{
				s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons,lhd , expr_2y, NULL);
			}
		}
		ap_linexpr0_free(expr_1y);
		ap_linexpr0_free(expr_2y);
	}

#ifndef NDEBUG1
	printf("  \n\t fold_with_closure_of_P21 returns : \n");
	ucons_fprint(stdout,pr->man,r,NULL);
	printf("\n");
#endif

	return r;
}


ucons_t * fold_without_closoure_of_P21(ucons_internal_t * pr, ucons_t * a, ap_dim_t *fold_dim,
		size_t size_fold_dim, bool update_length){

	ucons_t *r = ucons_copy_internal (pr, a);

#ifndef NDEBUG1
	fprintf(stdout,"  \n\t fold without closure: \n");
	ucons_fprint(stdout, pr->man, a, NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

	size_t i, j, k;
	ap_dim_t fst_ni, snd_ni;

	ap_dimchange_t dimchange;

	pattern_key_t * look;
	size_t u_seg, nr_y;
	unsigned keylen;
	pattern_t * n_dcons;

	ap_linexpr0_t *expr;
	ap_lincons0_t cons;

	ap_dim_t *snd_node_tdim = NULL; // the nodes to concatenate from dimensions in tdim
	ap_dim_t *fst_node_tdim = NULL; // the nodes to concatenate from dimensions in tdim

	ap_dim_t ddy1, ddy2, l01, l02, ly1, ly2 ,li, lj;

	size_t	size_fst_node_tdim = 0;
	size_t size_snd_node_tdim = 0;

	ap_abstract0_t * aux = NULL;
	ap_abstract0_t * old_pattern_ctr = NULL;

	size_t size = size_fold_dim/2;
	size_fold_dim = size_fold_dim/2;

	ap_dim_t *tdim = NULL;
	checked_realloc(tdim, ap_dim_t, sizeof(ap_dim_t),size,  return NULL;);

	for(i = 0; i < size; i++){
		tdim[i] = fold_dim[i+size];
	}
	checked_realloc(fold_dim, ap_dim_t,
			sizeof(ap_dim_t), size_fold_dim,return NULL;);

	for(i = 0; i < size; i++){
		if(tdim[i] >= r->datadim){
			size_snd_node_tdim++;
			checked_realloc(snd_node_tdim, ap_dim_t,
					sizeof(ap_dim_t),size_snd_node_tdim, return NULL;);
			snd_node_tdim[size_snd_node_tdim-1] = tdim[i] - r->datadim;
		}
	}
	for(i = 0; i < size_fold_dim; i++){
		if(fold_dim[i] >= r->datadim){
			size_fst_node_tdim++;
			checked_realloc(fst_node_tdim, ap_dim_t,
					sizeof(ap_dim_t),size_fst_node_tdim, return NULL;);
			fst_node_tdim[size_fst_node_tdim-1] = fold_dim[i] - r->datadim;
		}
	}

#ifndef NDEBUG1
	fprintf(stdout,"  \n\t fold_dim: ");
	for(size_t ll= 0 ;ll<size_fold_dim; ll++)
		fprintf(stdout,"  %zu", fold_dim[ll]);
	fprintf(stdout,"\n");
	fprintf(stdout,"  \n\t fst_node_tdim: ");
	for(size_t ll= 0 ;ll<size_fst_node_tdim; ll++)
		fprintf(stdout," %zu", fst_node_tdim[ll]);
	fprintf(stdout,"\n");
	for(size_t ll= 0 ;ll<size_snd_node_tdim; ll++)
		fprintf(stdout," %zu", snd_node_tdim[ll]);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

	/*
	 * fold singular segments
	 *
	 */
	ucons_t *rr = fold_with_closoure_of_P21(pr,r,fst_node_tdim, size_fst_node_tdim,snd_node_tdim[0],false);
	ucons_free_internal(pr,r);
	r = NULL;
	r = rr;
	if(pr->active_patterns[1]){
		rr = fold_with_closure_P11_or_P12(pr,r,fst_node_tdim,size_fst_node_tdim,false);
		ucons_free_internal(pr,r);
		r = NULL;
		r = rr;
	}

#ifndef NDEBUG2
	fprintf(stdout,"\t after fold with closure first segment : \n");
	ucons_fprint(stdout,pr->man,r, NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif
	rr = fold_with_closoure_of_P21(pr,r,snd_node_tdim, size_snd_node_tdim,fst_node_tdim[0],false);
	ucons_free_internal(pr,r);
	r = NULL;
	r = rr;
	if(pr->active_patterns[1]){
		rr = fold_with_closure_P11_or_P12(pr,r,snd_node_tdim,size_snd_node_tdim,false);
		ucons_free_internal(pr,r);
		r = NULL;
		r = rr;
	}
#ifndef NDEBUG2
	fprintf(stdout,"\t after fold with closure second segment: \n");
	ucons_fprint(stdout,pr->man,r, NULL);
	fprintf(stdout,"\n \t starting without closure \n ");
	fflush(stdout);
#endif



	if(fst_node_tdim[0] ==  snd_node_tdim[0])
	{
		/* fold snd_node_dim (fst_node_dim \subseteq snd_node_dim )*/

		ap_manager_raise_exception (pr->man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
				"not implemented");

	}
	if(!test_equal_length(pr->man_dcons,r->econs,r->datadim,r->segmentdim,fst_node_tdim[0], snd_node_tdim[0])){
		pattern_t *nd;
		checked_malloc(nd,pattern_t,1,(sizeof(pattern_t)+2*sizeof(size_t)),return NULL;);
		memset(nd, 0, (sizeof(pattern_t)+2*sizeof(size_t)));
		nd->dcons = ap_abstract0_top(pr->man_dcons, r->datadim + 2*r->segmentdim + 2*2, 0);
		nd->key.type = 1;
		if (fst_node_tdim[0] <  snd_node_tdim[0]){
			nd->key.segments[0] = fst_node_tdim[0];
			nd->key.segments[1] = snd_node_tdim[0];
		}
		else{
			nd->key.segments[1] = fst_node_tdim[0];
			nd->key.segments[0] = snd_node_tdim[0];
		}
		keylen = 2*sizeof(size_t) + sizeof(pattern_key_t);
		HASH_ADD(hh,r->udcons,key,keylen,nd);
		add_pattern_n2p(pr,r,&nd->key);
#ifndef NDEBUG2
		fprintf(stdout,"\t fold fst_node_tdim = ");
		for(size_t ll= 0 ;ll<size_fst_node_tdim; ll++)
			fprintf(stdout," %zu", fst_node_tdim[ll]);
		fprintf(stdout,"\n");
		for(size_t ll= 0 ;ll<size_snd_node_tdim; ll++)
			fprintf(stdout," %zu", snd_node_tdim[ll]);
		fprintf(stdout,"\n");
		fprintf(stdout,"with len[%zu]!=len[%zu] \n vi insufficient information \n fold_without_closure returns: \n",fst_node_tdim[0],snd_node_tdim[0]);
		ucons_fprint(stdout,pr->man,r, NULL);
		fflush(stdout);
#endif
		return r;
	}
	/* folding with pattern_type == 1 over two segments */
	u_seg = pr->PI[1].u_seg;
	nr_y = pr->PI[1].nr_y;
	checked_malloc(look, pattern_key_t, sizeof(pattern_key_t) + u_seg*sizeof(size_t), 1, return NULL;);
	memset(look, 0, sizeof(pattern_key_t) + u_seg*sizeof(size_t));
	look->type = 1;
	keylen = u_seg*sizeof(size_t) + sizeof(pattern_key_t);

	ap_abstract0_t * val = NULL;
	ap_abstract0_t * val_i = NULL;

	{
		if(fst_node_tdim[0] < snd_node_tdim[0]){
			look->segments[0] = fst_node_tdim[0];
			look->segments[1] = snd_node_tdim[0];
		}
		else{
			look->segments[1] = fst_node_tdim[0];
			look->segments[0] = snd_node_tdim[0];
		}
		n_dcons =NULL;
		HASH_FIND(hh, r->udcons, look, keylen, n_dcons);
		if(n_dcons) {
			old_pattern_ctr = ap_abstract0_copy(pr->man_dcons,n_dcons->dcons);
			ap_abstract0_free(pr->man_dcons,n_dcons->dcons);
			n_dcons->dcons = ap_abstract0_bottom(pr->man_dcons, r->datadim + 2*r->segmentdim + 2*nr_y, 0);
		}
		else {
			checked_malloc(n_dcons,pattern_t,1,(sizeof(pattern_t)+u_seg*sizeof(size_t)),return NULL;);
			memset(n_dcons, 0, (sizeof(pattern_t)+u_seg*sizeof(size_t)));

			n_dcons->dcons = ap_abstract0_bottom(pr->man_dcons, r->datadim + 2*r->segmentdim + 2*nr_y, 0);
			n_dcons->key.type= look->type;

			n_dcons->key.segments[0] = look->segments[0];
			n_dcons->key.segments[1] = look->segments[1];

			HASH_ADD(hh,r->udcons,key,keylen,n_dcons);
			add_pattern_n2p(pr,r,look);
		}

		if (old_pattern_ctr) aux = ap_abstract0_copy(pr->man_dcons,old_pattern_ctr);

#ifndef NDEBUG2
		fprintf(stdout,"\t initially aux from heads  : \n");
		if(aux) ap_abstract0_fprint(stdout,pr->man_dcons,aux, NULL);
		else fprintf(stdout,"\t  is empty heads are: %zu, %zu  \n",  look->segments[0], look->segments[1] );
		fflush(stdout);
#endif

		for(i = 1; i < size; i++){
			fst_ni = fst_node_tdim[i];
			snd_ni = snd_node_tdim[i];

			expr =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim );
			ap_linexpr0_set_cst_scalar_int (expr, 0);
			for(j = 0; j < i; j++){
				lj = r->datadim + r->segmentdim + fst_node_tdim[j] ;
				ap_linexpr0_set_coeff_scalar_int (expr, lj, 1);
				lj = r->datadim + r->segmentdim + snd_node_tdim[j] ;
				ap_linexpr0_set_coeff_scalar_int (expr, lj, -1);
			}
			cons = ap_lincons0_make(AP_CONS_EQ, expr, NULL);
			bool fst = false;
			// test length(fst_node_tdim[0], fst_node_tdim[i]) ==  length(snd_node_tdim[0], snd_node_tdim[i])
			if(ap_abstract0_sat_lincons(pr->man_dcons,r->econs, &cons)){

				if(fst_node_tdim[i] < snd_node_tdim[i]){
					look->segments[0] = fst_node_tdim[i];
					look->segments[1] = snd_node_tdim[i];
					fst = true;
				}
				else{
					look->segments[1] = fst_node_tdim[i];
					look->segments[0] = snd_node_tdim[i];
				}

				HASH_FIND(hh,r->udcons, look, keylen, n_dcons);
				if(n_dcons){

					ap_abstract0_t *transf_n_dcons = ap_abstract0_copy(pr->man_dcons,n_dcons->dcons);

					ap_linexpr0_t *fst_expr_y =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
					ap_linexpr0_t *snd_expr_y =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
					ap_linexpr0_set_cst_scalar_int (fst_expr_y, 0);/*TODO check the definition of the sustitution */
					ap_linexpr0_set_cst_scalar_int (snd_expr_y, 0);/*TODO check the definition of the sustitution */
					for(j = 0; j < i; j++){
						lj = r->datadim + r->segmentdim + fst_node_tdim[j] ;
						ap_linexpr0_set_coeff_scalar_int (fst_expr_y, lj, -1);
						lj = r->datadim + r->segmentdim + snd_node_tdim[j] ;
						ap_linexpr0_set_coeff_scalar_int (snd_expr_y, lj, -1);
					}
					if(fst){

						ly1 = r->datadim + 2*r->segmentdim;
						ly2 = r->datadim + 2*r->segmentdim + 1;
					}
					else{
						ly1 = r->datadim + 2*r->segmentdim +1;
						ly2 = r->datadim + 2*r->segmentdim;
					}
					ap_linexpr0_set_coeff_scalar_int (fst_expr_y, ly1, 1);
					ap_linexpr0_set_coeff_scalar_int (snd_expr_y, ly2, 1);

					transf_n_dcons = ap_abstract0_substitute_linexpr (pr->man_dcons,
							true, transf_n_dcons,ly1 , fst_expr_y, NULL);
					transf_n_dcons = ap_abstract0_substitute_linexpr (pr->man_dcons,
							true, transf_n_dcons,ly2 , snd_expr_y, NULL);

					if (aux) aux = ap_abstract0_join(pr->man_dcons, true, aux, transf_n_dcons);
					else {
						aux = ap_abstract0_copy(pr->man_dcons,transf_n_dcons);
						ap_abstract0_free(pr->man_dcons,transf_n_dcons);
						transf_n_dcons = NULL;
					}
#ifndef NDEBUG2
					fprintf(stdout,"\t intermidiate fold_without_closure returns aux right : \n");
					if(aux)ap_abstract0_fprint(stdout,pr->man_dcons,aux, NULL);
					fflush(stdout);
#endif
				}
				dimchange.intdim = 4;
				dimchange.realdim = 0;
				dimchange.dim = (ap_dim_t *)malloc(4*sizeof(ap_dim_t));
				dimchange.dim[0] = r->datadim + 2*r->segmentdim;
				dimchange.dim[1] = r->datadim + 2*r->segmentdim;
				dimchange.dim[2] = r->datadim + 2*r->segmentdim;
				dimchange.dim[3] = r->datadim + 2*r->segmentdim;


				val_i = ap_abstract0_add_dimensions(pr->man_dcons,false, r->econs, &dimchange, false);
				free(dimchange.dim);

				/*
				 * fac mai intai o substitutie ca sa tina poliedrele
				 *
				 * */




				ap_lincons0_array_t arr = ap_lincons0_array_make (4);
				// l[y1] == l[fst_node_dim[0]] + ... l[fst_node_dim[i-1]]
				arr.p[0].constyp = AP_CONS_EQ;
				arr.p[0].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4 );
				arr.p[0].scalar = NULL;

				// l[y2] == l[snd_node_dim[0]] + ... l[snd_node_dim[i-1]]
				arr.p[1].constyp = AP_CONS_EQ;
				arr.p[1].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4 );
				arr.p[1].scalar = NULL;

				for(j = 0; j < i; j++){
					lj = r->datadim + r->segmentdim + fst_node_tdim[j] ;
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, lj, 1);
					lj = r->datadim + r->segmentdim + snd_node_tdim[j] ;
					ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, lj, 1);
				}
				if(fst_node_tdim[0]<snd_node_tdim[0]){
					ly1 = r->datadim + 2*r->segmentdim;
					ly2 = r->datadim + 2*r->segmentdim + 1;
					ddy1 = fst_node_tdim[i];
					ddy2 = snd_node_tdim[i];
				}
				else{
					ly1 = r->datadim + 2*r->segmentdim +1;
					ly2 = r->datadim + 2*r->segmentdim;
					ddy2 = fst_node_tdim[i];
					ddy1 = snd_node_tdim[i];
				}
				ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,ly1, -1);
				ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

				ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,ly2, -1);
				ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);

				arr.p[2].constyp = AP_CONS_EQ;
				arr.p[2].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
				arr.p[2].scalar = NULL;
				ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0,
						r->datadim + ddy1, 1);
				ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0,
						r->datadim + 2*r->segmentdim + 2, -1);

				// d(y2) - d(snd_node_dim[i]) ==0
				arr.p[3].constyp = AP_CONS_EQ;
				arr.p[3].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
				arr.p[3].scalar = NULL;
				ap_linexpr0_set_coeff_scalar_int (arr.p[3].linexpr0,
						r->datadim + ddy2, 1);
				ap_linexpr0_set_coeff_scalar_int (arr.p[3].linexpr0,
						r->datadim + 2*r->segmentdim + 3, -1);

				ap_linexpr0_t *expr_y_tmp =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
				ap_linexpr0_set_cst_scalar_int (expr_y_tmp, 0);

				ap_linexpr0_set_coeff_scalar_int (expr_y_tmp,r->datadim + 2*r->segmentdim + 2, 1);
				ap_dim_t dn1 = r->datadim + fst_node_tdim[i];

				val_i = ap_abstract0_meet_lincons_array (pr->man_dcons,true, val_i,&arr);
				ap_lincons0_array_clear (&arr);


				//				ap_abtract0_t * val2_i = ap_abtract0_copy(pr->man_dcons,val_i);
				//
				//				val2_i = ap_abstract0_substitute_linexpr (pr->man_dcons, true, val2_i,dn1 , expr_y_tmp, NULL);
				//				ap_linexpr0_free(expr_y_tmp);
				//
				//				expr_y_tmp =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
				//				ap_linexpr0_set_cst_scalar_int (expr_y_tmp, 0);
				//
				//				ap_linexpr0_set_coeff_scalar_int (expr_y_tmp, r->datadim + 2*r->segmentdim + 3 , 1);
				//				ap_dim_t dn2 = r->datadim + snd_node_tdim[i];
				//
				//				val2_i = ap_abstract0_substitute_linexpr (pr->man_dcons, true, val2_i, dn2 , expr_y_tmp, NULL);
				//				ap_linexpr0_free(expr_y_tmp);

#ifndef NDEBUG2
				fprintf(stdout,"\t intermidiate fold_without_closure returns aux i : \n");
				if(val_i) ap_abstract0_fprint(stdout,pr->man_dcons,val_i, NULL);
				fflush(stdout);
#endif
				if(val_i){
					if (aux){
#ifndef NDEBUG2
						fprintf(stdout,"\t intermidiate join val_i cu aux %zu : \n", i);
						if(aux) ap_abstract0_fprint(stdout,pr->man_dcons,aux, NULL);
						fprintf(stdout,"\n");
						if(val_i) ap_abstract0_fprint(stdout,pr->man_dcons,val_i, NULL);
						fflush(stdout);
#endif
						aux = ap_abstract0_join(pr->man_dcons,true, aux, val_i);
					}
					else {
						aux = ap_abstract0_copy(pr->man_dcons, val_i);
						ap_abstract0_free(pr->man_dcons,val_i);
					}
					val_i = NULL;
				}
#ifndef NDEBUG2
				fprintf(stdout,"\t intermidiate fold_without_closure returns aux %zu : \n", i);
				if(aux) ap_abstract0_fprint(stdout,pr->man_dcons,aux, NULL);
				fflush(stdout);
#endif
			}
			else{
				aux = ap_abstract0_top(pr->man_dcons, r->datadim + 2*r->segmentdim + 4, 0);
			}
#ifndef NDEBUG2
			fprintf(stdout,"\t intermidiate fold_without_closure returns aux (right \\/ i) : \n");
			ap_abstract0_fprint(stdout,pr->man_dcons,aux, NULL);
			fflush(stdout);
#endif
			ap_linexpr0_free(expr);

		}//end for every pair of nodes to fold
		// join with constraints on first
		if(fst_node_tdim[0] < snd_node_tdim[0]){
			look->segments[0] = fst_node_tdim[0];
			look->segments[1] = snd_node_tdim[0];
		}
		else{
			look->segments[1] = fst_node_tdim[0];
			look->segments[0] = snd_node_tdim[0];
		}

		HASH_FIND(hh, r->udcons, look, keylen, n_dcons);
		if(n_dcons){
			if (n_dcons->dcons) n_dcons->dcons = ap_abstract0_join(pr->man_dcons,true, aux, n_dcons->dcons);
			else {
				n_dcons->dcons = ap_abstract0_copy(pr->man_dcons, aux);
				ap_abstract0_free(pr->man_dcons,aux);
			}
		}
		else{
			pattern_t *rn;

			u_seg = pr->PI[1].u_seg;
			keylen = u_seg*sizeof(size_t)+sizeof(pattern_key_t);

			checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+u_seg*sizeof(size_t),return NULL ;);
			memset(rn, 0, sizeof(pattern_t)+u_seg*sizeof(size_t));

			rn->dcons=ap_abstract0_copy(pr->man_dcons, aux);
			ap_abstract0_free(pr->man_dcons,aux);

			rn->key.type=1;

			rn->key.segments[0]=look->segments[0];
			rn->key.segments[1]=look->segments[1];

			HASH_ADD(hh,r->udcons,key,keylen,rn);
			add_pattern_n2p(pr,r,look);


		}

	}//end folding universals
	//update lengths

	if(update_length){
		/* recalculate length of tdim[0] */

		ap_dim_t l0;
		/* update length for the merged segments */
		/* l(cdim[0]) = l(cdim[0]) + l(cdim[1]) + ... + l(cdim[last])*/

		/* for econs */
		expr =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
		ap_linexpr0_set_cst_scalar_int (expr, 0);
		for(i=0;i<size;i++){
			li = r->datadim + r->segmentdim + fst_node_tdim[i] ;
			ap_linexpr0_set_coeff_scalar_int (expr,li, 1);
		}

		ap_linexpr0_t *expr2 =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
		ap_linexpr0_set_cst_scalar_int (expr2, 0);
		for(i = 0;i < size; i++){
			li = r->datadim + r->segmentdim + snd_node_tdim[i] ;
			ap_linexpr0_set_coeff_scalar_int (expr2, li, 1);
		}
		/* for udcons */

		ap_linexpr0_t *expr_11y =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
		ap_linexpr0_set_cst_scalar_int (expr_11y, 0);
		for(i=0;i<size;i++){
			li = r->datadim + r->segmentdim + fst_node_tdim[i] ;
			ap_linexpr0_set_coeff_scalar_int (expr_11y, li, 1);
		}
		ap_linexpr0_t *expr_12y =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
		ap_linexpr0_set_cst_scalar_int (expr_12y, 0);
		for(i=0;i<size;i++){
			li = r->datadim + r->segmentdim + snd_node_tdim[i] ;
			ap_linexpr0_set_coeff_scalar_int (expr_12y, li, 1);
		}

		ap_linexpr0_t *expr_21y =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
		ap_linexpr0_set_cst_scalar_int (expr_21y, 0);
		for(i=0;i<size;i++){
			li = r->datadim + r->segmentdim + fst_node_tdim[i] ;
			ap_linexpr0_set_coeff_scalar_int (expr_21y, li, 1);
		}
		ap_linexpr0_t *expr_22y =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
		ap_linexpr0_set_cst_scalar_int (expr_22y, 0);
		for(i=0;i<size;i++){
			li = r->datadim + r->segmentdim + snd_node_tdim[i] ;
			ap_linexpr0_set_coeff_scalar_int (expr_22y, li, 1);
		}


		l01 = r->datadim + r->segmentdim + fst_node_tdim[0] ;
		l02 = r->datadim + r->segmentdim + snd_node_tdim[0] ;

		r->econs =  ap_abstract0_assign_linexpr (pr->man_dcons, true, r->econs,l01 , expr, NULL);
		ap_linexpr0_free(expr);

		r->econs =  ap_abstract0_assign_linexpr (pr->man_dcons, true, r->econs,l02 , expr2, NULL);
		ap_linexpr0_free(expr2);

		pattern_t *s;
		for(s=r->udcons; s!=NULL; s=s->hh.next){
			if(pr->PI[s->key.type].nr_y == 2){
				s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons,l01 , expr_21y, NULL);
				s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons,l02 , expr_22y, NULL);
			}
			else{
				s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons,l01 , expr_11y, NULL);
				s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons,l02 , expr_12y, NULL);
			}
		}
		ap_linexpr0_free(expr_21y);
		ap_linexpr0_free(expr_22y);
		ap_linexpr0_free(expr_11y);
		ap_linexpr0_free(expr_12y);

	}
#ifndef NDEBUG1
	printf("  \n\t fold without closure returns : \n");
	ucons_fprint(stdout,pr->man, r, NULL);
	printf("\n");
#endif
	return r;
}


ucons_t * split_with_pattern_P21(ucons_internal_t* pr, ucons_t *r,
		pattern_key_t *pattern_j, ap_dim_t n1, ap_dim_t n2){

#ifndef NDEBUG1
	fprintf(stdout,"\t split P21 : \n");
	ucons_fprint(stdout,pr->man,r, NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif
	pattern_key_t * pat_trasf;
	unsigned keylen;
	size_t u_seg,e_seg,nr_y;
	size_t y_pos,dy_pos , dn2;

	pattern_t *aux_find2, *aux, *sub_pattern_data;
	u_seg = pr->PI[pattern_j->type].u_seg;
	e_seg = pr->PI[pattern_j->type].e_seg;
	keylen = (u_seg+e_seg)*sizeof(size_t) + sizeof(pattern_key_t);

	nr_y = pr->PI[pattern_j->type].nr_y;
	ap_lincons0_array_t arr, arr1 ;
	HASH_FIND(hh,r->udcons,pattern_j,keylen,aux);
	if(aux){
		/* 1. NO TRANSFER over the same pattern where n1 is replaced to n2.
		 * 2. search for the pattern y1\in n1 and y1 = l[n1]
		 * 3. pattern if y1 \in n1 y2 \in n3 y1 = y2
		 * 			 search for the pattern y1\in n2 y2\in n3 and y1 = y2 + l[n1]
		 * 4. pattern if y1 \in n3 y2 \in n1 y1 = y2
		 * 			 search for the pattern y1\in n3 y2\in n2 and y1 + l[n1] = y2
		 *
		 * 5. set aux the constraint on pattern_j to bottom  */


		/*2.***********************************************************************************************************/
		/* y_n2 = y_n1 + l[x]
		 * y_n2 = y_n1
		 * y_n1 + l[x] = y_n2
		 * */
		pattern_key_t * pat_y_lx;
		size_t u_seg_y_lx = u_seg - 1;
		size_t e_seg_y_lx = e_seg + 1;
		size_t nr_y_y_lx = nr_y -1;
		size_t pos_y=0;

		checked_malloc(pat_y_lx, pattern_key_t, 1,
				sizeof(pattern_key_t) + (u_seg_y_lx + e_seg_y_lx)*sizeof(size_t), return NULL;);

		pat_y_lx->type = get_pattern_type(pr,u_seg_y_lx ,e_seg_y_lx ,nr_y_y_lx , pattern_1_lx);

		if (pat_y_lx->type!=100){
			ap_abstract0_t* dcons_y_lx;
			dcons_y_lx = ap_abstract0_copy(pr->man_dcons,aux->dcons);



			unsigned keylen_y_lx = (u_seg_y_lx + e_seg_y_lx )*sizeof(size_t) + sizeof(pattern_key_t);
			for(size_t i=0;i<e_seg;i++){
				pat_y_lx->segments[u_seg_y_lx+i] = pattern_j->segments[u_seg+i];
			}
			pat_y_lx->segments[u_seg_y_lx+e_seg] = n1;

			for (size_t ii = 2; ii < 1+e_seg_y_lx; ii++)
			{
				size_t  jj = 1;
				while (jj != ii && pat_y_lx->segments[jj] <= pat_y_lx->segments[ii])
					jj++;
				if (jj < ii)
				{
					size_t d = pat_y_lx->segments[ii];
					for (size_t  kk = ii - 1; kk >= jj; kk--)
						pat_y_lx->segments[kk + 1] = pat_y_lx->segments[kk];
					pat_y_lx->segments[jj] = d;
					//sort = false;
				}
			}

			/* pos_y is the position of the universal to eliminate */
			if(pattern_j->type == 1){
				if(pattern_j->segments[0]!=n1){
					pos_y = 1;
					pat_y_lx->segments[0] = pattern_j->segments[0];
				}
				else {
					pos_y = 0;
					pat_y_lx->segments[0] = pattern_j->segments[1];
				}
			}
			else if(pr->PI[pattern_j->type].kind == pattern_2_1_lx){
				pos_y = 1;
				pat_y_lx->segments[0] = pattern_j->segments[0];
			}
			else {
				pos_y = 0;
				pat_y_lx->segments[0] = pattern_j->segments[1];
			}

#ifndef NDEBUG2
			fprintf(stdout,"\t split introduces the pattern last: \n");
			fprintf(stdout,"\n");
			fprintf(stdout,"\t pat_y_lx->type = %zu",pat_y_lx->type);
			for(size_t lk = 0 ;lk< u_seg_y_lx + e_seg_y_lx; lk ++)
				fprintf(stdout,"\t pat_y_lx->segments[%zu]=%zu \n",lk,pat_y_lx->segments[lk] );
			fflush(stdout);
#endif

			/* y_pos = 1 and d(y_pos) = d(n2) */

			y_pos = r->datadim + 2*r->segmentdim + pos_y;
			dy_pos = r->datadim + 2*r->segmentdim + pos_y+nr_y;
			dn2 = r->datadim + n2;

			arr = ap_lincons0_array_make (2);
			/*
			 *  y_pos-1 == 0
			 */
			arr.p[0].constyp = AP_CONS_EQ;
			arr.p[0].linexpr0 =
					ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2*nr_y);
			arr.p[0].scalar = NULL;
			ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,y_pos, 1);
			ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);

			/* dy_pos - dn2 == 0
			 */
			arr.p[1].constyp = AP_CONS_EQ;
			arr.p[1].linexpr0 =
					ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2*nr_y);
			arr.p[1].scalar = NULL;
			ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, (dy_pos), -1);
			ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dn2, 1);
			ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);

			dcons_y_lx = ap_abstract0_meet_lincons_array (pr->man_dcons, true, dcons_y_lx, &arr);
			ap_lincons0_array_clear (&arr);

			/*  remove dimension y_pos, dy_pos */
			ap_dim_t *tdim;
			checked_malloc(tdim,ap_dim_t,2,sizeof(ap_dim_t),return NULL;);
			tdim[0] = y_pos;
			tdim[1] = dy_pos;
			dcons_y_lx = ap_abstract0_forget_array(pr->man_dcons, true, dcons_y_lx, tdim, 1, true);

			ap_dimchange_t* dimchange=ap_dimchange_alloc(2,0);
			dimchange->dim[0] = y_pos;
			dimchange->dim[1] = dy_pos;
			dcons_y_lx = ap_abstract0_remove_dimensions (pr->man_dcons, true, dcons_y_lx, dimchange);
			ap_dimchange_free(dimchange);

			pattern_t * pattern_y_lx;
			HASH_FIND(hh,r->udcons,pat_y_lx,keylen,pattern_y_lx);
			if(pattern_y_lx){
				pattern_y_lx->dcons = ap_abstract0_meet(pr->man_dcons,true,pattern_y_lx->dcons, dcons_y_lx);
			}
			else{
				checked_malloc(pattern_y_lx,pattern_t,1,sizeof(pattern_t)+(u_seg_y_lx+e_seg_y_lx)*sizeof(size_t),return NULL;);
				memset(pattern_y_lx, 0, sizeof(pattern_t)+(u_seg_y_lx+e_seg_y_lx)*sizeof(size_t));
				pattern_y_lx->key.type = pat_y_lx->type;
				for (size_t i=0 ; i<(u_seg_y_lx+e_seg_y_lx); i++)
					pattern_y_lx->key.segments[i] = pat_y_lx->segments[i];
				pattern_y_lx->dcons = ap_abstract0_copy(pr->man_dcons, dcons_y_lx);
				HASH_ADD(hh,r->udcons,key,keylen,pattern_y_lx);
				r=add_pattern_n2p(pr,r,pat_y_lx);
				ap_abstract0_free(pr->man_dcons,dcons_y_lx);
			};
		}// end if pat_y_lx->type !== 100 am gasit un pattern cu suma
		else{
#ifndef NDEBUG1
			fprintf(stdout,"\t pattern not found: \n");
			fprintf(stdout,"\n");
			fflush(stdout);
#endif
		}
		free(pat_y_lx);

#ifndef NDEBUG1
		fprintf(stdout,"\t middle split P21 : \n");
		ucons_fprint(stdout,pr->man,r, NULL);
		fprintf(stdout,"\n");
		fflush(stdout);
#endif

		/* 3+4. ***********************************************************************************************************/
		pattern_key_t * pat_2_1_lx;
		size_t u_seg_2_lx = u_seg;
		size_t e_seg_2_lx = e_seg + 1;
		size_t nr_y_2_lx = nr_y;

		checked_malloc(pat_2_1_lx, pattern_key_t, 1,
				sizeof(pattern_key_t) + (u_seg_2_lx + e_seg_2_lx)*sizeof(size_t), return NULL;);


		/* pos_y is the position of the universal to eliminate */
		/* pattern_j->segments[] == {n01,n02} */
		bool perm_u_seg = false;
		if(pattern_j->type == 1){
			if(pattern_j->segments[0]!=n1){
				/* y1\in n01, y2\in n4 */
				pos_y = 1;
				pat_2_1_lx->segments[0] = pattern_j->segments[0];
				pat_2_1_lx->segments[1] = n2;
			}
			else {
				/* y1\in n02, y2\in n4*/
				pos_y = 0;
				pat_2_1_lx->segments[0] = pattern_j->segments[1];
				pat_2_1_lx->segments[1] = n2;
				perm_u_seg = true;
			}
		}
		else if (pr->PI[pattern_j->type].kind == pattern_2_1_lx){
			if(pattern_j->segments[0]!=n1){
				/* y1\in n01 y2\in n4*/
				pos_y = 1;
				pat_2_1_lx->segments[0] = pattern_j->segments[0];
				pat_2_1_lx->segments[1] = n2;
			}
			else {
				fprintf(stdout,"Error !! bad pattern call ");
			}
		}
		else {
			/* y1\in n02, y2\in n4*/
			pos_y = 0;
			pat_2_1_lx->segments[0] = pattern_j->segments[1];
			pat_2_1_lx->segments[1] = n2;
			perm_u_seg = true;
		}


		/* pattern_2_1 mlx becomes pattern_2_1 lx because
		 * split n2 in y_n2 + l[n1] = y_n3 is :
		 * y_n4 + l[n2] + l[n1] = y_n3 but n4 > n3 so they are permuted
		 * */
		pat_2_1_lx->type = get_pattern_type(pr,u_seg_2_lx ,e_seg_2_lx ,nr_y_2_lx , pattern_2_1_lx);

		if(pat_2_1_lx->type != 100 ){
			ap_abstract0_t* dcons_y1y2_lx;
			dcons_y1y2_lx = ap_abstract0_copy(pr->man_dcons, aux->dcons);

			unsigned keylen_2_lx = (u_seg_2_lx + e_seg_2_lx )*sizeof(size_t) + sizeof(pattern_key_t);
			for(size_t i=0;i<e_seg;i++){
				pat_2_1_lx->segments[u_seg_2_lx+i] = pattern_j->segments[u_seg+i];
			}

			pat_2_1_lx->segments[u_seg_2_lx+e_seg_2_lx-1] = n1;

			for (size_t ii = u_seg_2_lx+1; ii < u_seg_2_lx+e_seg_2_lx; ii++)
			{
				size_t  jj = u_seg_2_lx;
				while (jj != ii && pat_2_1_lx->segments[jj] <= pat_2_1_lx->segments[ii])
					jj++;
				if (jj < ii)
				{
					size_t d = pat_2_1_lx->segments[ii];
					size_t tmp;
					size_t kk;
					for (size_t  tmp = ii ; tmp > jj; tmp--){
						kk = tmp -1;
						pat_2_1_lx->segments[kk + 1] = pat_2_1_lx->segments[kk];
					}
					pat_2_1_lx->segments[jj] = d;
					//sort = false;
				}
			}

			/* y_pos = y_pos+1 */
			y_pos = r->datadim + 2*r->segmentdim + pos_y;
			dy_pos = r->datadim + 2*r->segmentdim + pos_y+nr_y;
			dn2 = r->datadim + n2;

			ap_linexpr0_t* expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
					r->datadim + 2 * r->segmentdim + 2*nr_y);
			ap_linexpr0_set_coeff_scalar_int (expr, y_pos, 1);
			ap_linexpr0_set_cst_scalar_int (expr, 1);

			dcons_y1y2_lx = ap_abstract0_substitute_linexpr( pr->man_dcons, true, dcons_y1y2_lx, y_pos,
					expr, NULL );
			ap_linexpr0_free(expr);


			if(perm_u_seg){
				ap_dimperm_t * perm = ap_dimperm_alloc (r->datadim + 2*r->segmentdim + 2*nr_y);
				ap_dimperm_set_id(perm);
				if(pos_y == 0) {
					perm->dim[y_pos] = y_pos + 1;
					perm->dim[y_pos+1] = y_pos;
					perm->dim[dy_pos] = dy_pos + 1;
					perm->dim[dy_pos + 1] = dy_pos ;
				}
				else{
					perm->dim[y_pos] = y_pos - 1;
					perm->dim[y_pos-1] = y_pos;
					perm->dim[dy_pos] = dy_pos - 1;
					perm->dim[dy_pos - 1] = dy_pos ;
				}
				dcons_y1y2_lx = ap_abstract0_permute_dimensions (pr->man_dcons, true, dcons_y1y2_lx,perm);
				ap_dimperm_free(perm);
			}

#ifndef NDEBUG1
			fprintf(stdout,"\t dcons_y1y2_lx : \n");
			if(dcons_y1y2_lx) ap_abstract0_fprint(stdout,pr->man_dcons,dcons_y1y2_lx, NULL);
			fprintf(stdout,"\n");
			fflush(stdout);
#endif
			//	ap_abstract0_fprint(stdout,pr->man_dcons,dcons_y1y2_lx, NULL);
			pattern_t * pattern_2_lx;
			HASH_FIND(hh,r->udcons,pat_2_1_lx,keylen_2_lx,pattern_2_lx);
			if(pattern_2_lx){
				pattern_2_lx->dcons = ap_abstract0_meet(pr->man_dcons,true,pattern_2_lx->dcons, dcons_y1y2_lx);
				/* !y_pos >= 1  */

				arr = ap_lincons0_array_make (1);
				arr.p[0].constyp = AP_CONS_SUPEQ;
				arr.p[0].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2*nr_y);
				arr.p[0].scalar = NULL;
				ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);

				ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
						(1+r->datadim + 2*r->segmentdim), 1);

				pattern_2_lx->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, pattern_2_lx->dcons, &arr);

				ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
						(1+r->datadim + 2*r->segmentdim), 0);
				ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
						(r->datadim + 2*r->segmentdim), 1);
				pattern_2_lx->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, pattern_2_lx->dcons, &arr);
				ap_lincons0_array_clear (&arr);
			}
			else{
				checked_malloc(pattern_2_lx,pattern_t,1,
						sizeof(pattern_t)+(u_seg_2_lx+e_seg_2_lx)*sizeof(size_t),return NULL;);
				memset(pattern_2_lx, 0, sizeof(pattern_t)+(u_seg_2_lx+e_seg_2_lx)*sizeof(size_t));
				pattern_2_lx->key.type = pat_2_1_lx->type;
				for (size_t i=0 ; i<(u_seg_2_lx+e_seg_2_lx); i++)
					pattern_2_lx->key.segments[i] = pat_2_1_lx->segments[i];
				pattern_2_lx->dcons = ap_abstract0_copy(pr->man_dcons, dcons_y1y2_lx);
				/* !y_pos >= 1  */

				arr = ap_lincons0_array_make (1);
				arr.p[0].constyp = AP_CONS_SUPEQ;
				arr.p[0].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2*nr_y);
				arr.p[0].scalar = NULL;
				ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);

				ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
						(1+r->datadim + 2*r->segmentdim), 1);

				pattern_2_lx->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, pattern_2_lx->dcons, &arr);

				ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
						(1+r->datadim + 2*r->segmentdim), 0);
				ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
						(r->datadim + 2*r->segmentdim), 1);
				pattern_2_lx->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, pattern_2_lx->dcons, &arr);
				ap_lincons0_array_clear (&arr);

				HASH_ADD(hh,r->udcons,key,keylen_2_lx,pattern_2_lx);
				add_pattern_n2p(pr,r,pat_2_1_lx);
				ap_abstract0_free(pr->man_dcons,dcons_y1y2_lx);
			};
		}//end if pat_2_1_lx->type!=100 am gasit pattern pe care sa pot exprima prelungirea propietatii
		else{
#ifndef NDEBUG1
			fprintf(stdout,"\t pattern not found: \n");
			fprintf(stdout,"\n");
			fflush(stdout);
#endif
		}
		free(pat_2_1_lx);

		/* set old constraint to bottom */
		ap_abstract0_free(pr->man_dcons,aux->dcons);
		aux->dcons = NULL;

		remove_pattern_n2p(pr,r,&aux->key);
		//		pattern_t *t = s;
		//		s = s->hh.next;
		HASH_DEL(r->udcons,aux);
		//		free(t);
		free(aux);


	}
	else{
#ifndef NDEBUG1
		fprintf(stdout,"\t pattern not found: \n");
		fprintf(stdout,"\n");
		fflush(stdout);
#endif
	}

#ifndef NDEBUG1
	fprintf(stdout,"\t split P21 returns : \n");
	ucons_fprint(stdout,pr->man,r, NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

	return r;
}

bool test_equal_length(ap_manager_t * man, ap_abstract0_t * dcons, size_t intdim, size_t realdim, ap_dim_t n, ap_dim_t m){

	ap_lincons0_t cons;
	ap_linexpr0_t *linexpr;

	bool eq=false;
	linexpr = ap_linexpr0_alloc ( AP_LINEXPR_DENSE, intdim + 2 * realdim);
	ap_linexpr0_set_coeff_scalar_int (linexpr, intdim + realdim + n , 1);
	ap_linexpr0_set_coeff_scalar_int (linexpr, intdim + realdim + m , -1);
	ap_linexpr0_set_cst_scalar_int (linexpr, 0);
	cons = ap_lincons0_make(AP_CONS_EQ, linexpr, NULL);
	//	ap_lincons0_fprint(stdout,&cons,NULL);

	eq = ap_abstract0_sat_lincons(man,dcons,&cons);

	ap_lincons0_clear(&cons);
	linexpr=NULL;

	return eq;
}

bool test_singleton(ap_manager_t * man, ap_abstract0_t * dcons, size_t intdim, size_t realdim, ap_dim_t n){

	ap_lincons0_t cons;
	ap_linexpr0_t *linexpr;

	bool singl=false;
	linexpr = ap_linexpr0_alloc ( AP_LINEXPR_DENSE, intdim + 2 * realdim);
	ap_linexpr0_set_coeff_scalar_int (linexpr,intdim + realdim + n , 1);
	ap_linexpr0_set_cst_scalar_int (linexpr, -1);
	cons = ap_lincons0_make(AP_CONS_EQ, linexpr, NULL);

	//	ap_lincons0_fprint(stdout,&cons,NULL);

	singl = ap_abstract0_sat_lincons(man,dcons,&cons);

	ap_lincons0_clear(&cons);
	linexpr=NULL;

	return singl;
}

bool test_l_leq_v(ap_manager_t * man, ap_abstract0_t * dcons,
		size_t intdim, size_t realdim, ap_dim_t n, int v){

	ap_lincons0_t cons;
	ap_linexpr0_t *linexpr;

	bool singl=false;
	linexpr = ap_linexpr0_alloc ( AP_LINEXPR_DENSE, intdim + 2 * realdim);
	ap_linexpr0_set_coeff_scalar_int (linexpr,intdim + realdim + n , -1);
	ap_linexpr0_set_cst_scalar_int (linexpr, v);
	cons = ap_lincons0_make(AP_CONS_SUPEQ, linexpr, NULL);
	//	ap_lincons0_fprint(stdout,&cons,NULL);

	singl = ap_abstract0_sat_lincons(man,dcons,&cons);

	ap_lincons0_clear(&cons);
	linexpr=NULL;

	return singl;
}

bool test_g_leq_v(ap_manager_t * man, ap_abstract0_t * dcons,
		size_t intdim, size_t realdim, ap_dim_t n, int v){

	ap_lincons0_t cons;
	ap_linexpr0_t *linexpr;

	bool singl=false;
	linexpr = ap_linexpr0_alloc ( AP_LINEXPR_DENSE, intdim + 2 * realdim);
	ap_linexpr0_set_coeff_scalar_int (linexpr,intdim + realdim + n , 1);
	ap_linexpr0_set_cst_scalar_int (linexpr, -v);
	cons = ap_lincons0_make(AP_CONS_SUPEQ, linexpr, NULL);
	//	ap_lincons0_fprint(stdout,&cons,NULL);

	singl = ap_abstract0_sat_lincons(man,dcons,&cons);

	ap_lincons0_clear(&cons);
	linexpr=NULL;

	return singl;
}


ap_dim_t * equal_lengths(ap_manager_t * man, ap_abstract0_t * dcons, size_t intdim,
		size_t realdim, ap_dim_t n, ap_dim_t *tdim, size_t size_tdim, size_t* size){

	ap_dim_t *len = NULL;
	size_t size_ps=realdim-1-size_tdim;
	ap_dim_t *ps = (ap_dim_t*)malloc((realdim-1-size_tdim)*sizeof(ap_dim_t));

	size_t kk = 0;
	for(size_t i=1;i<realdim;i++){
		bool found = false;
		for (size_t j=0; j<size_tdim && !found; j++)
			if (tdim[j] - intdim == i) found = true;
		if(!found)
			ps[kk++] = i;
	}
	size_t i_len=0;
	for(size_t i = 0; i<size_ps; i++)
		if(test_equal_length(man,dcons,intdim,realdim,n,ps[i])){
			len = (ap_dim_t*) malloc (sizeof(ap_dim_t));
			len[i_len] = ps[i];
			++i_len;
			*size = i_len;
			return len;
		}
	//todo multiple choices for length: ex search for l[n] = l[n1] + l[n2]
	return NULL;
}

ap_lincons0_t cons_dy0(size_t intdim, size_t realdim, ap_dim_t dy0, ap_dim_t n){

	ap_lincons0_t cons;
	ap_linexpr0_t *linexpr;

	linexpr = ap_linexpr0_alloc ( AP_LINEXPR_DENSE, intdim + 2 * realdim + 2);
	ap_linexpr0_set_coeff_scalar_int (linexpr,intdim + 2 * realdim + 1 , 1);
	ap_linexpr0_set_coeff_scalar_int (linexpr,n , -1);
	cons = ap_lincons0_make (AP_CONS_SUPEQ, linexpr, NULL);

	ap_lincons0_fprint(stdout,&cons,NULL);
	return cons;

}

ap_lincons0_t cons_k_m(size_t intdim, size_t realdim, ap_dim_t k, ap_dim_t m){

	ap_lincons0_t cons;
	ap_linexpr0_t *linexpr;

	linexpr = ap_linexpr0_alloc ( AP_LINEXPR_DENSE, intdim + 2 * realdim );
	ap_linexpr0_set_coeff_scalar_int (linexpr,k , -1);
	ap_linexpr0_set_coeff_scalar_int (linexpr,m , 1);
	cons = ap_lincons0_make (AP_CONS_SUPEQ, linexpr, NULL);

	//	ap_lincons0_fprint(stdout,&cons,NULL);
	return cons;

}

ap_lincons0_t cons_dy0_dy1(size_t intdim, size_t realdim, ap_dim_t dy0, ap_dim_t dy1){

	ap_lincons0_t cons;
	ap_linexpr0_t *linexpr;

	linexpr = ap_linexpr0_alloc ( AP_LINEXPR_DENSE, intdim + 2 * realdim + 4);
	ap_linexpr0_set_coeff_scalar_int (linexpr,intdim + 2 * realdim + 2 , -1);
	ap_linexpr0_set_coeff_scalar_int (linexpr,intdim + 2 * realdim + 3 , 1);
	cons = ap_lincons0_make (AP_CONS_SUPEQ, linexpr, NULL);

	//	ap_lincons0_fprint(stdout,&cons,NULL);

	return cons;

}

