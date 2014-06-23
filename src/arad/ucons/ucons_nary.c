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
 * N-ary ucons functions: meet, join, widening, narrowing & related.
 */

#include "ucons.h"
#include "ucons_internal.h"
#include "ucons_fun.h"
#include "apron2shape.h"
#include "shape_macros.h"
#include "ap_generic.h"

/* ============================================================ */
/* Meet and Join */
/* ============================================================ */

/* TODO: priority 3 */




/* TODO: priority 3 */
ucons_t *
ucons_meet (ap_manager_t * man, bool destructive, ucons_t * a1, ucons_t * a2)
{
	ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_MEET, 0);

	if (ucons_is_bottom (man, a1) || ucons_is_bottom (man, a2)){
		if(destructive){
			if(a1) ucons_free_internal(pr,a1);
			if(a2) ucons_free_internal(pr,a2);
		}
		return NULL;
	}



	arg_assert (a1->datadim == a2->datadim
			&& a1->segmentdim == a2->segmentdim, return NULL;);

#ifndef NDEBUG
	fprintf(stdout,"\n \t meet a1 \n ");
	ucons_fprint(stdout,man,a1,NULL);
	fprintf(stdout,"\n \t meet  a2 \n ");
	ucons_fprint(stdout,man,a2,NULL);
	fprintf(stdout,"\n");
#endif

	ucons_t * r  = ucons_alloc_internal(pr, a1->datadim,a1->segmentdim);

	r->econs = ap_abstract0_meet (pr->man_dcons, false, a1->econs, a2->econs);


	pattern_t *sa1;
	pattern_t *sa2;

	pattern_t *ru=NULL; //universal constraint to add to r
	size_t u_seg;
	size_t e_seg;


	size_t i;
	unsigned keylen;


	for(sa1=a1->udcons; sa1 != NULL; sa1=sa1->hh.next) {

		u_seg = pr->PI[sa1->key.type].u_seg; // number of universal segments
		e_seg = pr->PI[sa1->key.type].e_seg; // number of existential segments
		keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);



		HASH_FIND(hh,a2->udcons,&sa1->key,keylen,sa2);

		HASH_FIND(hh,r->udcons,&sa1->key,keylen,ru);

		//bool pat_sat = test_pattern_sat(pr, r, &sa1->key, 0);
		ap_dim_t n = sa1->key.segments[0];

		/* pat_sat == false <=> singula valoare pentru lungimea lui l este 1 in r->econs */
//		bool pat_sat = !(test_singleton(pr->man_dcons, r->econs, r->datadim, r->segmentdim ,n)&&
//				(!test_l_leq_2(pr->man_dcons, r->econs, r->datadim, r->segmentdim ,n)));

		bool pat_sat = !test_singleton(pr->man_dcons, r->econs, r->datadim, r->segmentdim ,n);
				//test_l_leq_2(pr->man_dcons, r->econs, r->datadim, r->segmentdim ,n);

		//checked_malloc(ru,pattern_t,1,sizeof(pattern_t)+u_seg*sizeof(size_t),return NULL;);
		//memset(ru, 0, sizeof(pattern_t)+u_seg*sizeof(size_t));

		ap_abstract0_t *dcons1, *dcons2;
		ap_dimchange_t dimadd;
		size_t nr_y, i;

		if(sa2){
			if (ru) {
				ru->dcons=ap_abstract0_meet(pr->man_dcons,false,sa1->dcons,sa2->dcons);
			}
			else{
				//if(pat_sat){
				pattern_t *rn;

				u_seg = pr->PI[sa1->key.type].u_seg;
				e_seg = pr->PI[sa1->key.type].e_seg;
				keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

				checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
				memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));

				rn->dcons=ap_abstract0_meet(pr->man_dcons,false,sa1->dcons,sa2->dcons);

				rn->key.type=sa1->key.type;

				for(size_t i=0; i<(u_seg+e_seg); i++)
					rn->key.segments[i]=sa1->key.segments[i];


				HASH_ADD(hh,r->udcons,key,keylen,rn);
				add_pattern_n2p(pr,r,&rn->key);
				//}

			}
		}
		else {
			//HASH_FIND(hh,r->udcons,&sa1->key,keylen,ru);
			if(pat_sat){

				if (ru) {
					ru->dcons=ap_abstract0_copy(pr->man_dcons,sa1->dcons);
				}
				else{
					pattern_t *rn;

					u_seg = pr->PI[sa1->key.type].u_seg;
					e_seg = pr->PI[sa1->key.type].e_seg;
					keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

					checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
					memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));

					rn->dcons=ap_abstract0_copy(pr->man_dcons,sa1->dcons);

					rn->key.type=sa1->key.type;

					for(size_t i=0; i<(u_seg+e_seg); i++)
						rn->key.segments[i]=sa1->key.segments[i];


					HASH_ADD(hh,r->udcons,key,keylen,rn);
					add_pattern_n2p(pr,r,&rn->key);
				}
			}
			else{
				if(ru){
					if(ru->dcons) ap_abstract0_free(pr->man_dcons,ru->dcons);
					ru = NULL;
				}
			}
		}
	}
	sa1 = NULL;

	for(sa2=a2->udcons; sa2 != NULL; sa2=sa2->hh.next) {

		u_seg = pr->PI[sa2->key.type].u_seg; // number of universal segments
		e_seg = pr->PI[sa2->key.type].e_seg; // number of existential segments
		keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);


		HASH_FIND(hh,a1->udcons,&sa2->key,keylen,sa1);

		HASH_FIND(hh,r->udcons,&sa2->key,keylen,ru);


		ap_dim_t n = sa2->key.segments[0];
		//		bool pat_sat = !(test_singleton(pr->man_dcons, r->econs, r->datadim, r->segmentdim ,n)&&
		//						(!test_l_leq_2(pr->man_dcons, r->econs, r->datadim, r->segmentdim ,n)));
		bool pat_sat = !test_singleton(pr->man_dcons, r->econs, r->datadim, r->segmentdim ,n);
				//test_l_leq_2(pr->man_dcons, r->econs, r->datadim, r->segmentdim ,n);


		//checked_malloc(ru,pattern_t,1,sizeof(pattern_t)+u_seg*sizeof(size_t),return NULL;);
		//memset(ru, 0, sizeof(pattern_t)+u_seg*sizeof(size_t));

		ap_abstract0_t *dcons1, *dcons2;
		ap_dimchange_t dimadd;
		size_t nr_y, i;

		if(!sa1){
			if (ru) {
				if(pat_sat)
					ru->dcons=ap_abstract0_copy(pr->man_dcons,sa2->dcons);
				else{
					ap_abstract0_free(pr->man_dcons,ru->dcons);
					ru->dcons=NULL;
				}
			}
			else{
				if(pat_sat){
					pattern_t *rn;

					u_seg = pr->PI[sa2->key.type].u_seg;
					e_seg = pr->PI[sa2->key.type].e_seg;

					keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

					checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
					memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));

					rn->dcons=ap_abstract0_copy(pr->man_dcons,sa2->dcons);

					rn->key.type=sa2->key.type;

					for(size_t i=0; i<(u_seg+e_seg); i++)
						rn->key.segments[i]=sa2->key.segments[i];


					HASH_ADD(hh,r->udcons,key,keylen,rn);
					add_pattern_n2p(pr,r,&rn->key);
				}
			}
		}
	}

	pattern_t * sr=r->udcons;
	while(sr) {
		if(sr->dcons == NULL){
			fprintf(stdout," \n Null in meet \n ");
			fflush(stdout);
			pattern_t * sr2 = sr->hh.next;
			remove_pattern_n2p(pr,r,&sr->key);
			HASH_DEL(r->udcons,sr);
			sr = sr2;
		}
		else {
			ap_abstract0_t * auxn = ap_abstract0_copy(pr->man_dcons,r->econs);
			ap_dimchange_t dimadd;
			if(pr->PI[sr->key.type].nr_y == 1){
				ap_dimchange_init (&dimadd, 2, 0);
				dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
				dimadd.dim[0] = r->datadim + 2*r->segmentdim ;
				dimadd.dim[1] = r->datadim + 2*r->segmentdim ;
				ap_abstract0_t * aux = ap_abstract0_add_dimensions(pr->man_dcons,true,auxn,&dimadd,false);
				ap_dimchange_clear (&dimadd);
				sr->dcons = ap_abstract0_meet(pr->man_dcons,true,aux,sr->dcons);
			}else{
				ap_dimchange_init (&dimadd, 4, 0);
				dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
				dimadd.dim[0] = r->datadim + 2*r->segmentdim ;
				dimadd.dim[1] = r->datadim + 2*r->segmentdim ;
				dimadd.dim[2] = r->datadim + 2*r->segmentdim ;
				dimadd.dim[3] = r->datadim + 2*r->segmentdim ;
				ap_abstract0_t * aux = ap_abstract0_add_dimensions(pr->man_dcons,true,auxn,&dimadd,false);
				ap_dimchange_clear (&dimadd);
				sr->dcons = ap_abstract0_meet(pr->man_dcons,true,aux,sr->dcons);
			}

			sr = sr->hh.next;
		}
	}

	if (destructive)
	{
		ucons_free_internal (pr, a1);
		ucons_free_internal (pr, a2);
	}

#ifndef NDEBUG
	fprintf(stdout,"\n \t meet returns  \n ");
	ucons_fprint(stdout,man,r,NULL);
	fprintf(stdout,"\n");
#endif

	return r;
}




//ucons_t *
//ucons_meet (ap_manager_t * man, bool destructive, ucons_t * a1, ucons_t * a2)
//{
//	ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_MEET, 0);
//
//	if (ucons_is_bottom (man, a1) && ucons_is_bottom (man, a2))
//		return (destructive) ? a1 : ucons_bottom (man, a1->datadim, a1->segmentdim);
//	if (ucons_is_bottom (man, a1))
//		return (destructive) ? a2 : ucons_copy (man, a2);
//	if (ucons_is_bottom (man, a2))
//		return (destructive) ? a1 : ucons_copy (man, a1);
//
//
//	arg_assert (a1->datadim == a2->datadim
//			&& a1->segmentdim == a2->segmentdim, return NULL;);
//
////#ifndef NDEBUG1
////		printf("  \n\t meet a1 : \n");
////		ucons_fprint(stdout,pr->man,a1,NULL);
////		printf("\n with a2 \n");
////		ucons_fprint(stdout,pr->man,a2,NULL);
////		printf("\n");
////#endif
//
//	ucons_t * r  = ucons_alloc_internal(pr,  a1->datadim,a1->segmentdim );
//
//	r->econs = ap_abstract0_meet (pr->man_dcons, destructive, a1->econs, a2->econs);
//
//	r->n2p  = pattern_key_set_copy (pr, a1->n2p ,a1->segmentdim);
//
//	pattern_t *sa1;
//	pattern_t *sa2;
//
//	pattern_t *ru=NULL; //universal constraint to add to r
//	size_t u_seg;
//
//	size_t i;
//	unsigned keylen;
//	for(sa1=a1->udcons; sa1 != NULL; sa1=sa1->hh.next) {
//
//		u_seg = pr->PI[sa1->key.type].u_seg; // number of universal segments
//		keylen = u_seg*sizeof(size_t)+sizeof(pattern_key_t);
//
//
//
//		HASH_FIND(hh,a2->udcons,&sa1->key,keylen,sa2);
//
//		HASH_FIND(hh,r->udcons,&sa1->key,keylen,ru);
//
//		if(sa2 && ru ){
//
//			ru->dcons=ap_abstract0_meet(pr->man_dcons,destructive,sa1->dcons,sa2->dcons);
//
//		}
//		else
//		{
//			ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
//					"not implemented for diffrent ucons object over differnt set of patterns ");
//		}
//
//	}
//	if (destructive)
//	{
//		ucons_free_internal (pr, a1);
//		ucons_free_internal (pr, a2);
//	}
////#ifndef NDEBUG1
////		printf("  \n\t meet returns : \n");
////		ucons_fprint(stdout,pr->man,r,NULL);
////		printf("\n");
////#endif
//	return r;
//
//}

/* TODO: priority 3 */
ucons_t *
ucons_join (ap_manager_t * man, bool destructive, ucons_t * a1, ucons_t * a2)
{
	ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_JOIN, 0);

	if (ucons_is_bottom (man, a1) && ucons_is_bottom (man, a2))
		return (destructive) ? a1 : ucons_bottom (man, a1->datadim, a1->segmentdim);
	if (ucons_is_bottom (man, a1))
		return (destructive) ? a2 : ucons_copy (man, a2);
	if (ucons_is_bottom (man, a2))
		return (destructive) ? a1 : ucons_copy (man, a1);


	arg_assert (a1->datadim == a2->datadim
			&& a1->segmentdim == a2->segmentdim, return NULL;);

	#ifndef NDEBUG1
		fprintf(stdout,"\n \t join a1 \n ");
		ucons_fprint(stdout,man,a1,NULL);
		fprintf(stdout,"\n \t with  a2 \n ");
		ucons_fprint(stdout,man,a2,NULL);
		fprintf(stdout,"\n");
	#endif

	ucons_t * r  = ucons_alloc_internal(pr, a1->datadim,a1->segmentdim);

	r->econs = ap_abstract0_join (pr->man_dcons, destructive, a1->econs, a2->econs);


	pattern_t *sa1;
	pattern_t *sa2;

	pattern_t *ru=NULL; //universal constraint to add to r
	size_t u_seg;
	size_t e_seg,nr_y;


	size_t i;
	unsigned keylen;


	ap_dimchange_t dimadd;
	ap_dimchange_init (&dimadd, 2, 0);
	dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
	dimadd.dim[0] = r->datadim + 2*r->segmentdim;
	dimadd.dim[1] = r->datadim + 2*r->segmentdim;

	ap_abstract0_t * aux_join = ap_abstract0_add_dimensions(pr->man_dcons,false,a2->econs,&dimadd,false);
	ap_abstract0_t * aux_join_a1 = ap_abstract0_add_dimensions(pr->man_dcons,false,a1->econs,&dimadd,false);

	ap_dimchange_clear (&dimadd);

	ap_dimchange_init (&dimadd, 4, 0);
	dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
	dimadd.dim[0] = r->datadim + 2*r->segmentdim;
	dimadd.dim[1] = r->datadim + 2*r->segmentdim;
	dimadd.dim[2] = r->datadim + 2*r->segmentdim;
	dimadd.dim[3] = r->datadim + 2*r->segmentdim;

	ap_abstract0_t * y2aux_join = ap_abstract0_add_dimensions(pr->man_dcons,false,a2->econs,&dimadd,false);
	ap_abstract0_t * y2aux_join_a1 = ap_abstract0_add_dimensions(pr->man_dcons,false,a1->econs,&dimadd,false);

	ap_dimchange_clear (&dimadd);

	for(sa1=a1->udcons; sa1 != NULL; sa1=sa1->hh.next) {

		u_seg = pr->PI[sa1->key.type].u_seg; // number of universal segments
		e_seg = pr->PI[sa1->key.type].e_seg; // number of existential segments
		nr_y = pr->PI[sa1->key.type].nr_y; // number of existential segments

		keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);


		ap_dim_t n = sa1->key.segments[0];
		bool pats1_sat_a2 =!test_singleton(pr->man_dcons, a2->econs, a2->datadim, a2->segmentdim ,n);
				//test_l_leq_2(pr->man_dcons, a2->econs, r->datadim, r->segmentdim ,n);

		HASH_FIND(hh,a2->udcons,&sa1->key,keylen,sa2);

		HASH_FIND(hh,r->udcons,&sa1->key,keylen,ru);

		//checked_malloc(ru,pattern_t,1,sizeof(pattern_t)+u_seg*sizeof(size_t),return NULL;);
		//memset(ru, 0, sizeof(pattern_t)+u_seg*sizeof(size_t));

		ap_abstract0_t *dcons1, *dcons2;
		ap_dimchange_t dimadd;
		size_t i;

		if(sa2){
			if (ru) {
				ru->dcons=ap_abstract0_join(pr->man_dcons,destructive,sa1->dcons,sa2->dcons);
			}
			else{

				pattern_t *rn;

				u_seg = pr->PI[sa1->key.type].u_seg;
				e_seg = pr->PI[sa1->key.type].e_seg;
				keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

				checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
				memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));

				rn->dcons=ap_abstract0_join(pr->man_dcons,destructive,sa1->dcons,sa2->dcons);

		
				ap_lincons0_array_t arr = ap_lincons0_array_make (1);

				arr.p[0].constyp = AP_CONS_SUPEQ;
				arr.p[0].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2*r->segmentdim + 2*nr_y );
				arr.p[0].scalar = NULL;
				bool ord, ord1, ord2;
				ord = false; ord1 = false; ord2 = false;

				if(nr_y==2){
					// d(y1) - d(y2) >=0
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 2, 1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 3, -1);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

					ord1 = ap_abstract0_sat_lincons(pr->man_dcons,sa1->dcons,&arr.p[0]);
					ord2 = ap_abstract0_sat_lincons(pr->man_dcons,sa2->dcons,&arr.p[0]);

					if(ord1 && ord2)
						rn->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true,rn->dcons, &arr);

					// d(y2) - d(y1) >=0
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 2, -1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 3, 1);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

					ord1 = ap_abstract0_sat_lincons(pr->man_dcons,sa1->dcons,&arr.p[0]);
					ord2 = ap_abstract0_sat_lincons(pr->man_dcons,sa2->dcons,&arr.p[0]);

					if(ord1 && ord2)
							rn->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true,rn->dcons, &arr);

					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 2, 0);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 3, 0);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

				}
				ap_dim_t dctr, dy1, dy2;
				for(size_t ctr = 1 ; ctr < a1->segmentdim; ctr ++){
					dctr = a1->datadim + ctr;
					if(nr_y == 1)
						dy1 =  r->datadim + 2*r->segmentdim + 1;
					else
						dy1 =  r->datadim + 2*r->segmentdim + 2;

					dy2 = r->datadim + 2*r->segmentdim + 3;

					 //d(ctr) - d(y1) >= 0
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy1, -1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dctr, 1);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

					ord1 = ap_abstract0_sat_lincons(pr->man_dcons,sa1->dcons,&arr.p[0]);
					ord2 = ap_abstract0_sat_lincons(pr->man_dcons,sa2->dcons,&arr.p[0]);

					if(ord1 && ord2)
							rn->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, rn->dcons, &arr);

					// d(y1) - d(ctr) >= 0
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy1, 1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dctr, -1);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

					ord1 = ap_abstract0_sat_lincons(pr->man_dcons,sa1->dcons,&arr.p[0]);
					ord2 = ap_abstract0_sat_lincons(pr->man_dcons,sa2->dcons,&arr.p[0]);

					if(ord1 && ord2)
						rn->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, rn->dcons, &arr);

					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy1, 0);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dctr, 0);

					if(nr_y == 2){
					//	 d(ctr) - d(y2) >= 0
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy2, -1);
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dctr, 1);
						ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

						ord1 = ap_abstract0_sat_lincons(pr->man_dcons,sa1->dcons,&arr.p[0]);
						ord2 = ap_abstract0_sat_lincons(pr->man_dcons,sa2->dcons,&arr.p[0]);

						if(ord1 && ord2)
							rn->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, rn->dcons, &arr);

					//	 d(y2) - d(ctr) >= 0
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy2, 1);
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dctr, -1);
						ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

						ord1 = ap_abstract0_sat_lincons(pr->man_dcons,sa1->dcons,&arr.p[0]);
						ord2 = ap_abstract0_sat_lincons(pr->man_dcons,sa2->dcons,&arr.p[0]);

						if(ord1 && ord2)
							rn->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, rn->dcons, &arr);

						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy2, 0);
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dctr, 0);
					}

				}


				ap_lincons0_array_clear (&arr);
				



				rn->key.type=sa1->key.type;

				for(size_t i=0; i<(u_seg+e_seg); i++)
					rn->key.segments[i]=sa1->key.segments[i];


				HASH_ADD(hh,r->udcons,key,keylen,rn);
				add_pattern_n2p(pr,r,&rn->key);

			}
		}
		else {
			//HASH_FIND(hh,r->udcons,&sa1->key,keylen,ru);
			if(!pats1_sat_a2){
				if (ru) {
					ru->dcons=ap_abstract0_copy(pr->man_dcons,sa1->dcons);
				}
				else{

					pattern_t *rn;

					u_seg = pr->PI[sa1->key.type].u_seg;
					e_seg = pr->PI[sa1->key.type].e_seg;
					keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

					checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
					memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));

					rn->dcons=ap_abstract0_copy(pr->man_dcons,sa1->dcons);

					rn->key.type=sa1->key.type;

					for(size_t i=0; i<(u_seg+e_seg); i++)
						rn->key.segments[i]=sa1->key.segments[i];


					HASH_ADD(hh,r->udcons,key,keylen,rn);
					add_pattern_n2p(pr,r,&rn->key);

				}
			}//end pattern not sat
			else{

#ifndef NDEBUG1
				fprintf(stdout,"\n \t join NOU  \n ");
				fprintf(stdout,"\n");
				fflush(stdout);
#endif


				if(ru){
					if(nr_y == 1)
						ru->dcons = ap_abstract0_join(pr->man_dcons,false,aux_join,sa1->dcons);
					else
						ru->dcons = ap_abstract0_join(pr->man_dcons,false,y2aux_join,sa1->dcons);
				}
				else{
					pattern_t *rn;

					u_seg = pr->PI[sa1->key.type].u_seg;
					e_seg = pr->PI[sa1->key.type].e_seg;
					keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

					checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
					memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));

					if(nr_y == 1)
						rn->dcons=ap_abstract0_join(pr->man_dcons,false,aux_join,sa1->dcons);
					else
						rn->dcons=ap_abstract0_join(pr->man_dcons,false,y2aux_join,sa1->dcons);

					rn->key.type=sa1->key.type;

					for(size_t i=0; i<(u_seg+e_seg); i++)
						rn->key.segments[i]=sa1->key.segments[i];


					HASH_ADD(hh,r->udcons,key,keylen,rn);
					add_pattern_n2p(pr,r,&rn->key);
				}
			}//end pattern sat
		}
	}

	sa1 = NULL;


	for(sa2=a2->udcons; sa2 != NULL; sa2=sa2->hh.next) {

		u_seg = pr->PI[sa2->key.type].u_seg; // number of universal segments
		e_seg = pr->PI[sa2->key.type].e_seg; // number of existential segments
		nr_y = pr->PI[sa2->key.type].nr_y;

		keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);


		HASH_FIND(hh,a1->udcons,&sa2->key,keylen,sa1);

		HASH_FIND(hh,r->udcons,&sa2->key,keylen,ru);

		//checked_malloc(ru,pattern_t,1,sizeof(pattern_t)+u_seg*sizeof(size_t),return NULL;);
		//memset(ru, 0, sizeof(pattern_t)+u_seg*sizeof(size_t));

		ap_dim_t n = sa2->key.segments[0];
		bool pats2_sat_a1 =!test_singleton(pr->man_dcons, a1->econs, a1->datadim, a1->segmentdim ,n);
				//test_l_leq_2(pr->man_dcons, a1->econs, r->datadim, r->segmentdim ,n);

		ap_abstract0_t *dcons1, *dcons2;
		ap_dimchange_t dimadd;
		size_t i;

		if(!sa1){
			if(!pats2_sat_a1){
				if (ru) {
					ru->dcons=ap_abstract0_copy(pr->man_dcons,sa2->dcons);
				}
				else{

					pattern_t *rn;

					u_seg = pr->PI[sa2->key.type].u_seg;
					e_seg = pr->PI[sa2->key.type].e_seg;

					keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

					checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
					memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));

					rn->dcons=ap_abstract0_copy(pr->man_dcons,sa2->dcons);

					rn->key.type=sa2->key.type;

					for(size_t i=0; i<(u_seg+e_seg); i++)
						rn->key.segments[i]=sa2->key.segments[i];


					HASH_ADD(hh,r->udcons,key,keylen,rn);
					add_pattern_n2p(pr,r,&rn->key);

				}
			}//pattern not sat
			else{

#ifndef NDEBUG1
				fprintf(stdout,"\n \t join2 NOU  \n ");
				fprintf(stdout,"\n");
				fflush(stdout);
#endif
				if(ru){
					if(nr_y == 1)
						ru->dcons = ap_abstract0_join(pr->man_dcons,false,aux_join_a1,sa2->dcons);
					else
						ru->dcons = ap_abstract0_join(pr->man_dcons,false,y2aux_join_a1,sa2->dcons);
				}
				else{
					pattern_t *rn;

					u_seg = pr->PI[sa2->key.type].u_seg;
					e_seg = pr->PI[sa2->key.type].e_seg;
					keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

					checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
					memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));

					if(nr_y == 1)
						rn->dcons=ap_abstract0_join(pr->man_dcons,false,aux_join_a1,sa2->dcons);
					else
						rn->dcons=ap_abstract0_join(pr->man_dcons,false,y2aux_join_a1,sa2->dcons);

					rn->key.type=sa2->key.type;

					for(size_t i=0; i<(u_seg+e_seg); i++)
						rn->key.segments[i]=sa2->key.segments[i];


					HASH_ADD(hh,r->udcons,key,keylen,rn);
					add_pattern_n2p(pr,r,&rn->key);
				}
			}//end pattern has instances
		}
	}
	if(aux_join)ap_abstract0_free(pr->man_dcons,aux_join);
	if(aux_join_a1)ap_abstract0_free(pr->man_dcons,aux_join_a1);

	if (destructive)
	{
		ucons_free_internal (pr, a1);
		ucons_free_internal (pr, a2);
	}

#ifndef NDEBUG1
		fprintf(stdout,"\n \t join returns  \n ");
		ucons_fprint(stdout,man,r,NULL);
	fprintf(stdout,"\n");
#endif

	return r;
}



/* TODO: priority 3 */
ucons_t *
ucons_meet_array (ap_manager_t * man, ucons_t ** tab, size_t size)
{
	ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_JOIN, 0);
	arg_assert (size > 0, return NULL;
	);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return tab[0];
}

/* TODO: priority 3 */
ucons_t *
ucons_join_array (ap_manager_t * man, ucons_t ** tab, size_t size)
{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_JOIN_ARRAY, 0);
	arg_assert (size > 0, return NULL;);
	ucons_t *r = ucons_copy_internal (pr, tab[0]);
	size_t i;
	for (i = 1; i < size; i++)
		r = ucons_join (man, false, tab[i], r);
	return r;

}


/* ============================================================ */
/* Widening, Narrowing */
/* ============================================================ */

/* TODO: priority 3 */
ucons_t *
ucons_widening (ap_manager_t * man, ucons_t * a1, ucons_t * a2)
{
	ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_WIDENING, 0);
	ucons_t * r  = ucons_alloc_internal(pr, a1->datadim, a1->segmentdim );

	arg_assert (a1->datadim == a2->datadim
			&& a1->segmentdim == a2->segmentdim, return NULL;);


#ifndef NDEBUG1
		fprintf(stdout,"\n \t wide a1 \n ");
		ucons_fprint(stdout,man,a1,NULL);
		fprintf(stdout,"\n \t wide  a2 \n ");
		ucons_fprint(stdout,man,a2,NULL);
		fprintf(stdout,"\n");
		fflush(stdout);
#endif


	r->econs = ap_abstract0_widening (pr->man_dcons, a1->econs, a2->econs);
	//	r->n2p  = pattern_key_set_copy (pr, a1->n2p ,a1->segmentdim);

	ap_lincons0_array_t arre = ap_lincons0_array_make (1);

	arre.p[0].constyp = AP_CONS_SUPEQ;
	arre.p[0].linexpr0 =
			ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2*r->segmentdim );
	arre.p[0].scalar = NULL;
	bool ord, ord1, ord2;
	ord = false; ord1 = false; ord2 = false;
	for(size_t yy = 1; yy<r->segmentdim ; yy++){
		for (size_t ll = 1 ; ll< r->segmentdim ; ll++)
			if(ll!=yy){


				ap_linexpr0_set_coeff_scalar_int (arre.p[0].linexpr0,
						r->datadim + yy, 1);
				ap_linexpr0_set_coeff_scalar_int (arre.p[0].linexpr0,
						r->datadim +ll, -1);
				ap_linexpr0_set_cst_scalar_int (arre.p[0].linexpr0, 0);

				ord1 = ap_abstract0_sat_lincons(pr->man_dcons,a1->econs,&arre.p[0]);
				ord2 = ap_abstract0_sat_lincons(pr->man_dcons,a2->econs,&arre.p[0]);

				if(ord1 && ord2){
					if(!ap_abstract0_sat_lincons(pr->man_dcons,r->econs,&arre.p[0])){
#ifndef NDEBUG2
						fprintf(stdout," helping econs\n");
						fflush(stdout);
#endif
					}
					r->econs = ap_abstract0_meet_lincons_array (pr->man_dcons, true,r->econs, &arre);
				}

				ap_linexpr0_set_coeff_scalar_int (arre.p[0].linexpr0,
						r->datadim + yy, 0);
				ap_linexpr0_set_coeff_scalar_int (arre.p[0].linexpr0,
						r->datadim +ll, 0);
			}
	}
	ap_lincons0_array_clear (&arre);

#ifndef NDEBUG1
		fprintf(stdout,"\n \t wide a1 \n ");
		ucons_fprint_econs(stdout,man,a1,NULL);
		fprintf(stdout,"\n");
		fflush(stdout);
#endif

	pattern_t *sa1;
	pattern_t *sa2;

	pattern_t *ru=NULL; //universal constraint to add to r
	size_t u_seg;
	size_t e_seg,nr_y;

	size_t i;
	unsigned keylen;


	ap_dimchange_t dimadd;
	ap_dimchange_init (&dimadd, 2, 0);
	dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
	dimadd.dim[0] = r->datadim + 2*r->segmentdim;
	dimadd.dim[1] = r->datadim + 2*r->segmentdim;

	ap_abstract0_t * aux_wide_a2 = ap_abstract0_add_dimensions(pr->man_dcons,false,a2->econs,&dimadd,false);
	ap_abstract0_t * aux_wide_a1 = ap_abstract0_add_dimensions(pr->man_dcons,false,a1->econs,&dimadd,false);

	ap_dimchange_clear (&dimadd);

	ap_dimchange_init (&dimadd, 4, 0);
	dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
	dimadd.dim[0] = r->datadim + 2*r->segmentdim;
	dimadd.dim[1] = r->datadim + 2*r->segmentdim;
	dimadd.dim[2] = r->datadim + 2*r->segmentdim;
	dimadd.dim[3] = r->datadim + 2*r->segmentdim;

	ap_abstract0_t * y2aux_wide_a2 = ap_abstract0_add_dimensions(pr->man_dcons,false,a2->econs,&dimadd,false);
	ap_abstract0_t * y2aux_wide_a1 = ap_abstract0_add_dimensions(pr->man_dcons,false,a1->econs,&dimadd,false);

	ap_dimchange_clear (&dimadd);

	//#if defined (UCONS_DCONS_OCT_P21) || defined (UCONS_DCONS_POLY_P21)
	for(sa1=a1->udcons; sa1 != NULL; sa1=sa1->hh.next) {

		u_seg = pr->PI[sa1->key.type].u_seg; // number of universal segments
		e_seg = pr->PI[sa1->key.type].e_seg;
		nr_y = pr->PI[sa1->key.type].nr_y;

		keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

		ap_dim_t n = sa1->key.segments[0];
		bool pats1_sat_a2 =!test_singleton(pr->man_dcons, a2->econs, a2->datadim, a2->segmentdim ,n);
				//test_l_leq_2(pr->man_dcons, a2->econs, r->datadim, r->segmentdim ,n);

		HASH_FIND(hh,a2->udcons,&sa1->key,keylen,sa2);

		HASH_FIND(hh,r->udcons,&sa1->key,keylen,ru);

		//checked_malloc(ru,pattern_t,1,sizeof(pattern_t)+u_seg*sizeof(size_t),return NULL;);
		//memset(ru, 0, sizeof(pattern_t)+u_seg*sizeof(size_t));

		ap_abstract0_t *dcons1, *dcons2;
		ap_dimchange_t dimadd;
		size_t i;
	//	ap_abstract0_closure(pr->man_dcons,false,sa2->dcons)
		if(sa2){
			if (ru) {
				ru->dcons=ap_abstract0_widening(pr->man_dcons,sa1->dcons,sa2->dcons);
#ifndef NDEBUG1
		fprintf(stdout,"\n \t wide neintarit  \n ");
		fprintf(stdout,"\n");
		fflush(stdout);
#endif
			}
			else{

				pattern_t *rn;

				u_seg = pr->PI[sa1->key.type].u_seg;
				e_seg = pr->PI[sa1->key.type].e_seg;
				keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

				checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
				memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));

				rn->dcons = ap_abstract0_widening(pr->man_dcons,sa1->dcons,sa2->dcons);

				ap_lincons0_array_t arr = ap_lincons0_array_make (1);

				arr.p[0].constyp = AP_CONS_SUPEQ;
				arr.p[0].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2*r->segmentdim + 2*nr_y );
				arr.p[0].scalar = NULL;
				bool ord, ord1, ord2;
				ord = false; ord1 = false; ord2 = false;

				if(nr_y==2){
					// d(y1) - d(y2) >=0
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 2, 1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 3, -1);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

					ord1 = ap_abstract0_sat_lincons(pr->man_dcons,sa1->dcons,&arr.p[0]);
					ord2 = ap_abstract0_sat_lincons(pr->man_dcons,sa2->dcons,&arr.p[0]);

					if(ord1 && ord2)
						rn->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true,rn->dcons, &arr);

					// d(y2) - d(y1) >=0
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 2, -1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 3, 1);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

					ord1 = ap_abstract0_sat_lincons(pr->man_dcons,sa1->dcons,&arr.p[0]);
					ord2 = ap_abstract0_sat_lincons(pr->man_dcons,sa2->dcons,&arr.p[0]);

					if(ord1 && ord2)
							rn->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true,rn->dcons, &arr);

					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 2, 0);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 3, 0);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

				}
				ap_dim_t dctr, dy1, dy2;
				for(size_t ctr = 1 ; ctr < a1->segmentdim; ctr ++){
					dctr = a1->datadim + ctr;
					if(nr_y == 1)
						dy1 =  r->datadim + 2*r->segmentdim + 1;
					else
						dy1 =  r->datadim + 2*r->segmentdim + 2;

					dy2 = r->datadim + 2*r->segmentdim + 3;

					// d(ctr) - d(y1) >= 0
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy1, -1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dctr, 1);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

					ord1 = ap_abstract0_sat_lincons(pr->man_dcons,sa1->dcons,&arr.p[0]);
					ord2 = ap_abstract0_sat_lincons(pr->man_dcons,sa2->dcons,&arr.p[0]);

					if(ord1 && ord2)
							rn->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, rn->dcons, &arr);

				//	 d(y1) - d(ctr) >= 0
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy1, 1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dctr, -1);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

					ord1 = ap_abstract0_sat_lincons(pr->man_dcons,sa1->dcons,&arr.p[0]);
					ord2 = ap_abstract0_sat_lincons(pr->man_dcons,sa2->dcons,&arr.p[0]);

					if(ord1 && ord2)
						rn->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, rn->dcons, &arr);

					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy1, 0);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dctr, 0);

					if(nr_y == 2){
					//	 d(ctr) - d(y2) >= 0
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy2, -1);
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dctr, 1);
						ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

						ord1 = ap_abstract0_sat_lincons(pr->man_dcons,sa1->dcons,&arr.p[0]);
						ord2 = ap_abstract0_sat_lincons(pr->man_dcons,sa2->dcons,&arr.p[0]);

						if(ord1 && ord2)
							rn->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, rn->dcons, &arr);

						// d(y2) - d(ctr) >= 0
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy2, 1);
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dctr, -1);
						ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

						ord1 = ap_abstract0_sat_lincons(pr->man_dcons,sa1->dcons,&arr.p[0]);
						ord2 = ap_abstract0_sat_lincons(pr->man_dcons,sa2->dcons,&arr.p[0]);

						if(ord1 && ord2)
							rn->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, rn->dcons, &arr);

						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy2, 0);
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dctr, 0);
					}

				}


				ap_lincons0_array_clear (&arr);

//				rn->dcons = ap_abstract0_meet(pr->man_dcons,true,rn->dcons,pt1);
#ifndef NDEBUG1
		fprintf(stdout,"\n \t in wide result after meet with join \n ");
		ap_abstract0_fprint(stdout,pr->man_dcons,rn->dcons,NULL);
		fprintf(stdout,"\n");
		fflush(stdout);
#endif

				rn->key.type=sa1->key.type;

				for(size_t i=0; i<(u_seg+e_seg); i++)
					rn->key.segments[i]=sa1->key.segments[i];


				HASH_ADD(hh,r->udcons,key,keylen,rn);
				add_pattern_n2p(pr,r,&rn->key);
#ifndef NDEBUG1
		fprintf(stdout,"\n \t wide s1 \n ");
		ucons_fprint_dcons(stdout,man,a1,NULL,&sa1->key);
		fprintf(stdout,"\n \t with s2 \n ");
		ucons_fprint_dcons(stdout,man,a2,NULL,&sa1->key);
		fprintf(stdout,"\n \t wide->ru->dcons is \n ");
		ucons_fprint_dcons(stdout,man,r,NULL,&sa1->key);
		fprintf(stdout,"\n");
		fflush(stdout);
#endif

			}
		}
		else {
			//HASH_FIND(hh,r->udcons,&sa1->key,keylen,ru);
			if(!pats1_sat_a2){
				if (ru) {
					ru->dcons=ap_abstract0_copy(pr->man_dcons,sa1->dcons);
				}
				else{

					pattern_t *rn;

					u_seg = pr->PI[sa1->key.type].u_seg;
					e_seg = pr->PI[sa1->key.type].e_seg;

					keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

					checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
					memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));

					rn->dcons=ap_abstract0_copy(pr->man_dcons,sa1->dcons);

					rn->key.type=sa1->key.type;

					for(size_t i=0; i<(u_seg+e_seg); i++)
						rn->key.segments[i]=sa1->key.segments[i];


					HASH_ADD(hh,r->udcons,key,keylen,rn);
					add_pattern_n2p(pr,r,&rn->key);

				}
			}//end no instances for the pattern length of n == 1 strict
			else{

				if(ru){
					if(nr_y == 1)
						ru->dcons = ap_abstract0_widening(pr->man_dcons,aux_wide_a2,sa1->dcons);
					else
						ru->dcons = ap_abstract0_widening(pr->man_dcons,y2aux_wide_a2,sa1->dcons);
				}
				else{
					pattern_t *rn;

					u_seg = pr->PI[sa1->key.type].u_seg;
					e_seg = pr->PI[sa1->key.type].e_seg;
					keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

					checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
					memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));
#ifndef NDEBUG1
		fprintf(stdout,"\n \t wide neintarit  \n ");
		fprintf(stdout,"\n");
		fflush(stdout);
#endif
					if(nr_y == 1)
						rn->dcons=ap_abstract0_widening(pr->man_dcons,aux_wide_a2,sa1->dcons);
					else
						rn->dcons=ap_abstract0_widening(pr->man_dcons,y2aux_wide_a2,sa1->dcons);

					rn->key.type=sa1->key.type;

					for(size_t i=0; i<(u_seg+e_seg); i++)
						rn->key.segments[i]=sa1->key.segments[i];


					HASH_ADD(hh,r->udcons,key,keylen,rn);
					add_pattern_n2p(pr,r,&rn->key);
				}
			}
		}
	}

	sa1 = NULL;

	for(sa2=a2->udcons; sa2 != NULL; sa2=sa2->hh.next) {

		u_seg = pr->PI[sa2->key.type].u_seg; // number of universal segments
		e_seg = pr->PI[sa2->key.type].e_seg;
		nr_y = pr->PI[sa2->key.type].nr_y;

		keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

		ap_dim_t n = sa2->key.segments[0];
		bool pats2_sat_a1 = !test_singleton(pr->man_dcons, a1->econs, a1->datadim, a1->segmentdim ,n);
				//test_l_leq_2(pr->man_dcons, a1->econs, r->datadim, r->segmentdim ,n);

		HASH_FIND(hh,a1->udcons,&sa2->key,keylen,sa1);

		HASH_FIND(hh,r->udcons,&sa2->key,keylen,ru);

		//checked_malloc(ru,pattern_t,1,sizeof(pattern_t)+u_seg*sizeof(size_t),return NULL;);
		//memset(ru, 0, sizeof(pattern_t)+u_seg*sizeof(size_t));

		ap_abstract0_t *dcons1, *dcons2;
		ap_dimchange_t dimadd;
		size_t i;

		if(!sa1){
			if(!pats2_sat_a1){
				if (ru) {
					ru->dcons=ap_abstract0_copy(pr->man_dcons,sa2->dcons);
				}
				else{

					pattern_t *rn;

					u_seg = pr->PI[sa2->key.type].u_seg;
					e_seg = pr->PI[sa2->key.type].e_seg;

					keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

					checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
					memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));

					rn->dcons=ap_abstract0_copy(pr->man_dcons,sa2->dcons);

					rn->key.type=sa2->key.type;

					for(size_t i=0; i<(u_seg+e_seg); i++)
						rn->key.segments[i]=sa2->key.segments[i];


					HASH_ADD(hh,r->udcons,key,keylen,rn);
					add_pattern_n2p(pr,r,&rn->key);

				}
			}// no instances for the pattern
			else{
				if(ru){
					if(nr_y == 1)
						ru->dcons = ap_abstract0_widening(pr->man_dcons,aux_wide_a1,sa2->dcons);
					else
						ru->dcons = ap_abstract0_widening(pr->man_dcons,y2aux_wide_a1,sa2->dcons);
				}
				else{
					pattern_t *rn;

					u_seg = pr->PI[sa2->key.type].u_seg;
					e_seg = pr->PI[sa2->key.type].e_seg;
					keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

					checked_malloc(rn,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL ;);
					memset(rn, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));

					if(nr_y == 1)
						rn->dcons=ap_abstract0_widening(pr->man_dcons,aux_wide_a1,sa2->dcons);
					else
						rn->dcons=ap_abstract0_widening(pr->man_dcons,y2aux_wide_a1,sa2->dcons);

					rn->key.type=sa2->key.type;

					for(size_t i=0; i<(u_seg+e_seg); i++)
						rn->key.segments[i]=sa2->key.segments[i];


					HASH_ADD(hh,r->udcons,key,keylen,rn);
					add_pattern_n2p(pr,r,&rn->key);
				}
			}
		}
	}

	if(aux_wide_a1) ap_abstract0_free(pr->man_dcons,aux_wide_a1);
   if(aux_wide_a2) ap_abstract0_free(pr->man_dcons,aux_wide_a2);


	return r;
}

/* TODO: priority 3 */
ucons_t *
ucons_widening_thresholds (ap_manager_t * man,
		ucons_t * a1, ucons_t * a2,
		ap_scalar_t ** array, size_t nb)
		{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_WIDENING, nb + 1);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return a2;
		}

/* NOT IMPLEMENTED */
ucons_t *
ucons_narrowing (ap_manager_t * man, ucons_t * a1, ucons_t * a2)
{
	ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_WIDENING, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return a2;
}
