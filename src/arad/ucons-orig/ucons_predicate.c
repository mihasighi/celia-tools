/*
 * ucons_predicate.c
 *
 * Predicates on uconss, extraction functions
 *
 * APRON Library / Shape Domain
 *
 * Copyright (C) LIAFA 2009
 */

/*
 * This file is part of the APRON Library, released under LGPL license.
 * Please read the COPYING file packaged in the distribution.
 */

#include "ucons.h"
#include "ucons_internal.h"
#include "ucons_fun.h"
#include "ap_generic.h"

#include <stdio.h>



/* ============================================================ */
/* Tests */
/* ============================================================ */

/* ********************************************************************** */
/* I. ucons_pattern_key_set_t */
/* ********************************************************************** */

/* ********************************************************************** */
/* III. pattern_t */
/* ********************************************************************** */

/* ********************************************************************** */
/* IV. ucons_pattern_t */
/* ********************************************************************** */
bool is_pattern_inst(ucons_internal_t * pr, ucons_t *r, pattern_key_t *p ){

	size_t type = p->type;
	size_t nry = pr->PI[type].nr_y;
	size_t useg = pr->PI[type].u_seg;

	size_t len = 1;
	ap_dim_t n = p->segments[0];

	bool res = true;

	switch(pr->PI[type].kind){
	case  pattern_2_1 :
	case pattern_1 :
	case pattern_1_l1:
	{
		res = !test_singleton(pr->man_dcons,r->econs, r->datadim, r->segmentdim , n);
		break;
	}
	case  pattern_1_lx_1 :
	{
		//pattern strengthens econs in split only if it's size is 2
		res = (!test_singleton(pr->man_dcons,r->econs, r->datadim, r->segmentdim , n)
				&& test_l_leq_v(pr->man_dcons,r->econs, r->datadim, r->segmentdim , n, 2));
				break;
	}
	case  pattern_succ_1_2 :
	{
		//res == false is fine also. not need to strengthen econs with this pattern
		res =  test_g_leq_v(pr->man_dcons,r->econs, r->datadim, r->segmentdim , n, 3);
		break;
	}

	default:{
		bool flag = true;
		//TODO test completely the patterns satisfiability when is defined on multiple segment
		for (size_t i = 0; i < (useg) && flag == true; i++){
			n =  p->segments[i];
			flag = !test_singleton(pr->man_dcons,r->econs, r->datadim, r->segmentdim , n);
		}
		res=flag;
	}
	}

	return res;
}
/* ********************************************************************** */
/* V. ucons_array_t */
/* ********************************************************************** */

/* The bottom value is the empty set. */
bool
ucons_is_bottom (ap_manager_t * man, ucons_t * a)
{
	ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_IS_BOTTOM, 0);
	return (!a || ap_abstract0_is_bottom (pr->man_dcons, a->econs));
}

/* The top value is top ucons in the first position. */
bool
ucons_is_top (ap_manager_t * man, ucons_t * a)
{
	size_t i;
	ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_IS_TOP, 0);

	pattern_t *s;
	for(s=a->udcons;s!=NULL;s=s->hh.next)
		if((s->dcons==NULL) || (ap_abstract0_is_top (pr->man_dcons, s->dcons)==false))
			return false;

	return (a && ap_abstract0_is_top (pr->man_dcons, a->econs));
}

/* TODO: priority 3 */
bool
ucons_is_leq (ap_manager_t * man, ucons_t * a1, ucons_t * a2)
{
	ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_IS_LEQ, 0);

//#ifndef NDEBUG2
//	fprintf(stdout,"\t a1 : \n");
//	ucons_fprint(stdout,pr->man,a1, NULL);
//	fprintf(stdout,"\n");
//	fprintf(stdout,"\t implies a2 : \n");
//	ucons_fprint(stdout,pr->man,a2, NULL);
//	fprintf(stdout,"\n");
//	fflush(stdout);
//#endif

	if (ucons_is_bottom (man, a1)){
		fprintf(stdout,"\n a1 bottom");
		fflush(stdout);
		return true;
	}
	if (ucons_is_bottom (man, a2)){
		fprintf(stdout,"\n a2 bottom");
		fflush(stdout);
		return false;
	}
	if (a1->datadim != a2->datadim || a1->segmentdim != a2->segmentdim)
		return false;
	bool is_leq = ap_abstract0_is_leq (pr->man_dcons, a1->econs, a2->econs);

	if (!is_leq ){
#ifndef NDEBUG2
		fprintf(stdout,"\n exists not leq ");
		fflush(stdout);
#endif
		return false;
	}
	pattern_t * s1, *s2;
	size_t u_seg;
	size_t e_seg;
	unsigned keylen;

	for(s1 = a1->udcons; s1!= NULL; s1=s1->hh.next){
		u_seg = pr->PI[s1->key.type].u_seg;
		e_seg = pr->PI[s1->key.type].e_seg;
		keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);
		HASH_FIND(hh,a2->udcons,&s1->key.type,keylen,s2);
		if(s2){
			if (!ap_abstract0_is_leq (pr->man_dcons, s1->dcons, s2->dcons)){
				#ifndef NDEBUG2
				fprintf(stdout,"\n ucons_leq returns false");
				#endif
				/*
	fprintf(stdout,"\n a1:");
	ucons_fprint(stdout, pr->man, a1, NULL);
	fprintf(stdout,"\n not leq: ");
	fprintf(stdout,"\n a2:");
	ucons_fprint(stdout, pr->man, a2, NULL);
				 */

				#ifndef NDEBUG2
				fprintf(stdout,"\n pattern_type := %zu", s1->key.type );

				ap_abstract0_fprint(stdout,pr->man_dcons,s1->dcons,NULL);
				fprintf(stdout,"\n not leq \n");
				ap_abstract0_fprint(stdout,pr->man_dcons,s2->dcons,NULL);
				fflush(stdout);
				#endif

				return false;
			}
		}
	}

	for(s2 = a2->udcons; s2!= NULL; s2=s2->hh.next){
		u_seg = pr->PI[s2->key.type].u_seg;
		e_seg = pr->PI[s2->key.type].e_seg;
		keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);
		HASH_FIND(hh,a1->udcons,&s2->key.type,keylen,s1);
		if(!s1){
			/* only if the patterns is instantialble in a1 consider the test
			 * the segments in s2 universal quantified, their lengths is greather then 1 in a1
			 * */
			bool sat_pattern = true;
			size_t useg = pr->PI[s2->key.type].u_seg;

			for(size_t i = 0 ; i<u_seg; i++){
				ap_linexpr0_t* linexpr = ap_linexpr0_alloc ( AP_LINEXPR_DENSE,
									a1->datadim + 2 * a1->segmentdim);
							ap_linexpr0_set_cst_scalar_int (linexpr, -1);
				ap_linexpr0_set_coeff_scalar_int (linexpr,
						a1->datadim + a1->segmentdim + s2->key.segments[i] , 1);
				ap_lincons0_t cons = ap_lincons0_make(AP_CONS_EQ, linexpr, NULL);
				if(ap_abstract0_sat_lincons(pr->man_dcons,a1->econs,&cons))
					sat_pattern=false;

				ap_lincons0_clear(&cons);
				linexpr = NULL;
			}

			if (sat_pattern&&!ap_abstract0_is_top (pr->man_dcons, s2->dcons)){
#ifndef NDEBUG2
	fprintf(stdout,"\n sat ucons_leq returns false");
	fflush(stdout);
#endif
				return false;
			}
		}
	}

#ifndef NDEBUG2
	fprintf(stdout,"\n ucons_leq returns true");
	fflush(stdout);
#endif
	return true;

	//	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
	//			"not implemented");
}

/* TODO: priority 3 */
bool
ucons_is_eq (ap_manager_t * man, ucons_t * a1, ucons_t * a2)
{
	ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_IS_EQ, 0);

	if ((ucons_is_bottom (man, a1) && ucons_is_bottom (man, a2)) ||
			(ucons_is_top (man, a1) && ucons_is_top (man, a2)))
		return true;
	if (a1->datadim != a2->datadim || a1->segmentdim != a2->segmentdim)
		return false;
	if(! ap_abstract0_is_eq (pr->man_dcons, a1->econs, a2->econs))
		return false;


	return ucons_is_leq(man,a1,a2) && ucons_is_leq(man,a2,a1);
//
//	pattern_t * s1, *s2;
//	size_t u_seg;
//	size_t e_seg;
//	unsigned keylen;
//
//	for(s1 = a1->udcons; s1!= NULL; s1=s1->hh.next){
//		u_seg = pr->PI[s1->key.type].u_seg;
//		e_seg = pr->PI[s1->key.type].e_seg;
//		keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);
//		HASH_FIND(hh,a2->udcons,&s1->key.type,keylen,s2);
//		if(s2){
//			if (!ap_abstract0_is_eq (pr->man_dcons, s1->dcons, s2->dcons))
//				return false;
//		}
//	}
//
//	return true;
//	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
//			"not implemented");
}

/* TODO: priority 1 */
/*
 * Possibly needed for checking linear constraints representing - aliasing
 * between pointer variables, e.g., x = y, - or constraints between the
 * program variables (lengths and data) in the assert statements.
 */
bool
ucons_sat_lincons (ap_manager_t * man, ucons_t * a, ap_lincons0_t * lincons)
{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_SAT_LINCONS, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return true;
}

/* TODO: priority 1 */
bool
ucons_sat_tcons (ap_manager_t * man, ucons_t * a, ap_tcons0_t * cons)
{
	ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_SAT_TCONS, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return true;
}

/* TODO: priority 1 */
/* Interval constraints are only between non-pointer variables */
bool
ucons_sat_interval (ap_manager_t * man, ucons_t * a,
		ap_dim_t dim, ap_interval_t * i)
{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_SAT_INTERVAL, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return true;
}

/* TODO: priority 1 */
bool
ucons_is_dimension_unconstrained (ap_manager_t * man, ucons_t * a,
		ap_dim_t dim)
{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_IS_DIMENSION_UNCONSTRAINED, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return false;
}

/* ============================================================ */
/* Extraction of properties */
/* ============================================================ */

/* NOT IMPLEMENTED */
ap_interval_t *
ucons_bound_linexpr (ap_manager_t * man, ucons_t * a, ap_linexpr0_t * expr)
{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_BOUND_LINEXPR, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return NULL;
}

/* NOT IMPLEMENTED */
ap_interval_t *
ucons_bound_texpr (ap_manager_t * man, ucons_t * a, ap_texpr0_t * expr)
{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_BOUND_TEXPR, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return NULL;
}

/* NOT IMPLEMENTED */
ap_interval_t *
ucons_bound_dimension (ap_manager_t * man, ucons_t * a, ap_dim_t dim)
{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_BOUND_DIMENSION, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return NULL;
}

/* NOT IMPLEMENTED */
ap_lincons0_array_t
ucons_to_lincons_array (ap_manager_t * man, ucons_t * a)
{
	ap_lincons0_array_t ar;
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_TO_LINCONS_ARRAY, 0);
	ar = ap_lincons0_array_make (1);
	ar.p[0] = ap_lincons0_make_unsat ();
	return ar;
}

/* NOT IMPLEMENTED */
ap_tcons0_array_t
ucons_to_tcons_array (ap_manager_t * man, ucons_t * a)
{
	return ap_generic_to_tcons_array (man, a);
}

/* NOT IMPLEMENTED */
ap_interval_t **
ucons_to_box (ap_manager_t * man, ucons_t * a)
{
	ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_TO_BOX, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return NULL;
}

/* NOT IMPLEMENTED */
ap_generator0_array_t
ucons_to_generator_array (ap_manager_t * man, ucons_t * a)
{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_TO_GENERATOR_ARRAY, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return ap_generator0_array_make (0);
}
