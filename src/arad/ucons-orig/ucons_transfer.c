/*
 * ucons_transfer.c
 *
 * Assignment, substitution, guard transfer functions.
 *
 * APRON Library / Shapes Domain
 *
 * Copyright (C) LIAFA 2009
 *
 */

/*
 * This file is part of the APRON Library, released under LGPL license.
 * Please read the COPYING file packaged in the distribution.
 */

#include <assert.h>
#include "ucons.h"
#include "ucons_internal.h"
#include "apron2shape.h"
#include "ap_generic.h"
#include "ap_linexpr0.h"
#include "ap_pcons0.h"


/* ============================================================ */
/* Meet constraints and Join generators */
/* ============================================================ */

bool * mset_nodes = NULL; // nodes to eliminate

size_t ** eq_mset_nodes = NULL;/*  array of mset inclusion between mset_nodes and
								 the other nodes
								 eq_mset_nodes[i][j][0] == 1 <=> mset(i) included in mset(j)
 */

ucons_t *
ucons_meet_lincons (ucons_internal_t * pr, bool destructive,
		ucons_t * a, ap_lincons0_t * lcons)
		{


	//	arg_assert (lcons!=NULL, return NULL;);
	//	arg_assert (a!=NULL, return NULL;);

	ucons_t *r = ucons_copy_internal (pr, a);

#ifndef NDEBUG1
	printf("@@@@ ucons_meet_lincons: with lincons=(");
	if(lcons) ap_lincons0_fprint(stdout, lcons, NULL);
	printf(") \n and offset (lcons->scalar): ");
	ap_scalar_fprint(stdout,lcons->scalar);
	printf("\n on a=(");
	ucons_fprint(stdout, pr->man, a,NULL);
	printf( ") and r=(");
	ucons_fprint(stdout, pr->man, r,NULL);
	printf( ")\n");
	fflush(stdout);
#endif

	ap_linexpr0_t *lexpr;
	bool eqstruct = false;
	size_t i, dim;
	offset_t  kind;
	ap_coeff_t *coeff = NULL;
	ap_constyp_t op = lcons->constyp;

	kind = OFFSET_OTHER;		/* unknown */
	if (lcons->scalar)
	{
		if (!ap_scalar_cmp_int (lcons->scalar, OFFSET_DATA))
			kind = OFFSET_DATA;	/* data */
		else if (!ap_scalar_cmp_int (lcons->scalar, OFFSET_LEN))
			kind = OFFSET_LEN;	/* len */
		else if(!ap_scalar_cmp_int (lcons->scalar, OFFSET_UCONS))
			kind = OFFSET_UCONS;
		else if(!ap_scalar_cmp_int (lcons->scalar, OFFSET_MSET))
			kind = OFFSET_MSET;	/* mset */
		else
			if(op == AP_CONS_EQ || op == AP_CONS_EQMOD)
			{
				eqstruct = true;
				kind = OFFSET_OTHER;
			}
			else
				assert (0);
	}

#ifndef NDEBUG1
	fprintf (stdout, "@@@@ ucons_meet_lincons: kind=%d\n", kind);
	fflush (stdout);
#endif

	lexpr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * a->segmentdim);
	ap_dim_t * useg = NULL; // for EQ_MOD: the segments with equality constraint
	size_t size_useg = 0;

	ap_dim_t n1, n2, n3; // Encoding to know when is
	// an eq constraint and when is just a node to eliminate
	n1 = 0; n2 = 0; n3 = 0;

	ap_linexpr0_ForeachLinterm (lcons->linexpr0, i, dim, coeff)
	{
		if (coeff && !ap_coeff_zero (coeff))
		{
			if(kind == OFFSET_MSET){
				/* transform the coefficient into scalar
				 * to types of constraints handled :
				 * 		n1 != 0 n2 != 0 n3==0 ==> mset(n1) == mset(n2)
				 * 	    n1 != 0 n2 != 0 n3!=0 ==> mset(n1)+mset(n2)= tl(n3)
				 * */
				ap_coeff_t * c1 = ap_coeff_alloc(AP_COEFF_SCALAR);
				ap_coeff_t * c2 = ap_coeff_alloc(AP_COEFF_SCALAR);
				ap_coeff_t * c3 = ap_coeff_alloc(AP_COEFF_SCALAR);
				ap_coeff_set_scalar_int(c1,1);
				ap_coeff_set_scalar_int(c2,-1);
				ap_coeff_set_scalar_int(c3,-2);
				if(ap_coeff_equal(coeff,c1)){
					if (n1==0) n1 = dim - a->datadim;
					else n2 = dim - a->datadim;
				}
				if(ap_coeff_equal(coeff,c2)){
					n2 = dim - a->datadim;
				}
				if(ap_coeff_equal(coeff,c3)){
					n3 = dim - a->datadim;
				}
			}
			else
				if (eqstruct){
					/*
					 * equal segments from call function
					 */
					checked_realloc(useg, ap_dim_t,
							(size_useg + 1),sizeof(ap_dim_t), return NULL;);
					useg[size_useg++] = dim-a->datadim;//+1 PB
				}
				else{
					if (dim < a->datadim)
						ap_coeff_set (&lexpr->p.coeff[dim], coeff);
					else			// dim >= a->datadim
					{
						if (kind == OFFSET_DATA || kind == OFFSET_UCONS)
							ap_coeff_set (&lexpr->p.coeff[dim], coeff);
						else
							if (kind == OFFSET_LEN)
								ap_coeff_set (&lexpr->p.coeff[dim +  a->segmentdim], coeff);
							else
								assert (0);
					}
				}
		}
	}


#ifndef NDEBUG1
	fprintf(stdout,"@@@@ ucons_meet_licons: lexpr=(");
	if(lexpr) ap_linexpr0_fprint(stdout, lexpr, NULL);
	fprintf(stdout,")\n");
	fflush(stdout);
#endif

	if(kind == OFFSET_UCONS){
                // Warning: ap_coeff_equal_int is not correct, use scalar
		ap_coeff_t * coeff = ap_linexpr0_cstref(lcons->linexpr0);
                ap_scalar_t* code = NULL;
                if (coeff && coeff->discr == AP_COEFF_SCALAR)
                    code = coeff->val.scalar;
		r = build_constraint(pr, r, lexpr, code);

		ap_linexpr0_free (lexpr);
	}
	else{

		if(kind == OFFSET_MSET){
#ifndef NDEBUG1
			fprintf(stdout,"====ucons_meet_mset_constraint \n ");
			if(lcons) ap_lincons0_fprint(stdout, lcons, NULL);
			fprintf(stdout,"\n n1 = %zu n2 = %zu n3=%zu \n ",n1,n2,n3);
			fflush(stdout);
#endif
			if(n2==0){
				if(mset_nodes==NULL){
					checked_realloc(mset_nodes, bool,
							(a->segmentdim),sizeof(bool), return NULL;);
					for(size_t i = 0;i<a->segmentdim; i++)
						mset_nodes[i]=false;
				}
				mset_nodes[n1] = true;
				n1 = 0;
				n3 = 0;
			}
			else {
				if(eq_mset_nodes == NULL){
					checked_realloc(eq_mset_nodes, size_t*,
							(a->segmentdim),sizeof(size_t*), return NULL;);
					for(size_t i = 0;i<a->segmentdim; i++){
						checked_malloc(eq_mset_nodes[i], size_t,
								(a->segmentdim),sizeof(size_t), return NULL;);
						for(size_t j=0 ; j<a->segmentdim; j++)
							eq_mset_nodes[i][j] = 0;
						eq_mset_nodes[i][i] = 1;
					}
				}
				if(n3==0){
					/*  ms(n1) == ms(n2) */
					eq_mset_nodes[n1][n2] = 1;
					eq_mset_nodes[n2][n1] = 1;

					n1 = 0;
					n2 = 0;
				}
				else{
					/*  ms(n1) \cup ms(n2) == ms(n3) */
					eq_mset_nodes[n1][n3] = 2;
					eq_mset_nodes[n2][n3] = 2;
					n1 = 0;
					n2 = 0;
					n3 = 0;
				}
			}
		}
		else
			if(eqstruct){
				/* EQ_MOD */
				ucons_meet_eq_structural_constraint(pr, a, r, useg, size_useg);

			}// end eqconstraint

			else{
				/* meet with normal constraint */
				ap_coeff_t *coeff = ap_coeff_alloc(AP_COEFF_SCALAR);
				ap_linexpr0_get_cst (coeff, lcons->linexpr0);
				ap_linexpr0_set_cst (lexpr, coeff);

				if(op == AP_CONS_SUP ){

					ap_coeff_t *oldcoeff = ap_coeff_alloc(AP_COEFF_SCALAR);
					double oldd;
					ap_linexpr0_get_cst(oldcoeff, lexpr);
					ap_double_set_scalar(&oldd, oldcoeff->val.scalar, GMP_RNDN);
					//       fprintf(stdout,"!!***!! scalar old %f \n", oldd);
					oldd -= 1;
					//        fprintf(stdout,"!!***!! scalar new %f \n", oldd);
					ap_linexpr0_set_cst_scalar_double(lexpr, oldd);
					ap_coeff_free(oldcoeff);
					lcons->constyp = AP_CONS_SUPEQ;

				}

				ap_lincons0_array_t arr = ap_lincons0_array_make (1);
				arr.p[0] = ap_lincons0_copy (lcons);

#ifndef NDEBUG1
				fprintf (stdout, "\n intersection with classic lincons :=");
				ap_lincons0_array_fprint (stdout, &arr, NULL);
				printf( "\n");
				ucons_fprint(stdout, pr->man, r,NULL);
				printf( "\n");
				fprintf (stdout, "\n");
#endif

				r->econs =
						ap_abstract0_meet_lincons_array (pr->man_dcons, false, a->econs, &arr);


				size_t u_seg, e_seg, nr_y;
				pattern_t * s, *ra,*rt;
				unsigned keylen;

				for(s=a->udcons;s!=NULL;s=s->hh.next){


					nr_y = pr->PI[s->key.type].nr_y;
					ap_linexpr0_realloc (lexpr, a->datadim + 2 * a->segmentdim + 2 * nr_y);

					ap_linexpr0_t * lexprr = ap_linexpr0_copy(lexpr);

					ap_lincons0_array_t arrr = ap_lincons0_array_make (1);
					arrr.p[0] = ap_lincons0_make (lcons->constyp, lexprr, NULL);

					u_seg=pr->PI[s->key.type].u_seg;
					e_seg=pr->PI[s->key.type].e_seg;

					keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);
					HASH_FIND(hh,r->udcons,&s->key,keylen,rt);


					if(rt){
						rt->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, false, s->dcons, &arrr);
					}

					else {


						checked_malloc(ra,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL;);
						memset(ra, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));
						//		ra->dcons =(s->dcons) ?
						//				ap_abstract0_copy (pr->man_dcons, s->dcons) : NULL;
						ra->dcons = NULL;
						ra->key.type=s->key.type;

						for(i=0;i<(u_seg+e_seg);i++)
							ra->key.segments[i]=s->key.segments[i];
						ra->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, false, s->dcons, &arrr);
						HASH_ADD(hh,r->udcons,key,keylen,ra);
						add_pattern_n2p(pr,r,&ra->key);
					}

					//ap_linexpr0_free (lexpr);
					ap_lincons0_array_clear (&arrr);

				}
#ifndef NDEBUG2
				fprintf(stdout, "====ucons_meet_lincons");
				fprintf(stdout, " builds arr=(");
				ap_lincons0_array_fprint(stdout, &arr, NULL);
				fprintf(stdout, ") returns: ");
				ucons_fprint(stdout, pr->man, r,NULL);
				fprintf(stdout, "\n");
				fflush(stdout);
#endif
				ap_lincons0_array_clear (&arr);
			}//end no EQ_MOD
	}

	if (destructive)
		ucons_free_internal (pr, a);
	return r;

		}

void ucons_meet_eq_structural_constraint(ucons_internal_t *pr, ucons_t *a, ucons_t *r, ap_dim_t *useg, size_t size_useg )
{
	if(useg){
#ifndef NDEBUG1
		printf("intersection with eq lincons between\n (\n ");
		for(size_t kk=0 ;kk<size_useg; kk++)
			printf("useg[%zu]=%zu\n",kk,useg[kk]);
		fprintf(stdout, " )\n");
		fflush(stdout);
#endif

		ap_lincons0_array_t arr_ex = ap_lincons0_array_make ((size_useg *(size_useg-1)));

		for(size_t i = 0; i<size_useg; i++)
			for(size_t j = i+1; j<size_useg; j++){

				size_t li = r->datadim+ useg[i] + r->segmentdim;
				size_t lj = r->datadim+ useg[j] + r->segmentdim;
				size_t di = r->datadim+ useg[i] ;
				size_t dj = r->datadim+ useg[j] ;

				/* data + len transfer */
				// li-lj == 0
				arr_ex.p[0].constyp = AP_CONS_EQ;
				arr_ex.p[0].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
				arr_ex.p[0].scalar = NULL;
				ap_linexpr0_set_coeff_scalar_int (arr_ex.p[0].linexpr0, li, 1);
				ap_linexpr0_set_cst_scalar_int (arr_ex.p[0].linexpr0, 0);
				ap_linexpr0_set_coeff_scalar_int (arr_ex.p[0].linexpr0, lj, -1);

				// di-dj == 0
				arr_ex.p[1].constyp = AP_CONS_EQ;
				arr_ex.p[1].linexpr0 =
						ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
				arr_ex.p[1].scalar = NULL;
				ap_linexpr0_set_coeff_scalar_int (arr_ex.p[1].linexpr0, di, 1);
				ap_linexpr0_set_cst_scalar_int (arr_ex.p[1].linexpr0, 0);
				ap_linexpr0_set_coeff_scalar_int (arr_ex.p[1].linexpr0, dj, -1);


				r->econs = ap_abstract0_meet_lincons_array (pr->man_dcons, false, a->econs, &arr_ex);




				size_t u_seg, e_seg, nr_y;
				pattern_t * s, *ra,*rt;
				ra = NULL;
				unsigned keylen;

				for(s=a->udcons;s!=NULL;s=s->hh.next){


					nr_y = pr->PI[s->key.type].nr_y;
					ap_linexpr0_realloc (arr_ex.p[0].linexpr0, a->datadim + 2 * a->segmentdim + 2 * nr_y);
					ap_linexpr0_realloc (arr_ex.p[1].linexpr0, a->datadim + 2 * a->segmentdim + 2 * nr_y);

					ap_linexpr0_t * lexprr0 = ap_linexpr0_copy(arr_ex.p[0].linexpr0);
					ap_linexpr0_t * lexprr1 = ap_linexpr0_copy(arr_ex.p[1].linexpr0);

					ap_lincons0_array_t arrr = ap_lincons0_array_make (2);
					arrr.p[0] = ap_lincons0_make (arr_ex.p[0].constyp, lexprr0, NULL);
					arrr.p[1] = ap_lincons0_make (arr_ex.p[1].constyp, lexprr1, NULL);

					u_seg=pr->PI[s->key.type].u_seg;
					e_seg=pr->PI[s->key.type].e_seg;

					keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);
					HASH_FIND(hh,r->udcons,&s->key,keylen,rt);


					if(rt){
						rt->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, false, s->dcons, &arrr);
					}

					else {
						//checked_malloc(ra,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL;);
						ra = (pattern_t*) malloc (sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));
						memset(ra, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));
						//		ra->dcons =(s->dcons) ?
						//				ap_abstract0_copy (pr->man_dcons, s->dcons) : NULL;
						ra->dcons = NULL;
						ra->key.type=s->key.type;

						for(i=0;i<(u_seg+e_seg);i++)
							ra->key.segments[i]=s->key.segments[i];
						ra->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, false, s->dcons, &arrr);
						HASH_ADD(hh,r->udcons,key,keylen,ra);
						add_pattern_n2p(pr,r,&ra->key);
					}
					ap_lincons0_array_clear (&arrr);
				}

#ifndef NDEBUG1
				printf("====ucons_meet_lincons: with lincons=(");
				printf( ") builds arr=(");
				ap_lincons0_array_fprint(stdout, &arr_ex, NULL);
				printf(") returns: ");
				ucons_fprint(stdout, pr->man, r,NULL);
				fprintf(stdout, "\n");
				fflush(stdout);
#endif
				ap_lincons0_array_clear (&arr_ex);


				/* pattern transfer */

				pattern_key_set_t psi = r->n2p[useg[i]];
				pattern_key_set_t psj = r->n2p[useg[j]];

				size_t si = r->n2p[useg[i]].size;
				size_t sj = r->n2p[useg[j]].size;

				if((si>0) || (sj>0)){
					/* transfer universals */
					for(size_t q=0; q<si; q++){
						pattern_key_t *pq = psi.p[q];
						u_seg = pr->PI[pq->type].u_seg;
						e_seg = pr->PI[pq->type].e_seg;
						keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);
						pattern_t *auxi, *auxj;
						HASH_FIND(hh,r->udcons,pq,keylen,auxi);

						pattern_key_t *lookj = NULL;
						//	checked_malloc(lookj, pattern_key_t, sizeof(pattern_key_t) + (u_seg+e_seg)*sizeof(size_t), 1, return NULL;);
						lookj = (pattern_key_t*) malloc (sizeof(pattern_key_t)+(u_seg+e_seg)*sizeof(size_t));
						memset(lookj, 0, sizeof(pattern_key_t) + (u_seg+e_seg)*sizeof(size_t));
						lookj->type = pq->type;
						/* copy + modif + sort */
						bool sort_eseg = false;
						for(size_t kk = 0 ; kk<u_seg+e_seg; kk++){
							if(pq->segments[kk] == useg[i]){
								if(kk>=u_seg) sort_eseg = true;
								lookj->segments[kk] = useg[j];
							}
						}
						if(sort_eseg){
							/*sorting exist segment */
							for (size_t ii = u_seg; ii < u_seg+e_seg; ii++)
							{
								size_t  jj = u_seg;
								while (jj != ii && lookj->segments[jj] <= lookj->segments[ii])
									jj++;
								if (jj < ii)
								{
									size_t d = lookj->segments[ii];
									size_t tmp;
									size_t kk;
									for (size_t  tmp = ii ; tmp > jj; tmp--){
										kk = tmp -1;
										lookj->segments[kk + 1] = lookj->segments[kk];
									}
									lookj->segments[jj] = d;
								}
							}
						}
						else
							if(u_seg==2){
								if(lookj->segments[0]>lookj->segments[1])
								{
									size_t jj = lookj->segments[0];
									lookj->segments[0] = lookj->segments[1];
									lookj->segments[0] = jj;
								}
							}
						//unsigned keylenj = (u_seg+e_seg)*sizeof(size_t) + sizeof(pattern_key_t);
						HASH_FIND(hh, r->udcons, lookj, keylen, auxj);

						if(auxi){
							if(auxj){
								ap_abstract0_t* dcons = ap_abstract0_meet(pr->man_dcons,true,
										auxi->dcons,auxj->dcons);
								auxi->dcons = ap_abstract0_copy(pr->man_dcons,dcons);
								auxj->dcons = ap_abstract0_copy(pr->man_dcons,dcons);
								ap_abstract0_free(pr->man_dcons,dcons);
							}
							else{
								//								checked_malloc(auxj,pattern_t,1,
								//										sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL;);
								auxj = (pattern_t*) malloc (sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));
								memset(auxj, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));
								auxj->key.type = lookj->type;
								for (size_t i=0 ; i<(u_seg+e_seg); i++)
									auxj->key.segments[i] = lookj->segments[i];
								auxj->dcons = ap_abstract0_copy(pr->man_dcons, auxi->dcons);

								HASH_ADD(hh,r->udcons,key,keylen,auxj);
								add_pattern_n2p(pr,r,lookj);
							}
						}

					}

					for(size_t q=0; q<sj; q++){
						pattern_key_t *pq = psj.p[q];
						size_t u_seg = pr->PI[pq->type].u_seg;
						size_t e_seg = pr->PI[pq->type].e_seg;
						keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);
						pattern_t *auxi, *auxj;
						HASH_FIND(hh,r->udcons,pq,keylen,auxj);

						pattern_key_t *looki = NULL;
						//checked_malloc(looki, pattern_key_t, sizeof(pattern_key_t) + (u_seg+e_seg)*sizeof(size_t), 1, return NULL;);
						looki = (pattern_key_t*) malloc (sizeof(pattern_key_t)+(u_seg+e_seg)*sizeof(size_t));
						memset(looki, 0, sizeof(pattern_key_t) + (u_seg+e_seg)*sizeof(size_t));
						looki->type = pq->type;
						/* copy + modif + sort */
						bool sort_eseg = false;
						for(size_t kk = 0 ; kk<u_seg+e_seg; kk++){
							if(pq->segments[kk] == useg[j]){
								if(kk>=u_seg) sort_eseg = true;
								looki->segments[kk] = useg[i];
							}
						}
						if(sort_eseg){
							/*sorting exist segment */
							for (size_t ii = u_seg; ii < u_seg+e_seg; ii++)
							{
								size_t  jj = u_seg;
								while (jj != ii && looki->segments[jj] <= looki->segments[ii])
									jj++;
								if (jj < ii)
								{
									size_t d = looki->segments[ii];
									size_t tmp;
									size_t kk;
									for (size_t  tmp = ii ; tmp > jj; tmp--){
										kk = tmp -1;
										looki->segments[kk + 1] = looki->segments[kk];
									}
									looki->segments[jj] = d;
								}
							}
						}
						else
							if(u_seg==2){
								if(looki->segments[0]>looki->segments[1])
								{
									size_t jj = looki->segments[0];
									looki->segments[0] = looki->segments[1];
									looki->segments[0] = jj;
								}
							}
						//unsigned keylenj = (u_seg+e_seg)*sizeof(size_t) + sizeof(pattern_key_t);
						HASH_FIND(hh, r->udcons, looki, keylen, auxi);

						if(auxj){
							if(!auxi){
								//								checked_malloc(auxi,pattern_t,1,
								//										sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL;);
								auxi = (pattern_t*) malloc (sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));
								memset(auxi, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));
								auxi->key.type = looki->type;
								for (size_t i=0 ; i<(u_seg+e_seg); i++)
									auxi->key.segments[i] = looki->segments[i];
								auxi->dcons = ap_abstract0_copy(pr->man_dcons, auxj->dcons);

								HASH_ADD(hh,r->udcons,key,keylen,auxi);
								add_pattern_n2p(pr,r,looki);
							}
						}

					}
				}

				/* equality pattern constraint */
				//#if defined(UCONS_DCONS_OCT_P21) || defined(UCONS_DCONS_POLY_P21)
				if(pr->active_patterns[2]){
#ifndef NDEBUG2
					fprintf (stdout, "\n intersection with eq constraint");
					fprintf (stdout, "\n");
					fflush(stdout);
#endif

					//	pattern_t *ra;
					if (ra) { free(ra); ra = NULL;}
					//checked_malloc(ra,pattern_t,1,sizeof(pattern_t)+(2)*sizeof(size_t),return NULL;);
					ra = (pattern_t*) malloc (sizeof(pattern_t)+(2)*sizeof(size_t));
					memset(ra, 0, sizeof(pattern_t)+(2*sizeof(size_t)));

					ra->dcons = NULL;
					ra->key.type=1;// y1=y2

					if(useg[i]<useg[j]){
						ra->key.segments[0]=useg[i];
						ra->key.segments[1]=useg[j];
					}
					else{
						ra->key.segments[0]=useg[j];
						ra->key.segments[1]=useg[i];
					}
					/* add dims to r->econs */
					ap_dimchange_t* dimchange = ap_dimchange_alloc(4,0);
					dimchange->dim[0] = r->datadim + 2*r->segmentdim;
					dimchange->dim[1] = r->datadim + 2*r->segmentdim;
					dimchange->dim[2] = r->datadim + 2*r->segmentdim;
					dimchange->dim[3] = r->datadim + 2*r->segmentdim;
					ap_dim_t ly1 = r->datadim + 2*r->segmentdim;
					ap_dim_t ly2 = r->datadim + 2*r->segmentdim + 1;
					ap_dim_t dy1 = r->datadim + 2*r->segmentdim + 2;
					ap_dim_t dy2 = r->datadim + 2*r->segmentdim + 3;

					ra->dcons = ap_abstract0_add_dimensions(pr->man_dcons,false, r->econs, dimchange, false);
					ap_dimchange_free(dimchange);

					ap_lincons0_array_t arr = ap_lincons0_array_make (2);

					// ly1-ly2 == 0
					arr.p[0].constyp = AP_CONS_EQ;
					arr.p[0].linexpr0 =
							ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4 );
					arr.p[0].scalar = NULL;
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ly1, 1);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ly2, -1);
					//dy1-dy2=0
					arr.p[1].constyp = AP_CONS_EQ;
					arr.p[1].linexpr0 =
							ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4 );
					arr.p[1].scalar = NULL;
					ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy2, 1);
					ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);
					ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy1, -1);

					ra->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, ra->dcons, &arr);
					ap_lincons0_array_clear (&arr);

					keylen=  2*sizeof(size_t) + sizeof(pattern_key_t);
					HASH_ADD(hh,r->udcons,key,keylen,ra);
					add_pattern_n2p(pr,r,&ra->key);
					//#endif
				}
			}
	}
}



/* TODO: priority 1 */
/* TODO: detect when is needed */
ucons_t *
ucons_meet_lincons_array (ap_manager_t * man,
		bool destructive, ucons_t * a,
		ap_lincons0_array_t * array)
		{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_MEET_LINCONS_ARRAY, 0);

#ifndef NDEBUG1
	fprintf(stdout,"\n ENTER ucons_meet_lincons_array \n:=");
	fflush(stdout);
#endif

	if (!a)
		return NULL;
	if (!array || array->size == 0 || !array->p)
		return (destructive) ? a : ucons_copy(man,a);
	size_t i;

#ifndef NDEBUG1
	fprintf(stdout,"\n @@ucons_meet_lincons_array \n:=");
	ap_lincons0_array_fprint (stdout, array, NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif


#ifndef NDEBUG1
	fprintf(stdout,"\n @@ucons_meet_lincons_array [0] \n:=");
	ap_lincons0_fprint (stdout, &array->p[0], NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

	ucons_t *r = ucons_meet_lincons (pr, false, a, &array->p[0]);

	for (i = 1; i < array->size; i++){
#ifndef NDEBUG1
		fprintf(stdout,"\n @@ucons_meet_lincons_array [%zu] \n:=",i);
		ap_lincons0_fprint (stdout, &array->p[i], NULL);
		fprintf(stdout,"\n");
		fflush(stdout);
#endif
		r = ucons_meet_lincons (pr, true, r, &array->p[i]);
	}

	if ((eq_mset_nodes != NULL )){
		/* if the array coresponded to an intersection with a multiset constraint */

		r = ucons_strengthen(pr,a, eq_mset_nodes);
		free(mset_nodes);
		free(eq_mset_nodes);
		mset_nodes = NULL;
		eq_mset_nodes = NULL;
	}
	if(mset_nodes != NULL){
		free(mset_nodes);
		mset_nodes = NULL;
	}

	if (destructive)
		ucons_free_internal (pr, a);
	return r;

		}

/* TODO: priority 1 */
/* TODO: detect when is needed */
ucons_t *
ucons_meet_tcons_array (ap_manager_t * man,
		bool destructive, ucons_t * a,
		ap_tcons0_array_t * array)
		{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_MEET_LINCONS_ARRAY, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return a;
		}

/* NOT IMPLEMENTED */
ucons_t *
ucons_add_ray_array (ap_manager_t * man,
		bool destructive, ucons_t * a,
		ap_generator0_array_t * array)
		{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_ADD_RAY_ARRAY, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return a;
		}


/* ============================================================ */
/* Assignement and Substitutions */
/* ============================================================ */

/* TODO: priority 3 */
/* Only done if real dimension (ptr vars) are used in constraints */

ucons_t *
ucons_assign_linexpr (ap_manager_t * man,
		bool destructive, ucons_t * a,
		ap_dim_t d, ap_linexpr0_t * expr, ucons_t * dest)
		{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_ASSIGN_LINEXPR_ARRAY, 0);
	if (!a)
		return NULL;
	assert (expr && ap_linexpr0_size (expr) == (a->datadim + a->segmentdim) &&
			d < (a->datadim + a->segmentdim) && !dest);
	ucons_t *r = ucons_alloc_internal (pr, a->datadim, a->segmentdim);

	ap_linexpr0_realloc (expr, a->datadim + 2 * a->segmentdim);

#ifndef NDEBUG1
	fprintf (stdout, "to dcons %d:=", d);
	ap_linexpr0_fprint (stdout, expr, NULL);
	fprintf (stdout, "\n");

	fprintf (stdout, "to ucons a");
	ucons_fprint (stdout,pr->man,a, NULL);
	fprintf (stdout, "\n");
#endif


	//pattern_t *s,*aux;
	ap_abstract0_t *aux_dcons;

	size_t i,u_seg,e_seg;
	unsigned keylen;


	pattern_t * s, *ra,*rt;

	for(s=a->udcons;s!=NULL;s=s->hh.next){
		/*
		 * TODO add dimensions to r->econs to cope with s->dcons; meet between r->econs and s->dcons ;
		 * */
		u_seg=pr->PI[s->key.type].u_seg;
		e_seg=pr->PI[s->key.type].e_seg;

		keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);

		HASH_FIND(hh,r->udcons,&s->key,keylen,rt);
		if(rt){
			rt->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, false, s->dcons, d, expr, NULL);
		}

		else {
			checked_malloc(ra,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL;);
			memset(ra, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));


			ra->dcons = NULL;
			ra->key.type=s->key.type;

			for(i=0;i<(u_seg+e_seg);i++)
				ra->key.segments[i]=s->key.segments[i];
			ra->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, false, s->dcons, d, expr, NULL);
			HASH_ADD(hh,r->udcons,key,keylen,ra);
			add_pattern_n2p(pr,r,&ra->key);
		}

	}

	r->econs = ap_abstract0_assign_linexpr (pr->man_dcons, false, a->econs, d, expr, NULL);

#ifndef NDEBUG1
	fprintf (stdout, "returns ");
	ucons_fprint (stdout, pr->man, r, NULL);
	fprintf (stdout, "\n");
#endif
	if (destructive)
		ucons_free_internal (pr, a);
	return r;
		}

/* TODO: priority 0 */
/* used for pre-image computation */
ucons_t *
ucons_substitute_linexpr (ap_manager_t * man,
		bool destructive, ucons_t * a,
		ap_dim_t d, ap_linexpr0_t * expr, ucons_t * dest)
		{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_SUBSTITUTE_LINEXPR_ARRAY, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return a;
		}

/* TODO: priority 3 */
/* Only done if real dimension (ptr vars) are used in constraints */
ucons_t *
ucons_assign_linexpr_array (ap_manager_t * man,
		bool destructive, ucons_t * a,
		ap_dim_t * tdim,
		ap_linexpr0_t ** texpr,
		size_t size, ucons_t * dest)
		{
	if (size == 1)
		return ucons_assign_linexpr (man, destructive, a, tdim[0], texpr[0],
				dest);

	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_ASSIGN_LINEXPR_ARRAY, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return a;
		}

/* TODO: priority 0 */
/* used for pre-image computation */
ucons_t *
ucons_substitute_linexpr_array (ap_manager_t * man,
		bool destructive, ucons_t * a,
		ap_dim_t * tdim,
		ap_linexpr0_t ** texpr,
		size_t size, ucons_t * dest)
		{
	if (size == 1)
		return ucons_substitute_linexpr (man, destructive, a, tdim[0], texpr[0],
				dest);

	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_SUBSTITUTE_TEXPR_ARRAY, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return a;
		}

/* TODO: priority 3 */
/* Only for constraints between realdim */
ucons_t *
ucons_assign_texpr_array (ap_manager_t * man,
		bool destructive, ucons_t * a,
		ap_dim_t * tdim,
		ap_texpr0_t ** texpr, size_t size, ucons_t * dest)
		{
	ucons_internal_t *pr =
			ucons_init_from_manager (man, AP_FUNID_ASSIGN_TEXPR_ARRAY, 0);
	ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			"not implemented");
	return a;
		}

/* TODO: priority 0 */
/* used only for pre-image computation */
ucons_t *
ucons_substitute_texpr_array (ap_manager_t * man,
		bool destructive, ucons_t * a,
		ap_dim_t * tdim,
		ap_texpr0_t ** texpr,
		size_t size, ucons_t * dest)
		{
	return ap_generic_substitute_texpr_array (man, destructive, a, tdim, texpr,
			size, dest);
		}

ucons_t*
ucons_strengthen(ucons_internal_t* pr, ucons_t* a, size_t **eqdim){


	ucons_t *r  = ucons_copy_internal(pr,a);
#ifndef NDEBUG1
	fprintf(stdout,"\n ucons_strenghten : \n");
	ucons_fprint(stdout, pr->man, r, NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

	close_eq(eqdim, a->segmentdim);

#ifndef NDEBUG1
	fprintf(stdout,"\n eqdim closed");
	for(size_t ll = 0 ; ll < a->segmentdim; ll++)
		for(size_t pp = 0 ; pp < a->segmentdim; pp++)
			if((pp!=ll) && eqdim[pp][ll]==1)
				fprintf(stdout," ms(%zu) <= ms(%zu)  ", pp,ll);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif

	//	while(!end){
	//		end = true;
	size_t i,j;
	pattern_key_t * look = NULL;
	// look for a pattern forall y
	checked_malloc(look,pattern_key_t, sizeof(pattern_key_t) + 1*sizeof(size_t), 1, return NULL;);
	memset(look, 0, sizeof(pattern_key_t) + 1*sizeof(size_t));
	look->type = 0;
	unsigned keylen = 1*sizeof(size_t) + sizeof(pattern_key_t);

	pattern_key_t * look2 = NULL;
	checked_malloc(look2,pattern_key_t, sizeof(pattern_key_t) + 1*sizeof(size_t), 1, return NULL;);
	memset(look2, 0, sizeof(pattern_key_t) + 1*sizeof(size_t));
	look2->type = 3;
	unsigned keylen2 = 1*sizeof(size_t) + sizeof(pattern_key_t);

	for(i=1;i<r->segmentdim ;i++)
		for(j=1; j<r->segmentdim; j++)
			//if ((i!=j) && ((eqdim[i][j]==2) || (eqdim[i][j]==1)) )
			if ((i!=j) && ((eqdim[i][j]==1)) ){

#ifndef NDEBUG1
				fprintf(stdout,"\n ucons_strenghten Nodul %zu cu Nodul %zu: \n",i,j);
				fflush(stdout);
#endif

				pattern_t  * ni_dcons = NULL;
				pattern_t  * ni2_dcons = NULL;
				pattern_t  * nj_dcons = NULL;
				ap_abstract0_t* eaux = NULL;

				look->segments[0] = i;
				HASH_FIND(hh, r->udcons, look, keylen, ni_dcons);
				look->segments[0] = j;
				HASH_FIND(hh, r->udcons, look, keylen, nj_dcons);

				if(nj_dcons && (nj_dcons->dcons!= NULL) && !ap_abstract0_is_bottom(pr->man_dcons,nj_dcons->dcons)){

					//end = false;

					ap_abstract0_t* aux = NULL;
					ap_abstract0_t* eaux = NULL;

					look2->segments[0] = i;
					HASH_FIND(hh, r->udcons, look2, keylen2, ni2_dcons);

					if(eqdim[i][j]==1){
						/* strenghten hd(ni)\cup tl(ni)\subseteq hd(nj)\cup tl(nj) w.r.t. the patten forall y */
						ap_dimchange_t dimadd;
						ap_dimchange_init (&dimadd, 2, 0);
						dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
						dimadd.dim[0] = r->datadim + 2*r->segmentdim;
						dimadd.dim[1] = r->datadim + 2*r->segmentdim;

						aux = ap_abstract0_add_dimensions(pr->man_dcons,false,r->econs,&dimadd,false);

						ap_dimchange_clear (&dimadd);

						ap_lincons0_array_t arr = ap_lincons0_array_make (1);
						// d(y) - d(nj) ==0
						arr.p[0].constyp = AP_CONS_EQ;
						arr.p[0].linexpr0 =
								ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
						arr.p[0].scalar = NULL;
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
								r->datadim + j, 1);
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
								r->datadim + 2*r->segmentdim + 1, -1);
						ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

						//								// ? de verif ca trebuie y == 0
						//								arr.p[1].constyp = AP_CONS_EQ;
						//								arr.p[1].linexpr0 =
						//										ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * r->segmentdim + 2);
						//								arr.p[1].scalar = NULL;
						//								ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
						//										r->datadim + 2* r->segmentdim, 1);
						//								ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);


						aux = ap_abstract0_meet_lincons_array (pr->man_dcons, true, aux, &arr);

						ap_lincons0_array_clear (&arr);

						ap_abstract0_t* aux2 = ap_abstract0_copy(pr->man_dcons,nj_dcons->dcons);

						//								ap_dim_t * fdim = (ap_dim_t *) malloc ((r->segmentdim +1) * sizeof (ap_dim_t));
						//								size_t kk;
						//								for(kk=0; kk<r->segmentdim ;kk++)
						//									fdim[kk] = r->datadim + r->segmentdim + kk;
						//								fdim[kk] = r->datadim + r->segmentdim + kk;

						ap_dim_t * fdim = (ap_dim_t *) malloc (1 * sizeof (ap_dim_t));
						fdim[0] = r->datadim + 2*r->segmentdim;
						/* forget the length constraints */
						aux2 = ap_abstract0_forget_array(pr->man_dcons,true,aux2,
								fdim,1,false);
						free(fdim);
						fdim = NULL;

						aux = ap_abstract0_join(pr->man_dcons,true, aux,aux2);

						arr = ap_lincons0_array_make (2);
						//  y - 1 >= 0
						arr.p[0].constyp = AP_CONS_SUPEQ;
						arr.p[0].linexpr0 =
								ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * r->segmentdim + 2);
						arr.p[0].scalar = NULL;
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
								r->datadim + 2* r->segmentdim, 1);
						ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);
						//  l[ni] - y - 1 >= 0
						arr.p[1].constyp = AP_CONS_SUPEQ;
						arr.p[1].linexpr0 =
								ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * r->segmentdim + 2);
						arr.p[1].scalar = NULL;
						ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
								r->datadim + r->segmentdim + i, 1);
						ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
								r->datadim + 2* r->segmentdim, -1);
						ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, -1);

						aux = ap_abstract0_meet_lincons_array (pr->man_dcons, true, aux, &arr);

						ap_lincons0_array_clear (&arr);

					}
					//						else{
					//							/* strenghten hd(ni)\cup tl(ni)\subseteq tl(nj) w.r.t. the patten forall y */
					//							aux = ap_abstract0_copy(pr->man_dcons,nj_dcons->dcons);
					//						}

					/* hd(ni)\subseteq tl(nj)\cup hd(nj)  */

					ap_lincons0_array_t arr = ap_lincons0_array_make (1);
					// d(y) - d(ni) ==0
					arr.p[0].constyp = AP_CONS_EQ;
					arr.p[0].linexpr0 =
							ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * r->segmentdim + 2);
					arr.p[0].scalar = NULL;
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + i, 1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + 2*r->segmentdim + 1, -1);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

					eaux = ap_abstract0_meet_lincons_array (pr->man_dcons, false, nj_dcons->dcons, &arr);

					ap_lincons0_array_clear (&arr);

					ap_linexpr0_t *expr_y =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
					ap_linexpr0_set_cst_scalar_int (expr_y, 0);

					ap_linexpr0_set_coeff_scalar_int (expr_y,r->datadim + i , 1);
					ap_dim_t dy = r->datadim + 2*r->segmentdim + 1;

					eaux = ap_abstract0_substitute_linexpr (pr->man_dcons, false, nj_dcons->dcons,dy , expr_y, NULL);
					ap_linexpr0_free(expr_y);

					ap_dimchange_t dimrm;
					dimrm.realdim = 0;
					dimrm.intdim = 2;
					dimrm.dim = (ap_dim_t *) malloc (dimrm.intdim * sizeof (ap_dim_t));
					dimrm.dim[0] =  r->datadim + 2*r->segmentdim;
					dimrm.dim[1] =  r->datadim + 2*r->segmentdim +1;


					eaux = ap_abstract0_remove_dimensions (pr->man_dcons, true, eaux, &dimrm);
					free(dimrm.dim);

					//						ap_dim_t * fdim = (ap_dim_t *) malloc (r->segmentdim * sizeof (ap_dim_t));
					//						for(size_t i=0; i<r->segmentdim ;i++)
					//							fdim[i] = r->datadim + r->segmentdim + i;
					//						/* forget the length constraints */
					//						eaux = ap_abstract0_forget_array(pr->man_dcons,true,eaux,
					//								fdim,r->segmentdim,true);
					//						free(fdim);
					//						fdim = NULL;


					if(eqdim[i][j]==1){
						/* hd(ni)\subseteq tl(nj) \CUP HD(nj)  */
						ap_abstract0_t *eaux2=NULL;

						ap_lincons0_array_t arr = ap_lincons0_array_make (1);
						// d(nj) - d(ni) ==0
						arr.p[0].constyp = AP_CONS_EQ;
						arr.p[0].linexpr0 =
								ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * r->segmentdim );
						arr.p[0].scalar = NULL;
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
								r->datadim + i, 1);
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
								r->datadim + j, -1);
						ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);
						eaux2 = ap_abstract0_meet_lincons_array (pr->man_dcons, false,r->econs, &arr);
						ap_lincons0_array_clear (&arr);



						ap_linexpr0_t *expr_y =  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
						ap_linexpr0_set_cst_scalar_int (expr_y, 0);

						ap_linexpr0_set_coeff_scalar_int (expr_y,r->datadim + i , 1);
						ap_dim_t dnj = r->datadim + j;
						eaux2 = ap_abstract0_substitute_linexpr (pr->man_dcons, false, r->econs,dnj , expr_y, NULL);
						ap_linexpr0_free(expr_y);
#ifndef NDEBUG1
						fprintf(stdout,"\n ucons_strenghten problem dupa  a->econs: \n");
						fflush(stdout);
#endif

						eaux = ap_abstract0_join(pr->man_dcons,true,eaux,eaux2);
#ifndef NDEBUG1
						fprintf(stdout,"\n ucons_strenghten problem dupa  join: \n");
						fflush(stdout);
#endif
						//insure equal lengths

						//ap_lincons0_array_t
						arr = ap_lincons0_array_make (1);
						// l(nj) - l(ni) >=0
						arr.p[0].constyp = AP_CONS_SUPEQ;
						arr.p[0].linexpr0 =
								ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * r->segmentdim );
						arr.p[0].scalar = NULL;
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
								r->datadim + r->segmentdim + j, 1);
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
								r->datadim + r->segmentdim + i, -1);
						ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

						eaux = ap_abstract0_meet_lincons_array (pr->man_dcons, true,eaux, &arr);

						ap_lincons0_array_clear (&arr);

						//							/*  ms(ni) inclus ms(nj)
						//							 * if (len(ni) == 1 and len(nj) == 1
						//							 * then hd(ni) == hd(nj)
						//							 * */
						//							if(test_singleton(pr->man_dcons,r->econs,r->datadim,r->segmentdim,i)
						//									&&test_singleton(pr->man_dcons,r->econs,r->datadim,r->segmentdim,j)){
						//								arr = ap_lincons0_array_make (1);
						//								// d(nj) - d(ni) ==0
						//								arr.p[0].constyp = AP_CONS_EQ;
						//								arr.p[0].linexpr0 =
						//										ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * r->segmentdim );
						//								arr.p[0].scalar = NULL;
						//								ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
						//										r->datadim + j, 1);
						//								ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
						//										r->datadim + i, -1);
						//								ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);
						//
						//								eaux = ap_abstract0_meet_lincons_array (pr->man_dcons, false,r->econs, &arr);
						//
						//								ap_lincons0_array_clear (&arr);
						//							}




					}
					//						else{
					//
					//#ifndef NDEBUG1
					//						fprintf(stdout,"\n ucons_strenghten problem a->econs: \n");
					//						ucons_fprint(stdout, pr->man, r, NULL);
					//						fprintf(stdout,"\n");
					//						fflush(stdout);
					//#endif
					//							eaux = ap_abstract0_join(pr->man_dcons,false,r->econs,eaux);
					//						}

#ifndef NDEBUG1
					fprintf(stdout,"\n ucons_strenghten before meet a->econs: \n");
					ucons_fprint_econs(stdout, pr->man, r, NULL);
					fprintf(stdout,"\n");
					fflush(stdout);
#endif
					if(eaux)
						r->econs = ap_abstract0_meet(pr->man_dcons,true,eaux,r->econs);
#ifndef NDEBUG1
					fprintf(stdout,"\n ucons_strenghten after meet a->econs: \n");
					ucons_fprint_econs(stdout, pr->man, r, NULL);
					fprintf(stdout,"\n");
					fflush(stdout);
#endif

					/* join the 2 universals */
					ap_abstract0_t *auxe1 = ap_abstract0_copy(pr->man_dcons,r->econs);
					ap_dimchange_t dimadd1;

					ap_dimchange_init (&dimadd1, 2, 0);
					dimadd1.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
					dimadd1.dim[0] = r->datadim + 2*r->segmentdim + 0;
					dimadd1.dim[1] = r->datadim + 2*r->segmentdim + 0;

					auxe1 = ap_abstract0_add_dimensions(pr->man_dcons,true,auxe1,&dimadd1,false);
					ap_dimchange_clear (&dimadd1);

					ap_abstract0_t *aux_y = ap_abstract0_copy(pr->man_dcons,aux);

					if(ni_dcons){
						//ap_abstract0_t * aux2 = ap_abstract0_copy(pr->man_dcons,aux);
						ni_dcons->dcons = ap_abstract0_meet(pr->man_dcons,true,
								ni_dcons->dcons,aux);
						ni_dcons->dcons = ap_abstract0_meet(pr->man_dcons,true,
								ni_dcons->dcons,auxe1);
					}
					else{
						//end = false;
						checked_malloc(ni_dcons,pattern_t,1,sizeof(pattern_t)+1*sizeof(size_t),return NULL;);
						memset(ni_dcons, 0, sizeof(pattern_t)+1*sizeof(size_t));
						ni_dcons->key.type = look->type;
						look->segments[0] = i;

						ni_dcons->key.segments[0] = look->segments[0];
						ni_dcons->dcons = ap_abstract0_copy(pr->man_dcons,aux);
						ni_dcons->dcons = ap_abstract0_meet(pr->man_dcons,true,ni_dcons->dcons,auxe1);

						ap_abstract0_free(pr->man_dcons,aux);
						HASH_ADD(hh,r->udcons,key,keylen,ni_dcons);
						r=add_pattern_n2p(pr,r,look);
					}



					/* computation for y1<=y2 constraint */
					ap_abstract0_t *auxe2 = ap_abstract0_copy(pr->man_dcons,r->econs);

					ap_dimchange_t dimadd2;
					ap_dimchange_init (&dimadd2, 4, 0);
					dimadd2.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
					dimadd2.dim[0] = r->datadim + 2*r->segmentdim + 0;
					dimadd2.dim[1] = r->datadim + 2*r->segmentdim + 0;
					dimadd2.dim[2] = r->datadim + 2*r->segmentdim + 0;
					dimadd2.dim[3] = r->datadim + 2*r->segmentdim + 0;

					auxe2 = ap_abstract0_add_dimensions(pr->man_dcons,true,auxe2,&dimadd2,false);
					ap_dimchange_clear (&dimadd2);


					//add dims  to aux before intersection
					//add y1 to renforce y2
					ap_dimchange_t dimadd;

					ap_dimchange_init (&dimadd, 2, 0);
					dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
					dimadd.dim[0] = r->datadim + 2*r->segmentdim + 1;
					dimadd.dim[1] = r->datadim + 2*r->segmentdim + 2;
					ap_abstract0_t * aux3 = ap_abstract0_add_dimensions(pr->man_dcons,false,aux_y,&dimadd,false);
					ap_dimchange_clear (&dimadd);

					//add y2 to renforce y1
					ap_dimchange_init (&dimadd, 2, 0);
					dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
					dimadd.dim[0] = r->datadim + 2*r->segmentdim + 0;
					dimadd.dim[1] = r->datadim + 2*r->segmentdim + 1;
					aux_y = ap_abstract0_add_dimensions(pr->man_dcons,true,aux_y,&dimadd,false);
					ap_dimchange_clear (&dimadd);

					if(ni2_dcons){

#ifndef NDEBUG1
						fprintf(stdout,"\n ucons_strenghten debug before strengthen P12: \n");
						fprintf(stdout," strengthen: \n");
						//						look->segments[0] = j;
						//						ucons_fprint_dcons(stdout, pr->man, r, NULL);
						fprintf(stdout,"\n");
						fflush(stdout);
#endif
						ni2_dcons->dcons = ap_abstract0_meet(pr->man_dcons,true,
								ni2_dcons->dcons,aux3);

						ni2_dcons->dcons = ap_abstract0_meet(pr->man_dcons,true,
								ni2_dcons->dcons,aux_y);

						ni2_dcons->dcons = ap_abstract0_meet(pr->man_dcons,true,
								ni2_dcons->dcons,auxe2);

#ifndef NDEBUG1
						fprintf(stdout,"\n results: \n");
						ap_abstract0_fprint(stdout, pr->man_dcons, ni2_dcons->dcons, NULL);
						fprintf(stdout,"\n");
						fflush(stdout);
#endif

					}
					else{
						//end = false;
						checked_malloc(ni2_dcons,pattern_t,1,sizeof(pattern_t)+1*sizeof(size_t),return NULL;);
						memset(ni2_dcons, 0, sizeof(pattern_t)+1*sizeof(size_t));
						ni2_dcons->key.type = look2->type;
						look2->segments[0] = i;

						ni2_dcons->key.segments[0] = look2->segments[0];

						aux_y = ap_abstract0_meet(pr->man_dcons,true, aux3, aux_y);

						ap_lincons0_array_t arr2 = ap_lincons0_array_make (1);
						// y2 - y1 >=0
						arr2.p[0].constyp = AP_CONS_SUPEQ;
						arr2.p[0].linexpr0 =
								ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * r->segmentdim +4);
						arr2.p[0].scalar = NULL;
						ap_linexpr0_set_coeff_scalar_int (arr2.p[0].linexpr0,
								r->datadim + 2*r->segmentdim , -1);
						ap_linexpr0_set_coeff_scalar_int (arr2.p[0].linexpr0,
								r->datadim + 2*r->segmentdim + 1, 1);
						ap_linexpr0_set_cst_scalar_int (arr2.p[0].linexpr0, 0);

						aux_y = ap_abstract0_meet_lincons_array (pr->man_dcons, true,aux_y, &arr2);

						ap_lincons0_array_clear (&arr2);


						ni2_dcons->dcons = ap_abstract0_copy(pr->man_dcons,aux_y);
						ni2_dcons->dcons = ap_abstract0_meet(pr->man_dcons,true, ni2_dcons->dcons,auxe2);

						HASH_ADD(hh,r->udcons,key,keylen,ni2_dcons);
						r=add_pattern_n2p(pr,r,look2);
						ap_abstract0_free(pr->man_dcons, aux_y);
					}

#ifndef NDEBUG1
					fprintf(stdout,"\n ucons_strenghten after meet a->udcons: \n");
					ucons_fprint(stdout, pr->man, r, NULL);
					fprintf(stdout,"\n");
					fflush(stdout);
#endif

				}
				else{
					//insure equal lengths

					ap_lincons0_array_t arr = ap_lincons0_array_make (1);
					// l(nj) - l(ni) >=0
					arr.p[0].constyp = AP_CONS_SUPEQ;
					arr.p[0].linexpr0 =
							ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * r->segmentdim );
					arr.p[0].scalar = NULL;
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + r->segmentdim + j, 1);
					ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
							r->datadim + r->segmentdim + i, -1);
					ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

					r->econs = ap_abstract0_meet_lincons_array (pr->man_dcons, true,r->econs, &arr);

					ap_lincons0_array_clear (&arr);

					/*  ms(ni) inclus ms(nj)
					 * if (len(ni) == 1 and len(nj) == 1
					 * then hd(ni) == hd(nj)
					 * */
					if(test_singleton(pr->man_dcons,r->econs,r->datadim,r->segmentdim,i)
							&&test_singleton(pr->man_dcons,r->econs,r->datadim,r->segmentdim,j)){

#ifndef NDEBUG1
						fprintf(stdout,"\n ucons_strenghten singletons i=%zu and j=%zu \n", i, j);
						fprintf(stdout,"\n");
						fflush(stdout);
#endif
						arr = ap_lincons0_array_make (1);
						// d(nj) - d(ni) ==0
						arr.p[0].constyp = AP_CONS_EQ;
						arr.p[0].linexpr0 =
								ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * r->segmentdim );
						arr.p[0].scalar = NULL;
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
								r->datadim + j, 1);
						ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
								r->datadim + i, -1);
						ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

						r->econs = ap_abstract0_meet_lincons_array (pr->man_dcons, true,r->econs, &arr);

						ap_lincons0_array_clear (&arr);
					}


				}
				/// TODO MODIFSSSS
				//	eqdim[i][j]=0;
			}
	if(look) free(look);
	look = NULL;
	if(look2) free(look2);
	look2 = NULL;
	//}
#ifndef NDEBUG1
	fprintf(stdout,"\n ucons_strenghten returns: \n");
	ucons_fprint(stdout, pr->man, r, NULL);
	fprintf(stdout,"\n");
	fflush(stdout);
#endif
	return r;
}

void close_eq(size_t ** tdim, size_t size){
	size_t i,j,k;
	size_t n = size;
	bool b = false;
	while(!b){
		b=true;
		for(i=1; i < n; i++)
			for(j=1; j < n; j++){
				if(tdim[i][j]==1){
					for(k=1; k < n; k++){
						if(tdim[j][k]==1) {
							if(tdim[i][k]==0) {
								tdim[i][k]= 1;
								b=false;
							}
							if(tdim[k][i]==0) {
								tdim[k][i]= 1;
								b=false;
							}
						}
					}
				}
			}
	}
}

ucons_t*
ucons_saturation(ucons_internal_t* pr, ucons_t* a){

	pattern_t *s =NULL;
	size_t u_seg,e_seg,nr_y;
	unsigned keylen;
	pattern_t *rt=NULL;

	for(s=a->udcons; s!=NULL; s=s->hh.next){

		//check definitions y1>=1....etc
		u_seg=pr->PI[s->key.type].u_seg;
		e_seg=pr->PI[s->key.type].e_seg;
		nr_y=pr->PI[s->key.type].nr_y;

		ap_lincons0_array_t arr = ap_lincons0_array_make (nr_y);
		for(size_t i=0; i<nr_y ;i++){
			// yi - 1 >= 0
			arr.p[i].constyp = AP_CONS_SUPEQ;
			arr.p[i].linexpr0 =
					ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim +
							2 * a->segmentdim + 2*nr_y);
			arr.p[i].scalar = NULL;
			ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, a->datadim +
					2 * a->segmentdim + i, 1);
			ap_linexpr0_set_cst_scalar_int (arr.p[i].linexpr0, -1);

		}
		s->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons,
				true, s->dcons, &arr);
		ap_lincons0_array_clear (&arr);

		if(pr->PI[s->key.type].kind==pattern_1_2){
			/* y2 >= y1 */
			arr = ap_lincons0_array_make (1);
			arr.p[0].constyp = AP_CONS_SUPEQ;
			arr.p[0].linexpr0 =
					ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim +
							2 * a->segmentdim + 2*nr_y);
			arr.p[0].scalar = NULL;
			ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, a->datadim +
					2 * a->segmentdim, -1);
			ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, a->datadim +
							2 * a->segmentdim+1, 1);
			ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

			s->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons,
							true, s->dcons, &arr);
			ap_lincons0_array_clear (&arr);

			/* meet with \forall y */
			keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);
			pattern_key_t *look=NULL;
			checked_malloc(look, pattern_key_t, sizeof(pattern_key_t) + u_seg*sizeof(size_t), 1, return NULL;);
			memset(look, 0, sizeof(pattern_key_t) + u_seg*sizeof(size_t));
			look->type = 0; /// PI[0] \forall y
			look->segments[0] = s->key.segments[0];
			HASH_FIND(hh, a->udcons, look, keylen, rt);

			if(rt!=NULL && rt->dcons!=NULL){

				ap_abstract0_t * aux = ap_abstract0_copy(pr->man_dcons,rt->dcons);

				ap_dimchange_t dimadd;
				ap_dimchange_init (&dimadd, 2, 0);
				dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
				dimadd.dim[0] = a->datadim + 2*a->segmentdim;
				dimadd.dim[1] = a->datadim + 2*a->segmentdim + 1;
				aux = ap_abstract0_add_dimensions(pr->man_dcons, true, aux, &dimadd, false);

				ap_dimchange_clear (&dimadd);

				s->dcons = ap_abstract0_meet(pr->man_dcons, true, s->dcons, aux);

				aux = ap_abstract0_copy(pr->man_dcons,rt->dcons);

				ap_dimchange_init (&dimadd, 2, 0);
				dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
				dimadd.dim[0] = a->datadim + 2*a->segmentdim + 1;
				dimadd.dim[1] = a->datadim + 2*a->segmentdim + 2 ;
				aux = ap_abstract0_add_dimensions(pr->man_dcons,true,aux,&dimadd,false);

				ap_dimchange_clear (&dimadd);

				s->dcons = ap_abstract0_meet(pr->man_dcons,true,s->dcons,aux);


			}
		}
		//TODO : all patterns
	}


	return a;
}
