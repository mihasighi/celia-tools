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


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

long int lrand48(void);
void srand48(long int seedval);

#include "shape.h"
#include "shape_fun.h"
#include "shape_internal.h"
#include "apron2shape.h"
#include "shape_macros.h"

ap_manager_t* ms; /* shape manager */

shape_internal_t* pr;

size_t realdim = 3; /* vars x,xi,y; add 3 for yi,z,zi */
size_t intdim = 3;  /* vars _l, _k, S; add 1 for m */
size_t SIZE = 4;     /* nb of nodes */

char* names_dim[9] = { "_l", "_k", "d", "x", "xi", "y", "yi", "z", "zi" };
#define DIM_X (intdim)
#define DIM_Y (intdim+2)

#define MAXN 50000
int N = 7;

char b1_[MAXN+4]= " [";
char b2_[MAXN+4];
int i_;

typedef enum { none = 0, best = 1, exact = 2 } exactness;

exactness flag;

#define LOOP								\
  {									\
    memset(b1_+2,' ',N); b1_[N+2]=']'; b1_[N+3]=0;			\
    memset(b2_,8,N+3); b2_[N+3]=0;					\
    for (i_=0;i_<N;i_++) {						\
      if (N<80) printf("%s%s",b1_,b2_);					\
      else if (i_) printf("%c",b1_[i_+1]);				\
      fflush(stdout);							\
      flag = exact;							\

#define FLAG(m)						    \
  if (m->result.flag_exact) ;				    \
  else if (m->result.flag_best && flag==exact) flag = best; \
  else flag = none;

#define RESULT(c) b1_[i_+2]=c

#define ERROR_RESULT(msg)	ERROR(msg,RESULT('!');)

#define ENDLOOP					\
  } }						\
  if (N<80) printf("%s%s\n",b1_,b2_);		\
  else printf("%c",b1_[i_+2])


void print_shape(const char* msg, shape_t* a)
{
  fprintf(stdout,"%s = ",msg);
  shape_fprint(stdout,ms,a,names_dim);
  //shape_fdump(stdout,ms,a);
  fprintf(stdout,"\n");
}




/* chose one linexpr within an interval linexpr */
/* TOD0: used? */
ap_linexpr0_t* random_from_linexpr(ap_linexpr0_t* a)
{
  size_t i;
  ap_linexpr0_t* l = ap_linexpr0_alloc(AP_LINEXPR_DENSE,a->size);
  assert(a->discr==AP_LINEXPR_DENSE);
  for (i=0;i<a->size;i++) {
    switch (a->p.coeff[i].discr) {
    case AP_COEFF_SCALAR:
      ap_coeff_set_scalar(l->p.coeff+i,a->p.coeff[i].val.scalar);
      break;
    case AP_COEFF_INTERVAL:
      if (lrand48()%10>8)
        ap_coeff_set_scalar(l->p.coeff+i,a->p.coeff[i].val.interval->inf);
      else
        ap_coeff_set_scalar(l->p.coeff+i,a->p.coeff[i].val.interval->sup);
      break;

    }
  }  switch (a->cst.discr) {
  case AP_COEFF_SCALAR:
    ap_coeff_set_scalar(&l->cst,a->cst.val.scalar);
    break;
  case AP_COEFF_INTERVAL:
    if (lrand48()%10>8) ap_coeff_set_scalar(&l->cst,a->cst.val.interval->inf);
    else ap_coeff_set_scalar(&l->cst,a->cst.val.interval->sup);
    break;

  }
 return l;
}

/* TODO: used? */
ap_lincons0_t random_from_lincons(ap_lincons0_t a)
{
  return ap_lincons0_make(a.constyp,random_from_linexpr(a.linexpr0),NULL);
}


/* ********************************* */
/*             infos                 */
/* ********************************* */



void info(void)
{
  printf("shapes:   %s (%s)\n",ms->library,ms->version);
  /*
  printf("nums:      %s (%s wto overflow,%s)\n",NUM_NAME,
	 num_safe ? "sound" : "unsound",
	 num_incomplete ? "incomplete" : "complete");
  */
}



/* ********************************* */
/*           various tests           */
/* ********************************* */

void test_misc(void)
{
  int i;
  shape_t* bot = shape_bottom(ms,intdim,realdim);
  shape_t* top = shape_top(ms,intdim,realdim);
  ap_dimension_t d1 = shape_dimension(ms,bot);
  ap_dimension_t d2 = shape_dimension(ms,top);
  printf("\n**********\n* test_misc: performing various tests\n**********\n");
  if (shape_check(pr,bot) != '.')	  printf("shape_bottom failed #0\n");
  if (shape_check(pr,top) != '.')	  printf("shape_top failed #0\n");
  if (d1.intdim!=intdim || d1.realdim!=realdim)      printf("shape_dimension failed #1\n");
  if (d2.intdim!=intdim ||
      d2.realdim!=realdim)         printf("shape_dimension failed #2\n");
  if (!shape_is_bottom(ms,bot))  printf("shape_is_bottom failed #3\n");
  if (shape_is_bottom(ms,top))   printf("shape_is_bottom failed #4\n");
  if (shape_is_top(ms,bot))      printf("shape_is_top failed #5\n");
  if (!shape_is_top(ms,top))     printf("shape_is_top failed #6\n");
  if (!shape_is_leq(ms,bot,top)) printf("shape_is_leq failed #7\n");
  if (shape_is_leq(ms,top,bot))  printf("shape_is_leq failed #8\n");
  if (!shape_is_eq(ms,bot,bot))  printf("shape_is_eq failed #9\n");
  if (!shape_is_eq(ms,top,top))  printf("shape_is_eq failed #10\n");
  if (shape_is_eq(ms,bot,top))   printf("shape_is_eq failed #11\n");
  if (shape_is_dimension_unconstrained(ms,bot,DIM_X))
    printf("shape_is_dimension_unconstrained #12\n");
  if (!shape_is_dimension_unconstrained(ms,top,DIM_X))
    printf("shape_is_dimension_unconstrained #13\n");
  /* Test of predefined shapes used in assume. */
  for (i=0;i<5;i++) {
    printf("***** shape_make %d\n", i);
    shape_t* o = shape_make(pr, i, intdim,realdim);
    shape_t* c = shape_copy(ms,o);
    shape_t* l = shape_closure(ms,false,o);
    ap_dimension_t d = shape_dimension(ms,o);
    if (d.intdim<intdim || d.realdim!=realdim) printf("shape_dimension failed #18\n");
    if (!shape_is_leq(ms,bot,o))  printf("shape_is_leq failed #19\n");
    if (!shape_is_leq(ms,o,top))  printf("shape_is_leq failed #20\n");
    if (!shape_is_eq(ms,o,c))     printf("shape_is_eq failed #21\n");
    shape_size(ms,o);
    shape_close(pr,o);
    // not implemented
    //shape_minimize(ms,o);
    //shape_canonicalize(ms,o);
    //shape_approximate(ms,o,0);
    //shape_is_minimal(ms,o);
    //shape_is_canonical(ms,o);
    shape_free(ms,o); shape_free(ms,c); shape_free(ms,l);
  }
  /* Test of random shapes (not used). */
  for (i=0;i<N;i++) {
    shape_t* o = shape_random(pr, SIZE, intdim,realdim);
    shape_t* c = shape_copy(ms,o);
    shape_t* l = shape_closure(ms,false,o);
    ap_dimension_t d = shape_dimension(ms,o);
    if (d.intdim<intdim || d.realdim!=realdim) printf("shape_dimension failed #14\n");
    if (!shape_is_leq(ms,bot,o))  printf("shape_is_leq failed #15\n");
    if (!shape_is_leq(ms,o,top))  printf("shape_is_leq failed #16\n");
    if (!shape_is_eq(ms,o,c))     printf("shape_is_eq failed #17\n");
    shape_size(ms,o);
    shape_close(pr,o);
    // not implemented
    //shape_minimize(ms,o);
    //shape_canonicalize(ms,o);
    //shape_approximate(ms,o,0);
    //shape_is_minimal(ms,o);
    //shape_is_canonical(ms,o);
    shape_free(ms,o); shape_free(ms,c); shape_free(ms,l);
  }
  shape_free(ms,bot); shape_free(ms,top);
}



/* ********************************* */
/*           closure                 */
/* ********************************* */


void test_closure(void)
{
  // printf("\nclosure %s\n",num_incomplete?"":"(c,o expected)");
  LOOP {
    shape_t* o = shape_random(pr, SIZE, intdim, realdim);
    shape_close(pr,o);
    RESULT(shape_check(pr,o));
    shape_free(ms,o);
  } ENDLOOP;
}



/* ********************************* */
/*            conversions            */
/* ********************************* */

void
test_lincons_conversion(exprmode_t mode)
{
  printf("\nconversion from %slincons\n",exprname[mode]);
  LOOP {
    int dim = 6, nb = 10, i;
    ap_abstract0_t *p, *p2;
    shape_t *o, *oo;
    ap_lincons0_array_t t = ap_lincons0_array_make(nb);
    ap_lincons0_array_t tt = ap_lincons0_array_make(nb);
    for (i=0;i<nb;i++) {
      t.p[i] = ap_lincons0_make((lrand48()%100>=90)?AP_CONS_EQ:
				(lrand48()%100>=90)?AP_CONS_SUP:
				AP_CONS_SUPEQ,
				shape_linexpr_random(pr, mode,intdim,realdim),
				NULL);
      tt.p[i] = random_from_lincons(t.p[i]);
    }
    o  = shape_top(ms,intdim,realdim);
    oo  = shape_meet_lincons_array(ms,true,o,&t); FLAG(ms);
    RESULT(shape_check(pr,o));
    if (shape_is_eq(ms,o,oo)) RESULT('x');
    if (shape_is_leq(ms,oo,o)) ERROR_RESULT("best flag");
    shape_free(ms,o); shape_free(ms,oo);
    ap_lincons0_array_clear(&t); ap_lincons0_array_clear(&tt);
  } ENDLOOP;
}


/* ********************************* */
/*               meet                */
/* ********************************* */

void test_meet(void)
{
  printf("\nmeet %s\n","(* expected)");
  LOOP {
    int dim = 6;
    shape_t *o1, *o2, *o;
    o1 = shape_random(pr, SIZE,intdim,realdim);
    o2 = shape_random(pr, SIZE,intdim,realdim);
    o  = shape_meet(ms,false,o1,o2); FLAG(ms);
    RESULT(shape_check(pr,o)); shape_check(pr,o1); shape_check(pr,o2);
    if (!shape_is_leq(ms,o,o1) || !shape_is_leq(ms,o,o2)) {
      ERROR_RESULT("not lower bound");
      print_shape("o1",o1); print_shape("o2",o2); print_shape("o",o);
    }
    shape_free(ms,o); shape_free(ms,o1); shape_free(ms,o2);
  } ENDLOOP;
  printf("\nmeet top %s\n","(* expected)");
  LOOP {
    int dim = 8;
    shape_t *o1, *o2, *o, *oo;
    o1 = shape_random(pr, SIZE, intdim,realdim);
    o2 = shape_top(ms,intdim,realdim);
    o  = shape_meet(ms,false,o1,o2);
    oo = shape_meet(ms,false,o,o1);
    shape_check(pr,o1); shape_check(pr,o2); shape_check(pr,o);
    if (!shape_is_eq(ms,o,o1)) {
      ERROR_RESULT("not eq #1");
      print_shape("o1",o1); print_shape("o",o);
    }
    else if (!shape_is_eq(ms,o,oo)) {
      ERROR_RESULT("not eq #2");
      print_shape("o1",o1); print_shape("o",o); print_shape("oo",oo);
    }
    else RESULT('*');
    shape_free(ms,o); shape_free(ms,o1); shape_free(ms,o2); shape_free(ms,oo);
  } ENDLOOP;
  printf("\nmeet bot %s\n","(* expected)");
  LOOP {
    int dim = 8;
    shape_t *o1, *o2, *o;
    o1 = shape_random(pr, SIZE,intdim,realdim);
    o2 = shape_bottom(ms,intdim,realdim);
    o  = shape_meet(ms,false,o1,o2);
    shape_check(pr,o1); shape_check(pr,o2); shape_check(pr,o);
    if (!shape_is_bottom(ms,o)) {
      ERROR_RESULT("not bottom");
      print_shape("o1",o1); print_shape("o",o);
    }
    else RESULT('*');
    shape_free(ms,o); shape_free(ms,o1); shape_free(ms,o2);
  } ENDLOOP;
}

#define NB_MEET 5
void test_meet_array(void)
{
  printf("\nmeet array %s\n","(* expected)");
  LOOP {
    int i, dim = 6;
    shape_t* o[NB_MEET], *oo;
    for (i=0;i<NB_MEET;i++)
      o[i] = shape_random(pr,SIZE,intdim,realdim);
    oo = shape_meet_array(ms,o,NB_MEET); FLAG(ms);
    RESULT(shape_check(pr,oo));
    for (i=0;i<NB_MEET;i++)
      if (!shape_is_leq(ms,oo,o[i])) ERROR_RESULT("not lower bound");
    for (i=0;i<NB_MEET;i++) shape_free(ms,o[i]);
    shape_free(ms,oo);
  } ENDLOOP;
}

shape_t*
test_add_pcons_once(exprmode_t mode, shape_t* o, pcons0_array_t* array)
{
    size_t i, o2size;
    shape_t* o1;
    shape_t* o2;
    printf("\n*** meet with pcons:\n");
#ifndef NDEBUG
    shape_fdump(stdout,pr->man,o);
    shape_pcons_array_fdump(stdout,array);
#endif
    /* apply domain function (underapproximation) */
    printf("\n*** (underapproximation):\n");
    o1 = shape_copy_internal(pr,o);
    o2 = shape_meet_pcons_array(pr,true,o1,array); FLAG(ms);
    shape_check(pr,o);
    if (o2) {
      RESULT(shape_check(pr,o2));
      if (!shape_is_leq(ms,o2,o)) {
        ERROR_RESULT("not included in");
        print_shape("o",o); print_shape("o2",o2);
          }
      shape_free(pr->man,o2);
    }

    /* apply internal function (exact) */
    printf("\n*** (exact):\n");
    o2 =  shape_meet_pcons_array(pr, false, o, array);
    if (o2) {
      RESULT(shape_check(pr,o2));
      if (!shape_is_leq(ms,o2,o)) {
        ERROR_RESULT("not included in");
        print_shape("o",o); print_shape("o2",o2);
          }
      // shape_free(pr->man,o2);
    }
    return o2;
}

void
test_add_lincons(exprmode_t mode)
{
  size_t i, nb = 4;
  ap_lincons0_array_t ar;
#ifdef PCONS
  pcons0_array_t* arr;
#endif
  shape_t* o, *or;
  printf("\n**********\n* test_add_lincons: meet with lincons %s (* expected)\n**********\n",exprname[mode]);
  /* predefined graphs and constraints */
  /* case (1): top /\ x --> null /\ x == y */
  o = shape_top(ms,intdim,realdim);
  ar = ap_lincons0_array_make(2);
  ar.p[0] = shape_lincons_reach_x_y(DIM_X, NULL_DIM, intdim,realdim);
  ar.p[1] = shape_lincons_same_x_y(DIM_X, DIM_Y,intdim,realdim);
#ifdef PCONS
  arr = shape_pcons_array_of_lincons_array (pr, &ar, intdim, realdim);
  or = test_add_pcons_once(mode,o,arr);
  shape_pcons0_array_clear(arr);
#else
  or = shape_meet_lincons_array(ms,false,o,&ar);
#endif
  shape_free(pr->man, o);
  shape_free(pr->man, or);

  /* case (2): x-->null, y-->null /\ y != null */
  o = shape_make(pr,2, intdim, realdim);
  ar = ap_lincons0_array_make(1);
  ar.p[0] = shape_lincons_diff_x_y(DIM_Y,NULL_DIM,intdim,realdim);
#ifdef PCONS
  arr = shape_pcons_array_of_lincons_array (pr, &ar, intdim, realdim);
  or = test_add_pcons_once(mode,o,arr);
#else
  or = shape_meet_lincons_array(ms,false,o,&ar);
#endif
  shape_free(pr->man, o);
  shape_free(pr->man, or);

  /* case (3): x-->y-->null /\ y != null */
  o = shape_make(pr,1, intdim, realdim);
#ifdef PCONS
  or = test_add_pcons_once(mode,o,arr);
  shape_pcons0_array_clear(arr);
#else
  or = shape_meet_lincons_array(ms,false,o,&ar);
#endif
  shape_free(pr->man, o);
  shape_free(pr->man, or);

  /* random graphs and constraints
  LOOP {
    ar = ap_lincons0_array_make(nb);
    o = shape_random(pr, SIZE,intdim,realdim);
    if (lrand48()%10>=8) shape_close(pr,o);
    for (i=0;i<nb;i++) {
      ar.p[i] = ap_lincons0_make((lrand48()%100>=80)?AP_CONS_EQ:
                 (lrand48()%100>=80)?AP_CONS_SUP:
                 AP_CONS_SUPEQ,
                 shape_linexpr_random(pr,mode,intdim,realdim),
                 NULL);
    }
    arr = shape_pcons_array_of_lincons_array (pr, &ar, intdim, realdim);
    test_add_pcons_once(mode,o,arr);
    shape_free(ms,o);
    shape_pcons0_array_clear(arr);
   } ENDLOOP;
   */
}


void
test_add_tcons(exprmode_t mode)
{
  size_t i, nb = 4;
  ap_tcons0_array_t ar;
#ifdef PCONS
  pcons0_array_t* arr;
#endif
  shape_t* o, *om, *or;
  printf("\n**********\n* test_add_tcons: meet with lincons %s (* expected *)\n**********\n",exprname[mode]);
  /* predefined graphs and constraints */
  /* case (0): top /\ assume code 0, i.e., x==0 */
  o = shape_top(ms,intdim,realdim);
  ar = ap_tcons0_array_make(1);
  ar.p[0] = shape_tcons_x_cst(AP_CONS_EQ,DIM_X,false,0,intdim,realdim);
#ifdef PCONS
  arr = shape_pcons_array_of_tcons_array (pr, &ar,intdim, realdim);
  or = test_add_pcons_once(mode,o,arr);
  shape_pcons0_array_clear(arr);
#else
  or = shape_meet_tcons_array(ms,false,o,&ar);
#endif
  om = shape_make(pr,0,intdim, realdim);
  if (!shape_is_eq(ms,or,om)) printf("shape_meet_tcons failed #0\n");
  shape_free(ms, o);
  shape_free(ms, or);
  shape_free(ms, om);

  /* case (1): top /\ assume code 1, i.e., x==1 */
  o = shape_top(ms,intdim,realdim);
  ar = ap_tcons0_array_make(1);
  ar.p[0] = shape_tcons_x_cst(AP_CONS_EQ,DIM_X,false,1,intdim,realdim);
#ifdef PCONS
  arr = shape_pcons_array_of_tcons_array (pr, &ar, intdim, realdim);
  or = test_add_pcons_once(mode,o,arr);
  shape_pcons0_array_clear(arr);
#else
  or = shape_meet_tcons_array(ms,false,o,&ar);
#endif
  om = shape_make(pr,1, intdim, realdim);
  if (!shape_is_eq(ms,or,om)) printf("shape_meet_tcons failed #1\n");
  shape_free(ms, o);
  shape_free(ms, om);
  shape_free(ms, or);

  /* case (2): top /\ assume code 2, i.e., x==2 */
  o = shape_top(ms,intdim,realdim);
  ar = ap_tcons0_array_make(1);
  ar.p[0] = shape_tcons_x_cst(AP_CONS_EQ,DIM_X,false,2,intdim,realdim);
#ifdef PCONS
  arr = shape_pcons_array_of_tcons_array (pr, &ar, intdim, realdim);
  or = test_add_pcons_once(mode,o,arr);
  shape_pcons0_array_clear(arr);
#else
  or = shape_meet_tcons_array(ms,false,o,&ar);
#endif
  om = shape_make(pr,2, intdim, realdim);
  if (!shape_is_eq(ms,or,om)) printf("shape_meet_tcons failed #2\n");
  shape_free(ms, o);
  shape_free(ms, om);
  shape_free(ms, or);

  /* case (3): top /\ assume code 3, i.e., x==3 */
  o = shape_top(ms,intdim,realdim);
  ar = ap_tcons0_array_make(1);
  ar.p[0] = shape_tcons_x_cst(AP_CONS_EQ,DIM_X,false,3,intdim,realdim);
#ifdef PCONS
  arr = shape_pcons_array_of_tcons_array (pr, &ar, intdim, realdim);
  or = test_add_pcons_once(mode,o,arr);
  shape_pcons0_array_clear(arr);
#else
  or = shape_meet_tcons_array(ms,false,o,&ar);
#endif
  om = shape_make(pr,3, intdim, realdim);
  if (!shape_is_eq(ms,or,om)) printf("shape_meet_tcons failed #3\n");
  shape_free(ms, o);
  shape_free(ms, om);
  shape_free(ms, or);

  /* case (4): top /\ assume code 2, i.e., x==2 */
  o = shape_top(ms,intdim,realdim);
  ar = ap_tcons0_array_make(1);
  ar.p[0] = shape_tcons_x_cst(AP_CONS_EQ,DIM_X,false,2,intdim,realdim);
#ifdef PCONS
  arr = shape_pcons_array_of_tcons_array (pr, &ar, intdim, realdim);
  or = test_add_pcons_once(mode,o,arr);
  shape_pcons0_array_clear(arr);
#else
  or = shape_meet_tcons_array(ms,false,o,&ar);
#endif
  om = shape_make(pr,2, intdim, realdim);
  if (!shape_is_eq(ms,or,om)) printf("shape_meet_tcons failed #4\n");
  shape_free(ms, o);
  shape_free(ms, om);
  shape_free(ms, or);

  /* case (5): code (2) /\ x->next != null */
  o = shape_make(pr,2,intdim, realdim);
  ar = ap_tcons0_array_make(1);
  ar.p[0] = shape_tcons_diff_x_y (DIM_X, true, NULL_DIM, false, intdim,
					  realdim);
#ifdef PCONS
  arr = shape_pcons_array_of_tcons_array (pr, &ar, intdim, realdim);
  or = test_add_pcons_once(mode,o,arr);
  shape_pcons0_array_clear(arr);
  shape_free(ms, or);
#else
 // TODO: or = shape_meet_tcons_array(ms,false,o,&ar);
#endif
  shape_free(ms, o);


  /* random graphs and constraints */
  /*
  LOOP {
    ar = ap_tcons0_array_make(nb);
    o = shape_random(pr, SIZE,intdim,realdim);
    if (lrand48()%10>=8) shape_close(pr,o);
    for (i=0;i<nb;i++)
      {
      if (mode==expr_data)
        ar.p[i] = ap_tcons0_make((lrand48()%100>=80)?AP_CONS_SUP:
        AP_CONS_SUPEQ,
        shape_texpr_random(pr,mode,intdim,realdim),
        NULL);
      else
        ar.p[i] = ap_tcons0_make((lrand48()%100>=80)?AP_CONS_EQ:
        AP_CONS_DISEQ,
        shape_texpr_random(pr,expr_reach,intdim,realdim),
        NULL);
      }
    ap_tcons0_array_fprint(stdout,&ar,NULL);
    arr = shape_pcons_array_of_tcons_array (pr, &ar, intdim, realdim);
    or = test_add_pcons_once(mode,o,arr);
    shape_pcons0_array_clear(arr);
    shape_free(ms,o);
    shape_free(ms,or);
   } ENDLOOP;
   */
}

/* ********************************* */
/*     assignments / substitutions   */
/* ********************************* */

void test_assign(int subst, exprmode_t mode)
{
  ap_linexpr0_t* l;
  ap_dim_t d;
  shape_t* o, *om, *or;

  printf("\n**********\n* test_%s_linexpr: with %s (* expected)\n**********\n",
      subst ? "subst" : "assign", exprname[mode]);

  /* identity assignments */
  /* shape 1 and xi := null */
  o = shape_make(pr, 1, intdim,realdim);
  d = DIM_X+1;
  l = shape_linexpr_null(intdim,realdim);
  or = shape_assign_linexpr_array(ms,false,o,&d,&l,1,NULL);
  FLAG(ms);
  RESULT(shape_check(pr,or));
  if (!shape_is_eq(ms,o,or)) printf("shape_assign_linexpr_array failed #1\n");
  ap_linexpr0_free(l);
  shape_free(ms,or);

  /* shape 1 and xi := x */
  l = shape_linexpr_x(DIM_X,intdim,realdim);
  or = shape_assign_linexpr_array(ms,false,o,&d,&l,1,NULL);
  FLAG(ms);
  RESULT(shape_check(pr,or));
  if (shape_is_eq(ms,o,or)) printf("shape_assign_linexpr_array failed #2\n");
  ap_linexpr0_free(l);
  shape_free(ms,o);

  /* random assignements */
  /*
  LOOP {
    size_t dim = 6;
    ap_dim_t d;
    shape_t *o, *o1;
    ap_linexpr0_t* l = shape_linexpr_random(pr, mode,intdim,realdim);
    if (mode != expr_lindata || mode != expr_data)
      d = RANDOM_PTRVAR(pr,intdim,realdim);
    else
      d = lrand48()%intdim;
    o = shape_random(pr, SIZE, intdim,realdim);
    if (lrand48()%10>=5) shape_close(pr,o);
    o1 = subst ? shape_substitute_linexpr_array(ms,false,o,&d,&l,1,NULL) :
                 shape_assign_linexpr_array(ms,false,o,&d,&l,1,NULL);
    FLAG(ms);
    RESULT(shape_check(pr,o1));
    if (shape_is_eq(ms,o,o1)) ERROR_RESULT("best flag");
    shape_free(ms,o); shape_free(ms,o1);
    ap_linexpr0_free(l);
  } ENDLOOP;
  */
}

void test_assign_texpr(int subst, exprmode_t mode)
{
  ap_texpr0_t* t1, *t2, *t3, *t4, *t5;
  ap_dim_t d;
  shape_t* o, *om, *or;

  printf("\n**********\n* test_%s_texpr: with %s (* expected)\n**********\n",
      subst ? "subst" : "assign", exprname[mode]);

  /* predefined graphs and constraints */
  /* case (1): code 1 and x := null */
  /* shape 1 and xi := null */
  o = shape_make(pr, 1, intdim,realdim);
  d = DIM_X+1;
  t1 = shape_texpr_null(intdim,realdim);
  or = shape_assign_texpr_array(ms,false,o,&d,&t1,1,NULL);
  FLAG(ms);
  RESULT(shape_check(pr,or));
  if (!shape_is_eq(ms,o,or)) printf("shape_assign_texpr_array failed #1\n");
  shape_free(ms,or);

  /* case (2): shape 1 and xi := x */
  t2 = shape_texpr_x_n_next(DIM_X,0,intdim,realdim);
  or = shape_assign_texpr_array(ms,false,o,&d,&t2,1,NULL);
  FLAG(ms);
  RESULT(shape_check(pr,or));
  if (shape_is_eq(ms,o,or)) printf("shape_assign_linexpr_array failed #2\n");
  //ap_texpr0_free(t2);
  shape_free(ms,o);
  o = or;
  or = NULL;

  /* case (3): shape 1 and xi->data := (xi->data+2) */
  t3 = shape_texpr_x_y_v_cst_data(1,DIM_X+1,0,0,0,0,2,true,intdim,realdim);
  or = shape_assign_texpr_array(ms,false,o,&d,&t3,1,NULL);
  FLAG(ms);
  RESULT(shape_check(pr,o));
  if (shape_is_eq(ms,o,or)) printf("shape_assign_texpr_array failed #3\n");
  //ap_texpr0_free(t3);
  shape_free(ms,o);
  o = or;
  or = NULL;

  /* case (4): case (3) and y := xi->next */
  d = DIM_Y;
  t4 = shape_texpr_x_n_next(DIM_X+1,1,intdim,realdim);
  or = shape_assign_texpr_array(ms,false,o,&d,&t4,1,NULL);
  FLAG(ms);
  RESULT(shape_check(pr,o));
  if (shape_is_eq(ms,o,or)) printf("shape_assign_texpr_array failed #4\n");
  //ap_texpr0_free(t4);
  shape_free(ms,o);
  o = or;
  or = NULL;

  /* case (5): case (4) and xi := null */
  d = DIM_X+1;
  or = shape_assign_texpr_array(ms,false,o,&d,&t1,1,NULL);
  FLAG(ms);
  RESULT(shape_check(pr,or));
  if (shape_is_eq(ms,o,or)) printf("shape_assign_texpr_array failed #5\n");
  shape_free(ms,o);
  o = or;
  or = NULL;

  /* case (6): case (5) and xi := y */
  t5 = shape_texpr_x_n_next(DIM_Y,0,intdim,realdim);
  or = shape_assign_texpr_array(ms,false,o,&d,&t5,1,NULL);
  FLAG(ms);
  RESULT(shape_check(pr,or));
  if (shape_is_eq(ms,o,or)) printf("shape_assign_texpr_array failed #6\n");
  shape_free(ms,o);
  o = or;
  or = NULL;

  /* case (7): case (6) and y := _null */
  d = DIM_Y;
  or = shape_assign_texpr_array(ms,false,o,&d,&t1,1,NULL);
  FLAG(ms);
  RESULT(shape_check(pr,or));
  if (shape_is_eq(ms,o,or)) printf("shape_assign_texpr_array failed #6\n");
  shape_free(ms,o);

  /*
  LOOP {
    size_t dim = 5;
    ap_dim_t d;
    shape_t *o, *o1;
    ap_texpr0_t* e = shape_texpr_random(pr, mode,intdim,realdim);
    if (mode != expr_lindata || mode != expr_data)
      d = RANDOM_PTRVAR(pr,intdim,realdim);
    else
      d = lrand48()%intdim;
    o = shape_random(pr, SIZE, intdim,realdim);
    if (lrand48()%10>=5) shape_close(pr,o);
    o1 = subst ? shape_substitute_texpr_array(ms,false,o,&d,&e,1,NULL) : shape_assign_texpr_array(ms,false,o,&d,&e,1,NULL);
    FLAG(ms);
    RESULT(shape_check(pr,o1));
    if (shape_is_eq(ms,o,o1)) RESULT('x');
    else RESULT('.');
    shape_free(ms,o); shape_free(ms,o1);
    ap_texpr0_free(e);
  } ENDLOOP;
  */
}

void test_prg_add2(void)
{
  shape_t *top, *l6, *l7, *l8, *l9, *l10, *l11, *l12, *l13, *l14, *l15, *ljoin, *lwiden, *l16;
  ap_tcons0_array_t whileCond, notWhileCond, init0, init;
  ap_dim_t xi, y, z;
  ap_texpr0_t* t7, *t10, *t11, *t12, *t13, *t14;
  size_t i;
  bool fin = false;

  printf("\n**********\n* test_prg_add2 (* expected)\n**********\n");

  pr->max_anon = 0;
  pr->segm_anon = 1;
  /*
   * init program
   */
  xi = DIM_X+1; /* x_i */
  y = DIM_Y; /* y */
  t7 = shape_texpr_x_n_next(DIM_X,0,intdim,realdim); /* x */
  whileCond = ap_tcons0_array_make(1);
  whileCond.p[0] = shape_tcons_diff_x_y (xi, false, NULL_DIM, false, intdim, realdim);
  notWhileCond = ap_tcons0_array_make(1);
  notWhileCond.p[0] = shape_tcons_same_x_y (xi, false, NULL_DIM, false, intdim, realdim);
  t10 = shape_texpr_x_y_v_cst_data(1,xi,0,0,0,0,2,true,intdim,realdim);
  t11 = shape_texpr_x_n_next(xi,1,intdim,realdim);
  t12 = shape_texpr_null(intdim,realdim);
  t13 = shape_texpr_x_n_next(y,0,intdim,realdim); /* y */
  t14 = t12; /* null */
  /*
   * init
   */
  // Solution 1
  //l7 = shape_make(pr,0,intdim,realdim); /* x --> null */
  // Solution 2
  top = shape_top(ms,intdim,realdim);
#ifndef NDEBUG
  printf("\ntop:");
  shape_fdump(stdout,ms,top);
  printf("\n");
#endif
  init0 = ap_tcons0_array_make(1);
  // init with 0=0
  init0.p[0] = ap_tcons0_make(AP_CONS_EQ,ap_texpr0_cst_scalar_int(0),NULL);
  shape_approximate(ms,top,1);
  l6 = shape_meet_tcons_array(ms,false,top,&init0);
  ap_tcons0_clear(init0.p);
  init = ap_tcons0_array_make(1);
  init.p[0] = shape_tcons_reach_x_y (DIM_X, false, NULL_DIM, false, intdim, realdim);
  shape_approximate(ms,top,1);
  l7 = shape_meet_tcons_array(ms,false,top, &init);
  l8 = shape_assign_texpr_array(ms,true,l7,&y,&t14,1,NULL);  /* y := null */
  l9 = shape_assign_texpr_array(ms,true,l8,&xi,&t7,1,NULL);  /* xi := null */
  for (i=0; i <= 4 && !fin; i++)
    {
      l10 = shape_meet_tcons_array(ms,false,l9,&whileCond); /* while (xi != null) */
      l11 = shape_assign_texpr_array(ms,true,l10,&xi,&t10,1,NULL); /* xi->data := xi->data+2 */
      l12 = shape_assign_texpr_array(ms,true,l11,&y,&t11,1,NULL);  /* y := xi->next */
      l13 = shape_assign_texpr_array(ms,true,l12,&xi,&t12,1,NULL);  /* xi := null */
      l14 = shape_assign_texpr_array(ms,true,l13,&xi,&t13,1,NULL);  /* xi := y */
      l15 = shape_assign_texpr_array(ms,true,l14,&y,&t14,1,NULL);  /* y := null */
      /* join l15 with l10 */
      ljoin = shape_join(ms,false,l15,l9);
      lwiden = shape_widening(ms,l9,ljoin);
      if (shape_is_bottom(ms,lwiden))
        {
          printf("test_prg_add2 failed #1\n");
          ERROR_RESULT("bottom reached");
          fin = true;
        }
      if (shape_is_leq(ms,lwiden,l9))
        {
          RESULT('.');
          fin = true;
        }
      shape_free(ms, l9);
      shape_free(ms, l15);
      shape_free(ms, ljoin);
      l9 = lwiden;
      lwiden = NULL;
    }
  if (!shape_is_bottom(ms,l9))
    {
      l16 = shape_meet_tcons_array(ms,true,l9,&notWhileCond); /* (xi == null) */
    }
}

void test_prg_pre_add2(void)
{
  shape_t *top, *l6, *l7, *l8, *l9, *l10, *l11, *l12, *l13, *l14, *l15, *ljoin, *lwiden, *l16;
  ap_tcons0_array_t whileCond, notWhileCond, init0, init;
  ap_dim_t xi, y, z;
  ap_texpr0_t* t7, *t10, *t11, *t12, *t13, *t14;
  size_t i;
  bool fin = false;
  /* for pre */
  shape_t *pre_l10, *pre_l11, *pre_l12;

  printf("\n**********\n* test_prg_pre_add2 (* expected)\n**********\n");

  pr->max_anon = 0;
  pr->segm_anon = 1;
  /*
   * init program
   */
  xi = DIM_X+1; /* x_i */
  y = DIM_Y; /* y */
  t7 = shape_texpr_x_n_next(DIM_X,0,intdim,realdim); /* x */
  whileCond = ap_tcons0_array_make(1);
  whileCond.p[0] = shape_tcons_diff_x_y (xi, false, NULL_DIM, false, intdim, realdim);
  notWhileCond = ap_tcons0_array_make(1);
  notWhileCond.p[0] = shape_tcons_same_x_y (xi, false, NULL_DIM, false, intdim, realdim);
  t10 = shape_texpr_x_y_v_cst_data(1,xi,0,0,0,0,2,true,intdim,realdim);
  t11 = shape_texpr_x_n_next(xi,1,intdim,realdim);
  t12 = shape_texpr_null(intdim,realdim);
  t13 = shape_texpr_x_n_next(y,0,intdim,realdim); /* y */
  t14 = t12; /* null */
  /*
   * init
   */
  top = shape_top(ms,intdim,realdim);
#ifndef NDEBUG
  printf("\ntop:");
  shape_fdump(stdout,ms,top);
  printf("\n");
#endif
  init0 = ap_tcons0_array_make(1);
  // init with 0=0
  init0.p[0] = ap_tcons0_make(AP_CONS_EQ,ap_texpr0_cst_scalar_int(0),NULL);
  shape_approximate(ms,top,1);
  l6 = shape_meet_tcons_array(ms,false,top,&init0);
  ap_tcons0_clear(init0.p);
  init = ap_tcons0_array_make(1);
  init.p[0] = shape_tcons_reach_x_y (DIM_X, false, NULL_DIM, false, intdim, realdim);
  shape_approximate(ms,top,1);
  l7 = shape_meet_tcons_array(ms,false,top, &init);
  l8 = shape_assign_texpr_array(ms,true,l7,&y,&t14,1,NULL);  /* y := null */
  l9 = shape_assign_texpr_array(ms,true,l8,&xi,&t7,1,NULL);  /* xi := null */
  for (i=0; i <= 4 && !fin; i++)
    {
      l10 = shape_meet_tcons_array(ms,false,l9,&whileCond); /* while (xi != null) */
      l11 = shape_assign_texpr_array(ms,true,l10,&xi,&t10,1,NULL); /* xi->data := xi->data+2 */
      // pre_l10 = shape_substitute_texpr_array(ms,false,l11,&xi,&t10,1,NULL);
      l12 = shape_assign_texpr_array(ms,true,l11,&y,&t11,1,NULL);  /* y := xi->next */
      // pre_l11 = shape_substitute_texpr_array(ms,false,l12,&y,&t11,1,NULL);
      l13 = shape_assign_texpr_array(ms,true,l12,&xi,&t12,1,NULL);  /* xi := null */
      // pre_l12 = shape_substitute_texpr_array(ms,false,l13,&xi,&t12,1,NULL);
      l14 = shape_assign_texpr_array(ms,true,l13,&xi,&t13,1,NULL);  /* xi := y */
      l15 = shape_assign_texpr_array(ms,true,l14,&y,&t14,1,NULL);  /* y := null */
      /* join l15 with l10 */
      ljoin = shape_join(ms,false,l15,l9);
      lwiden = shape_widening(ms,l9,ljoin);
      if (shape_is_bottom(ms,lwiden))
        {
          printf("test_prg_add2 failed #1\n");
          ERROR_RESULT("bottom reached");
          fin = true;
        }
      if (shape_is_leq(ms,lwiden,l9))
        {
          RESULT('.');
          fin = true;
        }
      shape_free(ms, l9);
      shape_free(ms, l15);
      shape_free(ms, ljoin);
      l9 = lwiden;
      lwiden = NULL;
    }
  if (!shape_is_bottom(ms,l9))
    {
      l16 = shape_meet_tcons_array(ms,true,l9,&notWhileCond); /* (xi == null) */
    }
}

void test_prg_add2copy(int code)
{
  shape_t *l7, *l9, *l10, *l11, *l12, *l13, *l14, *l15, *l15b, *l16, *l16b,
          *l17, *l17b, *ljoin, *lwiden, *l18;
  ap_tcons0_array_t whileCond, notWhileCond;
  ap_dim_t x, xi, y, yi, z1,z2;
  ap_texpr0_t* t9, *t10, *t11, *t12, *t13, *t14, *t15, *t16, *t16b;
  size_t i;
  bool fin = false;

  printf("\n**********\n* test_prg_add2copy (* expected)\n**********\n");

  pr->max_anon = 0;
  pr->segm_anon = 2;
  /*
   * init program
   */
  // all cases need three more pointer variables
  realdim += 3;
  x = DIM_X;
  xi = DIM_X+1; /* x_i */
  y = DIM_Y;
  yi = DIM_Y+1; /* yi */
  z1 = DIM_Y+2; /* z1 */
  z2 = DIM_Y+3; /* z2 */
  t9 = shape_texpr_x_n_next(DIM_X,0,intdim,realdim); /* x */
  t10 = shape_texpr_x_n_next(DIM_Y,0,intdim,realdim); /* y */
  if (code <=1)
    {
      whileCond = ap_tcons0_array_make(1);
      whileCond.p[0] = shape_tcons_diff_x_y (xi,0,NULL_DIM,0, intdim, realdim);
    }
  else
    {
      whileCond = ap_tcons0_array_make(2);
      whileCond.p[0] = shape_tcons_diff_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
      whileCond.p[1] = shape_tcons_diff_x_y (yi, 0, NULL_DIM, 0, intdim, realdim);
    }
  notWhileCond = ap_tcons0_array_make(1);
  notWhileCond.p[0] = shape_tcons_same_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
  t12 = shape_texpr_x_y_v_cst_data(1,xi,0,0,0,0,2,true,intdim,realdim);
  t13 = shape_texpr_x_n_next(xi,1,intdim,realdim); /* xi->next */
  t14 = shape_texpr_x_n_next(yi,1,intdim,realdim); /* yi->next */
  t15 = shape_texpr_null(intdim,realdim); /* null */
  t16 = shape_texpr_x_n_next(z1,0,intdim,realdim); /* z1 */
  t16b = shape_texpr_x_n_next(z2,0,intdim,realdim); /* z2 */
  /*
   * init
   */
  if (!code)
	  // x-->null and y-->null && length(x)=length(y)
    l7 = shape_make(pr,2,intdim,realdim);
  else
	  // x-->null and y-->null && length(x)>length(y)
    l7 = shape_make(pr,3,intdim,realdim);
  l9 = shape_assign_texpr_array(ms,true,l7,&xi,&t9,1,NULL);
  l10 = shape_assign_texpr_array(ms,true,l9,&yi,&t10,1,NULL);
  for (i=0; i <= 4 && !fin; i++)
    {
      l11 = shape_meet_tcons_array(ms,false,l10,&whileCond); /* while (xi != null) */
      l12 = shape_assign_texpr_array(ms,true,l11,&yi,&t12,1,NULL); /* yi->data := xi->data+2 */
      l13 = shape_assign_texpr_array(ms,true,l12,&z1,&t13,1,NULL); /* z1 := xi->next */
      l14 = shape_assign_texpr_array(ms,true,l13,&z2,&t14,1,NULL); /* z2 := yi->next */
      l15 = shape_assign_texpr_array(ms,true,l14,&xi,&t15,1,NULL); /* xi := null */
      l15b = shape_assign_texpr_array(ms,true,l15,&yi,&t15,1,NULL); /* yi := null */
      l16 = shape_assign_texpr_array(ms,true,l15b,&xi,&t16,1,NULL);  /* xi := z1 */
      l16b = shape_assign_texpr_array(ms,true,l16,&yi,&t16b,1,NULL);  /* yi := z2 */
      l17 = shape_assign_texpr_array(ms,true,l16b,&z1,&t15,1,NULL);  /* z1 := null */
      l17b = shape_assign_texpr_array(ms,true,l17,&z2,&t15,1,NULL);  /* z2 := null */
      /* join l15 with l10 */
      ljoin = shape_join(ms,false,l17b,l10);
      lwiden = shape_widening(ms,l10,ljoin);
      if (shape_is_bottom(ms,lwiden))
        {
          printf("test_prg_add2 failed #1\n");
          ERROR_RESULT("bottom reached");
          fin = true;
        }
      if (shape_is_leq(ms,lwiden,l10))
        {
          RESULT('.');
          fin = true;
        }
      shape_free(ms, l10);
      shape_free(ms, l17b);
      shape_free(ms, ljoin);
      l10 = lwiden;
      lwiden = NULL;
    }
  if (!shape_is_bottom(ms,l10))
    {
      l16 = shape_meet_tcons_array(ms,true,l10,&notWhileCond); /* (xi == null) */
    }
  realdim -= 3;
}

void test_prg_add2new(void)
{
  shape_t *l7, *l11, *l12, *l13, *l14, *l15, *l16, *l17, *l18, *l19, *l20,
          *l21, *l22, *l23, *l24, *l25, *l26, *l27, *l28, *ljoin, *lwiden, *l29;
  ap_tcons0_array_t whileCond, notWhileCond, ifCond, notIfCond;
  ap_dim_t x, xi, y, yi, z, n;
  ap_texpr0_t *t11, *t12, *t13, *t14, *t15, *t17, *t19, *t20, *t22, *t25;
  size_t i;
  bool fin = false;

  printf("\n**********\n* test_prg_add2new (* expected)\n**********\n");

  pr->max_anon = 0;
  pr->segm_anon = 2;
  /*
   * init program
   */
  realdim += 2;
  x = DIM_X;
  xi = DIM_X+1;
  y = DIM_Y;
  yi = DIM_Y+1;
  z = DIM_Y+2;
  t11 = shape_texpr_x_n_next(DIM_X,0,intdim,realdim); /* x */
  whileCond = ap_tcons0_array_make(1);
  whileCond.p[0] = shape_tcons_diff_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
  notWhileCond = ap_tcons0_array_make(1);
  notWhileCond.p[0] = shape_tcons_same_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
  ifCond = ap_tcons0_array_make(1);
  ifCond.p[0] = shape_tcons_same_x_y (yi, 0, NULL_DIM, 0, intdim, realdim);
  notIfCond = ap_tcons0_array_make(1);
  notIfCond.p[0] = shape_tcons_diff_x_y (yi, 0, NULL_DIM, 0, intdim, realdim);
  t13 = shape_texpr_x_alloc(z,intdim,realdim); /* new */
  t14 = shape_texpr_x_y_v_cst_data(1,xi,0,0,0,0,2,true,intdim,realdim);
  t15 = shape_texpr_deref_next(shape_texpr_null(intdim,realdim));
  t17 = shape_texpr_x_n_next(z,0,intdim,realdim); /* z */
  t20 = shape_texpr_x_n_next(z,-1,intdim,realdim); /* z / next */
  t22 = shape_texpr_null(intdim,realdim); /* null */
  t25 = shape_texpr_x_n_next(xi,1,intdim,realdim); /* xi->next */
  /*
   * init
   */
  l7 = shape_make(pr,1,intdim,realdim);
  l11 = shape_assign_texpr_array(ms,true,l7,&xi,&t11,1,NULL);
  for (i=0; i <= 4 && !fin; i++)
    {
      l12 = shape_meet_tcons_array(ms,false,l11,&whileCond); /* while (xi != null) */
      l13 = shape_assign_texpr_array(ms,true,l12,&z,&t13,1,NULL); /* z := new */
      l14 = shape_assign_texpr_array(ms,true,l13,&z,&t14,1,NULL); /* z->data := xi->data+2 */
      l15 = shape_assign_texpr_array(ms,true,l14,&z,&t15,1,NULL); /* z->next := null */
      l16 = shape_meet_tcons_array(ms,false,l15,&ifCond);  /* if (yi == _null) */
      l17 = shape_assign_texpr_array(ms,true,l16,&y,&t17,1,NULL);  /* y := z */
      l18 = shape_meet_tcons_array(ms,true,l15,&notIfCond); /* if !(yi == _null) */
      l19 = shape_assign_texpr_array(ms,true,l18,&yi,&t15,1,NULL);  /* yi->next := null */
      l20 = shape_assign_texpr_array(ms,true,l19,&yi,&t20,1,NULL);  /* yi->next := z */
      l21 = shape_join(ms,true,l17,l20);
      l22 = shape_assign_texpr_array(ms,true,l21,&yi,&t22,1,NULL); /* yi := null */
      l23 = shape_assign_texpr_array(ms,true,l22,&yi,&t17,1,NULL); /* yi := z */
      l24 = shape_assign_texpr_array(ms,true,l23,&z,&t22,1,NULL); /* z := null */
      l25 = shape_assign_texpr_array(ms,true,l24,&z,&t25,1,NULL); /* z := xi->next */
      l26 = shape_assign_texpr_array(ms,true,l25,&xi,&t22,1,NULL); /* xi := null */
      l27 = shape_assign_texpr_array(ms,true,l26,&xi,&t17,1,NULL); /* xi := z */
      l28 = shape_assign_texpr_array(ms,true,l27,&z,&t22,1,NULL); /* z := null */
      /* join l28 with l11 */
      ljoin = shape_join(ms,false,l28,l11);
      lwiden = shape_widening(ms,l11,ljoin);
      if (shape_is_bottom(ms,lwiden))
        {
          printf("test_prg_add2new failed #1\n");
          ERROR_RESULT("bottom reached");
          fin = true;
        }
      if (shape_is_leq(ms,lwiden,l11))
        {
          RESULT('.');
          fin = true;
        }
      shape_free(ms, l11);
      shape_free(ms, l28);
      shape_free(ms, ljoin);
      l11 = lwiden;
      lwiden = NULL;
    }
  if (!shape_is_bottom(ms,l11))
    {
      l29 = shape_meet_tcons_array(ms,true,l11,&notWhileCond); /* (xi == null) */
    }
  realdim -= 2;
}

void test_prg_fstNot0(void)
{
  shape_t *l5, *l7, *l8, *l9, *l10, *l11, *l11b, *l12, *l13, *l14, *l15, *ljoin, *lwiden, *l16;
  ap_tcons0_array_t whileCond, notWhileCond, ifCond, notIfCond;
  ap_dim_t x, xi, y, z;
  ap_texpr0_t* t7, *t10, *t12, *t13, *t14;
  size_t i;
  bool fin = false;

  printf("\n**********\n* test_prg_fstNot0 (* expected)\n**********\n");

  pr->max_anon = 0;
  pr->segm_anon = 1;
  /*
   * init program
   */
  realdim += 1;
  x = DIM_X;
  xi = DIM_X+1; /* x_i */
  y = DIM_Y;
  z = DIM_Y+1;
  t7 = shape_texpr_x_n_next(DIM_X,0,intdim,realdim); /* x */
  whileCond = ap_tcons0_array_make(2);
  whileCond.p[0] = shape_tcons_diff_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
  whileCond.p[1] = shape_tcons_same_x_y (y, 0, NULL_DIM, 0, intdim, realdim);
  ifCond = ap_tcons0_array_make(1);
  ifCond.p[0] = ap_tcons0_make (AP_CONS_DISEQ,
      shape_texpr_x_y_v_cst_data (1,xi,0,0,0,0,0,false,intdim,realdim),NULL);
  notWhileCond = ap_tcons0_array_make(1);
  notWhileCond.p[0] = shape_tcons_diff_x_y (y, 0, NULL_DIM, 0, intdim, realdim);
  notIfCond = ap_tcons0_array_make(1);
  notIfCond.p[0] = ap_tcons0_make (AP_CONS_EQ,
      shape_texpr_x_y_v_cst_data (1,xi,0,0,0,0,0,false,intdim,realdim),NULL);
  t10 = shape_texpr_x_n_next(xi,0,intdim,realdim); /* xi */
  t12 = shape_texpr_x_n_next(xi,1,intdim,realdim); /* xi->next */
  t13 = shape_texpr_null(intdim,realdim); /* null */
  t14 = shape_texpr_x_n_next(z,0,intdim,realdim); /* z */
  /*
   * init
   */
  l5 = shape_make(pr,1,intdim,realdim);
  l7 = shape_assign_texpr_array(ms,true,l5,&xi,&t7,1,NULL); /* xi := x */
  for (i=0; i <= 4 && !fin; i++)
    {
      l8 = shape_meet_tcons_array(ms,false,l7,&whileCond); /* while (xi != null && y==null) */
      l9 = shape_meet_tcons_array(ms,false,l8,&ifCond); /* if (xi->data != 0) */
      l10 = shape_assign_texpr_array(ms,true,l9,&y,&t10,1,NULL); /* y := xi */
      l11 = shape_meet_tcons_array(ms,true,l8,&notIfCond); /* if (xi->data == 0) */
      l11b = shape_join(ms,true,l10,l11);
      l12 = shape_assign_texpr_array(ms,true,l11b,&z,&t12,1,NULL); /* z := xi->next */
      l13 = shape_assign_texpr_array(ms,true,l12,&xi,&t13,1,NULL); /* xi := null */
      l14 = shape_assign_texpr_array(ms,true,l13,&xi,&t14,1,NULL); /* xi := z */
      l15 = shape_assign_texpr_array(ms,true,l14,&z,&t13,1,NULL);  /* z := null */
      /* join l15 with l7 */
      ljoin = shape_join(ms,false,l15,l7);
      lwiden = shape_widening(ms,l7,ljoin);
      if (shape_is_bottom(ms,lwiden))
        {
          printf("test_prg_fstNot0 failed #1\n");
          ERROR_RESULT("bottom reached");
          fin = true;
        }
      if (shape_is_leq(ms,lwiden,l7))
        {
          RESULT('.');
          fin = true;
        }
      shape_free(ms, l7);
      shape_free(ms, l15);
      shape_free(ms, ljoin);
      l7 = lwiden;
      lwiden = NULL;
    }
  if (!shape_is_bottom(ms,l7))
    {
      l16 = shape_meet_tcons_array(ms,true,l7,&notWhileCond); /* (y != null) */
    }
  realdim -= 1;
}

void test_prg_getMax(void)
{
  shape_t *l5, *l7, *l8, *l9, *l10, *l11, *l12b, *l12, *l13, *l14, *l15, *l16, *ljoin, *lwiden, *l17;
  ap_tcons0_array_t whileCond, notWhileCond, ifCond, notIfCond;
  ap_dim_t x, xi, y, m;
  ap_texpr0_t* t7, *t8, *t11, *t13, *t14, *t15;
  size_t i;
  bool fin = false;

  printf("\n**********\n* test_prg_getMax (* expected)\n**********\n");

  pr->max_anon = 0;
  pr->segm_anon = 1;
  /*
   * init program
   */
  intdim +=1;
  x = DIM_X;
  xi = DIM_X+1; /* x_i */
  y = DIM_Y;
  m = intdim-1;
  t7 = shape_texpr_x_y_v_cst_data (1,x,0,0,0,0,0,false,intdim,realdim);
  t8 = shape_texpr_x_n_next(DIM_X,0,intdim,realdim); /* x */
  whileCond = ap_tcons0_array_make(1);
  whileCond.p[0] = shape_tcons_diff_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
  ifCond = ap_tcons0_array_make(1);
  ifCond.p[0] = ap_tcons0_make (AP_CONS_SUPEQ,
      shape_texpr_x_y_v_cst_data (1,xi,0,0,-1,m,-1,false,intdim,realdim),NULL);
  notWhileCond = ap_tcons0_array_make(1);
  notWhileCond.p[0] = shape_tcons_same_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
  notIfCond = ap_tcons0_array_make(1);
  notIfCond.p[0] = ap_tcons0_make (AP_CONS_SUPEQ,
      shape_texpr_x_y_v_cst_data (-1,xi,0,0,1,m,0,false,intdim,realdim),NULL);
  t11 = shape_texpr_x_y_v_cst_data (1,xi,0,0,0,0,0,false,intdim,realdim); /* xi->data */
  t13 = shape_texpr_x_n_next(xi,1,intdim,realdim); /* xi->next */
  t14 = shape_texpr_null(intdim,realdim); /* null */
  t15 = shape_texpr_x_n_next(y,0,intdim,realdim); /* y */
  /*
   * init
   */
  l5 = shape_make(pr,0,intdim,realdim);
  l7 = shape_assign_texpr_array(ms,true,l5,&m,&t7,1,NULL); /* m := x->data */
  l8 = shape_assign_texpr_array(ms,true,l7,&xi,&t8,1,NULL); /* xi := x */
  for (i=0; i <= 4 && !fin; i++)
    {
      l9 = shape_meet_tcons_array(ms,false,l8,&whileCond); /* while (xi != null) */
      l10 = shape_meet_tcons_array(ms,false,l9,&ifCond); /* if (xi->data - m - 1 >= 0) */
      l11 = shape_assign_texpr_array(ms,true,l10,&m,&t11,1,NULL); /* m := xi->data */
      l12 = shape_meet_tcons_array(ms,true,l9,&notIfCond); /* if (m - xi->data >= 0) */
      l12b = shape_join(ms,true,l11,l12);
      l13 = shape_assign_texpr_array(ms,true,l12b,&y,&t13,1,NULL); /* y := xi->next */
      l14 = shape_assign_texpr_array(ms,true,l13,&xi,&t14,1,NULL); /* xi := null */
      l15 = shape_assign_texpr_array(ms,true,l14,&xi,&t15,1,NULL); /* xi := y */
      l16 = shape_assign_texpr_array(ms,true,l15,&y,&t14,1,NULL);  /* y := null */
      /* join l16 with l8 */
      ljoin = shape_join(ms,false,l16,l8);
      lwiden = shape_widening(ms,l8,ljoin);
      if (shape_is_bottom(ms,lwiden))
        {
          printf("test_prg_getMax failed #1\n");
          ERROR_RESULT("bottom reached");
          fin = true;
        }
      if (shape_is_leq(ms,lwiden,l8))
        {
          RESULT('.');
          fin = true;
        }
      shape_free(ms, l8);
      shape_free(ms, l16);
      shape_free(ms, ljoin);
      l8 = lwiden;
      lwiden = NULL;
    }
  if (!shape_is_bottom(ms,l8))
    {
      l17 = shape_meet_tcons_array(ms,true,l8,&notWhileCond); /* (xi == null) */
    }
  intdim -= 1;
}

void test_prg_initMod2(void)
{
  shape_t *l5, *l7, *l8, *l9, *l10, *l11, *l12b, *l12, *l13, *l14,
  *l15, *l16, *l17, *l18, *l19, *l20, *ljoin, *lwiden, *l21;
  ap_tcons0_array_t whileCond, notWhileCond, ifCond, notIfCond;
  ap_dim_t x, xi, y, k;
  ap_texpr0_t* t7, *t8, *t11, *t12, *t14, *t17, *t18, *t19;
  size_t i;
  bool fin = false;

  printf("\n**********\n* test_prg_initMod2 (* expected)\n**********\n");

  pr->max_anon = 1;
  pr->segm_anon = 1;
  /*
   * init program
   */
  x = DIM_X;
  xi = DIM_X+1; /* x_i */
  y = DIM_Y;
  k = 1;
  t7 = ap_texpr0_cst_scalar_int(0); /* value 0 */
  t8 = shape_texpr_x_n_next(DIM_X,0,intdim,realdim); /* x */
  whileCond = ap_tcons0_array_make(1);
  whileCond.p[0] = shape_tcons_diff_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
  ifCond = ap_tcons0_array_make(1);
  ifCond.p[0] = ap_tcons0_make (AP_CONS_SUPEQ,
      shape_texpr_x_y_v_cst_data (0,0,0,0,-1,k,0,false,intdim,realdim),NULL);
  notWhileCond = ap_tcons0_array_make(1);
  notWhileCond.p[0] = shape_tcons_same_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
  notIfCond = ap_tcons0_array_make(1);
  notIfCond.p[0] = ap_tcons0_make (AP_CONS_SUPEQ,
      shape_texpr_x_y_v_cst_data (0,0,0,0,1,k,-1,false,intdim,realdim),NULL);
  t11 = shape_texpr_x_y_v_cst_data (0,0,0,0,0,0,0,true,intdim,realdim); /* 0 / _data */
  t12 = ap_texpr0_cst_scalar_int(1); /* value 1 */
  t14 = shape_texpr_x_y_v_cst_data (0,0,0,0,0,0,1,true,intdim,realdim); /* 1 / _data */
  t17 = shape_texpr_x_n_next(xi,1,intdim,realdim); /* xi->next */
  t18 = shape_texpr_null(intdim,realdim); /* null */
  t19 = shape_texpr_x_n_next(y,0,intdim,realdim); /* y */
  /*
   * init
   */
  l5 = shape_make(pr,1,intdim,realdim);
  l7 = shape_assign_texpr_array(ms,true,l5,&k,&t7,1,NULL); /* k := 0 */
  l8 = shape_assign_texpr_array(ms,true,l7,&xi,&t8,1,NULL); /* xi := x */
  for (i=0; i <= 6 && !fin; i++)
    {
      l9 = shape_meet_tcons_array(ms,false,l8,&whileCond); /* while (xi != null) */
      l10 = shape_meet_tcons_array(ms,false,l9,&ifCond); /* if (k <= 0) */
      l11 = shape_assign_texpr_array(ms,true,l10,&xi,&t11,1,NULL); /* xi->data := 0 */
      l12 = shape_assign_texpr_array(ms,true,l11,&k,&t12,1,NULL); /* k := 1 */
      l13 = shape_meet_tcons_array(ms,true,l9,&notIfCond); /* if (k >= 1) */
      l14 = shape_assign_texpr_array(ms,true,l13,&xi,&t14,1,NULL); /* xi->data := 1 */
      l15 = shape_assign_texpr_array(ms,true,l14,&k,&t7,1,NULL); /* k := 0 */
      l16 = shape_join(ms,true,l12,l15);
      l17 = shape_assign_texpr_array(ms,true,l16,&y,&t17,1,NULL); /* y := xi->next */
      l18 = shape_assign_texpr_array(ms,true,l17,&xi,&t18,1,NULL); /* xi := null */
      l19 = shape_assign_texpr_array(ms,true,l18,&xi,&t19,1,NULL); /* xi := y */
      l20 = shape_assign_texpr_array(ms,true,l19,&y,&t18,1,NULL);  /* y := null */
      /* join l20 with l8 */
      ljoin = shape_join(ms,false,l20,l8);
      lwiden = shape_widening(ms,l8,ljoin);
      if (shape_is_bottom(ms,lwiden))
        {
          printf("test_prg_initMod2 failed #1\n");
          ERROR_RESULT("bottom reached");
          fin = true;
        }
      if (shape_is_leq(ms,lwiden,l8))
        {
          RESULT('.');
          fin = true;
        }
      shape_free(ms, l8);
      shape_free(ms, l20);
      shape_free(ms, ljoin);
      l8 = lwiden;
      lwiden = NULL;
    }
  if (!shape_is_bottom(ms,l8))
    {
      l21 = shape_meet_tcons_array(ms,true,l8,&notWhileCond); /* (xi == null) */
    }
}


void test_prg_delAllGeV(void)
{
  shape_t *l6, *l8, *l9, *l10, *l11, *l12, *l13, *l14, *l15, *l16,
  *l17, *l18, *l19, *l20, *l21, *l22, *l23, *l24, *l25, *l26, *l27, *l28, *l29, *l30,
  *ljoin, *lwiden, *l31;
  ap_tcons0_array_t whileCond, notWhileCond, ifCond1, notIfCond1, ifCond2, notIfCond2;
  ap_dim_t x, xi, y, z, v;
  ap_texpr0_t *t7, *t8, *t11, *t13, *t18, *t19,
              *t21, *t26, *t28;
  size_t i;
  bool fin = false;

  printf("\n**********\n* test_prg_delAllGeV (* expected)\n**********\n");

  pr->max_anon = 0;
  pr->segm_anon = 1;
  /*
   * init program
   */
  realdim+=1; /* for z */
  intdim+=1; /* for v */
  x = DIM_X;
  xi = DIM_X+1; /* x_i */
  y = DIM_Y;
  z = DIM_Y+1;
  v = 3;

  t7 = shape_texpr_null(intdim,realdim); /* null */
  t8 = shape_texpr_x_n_next(x,0,intdim,realdim); /* x */
  whileCond = ap_tcons0_array_make(1);
  whileCond.p[0] = shape_tcons_diff_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
  ifCond1 = ap_tcons0_array_make(1);
  ifCond1.p[0] = ap_tcons0_make (AP_CONS_SUPEQ,
      shape_texpr_x_y_v_cst_data (1,xi, 0,0,-1,v,0,false,intdim,realdim),NULL);
  ifCond2 = ap_tcons0_array_make(1);
  ifCond2.p[0] = shape_tcons_same_x_y (y, 0, NULL_DIM, 0, intdim, realdim);
  notWhileCond = ap_tcons0_array_make(1);
  notWhileCond.p[0] = shape_tcons_same_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
  notIfCond1 = ap_tcons0_array_make(1);
  notIfCond1.p[0] = ap_tcons0_make (AP_CONS_SUPEQ,
      shape_texpr_x_y_v_cst_data (-1,xi,0,0,1,v,-1,false,intdim,realdim),NULL);
  notIfCond2 = ap_tcons0_array_make(1);
  notIfCond2.p[0] = shape_tcons_diff_x_y (y, 0, NULL_DIM, 0, intdim, realdim);
  t11 = shape_texpr_x_n_next(xi,0,intdim,realdim); /* xi */
  t13 = shape_texpr_x_n_next(z,1,intdim,realdim); /* z->next */
  t18 = shape_texpr_deref_next(shape_texpr_null(intdim,realdim)); /* null / next */
  t19 = shape_texpr_x_n_next(xi,-1,intdim,realdim); /* xi / next */
  t21 = shape_texpr_x_free(z,intdim,realdim); /* free(z) */
  t26 = shape_texpr_x_n_next(xi,1,intdim,realdim); /* xi->next */
  t28 = shape_texpr_x_n_next(z,0,intdim,realdim); /* z */

  /*
   * init
   */
  l6 = shape_make(pr,0,intdim,realdim);
  l8 = shape_assign_texpr_array(ms,true,l6,&xi,&t8,1,NULL); /* y := x */
  for (i=0; i <= 6 && !fin; i++)
    {
      l9 = shape_meet_tcons_array(ms,false,l8,&whileCond); /* while (xi != null) */
      l10 = shape_meet_tcons_array(ms,false,l9,&ifCond1); /* if (xi->data >= v) */
      l11 = shape_assign_texpr_array(ms,true,l10,&z,&t11,1,NULL); /* z = xi */
      l12 = shape_assign_texpr_array(ms,true,l11,&xi,&t7,1,NULL); /* xi = null */
      l13 = shape_assign_texpr_array(ms,true,l12,&xi,&t13,1,NULL); /* xi = z->next */
      l14 = shape_meet_tcons_array(ms,false,l13,&ifCond2); /* y == null */
      l15 = shape_assign_texpr_array(ms,true,l14,&x,&t7,1,NULL); /* x = null */
      l16 = shape_assign_texpr_array(ms,true,l15,&x,&t11,1,NULL); /* x = xi */
      l17 = shape_meet_tcons_array(ms,true,l13,&notIfCond2); /* y != null */
      l18 = shape_assign_texpr_array(ms,true,l17,&y,&t18,1,NULL); /* y->next = null */
      l19 = shape_assign_texpr_array(ms,true,l18,&y,&t19,1,NULL); /* y->next = xi */
      l20 = shape_join(ms,true,l16,l19);
      l21 = shape_assign_texpr_array(ms,true,l20,&z,&t21,1,NULL); /* z = free */
      l22 = shape_assign_texpr_array(ms,true,l21,&z,&t7,1,NULL); /* z = null */

      l23 = shape_meet_tcons_array(ms,true,l9,&notIfCond1); /* if (xi->data +1<= v) */
      l24 = shape_assign_texpr_array(ms,true,l23,&y,&t7,1,NULL); /* y = null */
      l25 = shape_assign_texpr_array(ms,true,l24,&y,&t11,1,NULL); /* y = xi */
      l26 = shape_assign_texpr_array(ms,true,l25,&z,&t26,1,NULL); /* z = xi->next */
      l27 = shape_assign_texpr_array(ms,true,l26,&xi,&t7,1,NULL); /* xi = null */
      l28 = shape_assign_texpr_array(ms,true,l27,&xi,&t28,1,NULL); /* xi = z */
      l29 = shape_assign_texpr_array(ms,true,l28,&z,&t7,1,NULL); /* z = null */

      l30 = shape_join(ms,true,l22,l29);
      /* join l30 with l8 */
      ljoin = shape_join(ms,false,l30,l8);
      lwiden = shape_widening(ms,l8,ljoin);
      if (shape_is_bottom(ms,lwiden))
        {
          printf("test_prg_delAllGeV failed #1\n");
          ERROR_RESULT("bottom reached");
          fin = true;
        }
      if (shape_is_leq(ms,lwiden,l8))
        {
          RESULT('.');
          fin = true;
        }
      shape_free(ms, l8);
      shape_free(ms, l30);
      shape_free(ms, ljoin);
      l8 = lwiden;
      lwiden = NULL;
    }
  if (!shape_is_bottom(ms,l8))
    {
      l31 = shape_meet_tcons_array(ms,true,l8,&notWhileCond); /* (xi == null) */
    }
   realdim-=1;
   intdim-=1;
}


void test_prg_copyRev(void)
{
  shape_t *l6, *l8, *l9, *l10, *l11, *l12, *l13, *l14, *l15, *l16,
  *l17, *l18, *l19, *ljoin, *lwiden, *l20;
  ap_tcons0_array_t whileCond, notWhileCond;
  ap_dim_t x, xi, y, z;
  ap_texpr0_t *t7, *t8, *t10, *t11, *t12, *t14, *t16;
  size_t i;
  bool fin = false;

  printf("\n**********\n* test_prg_copyRev (* expected)\n**********\n");

  pr->max_anon = 0;
  pr->segm_anon = 1;
  /*
   * init program
   */
  realdim+=1; /* for z */
  x = DIM_X;
  xi = DIM_X+1; /* x_i */
  y = DIM_Y;
  z = DIM_Y+1;

  t7 = shape_texpr_null(intdim,realdim); /* null */
  t8 = shape_texpr_x_n_next(x,0,intdim,realdim); /* x */
  whileCond = ap_tcons0_array_make(1);
  whileCond.p[0] = shape_tcons_diff_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
  notWhileCond = ap_tcons0_array_make(1);
  notWhileCond.p[0] = shape_tcons_same_x_y (xi, 0, NULL_DIM, 0, intdim, realdim);
  t10 = shape_texpr_x_alloc(z,intdim,realdim); /* alloc(1,z) */
  t11 = shape_texpr_x_y_v_cst_data (1,xi,0,0,0,0,0,true,intdim,realdim); /* xi * data / _data */
  t12 = shape_texpr_x_n_next(y,-1,intdim,realdim); /* y / next */
  t14 = shape_texpr_x_n_next(z,0,intdim,realdim); /* z */
  t16 = shape_texpr_x_n_next(xi,1,intdim,realdim); /* xi->next */

  /*
   * init
   */
  l6 = shape_make(pr,0,intdim,realdim);
  l8 = shape_assign_texpr_array(ms,true,l6,&xi,&t8,1,NULL); /* xi := x */
  for (i=0; i <= 6 && !fin; i++)
    {
      l9 = shape_meet_tcons_array(ms,false,l8,&whileCond); /* while (xi != null) */
      l10 = shape_assign_texpr_array(ms,true,l9,&z,&t10,1,NULL); /* z = new() */
      l11 = shape_assign_texpr_array(ms,true,l10,&z,&t11,1,NULL); /* z->data = xi->data */
      l12 = shape_assign_texpr_array(ms,true,l11,&z,&t12,1,NULL);  /* z->next = y */
      l13 = shape_assign_texpr_array(ms,true,l12,&y,&t7,1,NULL); /* y = null */
      l14 = shape_assign_texpr_array(ms,true,l13,&y,&t14,1,NULL); /* y = z */
      l15 = shape_assign_texpr_array(ms,true,l14,&z,&t7,1,NULL); /* z = null */
      l16 = shape_assign_texpr_array(ms,true,l15,&z,&t16,1,NULL); /* z = xi->next */
      l17 = shape_assign_texpr_array(ms,true,l16,&xi,&t7,1,NULL); /* xi = null */
      l18 = shape_assign_texpr_array(ms,true,l17,&xi,&t14,1,NULL); /* xi = z */
      l19 = shape_assign_texpr_array(ms,true,l18,&z,&t7,1,NULL); /* z = null */

      /* join l19 with l8 */
      ljoin = shape_join(ms,false,l19,l8);
      lwiden = shape_widening(ms,l8,ljoin);
      if (shape_is_bottom(ms,lwiden))
        {
          printf("test_prg_copyRev failed #1\n");
          ERROR_RESULT("bottom reached");
          fin = true;
        }
      if (shape_is_leq(ms,lwiden,l8))
        {
          RESULT('.');
          fin = true;
        }
      shape_free(ms, l8);
      shape_free(ms, l19);
      shape_free(ms, ljoin);
      l8 = lwiden;
      lwiden = NULL;
    }
  if (!shape_is_bottom(ms,l8))
    {
      l20 = shape_meet_tcons_array(ms,true,l8,&notWhileCond); /* (xi == null) */
    }
   realdim-=1;
}

/* ********************************* */
/*           hash                    */
/* ********************************* */

void test_hash(void)
{
  printf("\nhash\n");
  LOOP {
    size_t dim = 5;
    shape_t *oa,*ob;
    int ra,rb;
    oa = shape_random(pr,SIZE,intdim,realdim);
    ra = shape_hash(ms,oa);
    ob = shape_copy(ms,oa);
    rb = shape_hash(ms,ob);
    RESULT('.');
    if (ra!=rb) ERROR_RESULT("different hash");
    shape_free(ms,oa); shape_free(ms,ob);
  } ENDLOOP;
}

/* ********************************* */
/*           main                    */
/* ********************************* */

void tests(int algo)
{
  int i;
  for (i=0;i<AP_FUNID_SIZE;i++) {
    ms->option.funopt[i].algorithm = algo;
  }
  printf("\nstarting tests with algo=%i\n",algo);
  /* tests */

  /*
  test_misc();

  test_add_lincons(expr_reach);
  test_add_lincons(expr_reachl);
  test_add_lincons(expr_lindata);
  test_add_tcons(expr_next);
  test_add_tcons(expr_data);

  test_assign(0,expr_ptr);
  test_assign(0,expr_next);
  test_assign(0,expr_lindata);
  */

  /*
  test_prg_add2(); 
  test_prg_add2copy(0);
  test_prg_add2copy(1);
  test_prg_add2copy(2);
  test_prg_fstNot0();
  test_prg_getMax();
  test_prg_initMod2();
  test_prg_delAllGeV();
  test_prg_copyRev();
  */
  test_prg_pre_add2(); 

}

int main(int argc, const char** argv)
{
  long int seed;
  int i;

  seed = time(NULL);
  for (i=1;i<argc;i++) {
    if (argv[i][0]=='+') N = atol(argv[i]+1);
    else seed = atol(argv[i]);
  }
  printf("seed = %ld, N= %i\n", seed, N);

  assert(N<MAXN);

  /* init */
  srand48(seed);
  ms = shape_manager_alloc();
  if (!ms) return 1;
  for (i=0;i<AP_EXC_SIZE;i++){
    ms->option.abort_if_exception[i] = true;
  }
  pr = shape_init_from_manager(ms,0,0);
  info();

  tests(0);

  /* quit */
  if (pr->error_) printf("\n%i error(s)!\n",pr->error_);
  else printf("\nall tests passed\n");
  ap_manager_free(ms);
  return 0;
}
