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

#include "ushape.h"
#include "ushape_fun.h"
#include "ushape_internal.h"
#include "apron2shape.h"

ap_manager_t* ms; /* shape manager */

shape_internal_t* pr;

size_t realdim = 6+3; /* vars x,y,z */
size_t intdim = 3;  /* vars _l, _k, d */
size_t SIZE = 4;     /* nb of nodes */

char* names_dim[12] = { "_l", "_k", "d", "_data", "_free", "_len",
                        "_next", "_new", "_null", "x", "xi", "y" };
#define DIM_NULL 8
#define DIM_X 9
#define DIM_Y 11

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


void print_ushape(const char* msg, ushape_t* a)
{
  fprintf(stdout,"%s = ",msg);
  ushape_fprint(stdout,ms,a,names_dim);
  //ushape_fdump(stdout,ms,a);
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
  printf("ushapes:   %s (%s)\n",ms->library,ms->version);
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
  ushape_t* bot = ushape_bottom(ms,intdim,realdim);
  ushape_t* top = ushape_top(ms,intdim,realdim);
  ap_dimension_t d1 = ushape_dimension(ms,bot);
  ap_dimension_t d2 = ushape_dimension(ms,top);
  printf("\nperforming various tests\n");
  if (ushape_check(pr,bot) != '.')	  printf("ushape_bottom failed #0\n");
  if (ushape_check(pr,top) != '.')	  printf("ushape_top failed #0\n");
  if (d1.intdim!=intdim ||
      d1.realdim!=realdim)         printf("ushape_dimension failed #1\n");
  if (d2.intdim!=intdim ||
      d2.realdim!=realdim)         printf("ushape_dimension failed #2\n");
  if (!ushape_is_bottom(ms,bot))  printf("ushape_is_bottom failed #3\n");
  if (ushape_is_bottom(ms,top))   printf("ushape_is_bottom failed #4\n");
  if (ushape_is_top(ms,bot))      printf("ushape_is_top failed #5\n");
  if (!ushape_is_top(ms,top))     printf("ushape_is_top failed #6\n");
  if (!ushape_is_leq(ms,bot,top)) printf("ushape_is_leq failed #7\n");
  if (ushape_is_leq(ms,top,bot))  printf("ushape_is_leq failed #8\n");
  if (!ushape_is_eq(ms,bot,bot))  printf("ushape_is_eq failed #9\n");
  if (!ushape_is_eq(ms,top,top))  printf("ushape_is_eq failed #10\n");
  if (ushape_is_eq(ms,bot,top))   printf("ushape_is_eq failed #11\n");
  if (ushape_is_dimension_unconstrained(ms,bot,pr->cst_ptr_null+1))
    printf("ushape_is_dimension_unconstrained #12\n");
  if (!ushape_is_dimension_unconstrained(ms,top,pr->cst_ptr_null+1))
    printf("ushape_is_dimension_unconstrained #13\n");
  for (i=0;i<N;i++) {
    ushape_t* o = ushape_random(pr, SIZE, intdim,realdim);
    ushape_t* c = ushape_copy(ms,o);
    ushape_t* l = ushape_closure(ms,false,o);
    ap_dimension_t d = ushape_dimension(ms,o);
    if (d.intdim!=intdim || d.realdim!=realdim) printf("ushape_dimension failed #14\n");
    if (!ushape_is_leq(ms,bot,o))  printf("ushape_is_leq failed #15\n");
    if (!ushape_is_leq(ms,o,top))  printf("ushape_is_leq failed #16\n");
    if (!ushape_is_eq(ms,o,c))     printf("ushape_is_eq failed #17\n");
    ushape_size(ms,o);
    ushape_close(pr,o);
    // not implemented
    //ushape_minimize(ms,o);
    //ushape_canonicalize(ms,o);
    //ushape_approximate(ms,o,0);
    //ushape_is_minimal(ms,o);
    //ushape_is_canonical(ms,o);
    ushape_free(ms,o); ushape_free(ms,c); ushape_free(ms,l);
  }
  ushape_free(ms,bot); ushape_free(ms,top);
}



/* ********************************* */
/*           closure                 */
/* ********************************* */


void test_closure(void)
{
  // printf("\nclosure %s\n",num_incomplete?"":"(c,o expected)");
  LOOP {
    ushape_t* o = ushape_random(pr, SIZE, intdim, realdim);
    ushape_close(pr,o);
    RESULT(ushape_check(pr,o));
    ushape_free(ms,o);
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
    ushape_t *o, *oo;
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
    o  = ushape_top(ms,intdim,realdim);
    oo  = ushape_meet_lincons_array(ms,true,o,&t); FLAG(ms);
    RESULT(ushape_check(pr,o));
    if (ushape_is_eq(ms,o,oo)) RESULT('x');
    if (ushape_is_leq(ms,oo,o)) ERROR_RESULT("best flag");
    ushape_free(ms,o); ushape_free(ms,oo);
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
    ushape_t *o1, *o2, *o;
    o1 = ushape_random(pr, SIZE,intdim,realdim);
    o2 = ushape_random(pr, SIZE,intdim,realdim);
    o  = ushape_meet(ms,false,o1,o2); FLAG(ms);
    RESULT(ushape_check(pr,o)); ushape_check(pr,o1); ushape_check(pr,o2);
    if (!ushape_is_leq(ms,o,o1) || !ushape_is_leq(ms,o,o2)) {
      ERROR_RESULT("not lower bound");
      print_ushape("o1",o1); print_ushape("o2",o2); print_ushape("o",o);
    }
    ushape_free(ms,o); ushape_free(ms,o1); ushape_free(ms,o2);
  } ENDLOOP;
  printf("\nmeet top %s\n","(* expected)");
  LOOP {
    int dim = 8;
    ushape_t *o1, *o2, *o, *oo;
    o1 = ushape_random(pr, SIZE, intdim,realdim);
    o2 = ushape_top(ms,intdim,realdim);
    o  = ushape_meet(ms,false,o1,o2);
    oo = ushape_meet(ms,false,o,o1);
    ushape_check(pr,o1); ushape_check(pr,o2); ushape_check(pr,o);
    if (!ushape_is_eq(ms,o,o1)) {
      ERROR_RESULT("not eq #1");
      print_ushape("o1",o1); print_ushape("o",o);
    }
    else if (!ushape_is_eq(ms,o,oo)) {
      ERROR_RESULT("not eq #2");
      print_ushape("o1",o1); print_ushape("o",o); print_ushape("oo",oo);
    }
    else RESULT('*');
    ushape_free(ms,o); ushape_free(ms,o1); ushape_free(ms,o2); ushape_free(ms,oo);
  } ENDLOOP;
  printf("\nmeet bot %s\n","(* expected)");
  LOOP {
    int dim = 8;
    ushape_t *o1, *o2, *o;
    o1 = ushape_random(pr, SIZE,intdim,realdim);
    o2 = ushape_bottom(ms,intdim,realdim);
    o  = ushape_meet(ms,false,o1,o2);
    ushape_check(pr,o1); ushape_check(pr,o2); ushape_check(pr,o);
    if (!ushape_is_bottom(ms,o)) {
      ERROR_RESULT("not bottom");
      print_ushape("o1",o1); print_ushape("o",o);
    }
    else RESULT('*');
    ushape_free(ms,o); ushape_free(ms,o1); ushape_free(ms,o2);
  } ENDLOOP;
}

#define NB_MEET 5
void test_meet_array(void)
{
  printf("\nmeet array %s\n","(* expected)");
  LOOP {
    int i, dim = 6;
    ushape_t* o[NB_MEET], *oo;
    for (i=0;i<NB_MEET;i++)
      o[i] = ushape_random(pr,SIZE,intdim,realdim);
    oo = ushape_meet_array(ms,o,NB_MEET); FLAG(ms);
    RESULT(ushape_check(pr,oo));
    for (i=0;i<NB_MEET;i++)
      if (!ushape_is_leq(ms,oo,o[i])) ERROR_RESULT("not lower bound");
    for (i=0;i<NB_MEET;i++) ushape_free(ms,o[i]);
    ushape_free(ms,oo);
  } ENDLOOP;
}

void test_add_pcons_once(exprmode_t mode, ushape_t* o, pcons0_array_t* array)
{
    size_t i, o2size;
    ushape_t* o1;
    ushape_array_t *o2;
    printf("\n*** meet with pcons:\n");
#ifndef NDEBUG
    ushape_fdump(stdout,pr->man,o);
    shape_pcons_array_fdump(stdout,array,o->datadim, o->ptrdim);
#endif
    /* apply domain function (underapproximation) */
    printf("\n*** (underapproximation):\n");
    o1 = ushape_copy_mem(pr,o);
    o2 = ushape_meet_pcons_array(pr,true,o1,array); FLAG(ms);
    ushape_check(pr,o);
    if (o2) {
    for (i=0; i < o2->size; i++) {
      RESULT(ushape_check(pr,o2->p[i]));
      if (!ushape_is_leq(ms,o2->p[i],o)) {
        ERROR_RESULT("not included in");
        print_ushape("o",o); print_ushape("o2",o2->p[i]);
          }
    }
    }
    ushape_array_clear(pr,o2);
    free(o2);

    /* apply internal function (exact) */
    printf("\n*** (exact):\n");
    o2 =  ushape_meet_pcons_array(pr, false, o, array);
    if (o2) {
    for (i=0; i < o2->size; i++) {
      RESULT(ushape_check(pr,o2->p[i]));
      if (!ushape_is_leq(ms,o2->p[i],o)) {
        ERROR_RESULT("not included in");
        print_ushape("o",o); print_ushape("o2",o2->p[i]);
          }
    }
    ushape_array_clear(pr,o2);
    }
}

void
test_add_lincons(exprmode_t mode)
{
  size_t i, nb = 4;
  ap_lincons0_array_t ar;
  pcons0_array_t* arr;
  ushape_t* o;
  printf("\nadd %slincons %s\n",exprname[mode],"(* expected)");
  /* predefined graphs and constraints */
  /* case (1): top /\ x --> null /\ x == y */
  o = ushape_top(ms,intdim,realdim);
  ar = ap_lincons0_array_make(2);
  ar.p[0] = shape_lincons_x_y_l(AP_CONS_SUPEQ,1,DIM_X,-1,DIM_NULL,0,0,intdim,realdim);
  ar.p[1] = shape_lincons_x_y_l(AP_CONS_EQ,1,DIM_X,-1,DIM_Y,0,0,intdim,realdim);
  arr = shape_pcons_array_of_lincons_array (pr, &ar, intdim, realdim);
  test_add_pcons_once(mode,o,arr);
  shape_pcons0_array_clear(arr);
  ushape_free(pr->man, o);
  /* case (2): x-->null, x==y  /\ y != null */
  o = ushape_make(pr,5, intdim, realdim);
  ar = ap_lincons0_array_make(1);
  ar.p[0] = shape_lincons_x_y(AP_CONS_DISEQ,DIM_Y,DIM_NULL,intdim,realdim);
  arr = shape_pcons_array_of_lincons_array (pr, &ar, intdim, realdim);
  test_add_pcons_once(mode,o,arr);
  ushape_free(pr->man, o);
  /* case (3): x-->y-->null /\ y != null */
  o = ushape_make(pr,6, intdim, realdim);
  test_add_pcons_once(mode,o,arr);
  ushape_free(pr->man, o);
  shape_pcons0_array_clear(arr);

  /* random graphs and constraints */
  LOOP {
    ar = ap_lincons0_array_make(nb);
    o = ushape_random(pr, SIZE,intdim,realdim);
    if (lrand48()%10>=8) ushape_close(pr,o);
    for (i=0;i<nb;i++) {
      ar.p[i] = ap_lincons0_make((lrand48()%100>=80)?AP_CONS_EQ:
                 (lrand48()%100>=80)?AP_CONS_SUP:
                 AP_CONS_SUPEQ,
                 shape_linexpr_random(pr,mode,intdim,realdim),
                 NULL);
    }
    arr = shape_pcons_array_of_lincons_array (pr, &ar, intdim, realdim);
    test_add_pcons_once(mode,o,arr);
    ushape_free(ms,o);
    shape_pcons0_array_clear(arr);
   } ENDLOOP;
}


void
test_add_tcons(exprmode_t mode)
{
  size_t i, nb = 4;
  ap_tcons0_array_t ar;
  pcons0_array_t* arr;
  ushape_t* o;
  printf("\nadd %s tcons %s\n",exprname[mode],"(* expected)");
  /* predefined graphs and constraints */
  /* case (1): x-->null, x==y /\ x->next != null */
  o = ushape_make(pr,5, intdim, realdim);
  ar = ap_tcons0_array_make(1);
  ar.p[0] = shape_tcons_x_y(AP_CONS_DISEQ,DIM_X,true,DIM_NULL,false,intdim,realdim);
  arr = shape_pcons_array_of_tcons_array (pr, &ar, intdim, realdim);
  test_add_pcons_once(mode,o,arr);
  ushape_free(pr->man, o);
  shape_pcons0_array_clear(arr);

  /* random graphs and constraints */
  LOOP {
    ar = ap_tcons0_array_make(nb);
    o = ushape_random(pr, SIZE,intdim,realdim);
    if (lrand48()%10>=8) ushape_close(pr,o);
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
    test_add_pcons_once(mode,o,arr);
    ushape_free(ms,o);
    shape_pcons0_array_clear(arr);
   } ENDLOOP;
}

/* ********************************* */
/*     assignments / substitutions   */
/* ********************************* */

void test_assign(int subst, exprmode_t mode)
{
  printf("\n%s %slinexpr %s\n", subst ? "subst" : "assign", exprname[mode],
	 "(* expected)");
  LOOP {
    size_t dim = 6;
    ap_dim_t d;
    ushape_t *o, *o1;
    ap_linexpr0_t* l = shape_linexpr_random(pr, mode,intdim,realdim);
    if (mode != expr_lindata || mode != expr_data)
      d = RANDOM_PTRDIM(pr,intdim,realdim);
    else
      d = lrand48()%intdim;
    o = ushape_random(pr, SIZE, intdim,realdim);
    if (lrand48()%10>=5) ushape_close(pr,o);
    o1 = subst ? ushape_substitute_linexpr_array(ms,false,o,&d,&l,1,NULL) :
                 ushape_assign_linexpr_array(ms,false,o,&d,&l,1,NULL);
    FLAG(ms);
    RESULT(ushape_check(pr,o1));
    if (ushape_is_eq(ms,o,o1)) ERROR_RESULT("best flag");
    ushape_free(ms,o); ushape_free(ms,o1);
    ap_linexpr0_free(l);
  } ENDLOOP;
}

void test_assign_texpr(int subst, exprmode_t mode)
{
  printf("\n%s texpr\n", subst ? "subst" : "assign");
  LOOP {
    size_t dim = 5;
    ap_dim_t d;
    ushape_t *o, *o1;
    ap_texpr0_t* e = shape_texpr_random(pr, mode,intdim,realdim);
    if (mode != expr_lindata || mode != expr_data)
      d = RANDOM_PTRDIM(pr,intdim,realdim);
    else
      d = lrand48()%intdim;
    o = ushape_random(pr, SIZE, intdim,realdim);
    if (lrand48()%10>=5) ushape_close(pr,o);
    o1 = subst ? ushape_substitute_texpr_array(ms,false,o,&d,&e,1,NULL) : ushape_assign_texpr_array(ms,false,o,&d,&e,1,NULL);
    FLAG(ms);
    RESULT(ushape_check(pr,o1));
    if (ushape_is_eq(ms,o,o1)) RESULT('x');
    else RESULT('.');
    ushape_free(ms,o); ushape_free(ms,o1);
    ap_texpr0_free(e);
  } ENDLOOP;
}

/* ********************************* */
/*           hash                    */
/* ********************************* */

void test_hash(void)
{
  printf("\nhash\n");
  LOOP {
    size_t dim = 5;
    ushape_t *oa,*ob;
    int ra,rb;
    oa = ushape_random(pr,SIZE,intdim,realdim);
    ra = ushape_hash(ms,oa);
    ob = ushape_copy(ms,oa);
    rb = ushape_hash(ms,ob);
    RESULT('.');
    if (ra!=rb) ERROR_RESULT("different hash");
    ushape_free(ms,oa); ushape_free(ms,ob);
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
  test_misc();
  /*
  test_serialize();
  test_closure();
  test_incremental_closure();
  test_lincons_conversion(expr_ushape);
  test_lincons_conversion(expr_lin);
  test_lincons_conversion(expr_interv);
  test_meet();
  test_meet_array();
  test_join();
  test_join_array();
  */
  test_add_lincons(expr_reach);
  test_add_lincons(expr_reachl);
  test_add_lincons(expr_lindata);
  test_add_tcons(expr_next);
  test_add_tcons(expr_data);
  /*
  test_sat_lincons(expr_reach);
  test_sat_lincons(expr_reachl);
  test_sat_lincons(expr_lindata);
    test_widening();
    test_widening_thrs();
    test_narrowing();
  */
  test_assign(0,expr_ptr);
  test_assign(0,expr_next);
  test_assign(0,expr_lindata);
  test_assign(0,expr_deref);
  test_assign(0,expr_data);
  test_assign(0,expr_deref_data);
  /*
    test_hash();
  */
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
  ms = ushape_manager_alloc();
  if (!ms) return 1;
  for (i=0;i<AP_EXC_SIZE;i++){
    ms->option.abort_if_exception[i] = true;
  }
  pr = ushape_init_from_manager(ms,0,0);
  info();

  tests(0);

  /* quit */
  if (pr->error_) printf("\n%i error(s)!\n",pr->error_);
  else printf("\nall tests passed\n");
  ap_manager_free(ms);
  return 0;
}
