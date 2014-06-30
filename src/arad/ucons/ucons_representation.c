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
 * Functions related to the representation of universal constraints.
 */

#include <stdio.h>
#include <stdlib.h>
#include "sh_macros.h"
#include "ucons.h"
#include "ucons_internal.h"
#include "ap_generic.h"
#include "apron2shape.h"
#include "ap_abstract0.h"
#include "box.h"
#include "oct.h"
#include "pk.h"
#include "ap_ppl.h"

void
initialize_PI (ucons_internal_t *pr)
{
  pattern_info_t* PI;
  size_t PI_size;

  checked_malloc (PI, pattern_info_t, 17, sizeof (pattern_info_t), return;);

  PI_size = 17;
  /* \forall y. */
  pattern_info_t p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;
  /*forall y \in n*/
  p1.kind = pattern_1;
  p1.e_seg = 0;
  p1.u_seg = 1;
  //p1.uvar=NULL;
  checked_malloc (p1.uvar, pattern_var_t, 1, sizeof (pattern_var_t), return;);
  p1.uvar[0].size = 1;
  p1.uvar[0].reach = NULL;
  p1.lcons = -1; // no constraints
  p1.nr_y = 1;

  /*\forall y1 in ?_1 \forall y2 in ?_2 y1=y2 */
  p2.kind = pattern_2_1;
  p2.e_seg = 0;
  p2.u_seg = 2;
  //p2.uvar=NULL;
  checked_malloc (p2.uvar, pattern_var_t, 2, sizeof (pattern_var_t), return;);
  p2.uvar[0].size = 1;
  p2.uvar[0].reach = NULL;
  p2.uvar[1].size = 1;
  p2.uvar[1].reach = NULL;
  p2.lcons = -1; // no constraints y1=y2
  p2.nr_y = 2;

  p3.kind = pattern_3_1;
  p3.e_seg = 0;
  p3.u_seg = 3;
  //p3.uvar=NULL;
  checked_malloc (p3.uvar, pattern_var_t, 3, sizeof (pattern_var_t), return;);
  p3.uvar[0].size = 1;
  p3.uvar[0].reach = NULL;
  p3.uvar[1].size = 1;
  p3.uvar[1].reach = NULL;
  p3.uvar[2].size = 1;
  p3.uvar[2].reach = NULL;
  p3.lcons = -1; // no constraints
  p3.nr_y = 3;

  /*\forall y1, y2. y1 <= y2  */
  p4.kind = pattern_1_2;
  p4.e_seg = 0;
  p4.u_seg = 1;
  checked_malloc (p4.uvar, pattern_var_t, 1, sizeof (pattern_var_t), return;);
  p4.uvar[0].size = 2;
  checked_malloc (p4.uvar[0].reach, bool, 1, sizeof (bool), return;);
  p4.uvar[0].reach[0] = true; //y1<y2
  p4.lcons = -1; // no constraints
  p4.nr_y = 2;


  p5.kind = pattern_succ_1_3;
  p5.e_seg = 0;
  p5.u_seg = 1;
  checked_malloc (p5.uvar, pattern_var_t, 1, sizeof (pattern_var_t), return;);
  p5.uvar[0].size = 3;
  checked_malloc (p5.uvar[0].reach, bool, 2, sizeof (bool), return;);
  p5.uvar[0].reach[0] = false; //y2 = y1 + 1;
  p5.uvar[0].reach[1] = false; //y3 = y2 + 1;
  p5.lcons = -1; // no constraints
  p5.nr_y = 3;

  /* \forall y\in ?. y=1*/
  p6.kind = pattern_1_l1;
  p6.e_seg = 0;
  p6.u_seg = 1;
  checked_malloc (p6.uvar, pattern_var_t, 1, sizeof (pattern_var_t), return;);
  p6.uvar[0].size = 1;
  p6.uvar[0].reach = NULL;
  p6.lcons = -1; // no constraints
  p6.nr_y = 1;

  /* \forall y\in ?. y=l[?]-1*/
  p7.kind = pattern_1_lx_1;
  p7.e_seg = 0;
  p7.u_seg = 1;
  //p1.uvar=NULL;
  checked_malloc (p7.uvar, pattern_var_t, 1, sizeof (pattern_var_t), return;);
  p7.uvar[0].size = 1;
  p7.uvar[0].reach = NULL;
  p7.lcons = -1; //  y=l[x]-1
  p7.nr_y = 1;

  //last added PI[6]

  //same name pattern_1_lx different nb of e_seg
  /* \forall y. y=l[?]*/
  p8.kind = pattern_1_lx;
  p8.e_seg = 1;
  p8.u_seg = 1;
  //p1.uvar=NULL;
  checked_malloc (p8.uvar, pattern_var_t, 1, sizeof (pattern_var_t), return;);
  p8.uvar[0].size = 1;
  p8.uvar[0].reach = NULL;
  p8.lcons = 0; // the first constaint y=l[x]
  p8.nr_y = 1;

  /* \forall y. y=l[?]+l[??]*/
  p9.kind = pattern_1_lx;
  p9.e_seg = 2;
  p9.u_seg = 1;
  //p1.uvar=NULL;
  checked_malloc (p9.uvar, pattern_var_t, 1, sizeof (pattern_var_t), return;);
  p9.uvar[0].size = 1;
  p9.uvar[0].reach = NULL;
  p9.lcons = 0; // the first constaint y=l[x]
  p9.nr_y = 1;

  /* \forall y. y=l[?]+l[??]+l[???]*/
  p10.kind = pattern_1_lx;
  p10.e_seg = 3;
  p10.u_seg = 1;
  //p1.uvar=NULL;
  checked_malloc (p10.uvar, pattern_var_t, 1, sizeof (pattern_var_t), return;);
  p10.uvar[0].size = 1;
  p10.uvar[0].reach = NULL;
  p10.lcons = 0; // the first constraint y=l[x]
  p10.nr_y = 1;


  //same name pattern_2_1_lx different nb of e_seg
  /* \forall y1\in ? y2\in ?? . y1=(l[???])+y2*/
  p11.kind = pattern_2_1_lx;
  p11.e_seg = 1;
  p11.u_seg = 2;
  //p1.uvar=NULL;
  checked_malloc (p11.uvar, pattern_var_t, 2, sizeof (pattern_var_t), return;);
  p11.uvar[0].size = 1;
  p11.uvar[0].reach = NULL;
  p11.uvar[1].size = 1;
  p11.uvar[1].reach = NULL;
  p11.lcons = 1; // the second constaint y1=l[x]+y2
  p11.nr_y = 2;

  /* \forall y1\in ? y2\in ?? . y1=(l[???]+l[???])+y2*/
  p12.kind = pattern_2_1_lx;
  p12.e_seg = 2;
  p12.u_seg = 2;
  //p1.uvar=NULL;
  checked_malloc (p12.uvar, pattern_var_t, 2, sizeof (pattern_var_t), return;);
  p12.uvar[0].size = 1;
  p12.uvar[0].reach = NULL;
  p12.uvar[1].size = 1;
  p12.uvar[1].reach = NULL;
  p12.lcons = 1; // the second constaint y1=l[x]+y2
  p12.nr_y = 2;

  /* \forall y1\in ? y2\in ?? . y1=l[!!!]+ l[!!] + l[!] +y2*/
  p13.kind = pattern_2_1_lx;
  p13.e_seg = 3;
  p13.u_seg = 2;
  //p1.uvar=NULL;
  checked_malloc (p13.uvar, pattern_var_t, 2, sizeof (pattern_var_t), return;);
  p13.uvar[0].size = 1;
  p13.uvar[0].reach = NULL;
  p13.uvar[1].size = 1;
  p13.uvar[1].reach = NULL;
  p13.lcons = 1; // the second constaint y1=l[x]+y2
  p13.nr_y = 2;


  //same name pattern_2_1_mlx different nb of e_seg
  /* \forall y1\in ? y2\in ?? . y1+ (l[???])= y2*/
  p14.kind = pattern_2_1_mlx;
  p14.e_seg = 1;
  p14.u_seg = 2;
  //p1.uvar=NULL;
  checked_malloc (p14.uvar, pattern_var_t, 2, sizeof (pattern_var_t), return;);
  p14.uvar[0].size = 1;
  p14.uvar[0].reach = NULL;
  p14.uvar[1].size = 1;
  p14.uvar[1].reach = NULL;
  p14.lcons = 1; // the second constaint y1=l[x]+y2
  p14.nr_y = 2;

  /* \forall y1\in ? y2\in ?? . y1+(l[???]+l[???])=y2*/
  p15.kind = pattern_2_1_mlx;
  p15.e_seg = 2;
  p15.u_seg = 2;
  //p1.uvar=NULL;
  checked_malloc (p15.uvar, pattern_var_t, 2, sizeof (pattern_var_t), return;);
  p15.uvar[0].size = 1;
  p15.uvar[0].reach = NULL;
  p15.uvar[1].size = 1;
  p15.uvar[1].reach = NULL;
  p15.lcons = 1; // the second constaint y1=l[x]+y2
  p15.nr_y = 2;

  /* \forall y1\in ? y2\in ?? . y1 + l[!!!]+ l[!!] + l[!]  =y2*/
  p17.kind = pattern_2_1_mlx;
  p17.e_seg = 3;
  p17.u_seg = 2;
  //p1.uvar=NULL;
  checked_malloc (p17.uvar, pattern_var_t, 2, sizeof (pattern_var_t), return;);
  p17.uvar[0].size = 1;
  p17.uvar[0].reach = NULL;
  p17.uvar[1].size = 1;
  p17.uvar[1].reach = NULL;
  p17.lcons = 1; // the second constaint y1=l[x]+y2
  p17.nr_y = 2;

  //	/* \forall y1\in ? y2\in ?? . y1 + l[!!!]+ l[!!] + l[!]  =y2*/
  //	p17.kind=pattern_2_1_mlx;
  //	p17.e_seg = 3;
  //	p17.u_seg = 2;
  //	//p1.uvar=NULL;
  //	checked_malloc(p17.uvar,pattern_var_t,2,sizeof(pattern_var_t),return ;);
  //	p17.uvar[0].size=1;
  //	p17.uvar[0].reach=NULL;
  //	p17.uvar[1].size=1;
  //	p17.uvar[1].reach=NULL;
  //	p17.lcons= 1;// the second constaint y1=l[x]+y2
  //	p17.nr_y = 2;

  /* \forall y1, y2 \in ?. y1 + 1 = y2 */
  p16.kind = pattern_succ_1_2;
  p16.e_seg = 0;
  p16.u_seg = 1;
  //p1.uvar=NULL;
  checked_malloc (p16.uvar, pattern_var_t, 1, sizeof (pattern_var_t), return;);
  p16.uvar[0].size = 1;
  p16.uvar[0].reach = NULL;
  p16.lcons = -1; // no constraint
  p16.nr_y = 2;


  PI[0] = p1;
  PI[1] = p2;
  PI[2] = p3;
  PI[3] = p4;
  PI[4] = p5;
  PI[5] = p6;
  PI[6] = p7;
  PI[7] = p8;
  PI[8] = p9;
  PI[9] = p10;
  PI[10] = p11;
  PI[11] = p12;
  PI[12] = p13;
  PI[13] = p14;
  PI[14] = p15;
  PI[15] = p16;
  PI[16] = p17;




  pr->PI = PI;
  pr->PI_size = PI_size;
  checked_malloc (pr->active_patterns, bool, PI_size, sizeof (bool), return;);
  for (size_t i = 0; i < PI_size; i++)
    {
      pr->active_patterns[i] = false;
    }

  //pr->active_patterns[1] = true;
  //pr->active_patterns[16] = true;
  pr->nr_active = 0;
}



/* ********************************************************************** */
/* V. ucons_t */
/* ********************************************************************** */

/* ============================================================ */

/* Internal management */


size_t
get_pattern_type (ucons_internal_t *pr, size_t u_seg, size_t e_seg, size_t nr_y, pattern_kind name)
{
  //todo transform encoding on patterns!!!!! for the same pattern.kind various combinations of existential segments
  for (size_t i = 0; i < pr->PI_size; i++)
    {
      if (pr->PI[i].kind == name && pr->PI[i].u_seg == u_seg &&
          pr->PI[i].e_seg == e_seg && pr->PI[i].nr_y == nr_y)
        {
#ifndef NDEBUG
          fprintf (stdout, "\n@@@@ ucons: get_pattern_type returns pattern %zu \n", i);
          fflush (stdout);
#endif
          return i;
        }
    }
  fprintf (stdout, "\nUcons domain: Pattern not found!\n");
  fflush (stdout);
  ap_manager_raise_exception (pr->man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return 100;
}

//size_t get_next_e_seg_pattern(size_t i){
//	//todo transform encoding on patterns!!!!! for the same pattern.kind various combinations of existential segments
//	if (i == 1)
//		return 10;
//	if (i==8 || i==9) return i+1;
//	if (i==10 || i==11) return i+1;
//	return 0;
//}


/* ============================================================ */

/* generic allocation routine, returns ... */

inline ucons_t *
ucons_alloc_internal (ucons_internal_t * pr, size_t intdim, size_t realdim)
{

  ucons_t *r;
  size_t i;

  checked_malloc (r, ucons_t, sizeof (ucons_t), 1, return NULL;);


  r->datadim = intdim;
  r->segmentdim = realdim;

  r->econs = NULL;

  checked_malloc (r->n2p, pattern_key_set_t, realdim, sizeof (pattern_key_set_t), return NULL;);

  // MS: produces a memory error because of bad typing
  // r->n2p->size = r->segmentdim;
  r->size = r->segmentdim;
  for (i = 0; i < realdim; i++)
    {
      r->n2p[i].size = 0;
      r->n2p[i].p = NULL;
    }

  /* initialization to null of the hash table*/
  r->udcons = NULL;

  /* instantiation of patterns */
  return r;
}

/* returns a top ucons */
ucons_t *
ucons_alloc_top (ucons_internal_t * pr, size_t intdim, size_t dim)
{
  ucons_t *r = ucons_alloc_internal (pr, intdim, dim);
  size_t i, nr_y, u_seg;

  r->econs = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim, 0);
  pattern_t *p;
  for (p = r->udcons; p != NULL; p = p->hh.next)
    {

      u_seg = pr->PI[p->key.type].u_seg;
      nr_y = 0;
      for (i = 0; i < u_seg; i++)
        {
          nr_y += pr->PI[p->key.type].uvar[i].size;
        }

      p->dcons = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 2 * nr_y, 0);
    }
#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_alloc_top: returns \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, "\n");
#endif
  return r;
}

inline void
ucons_free_internal (ucons_internal_t * pr, ucons_t * a)
{
  pattern_t *s;
  size_t i, j;

  if (a)
    {
      if (a->econs)
        {
          ap_abstract0_free (pr->man_dcons, a->econs);
          a->econs = NULL;
        }

      while (a->udcons)
        {
          s = a->udcons;
          if (s->dcons)
            {
              ap_abstract0_free (pr->man_dcons, s->dcons);
              s->dcons = NULL;
            }
          HASH_DEL (a->udcons, s);
          free (s);
        }
      if (a->n2p)
        {
          for (i = 0; i < a->segmentdim; i++)
            for (j = 0; j < a->n2p[i].size; j++)
              if (a->n2p[i].p[j])
                {
                  free (a->n2p[i].p[j]);
                  //	free(a->n2p[i].p[j]);
                  a->n2p[i].p[j] = NULL;
                }

          free (a->n2p);
        }
      free (a);
    }
  return;
}

inline ucons_t *
ucons_copy_internal (ucons_internal_t * pr, ucons_t * a)
{
  arg_assert (a, return NULL;);
  //	arg_assert (a && a->econs && a->udcons, return NULL;);

  ucons_t *r = ucons_alloc_internal (pr, a->datadim, a->segmentdim);

  pattern_t* s;
  pattern_t* ra, *rt;
  size_t u_seg;
  size_t e_seg;
  size_t i;
  unsigned keylen;

  r->econs = ap_abstract0_copy (pr->man_dcons, a->econs);

  for (s = a->udcons; s != NULL; s = s->hh.next)
    {
      u_seg = pr->PI[s->key.type].u_seg;
      e_seg = pr->PI[s->key.type].e_seg;

      checked_malloc (ra, pattern_t, 1, (sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t)), return NULL;);
      memset (ra, 0, (sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t)));

      ra->dcons = NULL;
      ra->key.type = s->key.type;

      for (i = 0; i < u_seg; i++)
        ra->key.segments[i] = s->key.segments[i];
      for (i = 0; i < e_seg; i++)
        ra->key.segments[i + u_seg] = s->key.segments[i + u_seg];

      keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);

      HASH_FIND (hh, r->udcons, &ra->key, keylen, rt);
      if (rt)
        //HASH_DEL(r->udcons,rt);
        rt->dcons = (s->dcons) ?
        ap_abstract0_copy (pr->man_dcons, s->dcons) : NULL;

      else
        {
          ra->dcons = (s->dcons) ?
                  ap_abstract0_copy (pr->man_dcons, s->dcons) : NULL;
          HASH_ADD (hh, r->udcons, key, keylen, ra);
        }
    }
  r->n2p = pattern_key_set_copy (pr, a->n2p, a->segmentdim);

  return r;
}


/* ============================================================ */
/* Memory */

/* ============================================================ */
pattern_key_set_t *
pattern_key_set_copy (ucons_internal_t * pr, pattern_key_set_t * a, size_t size)
{

  pattern_key_set_t *r;
  size_t u_seg, e_seg;
  size_t i, j, k;



  checked_malloc (r, pattern_key_set_t, size, sizeof (pattern_key_set_t), return NULL;);

  r->size = size;
  for (i = 0; i < size; i++)
    {
      r[i].size = a[i].size;
      r[i].p = NULL;
      //copy a[i].p into r[i].p
      checked_malloc (r[i].p, pattern_key_t*, (r[i].size), sizeof (pattern_key_t*), return NULL;);
      for (j = 0; j < r[i].size; j++)
        {
          u_seg = pr->PI[a[i].p[j]->type].u_seg;
          e_seg = pr->PI[a[i].p[j]->type].e_seg;

          checked_malloc (r[i].p[j], pattern_key_t, 1, (sizeof (pattern_key_t)+ (u_seg + e_seg) * sizeof (size_t)), return NULL;);

          r[i].p[j]->type = a[i].p[j]->type;
          for (k = 0; k < (u_seg + e_seg); k++)
            r[i].p[j]->segments[k] = a[i].p[j]->segments[k];

        }

    }
  return r;
}

ucons_t *
ucons_copy (ap_manager_t * man, ucons_t * a)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_COPY, 0);
  return (a) ? ucons_copy_internal (pr, a) : NULL;
}

void
ucons_free (ap_manager_t * man, ucons_t * a)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_FREE, 0);
  ucons_free_internal (pr, a);
}

size_t
ucons_size (ap_manager_t * man, ucons_t * a)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_ASIZE, 0);
  return a->size;
}


/* ============================================================ */
/* Control of internal representation */
/* ============================================================ */

/* TODO: priority 3 */

/* Return the set of minimal elements */
void
ucons_minimize (ap_manager_t * man, ucons_t * a)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_MINIMIZE, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
}

/* TODO: priority 3 */

/* Return the set of canonical elements */
void
ucons_canonicalize (ap_manager_t * man, ucons_t * a)
{
  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
}

/* TODO: priority 0 */
int
ucons_hash (ap_manager_t * man, ucons_t * a)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_HASH, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return 0;
}

/* NOT IMPLEMENTED: do nothing */
void
ucons_approximate (ap_manager_t * man, ucons_t * a, int algorithm)
{
  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_APPROXIMATE, 0);
  // analysis parameters
  // - algorithm is split in 3 parts: (anon, segm, patterns)
  int anon = (algorithm & 0xF);
  int segm = (algorithm & 0xF0) >> 4;
  int patt = (algorithm & 0xFF00) >> 8;
  pr->max_anon = anon;
  pr->segm_anon = segm;
#ifndef NDEBUG
  fprintf (stdout, "\n@@@@ ucons_approximate: algorithm=%zu set max_anon=%zu, segm_anon=%zu patt = %zu\n",
           algorithm, pr->max_anon, pr->segm_anon, patt);
  fflush (stdout);
#endif

  /*
   * patts encode a number on 5 bits
   * set active patterns
   * 1 = 0 0 0 0 1 -> P11 coresp to PI[0]
   * 2 = 0 0 0 1 0 -> P21 coresp to PI[1]
   * 4 = 0 0 1 0 0 -> P12 //always comes with closure coresp to PI[3]
   * 16 = 1 0 0 0 0 -> P12succ coresp to PI[15]
   */
  size_t ind, k;

  for (size_t i = 0; i < pr->PI_size; i++)
    pr->active_patterns[i] = false;

  k = 1;
  pr->nr_active = 0;
  while (patt > 0)
    {
      ind = patt % 2;
      if (ind != 0)
        {
          pr->active_patterns[ind * k] = true;
          pr->nr_active++;
        }
      k = k * 2;
      patt = patt / 2;
    }
  if (pr->active_patterns[4])
    {
      if (!pr->active_patterns[1])
        {
          pr->active_patterns[1] = true;
          pr->nr_active++;
        }
    }


  //	if(pr->active_patterns[16]){
  //		pr->nr_active++;
  //		//pr->active_patterns[7] = true; //pr->nr_active++;
  //		//pr->active_patterns[8] = true; //pr->nr_active++;
  //	}

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_approximate: active patterns ");
  size_t i;
  for (i = 0; i < pr->PI_size; i++)
    if (pr->active_patterns[i])
      fprintf (stdout, "active %zu", i);
  fflush (stdout);
#endif

  if (a)
    {

      /*remove the constraints with the patterns no longer available */
      pattern_t * s = a->udcons;
      size_t types;

      while (s != NULL)
        {
          types = s->key.type;
          /*  atentie la incodaj .
           * PI[0] = forall y incodat 1
           * PI[1] \forall y1 = y2 incodat 2
           * PI[3] \forall y1<= y2 incodat 4
           * */
          types++;
          if (!pr->active_patterns[types])
            {
#ifndef NDEBUG1
              fprintf (stdout, "\n@@@@ ucons_approximate: a!=NULL pattern %zu removed\n", types);
              fflush (stdout);
#endif
              /*remove pattern from n2p */
              remove_pattern_n2p (pr, a, &s->key);
              /* remove pattern from a->udcons */
              pattern_t *t = s;
              s = s->hh.next;
              HASH_DEL (a->udcons, t);
              free (t);
            }
          else s = s->hh.next;

        }
    }
}

/* TODO: priority 3 */
bool
ucons_is_minimal (ap_manager_t * man, ucons_t * a)
{
  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return true;
}

/* TODO: priority 3 */
bool
ucons_is_canonical (ap_manager_t * man, ucons_t * a)
{
  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return true;
}

/* ============================================================ */
/* Basic Constructors */

/* ============================================================ */

ucons_t *
ucons_bottom (ap_manager_t * man, size_t intdim, size_t realdim)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_BOTTOM, 0);

  return NULL;
}

ucons_t *
ucons_top (ap_manager_t * man, size_t intdim, size_t realdim)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_TOP, 0);
  ucons_t *r = ucons_alloc_top (pr, intdim, realdim);
  return r;
}

/* TODO: priority 3 */

/* put constraints on data variables */
ucons_t *
ucons_of_box (ap_manager_t * man, size_t intdim, size_t realdim,
              ap_interval_t ** t)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_OF_BOX, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return NULL;
}

/* NOT IMPLEMENTED */
ucons_t *
ucons_of_generator_array (ap_manager_t * man, size_t intdim, size_t realdim,
                          ap_generator0_array_t * ar)
{
  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_ADD_RAY_ARRAY, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return NULL;
}


/* ============================================================ */
/* Accessors */

/* ============================================================ */

ap_dimension_t
ucons_dimension (ap_manager_t * man, ucons_t * a)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_DIMENSION, 0);
  ap_dimension_t r;
  // if (!a->a || (a->a && !a->a[0]))
  if (!a)
    {
      r.realdim = 0;
      r.intdim = 0;
    }
  else
    {
      r.realdim = a->segmentdim;
      r.intdim = a->datadim;
    }
  return r;
}

/* ********************************************************************** */
/* VI. ucons_array_t */
/* ********************************************************************** */


/* ============================================================ */
/* Managers */

/* ============================================================ */

void
ucons_internal_free (ucons_internal_t * pr)
{
  /* TODO: free PI */

  for (size_t i = 0; i < pr->PI_size; i++)
    {
      for (size_t j = 0; j < pr->PI[i].u_seg; j++)
        free (pr->PI[i].uvar->reach);
      free (pr->PI[i].uvar);
    }
  free (pr->PI);
  pr->PI = NULL;
  /* TODO: free data manager */
  pr->man_dcons = NULL;
  free (pr);
}

ap_manager_t *
ucons_manager_alloc (void)
{
  size_t i;
  ap_manager_t *man;
  ucons_internal_t *pr;
  char domain[128];

  bool oct, poly, ppl;
  oct = false;
  poly = false;
  ppl = false;
  pr = (ucons_internal_t *) malloc (sizeof (ucons_internal_t));
  assert (pr);


  FILE * f = fopen ("./cinv-ucons.txt", "r");


  if (!f)
    {

      fprintf (stdout, "File cinv-ucons.txt not found, set default data domain for analysis with ucons\n");

      poly = true;
      ppl = false;
      oct = false;

    }
  else
    {

      fprintf (stdout, "File cinv-ucons.txt found, read input values for analysis\n");
      fflush (stdout);
      /* format %s */
      //while (fgetc(f) != EOF) {
      char* adom = (char *) malloc (10 * sizeof (char));
      memset (adom, 0, 10 * sizeof (char));
      int l = fscanf (f, "%s", adom);
      if (strcmp (adom, "oct") == 0)
        {
          oct = true;
        }
      else if (strcmp (adom, "poly") == 0)
        {
          poly = true;
        }
      else if (strcmp (adom, "ppl") == 0)
        {
          ppl = true;
        }
      else if (strlen (adom) != 0)
        fprintf (stdout, "*** unknown domain %s (choose in oct|poly)\n", adom);
      fflush (stdout);
      free (adom);
      //}
      fclose (f);
    }

  /* read parameters from input file */
  //	fprintf (stdout, "params anon %zu,%zu\n", pr->max_anon,
  //			 pr->segm_anon);

  //	if(ppl){
  //			fprintf(stdout,"alege ppl");
  //			fflush(stdout);
  //			pr->man_dcons = ap_ppl_poly_manager_alloc(false);
  //	}
  //	else{
  if (oct)
    pr->man_dcons = oct_manager_alloc ();
  else if (poly)
    pr->man_dcons = pk_manager_alloc (false);
  else
    pr->man_dcons = ap_ppl_poly_manager_alloc (false);
  //	}


  i = snprintf (domain, 127, "0.1 with (dcons=%s)",
                ((UCONS_DCONS_DOMAIN ==
                  DOM_BOX) ? "Box" : ((UCONS_DCONS_DOMAIN ==
                                       DOM_OCT) ? "Oct" : "Polka")));
  domain[i] = '\0';

  man = ap_manager_alloc ("ucons", domain, pr,
                          (void (*)(void *)) ucons_internal_free);

  pr->man = man;

  man->funptr[AP_FUNID_COPY] = &ucons_copy;
  man->funptr[AP_FUNID_FREE] = &ucons_free;
  man->funptr[AP_FUNID_ASIZE] = &ucons_size;
  man->funptr[AP_FUNID_MINIMIZE] = &ucons_minimize;
  man->funptr[AP_FUNID_CANONICALIZE] = &ucons_canonicalize;
  man->funptr[AP_FUNID_HASH] = &ucons_hash;
  man->funptr[AP_FUNID_APPROXIMATE] = &ucons_approximate;
  man->funptr[AP_FUNID_FPRINT] = &ucons_fprint;
  man->funptr[AP_FUNID_FPRINTDIFF] = &ucons_fprintdiff;
  man->funptr[AP_FUNID_FDUMP] = &ucons_fdump;
  man->funptr[AP_FUNID_SERIALIZE_RAW] = &ucons_serialize_raw;
  man->funptr[AP_FUNID_DESERIALIZE_RAW] = &ucons_deserialize_raw;
  man->funptr[AP_FUNID_BOTTOM] = &ucons_bottom;
  man->funptr[AP_FUNID_TOP] = &ucons_top;
  man->funptr[AP_FUNID_OF_BOX] = &ucons_of_box;
  man->funptr[AP_FUNID_DIMENSION] = &ucons_dimension;
  man->funptr[AP_FUNID_IS_BOTTOM] = &ucons_is_bottom;
  man->funptr[AP_FUNID_IS_TOP] = &ucons_is_top;
  man->funptr[AP_FUNID_IS_LEQ] = &ucons_is_leq;
  man->funptr[AP_FUNID_IS_EQ] = &ucons_is_eq;
  man->funptr[AP_FUNID_IS_DIMENSION_UNCONSTRAINED] =
          &ucons_is_dimension_unconstrained;
  man->funptr[AP_FUNID_SAT_INTERVAL] = &ucons_sat_interval;
  man->funptr[AP_FUNID_SAT_LINCONS] = &ucons_sat_lincons;
  man->funptr[AP_FUNID_SAT_TCONS] = &ucons_sat_tcons;
  man->funptr[AP_FUNID_BOUND_DIMENSION] = &ucons_bound_dimension;
  man->funptr[AP_FUNID_BOUND_LINEXPR] = &ucons_bound_linexpr;
  man->funptr[AP_FUNID_BOUND_TEXPR] = &ucons_bound_texpr;
  man->funptr[AP_FUNID_TO_BOX] = &ucons_to_box;
  man->funptr[AP_FUNID_TO_LINCONS_ARRAY] = &ucons_to_lincons_array;
  man->funptr[AP_FUNID_TO_TCONS_ARRAY] = &ucons_to_tcons_array;
  man->funptr[AP_FUNID_TO_GENERATOR_ARRAY] = &ucons_to_generator_array;
  man->funptr[AP_FUNID_MEET] = &ucons_meet;
  man->funptr[AP_FUNID_MEET_ARRAY] = &ucons_meet_array;
  man->funptr[AP_FUNID_MEET_LINCONS_ARRAY] = &ucons_meet_lincons_array;
  man->funptr[AP_FUNID_MEET_TCONS_ARRAY] = &ucons_meet_tcons_array;
  man->funptr[AP_FUNID_JOIN] = &ucons_join;
  man->funptr[AP_FUNID_JOIN_ARRAY] = &ucons_join_array;
  man->funptr[AP_FUNID_ADD_RAY_ARRAY] = &ucons_add_ray_array;
  man->funptr[AP_FUNID_ASSIGN_LINEXPR_ARRAY] = &ucons_assign_linexpr_array;
  man->funptr[AP_FUNID_SUBSTITUTE_LINEXPR_ARRAY] =
          &ucons_substitute_linexpr_array;
  man->funptr[AP_FUNID_ASSIGN_TEXPR_ARRAY] = &ucons_assign_texpr_array;
  man->funptr[AP_FUNID_SUBSTITUTE_TEXPR_ARRAY] =
          &ucons_substitute_texpr_array;
  man->funptr[AP_FUNID_ADD_DIMENSIONS] = &ucons_add_dimensions;
  man->funptr[AP_FUNID_REMOVE_DIMENSIONS] = &ucons_remove_dimensions;
  man->funptr[AP_FUNID_PERMUTE_DIMENSIONS] = &ucons_permute_dimensions;
  man->funptr[AP_FUNID_FORGET_ARRAY] = &ucons_forget_array;
  man->funptr[AP_FUNID_EXPAND] = &ucons_expand;
  man->funptr[AP_FUNID_FOLD] = &ucons_fold;
  man->funptr[AP_FUNID_WIDENING] = &ucons_widening;
  man->funptr[AP_FUNID_CLOSURE] = &ucons_closure;

  for (i = 0; i < AP_EXC_SIZE; i++)
    ap_manager_set_abort_if_exception (man, i, false);

  initialize_PI (pr);
  return man;
}

ucons_t *
ucons_of_abstract0 (ap_abstract0_t * a)
{
  return (ucons_t *) a->value;
}

ap_abstract0_t *
abstract0_of_ucons (ap_manager_t * man, ucons_t * a)
{
  ap_abstract0_t *r = malloc (sizeof (ap_abstract0_t));
  assert (r);
  r->value = a;
  r->man = ap_manager_copy (man);
  return r;
}
