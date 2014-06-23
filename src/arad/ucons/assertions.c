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

#include <assert.h>
#include "ucons.h"
#include "ucons_internal.h"
#include "apron2shape.h"
#include "ap_generic.h"
#include "ap_linexpr0.h"
#include "ap_pcons0.h"

ucons_t *
build_constraint (ucons_internal_t * pr, ucons_t * r,
                  ap_linexpr0_t *lexpr, ap_scalar_t *code)
{

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ [ucons_]build_constraint: with code=(");
  ap_scalar_print (code);
  fprintf (stdout, ")\n and lexpr=[");
  ap_linexpr0_fprint (stdout, lexpr, NULL);
  fprintf (stdout, "]\n");
  fflush (stdout);
#endif

  if (ap_scalar_equal_int (code, 1))
    {
      r = build_const_1 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 2))
    {
      r = build_const_2 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 3))
    {
      r = build_const_3 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 4))
    {
      r = build_const_4 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 5))
    {
      r = build_const_5 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 6))
    {
      r = build_const_6 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 7))
    {
      r = build_const_7 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 8))
    {
      r = build_const_8 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 10))
    {
      r = build_const_10 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 11))
    {
      r = build_const_11 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 20))
    {
      r = build_const_20 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 21))
    {
      r = build_const_21 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 22))
    {
      r = build_const_22 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 23))
    {
      r = build_const_23 (pr, r, lexpr);
    }
  else if (ap_scalar_equal_int (code, 24))
    {
      r = build_const_24 (pr, r, lexpr);
    }

#ifndef NDEBUG1
  fprintf (stdout, "\n build_constraint returns:");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif

  return r;
}

void
ucons_build_extract_dimensions (ap_linexpr0_t* lexpr,
                                ap_dim_t* nmain, size_t *nmain_size,
                                ap_dim_t* naux, size_t *naux_size)
{
  ap_coeff_t *coeff;
  size_t i, dim;
  size_t smain, saux;

  smain = saux = 0;

  ap_linexpr0_ForeachLinterm (lexpr, i, dim, coeff)
  {
    if (coeff && !ap_coeff_zero (coeff))
      {
        if (coeff->discr == AP_COEFF_SCALAR)
          {
            ap_scalar_t* scoeff = coeff->val.scalar;
            int cmp = ap_scalar_cmp_int (scoeff, 0);
            if (cmp > 0 && smain < *nmain_size)
              nmain[smain++] = dim;
            else if (cmp < 0 && saux < *naux_size)
              naux[saux++] = dim;
          }
      }
  }
  *nmain_size = smain;
  *naux_size = saux;
}

ucons_t *
build_const_2 (ucons_internal_t * pr, ucons_t * r,
               ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = 1;
  saux = 0;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  naux = NULL;
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1, return r;);

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_build \\forall y1,y2:[%zu]. y1+1=y2 ==> d(y1)+d(y2)=1\n",
           nmain[0] - r->datadim);
  fflush (stdout);
#endif
  /*  d(y1) + d(y2) = 1*/
  ap_abstract0_t *datad = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 4, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t y2 = r->datadim + 2 * r->segmentdim + 1;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 2;
  ap_dim_t dy2 = r->datadim + 2 * r->segmentdim + 3;

  ap_dim_t ln = nmain[0] + r->segmentdim; // nmain[0] contains r->datadim

  ap_lincons0_array_t arra = ap_lincons0_array_make (3);
  /*
   *  (d(y2)+ d(y1)) -1 == 0
   */
  arra.p[0].constyp = AP_CONS_EQ;
  arra.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy1, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy2, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[0].linexpr0, -1);
  /*
   *  y2 - y1 - 1 == 0
   */
  arra.p[1].constyp = AP_CONS_EQ;
  arra.p[1].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[1].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y2, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[1].linexpr0, -1);

  /* l[n] >= 2 */
  arra.p[2].constyp = AP_CONS_SUPEQ;
  arra.p[2].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[2].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[2].linexpr0, -2);

  datad = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                           true, datad, &arra);

  ap_lincons0_array_clear (&arra);


  /* meet with r->econs */

  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 4, 0);
  dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[2] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[3] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *auxa = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  datad = ap_abstract0_meet (pr->man_dcons, true, datad, auxa);

  /* add to r->udcons hash table*/
  pattern_key_t *looka = NULL;
  checked_malloc (looka, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  looka->type = get_pattern_type (pr, 1, 0, 2, pattern_succ_1_2);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  looka->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, looka, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = looka->type;
      //for (size_t i=0 ; i<(u_seg); i++)
      a->key.segments[0] = looka->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, looka);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);
    }

  free (looka);
  looka = NULL;
  ap_abstract0_free (pr->man_dcons, datad);

  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;
}

ucons_t*
build_const_3 (ucons_internal_t * pr, ucons_t * r,
               ap_linexpr0_t *lexpr)
{
  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = 1;
  saux = 1;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  checked_malloc (naux, ap_dim_t, sizeof (ap_dim_t), saux, return NULL;);
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1 && saux == 1, return r;);

  // k is given by the naux[0] dimension
#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ucons_build \\forall y = l[n%zu] - 1. d(y) = 1 - x%zu\n",
           nmain[0] - r->datadim, naux[0]);
  fflush (stdout);
#endif

  pattern_key_t *look = NULL;
  checked_malloc (look, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  look->type = get_pattern_type (pr, 1, 0, 1, pattern_1_lx_1);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  ap_abstract0_t * data = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 2, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 1;


  ap_dim_t ln = nmain[0] + r->segmentdim; // nmain dimensions includes r->datadim
  ap_dim_t dn = nmain[0];

  ap_lincons0_array_t arr;
  arr = ap_lincons0_array_make (3);
  /*
   *  y - l[n] + 1 == 0
   */
  arr.p[0].constyp = AP_CONS_EQ;
  arr.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arr.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ln, -1);
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, y1, 1);
  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 1);

  // P_1_lx_1 --> d(y) + k -1 ==0
  arr.p[1].constyp = AP_CONS_EQ;
  arr.p[1].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arr.p[1].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy1, 1);
  ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, naux[0], 1);
  ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, -1);

  /* l[n] >= 2 */
  arr.p[2].constyp = AP_CONS_SUPEQ;
  arr.p[2].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arr.p[2].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arr.p[2].linexpr0, -2);

  data = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                          true, data, &arr);

  ap_lincons0_array_clear (&arr);

  /* meet with r->econs */
  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 2, 0);
  dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *aux = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  data = ap_abstract0_meet (pr->man_dcons, true, data, aux);


  look->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, look, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = look->type;
      //for (size_t i=0 ; i<(u_seg); i++)
      a->key.segments[0] = look->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, look);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);
    }

  ap_abstract0_free (pr->man_dcons, data);
  data = NULL;

  free (look);
  look = NULL;

  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;

}

ucons_t*
build_const_4 (ucons_internal_t * pr, ucons_t * r,
               ap_linexpr0_t* lexpr)
{

#ifndef NDEBUG1
  fprintf (stdout, "\n NOT YET IMPLEMENTED!!!\n");
  fflush (stdout);
#endif
  return r;
}

ucons_t *
build_const_5 (ucons_internal_t * pr, ucons_t * r,
               ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = 2;
  saux = 1;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  checked_malloc (naux, ap_dim_t, sizeof (ap_dim_t), saux, return NULL;);
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 2 && (saux == 1 || saux == 0), return r;);

#ifndef NDEBUG1
  // c is naux[0]
  if (saux != 0)
    {
      fprintf (stdout, "\n@@@@ucons_build \\forall y1 in [n%zu], y2 in [n%zu]. y1 = y2 ==> d(y1) = d(y2)+x%zu\n",
               nmain[0] - r->datadim, nmain[1] - r->datadim, naux[0]);
    }
  else
    {
      fprintf (stdout, "\n@@@@ucons_build \\forall y1 in [n%zu], y2 in [n%zu]. y1 = y2 ==> d(y1) = d(y2)+2 \n",
               nmain[0] - r->datadim, nmain[1] - r->datadim);
    }
  fflush (stdout);
#endif

  //    d(y1) - d(y2) - c = 0
  ap_abstract0_t *datad = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 4, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t y2 = r->datadim + 2 * r->segmentdim + 1;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 2;
  ap_dim_t dy2 = r->datadim + 2 * r->segmentdim + 3;


  ap_lincons0_array_t arra = ap_lincons0_array_make (2);

  /*   d(y1) - d(y2) + c = 0*/

  arra.p[0].constyp = AP_CONS_EQ;
  arra.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy1, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy2, -1);
  if (saux != 0)
    {
      ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, naux[0], 1);
      ap_linexpr0_set_cst_scalar_int (arra.p[0].linexpr0, 0);
    }
  else
    {
      ap_linexpr0_set_cst_scalar_int (arra.p[0].linexpr0, 2);
    }

  /*  y2 == y1*/

  arra.p[1].constyp = AP_CONS_EQ;
  arra.p[1].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[1].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y2, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[1].linexpr0, 0);


  datad = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                           true, datad, &arra);

  ap_lincons0_array_clear (&arra);

  /* meet with r->econs */

  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 4, 0);
  dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[2] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[3] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *auxa = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  datad = ap_abstract0_meet (pr->man_dcons, true, datad, auxa);

  pattern_key_t *looka = NULL;
  checked_malloc (looka, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 2 * sizeof (size_t), return NULL;);
  looka->type = get_pattern_type (pr, 2, 0, 2, pattern_2_1);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  // sort the node dimensions
  if (nmain[0] > nmain[1])
    {
      size_t p = nmain[0];
      nmain[0] = nmain[1];
      nmain[1] = p;
    }

  looka->segments[0] = nmain[0] - r->datadim;
  looka->segments[1] = nmain[1] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, looka, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(2) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(2) * sizeof (size_t));
      a->key.type = looka->type;
      //for (size_t i=0 ; i<(u_seg); i++)
      a->key.segments[0] = looka->segments[0];
      a->key.segments[1] = looka->segments[1];
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, looka);

      //       ucons_fprint (stdout,pr->man,r, NULL);

    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);
      // 		ucons_fprint (stdout,pr->man,r, NULL);
    }


  free (looka);
  looka = NULL;
  ap_abstract0_free (pr->man_dcons, datad);

  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;
}

/*
ucons_t *
build_const_6(ucons_internal_t * pr, ucons_t * r,
        ap_linexpr0_t *lexpr) {

    ap_dim_t *nmain; // main dimensions for the pattern
    ap_dim_t *naux; // aux dimensions for the pattern
    size_t smain, saux;

    smain = 2;
    saux = 0;
    checked_malloc(nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
    checked_malloc(naux, ap_dim_t, sizeof (ap_dim_t), saux, return NULL;);
    ucons_build_extract_dimensions(lexpr, nmain, &smain, naux, &saux);

    arg_assert(smain == 2 && saux == 0, return r;);

    // c is naux[0]
#ifndef NDEBUG1
    fprintf(stdout, "\n@@@@ucons_build \\forall y1 in [n%zu], y2 in [n%zu]. y1 = y2 ==> d(y1) = d(y2) %zu\n",
            nmain[0] - r->datadim, nmain[1] - r->datadim, naux[0]);
    fflush(stdout);
#endif

      d(y1) - d(y2)= 0
    ap_abstract0_t *datad = ap_abstract0_top(pr->man_dcons, r->datadim + 2 * r->segmentdim + 4, 0);
    ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
    ap_dim_t y2 = r->datadim + 2 * r->segmentdim + 1;
    ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 2;
    ap_dim_t dy2 = r->datadim + 2 * r->segmentdim + 3;


    ap_lincons0_array_t arra = ap_lincons0_array_make(2);

 *   d(y1) - d(y2) = 0

    arra.p[0].constyp = AP_CONS_EQ;
    arra.p[0].linexpr0 =
            ap_linexpr0_alloc(AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
    arra.p[0].scalar = NULL;
    ap_linexpr0_set_coeff_scalar_int(arra.p[0].linexpr0, dy1, 1);
    ap_linexpr0_set_coeff_scalar_int(arra.p[0].linexpr0, dy2, -1);
    ap_linexpr0_set_cst_scalar_int(arra.p[0].linexpr0, 0);

 *  y2 == y1

    arra.p[1].constyp = AP_CONS_EQ;
    arra.p[1].linexpr0 =
            ap_linexpr0_alloc(AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
    arra.p[1].scalar = NULL;
    ap_linexpr0_set_coeff_scalar_int(arra.p[1].linexpr0, y1, -1);
    ap_linexpr0_set_coeff_scalar_int(arra.p[1].linexpr0, y2, 1);
    ap_linexpr0_set_cst_scalar_int(arra.p[1].linexpr0, 0);


    datad = ap_abstract0_meet_lincons_array(pr->man_dcons,
            true, datad, &arra);

    ap_lincons0_array_clear(&arra);

     meet with r->econs

    ap_dimchange_t dimadd;
    ap_dimchange_init(&dimadd, 4, 0);
    dimadd.dim = (ap_dim_t *) malloc(4 * sizeof (ap_dim_t));
    dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
    dimadd.dim[1] = r->datadim + 2 * r->segmentdim;
    dimadd.dim[2] = r->datadim + 2 * r->segmentdim;
    dimadd.dim[3] = r->datadim + 2 * r->segmentdim;

    ap_abstract0_t *auxa = ap_abstract0_add_dimensions(pr->man_dcons, false, r->econs, &dimadd, false);

    ap_dimchange_clear(&dimadd);

    datad = ap_abstract0_meet(pr->man_dcons, true, datad, auxa);

    pattern_key_t *looka = NULL;
    checked_malloc(looka, pattern_key_t, 1,
            sizeof (pattern_key_t) + 2 * sizeof (size_t), return NULL;);
    looka->type = get_pattern_type(pr, 2, 0, 2, pattern_2_1);
    unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

    // sort the node dimensions
    if (nmain[0] > nmain[1]) {
        size_t p = nmain[0];
        nmain[0] = nmain[1];
        nmain[1] = p;
    }

    looka->segments[0] = nmain[0] - r->datadim;
    looka->segments[1] = nmain[1] - r->datadim;
    pattern_t *a = NULL;
    HASH_FIND(hh, r->udcons, looka, keylen, a);
    if (!a) {
        checked_malloc(a, pattern_t, 1,
                sizeof (pattern_t)+(2) * sizeof (size_t), return NULL;);
        memset(a, 0, sizeof (pattern_t)+(2) * sizeof (size_t));
        a->key.type = looka->type;
        //for (size_t i=0 ; i<(u_seg); i++)
        a->key.segments[0] = looka->segments[0];
        a->key.segments[1] = looka->segments[1];
        a->dcons = ap_abstract0_copy(pr->man_dcons, datad);

        HASH_ADD(hh, r->udcons, key, keylen, a);
        r = add_pattern_n2p(pr, r, looka);
    } else {
        if (a->dcons != NULL)
            ap_abstract0_free(pr->man_dcons, a->dcons);
        a->dcons = ap_abstract0_copy(pr->man_dcons, datad);
    }


    free(looka);
    looka = NULL;
    ap_abstract0_free(pr->man_dcons, datad);

    if (nmain) free(nmain);
    if (naux) free(naux);

    return r;
}

 */




ucons_t *
build_const_8 (ucons_internal_t * pr, ucons_t * r,
               ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = 1;
  saux = 0;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  naux = NULL;
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1, return r;);

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_build \\forall y1,y2:[%zu]. y1+1=y2 ==> d(y1)+d(y2)=1\n",
           nmain[0] - r->datadim);
  fflush (stdout);
#endif
  /*  d(y1) + d(y2) = 1*/
  ap_abstract0_t *datad = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 4, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t y2 = r->datadim + 2 * r->segmentdim + 1;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 2;
  ap_dim_t dy2 = r->datadim + 2 * r->segmentdim + 3;

  ap_dim_t ln = nmain[0] + r->segmentdim; // nmain[0] contains r->datadim

  ap_lincons0_array_t arra = ap_lincons0_array_make (4);
  /*
   *  (d(y2)+ d(y1)) -1 == 0 !!!!!!! modifica la loc
   */
  arra.p[0].constyp = AP_CONS_SUPEQ;
  arra.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy2, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, y1, -1); // ??
  ap_linexpr0_set_cst_scalar_int (arra.p[0].linexpr0, 0); // -1
  /*
   *  y2 - y1 - 1 == 0
   */
  arra.p[1].constyp = AP_CONS_EQ;
  arra.p[1].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[1].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y2, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[1].linexpr0, -1);

  /* l[n] >= 3 */
  arra.p[2].constyp = AP_CONS_SUPEQ;
  arra.p[2].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[2].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[2].linexpr0, -3);

  /*
   *  y1 - 1 >= 0
   */
  arra.p[3].constyp = AP_CONS_SUPEQ;
  arra.p[3].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[3].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, y1, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[3].linexpr0, -1);

  //    /*
  //     *  d(y1)) - d(n) -1 >= 0 !!!!!!! modifica la loc
  //     */
  //    arra.p[4].constyp = AP_CONS_SUPEQ;
  //    arra.p[4].linexpr0 =
  //    		ap_linexpr0_alloc(AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  //    arra.p[4].scalar = NULL;
  //    ap_linexpr0_set_coeff_scalar_int(arra.p[4].linexpr0, dy1, 1);
  //    ap_linexpr0_set_coeff_scalar_int(arra.p[4].linexpr0, nmain[0], -1);
  //    ap_linexpr0_set_cst_scalar_int(arra.p[4].linexpr0, -1);


  datad = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                           true, datad, &arra);

  ap_lincons0_array_clear (&arra);


  /* meet with r->econs */

  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 4, 0);
  dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[2] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[3] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *auxa = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  datad = ap_abstract0_meet (pr->man_dcons, true, datad, auxa);

  /* add to r->udcons hash table*/
  pattern_key_t *looka = NULL;
  checked_malloc (looka, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  looka->type = get_pattern_type (pr, 1, 0, 2, pattern_succ_1_2);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  looka->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, looka, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = looka->type;
      //for (size_t i=0 ; i<(u_seg); i++)
      a->key.segments[0] = looka->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, looka);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);
    }

  free (looka);
  looka = NULL;
  ap_abstract0_free (pr->man_dcons, datad);


  //*** ************************************ build closure : *********************************************/////

  // \forall y in n. y=1 d(y)= d(n)+ 1

  datad = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 2, 0);

  y1 = r->datadim + 2 * r->segmentdim;
  dy1 = r->datadim + 2 * r->segmentdim + 1;
  ap_dim_t dn = nmain[0]; // includes r->datadim
  ln = nmain[0] + r->segmentdim;

  arra = ap_lincons0_array_make (3);
  /*
   *  d(n) = d(y) + 1
   */
  arra.p[0].constyp = AP_CONS_SUPEQ; //EQ
  arra.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy1, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dn, -1);
  ap_linexpr0_set_cst_scalar_int (arra.p[0].linexpr0, 0); //-1

  /*
   * y1 = 1
   */
  arra.p[1].constyp = AP_CONS_EQ;
  arra.p[1].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[1].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y1, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[1].linexpr0, -1);


  /*
   * y1 <= l[n] - 1 il elimin ln >= 2
   */
  arra.p[2].constyp = AP_CONS_SUPEQ;
  arra.p[2].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[2].scalar = NULL;
  //ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, y1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[2].linexpr0, -2);


  /*
   * ln >= 2

  arra.p[3].constyp = AP_CONS_SUPEQ;
  arra.p[3].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[3].scalar = NULL;
  //ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, y1, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[3].linexpr0, -2);

   */

  datad = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                           true, datad, &arra);
  ap_lincons0_array_clear (&arra);

  /* meet with r->econs */
  ap_dimchange_init (&dimadd, 2, 0);
  dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;

  auxa = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  datad = ap_abstract0_meet (pr->man_dcons, true, datad, auxa);

  looka = NULL;
  checked_malloc (looka, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  looka->type = get_pattern_type (pr, 1, 0, 1, pattern_1_l1);
  keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  looka->segments[0] = nmain[0] - r->datadim;
  a = NULL;
  HASH_FIND (hh, r->udcons, looka, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = looka->type;
      a->key.segments[0] = looka->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, looka);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);
    }

  free (looka);
  looka = NULL;
  ap_abstract0_free (pr->man_dcons, datad);

  //TODO add pattern_1_lx to express the prop on the last elem



  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;
}

ucons_t *
build_const_6 (ucons_internal_t * pr, ucons_t * r,
               ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = 1;
  saux = 0;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  naux = NULL;
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1, return r;);

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_build \\forall y1,y2:[%zu]. y1+1=y2 ==> d(y1)+d(y2)=1\n",
           nmain[0] - r->datadim);
  fflush (stdout);
#endif
  /*  d(y1) + d(y2) = 1*/
  ap_abstract0_t *datad = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 4, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t y2 = r->datadim + 2 * r->segmentdim + 1;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 2;
  ap_dim_t dy2 = r->datadim + 2 * r->segmentdim + 3;

  ap_dim_t ln = nmain[0] + r->segmentdim; // nmain[0] contains r->datadim

  ap_lincons0_array_t arra = ap_lincons0_array_make (5);
  /*
   *  (d(y2)+ d(y1)) -1 == 0
   */
  arra.p[0].constyp = AP_CONS_EQ;
  arra.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy2, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[0].linexpr0, -1); //
  /*
   *  y2 - y1 - 1 == 0
   */
  arra.p[1].constyp = AP_CONS_EQ;
  arra.p[1].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[1].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y2, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[1].linexpr0, -1);

  /* l[n] >= 3 */
  arra.p[2].constyp = AP_CONS_SUPEQ;
  arra.p[2].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[2].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[2].linexpr0, -3);

  /*
   *  y1 - 1 >= 0
   */
  arra.p[3].constyp = AP_CONS_SUPEQ;
  arra.p[3].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[3].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, y1, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[3].linexpr0, -1);

  /*
   *  d(y1)) - d(n) -1 >= 0
   */
  arra.p[4].constyp = AP_CONS_SUPEQ;
  arra.p[4].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[4].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[4].linexpr0, dy1, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[4].linexpr0, nmain[0], -1);
  ap_linexpr0_set_cst_scalar_int (arra.p[4].linexpr0, -1);


  datad = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                           true, datad, &arra);

  ap_lincons0_array_clear (&arra);


  /* meet with r->econs */

  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 4, 0);
  dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[2] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[3] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *auxa = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  datad = ap_abstract0_meet (pr->man_dcons, true, datad, auxa);

  /* add to r->udcons hash table*/
  pattern_key_t *looka = NULL;
  checked_malloc (looka, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  looka->type = get_pattern_type (pr, 1, 0, 2, pattern_succ_1_2);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  looka->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, looka, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = looka->type;
      //for (size_t i=0 ; i<(u_seg); i++)
      a->key.segments[0] = looka->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, looka);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);
    }

  free (looka);
  looka = NULL;
  ap_abstract0_free (pr->man_dcons, datad);


  //*** ************************************ build closure : *********************************************/////

  // \forall y in n. y=1 d(y)= d(n)+ 1

  datad = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 2, 0);

  y1 = r->datadim + 2 * r->segmentdim;
  dy1 = r->datadim + 2 * r->segmentdim + 1;
  ap_dim_t dn = nmain[0]; // includes r->datadim
  ln = nmain[0] + r->segmentdim;

  arra = ap_lincons0_array_make (3);
  /*
   *  d(n) = d(y) + 1
   */
  arra.p[0].constyp = AP_CONS_SUPEQ; //EQ
  arra.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy1, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dn, -1);
  ap_linexpr0_set_cst_scalar_int (arra.p[0].linexpr0, 0); //-1

  /*
   * y1 = 1
   */
  arra.p[1].constyp = AP_CONS_EQ;
  arra.p[1].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[1].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y1, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[1].linexpr0, -1);


  /*
   * y1 <= l[n] - 1 il elimin ln >= 2
   */
  arra.p[2].constyp = AP_CONS_SUPEQ;
  arra.p[2].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[2].scalar = NULL;
  //ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, y1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[2].linexpr0, -2);


  /*
   * ln >= 2

  arra.p[3].constyp = AP_CONS_SUPEQ;
  arra.p[3].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[3].scalar = NULL;
  //ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, y1, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[3].linexpr0, -2);

   */

  datad = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                           true, datad, &arra);
  ap_lincons0_array_clear (&arra);

  /* meet with r->econs */
  ap_dimchange_init (&dimadd, 2, 0);
  dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;

  auxa = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  datad = ap_abstract0_meet (pr->man_dcons, true, datad, auxa);

  looka = NULL;
  checked_malloc (looka, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  looka->type = get_pattern_type (pr, 1, 0, 1, pattern_1_l1);
  keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  looka->segments[0] = nmain[0] - r->datadim;
  a = NULL;
  HASH_FIND (hh, r->udcons, looka, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = looka->type;
      a->key.segments[0] = looka->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, looka);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);
    }

  free (looka);
  looka = NULL;
  ap_abstract0_free (pr->man_dcons, datad);

  //TODO add pattern_1_lx to express the prop on the last elem



  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;
}

ucons_t*
build_const_7 (ucons_internal_t * pr, ucons_t * r,
               ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = saux = 1;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  checked_malloc (naux, ap_dim_t, sizeof (ap_dim_t), saux, return NULL;);
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1 && saux == 1, return r;);

  // v is naux[0]
#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ucons_build \\forall y in [n%zu]. d(y) = x%zu\n",
           nmain[0] - r->datadim, naux[0]);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif

  ap_abstract0_t * data = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 2, 0);
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 1;

  ap_lincons0_array_t arr;
  arr = ap_lincons0_array_make (1);
  /*
   *  d(y) == v
   */
  arr.p[0].constyp = AP_CONS_EQ;
  arr.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arr.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, naux[0], 1);
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy1, -1);

  data = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                          true, data, &arr);

  ap_lincons0_array_clear (&arr);

  /* meet with r->econs */
  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 2, 0);
  dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *aux = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  data = ap_abstract0_meet (pr->man_dcons, true, data, aux);

  pattern_key_t *look = NULL;
  checked_malloc (look, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  look->type = get_pattern_type (pr, 1, 0, 1, pattern_1);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  look->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, look, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = look->type;
      //for (size_t i=0 ; i<(u_seg); i++)
      a->key.segments[0] = look->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, look);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);
    }

  free (look);
  look = NULL;
  ap_abstract0_free (pr->man_dcons, data);

  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;
}

ucons_t*
build_const_1 (ucons_internal_t * pr, ucons_t * r,
               ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = 1;
  saux = 0;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  naux = NULL;
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1, return r;);

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_build \\forall y in [n%zu]. d(y) = y\n",
           nmain[0] - r->datadim);
  fflush (stdout);
#endif

  ap_abstract0_t * data = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 2, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 1;

  ap_lincons0_array_t arr;
  arr = ap_lincons0_array_make (1);
  /*
   *  y - d(y) == 0
   */
  arr.p[0].constyp = AP_CONS_EQ;
  arr.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arr.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, y1, 1);
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy1, -1);
  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

  data = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                          true, data, &arr);

  ap_lincons0_array_clear (&arr);

  /* meet with r->econs */
  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 2, 0);
  dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *aux = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  data = ap_abstract0_meet (pr->man_dcons, true, data, aux);

  pattern_key_t *look = NULL;
  checked_malloc (look, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  look->type = get_pattern_type (pr, 1, 0, 1, pattern_1);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  look->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, look, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = look->type;
      //for (size_t i=0 ; i<(u_seg); i++)
      a->key.segments[0] = look->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, look);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);
    }


  free (look);
  look = NULL;
  ap_abstract0_free (pr->man_dcons, data);

  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;
}

ucons_t*
build_const_10 (ucons_internal_t * pr, ucons_t * r,
                ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = 1;
  saux = 0;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  naux = NULL;
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1, return r;);

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_build \\forall y1,y2 in [n%zu]. y2=y1+1 ==> d(y2)-d(y1)>=y1  \n",
           nmain[0] - r->datadim);
  fflush (stdout);
#endif

  /*  d(y2) - d(y1) >= y1 - 2 */
  ap_abstract0_t *datad = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 4, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t y2 = r->datadim + 2 * r->segmentdim + 1;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 2;
  ap_dim_t dy2 = r->datadim + 2 * r->segmentdim + 3;
  ap_dim_t ln = nmain[0] + r->segmentdim;

  ap_lincons0_array_t arra = ap_lincons0_array_make (5);
  /*
   *  (d(y2) - d(y1)) - y1 + 2 >= 0//(d(y2) - d(y1)) - y1 >= 0/
   */
  arra.p[0].constyp = AP_CONS_SUPEQ;
  arra.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy2, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, y1, -1);
  ap_linexpr0_set_cst_scalar_int (arra.p[0].linexpr0, 0);
  /*
   *  y2 - y1 -1 == 0
   */
  arra.p[1].constyp = AP_CONS_EQ;
  arra.p[1].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[1].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y2, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[1].linexpr0, -1);

  /*
   * y1 >= 1
   */
  arra.p[2].constyp = AP_CONS_SUPEQ;
  arra.p[2].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[2].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, y1, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[2].linexpr0, -1);

  /*
   * y2 <= l[n] - 1.
   */
  arra.p[3].constyp = AP_CONS_SUPEQ;
  arra.p[3].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[3].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, y2, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[3].linexpr0, -1);


  /*
   *  l[n] >= 3 .
   */
  arra.p[4].constyp = AP_CONS_SUPEQ;
  arra.p[4].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[4].scalar = NULL;
  //ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, y2, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[4].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[4].linexpr0, -3);


  datad = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                           true, datad, &arra);

  ap_lincons0_array_clear (&arra);

  /*

       meet with r->econs*/

  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 4, 0);
  dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[2] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[3] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *auxa = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  datad = ap_abstract0_meet (pr->man_dcons, true, datad, auxa);


  pattern_key_t *looka = NULL;
  checked_malloc (looka, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  looka->type = get_pattern_type (pr, 1, 0, 2, pattern_succ_1_2);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  looka->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, looka, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = looka->type;
      //for (size_t i=0 ; i<(u_seg); i++)
      a->key.segments[0] = looka->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, looka);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);
    }


  free (looka);
  looka = NULL;
  ap_abstract0_free (pr->man_dcons, datad);

  //*** ************************************ build closure : *********************************************/////

  // \forall y in n. y=1 d(y)= d(n)+ 1

  datad = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 2, 0);

  y1 = r->datadim + 2 * r->segmentdim;
  dy1 = r->datadim + 2 * r->segmentdim + 1;
  ap_dim_t dn = nmain[0]; // includes r->datadim
  ln = nmain[0] + r->segmentdim;

  arra = ap_lincons0_array_make (3);
  /*
   *  d(n) = d(y) + 1
   */
  arra.p[0].constyp = AP_CONS_EQ;
  arra.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy1, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dn, -1);
  ap_linexpr0_set_cst_scalar_int (arra.p[0].linexpr0, -1);

  /*
   * y1 = 1
   */
  arra.p[1].constyp = AP_CONS_EQ;
  arra.p[1].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[1].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y1, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[1].linexpr0, -1);


  /*
   * y1 <= l[n] - 1 il elimin ln >= 2
   */
  arra.p[2].constyp = AP_CONS_SUPEQ;
  arra.p[2].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[2].scalar = NULL;
  //ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, y1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[2].linexpr0, -2);


  /*
   * ln >= 2

  arra.p[3].constyp = AP_CONS_SUPEQ;
  arra.p[3].linexpr0 =
                  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[3].scalar = NULL;
  //ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, y1, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[3].linexpr0, -2);

   */

  datad = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                           true, datad, &arra);
  ap_lincons0_array_clear (&arra);

  /* meet with r->econs */
  ap_dimchange_init (&dimadd, 2, 0);
  dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;

  auxa = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  datad = ap_abstract0_meet (pr->man_dcons, true, datad, auxa);

  looka = NULL;
  checked_malloc (looka, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  looka->type = get_pattern_type (pr, 1, 0, 1, pattern_1_l1);
  keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  looka->segments[0] = nmain[0] - r->datadim;
  a = NULL;
  HASH_FIND (hh, r->udcons, looka, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = looka->type;
      a->key.segments[0] = looka->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, looka);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);
    }

  free (looka);
  looka = NULL;
  ap_abstract0_free (pr->man_dcons, datad);



  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;
}

ucons_t *
build_const_20 (ucons_internal_t * pr, ucons_t * r,
                ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = saux = 1;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  checked_malloc (naux, ap_dim_t, sizeof (ap_dim_t), saux, return NULL;);
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1 && saux <= 1, return r;);

  // sorted(n):
  //    \forall y1,y2. y1 <= y2 ==> d(y1) <= d(y2)
  // && \forall y. d(n) <= d(y)

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_build sorted(n%zu) && \\forall y:[n%zu]. d(y) <= d(%zu)\n",
           nmain[0] - r->datadim,
           nmain[0] - r->datadim, (saux == 1) ? naux[0] : 256);
  fflush (stdout);
#endif

  // \forall y1,y2. y1 <= y2 ==> d(y1) <= d(y2)
  /*  d(y1) <= d(y2) */
  ap_abstract0_t *datad = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 4, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t y2 = r->datadim + 2 * r->segmentdim + 1;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 2;
  ap_dim_t dy2 = r->datadim + 2 * r->segmentdim + 3;
  ap_dim_t ln = nmain[0] + r->segmentdim; // includes r->datadim
  ap_dim_t dn = nmain[0]; // includes r->datadim

  ap_lincons0_array_t arra = ap_lincons0_array_make (5);
  /*
   *  (d(y2) - d(y1)) >= 0
   */
  arra.p[0].constyp = AP_CONS_SUPEQ;
  arra.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy2, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[0].linexpr0, 0);
  /*
   *  y2 - y1 >= 0 // replace with y2 >=1
   */
  arra.p[1].constyp = AP_CONS_SUPEQ;
  arra.p[1].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[1].scalar = NULL;
  //ap_linexpr0_set_coeff_scalar_int(arra.p[1].linexpr0, y1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y2, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[1].linexpr0, -1);

  /*
   * y1 >= 1
   */
  arra.p[2].constyp = AP_CONS_SUPEQ;
  arra.p[2].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[2].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, y1, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[2].linexpr0, -1);

  /*
   * y2 <= l[n] - 1 //il elimin cu l[n] >= 2 .
   */
  arra.p[3].constyp = AP_CONS_SUPEQ;
  arra.p[3].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[3].scalar = NULL;
  //ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, y2, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[3].linexpr0, -2);
  /*
   *  (d(y1) - d(n)) >= 0
   */
  arra.p[4].constyp = AP_CONS_SUPEQ;
  arra.p[4].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[4].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[4].linexpr0, dn, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[4].linexpr0, dy1, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[4].linexpr0, 0);


  /*
   *  y2 - y1 >= 0

  arra.p[5].constyp = AP_CONS_SUPEQ;
  arra.p[5].linexpr0 =
              ap_linexpr0_alloc(AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arra.p[5].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int(arra.p[5].linexpr0, y1, -1);
  ap_linexpr0_set_coeff_scalar_int(arra.p[5].linexpr0, y2, 1);
  ap_linexpr0_set_cst_scalar_int(arra.p[5].linexpr0, 0);
   */

  datad = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                           true, datad, &arra);

  ap_lincons0_array_clear (&arra);

  /* meet with r->econs */
  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 4, 0);
  dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[2] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[3] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *auxa = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  datad = ap_abstract0_meet (pr->man_dcons, true, datad, auxa);

  pattern_key_t *looka = NULL;
  checked_malloc (looka, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  looka->type = get_pattern_type (pr, 1, 0, 2, pattern_1_2);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  looka->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, looka, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = looka->type;
      a->key.segments[0] = looka->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, looka);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);
    }


  free (looka);
  looka = NULL;
  ap_abstract0_free (pr->man_dcons, datad);


  // \forall y in n. d(n) <= d(y) <= naux
  /*  d(n) <= d(y) <= naux */
  datad = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 2, 0);

  y1 = r->datadim + 2 * r->segmentdim;
  dy1 = r->datadim + 2 * r->segmentdim + 1;
  dn = nmain[0]; // includes r->datadim
  ln = nmain[0] + r->segmentdim;

  arra = ap_lincons0_array_make ((saux == 0) ? 3 : 4);
  /*
   *  d(n) <= d(y) <= naux
   */
  arra.p[0].constyp = AP_CONS_SUPEQ;
  arra.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dy1, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[0].linexpr0, dn, -1);
  ap_linexpr0_set_cst_scalar_int (arra.p[0].linexpr0, 0);

  /*
   * y1 >= 1
   */
  arra.p[1].constyp = AP_CONS_SUPEQ;
  arra.p[1].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[1].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arra.p[1].linexpr0, y1, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[1].linexpr0, -1);


  /*
   * y1 <= l[n] - 1 il elimin ln >= 2
   */
  arra.p[2].constyp = AP_CONS_SUPEQ;
  arra.p[2].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[2].scalar = NULL;
  //ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, y1, -1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[2].linexpr0, -2);


  /*
   * ln >= 2

  arra.p[3].constyp = AP_CONS_SUPEQ;
  arra.p[3].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arra.p[3].scalar = NULL;
  //ap_linexpr0_set_coeff_scalar_int (arra.p[2].linexpr0, y1, 1);
  ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arra.p[3].linexpr0, -2);

   */
  if (saux == 1)
    {
      // v - d(y) >= 0
      arra.p[3].constyp = AP_CONS_SUPEQ;
      arra.p[3].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
      arra.p[3].scalar = NULL;
      ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, dy1, -1);
      ap_linexpr0_set_coeff_scalar_int (arra.p[3].linexpr0, naux[0], 1);
      ap_linexpr0_set_cst_scalar_int (arra.p[3].linexpr0, 0);
    }


  datad = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                           true, datad, &arra);

  ap_lincons0_array_clear (&arra);


  /* meet with r->econs */
  ap_dimchange_init (&dimadd, 2, 0);


  dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;

  auxa = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  datad = ap_abstract0_meet (pr->man_dcons, true, datad, auxa);

  looka = NULL;
  checked_malloc (looka, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  looka->type = get_pattern_type (pr, 1, 0, 1, pattern_1);
  keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  looka->segments[0] = nmain[0] - r->datadim;
  a = NULL;
  HASH_FIND (hh, r->udcons, looka, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = looka->type;
      a->key.segments[0] = looka->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, looka);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, datad);
    }


  free (looka);
  looka = NULL;
  ap_abstract0_free (pr->man_dcons, datad);

  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;
}

ucons_t*
build_const_11 (ucons_internal_t * pr, ucons_t * r,
                ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = saux = 1;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  checked_malloc (naux, ap_dim_t, sizeof (ap_dim_t), saux, return NULL;);
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1 && saux == 1, return r;);

  // m1 is the naux[0]
#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ucons_build \\forall y=l[%zu]-1. x%zu-d(y)>=l[n%zu]-2 && d(y) >= 1\n",
           nmain[0] - r->datadim, naux[0], nmain[0] - r->datadim);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif

  pattern_key_t *look = NULL;
  checked_malloc (look, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  look->type = get_pattern_type (pr, 1, 0, 1, pattern_1_lx_1);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  ap_abstract0_t * data = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 2, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 1;


  ap_dim_t ln = nmain[0] + r->segmentdim; // nmain contains r->datadim
  ap_dim_t dn = nmain[0];

  ap_lincons0_array_t arr;
  arr = ap_lincons0_array_make (3);
  /*
   *  y - l[n] + 1 == 0
   */
  arr.p[0].constyp = AP_CONS_EQ;
  arr.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arr.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ln, -1);
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, y1, 1);
  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 1);

  //  m1 - d(y) - l[n] + 2 >= 0 // m2 - d(y) - l[n] + 0 >= 0
  arr.p[1].constyp = AP_CONS_SUPEQ;
  arr.p[1].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arr.p[1].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy1, -1);
  ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, naux[0], 1);
  ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, ln, -1);
  ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);

  /* l[n] >= 2 */
  arr.p[2].constyp = AP_CONS_SUPEQ;
  arr.p[2].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arr.p[2].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0, ln, 1);
  ap_linexpr0_set_cst_scalar_int (arr.p[2].linexpr0, -2);

  data = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                          true, data, &arr);

  ap_lincons0_array_clear (&arr);

  /* meet with r->econs */
  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 2, 0);
  dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *aux = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  data = ap_abstract0_meet (pr->man_dcons, true, data, aux);


  look->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, look, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = look->type;
      //for (size_t i=0 ; i<(u_seg); i++)
      a->key.segments[0] = look->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, look);
    }
  else
    {
      if (a->dcons != NULL)
        ap_abstract0_free (pr->man_dcons, a->dcons);
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);
    }
  ap_abstract0_free (pr->man_dcons, data);
  data = NULL;

  free (look);
  look = NULL;

  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;

}

ucons_t*
build_const_21 (ucons_internal_t * pr, ucons_t * r,
                ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = saux = 1;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  checked_malloc (naux, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1 && saux == 1, return r;);

  printf ("\n@@@@ a == null ucons_build_21 \\forall y in [n%zu]. d(y) < x%zu\n",
          nmain[0] - r->datadim, naux[0]);
  fflush (stdout);

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_build \\forall y in [n%zu]. d(y) < x%zu\n",
           nmain[0] - r->datadim, naux[0]);
  fflush (stdout);
#endif

  ap_abstract0_t * data = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 2, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 1;

  ap_lincons0_array_t arr;
  arr = ap_lincons0_array_make (1);
  /*
   *  v - 1 - d(y) >= 0
   */
  arr.p[0].constyp = AP_CONS_SUPEQ;
  arr.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arr.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, naux[0], 1);
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy1, -1);
  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);

  data = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                          true, data, &arr);

  ap_lincons0_array_clear (&arr);

  /* meet with r->econs */
  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 2, 0);
  dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *aux = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  data = ap_abstract0_meet (pr->man_dcons, true, data, aux);

  pattern_key_t *look = NULL;
  checked_malloc (look, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  look->type = get_pattern_type (pr, 1, 0, 1, pattern_1);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  look->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, look, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = look->type;
      a->key.segments[0] = look->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, look);
    }
  else
    {
      if (a->dcons != NULL)
        a->dcons = ap_abstract0_meet (pr->man_dcons, true, a->dcons, data);
      else
        a->dcons = ap_abstract0_copy (pr->man_dcons, data);
    }


  free (look);
  look = NULL;
  ap_abstract0_free (pr->man_dcons, data);

  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;
}

ucons_t*
build_const_22 (ucons_internal_t * pr, ucons_t * r,
                ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = saux = 1;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  checked_malloc (naux, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1 && saux == 1, return r;);

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_build \\forall y in [n%zu]. d(y) < x%zu\n",
           nmain[0] - r->datadim, naux[0]);
  fflush (stdout);
#endif

  ap_abstract0_t * data = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 4, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t y2 = r->datadim + 2 * r->segmentdim + 1;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 2;
  ap_dim_t dy2 = r->datadim + 2 * r->segmentdim + 3;

  ap_lincons0_array_t arr;
  arr = ap_lincons0_array_make (1);
  /*
   *  v - 1 - d(y2) >= 0
   */
  arr.p[0].constyp = AP_CONS_SUPEQ;
  arr.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arr.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, naux[0], 1);
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy2, -1);
  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);

  data = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                          true, data, &arr);

  ap_lincons0_array_clear (&arr);

  /* meet with r->econs */
  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 4, 0);
  dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[2] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[3] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *aux = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  data = ap_abstract0_meet (pr->man_dcons, true, data, aux);

  pattern_key_t *look = NULL;
  checked_malloc (look, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  look->type = get_pattern_type (pr, 1, 0, 2, pattern_1_2);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  look->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, look, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = look->type;
      a->key.segments[0] = look->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, look);
      printf ("\n@@@@ a == null ucons_build_22 \\forall y1,y2 in [n%zu]. d(y2) < x%zu\n",
              nmain[0] - r->datadim, naux[0]);
      fflush (stdout);
    }
  else
    {
      if (a->dcons != NULL)
        {
          a->dcons = ap_abstract0_meet (pr->man_dcons, true, a->dcons, data);
          printf ("\n@@@@ a->dcons != null ucons_build_22 \\forall y1,y2 in [n%zu]. d(y2) < x%zu\n",
                  nmain[0] - r->datadim, naux[0]);
          fflush (stdout);

        }
      else
        {
          printf ("\n@@@@ a->dcons == null ucons_build_22 \\forall y1,y2 in [n%zu]. d(y2) < x%zu\n",
                  nmain[0] - r->datadim, naux[0]);
          fflush (stdout);
          a->dcons = ap_abstract0_copy (pr->man_dcons, data);
        }
    }


  free (look);
  look = NULL;
  ap_abstract0_free (pr->man_dcons, data);

  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;
}

ucons_t*
build_const_23 (ucons_internal_t * pr, ucons_t * r,
                ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = saux = 1;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  checked_malloc (naux, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1 && saux == 1, return r;);

  printf ("\n@@@@ a == null ucons_build_21 \\forall y in [n%zu]. d(y) < x%zu\n",
          nmain[0] - r->datadim, naux[0]);
  fflush (stdout);

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_build \\forall y in [n%zu]. d(y) < x%zu\n",
           nmain[0] - r->datadim, naux[0]);
  fflush (stdout);
#endif

  ap_abstract0_t * data = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 2, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 1;

  ap_lincons0_array_t arr;
  arr = ap_lincons0_array_make (1);
  /*
   *  v - 1 - d(y) >= 0
   */
  arr.p[0].constyp = AP_CONS_SUPEQ;
  arr.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
  arr.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, naux[0], 1);
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy1, -1);
  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0); // v - d(y) >= 0

  data = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                          true, data, &arr);

  ap_lincons0_array_clear (&arr);

  /* meet with r->econs */
  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 2, 0);
  dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *aux = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  data = ap_abstract0_meet (pr->man_dcons, true, data, aux);

  pattern_key_t *look = NULL;
  checked_malloc (look, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  look->type = get_pattern_type (pr, 1, 0, 1, pattern_1);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  look->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, look, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = look->type;
      a->key.segments[0] = look->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, look);
    }
  else
    {
      if (a->dcons != NULL)
        a->dcons = ap_abstract0_meet (pr->man_dcons, true, a->dcons, data);
      else
        a->dcons = ap_abstract0_copy (pr->man_dcons, data);
    }


  free (look);
  look = NULL;
  ap_abstract0_free (pr->man_dcons, data);

  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;
}

ucons_t*
build_const_24 (ucons_internal_t * pr, ucons_t * r,
                ap_linexpr0_t *lexpr)
{

  ap_dim_t *nmain; // main dimensions for the pattern
  ap_dim_t *naux; // aux dimensions for the pattern
  size_t smain, saux;

  smain = saux = 1;
  checked_malloc (nmain, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  checked_malloc (naux, ap_dim_t, sizeof (ap_dim_t), smain, return NULL;);
  ucons_build_extract_dimensions (lexpr, nmain, &smain, naux, &saux);

  arg_assert (smain == 1 && saux == 1, return r;);

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_build \\forall y1, y2 in [n%zu]. d(y2) < x%zu\n",
           nmain[0] - r->datadim, naux[0]);
  fflush (stdout);
#endif

  ap_abstract0_t * data = ap_abstract0_top (pr->man_dcons, r->datadim + 2 * r->segmentdim + 4, 0);
  ap_dim_t y1 = r->datadim + 2 * r->segmentdim;
  ap_dim_t y2 = r->datadim + 2 * r->segmentdim + 1;
  ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 2;
  ap_dim_t dy2 = r->datadim + 2 * r->segmentdim + 3;

  ap_lincons0_array_t arr;
  arr = ap_lincons0_array_make (1);
  /*
   *  v - 1 - d(y2) >= 0
   */
  arr.p[0].constyp = AP_CONS_SUPEQ;
  arr.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arr.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, naux[0], 1);
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, dy2, -1);
  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0); //v-d(y2)>=0

  data = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                          true, data, &arr);

  ap_lincons0_array_clear (&arr);

  /* meet with r->econs */
  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 4, 0);
  dimadd.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[2] = r->datadim + 2 * r->segmentdim;
  dimadd.dim[3] = r->datadim + 2 * r->segmentdim;

  ap_abstract0_t *aux = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimadd, false);

  ap_dimchange_clear (&dimadd);

  data = ap_abstract0_meet (pr->man_dcons, true, data, aux);

  pattern_key_t *look = NULL;
  checked_malloc (look, pattern_key_t, 1,
                  sizeof (pattern_key_t) + 1 * sizeof (size_t), return NULL;);
  look->type = get_pattern_type (pr, 1, 0, 2, pattern_1_2);
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  look->segments[0] = nmain[0] - r->datadim;
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, look, keylen, a);
  if (!a)
    {
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+(1) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(1) * sizeof (size_t));
      a->key.type = look->type;
      a->key.segments[0] = look->segments[0];
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      r = add_pattern_n2p (pr, r, look);
      printf ("\n@@@@ a == null ucons_build_22 \\forall y1,y2 in [n%zu]. d(y2) < x%zu\n",
              nmain[0] - r->datadim, naux[0]);
      fflush (stdout);
    }
  else
    {
      if (a->dcons != NULL)
        {
          a->dcons = ap_abstract0_meet (pr->man_dcons, true, a->dcons, data);
          printf ("\n@@@@ a->dcons != null ucons_build_22 \\forall y1,y2 in [n%zu]. d(y2) < x%zu\n",
                  nmain[0] - r->datadim, naux[0]);
          fflush (stdout);

        }
      else
        {
          printf ("\n@@@@ a->dcons == null ucons_build_22 \\forall y1,y2 in [n%zu]. d(y2) < x%zu\n",
                  nmain[0] - r->datadim, naux[0]);
          fflush (stdout);
          a->dcons = ap_abstract0_copy (pr->man_dcons, data);
        }
    }


  free (look);
  look = NULL;
  ap_abstract0_free (pr->man_dcons, data);

  if (nmain) free (nmain);
  if (naux) free (naux);

  return r;
}
