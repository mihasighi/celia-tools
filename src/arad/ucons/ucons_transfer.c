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
 * Assignment, substitution, guard transfer functions.
 */

#include <assert.h>
#include "ucons.h"
#include "ucons_internal.h"
#include "apron2shape.h"
#include "ap_generic.h"
#include "ap_linexpr0.h"
#include "ap_pcons0.h"
#include "shad.h"


/* ============================================================ */
/* Meet constraints and Join generators */
/* ============================================================ */

/* Internal functions used by ucons_meet_lincons.
 */
void ucons_meet_mset_build (ucons_internal_t * pr,
                            size_t datadim, size_t segmdim,
                            ap_lincons0_t * lcons);

ap_linexpr0_t* ucons_meet_ucons_build (ucons_internal_t * pr,
                                       size_t datadim, size_t segmdim,
                                       ap_lincons0_t * lcons);

ap_dim_t* ucons_meet_equals_build (ucons_internal_t * pr,
                                   size_t datadim, size_t segmdim,
                                   ap_lincons0_t * lcons,
                                   size_t* sz);
ucons_t* ucons_meet_equals_apply (ucons_internal_t *pr,
                                  ucons_t *a,
                                  ap_dim_t *useg, size_t size_useg);

ap_lincons0_t ucons_meet_econs_build (ucons_internal_t * pr,
                                      size_t datadim, size_t segmdim,
                                      ap_lincons0_t * lcons);
ucons_t* ucons_meet_econs_apply (ucons_internal_t * pr,
                                 ucons_t * a, ap_lincons0_t * lcons);

ucons_t* ucons_meet_formula (ucons_internal_t* pr, ucons_t* a,
                             sh_formula_t* f, size_t disj);

/**
 * Meet with a linear constraint encoding different types of constraints.
 * (Used also by ucons_of_lincons_array.)
 *
 * @param lcons  the lcons->scalar encodes the kind of constraint.
 * @return       the abstract value from met of a with the constraint
 */
ucons_t *
ucons_meet_lincons (ucons_internal_t * pr, bool destructive,
                    ucons_t * a, ap_lincons0_t * lcons)
{

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_lincons: with lincons=(");
  if (lcons) ap_lincons0_fprint (stdout, lcons, NULL);
  fprintf (stdout, ") \n and offset (lcons->scalar): ");
  ap_scalar_fprint (stdout, lcons->scalar);
  fprintf (stdout, "\n on a=(");
  ucons_fdump (stdout, pr->man, a);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif
  ap_constyp_t op = lcons->constyp;
  ucons_t* r = NULL;
  offset_t kind = OFFSET_NONE; /* unknown */
  bool eqstruct = false; /* is a structural equality? */

  if (lcons->scalar)
    {
      if (!ap_scalar_cmp_int (lcons->scalar, a->datadim + OFFSET_DATA))
        kind = OFFSET_DATA; /* data */
      else if (!ap_scalar_cmp_int (lcons->scalar, OFFSET_LEN))
        kind = OFFSET_LEN; /* len */
      else if (!ap_scalar_cmp_int (lcons->scalar, OFFSET_UCONS))
        kind = OFFSET_UCONS;
      else if (!ap_scalar_cmp_int (lcons->scalar, OFFSET_MSET))
        kind = OFFSET_MSET; /* mset */
      else if (!ap_scalar_cmp_int (lcons->scalar, OFFSET_SL3))
        kind = OFFSET_SL3; /* SL3 */
      else
        if (op == AP_CONS_EQ || op == AP_CONS_EQMOD)
        {
          eqstruct = true;
          kind = OFFSET_NONE;
        }
      else
        assert (0);
    }

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_lincons: kind=%d\n", kind);
  fflush (stdout);
#endif

  // dispatch treatment depending of kind and eqstruct
  if (kind == OFFSET_MSET)
    {
      /* interpret the constraint */
      ucons_meet_mset_build (pr, a->datadim, a->segmentdim, lcons);

      /* do the meet */
      r = ucons_copy_internal (pr, a);
      /* meet done for all lcons in ucons_strengthen */

      /* free allocated memory */
      /* NONE */
    }
  else if (kind == OFFSET_UCONS)
    {
      /* interpret the constraint */
      ap_linexpr0_t* lexpr = ucons_meet_ucons_build (pr, a->datadim, a->segmentdim,
                                                     lcons);
      // Warning: ap_coeff_equal_int is not correct, use scalar
      ap_coeff_t* coeff = ap_linexpr0_cstref (lcons->linexpr0);
      ap_scalar_t* code = NULL;
      if (coeff && coeff->discr == AP_COEFF_SCALAR)
        code = coeff->val.scalar;

      /* do the meet */
      r = ucons_copy_internal (pr, a); /* ucons_alloc_top is not working here!*/
      r = build_constraint (pr, r, lexpr, code);

      /* free allocated memory */
      ap_linexpr0_free (lexpr);
      // DO NOT free coeff which is inside lcons!
    }
  else if (kind == OFFSET_SL3)
    {
      /* interpret the constraint */
      ap_coeff_t* coeff = ap_linexpr0_cstref (lcons->linexpr0);
      double code;
      ap_double_set_scalar (&code, coeff->val.scalar, GMP_RNDN);

      /* do the meet */
      r = ucons_copy_internal (pr, a);
      r = ucons_meet_formula (pr, r, sh_crt, (size_t) code);

      /* free allocated memory */
      /* NONE */
      // DO NOT free coeff which is inside lcons!
    }
  else
    {
      size_t size_useg = 0;
      ap_dim_t* useg = NULL;
      if (eqstruct)
        {
          /* interpret the constraint */
          useg = ucons_meet_equals_build (pr, a->datadim, a->segmentdim,
                                          lcons, &size_useg);
        }
      /* do the meet */
      if (size_useg >= 1)
        {
          r = ucons_meet_equals_apply (pr, a, useg, size_useg);
          /* free allocated memory */
          free (useg);
        }
      else
        {
          // data, len, or other constraints
          /* interpret the constraint */
          ap_lincons0_t econs = ucons_meet_econs_build (pr,
                                                        a->datadim, a->segmentdim,
                                                        lcons);
          /* do the meet */
          r = ucons_copy_internal (pr, a);
          r = ucons_meet_econs_apply (pr, r, &econs);

          /* free allocated memory */
          ap_lincons0_clear (&econs);
        }
    }

  if (destructive)
    ucons_free_internal (pr, a);
  return r;

}

/**
 * Variables used to collect the full set of multiset constraints
 * used in ucons_meet_lincons_array to call
 * the reduction operator ucons_strengthen.
 */

/**
 *  Nodes to eliminate.
 *  MS: it seems to be filled but not used...
 */
bool * mset_nodes = NULL;

/**
 * Array of mset inclusion between mset_nodes and the other nodes
 * eq_mset_nodes[i][j] == 1 <=> mset(i) included in mset(j)
 */
size_t ** eq_mset_nodes = NULL;

/**
 * Method which collect informations above from a mset constraint.
 * Three types of constraints are handled:
 *          n1 != 0               ==> node n1 is involved!
 *          n1 != 0 n2 != 0 n3==0 ==> mset(n1) == mset(n2)
 * 	    n1 != 0 n2 != 0 n3!=0 ==> mset(n1)+mset(n2)=mstl(n3)
 */
void
ucons_meet_mset_build (ucons_internal_t * pr,
                       size_t datadim, size_t segmdim,
                       ap_lincons0_t * lcons)
{

  assert (lcons && !ap_scalar_cmp_int (lcons->scalar, OFFSET_MSET));

  ap_dim_t n1, n2, n3; // Encoding to know when is
  // an eq constraint and when is just a node to eliminate
  n1 = n2 = n3 = 0;

  size_t i, dim;
  ap_coeff_t *coeff = NULL;

  ap_linexpr0_ForeachLinterm (lcons->linexpr0, i, dim, coeff)
  {
    if (coeff && !ap_coeff_zero (coeff))
      {
        /* transform the coefficient into scalar
         * to types of constraints handled :
         *          n1 != 0 n2 != 0 n3==0 ==> mset(n1) == mset(n2)
         * 	    n1 != 0 n2 != 0 n3!=0 ==> mset(n1)+mset(n2)=mstl(n3)
         */
        ap_coeff_t * c1 = ap_coeff_alloc (AP_COEFF_SCALAR);
        ap_coeff_t * c2 = ap_coeff_alloc (AP_COEFF_SCALAR);
        ap_coeff_t * c3 = ap_coeff_alloc (AP_COEFF_SCALAR);
        ap_coeff_set_scalar_int (c1, 1);
        ap_coeff_set_scalar_int (c2, -1);
        ap_coeff_set_scalar_int (c3, -2);
        if (ap_coeff_equal (coeff, c1))
          {
            if (n1 == 0) n1 = dim - datadim;
            else n2 = dim - datadim;
          }
        if (ap_coeff_equal (coeff, c2))
          {
            n2 = dim - datadim;
          }
        if (ap_coeff_equal (coeff, c3))
          {
            n3 = dim - datadim;
          }
      }
  }
  coeff = NULL;

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_mset_build: constraint = (");
  if (lcons) ap_lincons0_fprint (stdout, lcons, NULL);
  fprintf (stdout, "\n) n1 = %zu n2 = %zu n3=%zu \n ", n1, n2, n3);
  fflush (stdout);
#endif
  if (n2 == 0)
    {
      if (mset_nodes == NULL)
        {
          checked_realloc (mset_nodes, bool,
                           (segmdim), sizeof (bool),;);
          for (i = 0; i < segmdim; i++)
            mset_nodes[i] = false;
        }
      mset_nodes[n1] = true;
      n1 = 0;
      n3 = 0;
    }
  else
    {
      if (eq_mset_nodes == NULL)
        {
          checked_realloc (eq_mset_nodes, size_t*,
                           (segmdim), sizeof (size_t*),;);
          for (i = 0; i < segmdim; i++)
            {
              checked_malloc (eq_mset_nodes[i], size_t,
                              (segmdim), sizeof (size_t),;);
              for (size_t j = 0; j < segmdim; j++)
                eq_mset_nodes[i][j] = 0;
              eq_mset_nodes[i][i] = 1;
            }
        }
      if (n3 == 0)
        {
          /*  ms(n1) == ms(n2) */
          eq_mset_nodes[n1][n2] = 1;
          eq_mset_nodes[n2][n1] = 1;
          n1 = 0;
          n2 = 0;
        }
      else
        {

          /*  ms(n1) \cup ms(n2) == ms(n3) */
          eq_mset_nodes[n1][n3] = 2;
          eq_mset_nodes[n2][n3] = 2;
          n1 = 0;
          n2 = 0;
          n3 = 0;
        }
    }
}

/**
 * Buils an abstract value from its encoding into a linear constraint.
 * @see build_constraint
 */
ap_linexpr0_t *
ucons_meet_ucons_build (ucons_internal_t * pr,
                        size_t datadim, size_t segmdim,
                        ap_lincons0_t * lcons)
{

  assert (lcons && !ap_scalar_cmp_int (lcons->scalar, OFFSET_UCONS));

  ap_linexpr0_t *lexpr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, datadim + 2 * segmdim);
  size_t i, dim;
  ap_coeff_t *coeff = NULL;

  ap_linexpr0_ForeachLinterm (lcons->linexpr0, i, dim, coeff)
  {
    if (coeff && !ap_coeff_zero (coeff))
      {
        if (dim < datadim)
          ap_coeff_set (&lexpr->p.coeff[dim], coeff);
        else // dim >= datadim
          ap_coeff_set (&lexpr->p.coeff[dim], coeff);
      }
  }
  coeff = NULL;

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_ucons_build: lexpr=(");
  if (lexpr) ap_linexpr0_fprint (stdout, lexpr, NULL);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  return lexpr;
}

/**
 * Meet with a structural equality of segments.
 * Nodes starting the segments are given in the lcons.
 */
ap_dim_t*
ucons_meet_equals_build (ucons_internal_t * pr,
                         size_t datadim, size_t segmdim,
                         ap_lincons0_t * lcons,
                         size_t* size_useg)
{
  assert (lcons &&
          ((lcons->constyp == AP_CONS_EQ) || (lcons->constyp == AP_CONS_EQMOD)));

  ap_dim_t *useg = NULL; // for EQ_MOD: the segments with equality constraint

  size_t i, dim;
  ap_coeff_t *coeff = NULL;

  ap_linexpr0_ForeachLinterm (lcons->linexpr0, i, dim, coeff)
  {
    if (coeff && !ap_coeff_zero (coeff) && dim >= datadim)
      {
        checked_realloc (useg, ap_dim_t,
                         ((*size_useg) + 1), sizeof (ap_dim_t), return NULL;);
        useg[*size_useg] = dim - datadim; //+1 PB
        (*size_useg) += 1;
      }
  }

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_equals_build: between ( ");
  for (i = 0; i < (*size_useg); i++)
    fprintf (stdout, "useg[%zu]=%zu, ", i, useg[i]);
  fprintf (stdout, " )\n");
  fflush (stdout);
#endif

  return useg;
}

/**
 * Build an existential constraint from lcons.
 */
ap_lincons0_t
ucons_meet_econs_build (ucons_internal_t * pr,
                        size_t datadim, size_t segmdim,
                        ap_lincons0_t * lcons)
{
  assert (lcons && lcons->scalar);
  offset_t kind = OFFSET_NONE;
  if (!ap_scalar_cmp_int (lcons->scalar, datadim + OFFSET_DATA))
    kind = OFFSET_DATA; /* data */
  else if (!ap_scalar_cmp_int (lcons->scalar, OFFSET_LEN))
    kind = OFFSET_LEN; /* len */
  else if (lcons->constyp == AP_CONS_EQ || lcons->constyp == AP_CONS_EQMOD)
    kind = OFFSET_DATA; /* eqstruct on data */
  else
    assert (0);

  // size of the existential constraint = datadim + segmdim (d(n)) + segmdim (l(n))
  ap_linexpr0_t *lexpr = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                            datadim + 2 * segmdim);

  size_t i, dim;
  ap_coeff_t *coeff = NULL;

  ap_linexpr0_ForeachLinterm (lcons->linexpr0, i, dim, coeff)
  {
    if (coeff && !ap_coeff_zero (coeff))
      {
        if (dim < datadim)
          ap_coeff_set (&lexpr->p.coeff[dim], coeff);
        else // dim >= a->datadim
          {
            if (kind == OFFSET_DATA)
              ap_coeff_set (&lexpr->p.coeff[dim], coeff);
            else
              ap_coeff_set (&lexpr->p.coeff[dim + segmdim], coeff);
          }
      }
  }
  coeff = NULL;

  /* coeff = ap_coeff_alloc (AP_COEFF_SCALAR);
   * ap_linexpr0_get_cst (coeff, lcons->linexpr0);
   */
  coeff = ap_linexpr0_cstref (lcons->linexpr0);
  ap_linexpr0_set_cst (lexpr, coeff);

  ap_constyp_t op = lcons->constyp;
  /* normalize the comparison to >= */
  if (op == AP_CONS_SUP)
    {
      ap_coeff_t *oldcoeff = ap_coeff_alloc (AP_COEFF_SCALAR);
      double oldd;
      ap_linexpr0_get_cst (oldcoeff, lexpr);
      ap_double_set_scalar (&oldd, oldcoeff->val.scalar, GMP_RNDN);
      //       fprintf(stdout,"!!***!! scalar old %f \n", oldd);
      oldd -= 1;
      //        fprintf(stdout,"!!***!! scalar new %f \n", oldd);
      ap_linexpr0_set_cst_scalar_double (lexpr, oldd);
      ap_coeff_free (oldcoeff);
      op = AP_CONS_SUPEQ;
    }

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_econs_build: lexpr=(");
  if (lexpr) ap_linexpr0_fprint (stdout, lexpr, NULL);
  fprintf (stdout, ") op=%d \n", op);
  fflush (stdout);
#endif

  return ap_lincons0_make (op, lexpr, NULL);
}

/**
 * Meet a with an existential constraint GIVEN (decoded) in lcons.
 * Works directly on a.
 */
ucons_t *
ucons_meet_econs_apply (ucons_internal_t * pr,
                        ucons_t * a, ap_lincons0_t * lcons)
{
  assert (lcons && !lcons->scalar);

  ucons_t * r = a;
  ap_lincons0_array_t arr = ap_lincons0_array_make (1);
  arr.p[0] = ap_lincons0_copy (lcons);
  ap_linexpr0_t* lexpr = arr.p[0].linexpr0;

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_econs_apply: meet with classic lincons = (");
  ap_lincons0_array_fprint (stdout, &arr, NULL);
  fprintf (stdout, "\n) and r = (");
  ucons_fdump (stdout, pr->man, r);
  fprintf (stdout, ")\n");
  fprintf (stdout, "\n");
#endif

  r->econs =
          ap_abstract0_meet_lincons_array (pr->man_dcons, false, a->econs, &arr);
  /* DO NOT: clear arr here because lexpr used afterthat */

  size_t i, u_seg, e_seg, nr_y;
  pattern_t * s, *ra, *rt;
  unsigned keylen;

  for (s = a->udcons; s != NULL; s = s->hh.next)
    {
      nr_y = pr->PI[s->key.type].nr_y;
      /*MS: old version which seems to me not very clean to modify lexpr
      ap_linexpr0_realloc (lexpr, a->datadim + 2 * a->segmentdim + 2 * nr_y);
      ap_linexpr0_t * lexprr = ap_linexpr0_copy (lexpr);
       */
      ap_linexpr0_t * lexprr = ap_linexpr0_copy (lexpr);
      ap_linexpr0_realloc (lexprr, a->datadim + 2 * a->segmentdim + 2 * nr_y);

      ap_lincons0_array_t arrr = ap_lincons0_array_make (1);
      arrr.p[0] = ap_lincons0_make (lcons->constyp, lexprr, NULL);

      u_seg = pr->PI[s->key.type].u_seg;
      e_seg = pr->PI[s->key.type].e_seg;

      keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);
      HASH_FIND (hh, r->udcons, &s->key, keylen, rt);

      if (rt)
        {
          rt->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, false,
                                                       s->dcons, &arrr);
        }
      else
        {
          checked_malloc (ra, pattern_t, 1,
                          sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t),
                          return NULL;);
          memset (ra, 0, sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t));
          //		ra->dcons =(s->dcons) ?
          //				ap_abstract0_copy (pr->man_dcons, s->dcons) : NULL;
          ra->dcons = NULL;
          ra->key.type = s->key.type;

          for (i = 0; i < (u_seg + e_seg); i++)
            ra->key.segments[i] = s->key.segments[i];
          ra->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, false,
                                                       s->dcons, &arrr);
          HASH_ADD (hh, r->udcons, key, keylen, ra);
          add_pattern_n2p (pr, r, &ra->key);
        }
      ap_lincons0_array_clear (&arrr); // free also lexprr
    }
#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_econs_apply: builds arr=(");
  ap_lincons0_array_fprint (stdout, &arr, NULL);
  fprintf (stdout, ")\n returns r=(");
  ucons_fdump (stdout, pr->man, r);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif
  ap_lincons0_array_clear (&arr);

  return r;
}

ucons_t*
ucons_meet_equals_apply (ucons_internal_t *pr,
                         ucons_t *a,
                         ap_dim_t *useg, size_t size_useg)
{

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_equals_apply: meet with eq lincons between = ( ");
  for (size_t kk = 0; kk < size_useg; kk++)
    fprintf (stdout, "useg[%zu]=%zu, ", kk, useg[kk]);
  fprintf (stdout, " )\n");
  fflush (stdout);
#endif

  ucons_t* r = ucons_copy_internal (pr, a);
  ap_lincons0_array_t arr_ex = ap_lincons0_array_make ((size_useg * (size_useg - 1)));

  for (size_t i = 0; i < size_useg; i++)
    for (size_t j = i + 1; j < size_useg; j++)
      {

        size_t li = r->datadim + useg[i] + r->segmentdim;
        size_t lj = r->datadim + useg[j] + r->segmentdim;
        size_t di = r->datadim + useg[i];
        size_t dj = r->datadim + useg[j];

        /* data + len transfer */
        // li-lj == 0
        arr_ex.p[0].constyp = AP_CONS_EQ;
        arr_ex.p[0].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                  r->datadim + 2 * r->segmentdim);
        arr_ex.p[0].scalar = NULL;
        ap_linexpr0_set_coeff_scalar_int (arr_ex.p[0].linexpr0, li, 1);
        ap_linexpr0_set_cst_scalar_int (arr_ex.p[0].linexpr0, 0);
        ap_linexpr0_set_coeff_scalar_int (arr_ex.p[0].linexpr0, lj, -1);

        // di-dj == 0
        arr_ex.p[1].constyp = AP_CONS_EQ;
        arr_ex.p[1].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                  r->datadim + 2 * r->segmentdim);
        arr_ex.p[1].scalar = NULL;
        ap_linexpr0_set_coeff_scalar_int (arr_ex.p[1].linexpr0, di, 1);
        ap_linexpr0_set_cst_scalar_int (arr_ex.p[1].linexpr0, 0);
        ap_linexpr0_set_coeff_scalar_int (arr_ex.p[1].linexpr0, dj, -1);


        r->econs = ap_abstract0_meet_lincons_array (pr->man_dcons, false,
                                                    a->econs, &arr_ex);


        size_t u_seg, e_seg, nr_y;
        pattern_t * s, *ra, *rt;
        ra = NULL;
        unsigned keylen;

        for (s = a->udcons; s != NULL; s = s->hh.next)
          {
            nr_y = pr->PI[s->key.type].nr_y;
            ap_linexpr0_realloc (arr_ex.p[0].linexpr0,
                                 a->datadim + 2 * a->segmentdim + 2 * nr_y);
            ap_linexpr0_realloc (arr_ex.p[1].linexpr0,
                                 a->datadim + 2 * a->segmentdim + 2 * nr_y);

            ap_linexpr0_t * lexprr0 = ap_linexpr0_copy (arr_ex.p[0].linexpr0);
            ap_linexpr0_t * lexprr1 = ap_linexpr0_copy (arr_ex.p[1].linexpr0);

            ap_lincons0_array_t arrr = ap_lincons0_array_make (2);
            arrr.p[0] = ap_lincons0_make (arr_ex.p[0].constyp, lexprr0, NULL);
            arrr.p[1] = ap_lincons0_make (arr_ex.p[1].constyp, lexprr1, NULL);

            u_seg = pr->PI[s->key.type].u_seg;
            e_seg = pr->PI[s->key.type].e_seg;

            keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);
            HASH_FIND (hh, r->udcons, &s->key, keylen, rt);

            if (rt)
              {
                rt->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, false,
                                                             s->dcons, &arrr);
              }
            else
              {
                //checked_malloc(ra,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL;);
                ra = (pattern_t*) malloc (sizeof (pattern_t)+(u_seg + e_seg) *
                                          sizeof (size_t));
                memset (ra, 0,
                        sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t));
                //		ra->dcons =(s->dcons) ?
                //				ap_abstract0_copy (pr->man_dcons, s->dcons) : NULL;
                ra->dcons = NULL;
                ra->key.type = s->key.type;

                for (i = 0; i < (u_seg + e_seg); i++)
                  ra->key.segments[i] = s->key.segments[i];
                ra->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, false,
                                                             s->dcons, &arrr);
                HASH_ADD (hh, r->udcons, key, keylen, ra);
                add_pattern_n2p (pr, r, &ra->key);
              }
            ap_lincons0_array_clear (&arrr);
          }

#ifndef NDEBUG1
        fprintf (stdout, "@@@@ ucons_meet_lincons: with lincons=(");
        fprintf (stdout, ") builds arr=(");
        ap_lincons0_array_fprint (stdout, &arr_ex, NULL);
        fprintf (stdout, ") returns r=( ");
        ucons_fdump (stdout, pr->man, r);
        fprintf (stdout, ")\n");
        fflush (stdout);
#endif
        ap_lincons0_array_clear (&arr_ex);


        /* pattern transfer */

        pattern_key_set_t psi = r->n2p[useg[i]];
        pattern_key_set_t psj = r->n2p[useg[j]];

        size_t si = r->n2p[useg[i]].size;
        size_t sj = r->n2p[useg[j]].size;

        if ((si > 0) || (sj > 0))
          {
            /* transfer universals */
            for (size_t q = 0; q < si; q++)
              {
                pattern_key_t *pq = psi.p[q];
                u_seg = pr->PI[pq->type].u_seg;
                e_seg = pr->PI[pq->type].e_seg;
                keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);
                pattern_t *auxi, *auxj;
                HASH_FIND (hh, r->udcons, pq, keylen, auxi);

                pattern_key_t *lookj = NULL;
                //	checked_malloc(lookj, pattern_key_t, sizeof(pattern_key_t) + (u_seg+e_seg)*sizeof(size_t), 1, return NULL;);
                lookj = (pattern_key_t*) malloc (sizeof (pattern_key_t)+(u_seg + e_seg) *
                                                 sizeof (size_t));
                memset (lookj, 0,
                        sizeof (pattern_key_t) + (u_seg + e_seg) * sizeof (size_t));
                lookj->type = pq->type;
                /* copy + modif + sort */
                bool sort_eseg = false;
                for (size_t kk = 0; kk < u_seg + e_seg; kk++)
                  {
                    if (pq->segments[kk] == useg[i])
                      {
                        if (kk >= u_seg) sort_eseg = true;
                        lookj->segments[kk] = useg[j];
                      }
                  }
                if (sort_eseg)
                  {
                    /* sorting existing segment */
                    for (size_t ii = u_seg; ii < u_seg + e_seg; ii++)
                      {
                        size_t jj = u_seg;
                        while (jj != ii &&
                               lookj->segments[jj] <= lookj->segments[ii])
                          jj++;
                        if (jj < ii)
                          {
                            size_t d = lookj->segments[ii];
                            size_t tmp;
                            size_t kk;
                            for (size_t tmp = ii; tmp > jj; tmp--)
                              {
                                kk = tmp - 1;
                                lookj->segments[kk + 1] = lookj->segments[kk];
                              }
                            lookj->segments[jj] = d;
                          }
                      }
                  }
                else
                  if (u_seg == 2)
                  {
                    if (lookj->segments[0] > lookj->segments[1])
                      {
                        size_t jj = lookj->segments[0];
                        lookj->segments[0] = lookj->segments[1];
                        lookj->segments[0] = jj;
                      }
                  }
                //unsigned keylenj = (u_seg+e_seg)*sizeof(size_t) + sizeof(pattern_key_t);
                HASH_FIND (hh, r->udcons, lookj, keylen, auxj);

                if (auxi)
                  {
                    if (auxj)
                      {
                        ap_abstract0_t* dcons = ap_abstract0_meet (pr->man_dcons, true,
                                                                   auxi->dcons, auxj->dcons);
                        auxi->dcons = ap_abstract0_copy (pr->man_dcons, dcons);
                        auxj->dcons = ap_abstract0_copy (pr->man_dcons, dcons);
                        ap_abstract0_free (pr->man_dcons, dcons);
                      }
                    else
                      {
                        //								checked_malloc(auxj,pattern_t,1,
                        //										sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL;);
                        auxj = (pattern_t*) malloc (sizeof (pattern_t)+(u_seg + e_seg) *
                                                    sizeof (size_t));
                        memset (auxj, 0,
                                sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t));
                        auxj->key.type = lookj->type;
                        for (size_t i = 0; i < (u_seg + e_seg); i++)
                          auxj->key.segments[i] = lookj->segments[i];
                        auxj->dcons = ap_abstract0_copy (pr->man_dcons, auxi->dcons);

                        HASH_ADD (hh, r->udcons, key, keylen, auxj);
                        add_pattern_n2p (pr, r, lookj);
                      }
                  }

              }

            for (size_t q = 0; q < sj; q++)
              {
                pattern_key_t *pq = psj.p[q];
                size_t u_seg = pr->PI[pq->type].u_seg;
                size_t e_seg = pr->PI[pq->type].e_seg;
                keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);
                pattern_t *auxi, *auxj;
                HASH_FIND (hh, r->udcons, pq, keylen, auxj);

                pattern_key_t *looki = NULL;
                //checked_malloc(looki, pattern_key_t, sizeof(pattern_key_t) + (u_seg+e_seg)*sizeof(size_t), 1, return NULL;);
                looki = (pattern_key_t*) malloc (sizeof (pattern_key_t)+(u_seg + e_seg) *
                                                 sizeof (size_t));
                memset (looki, 0,
                        sizeof (pattern_key_t) + (u_seg + e_seg) * sizeof (size_t));
                looki->type = pq->type;
                /* copy + modif + sort */
                bool sort_eseg = false;
                for (size_t kk = 0; kk < u_seg + e_seg; kk++)
                  {
                    if (pq->segments[kk] == useg[j])
                      {
                        if (kk >= u_seg) sort_eseg = true;
                        looki->segments[kk] = useg[i];
                      }
                  }
                if (sort_eseg)
                  {
                    /*sorting exist segment */
                    for (size_t ii = u_seg; ii < u_seg + e_seg; ii++)
                      {
                        size_t jj = u_seg;
                        while (jj != ii && looki->segments[jj] <= looki->segments[ii])
                          jj++;
                        if (jj < ii)
                          {
                            size_t d = looki->segments[ii];
                            size_t tmp;
                            size_t kk;
                            for (size_t tmp = ii; tmp > jj; tmp--)
                              {
                                kk = tmp - 1;
                                looki->segments[kk + 1] = looki->segments[kk];
                              }
                            looki->segments[jj] = d;
                          }
                      }
                  }
                else
                  if (u_seg == 2)
                  {
                    if (looki->segments[0] > looki->segments[1])
                      {
                        size_t jj = looki->segments[0];
                        looki->segments[0] = looki->segments[1];
                        looki->segments[0] = jj;
                      }
                  }
                //unsigned keylenj = (u_seg+e_seg)*sizeof(size_t) + sizeof(pattern_key_t);
                HASH_FIND (hh, r->udcons, looki, keylen, auxi);

                if (auxj)
                  {
                    if (!auxi)
                      {
                        //								checked_malloc(auxi,pattern_t,1,
                        //										sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL;);
                        auxi = (pattern_t*) malloc (sizeof (pattern_t)+(u_seg + e_seg) *
                                                    sizeof (size_t));
                        memset (auxi, 0,
                                sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t));
                        auxi->key.type = looki->type;
                        for (size_t i = 0; i < (u_seg + e_seg); i++)
                          auxi->key.segments[i] = looki->segments[i];
                        auxi->dcons = ap_abstract0_copy (pr->man_dcons, auxj->dcons);

                        HASH_ADD (hh, r->udcons, key, keylen, auxi);
                        add_pattern_n2p (pr, r, looki);
                      }
                  }

              }
          }

        /* equality pattern constraint */
        //#if defined(UCONS_DCONS_OCT_P21) || defined(UCONS_DCONS_POLY_P21)
        if (pr->active_patterns[2])
          {
#ifndef NDEBUG2
            fprintf (stdout, "@@@@ ucons_meet_equals_apply: meet with eq constraint\n");
            fflush (stdout);
#endif

            //	pattern_t *ra;
            if (ra)
              {
                free (ra);
                ra = NULL;
              }
            //checked_malloc(ra,pattern_t,1,sizeof(pattern_t)+(2)*sizeof(size_t),return NULL;);
            ra = (pattern_t*) malloc (sizeof (pattern_t)+(2) * sizeof (size_t));
            memset (ra, 0, sizeof (pattern_t)+(2 * sizeof (size_t)));

            ra->dcons = NULL;
            ra->key.type = 1; // y1=y2

            if (useg[i] < useg[j])
              {
                ra->key.segments[0] = useg[i];
                ra->key.segments[1] = useg[j];
              }
            else
              {

                ra->key.segments[0] = useg[j];
                ra->key.segments[1] = useg[i];
              }
            /* add dims to r->econs */
            ap_dimchange_t* dimchange = ap_dimchange_alloc (4, 0);
            dimchange->dim[0] = r->datadim + 2 * r->segmentdim;
            dimchange->dim[1] = r->datadim + 2 * r->segmentdim;
            dimchange->dim[2] = r->datadim + 2 * r->segmentdim;
            dimchange->dim[3] = r->datadim + 2 * r->segmentdim;
            ap_dim_t ly1 = r->datadim + 2 * r->segmentdim;
            ap_dim_t ly2 = r->datadim + 2 * r->segmentdim + 1;
            ap_dim_t dy1 = r->datadim + 2 * r->segmentdim + 2;
            ap_dim_t dy2 = r->datadim + 2 * r->segmentdim + 3;

            ra->dcons = ap_abstract0_add_dimensions (pr->man_dcons, false,
                                                     r->econs, dimchange, false);
            ap_dimchange_free (dimchange);

            ap_lincons0_array_t arr = ap_lincons0_array_make (2);

            // ly1-ly2 == 0
            arr.p[0].constyp = AP_CONS_EQ;
            arr.p[0].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                   r->datadim +
                                                   2 * r->segmentdim + 4);
            arr.p[0].scalar = NULL;
            ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ly1, 1);
            ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);
            ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ly2, -1);
            //dy1-dy2=0
            arr.p[1].constyp = AP_CONS_EQ;
            arr.p[1].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                   r->datadim +
                                                   2 * r->segmentdim + 4);
            arr.p[1].scalar = NULL;
            ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy2, 1);
            ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);
            ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy1, -1);

            ra->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true,
                                                         ra->dcons, &arr);
            ap_lincons0_array_clear (&arr);

            keylen = 2 * sizeof (size_t) + sizeof (pattern_key_t);
            HASH_ADD (hh, r->udcons, key, keylen, ra);
            add_pattern_n2p (pr, r, &ra->key);
            //#endif
          }
      }
  return r;
}

ucons_t *
ucons_meet_lincons_array (ap_manager_t * man,
                          bool destructive, ucons_t * a,
                          ap_lincons0_array_t * array)
{
  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_MEET_LINCONS_ARRAY, 0);

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_lincons_array: ENTER\n");
  fflush (stdout);
#endif

  if (!a)
    return NULL;
  if (!array || array->size == 0 || !array->p)
    return (destructive) ? a : ucons_copy (man, a);
  size_t i;

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_lincons_array: with array=(");
  ap_lincons0_array_fprint (stdout, array, NULL);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif


#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_lincons_array: with array[0]=(");
  ap_lincons0_fprint (stdout, &array->p[0], NULL);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  ucons_t * r = ucons_meet_lincons (pr, false, a, &array->p[0]);

  for (i = 1; i < array->size; i++)
    {
#ifndef NDEBUG1
      fprintf (stdout, "@@@@ ucons_meet_lincons_array: with array[%zu]=(", i);
      ap_lincons0_fprint (stdout, &array->p[i], NULL);
      fprintf (stdout, ")\n");
      fflush (stdout);
#endif
      r = ucons_meet_lincons (pr, true, r, &array->p[i]);
    }

  if ((eq_mset_nodes != NULL))
    {
      /* if the array corresponds to an intersection with a multiset constraint */

      r = ucons_strengthen (pr, a, eq_mset_nodes);
      free (mset_nodes);
      free (eq_mset_nodes);
      mset_nodes = NULL;
      eq_mset_nodes = NULL;
    }
  if (mset_nodes != NULL)
    {
      free (mset_nodes);
      mset_nodes = NULL;
    }

  if (destructive)
    ucons_free_internal (pr, a);
  return r;

}

/**
 * Translate a data dimension to separate universally quantified vars.
 * @returns either < datadim for data var or datadim + 2segmdim + i
 */
ap_dim_t
translate_dimension (ap_dim_t dim,
                     size_t datadim, size_t segmdim,
                     size_t ydim, ap_dim_t* y)
{
  ap_dim_t ndim;
  size_t i, less_y;
  less_y = 0;
  for (i = 0; i < ydim; i++)
    {
#ifndef NDEBUG1
      fprintf (stdout, "@@@@ translate_dimension: y[%zu]=%zu\n", i, y[i]);
      fflush (stdout);
#endif
      if (y[i] == dim) return (datadim + 2 * segmdim + i);
      if (y[i] < dim) less_y++;
    }
  // here is clear that it shall be less than datadim
  ndim = dim - less_y;
  return (ndim < datadim) ? ndim : 0;
}

/**
 * Buils a linear constraint supported by the ucons domain
 * from a dataformula.
 * Called in both existential and universal part.
 * In the universal part, datadim is computed excluding ydim!
 */
ap_lincons0_t
ucons_meet_formula_build_data (ucons_internal_t* pr,
                               size_t datadim, size_t segmdim,
                               size_t ydim, ap_dim_t* y,
                               sh_dataform_t* df)
{
  assert (df);
#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_formula_build_data: (datadim=%zu,segmdim=%zu,ydim=%zu)\n",
           datadim, segmdim, ydim);
  fflush (stdout);
#endif

  ap_linexpr0_t* lexpr = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                            datadim + 2 * segmdim + 2 * ydim);
  sh_linterm_t* ti = df->p;
  for (ti = df->p; ti != NULL; ti = ti->next)
    {
      // compute the dimension corresponding to this term
      ap_dim_t dim;
      switch (ti->t.funid)
        {
        case SH_F_SYMBOL: // may be intdim or ydim
          if (ydim == 0) dim = ti->t.p0; // ok dimension already computed
          else
            {
              // data variable or y
              dim = translate_dimension (ti->t.p0, datadim, segmdim, ydim, y);
              // returns < datadim or datadim + 2*segmdim + y
            }
#ifndef NDEBUG1
          fprintf (stdout, "@@@@ ucons_meet_formula_build_data: symbol(%zu)\n",
                   dim);
          fflush (stdout);
#endif
          break;
        case SH_F_LEN: // parameter is node (thus including datadim) but not nilNode!
          dim = segmdim + (ti->t.p0 + 1) - ydim;
          // TODO: ti->t.p0 is computed in the environment including datadim and thus ydim
#ifndef NDEBUG1
          fprintf (stdout, "@@@@ ucons_meet_formula_build_data: len(%zu) to %zu\n",
                   ti->t.p0, dim);
          fflush (stdout);
#endif
          break;
        case SH_F_DATA: // the second parameter is supposed to be 0!
          if (ydim == 0 && ti->t.p1 != AP_DIM_MAX)
            {
#ifndef NDEBUG1
              fprintf (stdout, "@@@@ ucons_meet_formula_build_data: error for data(%d,%d!=AP_DIM_MAX)\n",
                       ti->t.p0, ti->t.p1);
              fflush (stdout);
#endif
              goto end_of_dataform;
            }
          if (ydim == 0 || ti->t.p1 == AP_DIM_MAX)
            {
              // data(n,0)
              dim = ti->t.p0 + 1 - ydim; // TODO: ti->t.p0 includes datadim?
#ifndef NDEBUG1
              fprintf (stdout, "@@@@ ucons_meet_formula_build_data: data(n%d,y%d)=data(n%zu)\n",
                       ti->t.p0, ti->t.p1, dim);
              fflush (stdout);
#endif
            }
          else // data(n,y)
            {
              // TODO: check that ti->t.p1 corresponds to y position in the list!
              dim = ydim + translate_dimension (ti->t.p1, datadim, segmdim, ydim, y);
#ifndef NDEBUG1
              fprintf (stdout, "@@@@ ucons_meet_formula_build_data: data(n%d,y%d)=data(y%zu)\n",
                       ti->t.p0, ti->t.p1, dim);
              fflush (stdout);
#endif
            }
          break;
        default: /* not dealt: SUM, MSET, break the "for" */
#ifndef NDEBUG1
          fprintf (stdout, "@@@@ ucons_meet_formula_build_data: constraint not dealt (LSUM|MSET)\n");
          fflush (stdout);
#endif
          goto end_of_dataform;
          break;
        }
      if (ti->coeff)
        ap_coeff_set_scalar_int (&lexpr->p.coeff[dim], ti->coeff);
    }
  int cst;
  ap_constyp_t op;
end_of_dataform:
  // avoid the strict constraints over int
  cst = df->cst;
  op = df->constyp;
#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_formula_build_data: op = %d\n", op);
  fflush (stdout);
#endif
  if (df->constyp == AP_CONS_SUP)
    {
      cst -= 1;
      op = AP_CONS_SUPEQ;
    }
  ap_linexpr0_set_cst_scalar_int (lexpr, cst);
  return ap_lincons0_make (op, lexpr, NULL);
}

/**
 * Build an array of constraints from the data formula.
 */
ap_lincons0_array_t
ucons_meet_formula_build_univ (ucons_internal_t* pr,
                               size_t datadim, size_t segmdim,
                               sh_univform_t* du)
{
#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_formula_build_univ: (datadim=%zu,segmdim=%zu)\n",
           datadim, segmdim);
  fflush (stdout);
#endif
  // initialize the array of constraints with the sum of:
  // - size of du->length_data
  // - size of 2*du->length_y (1 <= y <= len(n)- (pos y) - 1)
  // - 1 additional constraint except for the Gall guard
  size_t arr_size = du->length_data + 2 * du->length_y + 1;
  ap_lincons0_array_t arr = ap_lincons0_array_make (arr_size);

  // fill with the constraints in du->data
  size_t i, j;
  for (i = 0; i < du->length_data; i++)
    {
      sh_dataform_t df = du->data[i];
      arr.p[i] = ucons_meet_formula_build_data (pr, datadim, segmdim,
                                                du->guard.length_y,
                                                du->guard.y,
                                                &df);
    }
  // i == du->length_data
  // fill with the constraints corresponding to conditions on du->y
  for (j = 0; j < du->length_y; j++)
    {
      // 2 <= len(n[j]) && y[j] <= len(n[j])-1
      ap_dim_t y = datadim + 2 * segmdim + j;
      ap_dim_t len_n = segmdim + du->guard.n[j] + 1 - du->length_y;
      // n[j] includes datadim  and length_y but not nilNode
#ifndef NDEBUG1
      fprintf (stdout, "@@@@ ucons_meet_formula_build_univ: 1 <= dim%zu(%zu) < dim%zu(%zu)\n",
               y, du->guard.y[j], len_n, du->guard.n[j]);
      fflush (stdout);
#endif
      // 2 <= len(n[j])
      int cst = (du->guard.guardtyp == SH_G_SUCC2) ? 3 : 2;
      arr.p[i].constyp = AP_CONS_SUPEQ;
      arr.p[i].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                             datadim + 2 * segmdim + 2 * du->guard.length_y);
      arr.p[i].scalar = NULL;
      ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, len_n, 1);
      ap_linexpr0_set_cst_scalar_int (arr.p[i].linexpr0, -cst);
      i++;
      // y[j] <= len(n[j])-1
      arr.p[i].constyp = AP_CONS_SUPEQ;
      arr.p[i].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                             datadim + 2 * segmdim + 2 * du->guard.length_y);
      arr.p[i].scalar = NULL;
      ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, len_n, 1);
      ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, y, -1);
      ap_linexpr0_set_cst_scalar_int (arr.p[i].linexpr0, -1);
      i++;
    }
  // fill with the additional constraint
  switch (du->guard.guardtyp)
    {
      /* Warning: this shall be included but does not work well with post
    case SH_G_FST:
      {
        ap_dim_t y = datadim + 2 * segmdim;
        // 1 = y
        arr.p[i].constyp = AP_CONS_EQ;
        arr.p[i].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                               datadim + 2 * segmdim + 2 * du->guard.length_y);
        arr.p[i].scalar = NULL;
        ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, y, 1);
        ap_linexpr0_set_cst_scalar_int (arr.p[i].linexpr0, -1);
        break;
      }
       */
    case SH_G_LST:
      {
        ap_dim_t y = datadim + 2 * segmdim;
        ap_dim_t len_n = segmdim + du->guard.n[0] + 1 - du->length_y;
        // du->guard.n[0] includes datadim and ydim but not nilNode
#ifndef NDEBUG1
        fprintf (stdout, "@@@@ ucons_meet_formula_build_univ: dim%zu(%zu) = dim%zu(%zu)\n",
                 y, du->guard.y[0], len_n, du->guard.n[0]);
        fflush (stdout);
#endif
        // y = len(n)-1
        arr.p[i].constyp = AP_CONS_EQ;
        arr.p[i].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                               datadim + 2 * segmdim + 2 * du->guard.length_y);
        arr.p[i].scalar = NULL;
        ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, len_n, 1);
        ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, y, -1);
        ap_linexpr0_set_cst_scalar_int (arr.p[i].linexpr0, -1);
        i++;
        break;
      }
    case SH_G_LE2:
      {
        ap_dim_t y1 = datadim + 2 * segmdim;
        ap_dim_t y2 = datadim + 2 * segmdim + 1;
        // -y1+y2 >= 0
        arr.p[i].constyp = AP_CONS_SUPEQ;
        arr.p[i].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                               datadim + 2 * segmdim + 2 * du->guard.length_y);
        arr.p[i].scalar = NULL;
        ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, y2, 1);
        ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, y1, -1);
        ap_linexpr0_set_cst_scalar_int (arr.p[i].linexpr0, 0);
        i++;
        break;
      }
    case SH_G_SUCC2:
      {
        ap_dim_t y1 = datadim + 2 * segmdim;
        ap_dim_t y2 = datadim + 2 * segmdim + 1;
        // -1-y1+y2 = 0
        arr.p[i].constyp = AP_CONS_EQ;
        arr.p[i].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                               datadim + 2 * segmdim + 2 * du->guard.length_y);
        arr.p[i].scalar = NULL;
        ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, y2, 1);
        ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, y1, -1);
        ap_linexpr0_set_cst_scalar_int (arr.p[i].linexpr0, -1);
        i++;
        break;
      }
    case SH_G_EQ2:
      {
        ap_dim_t y1 = datadim + 2 * segmdim;
        ap_dim_t y2 = datadim + 2 * segmdim + 1;
        // -y1+y2 = 0
        arr.p[i].constyp = AP_CONS_EQ;
        arr.p[i].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                               datadim + 2 * segmdim + 2 * du->guard.length_y);
        arr.p[i].scalar = NULL;
        ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, y2, 1);
        ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0, y1, -1);
        ap_linexpr0_set_cst_scalar_int (arr.p[i].linexpr0, 0);
        i++;
        break;
      }
    default:
      break;
    }
  ap_lincons0_array_resize (&arr, i);

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_formula_build_univ: constraints = (\n");
  ap_lincons0_array_fprint (stdout, &arr, NULL);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  return arr;
}

pattern_kind
translate_guard (sh_guardform_t* g, size_t* useg, size_t* nr_y)
{
  switch (g->guardtyp)
    {
    case SH_G_ALL:
      *useg = 1;
      *nr_y = 1;
      return pattern_1;
    case SH_G_LE2:
      *useg = 1;
      *nr_y = 2;
      return pattern_1_2;
    case SH_G_SUCC2:
      *useg = 1;
      *nr_y = 2;
      return pattern_succ_1_2;
    case SH_G_FST:
      *useg = 1;
      *nr_y = 1;
      return pattern_1_l1;
    case SH_G_LST:
      *useg = 1;
      *nr_y = 1;
      return pattern_1_lx_1;
    case SH_G_EQ2:
      *useg = 2;
      *nr_y = 2;
      return pattern_2_1;
    case SH_G_SL2:
      *useg = 2;
      *nr_y = 2;
      return pattern_2_1_mlx;
    case SH_G_SR2:
    default:
      *useg = 2;
      *nr_y = 2;
      return pattern_2_1_lx;
    }
}

/**
 * Meet the guard given by g in r with the array of constraints arr,
 * i.e., r->econs /\ r->ducons[g] /\ arr
 */
ucons_t*
ucons_meet_formula_apply_univ (ucons_internal_t* pr, ucons_t* r,
                               sh_guardform_t* g,
                               ap_lincons0_array_t* arr)
{
  if (!g)
    {
#ifndef NDEBUG1
      fprintf (stdout, "@@@@ ucons_meet_formula_apply_univ: guard null\n");
      fflush (stdout);
#endif
      return r;
    }

  // dimension size of r->ducons[g] and arr
  size_t dsize = r->datadim + 2 * r->segmentdim + g->length_y * 2;

  // Step 1: meet top /\ arr
  ap_abstract0_t * data = ap_abstract0_top (pr->man_dcons, dsize, 0);

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_formula_apply_univ: step 0 constraint ucons=(");
  ap_abstract0_fprint (stdout, pr->man_dcons, data, NULL);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  data = ap_abstract0_meet_lincons_array (pr->man_dcons, true, data, arr);

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_formula_apply_univ: step 1 constraint ucons=(");
  ap_abstract0_fprint (stdout, pr->man_dcons, data, NULL);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  // Step 2: meet with r->econs, add dimensions at the end of r->econs
  ap_dimchange_t dimadd;
  size_t szadd = g->length_y * 2;
  ap_dimchange_init (&dimadd, szadd, 0);
  dimadd.dim = (ap_dim_t *) malloc (szadd * sizeof (ap_dim_t));
  for (size_t i = 0; i < szadd; i++)
    dimadd.dim[i] = r->datadim + 2 * r->segmentdim;
  ap_abstract0_t *aux = ap_abstract0_add_dimensions (pr->man_dcons, false,
                                                     r->econs, &dimadd, false);
  ap_dimchange_clear (&dimadd);
  data = ap_abstract0_meet (pr->man_dcons, true, data, aux);
  ap_abstract0_free (pr->man_dcons, aux);

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_formula_apply_univ: step 2 constraint ucons=(");
  ap_abstract0_fprint (stdout, pr->man_dcons, data, NULL);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  // Step 3: meet with r->ducons[g]
  // - translate g in patern
  size_t u_seg, nr_y;
  pattern_kind p = translate_guard (g, &u_seg, &nr_y);
  // - build the pattern
  pattern_key_t *look = NULL;
  checked_malloc (look, pattern_key_t, 1,
                  sizeof (pattern_key_t) + (u_seg) * sizeof (size_t), return NULL;);
  look->type = get_pattern_type (pr, u_seg, 0, nr_y, p);
  unsigned keylen = (u_seg) * sizeof (size_t) + sizeof (pattern_key_t);
  // TODO: add e_seg to u_seg above iff guards with existential segments

  for (size_t i = 0; i < u_seg; i++)
    {
      // dim computed in environment with y, but not including nilNode
      look->segments[i] = g->n[i] + 1 - r->datadim - nr_y;
#ifndef NDEBUG1
      fprintf (stdout, "@@@@ ucons_meet_formula_apply_univ: segm[%zu]=%zu(%zu)\n",
               i, look->segments[i], g->n[i]);
      fflush (stdout);
#endif
    }
  pattern_t *a = NULL;
  HASH_FIND (hh, r->udcons, look, keylen, a);
  if (!a)
    {
      // pattern is not yet there, put (a copy of) the computed data
      checked_malloc (a, pattern_t, 1,
                      sizeof (pattern_t)+ (u_seg) * sizeof (size_t), return NULL;);
      memset (a, 0, sizeof (pattern_t)+(u_seg) * sizeof (size_t));
      a->key.type = look->type;
      for (size_t i = 0; i < u_seg; i++)
        a->key.segments[i] = look->segments[i];
      a->dcons = ap_abstract0_copy (pr->man_dcons, data);

      HASH_ADD (hh, r->udcons, key, keylen, a);
      if (add_pattern_n2p (pr, r, look) == NULL)
        {
#ifndef NDEBUG1
          fprintf (stderr, "@@@@ ucons_meet_formula_apply_univ: error add_pattern_n2p returns NULL.\n");
          fflush (stderr);
#endif
        }
    }
  else
    {
      // pattern is there, do the meet
      if (a->dcons != NULL)
        a->dcons = ap_abstract0_meet (pr->man_dcons, true, a->dcons, data);
      else
        a->dcons = ap_abstract0_copy (pr->man_dcons, data);
    }
  // free the allocated data
  ap_abstract0_free (pr->man_dcons, data);
  free (look);

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_formula_apply_univ: returns r=(");
  ucons_fdump (stdout, pr->man, r);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  return r;
}

/**
 * Meets a with the disjunct disj of formula f.
 * Works directly on a.
 */
ucons_t *
ucons_meet_formula (ucons_internal_t* pr, ucons_t* a,
                    sh_formula_t* f, size_t disj)
{

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_formula: disjunct %zu\n", disj);
  fflush (stdout);
#endif

  if (!f || !f->form[disj])
    return NULL;

  // check dimensions
  size_t datadim = f->env->intdim;
  size_t segmdim = f->form[disj]->nodes->realdim + 1;
  if (!a || datadim != a->datadim || segmdim != a->segmentdim)
    return NULL;

  ucons_t* r = a;

  // Step 1: meet with r->econs
  size_t i;
  for (i = 0; i < f->form[disj]->length_dform; i++)
    {
      // deal with the constraint f->form[disj]->dform[i]
      sh_dataform_t df = f->form[disj]->dform[i];
      // build a constraint compatible with econs
      ap_lincons0_t lcons = ucons_meet_formula_build_data (pr, datadim, segmdim,
                                                           0, NULL, &df);
#ifndef NDEBUG1
      fprintf (stdout, "@@@@ ucons_meet_formula_build_data: f[%zu].dataform[%zu]=(",
               disj, i);
      ap_lincons0_fprint (stdout, &lcons, NULL);
      fprintf (stdout, ")\n");
      fflush (stdout);
#endif
      // meet with this constraint
      r = ucons_meet_econs_apply (pr, r, &lcons);
      ap_lincons0_clear (&lcons);
    }

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_formula: after meet econs r=(");
  ucons_fdump (stdout, pr->man, r);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  // Step 2: meet with the corresponding values in r->udcons hashtable
  for (i = 0; i < f->form[disj]->length_uform; i++)
    {
      // deal with the constraint f->form[disj]->dform[i]
      sh_univform_t* du = &(f->form[disj]->uform[i]);

      // build the array of constraints from the guard and the left part
      ap_lincons0_array_t arr = ucons_meet_formula_build_univ (pr, datadim, segmdim,
                                                               du);

      // meet the array above with the constraint for guard in du
      r = ucons_meet_formula_apply_univ (pr, r, &(du->guard), &arr);

      // free the allocated memory
      ap_lincons0_array_clear (&arr);
    }

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_meet_formula: after meet ucons r=(");
  ucons_fdump (stdout, pr->man, r);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  return r;
}

/* NOT IMPLEMENTED */
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

/**
 *  Only done if real dimension (ptr vars) are used in constraints
 */
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
  ucons_t * r = ucons_alloc_internal (pr, a->datadim, a->segmentdim);

  ap_linexpr0_realloc (expr, a->datadim + 2 * a->segmentdim);

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_assign_linexpr: assign %d:=", d);
  ap_linexpr0_fprint (stdout, expr, NULL);
  fprintf (stdout, " on ucons a=(");
  ucons_fdump (stdout, pr->man, a);
  fprintf (stdout, ")\n");
#endif


  //pattern_t *s,*aux;
  ap_abstract0_t * aux_dcons;

  size_t i, u_seg, e_seg;
  unsigned keylen;


  pattern_t * s, *ra, *rt;

  for (s = a->udcons; s != NULL; s = s->hh.next)
    {
      /*
       * TODO add dimensions to r->econs to cope with s->dcons; meet between r->econs and s->dcons ;
       * */
      u_seg = pr->PI[s->key.type].u_seg;
      e_seg = pr->PI[s->key.type].e_seg;

      keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);

      HASH_FIND (hh, r->udcons, &s->key, keylen, rt);
      if (rt)
        {
          rt->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, false, s->dcons, d, expr, NULL);
        }

      else
        {
          checked_malloc (ra, pattern_t, 1, sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t), return NULL;);
          memset (ra, 0, sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t));


          ra->dcons = NULL;
          ra->key.type = s->key.type;

          for (i = 0; i < (u_seg + e_seg); i++)
            ra->key.segments[i] = s->key.segments[i];
          ra->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, false, s->dcons, d, expr, NULL);
          HASH_ADD (hh, r->udcons, key, keylen, ra);
          add_pattern_n2p (pr, r, &ra->key);
        }

    }

  r->econs = ap_abstract0_assign_linexpr (pr->man_dcons, false, a->econs, d, expr, NULL);

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_assign_linexpr: returns ");
  ucons_fdump (stdout, pr->man, r);
  fprintf (stdout, "\n");
#endif

  if (destructive)
    ucons_free_internal (pr, a);
  return r;
}

/**
 * Used for pre-image computation
 * NOT IMPLEMENTED
 */
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

/**
 *  Only done if real dimension (ptr vars) are used in constraints.
 */
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

  ucons_internal_t * pr =
          ucons_init_from_manager (man, AP_FUNID_ASSIGN_LINEXPR_ARRAY, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a;
}

/**
 * Used for pre-image computation
 * NOT IMPLEMENTED
 */
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

  ucons_internal_t * pr =
          ucons_init_from_manager (man, AP_FUNID_SUBSTITUTE_TEXPR_ARRAY, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a;
}

/* NOT IMPLEMENTED */
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

/* Used only for pre-image computation */
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

/**
 * Implements the reduced product between multisets and universals.
 */
ucons_t *
ucons_strengthen (ucons_internal_t* pr, ucons_t* a, size_t **eqdim)
{
  ucons_t *r = ucons_copy_internal (pr, a);
#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_strenghten: on r =(");
  ucons_fdump (stdout, pr->man, r);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  close_eq (eqdim, a->segmentdim);

#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_strengthen: eqdim closed=(");
  for (size_t ll = 0; ll < a->segmentdim; ll++)
    for (size_t pp = 0; pp < a->segmentdim; pp++)
      if ((pp != ll) && eqdim[pp][ll] == 1)
        fprintf (stdout, " ms(%zu) <= ms(%zu)  ", pp, ll);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  //	while(!end){
  //		end = true;
  size_t i, j;
  pattern_key_t * look = NULL;
  // look for a pattern forall y
  checked_malloc (look, pattern_key_t, sizeof (pattern_key_t) + 1 * sizeof (size_t),
                  1, return NULL;);
  memset (look, 0, sizeof (pattern_key_t) + 1 * sizeof (size_t));
  look->type = 0;
  unsigned keylen = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  pattern_key_t * look2 = NULL;
  checked_malloc (look2, pattern_key_t, sizeof (pattern_key_t) + 1 * sizeof (size_t),
                  1, return NULL;);
  memset (look2, 0, sizeof (pattern_key_t) + 1 * sizeof (size_t));
  look2->type = 3;
  unsigned keylen2 = 1 * sizeof (size_t) + sizeof (pattern_key_t);

  for (i = 1; i < r->segmentdim; i++)
    for (j = 1; j < r->segmentdim; j++)
      //if ((i!=j) && ((eqdim[i][j]==2) || (eqdim[i][j]==1)) )
      if ((i != j) && ((eqdim[i][j] == 1)))
        {

#ifndef NDEBUG1
          fprintf (stdout, "@@@@ ucons_strenghten: node i=%zu with node j=%zu\n", i, j);
          fflush (stdout);
#endif

          pattern_t * ni_dcons = NULL;
          pattern_t * ni2_dcons = NULL;
          pattern_t * nj_dcons = NULL;
          ap_abstract0_t* eaux = NULL;

          look->segments[0] = i;
          HASH_FIND (hh, r->udcons, look, keylen, ni_dcons);
          look->segments[0] = j;
          HASH_FIND (hh, r->udcons, look, keylen, nj_dcons);

          if (nj_dcons && (nj_dcons->dcons != NULL) &&
              !ap_abstract0_is_bottom (pr->man_dcons, nj_dcons->dcons))
            {

              //end = false;

              ap_abstract0_t* aux = NULL;
              ap_abstract0_t* eaux = NULL;

              look2->segments[0] = i;
              HASH_FIND (hh, r->udcons, look2, keylen2, ni2_dcons);

              if (eqdim[i][j] == 1)
                {
                  /* strenghten hd(ni)\cup tl(ni)\subseteq hd(nj)\cup tl(nj) w.r.t. the patten forall y */
                  ap_dimchange_t dimadd;
                  ap_dimchange_init (&dimadd, 2, 0);
                  dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
                  dimadd.dim[0] = r->datadim + 2 * r->segmentdim;
                  dimadd.dim[1] = r->datadim + 2 * r->segmentdim;

                  aux = ap_abstract0_add_dimensions (pr->man_dcons, false,
                                                     r->econs, &dimadd, false);

                  ap_dimchange_clear (&dimadd);

                  ap_lincons0_array_t arr = ap_lincons0_array_make (1);
                  // d(y) - d(nj) ==0
                  arr.p[0].constyp = AP_CONS_EQ;
                  arr.p[0].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                         r->datadim +
                                                         2 * r->segmentdim + 2);
                  arr.p[0].scalar = NULL;
                  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                    r->datadim + j, 1);
                  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                    r->datadim +
                                                    2 * r->segmentdim + 1, -1);
                  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

                  //								// ? de verif ca trebuie y == 0
                  //								arr.p[1].constyp = AP_CONS_EQ;
                  //								arr.p[1].linexpr0 =
                  //										ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * r->segmentdim + 2);
                  //								arr.p[1].scalar = NULL;
                  //								ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                  //										r->datadim + 2* r->segmentdim, 1);
                  //								ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);


                  aux = ap_abstract0_meet_lincons_array (pr->man_dcons, true,
                                                         aux, &arr);

                  ap_lincons0_array_clear (&arr);

                  ap_abstract0_t* aux2 = ap_abstract0_copy (pr->man_dcons,
                                                            nj_dcons->dcons);

                  //								ap_dim_t * fdim = (ap_dim_t *) malloc ((r->segmentdim +1) * sizeof (ap_dim_t));
                  //								size_t kk;
                  //								for(kk=0; kk<r->segmentdim ;kk++)
                  //									fdim[kk] = r->datadim + r->segmentdim + kk;
                  //								fdim[kk] = r->datadim + r->segmentdim + kk;

                  ap_dim_t * fdim = (ap_dim_t *) malloc (1 * sizeof (ap_dim_t));
                  fdim[0] = r->datadim + 2 * r->segmentdim;
                  /* forget the length constraints */
                  aux2 = ap_abstract0_forget_array (pr->man_dcons, true, aux2,
                                                    fdim, 1, false);
                  free (fdim);
                  fdim = NULL;

                  aux = ap_abstract0_join (pr->man_dcons, true, aux, aux2);

                  arr = ap_lincons0_array_make (2);
                  //  y - 1 >= 0
                  arr.p[0].constyp = AP_CONS_SUPEQ;
                  arr.p[0].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                         a->datadim +
                                                         2 * r->segmentdim + 2);
                  arr.p[0].scalar = NULL;
                  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                    r->datadim + 2 * r->segmentdim,
                                                    1);
                  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);
                  //  l[ni] - y - 1 >= 0
                  arr.p[1].constyp = AP_CONS_SUPEQ;
                  arr.p[1].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                         a->datadim +
                                                         2 * r->segmentdim + 2);
                  arr.p[1].scalar = NULL;
                  ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                                    r->datadim + r->segmentdim + i,
                                                    1);
                  ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                                    r->datadim + 2 * r->segmentdim,
                                                    -1);
                  ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, -1);

                  aux = ap_abstract0_meet_lincons_array (pr->man_dcons, true,
                                                         aux, &arr);

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
              arr.p[0].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                     a->datadim +
                                                     2 * r->segmentdim + 2);
              arr.p[0].scalar = NULL;
              ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                r->datadim + i,
                                                1);
              ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                r->datadim + 2 * r->segmentdim + 1,
                                                -1);
              ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

              eaux = ap_abstract0_meet_lincons_array (pr->man_dcons, false,
                                                      nj_dcons->dcons, &arr);

              ap_lincons0_array_clear (&arr);

              ap_linexpr0_t *expr_y = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                         r->datadim +
                                                         2 * r->segmentdim + 2);
              ap_linexpr0_set_cst_scalar_int (expr_y, 0);

              ap_linexpr0_set_coeff_scalar_int (expr_y, r->datadim + i, 1);
              ap_dim_t dy = r->datadim + 2 * r->segmentdim + 1;

              eaux = ap_abstract0_substitute_linexpr (pr->man_dcons, false,
                                                      nj_dcons->dcons, dy,
                                                      expr_y, NULL);
              ap_linexpr0_free (expr_y);

              ap_dimchange_t dimrm;
              dimrm.realdim = 0;
              dimrm.intdim = 2;
              dimrm.dim = (ap_dim_t *) malloc (dimrm.intdim * sizeof (ap_dim_t));
              dimrm.dim[0] = r->datadim + 2 * r->segmentdim;
              dimrm.dim[1] = r->datadim + 2 * r->segmentdim + 1;


              eaux = ap_abstract0_remove_dimensions (pr->man_dcons, true,
                                                     eaux, &dimrm);
              free (dimrm.dim);

              //						ap_dim_t * fdim = (ap_dim_t *) malloc (r->segmentdim * sizeof (ap_dim_t));
              //						for(size_t i=0; i<r->segmentdim ;i++)
              //							fdim[i] = r->datadim + r->segmentdim + i;
              //						/* forget the length constraints */
              //						eaux = ap_abstract0_forget_array(pr->man_dcons,true,eaux,
              //								fdim,r->segmentdim,true);
              //						free(fdim);
              //						fdim = NULL;


              if (eqdim[i][j] == 1)
                {
                  /* hd(ni)\subseteq tl(nj) \CUP HD(nj)  */
                  ap_abstract0_t *eaux2 = NULL;

                  ap_lincons0_array_t arr = ap_lincons0_array_make (1);
                  // d(nj) - d(ni) ==0
                  arr.p[0].constyp = AP_CONS_EQ;
                  arr.p[0].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                         a->datadim +
                                                         2 * r->segmentdim);
                  arr.p[0].scalar = NULL;
                  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                    r->datadim + i, 1);
                  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                    r->datadim + j, -1);
                  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);
                  eaux2 = ap_abstract0_meet_lincons_array (pr->man_dcons, false,
                                                           r->econs, &arr);
                  ap_lincons0_array_clear (&arr);



                  ap_linexpr0_t *expr_y = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                             r->datadim +
                                                             2 * r->segmentdim);
                  ap_linexpr0_set_cst_scalar_int (expr_y, 0);

                  ap_linexpr0_set_coeff_scalar_int (expr_y, r->datadim + i, 1);
                  ap_dim_t dnj = r->datadim + j;
                  eaux2 = ap_abstract0_substitute_linexpr (pr->man_dcons, false,
                                                           r->econs,
                                                           dnj, expr_y, NULL);
                  ap_linexpr0_free (expr_y);
#ifndef NDEBUG1
                  fprintf (stdout, "@@@@ ucons_strenghten: problem after a->econs\n");
                  fflush (stdout);
#endif

                  eaux = ap_abstract0_join (pr->man_dcons, true, eaux, eaux2);
#ifndef NDEBUG1
                  fprintf (stdout, "@@@@ ucons_strenghten: problem after join\n");
                  fflush (stdout);
#endif
                  //insure equal lengths

                  //ap_lincons0_array_t
                  arr = ap_lincons0_array_make (1);
                  // l(nj) - l(ni) >=0
                  arr.p[0].constyp = AP_CONS_SUPEQ;
                  arr.p[0].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                         a->datadim +
                                                         2 * r->segmentdim);
                  arr.p[0].scalar = NULL;
                  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                    r->datadim + r->segmentdim + j,
                                                    1);
                  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                    r->datadim + r->segmentdim + i,
                                                    -1);
                  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

                  eaux = ap_abstract0_meet_lincons_array (pr->man_dcons, true,
                                                          eaux, &arr);

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
              fprintf (stdout, "@@@@ ucons_strenghten: before meet a->econs, r=(");
              ucons_fprint_econs (stdout, pr->man, r, NULL);
              fprintf (stdout, ")\n");
              fflush (stdout);
#endif
              if (eaux)
                r->econs = ap_abstract0_meet (pr->man_dcons, true, eaux, r->econs);
#ifndef NDEBUG1
              fprintf (stdout, "@@@@ ucons_strenghten: after meet a->econs r=(");
              ucons_fprint_econs (stdout, pr->man, r, NULL);
              fprintf (stdout, ")\n");
              fflush (stdout);
#endif

              /* join the 2 universals */
              ap_abstract0_t * auxe1 = ap_abstract0_copy (pr->man_dcons, r->econs);
              ap_dimchange_t dimadd1;

              ap_dimchange_init (&dimadd1, 2, 0);
              dimadd1.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
              dimadd1.dim[0] = r->datadim + 2 * r->segmentdim + 0;
              dimadd1.dim[1] = r->datadim + 2 * r->segmentdim + 0;

              auxe1 = ap_abstract0_add_dimensions (pr->man_dcons, true, auxe1,
                                                   &dimadd1, false);
              ap_dimchange_clear (&dimadd1);

              ap_abstract0_t * aux_y = ap_abstract0_copy (pr->man_dcons, aux);

              if (ni_dcons)
                {
                  //ap_abstract0_t * aux2 = ap_abstract0_copy(pr->man_dcons,aux);
                  ni_dcons->dcons = ap_abstract0_meet (pr->man_dcons, true,
                                                       ni_dcons->dcons, aux);
                  ni_dcons->dcons = ap_abstract0_meet (pr->man_dcons, true,
                                                       ni_dcons->dcons, auxe1);
                }
              else
                {
                  //end = false;
                  checked_malloc (ni_dcons, pattern_t, 1,
                                  sizeof (pattern_t) + 1 * sizeof (size_t),
                                  return NULL;);
                  memset (ni_dcons, 0, sizeof (pattern_t) + 1 * sizeof (size_t));
                  ni_dcons->key.type = look->type;
                  look->segments[0] = i;

                  ni_dcons->key.segments[0] = look->segments[0];
                  ni_dcons->dcons = ap_abstract0_copy (pr->man_dcons, aux);
                  ni_dcons->dcons = ap_abstract0_meet (pr->man_dcons, true,
                                                       ni_dcons->dcons, auxe1);

                  ap_abstract0_free (pr->man_dcons, aux);
                  HASH_ADD (hh, r->udcons, key, keylen, ni_dcons);
                  r = add_pattern_n2p (pr, r, look);
                }

              /* computation for y1<=y2 constraint */
              ap_abstract0_t *auxe2 = ap_abstract0_copy (pr->man_dcons, r->econs);

              ap_dimchange_t dimadd2;
              ap_dimchange_init (&dimadd2, 4, 0);
              dimadd2.dim = (ap_dim_t *) malloc (4 * sizeof (ap_dim_t));
              dimadd2.dim[0] = r->datadim + 2 * r->segmentdim + 0;
              dimadd2.dim[1] = r->datadim + 2 * r->segmentdim + 0;
              dimadd2.dim[2] = r->datadim + 2 * r->segmentdim + 0;
              dimadd2.dim[3] = r->datadim + 2 * r->segmentdim + 0;

              auxe2 = ap_abstract0_add_dimensions (pr->man_dcons, true, auxe2,
                                                   &dimadd2, false);
              ap_dimchange_clear (&dimadd2);


              //add dims  to aux before intersection
              //add y1 to renforce y2
              ap_dimchange_t dimadd;

              ap_dimchange_init (&dimadd, 2, 0);
              dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
              dimadd.dim[0] = r->datadim + 2 * r->segmentdim + 1;
              dimadd.dim[1] = r->datadim + 2 * r->segmentdim + 2;
              ap_abstract0_t * aux3 = ap_abstract0_add_dimensions (pr->man_dcons, false,
                                                                   aux_y, &dimadd, false);
              ap_dimchange_clear (&dimadd);

              //add y2 to renforce y1
              ap_dimchange_init (&dimadd, 2, 0);
              dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
              dimadd.dim[0] = r->datadim + 2 * r->segmentdim + 0;
              dimadd.dim[1] = r->datadim + 2 * r->segmentdim + 1;
              aux_y = ap_abstract0_add_dimensions (pr->man_dcons, true,
                                                   aux_y, &dimadd, false);
              ap_dimchange_clear (&dimadd);

              if (ni2_dcons)
                {

#ifndef NDEBUG1
                  fprintf (stdout, "@@@@ ucons_strenghten: before strengthen P12: \n");
                  fprintf (stdout, " strengthen: \n");
                  //						look->segments[0] = j;
                  //						ucons_fprint_dcons(stdout, pr->man, r, NULL);
                  fprintf (stdout, "\n");
                  fflush (stdout);
#endif
                  ni2_dcons->dcons = ap_abstract0_meet (pr->man_dcons, true,
                                                        ni2_dcons->dcons, aux3);

                  ni2_dcons->dcons = ap_abstract0_meet (pr->man_dcons, true,
                                                        ni2_dcons->dcons, aux_y);

                  ni2_dcons->dcons = ap_abstract0_meet (pr->man_dcons, true,
                                                        ni2_dcons->dcons, auxe2);

#ifndef NDEBUG1
                  fprintf (stdout, "@@@@ ucons_strengthen: results (");
                  ap_abstract0_fprint (stdout, pr->man_dcons, ni2_dcons->dcons, NULL);
                  fprintf (stdout, ")\n");
                  fflush (stdout);
#endif

                }
              else
                {
                  //end = false;
                  checked_malloc (ni2_dcons, pattern_t, 1,
                                  sizeof (pattern_t) + 1 * sizeof (size_t),
                                  return NULL;);
                  memset (ni2_dcons, 0, sizeof (pattern_t) + 1 * sizeof (size_t));
                  ni2_dcons->key.type = look2->type;
                  look2->segments[0] = i;

                  ni2_dcons->key.segments[0] = look2->segments[0];

                  aux_y = ap_abstract0_meet (pr->man_dcons, true, aux3, aux_y);

                  ap_lincons0_array_t arr2 = ap_lincons0_array_make (1);
                  // y2 - y1 >=0
                  arr2.p[0].constyp = AP_CONS_SUPEQ;
                  arr2.p[0].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                          a->datadim +
                                                          2 * r->segmentdim + 4);
                  arr2.p[0].scalar = NULL;
                  ap_linexpr0_set_coeff_scalar_int (arr2.p[0].linexpr0,
                                                    r->datadim + 2 * r->segmentdim,
                                                    -1);
                  ap_linexpr0_set_coeff_scalar_int (arr2.p[0].linexpr0,
                                                    r->datadim + 2 * r->segmentdim + 1,
                                                    1);
                  ap_linexpr0_set_cst_scalar_int (arr2.p[0].linexpr0, 0);

                  aux_y = ap_abstract0_meet_lincons_array (pr->man_dcons, true,
                                                           aux_y, &arr2);

                  ap_lincons0_array_clear (&arr2);


                  ni2_dcons->dcons = ap_abstract0_copy (pr->man_dcons, aux_y);
                  ni2_dcons->dcons = ap_abstract0_meet (pr->man_dcons, true,
                                                        ni2_dcons->dcons, auxe2);

                  HASH_ADD (hh, r->udcons, key, keylen, ni2_dcons);
                  r = add_pattern_n2p (pr, r, look2);
                  ap_abstract0_free (pr->man_dcons, aux_y);
                }

#ifndef NDEBUG1
              fprintf (stdout, "@@@@ ucons_strenghten: after meet a->udcons, r=(");
              ucons_fdump (stdout, pr->man, r);
              fprintf (stdout, ")\n");
              fflush (stdout);
#endif

            }
          else
            {
              //insure equal lengths

              ap_lincons0_array_t arr = ap_lincons0_array_make (1);
              // l(nj) - l(ni) >=0
              arr.p[0].constyp = AP_CONS_SUPEQ;
              arr.p[0].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                     a->datadim + 2 * r->segmentdim);
              arr.p[0].scalar = NULL;
              ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                r->datadim + r->segmentdim + j,
                                                1);
              ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                r->datadim + r->segmentdim + i,
                                                -1);
              ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

              r->econs = ap_abstract0_meet_lincons_array (pr->man_dcons, true,
                                                          r->econs, &arr);

              ap_lincons0_array_clear (&arr);

              /*  ms(ni) inclus ms(nj)
               * if (len(ni) == 1 and len(nj) == 1
               * then hd(ni) == hd(nj)
               * */
              if (test_singleton (pr->man_dcons, r->econs, r->datadim,
                                  r->segmentdim, i)
                  && test_singleton (pr->man_dcons, r->econs, r->datadim,
                                     r->segmentdim, j))
                {

#ifndef NDEBUG1
                  fprintf (stdout, "@@@@ ucons_strenghten: singletons i=%zu and j=%zu \n", i, j);
                  fflush (stdout);
#endif
                  arr = ap_lincons0_array_make (1);
                  // d(nj) - d(ni) ==0
                  arr.p[0].constyp = AP_CONS_EQ;
                  arr.p[0].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                         a->datadim +
                                                         2 * r->segmentdim);
                  arr.p[0].scalar = NULL;
                  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                    r->datadim + j, 1);
                  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                                    r->datadim + i, -1);
                  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

                  r->econs = ap_abstract0_meet_lincons_array (pr->man_dcons, true,
                                                              r->econs, &arr);

                  ap_lincons0_array_clear (&arr);
                }


            }
          /// TODO MODIFSSSS
          //	eqdim[i][j]=0;
        }
  if (look) free (look);
  look = NULL;

  if (look2) free (look2);
  look2 = NULL;
  //}
#ifndef NDEBUG1
  fprintf (stdout, "@@@@ ucons_strenghten: returns r=(\n");
  ucons_fdump (stdout, pr->man, r);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif
  return r;
}

void
close_eq (size_t ** tdim, size_t size)
{
  size_t i, j, k;
  size_t n = size;
  bool b = false;
  while (!b)
    {
      b = true;
      for (i = 1; i < n; i++)
        for (j = 1; j < n; j++)
          {
            if (tdim[i][j] == 1)
              {
                for (k = 1; k < n; k++)
                  {
                    if (tdim[j][k] == 1)
                      {
                        if (tdim[i][k] == 0)
                          {
                            tdim[i][k] = 1;
                            b = false;
                          }
                        if (tdim[k][i] == 0)
                          {

                            tdim[k][i] = 1;
                            b = false;
                          }
                      }
                  }
              }
          }
    }
}

/**
 * Computes the abstract value from a 
 * which has saturated universal constraints for nodes in
 * nn (of size size).
 */
ucons_t *
ucons_saturation_all (ucons_internal_t* pr, ucons_t * a);
ucons_t *
ucons_saturation_one (ucons_internal_t* pr, ucons_t * a, size_t n);

ucons_t *
ucons_saturation (ucons_internal_t* pr, ucons_t * a, size_t* nn, size_t size)
{
  ucons_t* res = NULL;
  switch (size)
    {
    case 0:
      res = ucons_saturation_all (pr, a);
      break;
    case 1:
      res = ucons_saturation_one (pr, a, nn[0]);
      break;
    default:
#ifndef NDEBUG
      fprintf (stdout, "Ucons_saturation: case (2) not yet implemented.\n");
      fflush (stdout);
#endif
      res = ucons_copy_internal (pr, a);
      break;
    }
  return res;

}

/**
 * Saturates witn uniformely the universal constraints.
 * Works on a copy of a.
 */
ucons_t *
ucons_saturation_all (ucons_internal_t* pr, ucons_t * a)
{

  pattern_t *s = NULL;
  size_t u_seg, e_seg, nr_y;
  unsigned keylen;
  pattern_t *rt = NULL;
  ucons_t* r = ucons_copy_internal (pr, a);

#ifndef NDEBUG
  fprintf (stdout, "Ucons_saturation: case all, input a=(");
  ucons_fdump (stdout, pr->man, a);
  fprintf (stdout, ")\n");
#endif

  for (s = r->udcons; s != NULL; s = s->hh.next)
    {
      u_seg = pr->PI[s->key.type].u_seg;
      e_seg = pr->PI[s->key.type].e_seg;
      nr_y = pr->PI[s->key.type].nr_y;

      // add usual constraints y1>=1 ...
      ap_lincons0_array_t arr = ap_lincons0_array_make (nr_y);
      for (size_t i = 0; i < nr_y; i++)
        {

#ifndef NDEBUG
          fprintf (stdout, "Ucons_saturation: case all do y%zu>=1\n", i);
#endif

          // yi - 1 >= 0
          arr.p[i].constyp = AP_CONS_SUPEQ;
          arr.p[i].linexpr0 = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                 a->datadim +
                                                 2 * a->segmentdim + 2 * nr_y);
          arr.p[i].scalar = NULL;
          ap_linexpr0_set_coeff_scalar_int (arr.p[i].linexpr0,
                                            a->datadim +
                                            2 * a->segmentdim + i, 1);
          ap_linexpr0_set_cst_scalar_int (arr.p[i].linexpr0, -1);

        }
      s->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                                  true, s->dcons, &arr);
      ap_lincons0_array_clear (&arr);

      // pattern Gle2
      // - add constraint y1 <= y2
      if (pr->PI[s->key.type].kind == pattern_1_2)
        {
#ifndef NDEBUG
          fprintf (stdout, "Ucons_saturation: case all, guard Gle2, do y2>=y1\n");
#endif
          // y2 >= y1
          arr = ap_lincons0_array_make (1);
          arr.p[0].constyp = AP_CONS_SUPEQ;
          arr.p[0].linexpr0 =
                  ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim +
                                     2 * a->segmentdim + 2 * nr_y);
          arr.p[0].scalar = NULL;
          ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                            a->datadim + 2 * a->segmentdim,
                                            -1);
          ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                            a->datadim + 2 * a->segmentdim + 1,
                                            1);
          ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

          s->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                                      true, s->dcons, &arr);
          ap_lincons0_array_clear (&arr);
        }

      // pattern Gsucc2
      // - add constraint y1+1=y2
      // - add constraint from Gle2 if not empty
      if (pr->PI[s->key.type].kind == pattern_succ_1_2)
        {
          // y2=y1+1 
#ifndef NDEBUG
          fprintf (stdout, "Ucons_saturation: case all, guard Gsucc2, do y2=y1+1\n");
#endif
          arr = ap_lincons0_array_make (1);
          arr.p[0].constyp = AP_CONS_SUPEQ;
          arr.p[0].linexpr0 =
                  ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim +
                                     2 * a->segmentdim + 2 * nr_y);
          arr.p[0].scalar = NULL;
          ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                            a->datadim + 2 * a->segmentdim,
                                            -1);
          ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                            a->datadim + 2 * a->segmentdim + 1,
                                            1);
          ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);

          s->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                                      true, s->dcons, &arr);
          ap_lincons0_array_clear (&arr);

          // meet with U from \forall y1,y2 y1 <= y2 = U 
          keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);
          pattern_key_t *look = NULL;
          checked_malloc (look, pattern_key_t,
                          sizeof (pattern_key_t) + u_seg * sizeof (size_t), 1,
                          return NULL;);
          memset (look, 0, sizeof (pattern_key_t) + u_seg * sizeof (size_t));
          look->type = pattern_1_2; // PI[0] \forall y1,y2 y1 <= y2
          look->segments[0] = s->key.segments[0];
          HASH_FIND (hh, r->udcons, look, keylen, rt);

          if (rt != NULL && rt->dcons != NULL)
            {
              s->dcons = ap_abstract0_meet (pr->man_dcons, true, s->dcons, rt->dcons);
#ifndef NDEBUG
              fprintf (stdout, "Ucons_saturation: case all, Gsucc2, do meet Gle2\n");
#endif
            }
        }

      // patterns Gle2, Gsucc2
      // - add constraints in Gall
      if (pr->PI[s->key.type].kind == pattern_1_2 ||
          pr->PI[s->key.type].kind == pattern_succ_1_2)
        {
          // meet with \forall y 
          keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);
          pattern_key_t *look = NULL;
          checked_malloc (look, pattern_key_t,
                          sizeof (pattern_key_t) + u_seg * sizeof (size_t), 1,
                          return NULL;);
          memset (look, 0, sizeof (pattern_key_t) + u_seg * sizeof (size_t));
          look->type = pattern_1; // PI[0] \forall y
          look->segments[0] = s->key.segments[0];
          HASH_FIND (hh, r->udcons, look, keylen, rt);

          if (rt != NULL && rt->dcons != NULL)
            {
              ap_abstract0_t * aux = ap_abstract0_copy (pr->man_dcons, rt->dcons);

              ap_dimchange_t dimadd;
              ap_dimchange_init (&dimadd, 2, 0);
              dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
              dimadd.dim[0] = a->datadim + 2 * a->segmentdim;
              dimadd.dim[1] = a->datadim + 2 * a->segmentdim + 1;
              aux = ap_abstract0_add_dimensions (pr->man_dcons, true, aux,
                                                 &dimadd, false);
              ap_dimchange_clear (&dimadd);

              s->dcons = ap_abstract0_meet (pr->man_dcons, true, s->dcons, aux);

              aux = ap_abstract0_copy (pr->man_dcons, rt->dcons);

              ap_dimchange_init (&dimadd, 2, 0);
              dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
              dimadd.dim[0] = a->datadim + 2 * a->segmentdim + 1;
              dimadd.dim[1] = a->datadim + 2 * a->segmentdim + 2;
              aux = ap_abstract0_add_dimensions (pr->man_dcons, true, aux,
                                                 &dimadd, false);
              ap_dimchange_clear (&dimadd);

              s->dcons = ap_abstract0_meet (pr->man_dcons, true, s->dcons, aux);
#ifndef NDEBUG
              fprintf (stdout, "Ucons_saturation: case all, guards Gle2 and Gsucc2, meet Gall\n");
#endif
            }
        }
      // TODO: add the transfer relations from 
      // - Gle2, Gsucc2 to Gall, Gfst or Glst
      // - close by equality, if needed
    }

#ifndef NDEBUG
  fprintf (stdout, "Ucons_saturation: case all, result r=(");
  ucons_fdump (stdout, pr->man, r);
  fprintf (stdout, ")\n");
#endif

  return r;
}

ucons_t *
ucons_saturation_one (ucons_internal_t* pr, ucons_t * a, size_t n)
{
  ucons_t* res = ucons_copy_internal (pr, a);
  size_t nb_anon = pr->max_anon * pr->segm_anon;

  // Step 1: unfold until max_anon anonymous are generated, i.e., max_anon+1 times
  // array of dimensions added by split, stored to fold them again
  ap_dim_t* tdim = (ap_dim_t*) malloc (sizeof (ap_dim_t) * (nb_anon + 1));
  tdim[0] = n;

  for (size_t i = 1; i <= nb_anon; i++)
    {
      // after split, the node added has dimension res->segmentdim
      tdim[i] = res->segmentdim;
      res = ucons_split (pr, true, res, tdim[i - 1]);
    }
  // Split last anonymous node 
  res = ucons_split (pr, true, res, tdim[nb_anon]);

  // Step 2: fold all except last node
  res = ucons_fold (pr->man, true, res, tdim, nb_anon + 1);
  // remove dimensions folded
  ap_dimchange_t dimchange;
  dimchange.intdim = 0;
  dimchange.realdim = nb_anon;
  dimchange.dim = (ap_dim_t*) malloc (sizeof (ap_dim_t) * nb_anon);
  memcpy (dimchange.dim, &(tdim[1]), nb_anon * sizeof (ap_dim_t));
  res = ucons_remove_dimensions (pr->man, true, res, &dimchange);

  // Step 3: unfold again last node until max anon
  tdim[0] = n;
  tdim[1] = res->segmentdim - 1; // TODO: check that it is the last node generated above
  for (size_t i = 2; i <= nb_anon; i++)
    {
      // after split, the node added has dimension res->segmentdim
      tdim[i] = res->segmentdim;
      res = ucons_split (pr, true, res, tdim[i - 1]);
    }
  // set to 1 the last anonymous node 
  res = ucons_singleton (pr, true, res, tdim[nb_anon]);

  // Step 4: fold all
  res = ucons_fold (pr->man, true, res, tdim, nb_anon + 1);
  dimchange.dim = (ap_dim_t*) malloc (sizeof (ap_dim_t) * nb_anon);
  memcpy (dimchange.dim, &(tdim[1]), nb_anon * sizeof (ap_dim_t));
  res = ucons_remove_dimensions (pr->man, true, res, &dimchange);

  return res;
}
