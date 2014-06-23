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


#include "apron2shape.h"
#include "ushape.h"
#include "ushape_internal.h"
#include "shape.h"
#include "shape_internal.h"
#include "shape_macros.h"
#include "ap_generic.h"

/* ============================================================ */
/* Meet and Join */

/* ============================================================ */

shape_t *
shape_meet (ap_manager_t * man, bool destructive, shape_t * a1, shape_t * a2)
{
  shape_internal_t *pr = shape_init_from_manager (man, AP_FUNID_MEET, 0);
  size_t i, j, size;
  ushape_array_t arr;
  shape_t *r;
  /* simple cases */
  if (!a1 || !a2)
    return NULL;
  if (shape_is_bottom (man, a1))
    {
      if (destructive)
        {
          // destructive is only for a1, not for a2
          return a1;
        }
      else
        return shape_bottom (man, a2->intdim, a2->realdim);
    }
  if (shape_is_bottom (man, a2))
    {
      if (destructive)
        {
          shape_free_internal (pr, a1);
          a1 = shape_bottom (man, a2->intdim, a2->realdim);
        }
      else
        return shape_bottom (man, a2->intdim, a2->realdim);
    }
#ifndef NDEBUG1
  fprintf (stdout, "\n****shape_meet: on ");
  shape_fdump (stdout, man, a1);
  fprintf (stdout, " and ");
  shape_fdump (stdout, man, a2);
  fprintf (stdout, ")");
  fflush (stdout);
#endif
  /* general case: do meet between each pair of ushapes */
  arr.p = NULL;
  arr.size = 0;
  ushape_array_init (pr, &arr, a1->m.size);
  size = 0;
  for (i = 0; i < a1->m.size; i++)
    if (a1->m.p[i])
      {
        for (j = 0; j < a2->m.size; j++)
          if (a2->m.p[j])
            {
              ushape_t *ur = ushape_meet (man, false, a1->m.p[i], a2->m.p[j]);
              if (ur && !ushape_is_bottom (man, ur))
                size += ushape_array_add (pr, true, &arr, size, false, false, ur); /* do not copy and
											 * destroy */
              else if (!ur)
                ushape_free (man, ur);
            }
      }
  if (size == 0)
    r = shape_bottom (man, a1->intdim, a1->realdim);
  else
    {
      ushape_array_t *rr;
      r = shape_alloc_internal (pr, size);
      rr = ushape_array_copy (pr, &arr, size);
      r->m = *rr;
      rr->p = NULL;
      free (rr);
      r->set = true;
      r->closed = a1->closed && a2->closed;
      r->msize = size;
      r->intdim = a1->intdim;
      r->realdim = a1->realdim;
    }
#ifndef NDEBUG1
  fprintf (stdout, "\n****shape_meet returns: ");
  shape_fdump (stdout, man, r);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  if (destructive)
    shape_free (man, a1);
  // meet done, change the flag meet_algo
  shape_approximate (man, r, 0);
  return r;
}

shape_t *
shape_join (ap_manager_t * man, bool destructive, shape_t * a1, shape_t * a2)
{
  shape_internal_t *pr = shape_init_from_manager (man, AP_FUNID_JOIN, 0);
  shape_t *r;
  size_t i, j, size;
  bool *isin1, *isin2;
  /* simple cases */
  if (!a1 && !a2)
    return NULL;
  if (shape_is_bottom (man, a1))
    {
      if (destructive)
        {
          shape_free_internal (pr, a1);
          return a2;
        }
      else
        return shape_copy (man, a2);
    }
  if (shape_is_bottom (man, a2))
    {
      if (destructive)
        {
          shape_free_internal (pr, a2);
          return a1;
        }
      else
        return shape_copy (man, a1);
    }

#ifndef NDEBUG1
  fprintf (stdout, "\n****shape_join: on ");
  shape_fdump (stdout, man, a1);
  fprintf (stdout, " and ");
  shape_fdump (stdout, man, a2);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  /* general case, do the union */
  r = shape_alloc_internal (pr, a1->msize + a2->msize);
  isin1 = (bool *) malloc (a1->msize * sizeof (bool));
  memset (isin1, 0, a1->msize * sizeof (bool));
  isin2 = (bool *) malloc (a2->msize * sizeof (bool));
  memset (isin2, 0, a2->msize * sizeof (bool));
  size = 0;
  for (j = 0; j < a2->msize; j++)
    {
      for (i = 0; i < a1->msize; i++)
        if (!isin1[i])
          {
            ushape_t *ur = ushape_join (man, false, a1->m.p[i], a2->m.p[j]);
            if (ur && !ushape_is_bottom (man, ur))
              {
                isin1[i] = true;
                isin2[j] = true;
                size += ushape_array_add (pr, true, &r->m, size, false, false, ur); /* do not copy, nor
											 * destroy */
              }
          }
      if (!isin2[j])
        // nothing added for a2->m.p[i]
        size += ushape_array_add (pr, true, &r->m, size, true, false, a2->m.p[j]); /* copy but not destroy */
    }
  for (i = 0; i < a1->msize; i++)
    if (!isin1[i])
      size +=
            ushape_array_add (pr, true, &r->m, size, true, false, a1->m.p[i]);
  r->set = true;
  r->closed = false;
  r->msize = size;
  r->intdim = a1->intdim;
  r->realdim = a1->realdim;
  if (destructive)
    {
      shape_free_internal (pr, a1);
      shape_free_internal (pr, a2);
    }
#ifndef NDEBUG1
  fprintf (stdout, "\n****shape_join returns: ");
  shape_fdump (stdout, man, r);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  return r;
}

shape_t *
shape_meet_array (ap_manager_t * man, shape_t ** tab, size_t size)
{
  shape_internal_t *pr =
          shape_init_from_manager (man, AP_FUNID_MEET_ARRAY, 0);
  arg_assert (size > 0, return NULL;
              );
  shape_t *r = shape_copy_internal (pr, tab[0]);
  size_t i;
  for (i = 1; i < size && !shape_is_bottom (man, r); i++)
    {
      shape_t *rr = shape_meet (man, false, tab[i], r);
      shape_free (man, r);
      r = rr;
    }
  return r;
}

shape_t *
shape_join_array (ap_manager_t * man, shape_t ** tab, size_t size)
{
  shape_internal_t *pr =
          shape_init_from_manager (man, AP_FUNID_JOIN_ARRAY, 0);
  arg_assert (size > 0, return NULL;
              );
  shape_t *r = shape_copy_internal (pr, tab[0]);
  size_t i;
  for (i = 1; i < size && !shape_is_top (man, r); i++)
    {
      shape_t *rr = shape_join (man, false, tab[i], r);
      shape_free (man, r);
      r = rr;
    }
  return r;
}

/* ============================================================ */
/* Meet constraints and Join generators */
/* ============================================================ */

/*
 * Compute the effect of constraints in array (of size elements) on the shape
 * a. Better performances are obtained if array is sorted, @see
 * <code>shape_pcons_array_of_lincons_array</code>.
 * Modified to read SL3 constraints.
 */
bool
is_sl3_pred (pcons0_array_t* array)
{
  return array && (array->size == 1) &&
          (array->p[0] != NULL) &&
          (array->p[0]->type == SL3_CONS);
}

shape_t *
shape_meet_pcons_array (shape_internal_t * pr, bool destructive,
                        shape_t * a, pcons0_array_t * array)
{
  shape_t *rs = NULL;
  if (is_sl3_pred (array))
    {
#ifndef NDEBUG1
      fprintf (stdout, "\n****shape_meet_sl3: with constraint ");
      shape_pcons_array_fdump (stdout, array);
      fprintf (stdout, "\n");
      fflush (stdout);
#endif
     // Step 1: read the SL3 spec, transform it into a shape_t,
      // then push it with its number
      ap_coeff_t* c = ap_linexpr0_cstref (array->p[0]->info.data.cons.linexpr0);
      double spec_no;
      ap_double_set_scalar (&spec_no, c->val.scalar, GMP_RNDN);
      bool kind = (array->p[0]->info.data.cons.constyp == AP_CONS_EQ) ? true : false;
      shape_t* spec_f = shape_of_spec (pr, (int) spec_no, true);

      // Step 2a: if a is TOP (pre-condition) then 
      // return spec_f
      if (shape_is_top (pr->man, a))
        return spec_f;
      // Step 2b: otherwise (inv check or post-condition) then
      // (TODO: not the correct semantics)
      // if test inclusion of a in spec_f then
      //   rs=spec_f contraint par a->h
      // else bottom
      // if (kind==false) 
      //   inverser le resultat (si rs=bottom take a else return bottom)
      shape_approximate (pr->man, a, -2);
      rs = shape_meet (pr->man, false, a, spec_f);
      if (!kind)
        {
#ifndef NDEBUG1
      fprintf (stdout, "\n****shape_meet_sl3: with negative constraint\n");
      fflush (stdout);
#endif
          if (shape_is_bottom (pr->man, rs))
            rs = shape_copy (pr->man, a);
          else
            {
              shape_free (pr->man, rs);
              rs = NULL;
            }
        }
    }
  else
    {
      ushape_array_t *r = NULL;
      size_t i;
      // for TOP, one value exists in the array
      for (i = 0; i < a->msize; i++)
        {
          ushape_array_t *rr =
                  ushape_meet_pcons_array (pr, false, a->m.p[i], array);
          if (rr)
            {
              r = ushape_array_add_array (pr, true, r, rr);
              ushape_array_init (pr, rr, rr->size);
              free (rr);
            }
        }
      if (!r)
        return NULL;
      /* compute msize */
      for (i = 0; i < r->size; i++)
        if (!r->p[i])
          break;
      checked_malloc (rs, shape_t, sizeof (shape_t), 1, return NULL;);
      rs->m = *r;
      r->p = NULL;
      free (r);
      rs->msize = i;
      rs->set = rs->closed = false;
      rs->intdim = a->intdim;
      rs->realdim = a->realdim;
    }
  if (destructive)
    shape_free_internal (pr, a);

  return rs;
}

shape_t *
shape_meet_lincons_array (ap_manager_t * man,
                          bool destructive, shape_t * a,
                          ap_lincons0_array_t * array)
{
  shape_internal_t *pr =
          shape_init_from_manager (man, AP_FUNID_MEET_LINCONS_ARRAY, 0);
  if (shape_is_bottom (man, a) || !array || array->size == 0)
    /* nothing to do */
    return (destructive) ? a : shape_bottom (man, pr->intdim, pr->realdim);
  else
    {
      shape_t *b;
      pcons0_array_t *arr;
      shape_t *r;
      if (!destructive)
        b = shape_copy_internal (pr, a);
      else
        b = a;
      /* compute in arr the constraints sorted */
      arr =
              shape_pcons_array_of_lincons_array (pr, array, a->intdim, a->realdim);
#ifndef NDEBUG1
      fprintf (stdout, "\n****shape_meet_lincons_array: with constraint ");
      ap_lincons0_array_fprint (stdout, array, NULL);
      fprintf (stdout, "(i.e, ");
      shape_pcons_array_fdump (stdout, arr);
      fprintf (stdout, ")\n\t on  ");
      shape_fdump (stdout, man, b);
      fprintf (stdout, "\n");
      fflush (stdout);
#endif
      /* go */
      r = shape_meet_pcons_array (pr, false, b, arr);
      shape_free_internal (pr, b);
      if (!r) r = shape_bottom (man, pr->intdim, pr->realdim);
#ifndef NDEBUG1
      fprintf (stdout, "\n****shape_meet_lincons_array returns: ");
      shape_fdump (stdout, man, r);
      fprintf (stdout, "\n");
      fflush (stdout);
#endif
      return r;
    }
}

shape_t *
shape_meet_tcons_array (ap_manager_t * man,
                        bool destructive, shape_t * a,
                        ap_tcons0_array_t * array)
{

  if (shape_is_bottom (man, a) || !array || array->size == 0)
    /* nothing to do */
    return (destructive || !a) ? a : shape_copy (man, a);

  shape_t *b;
  pcons0_array_t *arr;
  shape_t *r;
  shape_internal_t *pr =
          shape_init_from_manager (man, AP_FUNID_MEET_TCONS_ARRAY, 0);
  if (!destructive)
    b = shape_copy_internal (pr, a);
  else
    b = a;
  /* compute in arr the constraints sorted */
  arr = shape_pcons_array_of_tcons_array (pr, array, a->intdim, a->realdim);
#ifndef NDEBUG1
  fprintf (stdout, "\n****shape_meet_tcons_array: with constraint ");
  //ap_tcons0_array_fprint (stdout, array,NULL);
  //fprintf (stdout, "\t\t (i.e, ");
  shape_pcons_array_fdump (stdout, arr);
  //fprintf (stdout, ")\n");
  fprintf (stdout, "\non ");
  shape_fdump (stdout, man, b);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  /* go */
  r = shape_meet_pcons_array (pr, false, b, arr);
  shape_free_internal (pr, b);
#ifndef NDEBUG1
  fprintf (stdout, "\n****shape_meet_tcons_array returns: ");
  shape_fdump (stdout, man, r);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  return r;
}

/* Abstract a conjunction of constraints. Based on meet, like in ap_abstract0. */
shape_t*
shape_of_lincons_array (ap_manager_t* man,
                        size_t intdim, size_t realdim,
                        ap_lincons0_array_t* array)
{
  shape_t* res = shape_top (man, intdim, realdim);
  res = shape_meet_lincons_array (man, true, res, array);
  shape_canonicalize (man, res);
  return res;
}

/* Abstract a conjunction of constraints. Based on meet, like in ap_abstract0. */
shape_t*
shape_of_tcons_array (ap_manager_t* man,
                      size_t intdim, size_t realdim,
                      ap_tcons0_array_t* array)
{
  shape_t* res = shape_top (man, intdim, realdim);
  res = shape_meet_tcons_array (man, true, res, array);
  shape_canonicalize (man, res);
  return res;
}

/* NOT IMPLEMENTED */
shape_t *
shape_add_ray_array (ap_manager_t * man,
                     bool destructive, shape_t * a,
                     ap_generator0_array_t * array)
{
  shape_internal_t *pr =
          shape_init_from_manager (man, AP_FUNID_ADD_RAY_ARRAY, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a;
}

/* ============================================================ */
/* Interface with sl3 */
/* ============================================================ */
/*
typedef struct shape_spec_t
{
  int id; // unique id given from code, use to sort the list
  shape_t* shape_pos;
  shape_t* shape_neg;
  struct shape_spec_t* next;
 *
} shape_spec_t;
 */
// stores specifications used during the analysis
shape_spec_t* shape_specs = NULL;

/** Search spec in shape_specs and returns 
 *  - the cell with the same spec if it exists
 *  - new cell inserted in shape_specs to be sorted, otherwise
 */
shape_spec_t*
shape_spec_get_spec (int spec)
{
  shape_spec_t* s = shape_specs;
  while ((s != NULL) && (s->id < spec))
    if (s->next != shape_specs) s = s->next;
    else s = NULL;
  if (s == NULL)
    { // insert before shape_specs, and take care if shape_specs is also NULL!
      s = (shape_spec_t*) malloc (sizeof (shape_spec_t));
      s->id = spec;
      s->shape_pos = NULL;
      s->shape_neg = NULL;
      s->next = (shape_specs != NULL) ? shape_specs : s;
      s->prev = (shape_specs != NULL) ? shape_specs->prev : s;
      if (shape_specs != NULL) shape_specs->prev = s;
      else shape_specs = s;
      if (s->prev->next != s) s->prev->next = s;
      return s;
    }
  else if (s->id == spec) return s;
  else // s !=NULL && s->id > spec
    {
      // insert before s != NULL
      shape_spec_t* ns = (shape_spec_t*) malloc (sizeof (shape_spec_t));
      ns->id = spec;
      ns->shape_pos = NULL;
      ns->shape_neg = NULL;
      ns->next = s;
      ns->prev = s->prev;
      s->prev = ns;
      if (ns->prev->next != ns) ns->prev->next = ns;
      return ns;
    }
}

/**
 *  Build shapes from sl3 spec, positive and negative versions.
 */
shape_t *
shape_of_spec (shape_internal_t* pr, int spec, bool version)
{
  // search spec in the list of specs,
  // returns its crt position or a new cell to be filled
  shape_spec_t* s = shape_spec_get_spec (spec);
  // assert (s != NULL)

  // Complex case: spec not in the list
  if (s->shape_pos == NULL && s->id == spec)
    {
      // init variable manager from apron
      ap_var_operations_t* ap_var_operations_old = ap_var_operations;
      ap_var_operations = &ap_var_operations_default;
      // - parse spec
      char* filename = (char*) malloc (sizeof (char) *(14 + 10)); // TODO: limit in file_no!
      snprintf (filename, 24, "pan/spec_%d.smt", spec);
      sh_fscan (filename);
      free (filename);
      // - translate SL3 to shapes
      s->id = spec;
      sh_crt = sh_pos;
      s->shape_pos = shape_of_formula (pr->man, sh_pos);
      sh_crt = sh_neg;
      s->shape_neg = shape_of_formula (pr->man, sh_neg);
      sh_crt = NULL;
      // restore manager for vars
      ap_var_operations = ap_var_operations_old;
      ap_var_operations_old = NULL;
    }
  return (version == true) ? s->shape_pos : s->shape_neg;
}

/* Transforming an SL3 formula into a shape value.
 */
shape_t*
shape_of_formula (ap_manager_t *man, sh_formula_t* f)
{

  shape_internal_t* pr =
          shape_init_from_manager (man, AP_FUNID_OF_BOX, 0);
  shape_t* rs;
  ushape_array_t *r = NULL;
  size_t rsize = 0;
  if (!f) return NULL;
  sh_crt = f;
  r = ushape_of_formula (pr, f, &rsize);
  if (!r) return NULL;
  checked_malloc (rs, shape_t, sizeof (shape_t), 1, return NULL;);
  rs->m = *r;
  rs->msize = rsize;
  rs->set = rs->closed = false;
  rs->intdim = rs->m.p[0]->datadim;
  rs->realdim = rs->m.p[0]->ptrdim;
  r->p = NULL; // remove link inside the result
  free (r);
  sh_crt = NULL;
  return rs;
}


/* ============================================================ */
/* Widening, Narrowing */
/* ============================================================ */

/* Windening is done between elements of shapes that have an
 * isomorphic hgraph; otherwise, these elements are joined to the shape.
 */
shape_t *
shape_widening (ap_manager_t * man, shape_t * a1, shape_t * a2)
{
  shape_internal_t *pr = shape_init_from_manager (man, AP_FUNID_WIDENING, 0);
  shape_t *r;
  size_t i, j, size;
  bool *isin1, *isin2;
  /* simple cases */
  if (!a1 || !a2 || shape_is_bottom (man, a1) || shape_is_bottom (man, a2))
    return shape_bottom (man, pr->intdim, pr->realdim);

#ifndef NDEBUG1
  fprintf (stdout, "\n****shape_widening: on ");
  shape_fdump (stdout, man, a1);
  fprintf (stdout, " and ");
  shape_fdump (stdout, man, a2);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  /* general case, apply widening on ushapes with isomorphic graphs */
  r = shape_alloc_internal (pr, a1->msize + a2->msize);
  isin1 = (bool *) malloc (a1->msize * sizeof (bool));
  memset (isin1, 0, a1->msize * sizeof (bool));
  isin2 = (bool *) malloc (a2->msize * sizeof (bool));
  memset (isin2, 0, a2->msize * sizeof (bool));
  size = 0;
  for (j = 0; j < a2->msize; j++)
    {
      for (i = 0; i < a1->msize; i++)
        if (!isin1[i])
          {
            ushape_t *ur = ushape_widening (man, a1->m.p[i], a2->m.p[j]);
            if (ur && !ushape_is_bottom (man, ur))
              {
                isin1[i] = true; /* only isomorphic graphs are widen */
                isin2[j] = true;
                size += ushape_array_add (pr, true, &r->m, size, false, false, ur); /* do not copy, nor
											 * destroy */
              }
          }
      if (!isin2[j])
        // nothing added for a2->m.p[i]
        size += ushape_array_add (pr, true, &r->m, size, true, false, a2->m.p[j]); /* copy but not destroy */
    }
  /*
     for (i = 0; i < a1->msize; i++)
     if (!isin1[i])
     size += ushape_array_add (pr, true, &r->m, size, true, false, a1->m.p[i]);
   */
  r->set = true;
  r->closed = false; /* TODO: widening conserves closed property */
  r->msize = size;
  r->intdim = a1->intdim;
  r->realdim = a1->realdim;
#ifndef NDEBUG1
  fprintf (stdout, "\n****shape_widening returns: ");
  shape_fdump (stdout, man, r);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  return r;
}

/* NOT IMPLEMENTED */
shape_t *
shape_widening_thresholds (ap_manager_t * man,
                           shape_t * a1, shape_t * a2,
                           ap_scalar_t ** array, size_t nb)
{
  shape_internal_t *pr =
          shape_init_from_manager (man, AP_FUNID_WIDENING, nb + 1);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a2;
}

/* NOT IMPLEMENTED */
shape_t *
shape_narrowing (ap_manager_t * man, shape_t * a1, shape_t * a2)
{
  shape_internal_t *pr = shape_init_from_manager (man, AP_FUNID_WIDENING, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a2;
}
