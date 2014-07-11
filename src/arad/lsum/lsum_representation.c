/**************************************************************************/
/*                                                                        */
/*  CELIA Tools / LSUM Abstract Domain                                    */
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
#include "lsum.h"
#include "lsum_internal.h"
#include "sh_macros.h"
#include "apron2shape.h"
#include "ap_abstract0.h"
#include "ap_generic.h"
#include "pk.h"
#include "ap_ppl.h"



/* ============================================================ */
/* Internal Management */
/* ============================================================ */

/* generic allocation routine */
inline lsum_t *
lsum_alloc_internal (lsum_internal_t * pr, size_t intdim, size_t realdim)
{
  lsum_t *r;
  if (pr != pr)
    return NULL;                /* to remove warning on unused parameter */
  checked_malloc (r, lsum_t, 1, sizeof (lsum_t), return NULL;
    );
  r->datadim = intdim;
  r->segmdim = realdim;
  r->dcons = NULL;
  return r;
}

/* returns a top lsum */
inline lsum_t *
lsum_alloc_top (lsum_internal_t * pr, size_t intdim, size_t realdim)
{
  lsum_t *r = lsum_alloc_internal (pr, intdim, realdim);
  r->dcons = ap_abstract0_top (pr->man_dcons, intdim + 3 * r->segmdim, 0);
  return r;
}

inline void
lsum_free_internal (lsum_internal_t * pr, lsum_t * a)
{
  if (a)
    {
      if (a->dcons)
        {
          ap_abstract0_free (pr->man_dcons, a->dcons);
          a->dcons = NULL;
        }
      free (a);
    }
  return;
}

inline lsum_t *
lsum_copy_internal (lsum_internal_t * pr, lsum_t * a)
{
  arg_assert (a, return NULL;
    );
  lsum_t *r = lsum_alloc_internal (pr, a->datadim, a->segmdim);
  r->dcons = ap_abstract0_copy (pr->man_dcons, a->dcons);
  return r;
}


/* ============================================================ */
/* Memory */

/* ============================================================ */

lsum_t *
lsum_copy (ap_manager_t * man, lsum_t * a)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_COPY, 0);
  return (a) ? lsum_copy_internal (pr, a) : NULL;
}

void
lsum_free (ap_manager_t * man, lsum_t * a)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_FREE, 0);
  lsum_free_internal (pr, a);
}

size_t
lsum_size (ap_manager_t * man, lsum_t * a)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_ASIZE, 0);
  if (!a)
    return 0;
  size_t s = ap_abstract0_size (pr->man_dcons, a->dcons);
  return s;
}


/* ============================================================ */
/* Control of internal representation */
/* ============================================================ */

/* Return the set of minimal elements */
void
lsum_minimize (ap_manager_t * man, lsum_t * a)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_MINIMIZE, 0);
  ap_abstract0_minimize (pr->man_dcons, a->dcons);
}

/* Return the set of canonical elements */
void
lsum_canonicalize (ap_manager_t * man, lsum_t * a)
{
  lsum_internal_t *pr =
    lsum_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);
  ap_abstract0_canonicalize (pr->man_dcons, a->dcons);
}

/* NOT IMPLEMENTED */
int
lsum_hash (ap_manager_t * man, lsum_t * a)
{
  if (a != a)
    return 1;                   /* to remove warning on unused parameter */
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_HASH, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return 0;
}

/* NOT IMPLEMENTED */
void
lsum_approximate (ap_manager_t * man, lsum_t * a, int algorithm)
{
  if ((a != a) || (algorithm != algorithm))
    return;                     /* to remove warning on unused parameter */

  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_APPROXIMATE, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
}

/* TODO: used? */
bool
lsum_is_minimal (ap_manager_t * man, lsum_t * a)
{
  if (a != a)
    return false;               /* to remove warning on unused parameter */
  lsum_internal_t *pr =
    lsum_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);
  return true;
}

/* TODO: used? */
bool
lsum_is_canonical (ap_manager_t * man, lsum_t * a)
{
  if (a != a)
    return false;               /* to remove warning on unused parameter */
  lsum_internal_t *pr =
    lsum_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);
  return true;
}

/* ============================================================ */
/* Basic Constructors */

/* ============================================================ */

lsum_t *
lsum_bottom (ap_manager_t * man, size_t intdim, size_t realdim)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_BOTTOM, 0);
  lsum_t *r = lsum_alloc_internal (pr, intdim, realdim);
  // all constraints are NULL, i.e., bottom
  return r;
}

lsum_t *
lsum_top (ap_manager_t * man, size_t intdim, size_t realdim)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_TOP, 0);
  lsum_t *r = lsum_alloc_top (pr, intdim, realdim);
  return r;
}

/* put constraints on data variables */
lsum_t *
lsum_of_box (ap_manager_t * man, size_t intdim, size_t realdim,
             ap_interval_t ** t)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_OF_BOX, 0);
  lsum_t *r = lsum_alloc_top (pr, intdim, realdim);
  ap_abstract0_free (pr->man_dcons, r->dcons);
  r->dcons = ap_abstract0_of_box (pr->man_dcons, intdim, 0, t); /* TODO */
  return r;
}

/* NOT IMPLEMENTED */
lsum_t *
lsum_of_generator_array (ap_manager_t * man, size_t intdim, size_t realdim,
                         ap_generator0_array_t * ar)
{
  if ((intdim != intdim) || (realdim != realdim) || (ar != ar))
    return NULL;                /* to remove warning on unused parameter */
  lsum_internal_t *pr =
    lsum_init_from_manager (man, AP_FUNID_ADD_RAY_ARRAY, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return NULL;
}


/* ============================================================ */
/* Accessors */

/* ============================================================ */

ap_dimension_t
lsum_dimension (ap_manager_t * man, lsum_t * a)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_DIMENSION, 0);
  ap_dimension_t r;
  r.intdim = 0;
  r.realdim = 0;
  if (a)
    {
      r.intdim = a->datadim;
      r.realdim = a->segmdim;
    }
  return r;
}

/* ============================================================ */
/* Managers */

/* ============================================================ */

void
lsum_internal_free (lsum_internal_t * pr)
{
  // ap_pkgrid_manager_free(pr->man_dcons);
  pr->man_dcons = NULL;
  free (pr);
}

ap_manager_t *
lsum_manager_alloc (void)
{
  ap_manager_t *man;
  lsum_internal_t *pr;
  size_t i;

  pr = (lsum_internal_t *) malloc (sizeof (lsum_internal_t));
  assert (pr);

  // pr->man_dcons = pk_manager_alloc(false);
  pr->man_dcons = ap_ppl_poly_manager_alloc (false);
  man = ap_manager_alloc ("lsum", "0.1 with (dcons=pk)", pr,
                          (void (*)(void *)) lsum_internal_free);

  pr->man = man;

  man->funptr[AP_FUNID_COPY] = &lsum_copy;
  man->funptr[AP_FUNID_FREE] = &lsum_free;
  man->funptr[AP_FUNID_ASIZE] = &lsum_size;
  man->funptr[AP_FUNID_MINIMIZE] = &lsum_minimize;
  man->funptr[AP_FUNID_CANONICALIZE] = &lsum_canonicalize;
  man->funptr[AP_FUNID_HASH] = &lsum_hash;
  man->funptr[AP_FUNID_APPROXIMATE] = &lsum_approximate;
  man->funptr[AP_FUNID_FPRINT] = &lsum_fprint;
  man->funptr[AP_FUNID_FPRINTDIFF] = &lsum_fprintdiff;
  man->funptr[AP_FUNID_FDUMP] = &lsum_fdump;
  man->funptr[AP_FUNID_SERIALIZE_RAW] = &lsum_serialize_raw;
  man->funptr[AP_FUNID_DESERIALIZE_RAW] = &lsum_deserialize_raw;
  man->funptr[AP_FUNID_BOTTOM] = &lsum_bottom;
  man->funptr[AP_FUNID_TOP] = &lsum_top;
  man->funptr[AP_FUNID_OF_BOX] = &lsum_of_box;
  man->funptr[AP_FUNID_DIMENSION] = &lsum_dimension;
  man->funptr[AP_FUNID_IS_BOTTOM] = &lsum_is_bottom;
  man->funptr[AP_FUNID_IS_TOP] = &lsum_is_top;
  man->funptr[AP_FUNID_IS_LEQ] = &lsum_is_leq;
  man->funptr[AP_FUNID_IS_EQ] = &lsum_is_eq;
  man->funptr[AP_FUNID_IS_DIMENSION_UNCONSTRAINED] =
    &lsum_is_dimension_unconstrained;
  man->funptr[AP_FUNID_SAT_INTERVAL] = &lsum_sat_interval;
  man->funptr[AP_FUNID_SAT_LINCONS] = &lsum_sat_lincons;
  man->funptr[AP_FUNID_SAT_TCONS] = &lsum_sat_tcons;
  man->funptr[AP_FUNID_BOUND_DIMENSION] = &lsum_bound_dimension;
  man->funptr[AP_FUNID_BOUND_LINEXPR] = &lsum_bound_linexpr;
  man->funptr[AP_FUNID_BOUND_TEXPR] = &lsum_bound_texpr;
  man->funptr[AP_FUNID_TO_BOX] = &lsum_to_box;
  man->funptr[AP_FUNID_TO_LINCONS_ARRAY] = &lsum_to_lincons_array;
  man->funptr[AP_FUNID_TO_TCONS_ARRAY] = &lsum_to_tcons_array;
  man->funptr[AP_FUNID_TO_GENERATOR_ARRAY] = &lsum_to_generator_array;
  man->funptr[AP_FUNID_MEET] = &lsum_meet;
  man->funptr[AP_FUNID_MEET_ARRAY] = &lsum_meet_array;
  man->funptr[AP_FUNID_MEET_LINCONS_ARRAY] = &lsum_meet_lincons_array;
  man->funptr[AP_FUNID_MEET_TCONS_ARRAY] = &lsum_meet_tcons_array;
  man->funptr[AP_FUNID_JOIN] = &lsum_join;
  man->funptr[AP_FUNID_JOIN_ARRAY] = &lsum_join_array;
  man->funptr[AP_FUNID_ADD_RAY_ARRAY] = &lsum_add_ray_array;
  man->funptr[AP_FUNID_ASSIGN_LINEXPR_ARRAY] = &lsum_assign_linexpr_array;
  man->funptr[AP_FUNID_SUBSTITUTE_LINEXPR_ARRAY] =
    &lsum_substitute_linexpr_array;
  man->funptr[AP_FUNID_ASSIGN_TEXPR_ARRAY] = &lsum_assign_texpr_array;
  man->funptr[AP_FUNID_SUBSTITUTE_TEXPR_ARRAY] = &lsum_substitute_texpr_array;
  man->funptr[AP_FUNID_ADD_DIMENSIONS] = &lsum_add_dimensions;
  man->funptr[AP_FUNID_REMOVE_DIMENSIONS] = &lsum_remove_dimensions;
  man->funptr[AP_FUNID_PERMUTE_DIMENSIONS] = &lsum_permute_dimensions;
  man->funptr[AP_FUNID_FORGET_ARRAY] = &lsum_forget_array;
  man->funptr[AP_FUNID_EXPAND] = &lsum_expand;
  man->funptr[AP_FUNID_FOLD] = &lsum_fold;
  man->funptr[AP_FUNID_WIDENING] = &lsum_widening;
  man->funptr[AP_FUNID_CLOSURE] = &lsum_closure;

  for (i = 0; i < AP_EXC_SIZE; i++)
    ap_manager_set_abort_if_exception (man, i, false);

  return man;
}

lsum_t *
lsum_of_abstract0 (ap_abstract0_t * a)
{
  return (lsum_t *) a->value;
}

ap_abstract0_t *
abstract0_of_lsum (ap_manager_t * man, lsum_t * a)
{
  ap_abstract0_t *r = malloc (sizeof (ap_abstract0_t));
  assert (r);
  r->value = a;
  r->man = ap_manager_copy (man);
  return r;
}
