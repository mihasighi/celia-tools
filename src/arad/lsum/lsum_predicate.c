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


#include "lsum.h"
#include "lsum_internal.h"
#include "ap_generic.h"



/* ============================================================ */
/* Tests */
/* ============================================================ */

/* The bottom value is bot on one constraint. */
bool
lsum_is_bottom (ap_manager_t * man, lsum_t * a)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_IS_BOTTOM, 0);
  if (!a)
    return true;
#ifndef NDEBUG2
  fprintf(stdout, "\n????? lsum_is_bottom: with a=(\n");
  ap_abstract0_fprint(stdout,pr->man_dcons, a->dcons, NULL);
  fprintf(stdout, ")\n");
  fflush(stdout);
#endif
  bool is_bot = ap_abstract0_is_bottom (pr->man_dcons, a->dcons);
#ifndef NDEBUG2
  fprintf(stdout, "\n????? lsum_is_bottom: returns %zu\n", (is_bot) ? 1: 0);
  fflush(stdout);
#endif
  return is_bot;
}

/* The top value is top in all constraints. */
bool
lsum_is_top (ap_manager_t * man, lsum_t * a)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_IS_TOP, 0);
  if (!a)
    return false;
#ifndef NDEBUG2
  fprintf(stdout, "\n????? lsum_is_top: with a=(\n");
  ap_abstract0_fprint(stdout,pr->man_dcons, a->dcons, NULL);
  fprintf(stdout, ")\n");
  fflush(stdout);
#endif
  bool is_top = ap_abstract0_is_top (pr->man_dcons, a->dcons);
#ifndef NDEBUG2
  fprintf(stdout, "\n????? lsum_is_top: returns %zu\n", (is_top) ? 1: 0);
  fflush(stdout);
#endif
  return is_top;
}

bool
lsum_is_leq (ap_manager_t * man, lsum_t * a1, lsum_t * a2)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_IS_LEQ, 0);
  if (lsum_is_bottom (man, a1))
    return true;
  if (lsum_is_bottom (man, a2))
    return false;
  if (a1->datadim != a2->datadim || a1->segmdim != a2->segmdim)
    return false;
  bool is_leq = ap_abstract0_is_leq (pr->man_dcons, a1->dcons, a2->dcons);
  return is_leq;
}

bool
lsum_is_eq (ap_manager_t * man, lsum_t * a1, lsum_t * a2)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_IS_EQ, 0);
  if ((lsum_is_bottom (man, a1) && lsum_is_bottom (man, a2)) ||
      (lsum_is_top (man, a1) && lsum_is_top (man, a2)))
    return true;
  if (a1->datadim != a2->datadim || a1->segmdim != a2->segmdim)
    return false;
  bool is_eq = ap_abstract0_is_eq (pr->man_dcons, a1->dcons, a2->dcons);
  return is_eq;
}

/* TODO: priority 1
 * Checking linear constraints between
 * program variables (pointer or scalar) in the assert statements.
 */
bool
lsum_sat_lincons (ap_manager_t * man, lsum_t * a, ap_lincons0_t * lincons)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_SAT_LINCONS, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			      "not implemented");
  return true;
}

/* TODO: priority 1
 * Checking constraints between expressions on pointer
 * program variables in the assert statements.
 */
bool
lsum_sat_tcons (ap_manager_t * man, lsum_t * a, ap_tcons0_t * cons)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_SAT_TCONS, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			      "not implemented");
  return true;
}

/* TODO: priority 1 */
/* Interval constraints are only between non-pointer variables */
bool
lsum_sat_interval (ap_manager_t * man, lsum_t * a,
		   ap_dim_t dim, ap_interval_t * i)
{
  lsum_internal_t *pr =
    lsum_init_from_manager (man, AP_FUNID_SAT_INTERVAL, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			      "not implemented");
  return true;
}

bool
lsum_is_dimension_unconstrained (ap_manager_t * man, lsum_t * a, ap_dim_t dim)
{
  lsum_internal_t *pr =
    lsum_init_from_manager (man, AP_FUNID_IS_DIMENSION_UNCONSTRAINED, 0);
  if (!a)
    return false;
  if (lsum_is_top (man, a) || dim >= (a->datadim + a->segmdim))
    return true;
  bool is_uncons =
    ap_abstract0_is_dimension_unconstrained (pr->man_dcons, a->dcons, dim);
  return is_uncons;
}

/* ============================================================ */
/* Extraction of properties */
/* ============================================================ */

/* NOT IMPLEMENTED */
ap_interval_t *
lsum_bound_linexpr (ap_manager_t * man, lsum_t * a, ap_linexpr0_t * expr)
{
  lsum_internal_t *pr =
    lsum_init_from_manager (man, AP_FUNID_BOUND_LINEXPR, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			      "not implemented");
  return NULL;
}

/* NOT IMPLEMENTED */
ap_interval_t *
lsum_bound_texpr (ap_manager_t * man, lsum_t * a, ap_texpr0_t * expr)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_BOUND_TEXPR, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			      "not implemented");
  return NULL;
}

/* NOT IMPLEMENTED */
ap_interval_t *
lsum_bound_dimension (ap_manager_t * man, lsum_t * a, ap_dim_t dim)
{
  lsum_internal_t *pr =
    lsum_init_from_manager (man, AP_FUNID_BOUND_DIMENSION, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			      "not implemented");
  return NULL;
}

/* NOT IMPLEMENTED */
ap_lincons0_array_t
lsum_to_lincons_array (ap_manager_t * man, lsum_t * a)
{
  ap_lincons0_array_t ar;
  lsum_internal_t *pr =
    lsum_init_from_manager (man, AP_FUNID_TO_LINCONS_ARRAY, 0);
  ar = ap_lincons0_array_make (1);
  ar.p[0] = ap_lincons0_make_unsat ();
  return ar;
}

/* NOT IMPLEMENTED */
ap_tcons0_array_t
lsum_to_tcons_array (ap_manager_t * man, lsum_t * a)
{
  return ap_generic_to_tcons_array (man, a);
}

/* NOT IMPLEMENTED */
ap_interval_t **
lsum_to_box (ap_manager_t * man, lsum_t * a)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_TO_BOX, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			      "not implemented");
  return NULL;
}

/* NOT IMPLEMENTED */
ap_generator0_array_t
lsum_to_generator_array (ap_manager_t * man, lsum_t * a)
{
  lsum_internal_t *pr =
    lsum_init_from_manager (man, AP_FUNID_TO_GENERATOR_ARRAY, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			      "not implemented");
  return ap_generator0_array_make (0);
}
