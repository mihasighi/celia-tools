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


#include "lsum.h"
#include "lsum_internal.h"


/* All closures are in-place. */


/* ============================================================ */
/* Full Closure */
/* ============================================================ */


/* NOT IMPLEMENTED */
bool
lsum_close (lsum_t * a, size_t dim)
{
  if ((a != a) || (dim != dim))
    return true;                /* to remove warning on unused parameter */

  return false;
}


/* ============================================================ */
/* Incremental Closure */
/* ============================================================ */

/* NOT IMPLEMENTED */
bool
lsum_close_incremental (lsum_t * a, size_t dim, size_t v)
{
  if ((a != a) || (v != v) || (dim != dim))
    return true;                /* to remove warning on unused parameter */

  return false;
}


/* ============================================================ */
/* Sanity Check */
/* ============================================================ */

/* NOT IMPLEMENTED */
bool
lsum_check_closed (lsum_t * a, size_t dim)
{
  if ((a != a) || (dim != dim))
    return true;                /* to remove warning on unused parameter */

  return false;
}


/* ============================================================ */
/* Topological closure operation */
/* ============================================================ */

lsum_t *
lsum_closure (ap_manager_t * man, bool destructive, lsum_t * a)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_CLOSURE, 0);
  if (!a)
    return NULL;
  lsum_t *r = lsum_alloc_internal (pr, a->datadim, a->segmdim);
  r->dcons = ap_abstract0_closure (pr->man_dcons, false, a->dcons);
  if (destructive)
    lsum_free_internal (pr, a);
  return r;
}
