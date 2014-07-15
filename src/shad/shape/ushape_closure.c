/**************************************************************************/
/*                                                                        */
/*  CELIA Tools / Shape Abstract Domain                                   */
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


#include "ushape.h"
#include "ushape_internal.h"


/* All closures are in-place. */


/* ============================================================ */
/* Full Closure */
/* ============================================================ */


/* NOT IMPLEMENTED */
bool
ushape_close (ushape_t * g, size_t dim)
{
  if ((g != g) || (dim != dim))
    return true;                /* to remove warning on unused parameter */
  return false;
}


/* ============================================================ */
/* Incremental Closure */
/* ============================================================ */

/* NOT IMPLEMENTED */
bool
ushape_close_incremental (ushape_t * m, size_t dim, size_t v)
{
  if ((m != m) || (dim != dim) || (v != v))
    return true;                /* to remove warning on unused parameter */
  return false;
}


/* ============================================================ */
/* Sanity Check */
/* ============================================================ */

/* NOT IMPLEMENTED */
bool
ushape_check_closed (ushape_t * m, size_t dim)
{
  if ((m != m) || (dim != dim))
    return true;                /* to remove warning on unused parameter */
  return false;
}


/* ============================================================ */
/* Topological closure operation */
/* ============================================================ */

/**
 * @brief Returns a copy of the input value depending on @p destructive.
 * 
 * The semantics of the APRON library is not implemented.
 * The copy of the value eliminates the anonymous nodes.
 */
ushape_t *
ushape_closure (ap_manager_t * man, bool destructive, ushape_t * a)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_CLOSURE, 0);
  if (destructive)
    return a;
  else
    return ushape_copy_internal (pr, a);
}
