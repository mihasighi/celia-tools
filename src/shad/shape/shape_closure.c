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


#include "shape.h"
#include "shape_internal.h"


/* All closures are in-place. */


/* ============================================================ */
/* Full Closure */
/* ============================================================ */


/* NOT IMPLEMENTED */
bool
shape_close (shape_internal_t * pr, shape_t * a)
{
  if ((pr != pr) || (a != a))
    return false;               /* to remove warning on unused parameter */
  return false;
}


/* ============================================================ */
/* Incremental Closure */
/* ============================================================ */

/* NOT IMPLEMENTED */
bool
shape_close_incremental (shape_t * m, size_t dim, size_t v)
{
  if ((m != m) || (dim != dim) || (v != v))
    return false;               /* to remove warning on unused parameter */
  return false;
}


/* ============================================================ */
/* Sanity Check */
/* ============================================================ */

/* NOT IMPLEMENTED */
bool
shape_check_closed (shape_t * m, size_t dim)
{
  if ((m != m) || (dim != dim))
    return false;               /* to remove warning on unused parameter */
  return false;
}


/* ============================================================ */
/* Topological closure operation */
/* ============================================================ */

/**
 * @brief Should eliminate anonymous nodes 
 * NOT IMPLEMENTED
 */
shape_t *
shape_closure (ap_manager_t * man, bool destructive, shape_t * a)
{
  if (destructive != destructive)
    return NULL;

  shape_internal_t *pr = shape_init_from_manager (man, AP_FUNID_CLOSURE, 0);
  if (destructive)
    return a;
  else
    return shape_copy_internal (pr, a);
}
