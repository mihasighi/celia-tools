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


#include "ushape.h"
#include "ushape_internal.h"


/* All closures are in-place. */


/* ============================================================ */
/* Full Closure */
/* ============================================================ */


/*
 * TODO: priority 1 What it is here closure? - close the hgraph and the lcons
 */

bool
ushape_close (ushape_t * g, size_t dim)
{
  return false;
}


/* ============================================================ */
/* Incremental Closure */
/* ============================================================ */

/*
 * TODO: priority 1
 *
 * ?? time. ?? space.
 */

bool
ushape_close_incremental (ushape_t * m, size_t dim, size_t v)
{
  return false;
}



/* ============================================================ */
/* Sanity Check */
/* ============================================================ */

/* TODO: priority 1 */
/* Check if closed field is not NULL if yes test if h field is closed */
bool
ushape_check_closed (ushape_t * m, size_t dim)
{
  return false;
}


/* ============================================================ */
/* Topological closure operation */
/* ============================================================ */

/* TODO: priority 3 */
/* Eliminate anonimous nodes */
ushape_t *
ushape_closure (ap_manager_t * man, bool destructive, ushape_t * a)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_CLOSURE, 0);
  if (destructive)
    return a;
  else
    return ushape_copy_internal (pr, a);
}
