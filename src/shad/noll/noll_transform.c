/* SHAD - Library of shape abstract domains
 * Copyright (C) 2012-2013 LIAFA
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * Please read the COPYING file packaged in the distribution.
 */

/* 
 * noll_representation.c: Implementation of functions related to the
 *                        data structure representation
 */

#include <assert.h>
#include "noll.h"
#include "noll_fun.h"
#include "noll_internal.h"


/* ============================================================ */
/* I.2 Control of the internal representation */
/* ============================================================ */

void 
noll_minimize (sh_manager_t * man, noll_val_t * a)
{
  return; // TODO
}

void noll_canonicalize (sh_manager_t * man, noll_val_t * a)
{
  return; // TODO
} 

int noll_hash (sh_manager_t * man, noll_val_t * a)
{
  return 1; // TODO
}

void 
noll_approximate (sh_manager_t * man, noll_val_t * a, int algorithm)
{
  return; // TODO
}

gboolean 
noll_is_minimal (sh_manager_t * man, noll_val_t * a)
{
  return TRUE; // TODO
}

gboolean 
noll_is_canonical (sh_manager_t * man, noll_val_t * a)
{
  return TRUE; // TODO
}

