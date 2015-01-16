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
 * sh_test.c: test of shad interface
 */

#include <assert.h>
#include "sh_manager.h"
#include "sh_utils.h"

/* ====================================================================== */
/* I. Types, Field, and predicates */
/* ====================================================================== */


/* ====================================================================== */
/* II. Variables, Frames, stack */
/* ====================================================================== */


/* ********************************************************************** */
/* III. Test functions */
/* ********************************************************************** */

int
main(void) 
{
  sh_manager_t* sman = sh_manager_alloc("shad", "0.1", NULL, NULL);

  fprintf (stdout, "Step 0: Initialized manager:");
  sh_manager_fdump(stdout, sman);

  sh_file_register (sman, "list.c");
  sh_type_register_ls (sman, "list.c", 3, "list", "next", "data");

  fprintf (stdout, "Step 1: Manager with types:");
  sh_manager_fdump(stdout, sman);
}
