/**************************************************************************/
/*                                                                        */
/*  CELIA Tools / SLL Abstract Domain                                     */
/*                                                                        */
/*  Copyright (C) 2009-2014                                               */
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


#ifndef HGRAPH_H_
#define HGRAPH_H_

/* dependencies */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "sh_manager.h"

/* *INDENT-OFF* */
#ifdef __cplusplus
extern          "C"
{
#endif
  /* *INDENT-ON* */

sh_manager_t *hgraph_manager_alloc (void);
  /* Creates a new manager for the hgraph library. 
   * @see hgraph_representation.c */

  /* ============================================================ */
  /* Supplementary functions & options, not in the APRON manager. */
  /* ============================================================ */

  /* TODO: add functions not in the APRON manager. */

  /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */

#endif /* HGRAPH_H_ */
