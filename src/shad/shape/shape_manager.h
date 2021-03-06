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


#ifndef SHAPE_MANAGER_H_
#define SHAPE_MANAGER_H_

#include "hgraph_fun.h"
#include "shape_fun.h"
#include "ap_pcons0.h"
#include "ap_passign0.h"
#include "ap_manager.h"

/* *INDENT-OFF* */
#ifdef __cplusplus
extern          "C"
{
#endif
  /* *INDENT-ON* */

  /* ********************************************************************** */
  /* Manager */
  /* ********************************************************************** */

  /* ============================================================ */
  /* Internal Representation */
  /* ============================================================ */

  /* manager-local data specific to shapes */
struct _shape_internal_t
{

  /* current function */
  ap_funid_t funid;

  /* local parameters for current function */
  ap_funopt_t *funopt;

  /* manager of the segment domains */
  size_t size_scons;
  ap_manager_t **man_scons;     /* array of scons_size managers */
  size_t man_mset;              /* position of the mset manager/constraint */
  size_t man_ucons;             /* position of the ucons manager/constraint */

  /* max number for anonymous nodes for closure */
  size_t max_anon;
  size_t segm_anon;

  /* count errors and files */
  int error_;
  long int filenum;

  /* default dimensions */
  size_t intdim;
  size_t realdim;

  /* approximate meet for hgraphs */
  int meet_algo;

  /* TODO: other data */

  /* back-pointer */
  ap_manager_t *man;
};

typedef struct _shape_internal_t shape_internal_t;
  /* Abstract data type of library-specific manager options. */

  /* ============================================================ */
  /* Basic management. */
  /* ============================================================ */

  /* called by each function to setup and get manager-local data */
static inline shape_internal_t *
shape_init_from_manager (ap_manager_t * man, ap_funid_t id, size_t size)
{
  shape_internal_t *pr = (shape_internal_t *) man->internal;
  pr->funid = id;
  pr->funopt = man->option.funopt + id;
  man->result.flag_exact = man->result.flag_best = true;

  /* TODO: set other internal data from manager */

  return pr;
}

  /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */


#endif /* SHAPE_MANAGER_H_ */
