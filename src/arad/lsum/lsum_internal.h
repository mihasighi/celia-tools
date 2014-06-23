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


#ifndef __LSUM_INTERNAL_H_
#define __LSUM_INTERNAL_H_

#include "shape_options.h"
#include "lsum_fun.h"
#include "shad.h"

/* *INDENT-OFF* */
#ifdef __cplusplus
extern "C"
{
#endif
  /* *INDENT-ON* */

  /* ********************************************************************** */
  /* Manager */
  /* ********************************************************************** */

  /* ============================================================ */
  /* Internal Representation */
  /* ============================================================ */

  /* manager-local data specific to universal constraints */
  struct _lsum_internal_t
  {
    /* current function */
    ap_funid_t funid;

    /* local parameters for current function */
    ap_funopt_t *funopt;

    /* manager of the length domain */
    ap_manager_t *man_lcons;

    /* manager of the data domain */
    ap_manager_t *man_dcons;

    /* count errors */
    int error_;

    /* back-pointer */
    ap_manager_t *man;
  };


  /* ============================================================ */
  /* Basic management. */
  /* ============================================================ */

  /* options for setting basic abstract domains */
  /* For the moment, the abstract domains used are:
   * - if relational then lcons: null, dcons: polyhedra
   * - if non relational then lcons: octagons, dcons: polyhedra
   */

  /* called by each function to setup and get manager-local data */
  static inline lsum_internal_t *
  lsum_init_from_manager (ap_manager_t * man, ap_funid_t id, size_t size)
  {
    lsum_internal_t *pr = (lsum_internal_t *) man->internal;
    pr->funid = id;
    pr->funopt = man->option.funopt + id;
    man->result.flag_exact = man->result.flag_best = true;

    return pr;
  }


  /* ********************************************************************** */
  /* List sum constraints */
  /* ********************************************************************** */

  /* ============================================================ */
  /* Internal Representation */
  /* ============================================================ */

  /** Listsum constraints maintain relations between the sums
   * of elements of list segments. For this, it also represents constraints
   * on data variables, length variables, and data of first nodes in
   * segments.
   * So we have a product domain on lcons /\ dcons
   */

  struct _lsum_t
  {
    size_t datadim; /* number of scalar data variables (length and data) */
    size_t segmdim; /* number of segments represented */
    void *dcons; /* constraint on data of dimension
				 * datadim + 3 segmdim */
  };


  /* ============================================================ */
  /* Internal Management */
  /* ============================================================ */

  lsum_t *lsum_alloc_internal (lsum_internal_t * pr, size_t intdim,
                               size_t realdim);
  lsum_t *lsum_alloc_top (lsum_internal_t * pr, size_t intdim, size_t realdim);
  void lsum_free_internal (lsum_internal_t * pr, lsum_t * a);
  lsum_t *lsum_copy_internal (lsum_internal_t * pr, lsum_t * a);


  /* ============================================================ */
  /* Extensions */
  /* ============================================================ */

  lsum_t* lsum_meet_formula (lsum_internal_t* pr, lsum_t* a,
                             sh_formula_t* f, size_t disj);
  /* Meets with a constraint from an SL3 disjunct. */

  /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */




#endif /* __LSUM_INTERNAL_H_ */
