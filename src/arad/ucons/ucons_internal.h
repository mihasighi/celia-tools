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

/*
 * Private definitions, access to internal structures and algorithms.
 * Use with care.
 */

#ifndef __UCONS_INTERNAL_H_
#define __UCONS_INTERNAL_H_

#include "uthash.h"
#include "shape_options.h"
#include "shape_macros.h"
#include "shape_manager.h"
#include "ucons_fun.h"


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



  typedef enum
  {
    pattern_1, /* y : PI[0]*/
    pattern_2_1, /* y1 in ?  /\\ y2 in ?? /\\ y1 = y2 : PI[1]*/
    pattern_3_1, /* y1 in ?  /\\ y2 in ??  /\\ y3 in ??? : PI[2] */
    pattern_1_2, /* y1 <= y2  : PI[3] */
    pattern_succ_1_3, /* y3 = y2 + 1 = y1 + 2  : PI[4] */
    pattern_1_l1, /* y /\\ y=1 : PI[5]*/
    pattern_1_lx_1, /* y /\\ y=l[?]-1 : PI[6]*/
    pattern_1_lx, /* y\in ! /\\ y = l[?]+.. : PI[7].PI[8],PI[9]*/
    //  pattern_1_lx2,            /* y\in! /\\ y = l[?]+l[??] : PI[8]*/
    //  pattern_1_lx3,            /* y\in! /\\ y = l[?]+l[??]+l[???] : PI[9]*/
    pattern_2_1_lx, /* y1\in ?, y2\in ?? /\ y1=l[!]...+y2 : PI[10], PI[11], PI[12]*/
    //  pattern_2_1_lx2,			/* y1\in ?, y2\in ?? /\ y1=l[!]+l[!!]+y2 : PI[11]*/
    //  pattern_2_1_lx3,			/* y1\in ?, y2\in ?? /\ y1= l[!] + l[!!] + l[!!!] + y2 : PI[12]*/
    pattern_2_1_mlx, /* y1\in ?, y2\in ?? /\ y1 + l[!]...= y2 : PI[10], PI[11], PI[12]*/
    //  pattern_2_1_mlx2,			/* y1\in ?, y2\in ?? /\ y1 + l[!]+l[!!]= y2 : PI[11]*/
    //  pattern_2_1_mlx3,			/* y1\in ?, y2\in ?? /\ y1+ l[!] + l[!!] + l[!!!] = y2 : PI[12]*/
    pattern_succ_1_2 /*y1,y2 \in n y2 = y1 + 1	*/

  } pattern_kind;

  typedef struct
  {
    size_t size; //number of universal vars on the segment
    bool *reach;
  } pattern_var_t;

  typedef struct
  {
    pattern_kind kind;
    size_t e_seg; // number of existential segments
    size_t u_seg; // number of universal segments
    pattern_var_t *uvar; //sizeof u_seg
    size_t lcons; // index in LI
    size_t nr_y; //number of universals on all segments
  } pattern_info_t;

  /* manager-local data specific to universal constraints */
  struct _ucons_internal_t
  {
    /* current function */
    ap_funid_t funid;

    /* local parameters for current function */
    ap_funopt_t *funopt;

    //              /* manager of the length domain */
    //              ap_manager_t *man_lcons;

    /* manager of the data domain */
    ap_manager_t *man_dcons;

    /* count errors */
    int error_;

    /* back-pointer */
    ap_manager_t *man;

    size_t segm_anon;
    size_t max_anon;
    bool *active_patterns;
    size_t nr_active;

    pattern_info_t *PI;
    /*  list of patterns PR /\ PL */

    size_t PI_size;
    /* number of patterns in the list */


  };


  /* ============================================================ */
  /* Basic management. */
  /* ============================================================ */

  /* options for setting basic abstract domains */

#if defined (UCONS_DCONS_OCT_P11)
#define UCONS_DCONS_DOMAIN DOM_OCT_P11
#elif defined (UCONS_DCONS_OCT_P12)
#define UCONS_DCONS_DOMAIN DOM_OCT_P12
#elif defined (UCONS_DCONS_OCT_P21)
#define UCONS_DCONS_DOMAIN DOM_OCT_P21
#elif defined (UCONS_DCONS_POLY_P11)
#define UCONS_DCONS_DOMAIN DOM_POLY_P11
#elif defined (UCONS_DCONS_POLY_P12)
#define UCONS_DCONS_DOMAIN DOM_POLY_P12
#elif defined (UCONS_DCONS_POLY_P21)
#define UCONS_DCONS_DOMAIN DOM_POLY_P21
#else
#define UCONS_DCONS_DOMAIN DOM_BOX
#endif

  /* called by each function to setup and get manager-local data */
  static inline ucons_internal_t *
  ucons_init_from_manager (ap_manager_t * man, ap_funid_t id, size_t size)
  {
    ucons_internal_t *pr = (ucons_internal_t *) man->internal;
    pr->funid = id;
    pr->funopt = man->option.funopt + id;
    man->result.flag_exact = man->result.flag_best = true;

    /* TODO: set other internal data from manager */
    /* DO NOT TOUCH TO THE hgraphs FIELD! */

    return pr;
  }


  /* ********************************************************************** */
  /* Universal constraints */
  /* ********************************************************************** */

  ucons_t * build_constraint (ucons_internal_t * pr, ucons_t * r,
                              ap_linexpr0_t *lexpr, ap_scalar_t *code);
  // Build the predefined constraint given by code.
  // The main nodes have positive coefficients.

  ucons_t* build_const_1 (ucons_internal_t * pr, ucons_t * r,
                          ap_linexpr0_t *lexpr);

  ucons_t* build_const_2 (ucons_internal_t * pr, ucons_t * r,
                          ap_linexpr0_t *lexpr);

  ucons_t* build_const_3 (ucons_internal_t * pr, ucons_t * r,
                          ap_linexpr0_t *lexpr);

  ucons_t* build_const_4 (ucons_internal_t * pr, ucons_t * r,
                          ap_linexpr0_t *lexpr);

  ucons_t* build_const_5 (ucons_internal_t * pr, ucons_t * r,
                          ap_linexpr0_t *lexpr);

  ucons_t* build_const_6 (ucons_internal_t * pr, ucons_t * r,
                          ap_linexpr0_t *lexpr);


  ucons_t* build_const_8 (ucons_internal_t * pr, ucons_t * r,
                          ap_linexpr0_t *lexpr);

  ucons_t* build_const_7 (ucons_internal_t * pr, ucons_t * r,
                          ap_linexpr0_t *lexpr);


  ucons_t* build_const_10 (ucons_internal_t * pr, ucons_t * r,
                           ap_linexpr0_t *lexpr);

  ucons_t* build_const_11 (ucons_internal_t * pr, ucons_t * r,
                           ap_linexpr0_t *lexpr);

  ucons_t* build_const_20 (ucons_internal_t * pr, ucons_t * r,
                           ap_linexpr0_t *lexpr);

  ucons_t* build_const_21 (ucons_internal_t * pr, ucons_t * r,
                           ap_linexpr0_t *lexpr);

  ucons_t* build_const_22 (ucons_internal_t * pr, ucons_t * r,
                           ap_linexpr0_t *lexpr);

  ucons_t* build_const_23 (ucons_internal_t * pr, ucons_t * r,
                           ap_linexpr0_t *lexpr);

  ucons_t*
  build_const_24 (ucons_internal_t * pr, ucons_t * r,
                  ap_linexpr0_t *lexpr);

  /* ============================================================ */
  /* Internal Representation */
  /* ============================================================ */

  /** Universal constraints are parameterized by a set of patterns and
   * the domains of data and length constraints. It represents
   * properties of segments of the heap having the form:
   * lcons /\ forall pattern in segment => /\i (lconsi ==> dconsi)
   */

  /** Element (lcons ==> dcons) of an universal constraint. */

  /** Data type for patterns.
   * The current implementation fixes the form of the patterns:
   * all patterns with at most 3 universals are considered.
   * The pivot is marqued by *.
   * Shall be used as an abstract data type.
   */


  typedef struct
  {
    size_t type; // index in PI
    size_t segments[]; //its size is PI[type]->u_seg
  } pattern_key_t;

  typedef struct
  {
    size_t size; //number of universals?
    pattern_key_t **p;
  } pattern_key_set_t;

  typedef struct
  {
    UT_hash_handle hh;
    ap_abstract0_t *dcons;
    pattern_key_t key; /* pattern description (reachability + length constraints) */
  } pattern_t;

  struct _ucons_t
  {
    size_t datadim; // number of program variables of type D
    size_t segmentdim; //number of segments

    size_t size;

    //bool *active_patterns; //active_patterns[0] = T <=> the pattern PI[0] is used with the universal formula

    ap_abstract0_t *econs;

    pattern_key_set_t *n2p; // size = segmentdim

    pattern_t *udcons;

  };

  /* ********************************************************************** */
  /* I. ucons_elm_t */
  /* ********************************************************************** */

  /* see ucons_predicate.c */

  //
  //      /* ********************************************************************** */
  //      /* II. ucons_set_t */
  //      /* ********************************************************************** */
  //
  //      /* see ucons_representation.c */
  //
  //      void *ucons_set_add_elm (ucons_internal_t * pr, ucons_elm_t * u, ucons_set_t * arr);
  //
  //      /* see ucons_predicate.c */
  //      bool ucons_set_isin (ucons_internal_t * pr, ucons_elm_t * u,
  //                      ucons_set_t * arr);
  //
  //      bool ucons_set_is_top (ucons_internal_t * pr, ucons_set_t * c);
  //
  //      bool ucons_set_is_bottom (ucons_internal_t * pr, ucons_set_t * c);
  //
  //
  //      /* ********************************************************************** */
  //      /* III. pattern_t */see ucons_representation.c
  //      /* ********************************************************************** */
  //

  pattern_key_set_t *pattern_key_set_copy (ucons_internal_t * pr,
                                           pattern_key_set_t * a, size_t size);

  void
  pattern_key_fprint (FILE * stream, ucons_internal_t *pr, pattern_key_t * a,
                      char **name_of_dim);
  void
  ucons_fprint_dcons (FILE * stream, ap_manager_t * man, ucons_t * a,
                      char **name_of_dim, pattern_key_t *key);
  void
  ucons_fprint_econs (FILE * stream, ap_manager_t * man, ucons_t * a,
                      char **name_of_dim);

  //      /* ********************************************************************** */
  //      /* IV. ucons_pattern_t */
  //      /* ********************************************************************** */
  //
  size_t get_pattern_type (ucons_internal_t *pr, size_t u_seg, size_t e_seg, size_t nr_y, pattern_kind name);


  /* used for ucons_expand */
  ucons_t *
  ucons_singleton (ucons_internal_t * pr, bool destructive,
                   ucons_t * a, ap_dim_t dim);
  ucons_t* ucons_split (ucons_internal_t* pr, bool destructive,
                        ucons_t* a, ap_dim_t dim);
  ucons_t * split_with_pattern_P12_P11 (ucons_internal_t* pr, ucons_t *r,
                                        pattern_key_t *pattern_j,
                                        ap_dim_t n1, ap_dim_t n2);

  ucons_t * split_with_pattern_succ_P12 (ucons_internal_t* pr, ucons_t *r,
                                         pattern_key_t *pattern_j,
                                         ap_dim_t n1, ap_dim_t n2);

  ucons_t *split_with_pattern_P21 (ucons_internal_t * pr, ucons_t * r,
                                   pattern_key_t* p,
                                   ap_dim_t node1, ap_dim_t node2);
  /* used for ucons fold */

  ucons_t * fold_with_closure_succ_P12 (ucons_internal_t *pr, ucons_t * a,
                                        ap_dim_t * tdim, size_t size,
                                        bool update_lenght);

  ucons_t * ucons_fold_pattern_1_lx (ucons_internal_t * pr, ucons_t *r, ucons_t * a,
                                     ap_dim_t * tdim, size_t size_tdim,
                                     pattern_key_t * pattern_j);
  ucons_t * ucons_fold_pattern_2_1 (ucons_internal_t * pr, ucons_t *r, ucons_t * a,
                                    ap_dim_t * tdim, size_t size_tdim,
                                    pattern_key_t * pattern_j);
  ucons_t * ucons_change_pattern_2_1 (ucons_internal_t * pr, ucons_t *r, ucons_t * a,
                                      ap_dim_t * tdim, size_t size_tdim,
                                      pattern_key_t * pattern_j);



  //      /* ********************************************************************** */
  //      /* V. ucons_t */see ucons_representation.c
  //      /* ********************************************************************** */

  ucons_t *ucons_alloc_internal (ucons_internal_t * pr, size_t intdim,
                                 size_t realdim);
  ucons_t *ucons_copy_internal (ucons_internal_t * pr, ucons_t * a);
  ucons_t *ucons_alloc_top (ucons_internal_t * pr, size_t intdim, size_t dim);


  void ucons_free_internal (ucons_internal_t * pr, ucons_t * a);

  /*replaced by is_pattern_inst */
  bool test_pattern_sat (ucons_internal_t * pr, ucons_t *r,
                         pattern_key_t *p, int s);

  bool is_pattern_inst (ucons_internal_t * pr, ucons_t *r,
                        pattern_key_t *p);

  void pattern_instance_pattern_1_lx (ucons_internal_t *pr, ucons_t *r);

  void pattern_instance_pattern_1_2 (ucons_internal_t *pr, ucons_t *r);

  void pattern_instance_pattern_2_1 (ucons_internal_t *pr, ucons_t *r);

  void pattern_instance_pattern_1 (ucons_internal_t *pr, ucons_t *r);

  ucons_t *add_pattern_n2p (ucons_internal_t * pr, ucons_t * r, pattern_key_t * key);
  ucons_t *remove_pattern_n2p (ucons_internal_t * pr, ucons_t * r, pattern_key_t * key);
  /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */




#endif /* UCONS_INTERNAL_H_ */
