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


#ifndef __APRON2SHAPE_H_
#define __APRON2SHAPE_H_

#include "shape_manager.h"

/* *INDENT-OFF* */
#ifdef __cplusplus
extern "C"
{
#endif
  /* *INDENT-ON* */


  /* ================================================================== */
  /* Dimensions */
  /* ================================================================== */

#define PTR2REAL_DIM(pr,d)   d
#define REAL2PTR_DIM(pr,d)   d
#define DIM2DATA(d,intdim)  d
#define DIM2PTR(d,intdim)   ((d==NULL_DIM)? NULL_DIM:(d-intdim))

#define IS_INTDIM(d,intdim,ptrdim) (d < intdim)
#define IS_PTRDIM(d,intdim,ptrdim) (intdim <= d && (d-intdim) < ptrdim)

#define IS_NULLDIM(d) (d == AP_DIM_MAX)
#define NULL_DIM      AP_DIM_MAX

#define DATA_DIM(intdim,ptrdim) intdim

  /* Encoding of operations on pointers */
  // with AP_TEXPR_CAST
#define IS_OFFSET(t,d)      \
     ((t==AP_RTYPE_INT || t==AP_RTYPE_REAL) && (d==AP_RDIR_ZERO || d==AP_RDIR_UP || d==AP_RDIR_NEAREST))
#define IS_FIELD_REF(t,d)   ((t==AP_RTYPE_INT || t==AP_RTYPE_REAL) && d==AP_RDIR_ZERO)
#define IS_FIELD_DEREF(t,d) (t==AP_RTYPE_REAL && d==AP_RDIR_ZERO)

#define IS_ARRAY_REF(t,d)   (t==AP_RTYPE_INT && d==AP_RDIR_UP)
#define IS_ARRAY_DEREF(t,d) (t==AP_RTYPE_REAL && d==AP_RDIR_UP)
#define IS_ARRAY_SPLIT(t,d) (t==AP_RTYPE_REAL && d==AP_RDIR_NEAREST)

  // with AP_TEXPR_CAST
#define IS_FUN(t,d)         (t==AP_RTYPE_QUAD)
#define IS_LENGTH_FUN(t,d)  (t==AP_RTYPE_QUAD && d==AP_RDIR_UP)
#define IS_SUM_FUN(t,d)     (t==AP_RTYPE_QUAD && d==AP_RDIR_ZERO)
#define IS_MSET_FUN(t,d)    (t==AP_RTYPE_QUAD && d==AP_RDIR_NEAREST)
#define IS_UCONS_FUN(t,d)   (t==AP_RTYPE_QUAD && d==AP_RDIR_DOWN)

  // with AP_TEXPR_SQRT
#define IS_FREE_STMT(t,d)   (t==AP_RTYPE_SINGLE && d==AP_RDIR_ZERO)
#define IS_ALLOC_STMT(t,d)  (t==AP_RTYPE_SINGLE && d==AP_RDIR_NEAREST)



  /* ================================================================== */
  /* Expressions and constraints */
  /* ================================================================== */

  pcons0_t *shape_pcons_of_lincons (shape_internal_t * pr, ap_lincons0_t * c,
                                    size_t intdim, size_t ptrdim);
  pcons0_t *shape_pcons_of_tcons (shape_internal_t * pr, ap_tcons0_t * c,
                                  size_t intdim, size_t ptrdim);
  pcons0_array_t *shape_pcons_array_of_lincons_array (shape_internal_t * pr,
                                                      ap_lincons0_array_t * c,
                                                      size_t intdim, size_t ptrdim);
  pcons0_array_t *shape_pcons_array_of_tcons_array (shape_internal_t * pr,
                                                    ap_tcons0_array_t * c,
                                                    size_t intdim, size_t ptrdim);
  /*
   * Convert Apron constraints to pointer constraints. In arrays, constraints
   * are sorted as follows: constraints involving NULL < EQ_CONS <
   * REACH*_CONS < NE_CONS < DATA_CONS
   */

  ap_lincons0_t shape_lincons_of_node (shape_internal_t * pr,
                                       ap_lincons0_t * c, int* offsets,
                                       size_t * v2n,
                                       size_t size, size_t intdim, size_t ptrdim);
  /*
   * Convert a linear constraint on ptr vars and data (intdim+ptrdim) to a
   * linear expression on nodes and data (size + intdim) using array v2n of
   * size (ptrdim)
   */


  /* ================================================================== */
  /* Assignments */
  /* ================================================================== */

  passign0_t *shape_passign_of_linexpr (shape_internal_t * pr,
                                        ap_dim_t lhs, ap_linexpr0_t * rhs,
                                        size_t intdim, size_t ptrdim);
  passign0_array_t *shape_passign_of_linexpr_array (shape_internal_t * pr,
                                                    ap_dim_t * lhs, ap_linexpr0_t ** rhs,
                                                    size_t size, size_t intdim, size_t ptrdim);
  passign0_t *shape_passign_of_texpr (shape_internal_t * pr, ap_dim_t lhs,
                                      ap_texpr0_t * rhs, size_t intdim, size_t ptrdim);
  passign0_array_t *
  shape_passign_of_texpr_array (shape_internal_t * pr,
                                ap_dim_t * lhs, ap_texpr0_t ** rhs,
                                size_t size, size_t intdim, size_t ptrdim);
  /*
   * Convert Apron assignments to pointer assignments
   */


  /* ================================================================== */
  /* Apron expressions for Shapes */

  /* ================================================================== */

  typedef enum
  {
    expr_ptr, /* x */
    expr_reach, /* x - y */
    expr_reachl, /* x - y - l */
    expr_lindata, /* data only */
    expr_next, /* x * next */
    expr_deref, /* x / next, x*next / next, */
    expr_data, /* v + x*data */
    expr_deref_data, /* expr_data / data */
  } exprmode_t;

  /* High level constructors */
  ap_linexpr0_t *shape_linexpr_null (size_t intdim, size_t ptrdim);
  ap_texpr0_t *shape_texpr_null (size_t intdim, size_t ptrdim);
  ap_linexpr0_t *shape_linexpr_x (size_t x, size_t intdim, size_t ptrdim);

  bool shape_tcons_array_is_true (ap_tcons0_array_t* array);

  ap_lincons0_t shape_lincons_same_x_y (size_t x, size_t y,
                                        size_t intdim, size_t ptrdim);
  ap_tcons0_t shape_tcons_same_x_y (size_t x, bool nx, size_t y, bool ny,
                                    size_t intdim, size_t ptrdim);

  ap_lincons0_t shape_lincons_diff_x_y (size_t x, size_t y,
                                        size_t intdim, size_t ptrdim);
  ap_tcons0_t shape_tcons_diff_x_y (size_t x, bool nx, size_t y, bool ny,
                                    size_t intdim, size_t ptrdim);

  ap_lincons0_t shape_lincons_reach_x_y (size_t x, size_t y,
                                         size_t intdim, size_t ptrdim);
  ap_tcons0_t shape_tcons_reach_x_y (size_t x, bool nx, size_t y, bool ny,
                                     size_t intdim, size_t ptrdim);

  ap_texpr0_t *shape_texpr_deref_data (ap_texpr0_t* t);
  ap_texpr0_t *shape_texpr_deref_next (ap_texpr0_t* t);
  ap_texpr0_t *shape_texpr_x_alloc (size_t x, size_t intdim, size_t ptrdim);
  ap_texpr0_t *shape_texpr_x_free (size_t x, size_t intdim, size_t ptrdim);
  ap_texpr0_t *shape_texpr_x_data (size_t x, size_t intdim, size_t ptrdim);
  ap_texpr0_t *shape_texpr_x_len (size_t x, size_t intdim, size_t ptrdim);
  ap_texpr0_t *shape_texpr_x_ucons (size_t x, size_t intdim, size_t ptrdim);
  ap_texpr0_t *shape_texpr_x_at_i (size_t x, size_t i, size_t intdim, size_t ptrdim);


  extern const char *exprname[];
#define RANDOM_PTRDIM(pr,intdim,ptrdim) (intdim + lrand48() % ptrdim)
#define RANDOM_PTRDIM_OR_NULL(pr,intdim,ptrdim) ((lrand48() % 2 == 0) ? NULL_DIM : RANDOM_PTRDIM(pr,intdim,ptrdim))
  ap_linexpr0_t *shape_linexpr_random (shape_internal_t * pr, exprmode_t mode,
                                       size_t intdim, size_t ptrdim);
  ap_texpr0_t *shape_texpr_random (shape_internal_t * pr, exprmode_t mode,
                                   size_t intdim, size_t ptrdim);
  /* Build radomly Apron expressions useful for Shapes */

  /* Low level constructors: take care with them! */
  ap_lincons0_t shape_lincons_x_cst (ap_constyp_t op, size_t x, int cst,
                                     size_t intdim, size_t ptrdim);
  ap_linexpr0_t *shape_linexpr_x_y_l (int cx, size_t x, int cy, size_t y,
                                      int cl, size_t l,
                                      size_t intdim, size_t ptrdim);
  ap_lincons0_t shape_lincons_x_y_l (ap_constyp_t op, int cx, size_t x, int cy,
                                     size_t y, int cl, size_t l,
                                     size_t intdim, size_t ptrdim);
  ap_lincons0_t shape_lincons_x_y_v_cst (ap_constyp_t op, size_t code, int cx,
                                         size_t x, int cy, size_t y, int cv,
                                         size_t v, int cst,
                                         size_t intdim, size_t ptrdim);
  ap_texpr0_t *shape_texpr_x_y_v_cst_data (int cx, size_t x, int cy, size_t y,
                                           int cv, size_t v, int cst, bool rhs,
                                           size_t intdim, size_t ptrdim);
  ap_tcons0_t shape_tcons_x_y (ap_constyp_t op, size_t x, bool nextx, size_t y,
                               bool nexty, size_t intdim, size_t ptrdim);
  ap_tcons0_t shape_tcons_x_cst (ap_constyp_t op, size_t x, bool nextx, int cst,
                                 size_t intdim, size_t ptrdim);
  /* Build Apron expressions useful for testing Shapes */

  ap_texpr0_t *shape_texpr_x_n_next (size_t x, int nextx,
                                     size_t intdim, size_t ptrdim);
  ap_texpr0_t *shape_texpr_x_cst (size_t x, bool nextx, int cst,
                                  size_t intdim, size_t ptrdim);
  ap_texpr0_t *shape_texpr_singleton (int n);
  ap_texpr0_t *shape_texpr_split (size_t src, size_t new);
  ap_texpr0_t **shape_texpr_merge (size_t ** nodes, size_t size, size_t max);
  /* Build Apron expressions on nodes for the segment abstract domain */

  ap_linexpr0_t *shape_linexpr_of_node (shape_internal_t * pr,
                                        ap_linexpr0_t * c,
                                        size_t * v2n, size_t size,
                                        size_t intdim, size_t ptrdim);
  /*
   * Convert a linear expression on ptr vars and data (intdim+ptrdim) to a
   * linear expression on nodes and data (size + intdim) using array v2n of
   * size (ptrdim)
   */

  void shape_dimperm_copy (ap_dimperm_t * dst, size_t begin,
                           ap_dimperm_t * src);
  /* Put the permutation of nodes from src in the dst starting with begin */
  void shape_dimchange_sort (ap_dimchange_t * p);
  /* Do insertion sort for dimensions */
  bool shape_dimperm_is_id (ap_dimperm_t * p);
  /* Test if a permutation is identity */

  /* ================================================================= */
  /* Printing through the Apron interface for Shapes */
  /* ================================================================= */

#define MAX_SEGM_DIM  20
  void shape_init_name_of_dim (size_t intdim, size_t segmdim);
  char* shape_name_of_dim (size_t dim);
  bool shape_get_print (void);
  void shape_set_print_dot (void);
  void shape_set_print_smt (void);

  /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */

#endif /* SHAPE_HGRAPH_H_ */
