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



/** Types for Shad: experimental 
 */

#ifndef _SHAD_H_
#define _SHAD_H_

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "ap_environment.h"
#include "ap_dimension.h"
#include "ap_lincons0.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* ====================================================================== */
  /* Datatypes */
  /* ====================================================================== */

  /* An edge formula = pair of nodes */
  typedef struct sh_edgeform_t
  {
    ap_dim_t src; /* Source node given dimension in nodes */
    ap_dim_t dst; /* Destination node, AP_DIM_MAX if nilNode */
  } sh_edgeform_t;

  /** A label formula = pair of ptrvar and node in the
   *  environment given by global env + node env
   */
  typedef struct sh_labelform_t
  {
    ap_dim_t var; /* Program variable in env */
    ap_dim_t node; /* Node, AP_DIM_MAX if nilNode, in nodes */
  } sh_labelform_t;

  /* An elementary term used in data formulas */
  typedef enum
  {
    SH_F_SYMBOL,
    SH_F_LEN,
    SH_F_DATA,
    SH_F_SUM,
    SH_F_MSET,
    /* PENDING: SH_F_ARRAY, */
    SH_F_OTHER /* not to be used */
  } sh_funid_t;

  typedef struct sh_term_t
  {
    sh_funid_t funid; /* index of function used for this term, @see formula */
    ap_dim_t p0; /* first parameter */
    ap_dim_t p1; /* optional 2nd parameter: @see funid */
  } sh_term_t;

  typedef struct sh_linterm_t
  {
    sh_term_t t;
    int coeff;
    struct sh_linterm_t* next;
  } sh_linterm_t;

  /* A data formula is a linear constraint over terms */
  typedef struct sh_dataform_t
  {
    sh_linterm_t* p; /* Sorted list (by t) of linear terms */
    int cst; /* constant */
    ap_constyp_t constyp; /* type of constraint */
  } sh_dataform_t;

  /* A guard is either a fixed constraint of a general data formula */
  typedef enum sh_guardtyp_t
  {
    SH_G_ALL, /* (y,n):           y\in tl(n) */
    SH_G_LE2, /* (y1,y2,n):       y1,y2\in tl(n) y1\le y2 */
    SH_G_SUCC2, /* (y1,y2,n):     y1,y2\in tl(n) y1+1\eq y2 */
    SH_G_FST, /* (y,n):           y\in tl(n) y\eq 1 */
    SH_G_LST, /* (y,n):           y\in tl(n) y\eq len(n)-1 */
    SH_G_EQ2, /* (y1,n1,y2,n2):   y1\in tl(n1), y2\in tl(n2) y1\eq y2 */
    SH_G_SL2, /* (y1,n1,y2,n2,c): y1\in tl(n1), y2\in tl(n2) y1\eq y2+c */
    SH_G_SR2, /* (y1,n1,y2,n2,c): y1\in tl(n1), y2\in tl(n2) y1\eq y2+c */
    SH_G_OTHER
  } sh_guardtyp_t;

  typedef struct sh_guardform_t
  {
    sh_guardtyp_t guardtyp;
    ap_dim_t * y; /* position variables */
    ap_dim_t * n; /* corresponding node variables */
    size_t length_y; /* length of y and n */
    sh_dataform_t* data; /* for SH_G_OTHER, else NULL */
    size_t length_data; /* for SH_G_OTHER, else NULL */
  } sh_guardform_t;

  /* An universal formula of the form forall y. guard => data */
  typedef struct sh_univform_t
  {
    ap_var_t* y; /* Universally quantified variables */
    size_t length_y;
    sh_guardform_t guard; /* Guard */
    sh_dataform_t* data; /* Data constraint over env, nodes, y */
    size_t length_data;
  } sh_univform_t;

  /* A SHAD formula = graph + pointer + word constraints */
  typedef struct sh_shadform_t
  {
    ap_environment_t* nodes; /* Declared nodes, nilNode is not there */
    sh_edgeform_t* eform; /* Graph edges over nodes */
    size_t length_eform;
    sh_labelform_t* pform; /* Pointer formula over env + nodes */
    size_t length_pform; /* may be less than ptrdim */
    sh_dataform_t* dform; /* Data formula over env(intdim) + nodes */
    size_t length_dform;
    sh_univform_t* uform; /* Segment formula over env(intdim) + nodes */
    size_t length_uform;
  } sh_shadform_t;

  /* General SHAD formulas */
  typedef struct sh_formula_t
  {
    ap_environment_t* env; /* Declared ptr and data variables */
    int dw; /* Data words logics used, a bitset, @see   */
    sh_shadform_t** form; /* Array of formulas, disjunct of its elements */
    size_t size; /* Array size */
  } sh_formula_t;

  /* ====================================================================== */
  /* Constructors/destructors */
  /* ====================================================================== */

  extern sh_formula_t* sh_pos;
  extern sh_formula_t* sh_neg;
  extern sh_formula_t* sh_crt;
  extern int sh_guards;
  /* Variables for building values in the abstract domains */

  void sh_init (void);
  /* Initializes the global vars above */
  void sh_push2 (sh_formula_t* a, sh_formula_t* b);
  void sh_push (sh_formula_t* f, bool sign);
  /* Pushes the formula in the global formula given by sign */
  void sh_formula_free (sh_formula_t* a);
  /* Free the environment and all disjuncts */

  sh_linterm_t* sh_linterm_alloc_symbol (ap_dim_t v);
  sh_linterm_t* sh_linterm_alloc_data (ap_dim_t n, ap_dim_t y);
  sh_linterm_t* sh_linterm_alloc_len (ap_dim_t v);
  sh_linterm_t* sh_linterm_alloc_sum (ap_dim_t v);
  sh_linterm_t* sh_linterm_alloc_mset (ap_dim_t v);
  void sh_linterm_free (sh_linterm_t* a);
  /* Free the list of linterms. */
  sh_linterm_t* sh_linterm_merge (sh_linterm_t* a, int coeffa, sh_linterm_t* b, int coeffb);
  /* Computes coeffa * a + coeffb * b in a and free common terms */

  /* ====================================================================== */
  /* Testing */
  /* ====================================================================== */

  int sh_check (void);
  /* Checks implication sh_pos => sh_neg
   * corresponding to satisfaction of sh_pos /\ not(sh_neg)
   */

  /* ====================================================================== */
  /* Printing and scaning */
  /* ====================================================================== */

  void sh_env_fprint (FILE* stream, ap_environment_t* env);
  void sh_fprint (FILE* stream, sh_formula_t * f);
  /*
   * Print the formula in a pretty way.
   */
  int sh_fscan (char* filename);
  /* 
   * Read the formulas above from the file filename.
   */


 
#ifdef __cplusplus
}
#endif

#endif

