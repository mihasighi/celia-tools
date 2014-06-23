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



/** Types for SL3: experimental 
 */

#ifndef _SL3_H_
#define _SL3_H_

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "ap_var.h"
#include "ap_environment.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* ====================================================================== */
  /* Datatypes */
  /* ====================================================================== */

  /* AST */

  /* Logics supported in SL3 */
  typedef enum
  {
    SL3_LOG_SL3 = 0, /* general: graph + segments */
    SL3_LOG_AUFLIA, /* only on segments */
    SL3_LOG_OTHER /* NOT TO BE USED */
  } sl3_logic_t;

  /** Basic types used in SL3 */
  typedef enum
  {
    SL3_TYP_BOOL,
    SL3_TYP_INT,
    SL3_TYP_NODE,
    SL3_TYP_PTR,
    SL3_TYP_OTHER
  } sl3_typ_t;

  /** Type expressions used in SL3 are
   *  - basic type
   *  - tuple of basic types
   *  - function with basic type as codomain
   */
  typedef struct sl3_type_t
  {
    sl3_typ_t* arr; // used in tuples and fun types
    size_t size;
    sl3_typ_t ty; // used in basic and fun types
  } *sl3_type_t;

  /* Valid term builder in SL3 */
  typedef enum
  {
    SL3_F_FALSE = 0,
    SL3_F_TRUE,
    SL3_F_OR,
    SL3_F_AND,
    SL3_F_NOT,
    SL3_F_IMPLIES,
    SL3_F_NUMBER,
    SL3_F_LT,
    SL3_F_LE,
    SL3_F_GT,
    SL3_F_GE,
    SL3_F_EQ,
    SL3_F_PLUS,
    SL3_F_MINUS,
    SL3_F_TIMES,
    SL3_F_SYMBOL,
    SL3_F_FORALL,
    SL3_F_EXISTS,
    SL3_F_DATA,
    SL3_F_SELECT,
    SL3_F_LEN,
    SL3_F_SUM,
    SL3_F_MSET,
    SL3_F_NIL,
    SL3_F_LS,
    SL3_F_LABEL,
    SL3_F_GALL,
    SL3_F_GEQ2,
    SL3_F_GLE2,
    SL3_F_GSUCC2,
    SL3_F_GFST,
    SL3_F_GLST,
    SL3_F_GSL2,
    SL3_F_GSR2,
    SL3_F_OTHER
  } sl3_funkind_t;

  typedef struct sl3_expr_t
  {
    sl3_funkind_t discr;

    union
    {

      /* normal function */
      struct
      {
        struct sl3_expr_t** arr;
        size_t size;
      } args;

      /* quantifiers */
      struct
      {
        ap_var_t* qarr;
        size_t qsize;
        struct sl3_expr_t* arg;
      } quant;
      /* numbers */
      int ivalue;
      /* symbols */
      char* name;
    } p;
  } * sl3_expr_t;

  /* Context used to parse smtlib2 formulas */
  typedef struct sl3_context_t
  {
    /* logic of work */
    sl3_logic_t logic;
    /* bitset of guards, see ucons_internal.h */
    int guards_bs;
    /* global variables */
    ap_environment_t* g_env;
    /* quantifier level */
    size_t qstack;
    /* local variables */
    ap_environment_t* l_env;
  } * sl3_context_t;


  /* ====================================================================== */
  /* Constructors/destructors */
  /* ====================================================================== */

  /* Parsing context */
  sl3_context_t sl3_mk_context (void);
  void sl3_del_context (sl3_context_t ctx);

  /* Parsing logic */
  sl3_logic_t sl3_set_logic (sl3_context_t ctx, const char* logic);

  /* Types */
  sl3_type_t sl3_mk_type (sl3_context_t ctx, const char* name, int arity);
  sl3_type_t sl3_mk_type_lst (sl3_context_t ctx, const char* name,
                              sl3_type_t* ty, size_t size);
  void sl3_del_type (sl3_type_t ty);
  void sl3_del_type_lst (sl3_type_t* ty, size_t size);

  /* Functions */
  sl3_type_t sl3_mk_fun_type (sl3_context_t ctx,
                              sl3_type_t* domain, size_t domain_size, sl3_type_t range);
  sl3_type_t sl3_mk_fun_decl (sl3_context_t ctx,
                              const char* name, sl3_type_t ty);
  sl3_type_t sl3_mk_fun_data (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_select (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_len (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_fsum (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_fmset (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_nilnode (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_ls (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_sep (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_label (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_gall (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_gle2 (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_gsucc2 (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_gfst (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_glst (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_geq2 (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_gsr2 (sl3_type_t ty);
  sl3_type_t sl3_mk_fun_gsl2 (sl3_type_t ty);

  /* Commands */
  bool sl3_assert (sl3_context_t ctx, sl3_expr_t term);
  int sl3_check (sl3_context_t ctx);

  /* Terms */
  void sl3_push_var (sl3_context_t ctx, const char* name, sl3_typ_t vty);
  bool sl3_push_quant (sl3_context_t ctx);
  bool sl3_pop_quant (sl3_context_t ctx);
  sl3_expr_t sl3_mk_forall (sl3_context_t ctx, sl3_expr_t term);
  sl3_expr_t sl3_mk_exists (sl3_context_t ctx, sl3_expr_t term);
  sl3_expr_t sl3_mk_num_from_string (sl3_context_t ctx, const char* num);
  sl3_expr_t sl3_mk_app (sl3_context_t ctx, const char* name,
                         sl3_type_t ty, sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_var (sl3_context_t ctx, const char* name, sl3_type_t ty);
  sl3_expr_t sl3_mk_true (sl3_context_t ctx);
  sl3_expr_t sl3_mk_false (sl3_context_t ctx);
  sl3_expr_t sl3_mk_and (sl3_context_t ctx,
                         sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_or (sl3_context_t ctx,
                        sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_not (sl3_context_t ctx,
                         sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_implies (sl3_context_t ctx,
                             sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_num (sl3_context_t ctx, int n);
  sl3_expr_t sl3_mk_sum (sl3_context_t ctx,
                         sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_sub (sl3_context_t ctx,
                         sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_mul (sl3_context_t ctx,
                         sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_eq (sl3_context_t ctx,
                        sl3_expr_t arg1, sl3_expr_t arg2);
  sl3_expr_t sl3_mk_le (sl3_context_t ctx,
                        sl3_expr_t arg1, sl3_expr_t arg2);
  sl3_expr_t sl3_mk_lt (sl3_context_t ctx,
                        sl3_expr_t arg1, sl3_expr_t arg2);
  sl3_expr_t sl3_mk_ge (sl3_context_t ctx,
                        sl3_expr_t arg1, sl3_expr_t arg2);
  sl3_expr_t sl3_mk_gt (sl3_context_t ctx,
                        sl3_expr_t arg1, sl3_expr_t arg2);
  sl3_expr_t sl3_mk_data (sl3_context_t ctx,
                          sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_select (sl3_context_t ctx,
                            sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_len (sl3_context_t ctx,
                         sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_fsum (sl3_context_t ctx,
                          sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_fmset (sl3_context_t ctx,
                           sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_nilnode (sl3_context_t ctx);
  sl3_expr_t sl3_mk_ls (sl3_context_t ctx,
                        sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_label (sl3_context_t ctx,
                           sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_gall (sl3_context_t ctx,
                          sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_gle2 (sl3_context_t ctx,
                          sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_gsucc2 (sl3_context_t ctx,
                            sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_gfst (sl3_context_t ctx,
                          sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_glst (sl3_context_t ctx,
                          sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_geq2 (sl3_context_t ctx,
                          sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_gsr2 (sl3_context_t ctx,
                          sl3_expr_t args[], size_t size);
  sl3_expr_t sl3_mk_gsl2 (sl3_context_t ctx,
                          sl3_expr_t args[], size_t size);

#ifdef __cplusplus
}
#endif

#endif
