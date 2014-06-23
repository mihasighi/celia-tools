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



#ifndef _SL3SHAD_H_
#define _SL3SHAD_H_

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "sl3.h"
#include "shad.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* ======================================================================
   * Translation to shad
   * ====================================================================== */

  sh_formula_t* sl3shad (sl3_context_t ctx, sl3_expr_t term, bool sign);

  sh_formula_t* sl3shad_or (sl3_context_t ctx, sl3_expr_t term,
                            ap_var_t* qarr, size_t qsize,
                            sh_formula_t* f, size_t *n);
  size_t sl3shad_or_size (sl3_expr_t term);

  sh_formula_t* sl3shad_exists (sl3_context_t ctx, sl3_expr_t term,
                                ap_var_t* qarr, size_t qsize,
                                sh_formula_t* f, size_t n);

  void sl3shad_and (sl3_context_t ctx, sl3_expr_t term,
                    sh_shadform_t* shf);
  void sl3shad_and_size (sl3_expr_t term,
                         size_t *nls, size_t* nlab, size_t* ndt, size_t *nfor);
  void sl3shad_and_aux (sl3_context_t ctx, sl3_expr_t term,
                        sh_shadform_t* shf,
                        size_t *nls, size_t* nlab, size_t* ndt, size_t *nfor);

  void sl3shad_edge (sl3_context_t ctx, sl3_expr_t term,
                     sh_shadform_t* shf, size_t* n);
  void sl3shad_label (sl3_context_t ctx, sl3_expr_t term,
                      sh_shadform_t* shf, size_t* n);

  void sl3shad_data (ap_environment_t* env, sl3_expr_t term,
                     sh_dataform_t* df, size_t len, size_t* n);
  sh_linterm_t* sl3shad_linterm (ap_environment_t* env, sl3_expr_t term,
                                 int *cst);

  void sl3shad_univ (ap_environment_t* env, sl3_expr_t term,
                     sh_univform_t* uf, size_t len, size_t* n);
  void sl3shad_univ_data (ap_environment_t* env, sl3_expr_t term,
                          sh_dataform_t* uf, size_t len, size_t* n);
  void sl3shad_guard (ap_environment_t* env, sl3_expr_t term,
                      sh_guardform_t* gf);

#ifdef __cplusplus
}
#endif

#endif

