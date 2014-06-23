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


#include <stdlib.h>
#include "sl3.h"
#include "sl3shad.h"

/* ======================================================================
 * Messages
 * ====================================================================== */

void
sl3_error (int level, char* fun, char* msg)
{
  fprintf (stderr, "SL3 Error of level %d in %s: %s.\n", level, fun, msg);
  if (level == 0)
    // terminate
    exit (0);
}

void
sl3_error_args (int level, char* fun, size_t size)
{
  fprintf (stderr, "SL3 Error of level %d in %s: bad number (%zu) of arguments.\n",
           level, fun, size);
  if (level == 0)
    // terminate
    exit (0);
}

void
sl3_error_id (int level, char* fun, const char* name)
{
  fprintf (stderr, "SL3 Error of level %d in %s: identifier %s not declared.\n",
           level, fun, name);
  if (level == 0)
    // terminate
    exit (0);
}

/* ======================================================================
 * Context
 * ====================================================================== */

sl3_context_t
sl3_mk_context (void)
{
  sl3_context_t r = (sl3_context_t) malloc (sizeof (struct sl3_context_t));
  r->g_env = NULL;
  r->guards_bs = 0;
  r->qstack = 0;
  r->l_env = NULL;
  sh_init ();
  return r;
}

void
sl3_del_context (sl3_context_t ctx)
{
  if (ctx->g_env) ap_environment_free (ctx->g_env);
  if (ctx->l_env) ap_environment_free (ctx->l_env);
  free (ctx);
}

/* ====================================================================== 
 * Logic
 * ====================================================================== */

/** Checks that the logic is supported and set it in the context.
 */
sl3_logic_t
sl3_set_logic (sl3_context_t ctx, const char* logic)
{
  if (strcmp (logic, "SL3") == 0)
    ctx->logic = SL3_LOG_SL3;
  else if (strcmp (logic, "AUFLIA") == 0)
    ctx->logic = SL3_LOG_AUFLIA;
  else
    ctx->logic = SL3_LOG_OTHER;
  return ctx->logic;
}

/* ====================================================================== 
 * Types 
 * ====================================================================== */

sl3_type_t
sl3_new_type (sl3_typ_t v)
{
  sl3_type_t res = (sl3_type_t) malloc (sizeof (struct sl3_type_t));
  res->arr = NULL;
  res->size = 0;
  res->ty = v;
  return res;
}

sl3_type_t
sl3_new_type_2 (sl3_typ_t v1, sl3_typ_t v2)
{
  sl3_type_t res = (sl3_type_t) malloc (sizeof (struct sl3_type_t));
  res->arr = (sl3_typ_t*) malloc (2 * sizeof (sl3_typ_t));
  res->size = 2;
  res->arr[0] = v1;
  res->arr[1] = v2;
  res->ty = SL3_TYP_OTHER;
  return res;
}

sl3_type_t
sl3_new_type_3 (sl3_typ_t v1, sl3_typ_t v2, sl3_typ_t v3)
{
  sl3_type_t res = (sl3_type_t) malloc (sizeof (struct sl3_type_t));
  res->arr = (sl3_typ_t*) malloc (3 * sizeof (sl3_typ_t));
  res->size = 3;
  res->arr[0] = v1;
  res->arr[1] = v2;
  res->arr[2] = v3;
  res->ty = SL3_TYP_OTHER;
  return res;
}

sl3_type_t
sl3_new_type_5 (sl3_typ_t v1, sl3_typ_t v2, sl3_typ_t v3,
                sl3_typ_t v4, sl3_typ_t v5)
{
  sl3_type_t res = (sl3_type_t) malloc (sizeof (struct sl3_type_t));
  res->arr = (sl3_typ_t*) malloc (5 * sizeof (sl3_typ_t));
  res->size = 5;
  res->arr[0] = v1;
  res->arr[1] = v2;
  res->arr[2] = v3;
  res->arr[3] = v4;
  res->arr[4] = v5;
  res->ty = SL3_TYP_OTHER;
  return res;
}

sl3_type_t
sl3_mk_type (sl3_context_t ctx, const char* name, int arity)
{
  sl3_type_t res = NULL;
  if (arity == 0)
    {
      if (strcmp (name, "Bool") == 0)
        res = sl3_new_type (SL3_TYP_BOOL);
      else if (strcmp (name, "Int") == 0)
        res = sl3_new_type (SL3_TYP_INT);
      else if (ctx->logic == SL3_LOG_SL3
               && strcmp (name, "Ptr") == 0)
        res = sl3_new_type (SL3_TYP_PTR);
      else if (strcmp (name, "Node") == 0)
        res = sl3_new_type (SL3_TYP_NODE);
    }
  return res;
}

sl3_type_t
sl3_mk_type_lst (sl3_context_t ctx,
                 const char* name, sl3_type_t* ty, size_t size)
{
  sl3_type_t res = NULL;
  switch (size)
    {
    case 0:
      res = sl3_mk_type (ctx, name, 0);
      break;

    case 1:
      {
        // may be (Bool Bool), (Int Int), (Node Int), (Ptr Node), (Node Node),
        // (Int Node)
        sl3_type_t ty0 = sl3_mk_type (ctx, name, 0);
        sl3_type_t ty1 = ty[0];
        if (ty0 && ty1 &&
            !ty0->arr && !ty0->size &&
            !ty1->arr && !ty1->size &&
            ((ty0->ty == SL3_TYP_BOOL && ty1->ty == SL3_TYP_BOOL) ||
             (ty0->ty == SL3_TYP_INT && ty1->ty == SL3_TYP_INT) ||
             (ty0->ty == SL3_TYP_NODE && ty1->ty == SL3_TYP_INT) ||
             (ty0->ty == SL3_TYP_PTR && ty1->ty == SL3_TYP_NODE) ||
             (ty0->ty == SL3_TYP_INT && ty1->ty == SL3_TYP_NODE)))
          {
            res = sl3_new_type_2 (ty0->ty, ty1->ty);
          }
        break;
      }
    case 2:
      {
        // may be (Array Int Int), (Int Int Node)
        sl3_type_t ty0 = sl3_mk_type (ctx, name, 0);
        sl3_type_t ty1 = ty[0];
        sl3_type_t ty2 = ty[1];
        if (ctx->logic == SL3_LOG_AUFLIA &&
            !strcmp (name, "Array") &&
            !ty0 && ty1 && ty2 &&
            !ty1->arr && !ty1->size &&
            !ty2->arr && !ty2->size &&
            ty0->ty == SL3_TYP_INT && ty1->ty == SL3_TYP_INT)
          res = sl3_new_type (SL3_TYP_NODE);
        else if (ty0 && ty1 && ty2 &&
                 !ty0->arr && !ty0->size && ty0->ty == SL3_TYP_INT &&
                 !ty1->arr && !ty1->size && ty1->ty == SL3_TYP_INT &&
                 !ty2->arr && !ty2->size && ty2->ty == SL3_TYP_NODE)
          {
            res = sl3_new_type_3 (ty0->ty, ty1->ty, ty2->ty);
          }
        break;
      }
    case 4:
      {
        // may be (Int Node Int Node Int)
        sl3_type_t ty0 = sl3_mk_type (ctx, name, 0);
        sl3_type_t ty1 = ty[0];
        sl3_type_t ty2 = ty[1];
        sl3_type_t ty3 = ty[2];
        sl3_type_t ty4 = ty[3];
        if (ty0 && ty1 && ty2 && ty3 && ty4 &&
            !ty0->arr && !ty0->size && ty0->ty == SL3_TYP_INT &&
            !ty1->arr && !ty1->size && ty1->ty == SL3_TYP_NODE &&
            !ty2->arr && !ty2->size && ty2->ty == SL3_TYP_INT &&
            !ty3->arr && !ty3->size && ty3->ty == SL3_TYP_NODE &&
            !ty4->arr && !ty4->size && ty4->ty == SL3_TYP_INT)
          {
            res = sl3_new_type_5 (ty0->ty, ty1->ty, ty2->ty, ty3->ty, ty4->ty);
          }
        break;
      }
    default: break;
    }
  return res;
}

void
sl3_del_type (sl3_type_t ty)
{
  if (ty && ty->arr && ty->size > 0)
    {
      free (ty->arr);
      ty->arr = NULL;
    }
  free (ty);
}

void
sl3_del_type_lst (sl3_type_t* ty, size_t size)
{
  if (ty)
    {
      size_t i;
      for (i = 0; i < size; i++)
        {
          sl3_del_type (ty[i]);
          ty[i] = NULL;
        }
      free (ty);
    }
}

/* ====================================================================== 
 * Functions
 * ====================================================================== */

sl3_type_t
sl3_mk_fun_type (sl3_context_t ctx,
                 sl3_type_t* domain, size_t domain_size, sl3_type_t range)
{
  sl3_type_t res = NULL;
  size_t i;
  // the range shall be a basic type
  if (range && !range->arr && !range->size && range->ty != SL3_TYP_OTHER)
    {
      if (domain && domain_size > 0)
        { // the domain shall be built only from basic types
          for (i = 0; i < domain_size; i++)
            if (domain[i] && domain[i]->arr)
              return res;
          res = (sl3_type_t) malloc (sizeof (struct sl3_type_t));
          res->arr = (sl3_typ_t *) malloc (domain_size * sizeof (sl3_typ_t));
          for (i = 0; i < domain_size; i++)
            res->arr[i] = domain[i]->ty;
          res->size = domain_size;
          res->ty = range->ty;
        }
      else
        {
          res = sl3_new_type (range->ty);
        }
    }
  return res;
}

sl3_type_t
sl3_mk_fun_decl (sl3_context_t ctx, const char* name, sl3_type_t ty)
{
  // recognize pre-defined functions for SL3
  if (strcmp (name, "nilNode") == 0)
    return sl3_mk_fun_nilnode (ty);
  if (strcmp (name, "len") == 0)
    return sl3_mk_fun_len (ty);
  if (strcmp (name, "data") == 0)
    return sl3_mk_fun_data (ty);
  if (strcmp (name, "sum") == 0)
    return sl3_mk_fun_fsum (ty);
  if (strcmp (name, "mset") == 0)
    return sl3_mk_fun_fmset (ty);
  if (strcmp (name, "sep") == 0)
    return sl3_mk_fun_sep (ty);
  if (strcmp (name, "ls") == 0)
    return sl3_mk_fun_ls (ty);
  if (strcmp (name, "label") == 0)
    return sl3_mk_fun_label (ty);
  if (strcmp (name, "Gall") == 0)
    return sl3_mk_fun_gall (ty);
  if (strcmp (name, "Gle2") == 0)
    return sl3_mk_fun_gle2 (ty);
  if (strcmp (name, "Gsucc2") == 0)
    return sl3_mk_fun_gsucc2 (ty);
  if (strcmp (name, "Gfst") == 0)
    return sl3_mk_fun_gfst (ty);
  if (strcmp (name, "Glst") == 0)
    return sl3_mk_fun_glst (ty);
  if (strcmp (name, "Geq2") == 0)
    return sl3_mk_fun_geq2 (ty);
  if (strcmp (name, "Gsr2") == 0)
    return sl3_mk_fun_gsr2 (ty);
  if (strcmp (name, "Gsl2") == 0)
    return sl3_mk_fun_gsl2 (ty);
  // other function declaration that those are global variables
  // check that the type is basic
  if (ty && ty->size == 0 &&
      (ctx->logic != SL3_LOG_SL3 ||
       ty->ty == SL3_TYP_PTR ||
       ty->ty == SL3_TYP_INT) &&
      (ctx->logic != SL3_LOG_AUFLIA ||
       ty->ty == SL3_TYP_NODE ||
       ty->ty == SL3_TYP_INT))
    {
      // check that the variable is not already in the environment
      ap_environment_t* nenv;
      ap_var_t v = (ap_var_t) ap_var_operations->copy((ap_var_t) name);
      // (ap_var_t) strdup(name);
      if (!ctx->g_env)
        ctx->g_env = ap_environment_alloc_empty ();
      if (ap_environment_mem_var (ctx->g_env, v))
        /* error */ return NULL;
      if (ty->ty == SL3_TYP_INT)
        nenv = ap_environment_add (ctx->g_env, &v, 1, NULL, 0);
      else
        nenv = ap_environment_add (ctx->g_env, NULL, 0, &v, 1);
      ap_environment_free (ctx->g_env);
      ap_var_operations->free (v);
      ctx->g_env = nenv;
      return ty;
    }
  return NULL;
}

sl3_type_t
sl3_mk_fun_data (sl3_type_t ty)
{
  if (ty && ty->size == 2 && ty->ty == SL3_TYP_INT &&
      ty->arr[0] == SL3_TYP_NODE &&
      ty->arr[1] == SL3_TYP_INT)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_select (sl3_type_t ty)
{
  if (ty && ty->size == 2 && ty->ty == SL3_TYP_INT &&
      ty->arr[0] == SL3_TYP_NODE &&
      ty->arr[1] == SL3_TYP_INT)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_len (sl3_type_t ty)
{
  if (ty && ty->size == 1 && ty->ty == SL3_TYP_INT &&
      ty->arr[0] == SL3_TYP_NODE)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_fsum (sl3_type_t ty)
{
  if (ty && ty->size == 1 && ty->ty == SL3_TYP_INT &&
      ty->arr[0] == SL3_TYP_NODE)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_fmset (sl3_type_t ty)
{
  if (ty && ty->size == 1 && ty->ty == SL3_TYP_INT &&
      ty->arr[0] == SL3_TYP_NODE)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_nilnode (sl3_type_t ty)
{
  if (ty && ty->size == 0 && ty->ty == SL3_TYP_NODE)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_ls (sl3_type_t ty)
{
  if (ty && ty->size == 2 && ty->ty == SL3_TYP_BOOL &&
      ty->arr[0] == SL3_TYP_NODE &&
      ty->arr[1] == SL3_TYP_NODE)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_sep (sl3_type_t ty)
{
  if (ty && ty->size == 2 && ty->ty == SL3_TYP_BOOL &&
      ty->arr[0] == SL3_TYP_BOOL &&
      ty->arr[1] == SL3_TYP_BOOL)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_label (sl3_type_t ty)
{
  if (ty && ty->size == 2 && ty->ty == SL3_TYP_BOOL &&
      ty->arr[0] == SL3_TYP_PTR &&
      ty->arr[1] == SL3_TYP_NODE)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_gall (sl3_type_t ty)
{
  if (ty && ty->size == 2 && ty->ty == SL3_TYP_BOOL &&
      ty->arr[0] == SL3_TYP_INT &&
      ty->arr[1] == SL3_TYP_NODE)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_gle2 (sl3_type_t ty)
{
  if (ty && ty->size == 3 && ty->ty == SL3_TYP_BOOL &&
      ty->arr[0] == SL3_TYP_INT &&
      ty->arr[1] == SL3_TYP_INT &&
      ty->arr[2] == SL3_TYP_NODE)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_gsucc2 (sl3_type_t ty)
{
  if (ty && ty->size == 3 && ty->ty == SL3_TYP_BOOL &&
      ty->arr[0] == SL3_TYP_INT &&
      ty->arr[1] == SL3_TYP_INT &&
      ty->arr[2] == SL3_TYP_NODE)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_gfst (sl3_type_t ty)
{
  if (ty && ty->size == 2 && ty->ty == SL3_TYP_BOOL &&
      ty->arr[0] == SL3_TYP_INT &&
      ty->arr[1] == SL3_TYP_NODE)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_glst (sl3_type_t ty)
{
  if (ty && ty->size == 2 && ty->ty == SL3_TYP_BOOL &&
      ty->arr[0] == SL3_TYP_INT &&
      ty->arr[1] == SL3_TYP_NODE)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_geq2 (sl3_type_t ty)
{
  if (ty && ty->size == 4 && ty->ty == SL3_TYP_BOOL &&
      ty->arr[0] == SL3_TYP_INT &&
      ty->arr[1] == SL3_TYP_NODE &&
      ty->arr[2] == SL3_TYP_INT &&
      ty->arr[3] == SL3_TYP_NODE)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_gsr2 (sl3_type_t ty)
{
  if (ty && ty->size == 5 && ty->ty == SL3_TYP_BOOL &&
      ty->arr[0] == SL3_TYP_INT &&
      ty->arr[1] == SL3_TYP_NODE &&
      ty->arr[2] == SL3_TYP_INT &&
      ty->arr[3] == SL3_TYP_NODE &&
      ty->arr[4] == SL3_TYP_INT)
    return ty;
  return NULL;
}

sl3_type_t
sl3_mk_fun_gsl2 (sl3_type_t ty)
{
  if (ty && ty->size == 5 && ty->ty == SL3_TYP_BOOL &&
      ty->arr[0] == SL3_TYP_INT &&
      ty->arr[1] == SL3_TYP_NODE &&
      ty->arr[2] == SL3_TYP_INT &&
      ty->arr[3] == SL3_TYP_NODE &&
      ty->arr[4] == SL3_TYP_INT)
    return ty;
  return NULL;
}

/* ====================================================================== 
 * Commands
 * ====================================================================== */

bool
sl3_assert (sl3_context_t ctx, sl3_expr_t term)
{
  sh_formula_t* r = NULL;
  /* check that l_env is empty (due to push and pop). */
  if (ctx->l_env != NULL && 
      (ctx->l_env->intdim != 0 || ctx->l_env->realdim != 0))
    {
#ifndef NDEBUG
      fprintf (stdout, "Sl3_assert: local environment: ");
      ap_environment_fdump (stdout, ctx->l_env);
      fflush (stdout);
#endif
      sl3_error (1, "sl3_assert", "non empty local environment");
    }
  else if (!term)
    {
      sl3_error (1, "sl3_assert", "empty term");
    }
  else if (term->discr <= SL3_F_AND || term->discr == SL3_F_EXISTS)
    {
      r = sl3shad (ctx, term, true);
      sh_push (r, true);
    }
  else if (term->discr == SL3_F_NOT)
    {
      r = sl3shad (ctx, term, false); // remove NOT
      // r is term without NOT
      sh_push (r, false);
    }
  else
    {
      sl3_error (1, "sl3_assert", "bad formula");
    }
  if (r)
    {
      r->dw = ctx->logic;
#ifndef NDEBUG1
      fprintf (stdout, "Sl3_assert: sh_formula result\n");
      sh_fprint (stdout, r);
      fflush (stdout);
#endif
      // TODO: free r
      return true;
    }
  return false;
}

int
sl3_check (sl3_context_t ctx)
{
  return sh_check ();
}

/* ======================================================================
 * Terms
 * ====================================================================== */

void
sl3_push_var (sl3_context_t ctx, const char* name, sl3_typ_t vty)
{
  ap_var_t v = (ap_var_t) ap_var_operations->copy((ap_var_t) name);
  ap_environment_t* nenv;
  if (vty == SL3_TYP_INT)
    nenv = ap_environment_add (ctx->l_env, &v, 1, NULL, 0);
  else
    nenv = ap_environment_add (ctx->l_env, NULL, 0, &v, 1);
  ap_environment_free (ctx->l_env);
  ap_var_operations->free (v);
  ctx->l_env = nenv;
}

bool
sl3_push_quant (sl3_context_t ctx)
{
  // create a lists of quantifiers iff it is NULL
  if (!ctx->l_env)
    ctx->l_env = ap_environment_alloc_empty ();
  ctx->qstack++;
  return true;
}

bool
sl3_pop_quant (sl3_context_t ctx)
{
  if ((ctx->logic == SL3_LOG_SL3 && ctx->qstack == 2) ||
      (ctx->logic == SL3_LOG_AUFLIA && ctx->qstack == 1))
    {
      // remove all integer variables from the local environment
      ap_environment_t* nenv = ap_environment_remove (ctx->l_env,
                                                      ctx->l_env->var_of_dim,
                                                      ctx->l_env->intdim);
      if (ctx->l_env) ap_environment_free (ctx->l_env);
      ctx->l_env = nenv;
    }
  if (ctx->logic == SL3_LOG_SL3 && ctx->qstack == 1)
    {
      // remove all variables from the environment
      if (ctx->l_env) ap_environment_free (ctx->l_env);
      ctx->l_env = NULL;
    }
  ctx->qstack--;
  return true;
}

sl3_expr_t
sl3_mk_op (sl3_funkind_t f, sl3_expr_t args[], size_t size)
{
  size_t i;
  sl3_expr_t res = (sl3_expr_t) malloc (sizeof (struct sl3_expr_t));
  res->discr = f;
  res->p.args.size = size;
  res->p.args.arr = NULL;
  if (size)
    {
      res->p.args.arr = (sl3_expr_t*) malloc (size * sizeof (struct sl3_expr_t*));
      for (i = 0; i < size; i++)
        res->p.args.arr[i] = args[i];
    }
  return res;
}

sl3_expr_t
sl3_mk_forall (sl3_context_t ctx, sl3_expr_t term)
{
  // extract the int variables from the local env to put as quantifiers
  if (ctx->l_env->intdim > 0)
    {
      size_t i;
      sl3_expr_t res = sl3_mk_op (SL3_F_FORALL, NULL, 0);
      res->p.quant.qsize = ctx->l_env->intdim;
      res->p.quant.qarr = NULL;
      res->p.quant.arg = term;
      res->p.quant.qarr = (ap_var_t*) malloc (res->p.quant.qsize *
                                              sizeof (ap_var_t));
      for (i = 0; i < res->p.quant.qsize; i++)
        res->p.quant.qarr[i] =
              ap_var_operations->copy (ctx->l_env->var_of_dim[i]);
      return res;
    }
  else
    return term;
}

sl3_expr_t
sl3_mk_exists (sl3_context_t ctx, sl3_expr_t term)
{
  // extract the real variables from the local env to put as quantifiers
  if (ctx->l_env->realdim > 0)
    {
      size_t i;
      size_t intdim = ctx->l_env->intdim;
      sl3_expr_t res = sl3_mk_op (SL3_F_EXISTS, NULL, 0);
      res->p.quant.qsize = ctx->l_env->realdim;
      res->p.quant.qarr = NULL;
      res->p.quant.arg = term;
      res->p.quant.qarr = (ap_var_t*) malloc (res->p.quant.qsize *
                                              sizeof (ap_var_t));
      for (i = 0; i < res->p.quant.qsize; i++)
        res->p.quant.qarr[i] =
              ap_var_operations->copy (ctx->l_env->var_of_dim[intdim + i]);
      return res;
    }
  else
    return term;
}

sl3_expr_t
sl3_mk_num_from_string (sl3_context_t ctx, const char* num)
{
  sl3_expr_t res = (sl3_expr_t) malloc (sizeof (struct sl3_expr_t));
  res->discr = SL3_F_NUMBER;
  res->p.ivalue = atoi (num);
  return res;
}

sl3_expr_t
sl3_mk_app (sl3_context_t ctx, const char* name,
            sl3_type_t ty, sl3_expr_t args[], size_t size)
{
  /* combined with mk_op below, it can be called for:
   * - function not defined
   * - variable use (ty is NULL, size is 0)
   * - quantified variable declaration (ty is not NULL, size id 0)
   */
  if (size > 0)
    {
      sl3_error_id (1, "sl3_mk_app", name);
      return NULL;
    }
  // size == 0
  // maybe true or false or variable
  if (strcmp (name, "true") == 0)
    return sl3_mk_true (ctx);
  if (strcmp (name, "false") == 0)
    return sl3_mk_false (ctx);
  return sl3_mk_var (ctx, name, ty);
}

sl3_expr_t
sl3_mk_var (sl3_context_t ctx, const char* name, sl3_type_t ty)
{
  sl3_expr_t res = NULL;
  ap_var_t v = NULL;
  // check that the variable is declared in some environment
  if (name[0] == '?')
    {
      // local variable, check in the local context
      v = (ap_var_t) ap_var_operations->copy((ap_var_t) (name+1));
      // problem with strdup (&(name[1]));
      if (!ap_environment_mem_var (ctx->l_env, v))
        /* for some printing problems in shad.ml, not
         * all local variables are prefixed with `?' */
        {
          /* fprintf (stderr, "Sl3_mk_var: local variable `%s' not declared.\n",
                   name);
           */
          ap_var_operations->free (v);
          v = NULL;
        }
    }
  if (!v)
    {
      // check in the global and local context
      v = (ap_var_t) ap_var_operations->copy((ap_var_t) name);
      // problem with strdup (name);
      if ((ctx->l_env == NULL || !ap_environment_mem_var (ctx->l_env, v))
          &&
          (ctx->g_env == NULL || !ap_environment_mem_var (ctx->g_env, v)))
        {
          sl3_error_id (1, "sl3_mk_var", name);
          ap_var_operations->free (v);
          return NULL;
        }
    }
  res = (sl3_expr_t) malloc (sizeof (struct sl3_expr_t));
  res->discr = SL3_F_SYMBOL;
  res->p.name = ap_var_operations->to_string (v);
  ap_var_operations->free (v);
  return res;
}

sl3_expr_t
sl3_mk_true (sl3_context_t ctx)
{
  sl3_expr_t res = (sl3_expr_t) malloc (sizeof (struct sl3_expr_t));
  res->discr = SL3_F_TRUE;
  return res;
}

sl3_expr_t
sl3_mk_false (sl3_context_t ctx)
{
  sl3_expr_t res = (sl3_expr_t) malloc (sizeof (struct sl3_expr_t));
  res->discr = SL3_F_FALSE;
  return res;
}

sl3_expr_t
sl3_mk_and (sl3_context_t ctx,
            sl3_expr_t args[], size_t size)
{
  /* do typechecking and return tree term */
  if (size <= 0)
  /* 0 arguments is false */
    return sl3_mk_false(ctx);
  else if (size == 1)
    return args[0];
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_AND, args, size);
}

sl3_expr_t
sl3_mk_or (sl3_context_t ctx,
           sl3_expr_t args[], size_t size)
{
  /* do typechecking and return tree term */
  if (size <= 0)
  /* 0 arguments is true */
    return sl3_mk_true(ctx);
  else if (size == 1)
    return args[0];
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_OR, args, size);
}

sl3_expr_t
sl3_mk_not (sl3_context_t ctx,
            sl3_expr_t args[], size_t size)
{
  /* do typechecking and return tree term */
  if (size != 1)
    sl3_error_args (1, "sl3_mk_not", size);
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_NOT, args, size);
}

sl3_expr_t
sl3_mk_implies (sl3_context_t ctx,
                sl3_expr_t args[], size_t size)
{
  /* do typechecking and return tree term */
  if (size != 2)
    sl3_error_args (1, "sl3_mk_implies", size);
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_IMPLIES, args, size);
}

sl3_expr_t
sl3_mk_num (sl3_context_t ctx, int n)
{
  sl3_expr_t res = (sl3_expr_t) malloc (sizeof (struct sl3_expr_t));
  res->discr = SL3_F_NUMBER;
  res->p.ivalue = n;
  return res;
}

sl3_expr_t
sl3_mk_sum (sl3_context_t ctx,
            sl3_expr_t args[], size_t size)
{
  /* do typechecking and return tree term */
  if (size != 2)
    sl3_error_args (1, "sl3_mk_plus", size);
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_PLUS, args, size);
}

sl3_expr_t
sl3_mk_sub (sl3_context_t ctx,
            sl3_expr_t args[], size_t size)
{
  /* do typechecking and return tree term */
  if (size != 2)
    sl3_error_args (1, "sl3_mk_minus", size);
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_MINUS, args, size);
}

sl3_expr_t
sl3_mk_mul (sl3_context_t ctx,
            sl3_expr_t args[], size_t size)
{
  /* do typechecking and return tree term */
  if (size != 2)
    sl3_error_args (1, "sl3_mk_mul", size);
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_MINUS, args, size);
}

sl3_expr_t
sl3_mk_eq (sl3_context_t ctx,
           sl3_expr_t arg1, sl3_expr_t arg2)
{
  /* do typechecking and return tree term */
  sl3_expr_t a[2] = {arg1, arg2};
  return sl3_mk_op (SL3_F_EQ, a, 2);
}

sl3_expr_t
sl3_mk_le (sl3_context_t ctx,
           sl3_expr_t arg1, sl3_expr_t arg2)
{
  /* do typechecking and return tree term */
  sl3_expr_t a[2] = {arg1, arg2};
  return sl3_mk_op (SL3_F_LE, a, 2);
}

sl3_expr_t
sl3_mk_lt (sl3_context_t ctx,
           sl3_expr_t arg1, sl3_expr_t arg2)
{
  sl3_expr_t a[2] = {arg1, arg2};
  return sl3_mk_op (SL3_F_LT, a, 2);
}

sl3_expr_t
sl3_mk_ge (sl3_context_t ctx,
           sl3_expr_t arg1, sl3_expr_t arg2)
{
  sl3_expr_t a[2] = {arg1, arg2};
  return sl3_mk_op (SL3_F_GE, a, 2);
}

sl3_expr_t
sl3_mk_gt (sl3_context_t ctx,
           sl3_expr_t arg1, sl3_expr_t arg2)
{
  sl3_expr_t a[2] = {arg1, arg2};
  return sl3_mk_op (SL3_F_GT, a, 2);
}

sl3_expr_t
sl3_mk_data (sl3_context_t ctx,
             sl3_expr_t args[], size_t size)
{
  if (size != 2)
    sl3_error_args (1, "sl3_mk_data", size);
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_DATA, args, size);
}

sl3_expr_t
sl3_mk_select (sl3_context_t ctx,
               sl3_expr_t args[], size_t size)
{
  if (size != 2)
    sl3_error_args (1, "sl3_mk_select", size);
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_SELECT, args, size);
}

sl3_expr_t
sl3_mk_len (sl3_context_t ctx,
            sl3_expr_t args[], size_t size)
{
  if (size != 1)
    sl3_error_args (1, "sl3_mk_len", size);
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_LEN, args, size);
}

sl3_expr_t
sl3_mk_fsum (sl3_context_t ctx,
             sl3_expr_t args[], size_t size)
{
  if (size != 1)
    sl3_error_args (1, "sl3_mk_sum", size);
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_SUM, args, size);
}

sl3_expr_t
sl3_mk_fmset (sl3_context_t ctx,
              sl3_expr_t args[], size_t size)
{
  if (size != 1)
    sl3_error_args (1, "sl3_mk_mset", size);
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_MSET, args, size);
}

sl3_expr_t
sl3_mk_nilnode (sl3_context_t ctx)
{
  return sl3_mk_op (SL3_F_NIL, NULL, 0);
}

sl3_expr_t
sl3_mk_ls (sl3_context_t ctx,
           sl3_expr_t args[], size_t size)
{
  if (size != 2)
    sl3_error_args (1, "sl3_mk_ls", size);
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_LS, args, size);
}

sl3_expr_t
sl3_mk_label (sl3_context_t ctx,
              sl3_expr_t args[], size_t size)
{
  if (size != 2)
    sl3_error_args (1, "sl3_mk_label", size);
  /* TODO: typechecking */
  return sl3_mk_op (SL3_F_LABEL, args, size);
}

sl3_expr_t
sl3_mk_gall (sl3_context_t ctx,
             sl3_expr_t args[], size_t size)
{
  if (size != 2)
    sl3_error_args (1, "sl3_mk_gall", size);
  /* TODO: typechecking */
  // update set of guards in the context
  ctx->guards_bs |= 1 << 8;
  return sl3_mk_op (SL3_F_GALL, args, size);
}

sl3_expr_t
sl3_mk_gle2 (sl3_context_t ctx,
             sl3_expr_t args[], size_t size)
{
  if (size != 3)
    sl3_error_args (1, "sl3_mk_gle2", size);
  /* TODO: typechecking */
  // update set of guards in the context
  ctx->guards_bs |= 5 << 8;
  return sl3_mk_op (SL3_F_GLE2, args, size);
}

sl3_expr_t
sl3_mk_gsucc2 (sl3_context_t ctx,
               sl3_expr_t args[], size_t size)
{
  if (size != 3)
    sl3_error_args (1, "sl3_mk_gsucc2", size);
  /* TODO: typechecking */
  // update set of guards in the context
  ctx->guards_bs |= 1 << 12; // 16 << 8
  return sl3_mk_op (SL3_F_GSUCC2, args, size);
}

sl3_expr_t
sl3_mk_gfst (sl3_context_t ctx,
             sl3_expr_t args[], size_t size)
{
  if (size != 2)
    sl3_error_args (1, "sl3_mk_gfst", size);
  /* TODO: typechecking */
  // update set of guards in the context
  ctx->guards_bs |= 1 << 12; // add <_1
  return sl3_mk_op (SL3_F_GFST, args, size);
}

sl3_expr_t
sl3_mk_glst (sl3_context_t ctx,
             sl3_expr_t args[], size_t size)
{
  if (size != 2)
    sl3_error_args (1, "sl3_mk_glst", size);
  /* TODO: typechecking */
  // update set of guards in the context
  ctx->guards_bs |= 1 << 12; // add <_1
  return sl3_mk_op (SL3_F_GLST, args, size);
}

sl3_expr_t
sl3_mk_geq2 (sl3_context_t ctx,
             sl3_expr_t args[], size_t size)
{
  if (size != 4)
    sl3_error_args (1, "sl3_mk_geq2", size);
  /* TODO: typechecking */
  // update set of guards in the context
  ctx->guards_bs |= 1 << 9; // 2 << 8
  return sl3_mk_op (SL3_F_GEQ2, args, size);
}

sl3_expr_t
sl3_mk_gsr2 (sl3_context_t ctx,
             sl3_expr_t args[], size_t size)
{
  if (size != 5)
    sl3_error_args (1, "sl3_mk_gsr2", size);
  /* TODO: typechecking */
  // update set of guards in the context
  ctx->guards_bs |= 1 << 12; // TODO
  return sl3_mk_op (SL3_F_GSR2, args, size);
}

sl3_expr_t
sl3_mk_gsl2 (sl3_context_t ctx,
             sl3_expr_t args[], size_t size)
{
  if (size != 5)
    sl3_error_args (1, "sl3_mk_gsl2", size);
  /* TODO: typechecking */
  // update set of guards in the context
  ctx->guards_bs |= 1 << 12; // TODO
  return sl3_mk_op (SL3_F_GSL2, args, size);
}

