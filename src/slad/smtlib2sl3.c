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

/* -*- C -*-
 *
 * SMT-LIB v2 interface to Yices 1
 *
 * Author: Alberto Griggio <griggio@fbk.eu>
 *
 * Copyright (C) 2010 Alberto Griggio
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#include "smtlib2sl3.h"
#include <stdlib.h>
#include <string.h>


static void smtlib2_sl3_parser_set_logic (smtlib2_parser_interface *p,
                                          const char *logic);
static void smtlib2_sl3_parser_declare_sort (smtlib2_parser_interface *p,
                                             const char *sortname,
                                             int arity);
static void smtlib2_sl3_parser_define_sort (smtlib2_parser_interface *p,
                                            const char *sortname,
                                            smtlib2_vector *params,
                                            smtlib2_sort sort);
static void smtlib2_sl3_parser_declare_function (smtlib2_parser_interface *p,
                                                 const char *name,
                                                 smtlib2_sort sort);
static void smtlib2_sl3_parser_declare_variable (
                                                 smtlib2_parser_interface *parser,
                                                 const char *name, smtlib2_sort sort);
static void smtlib2_sl3_parser_define_function (
                                                smtlib2_parser_interface *p,
                                                const char *name,
                                                smtlib2_vector *params,
                                                smtlib2_sort sort,
                                                smtlib2_term term);
static void smtlib2_sl3_parser_assert_formula (smtlib2_parser_interface *p,
                                               smtlib2_term term);
static void smtlib2_sl3_parser_check_sat (smtlib2_parser_interface *p);
static smtlib2_term smtlib2_sl3_parser_push_quantifier_scope (
                                                              smtlib2_parser_interface *p);
static smtlib2_term smtlib2_sl3_parser_pop_quantifier_scope (
                                                             smtlib2_parser_interface *p);
/*
 static smtlib2_term smtlib2_sl3_parser_make_term (
                                                  smtlib2_parser_interface *p,
                                                  const char *symbol, smtlib2_sort sort,
                                                  smtlib2_vector *index,
                                                  smtlib2_vector *args);
static smtlib2_term smtlib2_sl3_parser_make_number_term (
                                                         smtlib2_parser_interface *p,
                                                         const char *numval, int width, int base);
 */
// function above are rewritten to call function below
static smtlib2_term smtlib2_sl3_parser_mk_function (smtlib2_context ctx,
                                                    const char *symbol,
                                                    smtlib2_sort sort,
                                                    smtlib2_vector *index,
                                                    smtlib2_vector *args);
static smtlib2_term smtlib2_sl3_parser_mk_number (smtlib2_context ctx,
                                                  const char *rep,
                                                  unsigned int width,
                                                  unsigned int base);
static smtlib2_term smtlib2_sl3_parser_make_forall_term (
                                                         smtlib2_parser_interface *p,
                                                         smtlib2_term term);
static smtlib2_term smtlib2_sl3_parser_make_exists_term (
                                                         smtlib2_parser_interface *p,
                                                         smtlib2_term term);
static smtlib2_sort smtlib2_sl3_parser_make_sort (
                                                  smtlib2_parser_interface *p,
                                                  const char *sortname,
                                                  smtlib2_vector *index);
static smtlib2_sort smtlib2_sl3_parser_make_function_sort (
                                                           smtlib2_parser_interface *p,
                                                           smtlib2_vector *tps);

#define SMTLIB2_SL3_DECLHANDLER(name)				      \
    static smtlib2_term smtlib2_sl3_parser_mk_ ## name (              \
        smtlib2_context ctx,                                            \
        const char *symbol,                                             \
        smtlib2_sort sort,                                              \
        smtlib2_vector *idx,                                            \
        smtlib2_vector *args)

SMTLIB2_SL3_DECLHANDLER (and);
SMTLIB2_SL3_DECLHANDLER (or);
SMTLIB2_SL3_DECLHANDLER (not);
SMTLIB2_SL3_DECLHANDLER (implies);
SMTLIB2_SL3_DECLHANDLER (leq);
SMTLIB2_SL3_DECLHANDLER (lt);
SMTLIB2_SL3_DECLHANDLER (geq);
SMTLIB2_SL3_DECLHANDLER (gt);
SMTLIB2_SL3_DECLHANDLER (eq);
SMTLIB2_SL3_DECLHANDLER (plus);
SMTLIB2_SL3_DECLHANDLER (minus);
SMTLIB2_SL3_DECLHANDLER (times);
SMTLIB2_SL3_DECLHANDLER (data);
SMTLIB2_SL3_DECLHANDLER (select);
SMTLIB2_SL3_DECLHANDLER (len);
SMTLIB2_SL3_DECLHANDLER (sum);
SMTLIB2_SL3_DECLHANDLER (mset);
SMTLIB2_SL3_DECLHANDLER (nilnode);
SMTLIB2_SL3_DECLHANDLER (ls);
SMTLIB2_SL3_DECLHANDLER (label);
SMTLIB2_SL3_DECLHANDLER (gall);
SMTLIB2_SL3_DECLHANDLER (gle2);
SMTLIB2_SL3_DECLHANDLER (gsucc2);
SMTLIB2_SL3_DECLHANDLER (gfst);
SMTLIB2_SL3_DECLHANDLER (glst);
SMTLIB2_SL3_DECLHANDLER (geq2);
SMTLIB2_SL3_DECLHANDLER (gsl2);
SMTLIB2_SL3_DECLHANDLER (gsr2);

#define SMTLIB2_SL3_SETHANDLER(tp, s, name) \
    smtlib2_term_parser_set_handler(tp, s, smtlib2_sl3_parser_mk_ ## name)

/**
 * Initialization of the parser.
 */
smtlib2_sl3_parser *
smtlib2_sl3_parser_new (void)
{
  smtlib2_sl3_parser *ret =
          (smtlib2_sl3_parser *) malloc (sizeof (smtlib2_sl3_parser));
  smtlib2_parser_interface *pi;
  smtlib2_term_parser *tp;
  smtlib2_abstract_parser *ap;

  ret->ctx_ = sl3_mk_context ();
  smtlib2_abstract_parser_init ((smtlib2_abstract_parser *) ret,
                                (smtlib2_context) ret);
  ret->sorts_ = smtlib2_hashtable_new (smtlib2_hashfun_str,
                                       smtlib2_eqfun_str);
  ret->funs_ = smtlib2_hashtable_new (smtlib2_hashfun_str,
                                      smtlib2_eqfun_str);

  /* initialize the term parser and override virtual methods */
  pi = SMTLIB2_PARSER_INTERFACE (ret);
  pi->set_logic = smtlib2_sl3_parser_set_logic;
  pi->declare_sort = smtlib2_sl3_parser_declare_sort;
  pi->define_sort = smtlib2_sl3_parser_define_sort;
  pi->declare_function = smtlib2_sl3_parser_declare_function;
  pi->declare_variable = smtlib2_sl3_parser_declare_variable;
  pi->assert_formula = smtlib2_sl3_parser_assert_formula;
  pi->check_sat = smtlib2_sl3_parser_check_sat;
  pi->push_quantifier_scope = smtlib2_sl3_parser_push_quantifier_scope;
  pi->pop_quantifier_scope = smtlib2_sl3_parser_pop_quantifier_scope;
  // parsing of basic terms set in termparser
  // pi->make_term = smtlib2_sl3_parser_make_term;
  // pi->make_number_term = smtlib2_sl3_parser_make_number_term;
  pi->make_forall_term = smtlib2_sl3_parser_make_forall_term;
  pi->make_exists_term = smtlib2_sl3_parser_make_exists_term;
  pi->make_sort = smtlib2_sl3_parser_make_sort;
  pi->make_function_sort = smtlib2_sl3_parser_make_function_sort;

  tp = ((smtlib2_abstract_parser *) ret)->termparser_;
  smtlib2_term_parser_set_function_handler (tp,
                                            smtlib2_sl3_parser_mk_function);
  smtlib2_term_parser_set_number_handler (tp,
                                          smtlib2_sl3_parser_mk_number);

  SMTLIB2_SL3_SETHANDLER (tp, "or", or);
  SMTLIB2_SL3_SETHANDLER (tp, "and", and);
  SMTLIB2_SL3_SETHANDLER (tp, "not", not);
  SMTLIB2_SL3_SETHANDLER (tp, "=>", implies);
  SMTLIB2_SL3_SETHANDLER (tp, "<=", leq);
  SMTLIB2_SL3_SETHANDLER (tp, "<", lt);
  SMTLIB2_SL3_SETHANDLER (tp, ">=", geq);
  SMTLIB2_SL3_SETHANDLER (tp, ">", gt);
  SMTLIB2_SL3_SETHANDLER (tp, "=", eq);
  SMTLIB2_SL3_SETHANDLER (tp, "+", plus);
  SMTLIB2_SL3_SETHANDLER (tp, "-", minus);
  SMTLIB2_SL3_SETHANDLER (tp, "*", times);
  SMTLIB2_SL3_SETHANDLER (tp, "data", data);
  SMTLIB2_SL3_SETHANDLER (tp, "select", select);
  SMTLIB2_SL3_SETHANDLER (tp, "len", len);
  SMTLIB2_SL3_SETHANDLER (tp, "sum", sum);
  SMTLIB2_SL3_SETHANDLER (tp, "mset", mset);
  SMTLIB2_SL3_SETHANDLER (tp, "nilNode", nilnode);
  SMTLIB2_SL3_SETHANDLER (tp, "ls", ls);
  SMTLIB2_SL3_SETHANDLER (tp, "label", label);
  SMTLIB2_SL3_SETHANDLER (tp, "Gall", gall);
  SMTLIB2_SL3_SETHANDLER (tp, "Gle2", gle2);
  SMTLIB2_SL3_SETHANDLER (tp, "Gsucc2", gsucc2);
  SMTLIB2_SL3_SETHANDLER (tp, "Gfst", gfst);
  SMTLIB2_SL3_SETHANDLER (tp, "Glst", glst);
  SMTLIB2_SL3_SETHANDLER (tp, "Geq2", geq2);
  SMTLIB2_SL3_SETHANDLER (tp, "Gsl2", gsl2);
  SMTLIB2_SL3_SETHANDLER (tp, "Gsr2", gsr2);

  /* the built-in sorts */
  smtlib2_hashtable_set (ret->sorts_,
                         (intptr_t) smtlib2_strdup ("Bool"), (intptr_t) 0);
  smtlib2_hashtable_set (ret->sorts_,
                         (intptr_t) smtlib2_strdup ("Int"), (intptr_t) 0);

  /* set options in the abstract parser */

  ap = (smtlib2_abstract_parser *) pi;
  ap->print_success_ = false;

  return ret;
}

#define SL3CTX(p) (((smtlib2_sl3_parser *)(p))->ctx_)

void
smtlib2_sl3_parser_delete (smtlib2_sl3_parser *p)
{
  smtlib2_hashtable_delete (p->sorts_, NULL, NULL);
  smtlib2_hashtable_delete (p->funs_, NULL, NULL);
  smtlib2_abstract_parser_deinit (&(p->parent_));
  sl3_del_context (p->ctx_);
  free (p);
}

static void
smtlib2_sl3_parser_set_logic (smtlib2_parser_interface *p,
                              const char *logic)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  /* fix logic only one time */
  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      /* check that the logic is supported, i.e., SL3 or AUFLIA */
      sp->logic_ = sl3_set_logic (sp->ctx_, logic);
      if (sp->logic_ == SL3_LOG_OTHER)
        {
          ap->response_ == SMTLIB2_RESPONSE_ERROR;
          ap->errmsg_ = smtlib2_sprintf ("logic `%s' is not supported", logic);
        }
    }
}

/** Sort declaration is supported only for SL3
 *  for sorts `Node' and `Ptr'.
 */
static void
smtlib2_sl3_parser_declare_sort (smtlib2_parser_interface *p,
                                 const char *sortname, int arity)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      /* check that the sort is supported */
      sl3_type_t tp = sl3_mk_type (sp->ctx_, sortname, arity);
      if (sp->logic_ == SL3_LOG_SL3 &&
          tp != NULL)
        { /* check that the sort is not already declared */
          intptr_t k;
          if (smtlib2_hashtable_find (sp->sorts_, (intptr_t) sortname, &k))
            {
              ap->response_ = SMTLIB2_RESPONSE_ERROR;
              ap->errmsg_ = smtlib2_sprintf ("sort `%s' already declared or defined",
                                             sortname);
            }
          else
            {
              smtlib2_hashtable_set (sp->sorts_,
                                     (intptr_t) smtlib2_strdup (sortname),
                                     (intptr_t) tp);
              ap->response_ = SMTLIB2_RESPONSE_SUCCESS;
            }
          // sl3_del_type (tp);
        }
      else
        {
          ap->response_ = SMTLIB2_RESPONSE_ERROR;
          ap->errmsg_ = smtlib2_sprintf ("sort declaration `%s' not supported (logic %d)",
                                         sortname, sp->logic_);
        }
    }
}

/** Sort definition is supported only in AUFLIA 
 *  for sort Node defined as an (Array Int Int).
 */
static void
smtlib2_sl3_parser_define_sort (smtlib2_parser_interface *p,
                                const char *sortname,
                                smtlib2_vector *params,
                                smtlib2_sort sort)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      sl3_type_t tp = (sl3_type_t) sort;
      if (params != NULL
          || sp->logic_ != SL3_LOG_AUFLIA
          || tp->ty != SL3_TYP_NODE
          || !strcmp (sortname, "Node"))
        {
          ap->response_ = SMTLIB2_RESPONSE_ERROR;
          ap->errmsg_ = smtlib2_sprintf ("define-sort `%s' with parameters unsupported (logic %d)",
                                         sortname, sp->logic_);
        }
      else
        {
          intptr_t k;
          if (smtlib2_hashtable_find (sp->sorts_, (intptr_t) sortname, &k))
            {
              ap->response_ = SMTLIB2_RESPONSE_ERROR;
              ap->errmsg_ = smtlib2_sprintf ("sort `%s' already declared or defined",
                                             sortname);
            }
          else
            {
              smtlib2_hashtable_set (sp->sorts_,
                                     (intptr_t) smtlib2_strdup (sortname),
                                     (intptr_t) 0);
              ap->response_ = SMTLIB2_RESPONSE_SUCCESS;
            }
        }
    }
  else if (sp->logic_ == SL3_LOG_AUFLIA &&
           !strcmp (sortname, "Node"))
    {
      smtlib2_hashtable_set (sp->sorts_,
                             (intptr_t) smtlib2_strdup (sortname),
                             (intptr_t) 0);
      ap->response_ = SMTLIB2_RESPONSE_SUCCESS;
    }
}

/** Function declaration supported for both logics
 *  for different sets of predefined functions:
 *  - SL3: functions and predicates on Ptr and Node;
 *  - AUFLIA: function len and guard predicates.
 */
static void
smtlib2_sl3_parser_declare_function (smtlib2_parser_interface *p,
                                     const char *name,
                                     smtlib2_sort sort)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      sl3_type_t f = sl3_mk_fun_decl (sp->ctx_, name, (sl3_type_t) sort);
      if (f)
        {
          intptr_t k;
          // check that it is not already defined
          if (smtlib2_hashtable_find (sp->funs_, (intptr_t) name, &k))
            {
              ap->response_ = SMTLIB2_RESPONSE_ERROR;
              ap->errmsg_ = smtlib2_sprintf ("function `%s' already declared",
                                             name);
            }
          else
            {
              smtlib2_hashtable_set (sp->funs_, (intptr_t) smtlib2_strdup (name),
                                     (intptr_t) f);
              ap->response_ = SMTLIB2_RESPONSE_SUCCESS;
            }
        }
      else
        {
          ap->response_ = SMTLIB2_RESPONSE_ERROR;
          ap->errmsg_ = smtlib2_sprintf ("error declaring function `%s'",
                                         name);
        }
    }
}

/** Called for quantified variable list.
 */
static void
smtlib2_sl3_parser_declare_variable (smtlib2_parser_interface *p,
                                     const char *name,
                                     smtlib2_sort sort)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      sl3_type_t ty = (sl3_type_t) sort;
      if (ty != NULL)
        {
          // variable declaration
          // check that this name starts with  ?
          if (name[0] != '?')
            {
              ap->response_ = SMTLIB2_RESPONSE_ERROR;
              ap->errmsg_ = smtlib2_sprintf ("local variable `%s' is not a constant.",
                                             name);
            }
          else
            {
              // select type to include in environment
              if (ty->size == 0 && ty->ty == SL3_TYP_INT)
                {
                  sl3_push_var (sp->ctx_, (name+1), SL3_TYP_INT);
                }
              if (ty->size == 0 && ty->ty == SL3_TYP_NODE)
                {
                  sl3_push_var (sp->ctx_, (name+1), SL3_TYP_NODE);
                }
            }
        }
    }
}

/** Not called at present.
 *  TODO: called for predicates corresponding to guards.
 */
static void
smtlib2_sl3_parser_define_function (smtlib2_parser_interface *p,
                                    const char *name,
                                    smtlib2_vector *params,
                                    smtlib2_sort sort,
                                    smtlib2_term term)
{
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  // smtlib2_abstract_parser_define_function (p, name, params, sort, term);

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      ap->response_ = SMTLIB2_RESPONSE_ERROR;
      ap->errmsg_ = smtlib2_strdup ("function definition not supported");
    }
}

/** Called for sort expressions built from (sortname index).
 */
static smtlib2_sort
smtlib2_sl3_parser_make_sort (smtlib2_parser_interface *p,
                              const char *sortname,
                              smtlib2_vector * index)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  smtlib2_sort res = NULL;

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      if (index != NULL)
        {
          size_t size = smtlib2_vector_size (index);
          sl3_type_t* tlst = (sl3_type_t *)&(smtlib2_vector_at (index, 0));
          res = (smtlib2_sort) sl3_mk_type_lst (sp->ctx_, sortname, tlst, size);
          if (!res)
            {
              ap->response_ = SMTLIB2_RESPONSE_ERROR;
              ap->errmsg_ = smtlib2_sprintf ("error creating sort starting with %s",
                                             sortname);
            }
        }
      else
        {
          intptr_t v;
          res = (smtlib2_sort) sl3_mk_type (sp->ctx_, sortname, 0);
          ap->response_ = SMTLIB2_RESPONSE_SUCCESS;
          if (!smtlib2_hashtable_find (sp->sorts_, (intptr_t) sortname, &v))
            {
              ap->response_ = SMTLIB2_RESPONSE_ERROR;
              ap->errmsg_ = smtlib2_sprintf ("unknown sort `%s'", sortname);
            }
        }
    }
  return res;
}

/** Called in declare-fun or define-fun to collect typing information.
 */
static smtlib2_sort
smtlib2_sl3_parser_make_function_sort (smtlib2_parser_interface *p,
                                       smtlib2_vector * tps)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  smtlib2_sort ret = NULL;
  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      sl3_type_t *domain;
      sl3_type_t range;
      sl3_type_t tp;
      size_t domain_size;

      domain_size = smtlib2_vector_size (tps) - 1;
      domain = (sl3_type_t *)&(smtlib2_vector_at (tps, 0));
      range = (sl3_type_t) smtlib2_vector_last (tps);

      tp = sl3_mk_fun_type (sp->ctx_, domain, domain_size, range);
      ret = (smtlib2_sort) tp;
    }
  return ret;
}

static smtlib2_sort
smtlib2_sl3_parser_make_parametric_sort (smtlib2_parser_interface *p,
                                         const char *name,
                                         smtlib2_vector * tps)
{
  smtlib2_sl3_parser *yp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  smtlib2_sort ret = NULL;

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      ap->response_ = SMTLIB2_RESPONSE_ERROR;
      ap->errmsg_ = smtlib2_sprintf ("parametric sort `%s' not supported",
                                     name);
    }
  return ret;
}

static void
smtlib2_sl3_parser_assert_formula (smtlib2_parser_interface *p,
                                   smtlib2_term term)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      // check also that assert build a good formula
      if (sl3_assert (sp->ctx_, (sl3_expr_t) term) == false)
        {
          ap->response_ = SMTLIB2_RESPONSE_ERROR;
          ap->errmsg_ = smtlib2_strdup ("assert not an SL3 formula");
        }
      else
        ap->response_ = SMTLIB2_RESPONSE_SUCCESS;
    }
}

static void
smtlib2_sl3_parser_check_sat (smtlib2_parser_interface * p)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      int s = sl3_check (sp->ctx_);
      // returns 1 if phi1 ==> phi2 valid, 
      // i.e., phi1 /\ not(phi2) unsat
      ap->response_ = SMTLIB2_RESPONSE_STATUS;
      switch (s)
        {
        case 1: ap->status_ = SMTLIB2_STATUS_UNSAT;
          break;
        case 0: ap->status_ = SMTLIB2_STATUS_SAT;
          break;
        default: ap->status_ = SMTLIB2_STATUS_UNKNOWN;
        }
    }
}

static smtlib2_term
smtlib2_sl3_parser_push_quantifier_scope (smtlib2_parser_interface * p)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      if (!sl3_push_quant (sp->ctx_))
        {
          ap->response_ = SMTLIB2_RESPONSE_ERROR;
          ap->errmsg_ = smtlib2_strdup ("error in quantifiers");
        }
    }
  return NULL;
}

static smtlib2_term
smtlib2_sl3_parser_pop_quantifier_scope (smtlib2_parser_interface * p)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      if (!sl3_pop_quant (sp->ctx_))
        {
          ap->response_ = SMTLIB2_RESPONSE_ERROR;
          ap->errmsg_ = smtlib2_strdup ("error in quantifiers");
        }
    }
  return NULL;
}

static smtlib2_term
smtlib2_sl3_parser_make_forall_term (smtlib2_parser_interface *p,
                                     smtlib2_term term)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;
  sl3_expr_t res = NULL;

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      res = sl3_mk_forall (sp->ctx_, (sl3_expr_t) term);
      if (!res)
        {
          ap->response_ = SMTLIB2_RESPONSE_ERROR;
          ap->errmsg_ = smtlib2_strdup ("error in quantifiers");
        }
    }
  return res;
}

static smtlib2_term
smtlib2_sl3_parser_make_exists_term (smtlib2_parser_interface *p,
                                     smtlib2_term term)
{
  smtlib2_sl3_parser *sp = (smtlib2_sl3_parser *) p;
  smtlib2_abstract_parser *ap = (smtlib2_abstract_parser *) p;
  sl3_expr_t res = NULL;

  if (ap->response_ != SMTLIB2_RESPONSE_ERROR)
    {
      res = sl3_mk_exists (sp->ctx_, (sl3_expr_t) term);
      if (!res)
        {
          ap->response_ = SMTLIB2_RESPONSE_ERROR;
          ap->errmsg_ = smtlib2_strdup ("error in quantifiers");
        }
    }
  return res;
}

static smtlib2_term
smtlib2_sl3_parser_mk_function (smtlib2_context ctx,
                                const char *symbol,
                                smtlib2_sort sort,
                                smtlib2_vector *index,
                                smtlib2_vector * args)
{
  sl3_context_t sctx = SL3CTX (ctx);
  sl3_expr_t res = NULL;
  if (index)
    { // indexed terms not supported
      return (smtlib2_term) res;
    }
  if (args)
    { // n-ary function, n>0
      res = sl3_mk_app (sctx, symbol, (sl3_type_t) sort,
                        (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                        smtlib2_vector_size (args));
    }
  else
    { // constant, variable or quantified variable
      res = sl3_mk_app (sctx, symbol, (sl3_type_t) sort, NULL, 0);
    }
  // no way to return a message using the context
  return (smtlib2_term) res;
}

static smtlib2_term
smtlib2_sl3_parser_mk_number (smtlib2_context ctx,
                              const char *rep,
                              unsigned int width,
                              unsigned int base)
{
  smtlib2_term ret = NULL;
  if (width == 0 && base == 10)
    {
      ret = (smtlib2_term) sl3_mk_num_from_string (SL3CTX (ctx), rep);
    }
  return (smtlib2_term) ret;
}

SMTLIB2_SL3_DECLHANDLER (and)
{

  return sl3_mk_and (SL3CTX (ctx),
                     (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                     smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (or)
{

  return sl3_mk_or (SL3CTX (ctx),
                    (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                    smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (not)
{

  return sl3_mk_not (SL3CTX (ctx),
                     (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                     smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (implies)
{
  return sl3_mk_implies (SL3CTX (ctx),
                         (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                         smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (eq)
{
  sl3_context_t yctx = SL3CTX (ctx);
  sl3_expr_t ret = sl3_mk_eq (yctx, (sl3_expr_t) smtlib2_vector_at (args, 0),
                              (sl3_expr_t) smtlib2_vector_at (args, 1));
  size_t i;

  for (i = 2; i < smtlib2_vector_size (args); ++i)
    {

      sl3_expr_t prev = (sl3_expr_t) smtlib2_vector_at (args, i - 1);
      sl3_expr_t cur = (sl3_expr_t) smtlib2_vector_at (args, i);
      sl3_expr_t aa[2] = {ret, sl3_mk_eq (yctx, prev, cur)};
      ret = sl3_mk_and (yctx, aa, 2);
    }
  return ret;
}

SMTLIB2_SL3_DECLHANDLER (plus)
{
  return sl3_mk_sum (SL3CTX (ctx),
                     (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                     smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (times)
{
  return sl3_mk_mul (SL3CTX (ctx),
                     (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                     smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (minus)
{
  sl3_expr_t *aa = (sl3_expr_t *) smtlib2_vector_array (args);
  if (smtlib2_vector_size (args) == 1)
    {

      sl3_expr_t a0[2] = {sl3_mk_num (SL3CTX (ctx), -1), aa[0]};
      return sl3_mk_mul (SL3CTX (ctx), a0, 2);
    }
  return sl3_mk_sub (SL3CTX (ctx), aa, smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (leq)
{
  return sl3_mk_le (SL3CTX (ctx),
                    (sl3_expr_t) smtlib2_vector_at (args, 0),
                    (sl3_expr_t) smtlib2_vector_at (args, 1));
}

SMTLIB2_SL3_DECLHANDLER (lt)
{
  return sl3_mk_lt (SL3CTX (ctx),
                    (sl3_expr_t) smtlib2_vector_at (args, 0),
                    (sl3_expr_t) smtlib2_vector_at (args, 1));
}

SMTLIB2_SL3_DECLHANDLER (geq)
{
  return sl3_mk_ge (SL3CTX (ctx),
                    (sl3_expr_t) smtlib2_vector_at (args, 0),
                    (sl3_expr_t) smtlib2_vector_at (args, 1));
}

SMTLIB2_SL3_DECLHANDLER (gt)
{
  return sl3_mk_gt (SL3CTX (ctx),
                    (sl3_expr_t) smtlib2_vector_at (args, 0),
                    (sl3_expr_t) smtlib2_vector_at (args, 1));
}

SMTLIB2_SL3_DECLHANDLER (data)
{
  return sl3_mk_data (SL3CTX (ctx),
                      (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                      smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (select)
{
  return sl3_mk_select (SL3CTX (ctx),
                        (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                        smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (len)
{
  return sl3_mk_len (SL3CTX (ctx),
                     (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                     smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (sum)
{
  return sl3_mk_fsum (SL3CTX (ctx),
                      (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                      smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (mset)
{
  return sl3_mk_fmset (SL3CTX (ctx),
                       (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                       smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (nilnode)
{
  if (!args)
    return sl3_mk_nilnode (SL3CTX (ctx));
  else
    return NULL;
}

SMTLIB2_SL3_DECLHANDLER (ls)
{
  return sl3_mk_ls (SL3CTX (ctx),
                    (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                    smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (label)
{
  return sl3_mk_label (SL3CTX (ctx),
                       (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                       smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (gall)
{
  return sl3_mk_gall (SL3CTX (ctx),
                      (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                      smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (gle2)
{
  return sl3_mk_gle2 (SL3CTX (ctx),
                      (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                      smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (gsucc2)
{
  return sl3_mk_gsucc2 (SL3CTX (ctx),
                        (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                        smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (gfst)
{
  return sl3_mk_gfst (SL3CTX (ctx),
                      (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                      smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (glst)
{
  return sl3_mk_glst (SL3CTX (ctx),
                      (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                      smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (geq2)
{
  return sl3_mk_geq2 (SL3CTX (ctx),
                      (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                      smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (gsl2)
{
  return sl3_mk_gsl2 (SL3CTX (ctx),
                      (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                      smtlib2_vector_size (args));
}

SMTLIB2_SL3_DECLHANDLER (gsr2)
{
  return sl3_mk_gsr2 (SL3CTX (ctx),
                      (sl3_expr_t *)&(smtlib2_vector_at (args, 0)),
                      smtlib2_vector_size (args));
}
