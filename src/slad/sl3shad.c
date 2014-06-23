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


#include "sl3.h"
#include "sl3shad.h"
#include "shad.h"

/* ====================================================================== 
 * Translation to shad
 * ====================================================================== */

/** Entry point of the translation.
 *  Deals differently with SL3 logic or AUFLIA logic.
 */
sh_formula_t*
sl3shad (sl3_context_t ctx, sl3_expr_t term, bool sign)
{
  sh_formula_t* res = NULL; // resulting formula

  ap_var_t* qarr = NULL; // used to store quantified nodes for logic AUFLIA
  size_t qsize = 0;

  // set also the guards used
  sh_guards |= ctx->guards_bs;

  // error case
  if (!term) return res;

  // aloc the result
  res = (sh_formula_t*) malloc (sizeof (sh_formula_t));

  // Build the global environement = program pointer and integer vars
  res->env = ap_environment_copy (ctx->g_env);
  if (ctx->logic == SL3_LOG_AUFLIA)
    {
      // For logic AUFLIA, nodes are declared globals
      // - replace each by a ptr var declaration
      // - put nodes in the qarr
      size_t realdim = ctx->g_env->realdim;
      size_t intdim = ctx->g_env->intdim;
      // initialize the array of ptr variables from the one of node variables
      qarr = (ap_var_t*) malloc (realdim * sizeof (ap_var_t));
      qsize = realdim;
      ap_var_t* parr = (ap_var_t*) malloc (realdim * sizeof (ap_var_t));
      for (size_t i = 0; i < realdim; i++)
        {
          char* nname = (char*) ctx->g_env->var_of_dim[intdim + i];
          size_t lname = strlen (nname);
          char* vname = (char*) malloc ((2 + lname) * sizeof (char));
          // generate a name such that ptr vars are befor all node vars
          snprintf (vname, lname + 2, "_%s", nname);
          qarr[i] = (ap_var_t) strdup (nname);
          parr[i] = (ap_var_t) vname;
        }
      // remove nodes from the global environment
      res->env = ap_environment_remove (res->env, &(ctx->g_env->var_of_dim[intdim]),
                                        realdim);
      // add ptr variables
      res->env = ap_environment_add (res->env, NULL, 0, parr, realdim);
#ifndef NDEBUG
      fprintf (stdout, "Sl3shad: formula global env\n");
      ap_environment_fdump (stdout, res->env);
      fflush (stdout);
#endif
      // free ptr variables
      for (size_t i = 0; i < realdim; i++)
        free (parr[i]);
      free (parr);
    }

  if (sign == true)
    {
      size_t n = sl3shad_or_size (term);
      res->form = (sh_shadform_t**) malloc (n * sizeof (sh_shadform_t*));
      res->size = n;
      res = sl3shad_or (ctx, term, qarr, qsize, res, &n);
    }
  else
    {
      size_t n = sl3shad_or_size (term->p.args.arr[0]);
      res->form = (sh_shadform_t**) malloc (n * sizeof (sh_shadform_t*));
      res->size = n;
      res = sl3shad_or (ctx, term->p.args.arr[0], qarr, qsize, res, &n);
    }
  // free qarr
  if (qarr)
    {
      for (size_t i = 0; i < qsize; i++)
        free (qarr[i]);
      free (qarr);
    }

  return res;
}

/** Creates a new <code>shadform</code> and fills <code>f->form</code>
 *  at poisition <code>n</code> with it.
 *  The term may be "exists" or "and" since at least some existential
 *  part shall be given.
 *  Parameters qarr and qsize gives the additional nodes to be inserted
 *  in the node declaration.
 */
sh_formula_t*
sl3shad_exists (sl3_context_t ctx, sl3_expr_t term,
                ap_var_t* qarr, size_t qsize,
                sh_formula_t* f, size_t n)
{
  sh_shadform_t* shf = (sh_shadform_t*) malloc (sizeof (sh_shadform_t));
  if (term && term->discr == SL3_F_EXISTS)
    {
      // build the list of node variables from qarr
      shf->nodes = ap_environment_alloc (NULL, 0,
                                         term->p.quant.qarr,
                                         term->p.quant.qsize);
      shf->nodes = ap_environment_add (shf->nodes, NULL, 0, qarr, qsize);
      // Warning: do not free qarr because used afterwards
      sl3shad_and (ctx, term->p.quant.arg, shf);
    }
  else if (ctx->logic == SL3_LOG_AUFLIA)
    {
      size_t intdim = (f->env) ? f->env->intdim : 0;
      // Nodes are given in qarr and qsize and removed from f->env.
      shf->nodes = ap_environment_alloc (NULL, 0,
                                         qarr, qsize);
      // Add ptr formulas and edge formulas to shf
      shf->length_pform = qsize;
      shf->pform = (sh_labelform_t*) malloc (qsize * sizeof (sh_labelform_t));
      shf->length_eform = qsize;
      shf->eform = (sh_edgeform_t*) malloc (qsize * sizeof (sh_edgeform_t));
      for (size_t i = 0; i < qsize; i++)
        {
          shf->pform[i].node = i;
          shf->pform[i].var = intdim + i;
          shf->eform[i].src = i;
          shf->eform[i].dst = AP_DIM_MAX; // nilNode
#ifndef NDEBUG
          fprintf (stdout, "Sl3shad_exists: (it %zu) label(%zu,%zu) and %zu --> %zu\n",
                   i,
                   shf->pform[i].var, shf->pform[i].node,
                   shf->eform[i].src, shf->eform[i].dst);
#endif
        }
      sl3shad_and (ctx, term, shf);
    }
  else
    {
      // empty list of nodes for this SL3 formula
      shf->nodes = ap_environment_alloc_empty ();
      sl3shad_and (ctx, term, shf);
    }
  // add the built formula to the list
  f->form[f->size - n] = shf;
  return f;
}

/** Disjuncts are collected in <code>f</code> one by one.
 */
sh_formula_t*
sl3shad_or (sl3_context_t ctx, sl3_expr_t term,
            ap_var_t* qarr, size_t qsize,
            sh_formula_t* f, size_t* n)
{
  if (f && n && (f->size >= (*n)))
    {
      switch (term->discr)
        {
        case SL3_F_FALSE:
          sh_formula_free (f);
          return NULL;
        case SL3_F_TRUE:
          free (f->form);
          *n = f->size + 1;
          f->form = NULL;
          f->size = 0;
          break;
        case SL3_F_OR:
          {
            // or is n-ary
            size_t i;
            for (i = 0; i < term->p.args.size; i++)
              sl3shad_or (ctx, term->p.args.arr[i], qarr, qsize, f, n);
            break;
          }
        case SL3_F_EXISTS:
        case SL3_F_AND: /* no nodes declared */
        default:
          sl3shad_exists (ctx, term, qarr, qsize, f, *n);
          *n -= 1;
          break;
        }
    }
  return f;
}

size_t
sl3shad_or_size (sl3_expr_t term)
{
  size_t res = 0;
  if (term)
    {
      switch (term->discr)
        {
        case SL3_F_OR:
          {
            // or is n-ary
            size_t i;
            for (i = 0; i < term->p.args.size; i++)
              res += sl3shad_or_size (term->p.args.arr[i]);
            break;
          }
        case SL3_F_EXISTS:
        case SL3_F_AND:
          res = 1;
          break;
        default:
          break;
        }
    }
  return res;
}

void
sl3shad_and (sl3_context_t ctx, sl3_expr_t term, sh_shadform_t* shf)
{
  size_t nls, nlab, ndt, nfor;
  // compute the number of different kind of formulas 
  nls = nlab = ndt = nfor = 0;
  sl3shad_and_size (term, &nls, &nlab, &ndt, &nfor);
  // except for edge and label formula already computed in logic AUFLIA 
  // (@see sl3shad_exists)
  if (ctx->logic == SL3_LOG_SL3)
    {
      shf->length_eform = nls;
      shf->eform = (nls == 0) ? NULL
              : (sh_edgeform_t*) malloc (nls * sizeof (sh_edgeform_t));
      shf->length_pform = nlab;
      shf->pform = (nlab == 0) ? NULL
              : (sh_labelform_t*) malloc (nlab * sizeof (sh_labelform_t));
    }
  shf->length_dform = ndt;
  shf->dform = (ndt == 0) ? NULL
          : (sh_dataform_t*) malloc (ndt * sizeof (sh_dataform_t));
  shf->length_uform = nfor;
  shf->uform = (nfor == 0) ? NULL
          : (sh_univform_t*) malloc (nfor * sizeof (sh_univform_t));
  sl3shad_and_aux (ctx, term, shf, &nls, &nlab, &ndt, &nfor);
}

void
sl3shad_and_aux (sl3_context_t ctx,
                 sl3_expr_t term, sh_shadform_t* shf,
                 size_t *nls, size_t* nlab, size_t* ndt, size_t *nfor)
{
  if (!term) return;
  switch (term->discr)
    {
    case SL3_F_AND:
      { // and is n-ary
        size_t i;
        for (i = 0; i < term->p.args.size; i++)
          sl3shad_and_aux (ctx, term->p.args.arr[i], shf, nls, nlab, ndt, nfor);
        break;
      }
    case SL3_F_TRUE:
      break;
    case SL3_F_FALSE:
      break;
    case SL3_F_LT:
    case SL3_F_LE:
    case SL3_F_GT:
    case SL3_F_GE:
    case SL3_F_EQ:
      {
        ap_environment_t* nenv = ap_environment_add (shf->nodes,
                                                     ctx->g_env->var_of_dim,
                                                     ctx->g_env->intdim,
                                                     NULL, 0);
        sl3shad_data (nenv, term, shf->dform, shf->length_dform, ndt);
        ap_environment_free (nenv);
        break;
      }
    case SL3_F_FORALL:
      {
        ap_environment_t* nenv = ap_environment_add (shf->nodes,
                                                     ctx->g_env->var_of_dim,
                                                     ctx->g_env->intdim,
                                                     NULL, 0);
        sl3shad_univ (nenv, term, shf->uform, shf->length_uform, nfor);
        ap_environment_free (nenv);
        break;
      }
    case SL3_F_LS:
      sl3shad_edge (ctx, term, shf, nls);
      break;

    case SL3_F_LABEL:
      sl3shad_label (ctx, term, shf, nlab);
      break;

    default: /* error */
      fprintf (stderr, "Sl3shad_and_aux: not an SL3 formula.\n");
    }
}

void
sl3shad_and_size (sl3_expr_t term,
                  size_t *nls, size_t* nlab, size_t* ndt, size_t * nfor)
{
  if (!term) return;
  switch (term->discr)
    {
    case SL3_F_AND:
      {
        // and is n-ary
        size_t i;
        for (i = 0; i < term->p.args.size; i++)
          sl3shad_and_size (term->p.args.arr[i], nls, nlab, ndt, nfor);
        break;
      }
    case SL3_F_TRUE:
      break;
    case SL3_F_FALSE:
      break;
    case SL3_F_LT:
    case SL3_F_LE:
    case SL3_F_GT:
    case SL3_F_GE:
    case SL3_F_EQ:
      (*ndt) += 1;
      break;

    case SL3_F_FORALL:
      (*nfor) += 1;
      break;

    case SL3_F_LS:
      (*nls) += 1;
      break;

    case SL3_F_LABEL:
      (*nlab) += 1;
      break;

    default: /* error */
      fprintf (stderr, "Sl3shad_and_aux: not an SL3 formula.\n");
    }
}

/** Build an edge formula and fills shf->eform[shf->length_eform-*n].
 */
void
sl3shad_edge (sl3_context_t ctx,
              sl3_expr_t term, sh_shadform_t* shf, size_t * n)
{
  if (term && term->discr == SL3_F_LS &&
      term->p.args.arr[0]->discr == SL3_F_SYMBOL &&
      (term->p.args.arr[1]->discr == SL3_F_SYMBOL ||
       term->p.args.arr[1]->discr == SL3_F_NIL))
    {
      char* nn1 = term->p.args.arr[0]->p.name;
      char* nn2 = term->p.args.arr[1]->p.name;
      ap_dim_t n1 = ap_environment_dim_of_var (shf->nodes, (ap_var_t) nn1);
      ap_dim_t n2 = (!nn2) ? AP_DIM_MAX :
              ap_environment_dim_of_var (shf->nodes, (ap_var_t) nn2);
      if (n1 != AP_DIM_MAX)
        {
          sh_edgeform_t ef;
          ef.src = n1;
          ef.dst = n2;
          shf->eform[shf->length_eform - (*n)] = ef;
          (*n) -= 1;
          return;
        }
    }
  // error cases
  fprintf (stderr, "Sl3shad_edge: error in edge formula.\n");
  shf->eform[shf->length_eform - (*n)].src = AP_DIM_MAX;
  shf->eform[shf->length_eform - (*n)].dst = AP_DIM_MAX;
  (*n) -= 1;
  return;
}

void
sl3shad_label (sl3_context_t ctx,
               sl3_expr_t term, sh_shadform_t* shf, size_t * n)
{
  if (term && term->discr == SL3_F_LABEL &&
      term->p.args.arr[0]->discr == SL3_F_SYMBOL &&
      (term->p.args.arr[1]->discr == SL3_F_SYMBOL ||
       term->p.args.arr[1]->discr == SL3_F_NIL))
    {
      char* nv = term->p.args.arr[0]->p.name;
      char* nn = term->p.args.arr[1]->p.name;
      ap_dim_t v = ap_environment_dim_of_var (ctx->g_env, (ap_var_t) nv);
      ap_dim_t n1 = (!nn) ? AP_DIM_MAX :
              ap_environment_dim_of_var (shf->nodes, (ap_var_t) nn);
      if (v != AP_DIM_MAX)
        {
          sh_labelform_t lf;
          lf.var = v;
          lf.node = n1;
          shf->pform[shf->length_pform - (*n)] = lf;
          (*n) -= 1;
          return;
        }
    }
  // error cases
  fprintf (stderr, "Sl3shad_label: error in label formula.\n");
  shf->pform[shf->length_pform - (*n)].var = AP_DIM_MAX;
  shf->pform[shf->length_pform - (*n)].node = AP_DIM_MAX;
  (*n) -= 1;
  return;
}

/** Builds the (shf->length_dform - n)-th existential constraint.
 *  The environment is given by the int part of ctx->g_env
 *  and by shf->nodes.
 */
void
sl3shad_data (ap_environment_t *env,
              sl3_expr_t term, sh_dataform_t* dform, size_t len, size_t * n)
{
  // constant coefficients
  int cst, acst, bcst;
  // left resp. right terms
  sh_linterm_t* a, *b;
  size_t p = len - (*n);
  if (p >= len)
    {
      fprintf (stderr, "Sl3shad_data: internal error.\n");
      return;
    }
  // empty constraint
  dform[p].p = NULL;
  dform[p].cst = 0;
  dform[p].constyp = AP_CONS_EQ;
  (*n) -= 1;
  if (!term ||
      term->discr < SL3_F_LT ||
      term->discr > SL3_F_EQ)
    {
      // error case
      fprintf (stderr, "Sl3shad_data: error in data formula.\n");
      return;
    }
  // < a b   --> b - a -1 >= 0
  // <= a b  --> b - a >= 0
  // > a b   --> a - b -1 >= 0
  // >= a b  --> a - b >= 0
  // = a b   --> a - b = 0
  cst = acst = bcst = 0;
  a = b = NULL;
  if (term->discr == SL3_F_EQ)
    dform[p].constyp = AP_CONS_EQ;
  else
    dform[p].constyp = AP_CONS_SUPEQ;
  if (term->discr == SL3_F_LT || term->discr == SL3_F_GT)
    cst = -1;
  a = sl3shad_linterm (env, term->p.args.arr[0], &acst);
  b = sl3shad_linterm (env, term->p.args.arr[1], &bcst);
  if (term->discr == SL3_F_LT || term->discr == SL3_F_LE)
    {
      a = sh_linterm_merge (b, 1, a, -1); // b - a, free a
      cst += bcst - acst;
    }
  else
    {
      a = sh_linterm_merge (a, 1, b, -1); // a - b, free b
      cst += acst - bcst;
    }
  dform[p].p = a;
  dform[p].cst = cst;
  return;
}

/** Builds the linear term corresponding to term and
 * collects the constant in cst.
 * Dimensions are given wrt to the crt environment env.
 */
sh_linterm_t*
sl3shad_linterm (ap_environment_t* env, sl3_expr_t term,
                 int *cst)
{
  if (term)
    {
      switch (term->discr)
        {
        case SL3_F_NUMBER:
          (*cst) += term->p.ivalue;
          return NULL;

        case SL3_F_PLUS:
          {
            size_t i;
            int clt = 0;
            sh_linterm_t* lt = sl3shad_linterm (env, term->p.args.arr[0], &clt);
            for (i = 1; i < term->p.args.size; i++)
              {
                int crt = 0;
                sh_linterm_t* rt = sl3shad_linterm (env, term->p.args.arr[i], &crt);
                lt = sh_linterm_merge (lt, 1, rt, 1); // in place, free rt
                clt += crt;
              }
            (*cst) = clt;
            return lt;
          }
        case SL3_F_MINUS:
          { // only binary !!
            int clt = 0;
            int crt = 0;
            sh_linterm_t* lt = sl3shad_linterm (env, term->p.args.arr[0], &clt);
            sh_linterm_t* rt = sl3shad_linterm (env, term->p.args.arr[1], &crt);
            lt = sh_linterm_merge (lt, 1, rt, -1); // in place, free rt
            (*cst) = clt - crt;
            return lt;
          }
        case SL3_F_TIMES:
          { // only binary and with a constant
            int clt = 0;
            int crt = 0;
            sh_linterm_t* lt = sl3shad_linterm (env, term->p.args.arr[0], &clt);
            sh_linterm_t* rt = sl3shad_linterm (env, term->p.args.arr[1], &crt);
            if (lt && rt)
              {
                fprintf (stderr, "Sl3shad_linterm: non linear term.\n");
                return NULL;
              }
            else
              {
                int c = (lt) ? crt : clt;
                lt = sh_linterm_merge (lt, crt, rt, clt);
                // multiply the non-zero term by the constant of the zero term
                (*cst) = c * ((lt) ? clt : crt);
                return lt;
              }
          }

        case SL3_F_SYMBOL:
          {
            // can be only an integer data variable since other cases are dealt below
            ap_dim_t d = ap_environment_dim_of_var (env, (ap_var_t) term->p.name);
            sh_linterm_t* lt = sh_linterm_alloc_symbol (d);
            (*cst) += 0;
            return lt;
          }
        case SL3_F_SELECT: /* PENDING, for the moment the same as data */
        case SL3_F_DATA:
          {
            // data is for data(n,0) or data(n,y)
            // some typechecking
            ap_dim_t n, y;
            n = y = AP_DIM_MAX;
            if (term->p.args.size != 2)
              {
                fprintf (stderr, "Sl3shad_linterm: data term with `%zu' arguments.\n",
                         term->p.args.size);
                return NULL;
              }
            if (term->p.args.arr[0]->discr != SL3_F_SYMBOL)
              {
                fprintf (stderr, "Sl3shad_linterm: error on type of the 1st argument of data.\n");
                return NULL;
              }

            if (term->p.args.arr[1]->discr != SL3_F_SYMBOL &&
                term->p.args.arr[1]->discr != SL3_F_NUMBER)
              {
                fprintf (stderr, "Sl3shad_linterm: error on type of the 2nd argument of data.\n");
                return NULL;
              }
            n = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[0]->p.name);
            if (n == AP_DIM_MAX)
              {
                fprintf (stderr, "Sl3shad_linterm: unknown node as 1st argument of data.\n");
                return NULL;
              }
            if (term->p.args.arr[1]->discr == SL3_F_NUMBER &&
                term->p.args.arr[1]->p.ivalue != 0)
              {
                fprintf (stderr, "Sl3shad_linterm: only 0 constant for the 2nd argument of data.\n");
                return NULL;
              }
            if (term->p.args.arr[1]->discr == SL3_F_SYMBOL &&
                (y = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[1]->p.name))
                == AP_DIM_MAX)
              {
                fprintf (stderr, "Sl3shad_linterm: unknown position variable as 2nd argument of data.\n");
                return NULL;
              }
#ifndef NDEBUG1
            fprintf (stderr, "Sl3shad_linterm: build term data(%zu,%zu) in env:\n", n, y);
            sh_env_fprint (stderr, env);
            fflush (stderr);
#endif
            if (term->p.args.arr[1]->discr == SL3_F_NUMBER)
              return sh_linterm_alloc_data (n, AP_DIM_MAX);
            else
              return sh_linterm_alloc_data (n, y);
          }
        case SL3_F_LEN:
          {
            // some typechecking
            ap_dim_t n;
            n = AP_DIM_MAX;
            if (term->p.args.size != 1)
              {
                fprintf (stderr, "Sl3shad_linterm: len term with `%zu' arguments.\n",
                         term->p.args.size);
                return NULL;
              }
            if (term->p.args.arr[0]->discr != SL3_F_SYMBOL)
              {
                fprintf (stderr, "Sl3shad_linterm: error on type of the argument of len.\n");
                return NULL;
              }
            n = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[0]->p.name);
            if (n == AP_DIM_MAX)
              {
                fprintf (stderr, "Sl3shad_linterm: unknown node as argument of len.\n");
                return NULL;
              }
            return sh_linterm_alloc_len (n);
          }
        case SL3_F_SUM:
          {
            // some typechecking
            ap_dim_t n;
            n = AP_DIM_MAX;
            if (term->p.args.size != 1)
              {
                fprintf (stderr, "Sl3shad_linterm: sum term with `%zu' arguments.\n",
                         term->p.args.size);
                return NULL;
              }
            if (term->p.args.arr[0]->discr != SL3_F_SYMBOL)
              {
                fprintf (stderr, "Sl3shad_linterm: error on type of the argument of sum.\n");
                return NULL;
              }
            n = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[0]->p.name);
            if (n == AP_DIM_MAX)
              {
                fprintf (stderr, "Sl3shad_linterm: unknown node as argument of sum.\n");
                return NULL;
              }
            return sh_linterm_alloc_sum (n);
          }
        case SL3_F_MSET:
          {
            // some typechecking
            ap_dim_t n;
            n = AP_DIM_MAX;
            if (term->p.args.size != 1)
              {
                fprintf (stderr, "Sl3shad_linterm: mset term with `%zu' arguments.\n",
                         term->p.args.size);
                return NULL;
              }
            if (term->p.args.arr[0]->discr != SL3_F_SYMBOL)
              {
                fprintf (stderr, "Sl3shad_linterm: error on type of the argument of mset.\n");
                return NULL;
              }
            n = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[0]->p.name);
            if (n == AP_DIM_MAX)
              {
                fprintf (stderr, "Sl3shad_linterm: unknown node as argument of mset.\n");
                return NULL;
              }
            return sh_linterm_alloc_mset (n);
          }

        default:
          fprintf (stderr, "Sl3shad_linterm: data term not supported.\n");
          return NULL;
        }
    }
  return NULL;
}

/** Builds the (shf->length_dform - n)-th universal constraint.
 *  The env is built from nodes and integer global vars.
 */
void
sl3shad_univ (ap_environment_t* env, sl3_expr_t term,
              sh_univform_t* uform, size_t len, size_t * n)
{
  ap_environment_t* nenv = NULL;
  // assert (term && term->discr != SL3_F_FORALL);
  size_t p = len - (*n);
  // put the empty constraint
  uform[p].y = NULL;
  uform[p].length_y = 0;
  uform[p].guard.guardtyp = SH_G_OTHER;
  uform[p].guard.data = NULL;
  uform[p].guard.length_data = 0;
  uform[p].data = NULL;
  uform[p].length_data = 0;
  (*n) -= 1;
  // some well formness test
  if (term->p.quant.arg->discr != SL3_F_IMPLIES &&
      term->p.quant.arg->p.args.size != 2)
    {
      fprintf (stderr, "Sl3shad_univ: not an SL3 formula.\n");
      return;
    }
  // fill y
  uform[p].y = term->p.quant.qarr;
  uform[p].length_y = term->p.quant.qsize;
  // build the environement for this formula
  nenv = ap_environment_add (env, term->p.quant.qarr, term->p.quant.qsize, NULL, 0);
  // fill the guard
  sl3shad_guard (nenv, term->p.quant.arg->p.args.arr[0], &(uform[p].guard));
  // fill the data part
  {
    // count the number of data expressions
    size_t nls, nlab, ndt, nfor;
    nls = nlab = ndt = nfor = 0;
    sl3shad_and_size (term->p.quant.arg->p.args.arr[1], &nls, &nlab, &ndt, &nfor);
    if (nls != 0 || nlab != 0 || nfor != 0 || ndt == 0)
      {
        fprintf (stderr, "Sl3shad_univ: error inside universally quantified SL3 formula.\n");
        return;
      }
    uform[p].length_data = ndt;
    // fill the data part
    uform[p].data = (sh_dataform_t*) malloc (ndt * sizeof (sh_dataform_t));
    sl3shad_univ_data (nenv, term->p.quant.arg->p.args.arr[1], uform[p].data, ndt, &ndt);
  }
  ap_environment_free (nenv);
  return;
}

void
sl3shad_univ_data (ap_environment_t* env, sl3_expr_t term,
                   sh_dataform_t* df, size_t len, size_t* n)
{
  size_t p = len - (*n);
  if (p >= len || !term)
    {
      fprintf (stderr, "Sl3shad_univ_data: internal error.\n");
      return;
    }
  switch (term->discr)
    {
    case SL3_F_AND:
      { // and is n-ary
        size_t i;
        for (i = 0; i < term->p.args.size; i++)
          sl3shad_univ_data (env, term->p.args.arr[i], df, len, n);
        break;
      }
    case SL3_F_TRUE:
      break;
    case SL3_F_FALSE:
      break;
    case SL3_F_LT:
    case SL3_F_LE:
    case SL3_F_GT:
    case SL3_F_GE:
    case SL3_F_EQ:
      sl3shad_data (env, term, df, len, n);
      break;
    default:
      fprintf (stderr, "Sl3shad_univ_data: not an SL3 formula.\n");
      break;
    }
  return;
}

void
sl3shad_guard (ap_environment_t* env, sl3_expr_t term,
               sh_guardform_t * gf)
{
  // assert(term)
  switch (term->discr)
    {
    case SL3_F_GALL: gf->guardtyp = SH_G_ALL;
      break;
    case SL3_F_GFST: gf->guardtyp = SH_G_FST;
      break;
    case SL3_F_GLST: gf->guardtyp = SH_G_LST;
      break;
    case SL3_F_GLE2: gf->guardtyp = SH_G_LE2;
      break;
    case SL3_F_GSUCC2: gf->guardtyp = SH_G_SUCC2;
      break;
    case SL3_F_GEQ2: gf->guardtyp = SH_G_EQ2;
      break;
    case SL3_F_GSL2: gf->guardtyp = SH_G_SL2;
      break;
    case SL3_F_GSR2: gf->guardtyp = SH_G_SR2;
      break;
    case SL3_F_TRUE:
    case SL3_F_FALSE:
    case SL3_F_AND:
    case SL3_F_LT:
    case SL3_F_LE:
    case SL3_F_GT:
    case SL3_F_GE:
    case SL3_F_EQ: /* other guard, ignore */
      fprintf (stderr, "Sl3shad_guard: unsupported guard.\n");
      return;
    default: /* error */
      fprintf (stderr, "Sl3shad_guard: not an SL3 formula.\n");
      return;
    }

  switch (term->discr)
    {
    case SL3_F_GALL:
    case SL3_F_GFST:
    case SL3_F_GLST: // binary guards (y,n)
      {
        ap_dim_t dy, dn;
        dy = dn = AP_DIM_MAX;
        if (term->p.args.arr[0]->discr != SL3_F_SYMBOL)
          {
            fprintf (stderr, "Sl3shad_guard: not a variable for the 1st argument of GALL|GFST|GLST guard.\n");
            return;
          }
        dy = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[0]->p.name);
        if (dy == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: unknown variable `%s' for the 1st argument of GALL|GFST|GLST guard.\n",
                     term->p.args.arr[0]->p.name);
            return;
          }
        if (term->p.args.arr[1]->discr != SL3_F_SYMBOL)
          {
            fprintf (stderr, "Sl3shad_guard: not a variable for the 2nd argument of GALL|GFST|GLST guard.\n");
            return;
          }
        dn = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[1]->p.name);
        if (dn == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: unknown variable `%s' for the 2nd argument of GALL|GFST|GLST guard.\n",
                     term->p.args.arr[1]->p.name);
            return;
          }
        gf->y = (ap_dim_t*) malloc (sizeof (ap_dim_t));
        gf->n = (ap_dim_t*) malloc (sizeof (ap_dim_t));
        gf->y[0] = dy;
        gf->n[0] = dn;
        gf->length_y = 1;
        break;
      }
    case SL3_F_GLE2:
    case SL3_F_GSUCC2: // ternary guards (y1, y2, n)
      {
        ap_dim_t dy1, dy2, dn;
        dy1 = dy2 = dn = AP_DIM_MAX;
        if (term->p.args.arr[0]->discr != SL3_F_SYMBOL)
          {
            fprintf (stderr, "Sl3shad_guard: not a variable for the 1st argument of GLE2|GSUCC2 guard.\n");
            return;
          }
        dy1 = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[0]->p.name);
        if (dy1 == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: unknown variable `%s' for the 1st argument of GLE2|GSUCC2 guard.\n",
                     term->p.args.arr[0]->p.name);
            return;
          }
        if (term->p.args.arr[1]->discr != SL3_F_SYMBOL)
          {
            fprintf (stderr, "Sl3shad_guard: not a variable for the 2nd argument of GLE2|GSUCC2 guard.\n");
            return;
          }
        dy2 = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[1]->p.name);
        if (dy2 == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: unknown variable  `%s' for the 2nd argument of GLE2|GSUCC2 guard.\n",
                     term->p.args.arr[1]->p.name);
            return;
          }
        if (term->p.args.arr[2]->discr != SL3_F_SYMBOL)
          {
            fprintf (stderr, "Sl3shad_guard: not a variable for the 3rd argument of GLE2|GSUCC2 guard.\n");
            return;
          }
        dn = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[2]->p.name);
        if (dn == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: unknown variable `%s' for the 3rd argument of GLE2|GSUCC2 guard.\n",
                     term->p.args.arr[2]->p.name);
            return;
          }
        gf->y = (ap_dim_t*) malloc (2 * sizeof (ap_dim_t));
        gf->n = (ap_dim_t*) malloc (2 * sizeof (ap_dim_t));
        gf->y[0] = dy1;
        gf->y[1] = dy2;
        gf->n[0] = dn;
        gf->n[1] = dn;
        gf->length_y = 2;
        break;
      }


    case SL3_F_GEQ2: // 4-ary guards
      {
        ap_dim_t dy1, dy2;
        ap_dim_t dn1, dn2;
        dy1 = dy2 = dn1 = dn2 = AP_DIM_MAX;
        if (term->p.args.arr[0]->discr != SL3_F_SYMBOL ||
            (dy1 = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[0]->p.name)) == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: error for the 1st argument of GEQ2 guard.\n");
            return;
          }
        if (term->p.args.arr[1]->discr != SL3_F_SYMBOL ||
            (dn1 = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[1]->p.name)) == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: error for the 2nd argument of GEQ2 guard.\n");
            return;
          }
        if (term->p.args.arr[2]->discr != SL3_F_SYMBOL ||
            (dy2 = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[2]->p.name)) == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: error for the 3rd argument of GEQ2 guard.\n");
            return;
          }
        if (term->p.args.arr[3]->discr != SL3_F_SYMBOL ||
            (dn2 = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[3]->p.name)) == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: error for the 4th argument of GEQ2 guard.\n");
            return;
          }
        gf->y = (ap_dim_t*) malloc (2 * sizeof (ap_dim_t));
        gf->n = (ap_dim_t*) malloc (2 * sizeof (ap_dim_t));
        gf->y[0] = dy1;
        gf->y[1] = dy2;
        gf->n[0] = dn1;
        gf->n[1] = dn2;
        gf->length_y = 2;
        break;
      }

    case SL3_F_GSL2:
    case SL3_F_GSR2: // 5-ary guards
      {
        ap_dim_t dy1, dy2;
        ap_dim_t dn1, dn2, dn3; /* PENDING: support an expresions for the 5th argument*/
        dy1 = dy2 = dn1 = dn2 = dn3 = AP_DIM_MAX;
        if (term->p.args.arr[0]->discr != SL3_F_SYMBOL ||
            (dy1 = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[0]->p.name)) == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: error for the 1st argument of GSL2|GSR2 guard.\n");
            return;
          }
        if (term->p.args.arr[1]->discr != SL3_F_SYMBOL ||
            (dn1 = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[1]->p.name)) == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: error for the 2nd argument of GSL2|GSR2 guard.\n");
            return;
          }
        if (term->p.args.arr[2]->discr != SL3_F_SYMBOL ||
            (dy2 = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[2]->p.name)) == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: error for the 3rd argument of GSL2|GSR2 guard.\n");
            return;
          }
        if (term->p.args.arr[3]->discr != SL3_F_SYMBOL ||
            (dn2 = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[3]->p.name)) == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: error for the 4th argument of GSL2|GSR2 guard.\n");
            return;
          }
        if (term->p.args.arr[4]->discr != SL3_F_LEN ||
            term->p.args.arr[4]->p.args.arr[0]->discr != SL3_F_SYMBOL ||
            (dn3 = ap_environment_dim_of_var (env, (ap_var_t) term->p.args.arr[4]->p.args.arr[0]->p.name)) == AP_DIM_MAX)
          {
            fprintf (stderr, "Sl3shad_guard: error for the 5th argument of GSL2|GSR2 guard.\n");
            return;
          }
        gf->y = (ap_dim_t*) malloc (2 * sizeof (ap_dim_t));
        gf->n = (ap_dim_t*) malloc (2 * sizeof (ap_dim_t));
        gf->y[0] = dy1;
        gf->y[1] = dy2;
        gf->n[0] = dn1;
        gf->n[1] = dn2; /* PENDING: store dn3 */
        gf->length_y = 2;
        break;
      }
    default:
      break;
    }
  return;
}
