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


#include "shad.h"
#include "shape.h"
#include "shape_fun.h"
#include "smtlib2sl3.h"

/* variable for exchanging between domains */
sh_formula_t* sh_pos = NULL;
sh_formula_t* sh_neg = NULL;
sh_formula_t* sh_crt = NULL;
int sh_guards = 0;

/* ====================================================================== */
/* Constructors/destructors */

/* ====================================================================== */

void
sh_init (void)
{
  if (sh_pos != NULL) sh_formula_free (sh_pos);
  sh_pos = NULL;
  if (sh_neg != NULL) sh_formula_free (sh_neg);
  sh_neg = NULL;
  sh_crt = NULL;
  sh_guards = 0; // no guards
}

void
sh_push2 (sh_formula_t* a, sh_formula_t* b)
{
  size_t i, j;
  if (b->size == 0)
    return;
  i = a->size;
  a->size += b->size;
  a->form = (sh_shadform_t**) realloc (a->form,
                                       a->size * sizeof (sh_shadform_t*));
  for (j = 0; j < b->size; j++)
    {
      a->form[i + j] = b->form[j];
      b->form[j] = NULL;
    }
  b->size = 0;
  free (b->form);
}

void
sh_push (sh_formula_t* f, bool sign)
{
  sh_formula_t** sh = (sign == true) ? &sh_pos : &sh_neg;
  if (!f)
    {
      fprintf (stderr, "Sh_push: NULL formula of sign %d pushed.\n", sign);
      return;
    }

  if (!(*sh))
    {
      (*sh) = f;
    }
  else if (ap_environment_is_eq ((*sh)->env, f->env))
    {
      /* TODO: check f->dw */
      sh_push2 ((*sh), f); /* free f->form */
      ap_environment_free (f->env);
      free (f);
    }
  else
    fprintf (stderr, "Sh_push: formula with different environments.\n");
  return;
}

/** Free array a containing len data formula
 */
void
sh_dataform_free (sh_dataform_t* a, size_t len)
{
  return; /* PENDING */
}

/** Free array a containing len univ formula
 */
void
sh_univform_free (sh_univform_t* a, size_t len)
{
  return; /* PENDING */
}

void
sh_shadform_free (sh_shadform_t* a)
{
  ap_environment_free (a->nodes);
  if (a->eform) free (a->eform);
  if (a->pform) free (a->pform);
  if (a->dform) sh_dataform_free (a->dform, a->length_dform);
  if (a->uform) sh_univform_free (a->uform, a->length_uform);
  return;
}

void
sh_formula_free (sh_formula_t* a)
{
  if (a->env) ap_environment_free (a->env);
  a->env = NULL;
  if (a->size > 0)
    {
      size_t i;
      for (i = 0; i < a->size; i++)
        {
          sh_shadform_free (a->form[i]);
        }
      free (a->form);
      a->form = NULL;
    }
  free (a);
}

sh_linterm_t*
sh_linterm_alloc_1 (sh_funid_t ty, ap_dim_t v)
{
  sh_linterm_t* lt = (sh_linterm_t*) malloc (sizeof (sh_linterm_t));
  lt->t.funid = ty;
  lt->t.p0 = v;
  lt->t.p1 = 0;
  lt->next = NULL;
  lt->coeff = 1;
  return lt;
}

sh_linterm_t*
sh_linterm_alloc_symbol (ap_dim_t v)
{
  return sh_linterm_alloc_1 (SH_F_SYMBOL, v);
}

sh_linterm_t*
sh_linterm_alloc_data (ap_dim_t n, ap_dim_t y)
{
  sh_linterm_t* lt = (sh_linterm_t*) malloc (sizeof (sh_linterm_t));
  lt->t.funid = SH_F_DATA;
  lt->t.p0 = n;
  lt->t.p1 = y; /* MAYBE AP_DIM_MAX for 0 */
  lt->next = NULL;
  lt->coeff = 1;
  return lt;
}

sh_linterm_t*
sh_linterm_alloc_len (ap_dim_t v)
{
  return sh_linterm_alloc_1 (SH_F_LEN, v);
}

sh_linterm_t*
sh_linterm_alloc_sum (ap_dim_t v)
{
  return sh_linterm_alloc_1 (SH_F_SUM, v);
}

sh_linterm_t*
sh_linterm_alloc_mset (ap_dim_t v)
{
  return sh_linterm_alloc_1 (SH_F_MSET, v);
}

void
sh_linterm_free (sh_linterm_t* a)
{
  sh_linterm_t* ai = a;
  while (ai)
    {
      sh_linterm_t* an = ai->next;
      free (ai);
      ai = NULL;
      ai = an;
    }
  return;
}

/**
 * Compares using the lexicographic order on (funid,p0,p1)
 */
int
sh_term_cmp (sh_term_t* t1, sh_term_t* t2)
{
  if (t1->funid < t2->funid)
    return 1;
  else if (t1->funid > t2->funid)
    return -1;
  else // t1->funid == t2->funid
    if (t1->p0 < t2->p0)
    return 1;
  else if (t1->p0 > t2->p0)
    return -1;
  else // t1->p0 == t2->p0
    return ((int) t2->p1 - (int) t1->p1);
}

/** Computes coeffa * a + coeffb * b in a and free common terms
 * Works only on sorted (via sh_term_cmp) lists.
 */
sh_linterm_t*
sh_linterm_merge (sh_linterm_t* a, int coeffa, sh_linterm_t* b, int coeffb)

{
  sh_linterm_t* ai, *bi, *res, *resi;
  ai = a;
  bi = b;
  res = resi = NULL;
  while (ai || bi)
    {
      // compute the next term to be inserted in res
      sh_linterm_t* tmp = NULL;
      int coeff;
      if (ai && bi && sh_term_cmp (&ai->t, &bi->t) == 0)
        {
          coeff = ai->coeff * coeffa + bi->coeff * coeffb;
          // free the cell bi
          tmp = bi;
          bi = bi->next;
          free (tmp);
          tmp = ai;
          ai = ai->next;
          if (coeff != 0)
            {
              // use the cell ai
              tmp->coeff = coeff;
            }
          else // coeff is 0
            { // this term disapear
              // free cell ai
              free (tmp);
              tmp = NULL;
            }
        }
      else if (!bi ||
               (ai && sh_term_cmp (&ai->t, &bi->t) > 0))
        {
          // use cell ai
          coeff = ai->coeff * coeffa;
          tmp = ai;
          ai = ai->next;
          if (coeff != 0)
            {
              tmp->coeff = coeff;
            }
          else
            {
              free (tmp);
              tmp = NULL;
            }
        }
      else if (!ai ||
               (bi && sh_term_cmp (&ai->t, &bi->t) < 0))
        {
          // use cell bi
          coeff = bi->coeff * coeffb;
          tmp = bi;
          bi = bi->next;
          if (coeff != 0)
            {
              tmp->coeff = coeff;
            }
          else
            {
              free (tmp);
              tmp = NULL;
            }
        }
      // insert tmp after resi in res
      if (tmp)
        {
          if (resi)
            {
              resi->next = tmp;
            }
          else
            {
              // first element in the list
              res = tmp;
            }
          resi = tmp;
        }
    }
  return res;
}


/* ====================================================================== */
/* Testing */
/* ====================================================================== */

/**
 * Checks the entailment sh_pos => sh_neg
 * using the sound procedure of shad.
 */
int
sh_check (void)
{
#ifndef NDEBUG1
  fprintf (stdout, "Sh_check: sh_pos = \n");
  sh_fprint (stdout, sh_pos);
  fflush (stdout);
  fprintf (stdout, "Sh_check: sh_neg = \n");
  sh_fprint (stdout, sh_neg);
  fflush (stdout);
#endif
  // Step 1: deal with the trivial cases
  if (!sh_pos)
    {
      // consider sh_pos as true
      fprintf (stdout, "Sh_check: trivial entailment true => phi returns false.\n");
      return 0;
    }
  if (!sh_neg)
    {
      // consider sh_neg as true
      fprintf (stdout, "Sh_check: trivial entailment phi => true returns true.\n");
      return 1;
    }
  // Step 2: collect informations about the dataword domains
  // PENDING: generate files cinv.txt and, for guards, cinv-ucons.txt

  // Step 3: build the abstract values
  // create a manager for shape, parameters are read from cinv files above
  ap_manager_t* shm = shape_manager_alloc ();
  if (!shm) return 0;
  shape_t * top = shape_top (shm, 1, 1);
  for (int i = 0; i < AP_EXC_SIZE; i++)
    {
      shm->option.abort_if_exception[i] = true;
    }
  int max_anon = (sh_guards >= (4 << 8)) ? 1 : 0;
  int segm_anon = (sh_guards >= (2 << 8)) ? 2 : 1;
  shape_approximate (shm, top, max_anon | (segm_anon << 4) | sh_guards);
  shape_free (shm, top);

  // build abstract values from shad formula
  shape_t *shp = shape_of_formula (shm, sh_pos); // assigns shp_crt
  shape_t *shn = shape_of_formula (shm, sh_neg);

  // Step 4: do the inclusion test
  // The inclusion test call saturation if the test does not follow from
  // the shape form
  bool r = shape_is_leq (shm, shp, shn);

  return (r == true) ? 1 : -1;
  ;
}


/* ====================================================================== */
/* Printing */

/* ====================================================================== */

void sh_env_fprint (FILE* stream, ap_environment_t* env);
void sh_shadform_fprint (FILE* stream, ap_environment_t* env, sh_shadform_t* f);
void sh_edgeform_fprint (FILE* stream, ap_environment_t* nodes,
                         sh_edgeform_t* f, size_t size);
void sh_labelform_fprint (FILE* stream, ap_environment_t* env, ap_environment_t* nodes,
                          sh_labelform_t* f, size_t size);
void sh_dataform_fprint (FILE* stream, ap_environment_t* env,
                         sh_dataform_t* f, size_t size);
void sh_univform_fprint (FILE* stream, ap_environment_t* env,
                         sh_univform_t* f, size_t size);

void
sh_fprint (FILE * stream, sh_formula_t * f)
{
  fprintf (stream, "shad formula { \n");
  if (f)
    {
      sh_env_fprint (stream, f->env);
      fprintf (stream, "\tlogic %d;\n", f->dw);
      if (f->form)
        {
          size_t i;
          for (i = 0; i < f->size; i++)
            sh_shadform_fprint (stream, f->env, f->form[i]);
        }

      else
        fprintf (stream, "\tempty shadform;\n");
    }
  fprintf (stream, "\n}\n");

}

void
sh_env_fprint (FILE* stream, ap_environment_t* env)
{
  fprintf (stream, "\tenv: ");
  if (env)
    {
      size_t i, size;
      size = env->intdim + env->realdim;
      fprintf (stream, "(%d,%d) ", env->intdim, env->realdim);
      for (i = 0; i < size; i++)
        {
          fprintf (stream, "%s:%s ",
                   ap_var_operations->to_string (env->var_of_dim[i]),
                   (i < env->intdim) ? "int" : "ptr");
        }
    }

  else
    fprintf (stream, "empty");
  fprintf (stream, ";\n");
  fflush (stream);
}

void
sh_shadform_fprint (FILE* stream, ap_environment_t* env, sh_shadform_t* f)
{
  fprintf (stream, "\tor {\n");
  if (f)
    {

      size_t i;
      ap_environment_t* nnenv;
      fprintf (stream, "\t\tnodes ");
      sh_env_fprint (stream, f->nodes);
      nnenv = ap_environment_alloc (env->var_of_dim, env->intdim,
                                    f->nodes->var_of_dim, f->nodes->realdim);

      fprintf (stream, "\t\tedges ");
      sh_edgeform_fprint (stream, f->nodes, f->eform, f->length_eform);
      fprintf (stream, ";\n");
      fprintf (stream, "\t\tlabels ");
      sh_labelform_fprint (stream, env, f->nodes, f->pform, f->length_pform);
      fprintf (stream, ";\n");
      fprintf (stream, "\t\tdata \t");
      sh_dataform_fprint (stream, nnenv, f->dform, f->length_dform);
      fprintf (stream, ";\n");
      fprintf (stream, "\t\twords ");
      sh_univform_fprint (stream, nnenv, f->uform, f->length_uform);
      fprintf (stream, ";\n");
    }
  fprintf (stream, "\t};\n");
  fflush (stream);
}

void
sh_edgeform_fprint (FILE* stream, ap_environment_t* env,
                    sh_edgeform_t* f, size_t size)
{
  if (f && size > 0)
    {
      size_t i;
      for (i = 0; i < size; i++)
        {
#ifndef NDEBUG
          fprintf (stdout, "Sh_edgeform_fprint: src=%zu --> dst=%zu\n",
                   f[i].src, f[i].dst);
          fflush (stdout);
#endif
          char* n1 = (f[i].src == AP_DIM_MAX) ? "nilNode" :
                  ap_var_operations->to_string (ap_environment_var_of_dim (env, f[i].src));
          char* n2 = (f[i].dst == AP_DIM_MAX) ? "nilNode" :
                  ap_var_operations->to_string (ap_environment_var_of_dim (env, f[i].dst));
          fprintf (stream, "(%s --> %s) ", n1, n2);
        }
    }

  else
    fprintf (stream, "empty");
}

void
sh_labelform_fprint (FILE* stream, ap_environment_t* env, ap_environment_t* nodes,
                     sh_labelform_t* f, size_t size)
{
  if (f && size > 0)
    {
      size_t i;
      for (i = 0; i < size; i++)
        {
          char* nv = (f[i].var == AP_DIM_MAX) ? "??" :
                  ap_var_operations->to_string (ap_environment_var_of_dim (env, f[i].var));
          char* nn = (f[i].node == AP_DIM_MAX) ? "nilNode" :
                  ap_var_operations->to_string (ap_environment_var_of_dim (nodes, f[i].node));
          fprintf (stream, "%s(%s) ", nv, nn);
        }
    }

  else
    fprintf (stream, "empty");
}

void
sh_constyp_fprint (FILE* stream, ap_constyp_t constyp)
{
  switch (constyp)
    {

    case AP_CONS_EQ: fprintf (stream, " == ");
      break;
    case AP_CONS_SUPEQ: fprintf (stream, " >= ");
      break;
    case AP_CONS_SUP: fprintf (stream, " > ");
      break;
    case AP_CONS_EQMOD: fprintf (stream, " =mod ");
      break;
    case AP_CONS_DISEQ:
    default: fprintf (stream, " != ");
      break;
    }
}

void
sh_term_fprint (FILE* stream, ap_environment_t* env,
                sh_funid_t f, ap_dim_t p0, ap_dim_t p1)
{
  switch (f)
    {
    case SH_F_SYMBOL:
      {
        fprintf (stream, "%s",
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, p0)));
        break;
      }
    case SH_F_LEN:
      {
        fprintf (stream, "len(%s)",
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, p0)));
        break;
      }
    case SH_F_DATA:
      {
        fprintf (stream, "data(%s,%s)",
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, p0)),
                 (p1 == AP_DIM_MAX) ? 0 :
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, p1)));
        break;
      }
    case SH_F_SUM:
      {
        fprintf (stream, "sum(%s)",
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, p0)));
        break;
      }
    case SH_F_MSET:
      {
        fprintf (stream, "mset(%s)",
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, p0)));
        break;
      }
    default:
      {

        fprintf (stream, "unknown()");
        break;
      }
    }
}

void
sh_dataform_fprint (FILE* stream, ap_environment_t* env,
                    sh_dataform_t* f, size_t size)
{
  if (f && size > 0)
    {
      size_t i;
      for (i = 0; i < size; i++)
        {
          // print f[i]
          sh_linterm_t* ti;
          ti = f[i].p;
          while (ti)
            {
              fprintf (stream, " + ");
              if (ti->coeff != 1)
                fprintf (stream, "(%d) * ", ti->coeff);
              sh_term_fprint (stream, env, ti->t.funid, ti->t.p0, ti->t.p1);
              ti = ti->next;
            }
          if (f[i].cst != 0)
            fprintf (stream, " + (%d)", f[i].cst);
          sh_constyp_fprint (stream, f[i].constyp);
          fprintf (stream, "0\n\t\t\t");
        }

    }

  else
    fprintf (stream, "empty");
}

void
sh_guardform_fprint (FILE* stream, ap_environment_t* env,
                     sh_guardform_t f)
{
  switch (f.guardtyp)
    {
    case SH_G_ALL:
      {
        fprintf (stream, "Gall(%s,%s)",
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.y[0])),
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.n[0]))
                 );
        break;
      }
    case SH_G_LE2:
    case SH_G_SUCC2:
      {
        fprintf (stream, "G%s2(%s,%s,%s)",
                 (f.guardtyp == SH_G_LE2) ? "le" : "succ",
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.y[0])),
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.y[1])),
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.n[0]))
                 );
#ifndef NDEBUG1
        fprintf (stderr, "Sh_guardform: GLE|SUCC2(%s,%s,%s,%s)\n",
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.y[0])),
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.y[1])),
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.n[0])),
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.n[1]))
                 );
        fflush (stderr);
#endif
        break;
      }
    case SH_G_FST:
    case SH_G_LST:
      {
        fprintf (stream, "G%s(%s,%s)",
                 (f.guardtyp == SH_G_FST) ? "fst" : "lst",
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.y[0])),
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.n[0]))
                 );
        break;
      }
    case SH_G_EQ2:
    case SH_G_SL2:
    case SH_G_SR2:
      {
        fprintf (stream, "G%s2(%s,%s,%s,%s)",
                 (f.guardtyp == SH_G_EQ2) ? "eq" : "slr",
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.y[0])),
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.n[0])),
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.y[1])),
                 ap_var_operations->to_string (ap_environment_var_of_dim (env, f.n[1]))
                 );
        break;
      }
    default:
      {

        fprintf (stream, "G-unknown(%d)", f.length_y);
        break;
      }
    }
}

void
sh_univform_fprint (FILE* stream, ap_environment_t* env,
                    sh_univform_t* f, size_t size)
{
  if (f && size > 0)
    {
      size_t i, j;
      for (i = 0; i < size; i++)
        {
          // print f[i]
          // the new environement contains f[i].y
          ap_environment_t* nenv = ap_environment_add (env,
                                                       f[i].y, f[i].length_y,
                                                       NULL, 0);
          fprintf (stream, "\n\t\t(forall ");
          for (j = 0; j < f[i].length_y; j++)
            fprintf (stream, "%s, ",
                     ap_var_operations->to_string (f[i].y[j]));
          fprintf (stream, " .\n\t\t\t(");
          sh_guardform_fprint (stream, nenv, f[i].guard);
          fprintf (stream, "\n\t\t\t==>\n\t\t\t");
          sh_dataform_fprint (stream, nenv, f[i].data, f[i].length_data);
          fprintf (stream, "\n\t\t\t)\n\t\t)");
          ap_environment_free (nenv);
        }
    }

  else
    fprintf (stream, "empty");
}

/**
 * Reads the formulas @code{sh_pos} and @code{sh_neg}
 * from the file filename.
 * @param filename
 * @return 0 the number of formulas read, 0 if an eror occurs
 */
int
sh_fscan (char* filename)
{
  // open the file
  FILE* f = fopen (filename, "r");
  // create the parser
  smtlib2_sl3_parser *sp = smtlib2_sl3_parser_new ();
  // init formulas
  sh_init ();
  // call the parser
  smtlib2_abstract_parser_parse ((smtlib2_abstract_parser *) sp, f);
  smtlib2_sl3_parser_delete (sp);
  // close the file
  fclose (f);
  // either read sh_pos or (sh_pos and sh_neg)
  if (sh_pos != NULL && sh_neg != NULL) return 2;
  else if (sh_pos != NULL) return 1;
  // sh_neg cannot be read alone
  if (sh_neg) sh_formula_free (sh_neg);
  return 0;
}

