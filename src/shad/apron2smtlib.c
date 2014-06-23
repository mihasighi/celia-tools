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


#include "apron2smtlib.h"

void
ap_linexpr0_fprint_smtlib (FILE* stream, ap_linexpr0_t* a, char** name_of_dim)
{
  size_t i;
  ap_scalar_t* pscalar = 0;
  ap_scalar_t* scalar;
  ap_coeff_t* coeff;
  ap_dim_t dim;
  int sgn;
  int nargs = 0;

  fprintf (stream, "(- ");

  scalar = ap_scalar_alloc ();
  // put only positive terms 
  fprintf (stream, "(+ ");

  ap_linexpr0_ForeachLinterm (a, i, dim, coeff)
  {
    if (!ap_coeff_zero (coeff))
      {
        switch (coeff->discr)
          {
          case AP_COEFF_SCALAR:
            pscalar = coeff->val.scalar;
            sgn = ap_scalar_sgn (pscalar);
            if (sgn >= 0)
              {
                ap_scalar_set (scalar, pscalar);
                if (!ap_scalar_equal_int (scalar, 1))
                  {
                    fprintf (stream, "(* ");
                    ap_scalar_fprint (stream, scalar);
                  }
                if (name_of_dim)
                  fprintf (stream, " %s ", name_of_dim[dim]);
                else
                  fprintf (stream, " x%lu ", (unsigned long) dim);
                if (!ap_scalar_equal_int (scalar, 1))
                  fprintf (stream, ") ");
                nargs++;
              }
            break;
          case AP_COEFF_INTERVAL:
            assert (false);
            break;
          }
      }
  }
  /* Constant */
  if (!ap_coeff_zero (&a->cst))
    {
      switch (a->cst.discr)
        {
        case AP_COEFF_SCALAR:
          pscalar = a->cst.val.scalar;
          sgn = ap_scalar_sgn (pscalar);
          if (sgn >= 0)
            {
              ap_scalar_set (scalar, pscalar);
              ap_scalar_fprint (stream, scalar);
              nargs++;
            }
          break;
        case AP_COEFF_INTERVAL:
          assert (false);
          break;
        }
    }
  if (nargs == 0)
    fprintf (stream, " 0 0 ");
  if (nargs == 1)
    fprintf (stream, " 0 ");
  fprintf (stream, ") ");

  // put only negative terms 
  fprintf (stream, " (+ ");
  nargs = 0;

  ap_linexpr0_ForeachLinterm (a, i, dim, coeff)
  {
    if (!ap_coeff_zero (coeff))
      {
        switch (coeff->discr)
          {
          case AP_COEFF_SCALAR:
            pscalar = coeff->val.scalar;
            sgn = ap_scalar_sgn (pscalar);
            if (sgn < 0)
              {
                ap_scalar_neg (scalar, pscalar);
                if (!ap_scalar_equal_int (scalar, 1))
                  {
                    fprintf (stream, "(* ");
                    ap_scalar_fprint (stream, scalar);
                  }
                if (name_of_dim)
                  fprintf (stream, " %s ", name_of_dim[dim]);
                else
                  fprintf (stream, " x%lu ", (unsigned long) dim);
                if (!ap_scalar_equal_int (scalar, 1))
                  fprintf (stream, ") ");
                nargs++;
              }
            break;
          case AP_COEFF_INTERVAL:
            assert (false);
            break;
          }
      }
  }
  /* Constant */
  if (!ap_coeff_zero (&a->cst))
    {
      switch (a->cst.discr)
        {
        case AP_COEFF_SCALAR:
          pscalar = a->cst.val.scalar;
          sgn = ap_scalar_sgn (pscalar);
          if (sgn < 0)
            {
              ap_scalar_neg (scalar, pscalar);
              ap_scalar_fprint (stream, scalar);
              nargs++;
            }
          break;
        case AP_COEFF_INTERVAL:
          assert (false);
          break;
        }
    }
  if (nargs == 0)
    fprintf (stream, " 0 0 ");
  if (nargs == 1)
    fprintf (stream, " 0 ");
  fprintf (stream, ") ");

  fprintf (stream, ") "); // end of -

  ap_scalar_free (scalar);
}

void
ap_lincons0_fprint_smtlib (FILE* stream, ap_lincons0_t* cons, char** name_of_dim)
{

  fprintf (stream, "\t\t");
  switch (cons->constyp)
    {
    case AP_CONS_EQ:
    case AP_CONS_EQMOD:
      fprintf (stream, "(= 0 ");
      break;
    case AP_CONS_SUPEQ:
      fprintf (stream, "(<= 0 ");
      break;
    case AP_CONS_SUP:
      fprintf (stream, "(< 0 ");
      break;
    default:
      fprintf (stream, "\"ERROR in ap_lincons0_fprint\"");
      break;
    }
  if (cons->constyp == AP_CONS_EQMOD)
    {
      assert (cons->scalar != NULL);
      fprintf (stream, "(mod ");
    }
  ap_linexpr0_fprint_smtlib (stream, cons->linexpr0, name_of_dim);

  if (cons->constyp == AP_CONS_EQMOD)
    {
      assert (cons->scalar != NULL);
      ap_scalar_fprint (stream, cons->scalar);
      fprintf (stream, ") ");
    }
  fprintf (stream, ")\n");
}

void
ap_lincons0_array_fprint_smtlib (FILE* stream, ap_lincons0_array_t* array, char** name_of_dim)
{
  size_t i;

  if (array->size == 0)
    {
      fprintf (stream, "true\n");
    }
  else
    {
      for (i = 0; i < array->size; i++)
        {
          ap_lincons0_fprint_smtlib (stream, &array->p[i], name_of_dim);
        }
    }
}

static const char* smtlib_texpr_op_name[] = {"+", "-", "*", "/", "%", /* binary */
  "-", "cast", "sqrt", /* unary */};

void
ap_texpr0_node_fprint_smtlib (FILE* stream, ap_texpr0_node_t* a, char** name_of_dim)
{
  fprintf (stream, "(%s ", smtlib_texpr_op_name[a->op]);

  /* left argument (if binary) */
  ap_texpr0_fprint_smtlib (stream, a->exprA, name_of_dim);

  /* right argument */
  if (a->exprB)
    {
      fprintf (stream, " ");
      ap_texpr0_fprint_smtlib (stream, a->exprB, name_of_dim);
    }
  fprintf (stream, ") ");
}

void
ap_texpr0_fprint_smtlib (FILE* stream, ap_texpr0_t* a, char** name_of_dim)
{
  if (!a) return;

  ap_scalar_t * pscalar = 0;
  ap_scalar_t * scalar = ap_scalar_alloc ();
  int sgn;
  switch (a->discr)
    {
    case AP_TEXPR_CST:
      // deal with negative constants since not allowed by SMTLIB numerals
      switch (a->val.cst.discr)
        {
        case AP_COEFF_SCALAR:
          pscalar = a->val.cst.val.scalar;
          sgn = ap_scalar_sgn (pscalar);
          if (sgn >= 0)
            {
              ap_scalar_set (scalar, pscalar);
              ap_scalar_fprint (stream, scalar);
            }
          else
            {
              fprintf (stream, "(- 0 ");
              ap_scalar_neg (scalar, pscalar);
              fprintf (stream, ") ");
            }
          break;
        default: assert (false);
        }
      break;
    case AP_TEXPR_DIM:
      if (name_of_dim) fprintf (stream, " %s ", name_of_dim[a->val.dim]);
      else fprintf (stream, " x%lu ", (unsigned long) a->val.dim);
      break;
    case AP_TEXPR_NODE:
      ap_texpr0_node_fprint_smtlib (stream, a->val.node, name_of_dim);
      break;
    default:
      assert (false);
    }

}

void
ap_tcons0_fprint_smtlib (FILE* stream, ap_tcons0_t* cons, char** name_of_dim)
{

  fprintf (stream, "\t\t");
  switch (cons->constyp)
    {
    case AP_CONS_EQ:
    case AP_CONS_EQMOD:
      fprintf (stream, "(= 0 ");
      break;
    case AP_CONS_SUPEQ:
      fprintf (stream, "(<= 0 ");
      break;
    case AP_CONS_SUP:
      fprintf (stream, "(< 0 ");
      break;
    default:
      fprintf (stream, "\"ERROR in ap_tcons0_fprint\"");
      break;
    }

  if (cons->constyp == AP_CONS_EQMOD)
    {
      assert (cons->scalar != NULL);
      fprintf (stream, "(mod ");
    }
  ap_texpr0_fprint_smtlib (stream, cons->texpr0, name_of_dim);

  if (cons->constyp == AP_CONS_EQMOD)
    {
      assert (cons->scalar != NULL);
      ap_scalar_fprint (stream, cons->scalar);
      fprintf (stream, ") ");
    }
  fprintf (stream, ")\n");
}

void
ap_tcons0_array_fprint_smtlib (FILE* stream, ap_tcons0_array_t* array, char** name_of_dim)
{
  size_t i;

  if (array->size == 0)
    {
      fprintf (stream, "true\n");
    }
  else
    {
      for (i = 0; i < array->size; i++)
        {
          ap_tcons0_fprint_smtlib (stream, &array->p[i], name_of_dim);
        }
    }
}

char**
ap_dimension_to_smt (size_t size, size_t datadim, size_t segmdim, size_t offset,
                     char** name_of_dim)
{
  size_t i;
  if (size < (datadim + 2 * segmdim + offset))
    size = datadim + 2 * segmdim + offset;

  char** dname = (char **) malloc (size * sizeof (char *));
  for (i = 0; i < datadim; i++)
    {
      char * n = (name_of_dim) ? name_of_dim[i] : shape_name_of_dim (i);
      size_t lsize = (1 + strlen (n));
      dname[i] = (char *) malloc (lsize * sizeof (char));
      snprintf (dname[i], lsize, "%s", n);
    }
  for (i = 0; i < segmdim; i++)
    {
      char *n = (name_of_dim) ? name_of_dim[datadim + i] : shape_name_of_dim (datadim + i);
      size_t lsize = strlen (n);
      // data of n
      dname[datadim + i] = (char *) malloc ((lsize + 11) * sizeof (char));
      snprintf (dname[datadim + i], (lsize + 11), "(data ?%s 0)", n);
      // length of n, after offset
      dname[datadim + segmdim + offset + i] = (char *) malloc ((lsize + 8) * sizeof (char));
      snprintf (dname[datadim + segmdim + offset + i],
                (lsize + 8), "(len ?%s)", n);
    }
  return dname;
}
