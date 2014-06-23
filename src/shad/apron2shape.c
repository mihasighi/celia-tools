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



#include "apron2shape.h"
#include "shape_macros.h"

/* ====================================================================== */
/* Manipulation of the Htable for constraints.                            */
/* ====================================================================== */

/* Search an entry */
pcons0_t *
shape_pcons_search (shape_internal_t * pr, ap_lincons0_t * lcons,
                    ap_tcons0_t * tcons)
{
  pcons0_t *r, rr;
  unsigned keylen;
  /* search in the htable */
  keylen = offsetof (pcons0_t, hh) - offsetof (pcons0_t, lcons);
  r = NULL;
  rr.lcons = lcons;
  rr.tcons = tcons;
  HASH_FIND (hh, pr->pcons, &rr.lcons, keylen, r);
  return r;
}

/*
 * Add an entry return a pointer to it
 */
pcons0_t *
shape_pcons_add (shape_internal_t * pr, pcons0_t * cons)
{
  pcons0_t *r;
  unsigned keylen;
  /* search in the htable */
  keylen = offsetof (pcons0_t, hh) - offsetof (pcons0_t, lcons);
  r = NULL;
  HASH_FIND (hh, pr->pcons, &cons->lcons, keylen, r);
  if (!r)
    {
      HASH_ADD (hh, pr->pcons, lcons, keylen, cons);
      HASH_FIND (hh, pr->pcons, &cons->lcons, keylen, r);
    }
#ifndef NDEBUG2
  fprintf (stdout, "\n====shape_pcons_add: ");
  shape_pcons_fdump (stdout, cons);
  fprintf (stdout, "\n===returns: ");
  shape_pcons_fdump (stdout, r);
  fprintf (stdout, "\n");
#endif

  return r;
}

/* ====================================================================== */
/* Translation of linear expressions.                                     */
/* ====================================================================== */

/* Extract a (partial) constraint from a linear expression */
void
shape_pcons_of_linexpr (shape_internal_t * pr, ap_linexpr0_t * e,
                        pcons0_t * r, size_t intdim, size_t ptrdim)
{
  size_t i;
  ap_coeff_t *coeff;
  ap_dim_t d;
  arg_assert (r, return;
              );
#ifndef NDEBUG2
  fprintf (stdout, "\n====shape_pcons_of_linexpr: with e=(");
  ap_linexpr0_fprint (stdout, e, NULL);
  fprintf (stdout, ")\n");
#endif
  switch (r->type)
    {
    case EQ_CONS:
    case NE_CONS:
    case ISO_CONS:
    case SEP_CONS:
    case REACH_CONS:
      /* linear expression shall be x - y or x */
      r->info.ptr.x = NULL_DIM;
      r->info.ptr.y = NULL_DIM;
      r->info.ptr.offx = 0;
      r->info.ptr.offy = 0;

      ap_linexpr0_ForeachLinterm (e, i, d, coeff)
      {
        if (!ap_coeff_zero (coeff))
          {
            int coef = ap_scalar_sgn (coeff->val.scalar);
            arg_assert ((coeff->discr == AP_COEFF_SCALAR), return;
                        );
            arg_assert (IS_PTRDIM (d, intdim, ptrdim), return;
                        );
            if (coef == -1)
              r->info.ptr.y = d;
            else if (coef == 1)
              r->info.ptr.x = d;
            else
              arg_assert (0, return;
                          );
          }
      }
      if (r->info.ptr.y == 0 && r->type == REACH_CONS)
        r->type = ACYCLIC_CONS;

      break;

    case ACYCLIC_CONS:
    case CYCLIC_CONS:
      /* linear expression shall be x or x->next - x, search x */
      r->info.ptr.x = NULL_DIM;
      r->info.ptr.y = NULL_DIM;
      r->info.ptr.offx = 0;
      r->info.ptr.offy = 0;

      ap_linexpr0_ForeachLinterm (e, i, d, coeff)
      {
        if (!ap_coeff_zero (coeff))
          {
            arg_assert ((coeff->discr == AP_COEFF_SCALAR), return;
                        );
            arg_assert (IS_PTRDIM (d, intdim, ptrdim), return;
                        );
            /* TODO: be more precise in checking */
            r->info.ptr.x = d;
          }
      }
      break;

    case DATA_CONS:
    case SL3_CONS:
      r->info.data.cons.linexpr0 = ap_linexpr0_copy (e);
      r->info.data.offsets = NULL;
      break;


    case OTHER_CONS:
      {
        /* linear constraint shall be only one ptr variable or NULL */
        bool found = false;
        r->info.ptr.x = NULL_DIM;
        r->info.ptr.offx = 0;

        ap_linexpr0_ForeachLinterm (e, i, d, coeff)
        {
          if (!ap_coeff_zero (coeff))
            {
              int coef = ap_scalar_sgn (coeff->val.scalar);
              arg_assert (!found && (coeff->discr == AP_COEFF_SCALAR), return;
                          );
              found = true;
              if (IS_PTRDIM (d, intdim, ptrdim) && (coef == 1))
                r->info.ptr.x = d;
              else
                arg_assert (0, return;
                            );
              /* TODO: be more precise in checking */
            }
        }
        break;
      }
    }
}

/* ====================================================================== */
/* Translation of linear constraints                                      */

/* ====================================================================== */

pcons0_t *
shape_pcons_of_lincons (shape_internal_t * pr, ap_lincons0_t * c,
                        size_t intdim, size_t ptrdim)
{
  pcons0_t *r;
#ifndef NDEBUG2
  fprintf (stdout, "\n====shape_pcons_of_lincons: with c=(");
  ap_lincons0_fprint (stdout, c, NULL);
  fprintf (stdout, ")\n");
#endif
  /*
   * lincons may be:
   * - only reals for ptr constraints
   * - only integers
   */
  if (!c)
    return NULL;
  else if ((r = shape_pcons_search (pr, c, NULL)) != NULL)
    /* already in the htable */
    return r;
  else
    {
      ap_linexpr0_t *e = c->linexpr0;
      checked_malloc (r, pcons0_t, sizeof (pcons0_t), 1, return NULL;
                      );
      if (ap_linexpr0_is_real (e, intdim))
        {
          /* only constraints */
          /* x == 0.0, x > 0.0, x <> 0.0, x - y == 0.0, x - y > 0 or x <> 0.0 */
          switch (c->constyp)
            {
            case AP_CONS_EQ:
              r->type = EQ_CONS;
              break;
            case AP_CONS_DISEQ:
              r->type = NE_CONS;
              break;

            case AP_CONS_SUP:
              r->type = REACH_CONS;
              break;

            case AP_CONS_EQMOD:
              r->type = ISO_CONS;
              break;
            default: /* ERROR */
              arg_assert (0, return r;
                          );
              break;
            }
        }
      else if (ap_linexpr0_is_integer (e, intdim))
        { /* only constraints on integer variables */
          r->type = DATA_CONS;
          r->info.data.cons.constyp = c->constyp;
          r->info.data.cons.scalar = c->scalar;
        }
      else /* ERROR */
        arg_assert (0, return r;
                    );
      /* fill from expression */
      shape_pcons_of_linexpr (pr, c->linexpr0, r, intdim, ptrdim);
      r->lcons = c;
      r->tcons = NULL;
      //TODO:free r ?
      return shape_pcons_add (pr, r);
    }
}

/* ===================================================================== */
/* Translation of tree expressions.                                      */
/* ===================================================================== */

/* Functions needed for recovering informations from tree expressions */

/* Rewritten because the Apron function is not corect */
bool
shape_coeff_equal_int (ap_coeff_t* coeff, int i)
{
  switch (coeff->discr)
    {
    case AP_COEFF_SCALAR:
      return ap_scalar_equal_int (coeff->val.scalar, i);
    case AP_COEFF_INTERVAL:
      return (ap_scalar_cmp_int (coeff->val.interval->inf, i) <= 0)
              && (ap_scalar_cmp_int (coeff->val.interval->sup, i) >= 0);
    default:
      ;
    }
  return false;
}

/* Texpr e represents a variable of type integer */
bool
shape_is_int_var (shape_internal_t * pr, ap_texpr0_t * e,
                  size_t intdim, size_t ptrdim)
{
  return (e->discr == AP_TEXPR_DIM) && IS_INTDIM (e->val.dim, intdim, ptrdim);
}

/* Texpr e represents a variable of type pointer */
bool
shape_is_ptr_var (shape_internal_t * pr, ap_texpr0_t * e,
                  size_t intdim, size_t ptrdim)
{
  return (e->discr == AP_TEXPR_DIM)
          && IS_PTRDIM (e->val.dim, intdim, ptrdim);
}

/* Texpr e represents a free statement */
bool
shape_is_free_stmt (ap_texpr0_t * e)
{
  return (e->discr == AP_TEXPR_NODE) &&
          (e->val.node->op == AP_TEXPR_SQRT) &&
          IS_FREE_STMT (e->val.node->type, e->val.node->dir);
}

/* Texpr e represents a dereferencing to NULL, i.e., x->f = NULL
 * (encodes in x = cast_r,0(NULL / f)
 * The field f is returned as result (or 0) if not.
 */
int
shape_is_field_null_stmt (ap_texpr0_t *e, size_t intdim, size_t ptrdim)
{

  if ((e->discr == AP_TEXPR_NODE) &&
      (e->val.node->op == AP_TEXPR_CAST) &&
      IS_FIELD_DEREF (e->val.node->type, e->val.node->dir) &&
      (e->val.node->exprA->discr == AP_TEXPR_NODE) &&
      (e->val.node->exprA->val.node->op == AP_TEXPR_DIV) &&
      (e->val.node->exprA->val.node->exprA->discr == AP_TEXPR_CST) &&
      (e->val.node->exprA->val.node->exprB->discr == AP_TEXPR_CST))
    {
      ap_coeff_t e_null = e->val.node->exprA->val.node->exprA->val.cst;
      ap_coeff_t e_fld = e->val.node->exprA->val.node->exprB->val.cst;
#ifndef NDEBUG2
      fprintf (stdout, "+++++ shape_is_field_null_stmt: field=(");
      ap_coeff_fprint (stdout, &e_fld);
      fprintf (stdout, ")\n");
#endif
      if (ap_coeff_zero (&e_null))
        {
          int i;
          for (i = OFFSET_NEXT; i < OFFSET_FLD_MAX; i++)
            if (shape_coeff_equal_int (&e_fld, i))
              {
#ifndef NDEBUG2
                fprintf (stdout, "+++++ shape_is_field_null_stmt: with dim (%zu,%zu) field=(%d)\n",
                         intdim, ptrdim, i);
#endif
                return intdim + OFFSET_DATA + i;
              }
        }
    }
  return OFFSET_NONE;
}

/* Texpr e represents an alloc statement */
bool
shape_is_alloc_stmt (ap_texpr0_t * e)
{

  return (e->discr == AP_TEXPR_NODE) &&
          (e->val.node->op == AP_TEXPR_SQRT) &&
          IS_ALLOC_STMT (e->val.node->type, e->val.node->dir);
}

/* Texpr e represents a data field access or dereferencing. */
bool
shape_is_ptrfield_texpr (shape_internal_t * pr, ap_texpr0_t* e,
                         size_t intdim, size_t ptrdim)
{
  if ((e->discr == AP_TEXPR_NODE) &&
      (e->val.node->op == AP_TEXPR_CAST) &&
      (IS_FIELD_REF (e->val.node->type, e->val.node->dir) ||
       IS_FIELD_DEREF (e->val.node->type, e->val.node->dir)) &&
      (e->val.node->exprA->discr == AP_TEXPR_NODE) &&
      ((e->val.node->exprA->val.node->op == AP_TEXPR_DIV) ||
       (e->val.node->exprA->val.node->op == AP_TEXPR_MUL)) &&
      (e->val.node->exprA->val.node->exprB->discr == AP_TEXPR_CST))
    {
      ap_coeff_t e_fld = e->val.node->exprA->val.node->exprB->val.cst;
      assert (OFFSET_DATA == 0);
      if (!ap_coeff_zero (&e_fld))
        return true;
    }
  return false;
}

/* Texpr e represents a data field access or dereferencing. */
bool
shape_is_datafield_texpr (shape_internal_t * pr, ap_texpr0_t* e,
                          size_t intdim, size_t ptrdim)
{

  if ((e->discr == AP_TEXPR_NODE) &&
      (e->val.node->op == AP_TEXPR_CAST))
    {
      if (IS_ARRAY_REF (e->val.node->type, e->val.node->dir) ||
          IS_ARRAY_DEREF (e->val.node->type, e->val.node->dir))
        return true;
      if ((IS_FIELD_REF (e->val.node->type, e->val.node->dir) ||
           IS_FIELD_DEREF (e->val.node->type, e->val.node->dir)) &&
          (e->val.node->exprA->discr == AP_TEXPR_NODE) &&
          ((e->val.node->exprA->val.node->op == AP_TEXPR_DIV) ||
           (e->val.node->exprA->val.node->op == AP_TEXPR_MUL)) &&
          (e->val.node->exprA->val.node->exprB->discr == AP_TEXPR_CST))
        {
          ap_coeff_t e_fld = e->val.node->exprA->val.node->exprB->val.cst;
          assert (OFFSET_DATA == 0);
#ifndef NDEBUG2
          fprintf (stdout, "+++++shape_is_datafield_texpr: field = (");
          ap_coeff_fprint (stdout, &e_fld);
          fprintf (stdout, ")\n");
#endif
          return ap_coeff_zero (&e_fld);
        }
    }
  return false;
}

/* Texpr e contains a function application (length, sum, multi-set). */
bool
shape_is_fun_texpr (shape_internal_t * pr, ap_texpr0_t* e,
                    size_t intdim, size_t ptrdim)
{
  if (e->discr == AP_TEXPR_NODE &&
      e->val.node->op == AP_TEXPR_CAST)
    return IS_LENGTH_FUN (e->val.node->type, e->val.node->dir) ||
    IS_SUM_FUN (e->val.node->type, e->val.node->dir) ||
    IS_MSET_FUN (e->val.node->type, e->val.node->dir) ||
    IS_UCONS_FUN (e->val.node->type, e->val.node->dir);

  return false;
}

bool
shape_has_field_texpr (shape_internal_t * pr, ap_texpr0_t* e,
                       size_t intdim, size_t ptrdim,
                       bool (*is_field)(shape_internal_t*, ap_texpr0_t*, size_t, size_t))
{
  if (e->discr == AP_TEXPR_NODE)
    {
      if ((*is_field)(pr, e, intdim, ptrdim))
        return true;
      switch (e->val.node->op)
        {
        case AP_TEXPR_ADD:
        case AP_TEXPR_SUB:
        case AP_TEXPR_MUL:
        case AP_TEXPR_DIV:
        case AP_TEXPR_MOD:
          return shape_has_field_texpr (pr, e->val.node->exprA, intdim, ptrdim, is_field) ||
                  shape_has_field_texpr (pr, e->val.node->exprB, intdim, ptrdim, is_field);

        case AP_TEXPR_NEG:
          return shape_has_field_texpr (pr, e->val.node->exprA, intdim, ptrdim, is_field);

        default:
          return false;
        }
    }
  return false;
}

/* Set coefficient in l and offset in offsets for a pointer dimension pdim.
 * If offset can be set (no other offset) return true, otherwise false.
 */
bool
shape_set_pexpr (ap_linexpr0_t * l, int* offsets, ap_dim_t pdim, double coef,
                 int offset, size_t intdim, size_t ptrdim)
{
#ifndef NDEBUG
  fprintf (stdout, "\n++++shape_set_pexpr: for dim=%zu, offset=%d\n",
           pdim, offset);
  fflush (stdout);
#endif
  if (IS_PTRDIM (pdim, intdim, ptrdim)
      && offsets[DIM2PTR (pdim, intdim)] != OFFSET_NONE)
    return false;
  ap_linexpr0_set_coeff_scalar_double (l, pdim, coef);
  offsets[DIM2PTR (pdim, intdim)] = offset; /* TODO: check the offset encoding */

  return true;
}

/* Computes the offset encoded in the expression
 *     cast_t,d(e op offset)
 * where op may be * or / and offset may be field or array.
 */
int
shape_offset_from_cast (ap_texpr_rtype_t t, ap_texpr_rdir_t d,
                        ap_texpr0_t* e, size_t intdim, size_t ptrdim)
{
  if (IS_OFFSET (t, d) &&
      e->discr == AP_TEXPR_NODE &&
      e->val.node->op == AP_TEXPR_CAST)
    {
      ap_texpr0_t* e_op = e->val.node->exprA;
      if ((IS_FIELD_REF (t, d) || IS_FIELD_DEREF (t, d)) &&
          e_op->discr == AP_TEXPR_NODE &&
          (e_op->val.node->op == AP_TEXPR_DIV ||
           e_op->val.node->op == AP_TEXPR_MUL) &&
          e_op->val.node->exprB->discr == AP_TEXPR_CST)
        {
          // extract field no
          size_t i;
#ifndef NDEBUG2
          fprintf (stdout, "===== shape_offset_from_cast: dim=(%zu,%zu), field = (",
                   intdim, ptrdim);
          ap_coeff_fprint (stdout, &e_op->val.node->exprB->val.cst);
          fprintf (stdout, ")\n");
#endif
          assert (OFFSET_DATA == 0);
          for (i = 0; i < OFFSET_FLD_MAX; i++)
            if (shape_coeff_equal_int (&(e_op->val.node->exprB->val.cst), i))
              {
#ifndef NDEBUG2
                fprintf (stdout, "===== shape_offset_from_cast: i = (%zu)\n", i);
#endif
                return intdim + OFFSET_DATA + i;
              }

        }
      else if ((IS_ARRAY_REF (t, d) || IS_ARRAY_DEREF (t, d)) &&
               e_op->discr == AP_TEXPR_NODE &&
               (e_op->val.node->op == AP_TEXPR_DIV ||
                e_op->val.node->op == AP_TEXPR_MUL) &&
               e_op->val.node->exprB->discr == AP_TEXPR_CST)
        {
          // extract constant : only 0 allowed !
          if (ap_coeff_zero (&e_op->val.node->exprB->val.cst))
            return intdim + OFFSET_DATA;
        }
      else if ((IS_ARRAY_REF (t, d) || IS_ARRAY_DEREF (t, d)) &&
               e_op->discr == AP_TEXPR_NODE &&
               (e_op->val.node->op == AP_TEXPR_DIV ||
                e_op->val.node->op == AP_TEXPR_MUL) &&
               e_op->val.node->exprB->discr == AP_TEXPR_DIM)
        {
          // extract dimension of index
          size_t i;
          for (i = 0; i < intdim; i++)
            if (e_op->val.node->exprB->val.dim == i)
              return i;
        }
      else if (IS_ARRAY_SPLIT (t, d) &&
               e_op->discr == AP_TEXPR_NODE &&
               e_op->val.node->op == AP_TEXPR_MUL &&
               e_op->val.node->exprB->discr == AP_TEXPR_DIM)
        {
          // extract dimension of index
          size_t i;
          for (i = 0; i < intdim; i++)
            if (e_op->val.node->exprB->val.dim == i)
              return i;
        }
      return OFFSET_NONE;
    }
  if (IS_LENGTH_FUN (t, d)) return OFFSET_LEN;
  if (IS_SUM_FUN (t, d)) return OFFSET_SUM;
  if (IS_MSET_FUN (t, d)) return OFFSET_MSET;
  if (IS_UCONS_FUN (t, d)) return OFFSET_UCONS;

  return OFFSET_NONE;
}

int*
shape_init_offsets (size_t intdim, size_t ptrdim)
{
  size_t size = ptrdim;
  int* r = (int*) malloc (size * sizeof (int));
  size_t i;
  for (i = 0; i < size; i++)
    r[i] = OFFSET_NONE;

  return r;
}

/*
 * Extract in @code{l} a data constraint from tree expression @code{t};
 * @code{coef} keeps the final coefficient for the expression.
 * We suppose that t is a linear expression with terms 
 *      varint, nb * varint, nb * unary(dataexpr),
 *      nb * cast(varptr * offset)
 * The @code{offsets} array shall be allocate with size @code{ptrdim}.
 *
 * @return the scalar if @code{t} is scalar or NULL otherwise
 */
ap_scalar_t *
shape_dataexpr_of_texpr (shape_internal_t * pr, ap_texpr0_t * t,
                         ap_scalar_t * coef,
                         ap_linexpr0_t * l, int* offsets, size_t intdim, size_t ptrdim)
{

  ap_scalar_t *r = NULL;
  double dcoef = 0.0;

  if (!t)
    return r;
  if (coef)
    ap_double_set_scalar (&dcoef, coef, GMP_RNDN);
  else
    dcoef = 1.0;
  switch (t->discr)
    {
    case AP_TEXPR_CST:
      /* case @code{t} is a scalar */
      if (t->val.cst.discr == AP_COEFF_SCALAR)
        {
          double ncoef;
          ap_double_set_scalar (&ncoef, t->val.cst.val.scalar, GMP_RNDN);
          ncoef = dcoef * ncoef;
          r = ap_scalar_alloc_set_double (ncoef);
        }
      else
        arg_assert (0, return r;
                    );
      break;

    case AP_TEXPR_DIM:
      /* case @code{t} is varint */
      if (IS_INTDIM (t->val.dim, intdim, ptrdim))
        {
          /* compute the final coefficient of this integer var */
          ap_coeff_t *dimcoef = ap_linexpr0_coeffref (l, t->val.dim);
          if (dimcoef && dimcoef->discr == AP_COEFF_SCALAR)
            {
              double oldcoef;
              ap_double_set_scalar (&oldcoef, dimcoef->val.scalar, GMP_RNDN);
              dcoef = (oldcoef == 0.0) ? dcoef : (dcoef * oldcoef);
              ap_linexpr0_set_coeff_scalar_double (l, t->val.dim, dcoef);
            }
          else
            arg_assert (0, return r;
                        );
        }
      else /* ptr vars can be seen only above a cast node texpr */
        arg_assert (0, return r;
                    );
      break;

    case AP_TEXPR_NODE:
      /* all other cases, i.e., linear combination of */
      /* nb, varint, unary(dataexpr),
       * cast(dataexpr / offset), cast(varptr * offset) */
      switch (t->val.node->op)
        {
        case AP_TEXPR_ADD:
        case AP_TEXPR_SUB:
          {
            ap_scalar_t *ncoef = (t->val.node->op == AP_TEXPR_SUB) ?
                    ap_scalar_alloc_set_double (-dcoef) : coef;
            /* recursively apply on subexpressions */
            ap_scalar_t *r1 =
                    shape_dataexpr_of_texpr (pr, t->val.node->exprA,
                                             coef, l, offsets,
                                             intdim, ptrdim);
            ap_scalar_t *r2 =
                    shape_dataexpr_of_texpr (pr, t->val.node->exprB,
                                             ncoef, l, offsets,
                                             intdim, ptrdim);
            /* compute the final value of coefficients for each term */
            double d1, d2;
            if (r1)
              ap_double_set_scalar (&d1, r1, GMP_RNDN);
            else
              d1 = 0.0;
            if (r2)
              ap_double_set_scalar (&d2, r2, GMP_RNDN);
            else
              d2 = 0.0;
            d1 += d2;
            if (r1 || r2)
              {
                ap_coeff_t *oldcoeff = ap_coeff_alloc (AP_COEFF_SCALAR);
                double oldd;
                ap_linexpr0_get_cst (oldcoeff, l);
                ap_double_set_scalar (&oldd, oldcoeff->val.scalar, GMP_RNDN);
                d1 += oldd;
                ap_linexpr0_set_cst_scalar_double (l, d1);
                ap_coeff_free (oldcoeff);
              }
            if (r1)
              ap_scalar_free (r1);
            if (r2)
              ap_scalar_free (r2);
            if (ncoef != coef)
              ap_scalar_free (ncoef);
            break;
          }

        case AP_TEXPR_MUL:
          /*
           * can be only nb * texpr
           */
          {
            if (t->val.node->exprA->discr == AP_TEXPR_CST)
              {
                ap_scalar_t *newcoef =
                        shape_dataexpr_of_texpr (pr, t->val.node->exprA,
                                                 coef, l, offsets,
                                                 intdim, ptrdim);
                shape_dataexpr_of_texpr (pr, t->val.node->exprB,
                                         newcoef, l, offsets,
                                         intdim, ptrdim);
              }
            else
              arg_assert (0, return r;
                          );

            break;
          }

        case AP_TEXPR_NEG:
          {
            ap_scalar_t *newcoef = ap_scalar_alloc_set_double (-dcoef);
            shape_dataexpr_of_texpr (pr, t->val.node->exprA,
                                     newcoef, l, offsets,
                                     intdim, ptrdim);
            ap_scalar_free (newcoef);
            break;
          }

        case AP_TEXPR_CAST:
          /*
           * can be cast(dataexpr) or cast(varptr) or cast(varptr * offset)
           */
          {
            if (IS_OFFSET (t->val.node->type, t->val.node->dir))
              {
                /* the operand shall be a node * */
                arg_assert (t->val.node->exprA->discr == AP_TEXPR_NODE &&
                            t->val.node->exprA->val.node->op == AP_TEXPR_MUL &&
                            t->val.node->exprA->val.node->exprA->discr == AP_TEXPR_DIM &&
                            IS_PTRDIM (t->val.node->exprA->val.node->exprA->val.dim,
                                       intdim, ptrdim), return r;);
                /* set coefficient and offset */
                int ofs = shape_offset_from_cast (
                                                  t->val.node->type, t->val.node->dir,
                                                  t->val.node->exprA,
                                                  intdim, ptrdim);
                ap_dim_t pdim = t->val.node->exprA->val.node->exprA->val.dim;
                arg_assert (shape_set_pexpr (l, offsets, pdim, dcoef, ofs,
                                             intdim, ptrdim),
                            return r;
                            );
              }
            else if (IS_FUN (t->val.node->type, t->val.node->dir))
              {
                /* the operand shall be a varptr */
                arg_assert (t->val.node->exprA->discr == AP_TEXPR_DIM &&
                            IS_PTRDIM (t->val.node->exprA->val.dim,
                                       intdim, ptrdim), return r;);
                /* set ptrvar and offset */
                int ofs = shape_offset_from_cast (
                                                  t->val.node->type, t->val.node->dir,
                                                  t->val.node->exprA,
                                                  intdim, ptrdim);
                ap_dim_t pdim = t->val.node->exprA->val.dim;
                arg_assert (shape_set_pexpr (l, offsets, pdim, dcoef, ofs,
                                             intdim, ptrdim),
                            return r;
                            );
              }
            else
              /* the expression shall be a data expression */
              shape_dataexpr_of_texpr (pr, t->val.node->exprA,
                                       coef, l, offsets,
                                       intdim, ptrdim);
            break;
          }

        default: /* no other operators are valid */
#ifndef NDEBUG1
          fprintf (stdout, "+++++shape_dataexpr_of_texpr: operator %d not known!\n", t->val.node->op);
#endif
          arg_assert (0, return r;
                      );
        }
    }
  return r;
}

/*
 * Extract a partial constraint from a tree expression.
 * @return offset if a dereferencing (cast field or array) is used
 */
int
shape_pcons_of_texpr (shape_internal_t * pr, ap_texpr0_t * e,
                      pcons0_t * r, size_t intdim, size_t ptrdim)
{

  ap_dim_t *dim_array = ap_texpr0_dimlist (e);
  bool hasptr, hasdata, hasdatafield, hasptrfield;
  int deref_offset;
  int i;
#ifndef NDEBUG2
  fprintf (stdout, "\n====shape_pcons_of_texpr: with dim=(%d,%d), e=(",
           intdim, ptrdim);
  ap_texpr0_fprint (stdout, e, NULL);
  fprintf (stdout, "), type=(%d)\n", r->type);
#endif
  /* 1. Collect informations about the kind of expression at hand */
  /* - by looking at the dimensions of variables used */
  for (hasptr = hasdata = false, i = 0; dim_array[i] != AP_DIM_MAX; i++)
    if (IS_INTDIM (dim_array[i], intdim, ptrdim))
      hasdata = true;
    else if (IS_PTRDIM (dim_array[i], intdim, ptrdim))
      hasptr = true;
    else
      arg_assert (0, return false;
                  ); /* bad dimension */
  free (dim_array);
  /* do not assign i, since used below to build linear constraints */

  /* - by traversal of the tree of expression */
  hasdatafield = shape_has_field_texpr (pr, e, intdim, ptrdim,
                                        shape_is_datafield_texpr);
  hasdatafield = hasdatafield ||
          shape_has_field_texpr (pr, e, intdim, ptrdim, shape_is_fun_texpr);
  hasptrfield = shape_has_field_texpr (pr, e, intdim, ptrdim,
                                       shape_is_ptrfield_texpr);

#ifndef NDEBUG2
  fprintf (stdout, "\n====shape_pcons_of_texpr: with info hasdatafield=%zu, hasptrfield=%zu\n",
           hasdatafield, hasptrfield);
#endif

  deref_offset = OFFSET_NONE;

  /* 2. Analyze the expression using the information collected */
  if (hasptr && !hasdata && !hasdatafield)
    {
      /*
       * pure pointer constraints, only x[->f] - y[->f] or
       * x, x->f, or x/f for assign
       */
      arg_assert (r->type != DATA_CONS, return OFFSET_NONE;
                  );

      if (e->discr == AP_TEXPR_DIM)
        {
          r->info.ptr.x = e->val.dim;
          r->info.ptr.offx = OFFSET_NONE;
          r->info.ptr.y = NULL_DIM;
          r->info.ptr.offy = OFFSET_NONE;
        }
      else if (e->discr == AP_TEXPR_NODE)
        {
          /* operand here can be only -, + (with neg), unary or 1*unary */
          /* fill x, y */
          ap_texpr0_t *x, *y;
          x = y = NULL;
          switch (e->val.node->op)
            {
            case AP_TEXPR_ADD:
              arg_assert (r->type == OTHER_CONS, return false;);
              /* one of subnodes shall be NEG */
              /* first operand */
              if (e->val.node->exprA &&
                  e->val.node->exprA->discr == AP_TEXPR_NODE)
                {
                  if (e->val.node->exprA->val.node->op != AP_TEXPR_NEG
                      && x == NULL)
                    x = e->val.node->exprA;
                  else
                    {
                      if (e->val.node->exprA->val.node->op == AP_TEXPR_NEG
                          && y == NULL)
                        y = e->val.node->exprA->val.node->exprA;
                      else
                        arg_assert (0, return OFFSET_NONE;
                                    );
                    }
                }
              else
                {
                  if (e->val.node->exprA &&
                      e->val.node->exprA->discr == AP_TEXPR_DIM &&
                      x == NULL)
                    x = e->val.node->exprA;
                  else
                    arg_assert (0, return OFFSET_NONE;
                                );
                }
              /* second operand */
              if (e->val.node->exprB
                  && e->val.node->exprB->discr == AP_TEXPR_NODE)
                {
                  if (e->val.node->exprB->val.node->op != AP_TEXPR_NEG
                      && x == NULL)
                    x = e->val.node->exprB;
                  else
                    {
                      if (e->val.node->exprB->val.node->op == AP_TEXPR_NEG
                          && y == NULL)
                        y = e->val.node->exprB->val.node->exprA;
                      else
                        arg_assert (0, return OFFSET_NONE;
                                    );
                    }
                }
              else
                {
                  if (e->val.node->exprB &&
                      e->val.node->exprB->discr == AP_TEXPR_DIM &&
                      x == NULL)
                    x = e->val.node->exprB;
                  else /* for ADD, both operands exist */
                    arg_assert (0, return OFFSET_NONE;
                                );
                }
              break;

            case AP_TEXPR_SUB:
              //arg_assert (r->type == OTHER_CONS, return false;);
              /* first operand */
              if (e->val.node->exprA && x == NULL &&
                  (e->val.node->exprA->discr == AP_TEXPR_NODE ||
                   e->val.node->exprA->discr == AP_TEXPR_DIM))
                x = e->val.node->exprA;
              else
                arg_assert (0, return OFFSET_NONE;
                            );

              /* second operand */
              if (e->val.node->exprB && y == NULL)
                y = e->val.node->exprB;
                // TODO: e->val.node->exprB->discr == AP_TEXPR_CST
                // accepted only used for coding INIT
              else
                arg_assert (0, return OFFSET_NONE;
                            );
              break;

            case AP_TEXPR_MUL: /* shall be multiplication by constant 1 */
              {
                ap_coeff_t* coeff1 = ap_coeff_alloc (AP_COEFF_SCALAR);
                ap_coeff_set_scalar_int (coeff1, 1);
                if (e->val.node->exprA->discr == AP_TEXPR_CST &&
                    ap_coeff_equal (coeff1, &e->val.node->exprA->val.cst))
                  x = e->val.node->exprB;
                else
                  if (e->val.node->exprB->discr == AP_TEXPR_CST &&
                      ap_coeff_equal (coeff1, &e->val.node->exprB->val.cst))
                  x = e->val.node->exprA;
                else
                  arg_assert (0, return OFFSET_NONE;);
                ap_coeff_free (coeff1);
                break;
              }

            case AP_TEXPR_CAST: /* uses ptr field */
              arg_assert (r->type == OTHER_CONS, return false;
                          );
              if (IS_FIELD_REF (e->val.node->type, e->val.node->dir) &&
                  e->val.node->exprA->discr == AP_TEXPR_NODE &&
                  e->val.node->exprA->val.node->op == AP_TEXPR_MUL)
                /* ...->next */
                x = e;
              else if (IS_FIELD_DEREF (e->val.node->type, e->val.node->dir) &&
                       e->val.node->exprA->discr == AP_TEXPR_NODE &&
                       e->val.node->exprA->val.node->op == AP_TEXPR_DIV)
                /* ->next = x */
                {
                  x = e->val.node->exprA->val.node->exprA;
                  deref_offset = shape_offset_from_cast (e->val.node->type,
                                                         e->val.node->dir, e, intdim, ptrdim);
                }
              else
                arg_assert (0, return OFFSET_NONE;
                            );
              break;

            default:
              arg_assert (0, return OFFSET_NONE;
                          );
            }
          /* fill r from x */
          r->info.ptr.offx = OFFSET_NONE;
          while (x && x->discr != AP_TEXPR_DIM)
            if (x->discr == AP_TEXPR_NODE &&
                x->val.node->op == AP_TEXPR_CAST &&
                IS_OFFSET (x->val.node->type, x->val.node->dir) &&
                x->val.node->exprA->discr == AP_TEXPR_NODE &&
                x->val.node->exprA->val.node->op == AP_TEXPR_MUL)
              {
                r->info.ptr.offx = shape_offset_from_cast (x->val.node->type,
                                                           x->val.node->dir, x, intdim, ptrdim);
                x = x->val.node->exprA->val.node->exprA;
              }
            else
              arg_assert (0, return OFFSET_NONE;
                          );
          if (x)
            r->info.ptr.x = x->val.dim;

          /* fill r from y */
          r->info.ptr.offy = OFFSET_NONE;
          while (y && y->discr == AP_TEXPR_NODE)
            if (y->val.node->op == AP_TEXPR_CAST &&
                IS_OFFSET (y->val.node->type, y->val.node->dir) &&
                y->val.node->exprA->discr == AP_TEXPR_NODE &&
                y->val.node->exprA->val.node->op == AP_TEXPR_MUL)
              {
                r->info.ptr.offy = shape_offset_from_cast (x->val.node->type,
                                                           x->val.node->dir, x, intdim, ptrdim);
                y = y->val.node->exprA->val.node->exprA;
              }
            else
              arg_assert (0, return OFFSET_NONE;
                          );
          if (y)
            {
              if (y->discr == AP_TEXPR_DIM)
                r->info.ptr.y = y->val.dim;
              else if (y->discr == AP_TEXPR_CST
                       && ap_coeff_zero (&y->val.cst))
                // shall be NULL
                r->info.ptr.y = NULL_DIM;
              else
                arg_assert (0, return OFFSET_NONE;
                            );
            }
          else
            r->info.ptr.y = NULL_DIM;
        }
      else
        arg_assert (0, return OFFSET_NONE;
                    );
    }
  else if (hasptr && hasdata && !hasdatafield && !hasptrfield)
    {
      /* split constraint x * l */
      arg_assert (e->discr == AP_TEXPR_NODE, return OFFSET_NONE;
                  );
      if (e->val.node->op == AP_TEXPR_CAST &&
          IS_ARRAY_SPLIT (e->val.node->type, e->val.node->dir))
	{
          r->info.ptr.x = e->val.node->exprA->val.node->exprA->val.dim;
          r->info.ptr.offx = shape_offset_from_cast (e->val.node->type, 
						     e->val.node->dir,
						     e, intdim, ptrdim);
#ifndef NDEBUG2
	  fprintf (stdout, "\n====shape_pcons_of_texpr: x%d[x%d]\n", 
		   r->info.ptr.x, r->info.ptr.offx);
#endif
          r->info.ptr.y = NULL_DIM;
          r->info.ptr.offy = OFFSET_NONE;
        }
      else
	arg_assert (0, return OFFSET_NONE;
		    );
    }
  else if (hasdatafield && !hasptrfield)
    {
      /*
       * mixed constraints, transform into a linear constraint
       * and an array of offsets for the ptr vars.
       */
      ap_texpr0_t *linexpr;
      arg_assert (e->discr == AP_TEXPR_NODE, return OFFSET_NONE;
                  );
      /* for assign expression, fill dereferencing field */
      if (e->val.node->op == AP_TEXPR_CAST &&
          IS_OFFSET (e->val.node->type, e->val.node->dir))
        {
          /* ->data = texpr or  texpr->data */
          // TODO: assertion no more needed if complex constraints on data
          // arg_assert(r->type == OTHER_CONS, return OFFSET_OTHER;);
          if ((IS_FIELD_DEREF (e->val.node->type, e->val.node->dir) ||
               IS_ARRAY_DEREF (e->val.node->type, e->val.node->dir)) &&
              e->val.node->exprA->discr == AP_TEXPR_NODE &&
              e->val.node->exprA->val.node->op == AP_TEXPR_DIV)
            {
              deref_offset = shape_offset_from_cast (e->val.node->type, e->val.node->dir,
                                                     e, intdim, ptrdim);
#ifndef NDEBUG2
              fprintf (stdout, "\n====shape_pcons_of_texpr: deref_offset=%d\n", deref_offset);
#endif
              linexpr = e->val.node->exprA->val.node->exprA;
            }
          else // field or array reference, normal cast, function application, etc
            linexpr = e;
#ifndef NDEBUG2
          fprintf (stdout, "\n====shape_pcons_of_texpr: offset=%d, op=%d\n", deref_offset, e->val.node->exprA->val.node->op);
#endif
        }
      else
        linexpr = e;
      r->type = DATA_CONS;

      /* fill r->info.data.cons.linexpr0 field from linexpr */
      r->info.data.cons.linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, intdim + ptrdim);
      /* fill r->info.data.offsets also */
      r->info.data.offsets = shape_init_offsets (intdim, ptrdim);
      // call analysis of texpr
      ap_scalar_t *cst =
              shape_dataexpr_of_texpr (pr, linexpr, NULL,
                                       r->info.data.cons.linexpr0, r->info.data.offsets,
                                       intdim, ptrdim);
      if (cst)
        {
          ap_linexpr0_set_cst_scalar (r->info.data.cons.linexpr0, cst);
          ap_scalar_free (cst);
        }
#ifndef NDEBUG
      fprintf (stdout, "\n====shape_pcons_of_texpr: for dim=(%d,%d) returns (", intdim, ptrdim);
      // produces a memory access error because the pcons has no final scalar
      // shape_pcons_fdump(stdout, r);
      fprintf (stdout, ")\n");
      fflush (stdout);
#endif
    }
  else if (hasdata && !hasptr)
    {
      /* shall be a linear constraint on data vars only */
      r->type = DATA_CONS;
      r->info.data.cons.linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, intdim + ptrdim);
      r->info.data.offsets = NULL;
      shape_dataexpr_of_texpr (pr, e, NULL, r->info.data.cons.linexpr0, r->info.data.offsets,
                               intdim, ptrdim);
    }
  else if (!hasptr && !hasdata && !hasdatafield && !hasptrfield)
    {
      /* shall be a constant or other constant expressions... */
      /*MS: commented condition to support SL3 constraints */
      /* if (e->discr == AP_TEXPR_CST) { */
      r->type = (e->discr == AP_TEXPR_CST) ? DATA_CONS : SL3_CONS;
      r->info.data.cons.linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, intdim + ptrdim);
      r->info.data.offsets = NULL;
      ap_scalar_t *cst =
              shape_dataexpr_of_texpr (pr, e, NULL, r->info.data.cons.linexpr0, r->info.data.offsets,
                                       intdim, ptrdim);
      if (cst)
        { /* TODO: test that ths computations is corectly done */
          ap_linexpr0_set_cst_scalar (r->info.data.cons.linexpr0, cst);
          ap_scalar_free (cst);
        }
      /*  } else
          arg_assert(0, return OFFSET_NONE;
          );
       */
    }
  else /* combination of data and next dereferencing is not allowed */
    arg_assert (0, return OFFSET_NONE;
                );

  return deref_offset;
}


/* ===================================================================== */
/* Translation of tree constraints.                                      */

/* ===================================================================== */

pcons0_t *
shape_pcons_of_tcons (shape_internal_t * pr, ap_tcons0_t * c,
                      size_t intdim, size_t ptrdim)
{
  pcons0_t *r;
#ifndef NDEBUG2
  fprintf (stdout, "\n====shape_pcons_of_tcons: with c=(");
  ap_tcons0_fprint (stdout, c, NULL);
  fprintf (stdout, ")\n");
#endif
  /* tcons allowed: */
  /* - trivial constraint 0 = 0 */
  /* - constraints over prg data variables */
  /* - only ptr vars for x[->f] == y[->f] or x[->f] <> y[->f] or x[->f] > y */
  /* - ptr and int vars for constraints using x->data or x[i]  */
  /* Only scalar coefficients are allowed.
   */
  if (!c || !c->texpr0 || !ap_tcons0_is_scalar (c))
    return NULL;
  else if ((r = shape_pcons_search (pr, NULL, c)) != NULL)
    /* already in the htable */
    return r;
  else
    {
      ap_texpr0_t *e = c->texpr0;
      checked_malloc (r, pcons0_t, sizeof (pcons0_t), 1, return NULL;
                      );
      memset (r, 0, sizeof (pcons0_t));
      r->intdim = intdim;
      r->ptrdim = ptrdim;
      switch (c->constyp)
        {
        case AP_CONS_EQ:
          r->type = EQ_CONS;
          break;
        case AP_CONS_DISEQ:
          r->type = NE_CONS;
          break;
        case AP_CONS_EQMOD:
          if (ap_scalar_equal_int (c->scalar, 0))
            r->type = ISO_CONS;
          else
            r->type = SEP_CONS;
          break;

        case AP_CONS_SUPEQ:
        case AP_CONS_SUP:
          r->type = OTHER_CONS;
          break;
        default: /* ERROR */
          arg_assert (0, return r;
                      );
          break;
        }

      /* fill expression, the offset is not important, ignored */
      shape_pcons_of_texpr (pr, e, r, intdim, ptrdim);

      /* put comparison op for data constraints */
      if (r->type == DATA_CONS || r->type == SL3_CONS)
        {
          // special tratment of AP_CONS_SUP constraints to
          // obtain better precision with octagons
          // e > 0 ---> e - 1 >= 0
          if (c->constyp == AP_CONS_SUP)
            {
              ap_coeff_t *oldcoeff = ap_coeff_alloc (AP_COEFF_SCALAR);
              double oldd;
              ap_linexpr0_get_cst (oldcoeff, r->info.data.cons.linexpr0);
              ap_double_set_scalar (&oldd, oldcoeff->val.scalar, GMP_RNDN);
              //fprintf(stdout,"!!***!! scalar old %f \n", oldd);
              oldd -= 1;
              //fprintf(stdout,"!!***!! scalar new %f \n", oldd);
              ap_linexpr0_set_cst_scalar_double (r->info.data.cons.linexpr0,
                                                 oldd);
              ap_coeff_free (oldcoeff);
              r->info.data.cons.constyp = AP_CONS_SUPEQ;
            }
          else
            r->info.data.cons.constyp = c->constyp;

          r->info.data.cons.scalar = NULL;
        }
      else /* pointer constraints */
        {
          arg_assert ((c->constyp == AP_CONS_EQ
                       || c->constyp == AP_CONS_DISEQ
                       || c->constyp == AP_CONS_SUP
                       || c->constyp == AP_CONS_SUPEQ
                       || c->constyp == AP_CONS_EQMOD),
                      return r;
                      );
          if (c->constyp == AP_CONS_SUP || c->constyp == AP_CONS_SUPEQ)
            {
              // may be REACH, ACYCLIC, CYCLIC
              if (r->info.ptr.y == 0)
                r->type = ACYCLIC_CONS;
              else if (r->info.ptr.x == r->info.ptr.y)
                r->type = ACYCLIC_CONS;
              else
                r->type = REACH_CONS;
            }
          else if (c->constyp == AP_CONS_EQMOD)
            {
              if (ap_scalar_equal_int (c->scalar, 0))
                r->type = ISO_CONS;
              else
                r->type = SEP_CONS;
            }
        }
      //TODO:free r ?
      r->lcons = NULL;
      r->tcons = c;

      return shape_pcons_add (pr, r);
    }
}

/* Code for the following functions */
#define PCONS_ARRAY_OF_CONS_ARRAY(pr,array,f,intdim,ptrdim)      \
   pcons0_array_t *arr1, *arr2, *arr3, *arr4;                     \
   size_t i, s1, s2, s3, s4;                   \
   pcons0_t *cons;                       \
   checked_malloc(arr, pcons0_array_t, sizeof(pcons0_array_t), 1, return NULL;);      \
   checked_malloc(arr1, pcons0_array_t, sizeof(pcons0_array_t), 1, return NULL;);      \
   checked_malloc(arr2, pcons0_array_t, sizeof(pcons0_array_t), 1, return NULL;);      \
   checked_malloc(arr3, pcons0_array_t, sizeof(pcons0_array_t), 1, return NULL;);      \
   checked_malloc(arr4, pcons0_array_t, sizeof(pcons0_array_t), 1, return NULL;);      \
   arr->size = arr1->size = arr2->size = arr3->size = arr4->size = array->size;      \
   s = s1 = s2 = s3 = s4 = 0;                \
   checked_malloc(arr->p, pcons0_t*, sizeof(pcons0_t*), array->size, return NULL;);      \
   checked_malloc(arr1->p, pcons0_t*, sizeof(pcons0_t*), array->size, return NULL;);      \
   checked_malloc(arr2->p, pcons0_t*, sizeof(pcons0_t*), array->size, return NULL;);      \
   checked_malloc(arr3->p, pcons0_t*, sizeof(pcons0_t*), array->size, return NULL;);      \
   checked_malloc(arr4->p, pcons0_t*, sizeof(pcons0_t*), array->size, return NULL;);      \
   for (i = 0; i < array->size; i++)    \
     {                                  \
     cons = f(pr, &array->p[i], intdim, ptrdim);       \
     if (cons != NULL && cons->type < OTHER_CONS)       \
     {  \
       if (IS_NULLDIM(cons->info.ptr.x) ||      \
           IS_NULLDIM(cons->info.ptr.y) ||      \
           cons->type == ACYCLIC_CONS ||        \
           cons->type == CYCLIC_CONS)      \
         { arr->p[s] = cons;    \
         s++;   \
         }      \
       else if (cons->type == REACH_CONS)  \
         { arr1->p[s1] = cons;  \
         s1++;  \
         }      \
       else if (cons->type == EQ_CONS || cons->type == ISO_CONS || cons->type == SEP_CONS)  \
         { arr2->p[s2] = cons;  \
         s2++;  \
         }      \
       else if (cons->type == NE_CONS) \
       { arr3->p[s3] = cons;  \
       s3++;  \
       }      \
     else { arr4->p[s4] = cons;  \
     s4++;  \
     }      \
     } \
     } \
     for (i=0; i < s1; i++)     \
       { arr->p[s] = arr1->p[i];        \
       arr1->p[i] = NULL;       \
       s++;     \
       }        \
     free (arr1->p);   free(arr1); \
     for (i=0; i < s2; i++)     \
       { arr->p[s] = arr2->p[i];        \
       arr2->p[i] = NULL;       \
       s++;     \
       }        \
       free(arr2->p); free (arr2);    \
     for (i=0; i < s3; i++)     \
       { arr->p[s] = arr3->p[i];      \
       arr3->p[i] = NULL;       \
       s++;    \
       }        \
       free(arr3->p); free (arr3);    \
     for (i=0; i < s4; i++)     \
       { arr->p[s] = arr4->p[i];        \
       arr4->p[i] = NULL;       \
       s++;     \
       }        \
       free (arr4->p); free (arr4);

pcons0_array_t *
shape_pcons_array_of_lincons_array (shape_internal_t * pr,
                                    ap_lincons0_array_t * array,
                                    size_t intdim, size_t ptrdim)
{
  pcons0_array_t *arr;
  size_t s;
  PCONS_ARRAY_OF_CONS_ARRAY (pr, array, shape_pcons_of_lincons, intdim,
                             ptrdim);

  return arr;
}

pcons0_array_t *
shape_pcons_array_of_tcons_array (shape_internal_t * pr,
                                  ap_tcons0_array_t * array, size_t intdim,
                                  size_t ptrdim)
{
  pcons0_array_t *arr;
  size_t s;
  PCONS_ARRAY_OF_CONS_ARRAY (pr, array, shape_pcons_of_tcons, intdim,
                             ptrdim);

  return arr;
}

/* Fill a linear constraint on data on nodes starting
 * from a linear constraint on ptrvar.
 */
ap_lincons0_t
shape_lincons_of_node (shape_internal_t * pr,
                       ap_lincons0_t * c, int* offsets,
                       size_t * v2n,
                       size_t size, size_t intdim, size_t ptrdim)
{
  ap_linexpr0_t *l =
          shape_linexpr_of_node (pr, c->linexpr0, v2n, size, intdim, ptrdim);
  ap_scalar_t * zero = ap_scalar_alloc ();
  /* TODO: the code below works only for uniforms offsets and not for arrays */
  int kind = intdim + OFFSET_DATA;
  if (offsets != NULL)
    {
      size_t i;
      for (i = 0; i < ptrdim; i++)
        if (offsets[i] != OFFSET_NONE &&
            offsets[i] != kind)
          kind = offsets[i];
    }
  ap_scalar_set_int (zero, kind);

  return ap_lincons0_make (c->constyp, l, zero);
}

/* ===================================================================== */
/* Manipulation of the Htable for assignments.                           */
/* ===================================================================== */

/* Search a value in the htable using keys lhs, lexpr, and texpr */
passign0_t *
shape_passign_search (shape_internal_t * pr, ap_dim_t lhs,
                      size_t intdim, size_t ptrdim,
                      ap_linexpr0_t * lexpr, ap_texpr0_t * texpr)
{
  passign0_t *r, rr;
  unsigned keylen;
  /* search in the htable */
  keylen = offsetof (passign0_t, hh) - offsetof (passign0_t, lhs);
  r = NULL;
  rr.lhs = lhs;
  rr.intdim = intdim;
  rr.ptrdim = ptrdim;
  rr.lexpr = lexpr;
  rr.texpr = texpr;
  HASH_FIND (hh, pr->passigns, &rr.lhs, keylen, r);

  return r;
}

/* Add a value to the htable */
passign0_t *
shape_passign_add (shape_internal_t * pr, passign0_t * a)
{
  passign0_t *r;
  unsigned keylen;
  /* search in the htable */
  keylen = offsetof (passign0_t, hh) - offsetof (passign0_t, lhs);
  r = NULL;
  HASH_FIND (hh, pr->passigns, &a->lhs, keylen, r);
  if (!r)
    {
      HASH_ADD (hh, pr->passigns, lhs, keylen, a);
      HASH_FIND (hh, pr->passigns, &a->lhs, keylen, r);
    }
  return r;
}

/* ===================================================================== */
/* Translation of linear assignments.                                    */

/* ===================================================================== */

passign0_t *
shape_passign_of_linexpr (shape_internal_t * pr,
                          ap_dim_t lhs, ap_linexpr0_t * rhs,
                          size_t intdim, size_t ptrdim)
{
  passign0_t *r;

  arg_assert (rhs && lhs < intdim + ptrdim, return NULL;
              );

#ifndef NDEBUG
  fprintf (stdout, "****shape_passign_of_linexpr: x%zu = ", lhs);
  ap_linexpr0_fprint (stdout, rhs, NULL);
  fprintf (stdout, "\n");
#endif

  if ((r = shape_passign_search (pr, lhs, intdim, ptrdim, rhs, NULL)) != NULL)
    return r;

  checked_malloc (r, passign0_t, sizeof (passign0_t), 1, return NULL;
                  );

  /* Fill the lhs part */
  r->x = lhs;
  r->intdim = intdim;
  r->ptrdim = ptrdim;
  r->offx = OFFSET_NONE;
  if (IS_PTRDIM (lhs, intdim, ptrdim) &&
      ap_linexpr0_is_real (rhs, intdim))
    {
      /* x = y ptr assignment, the expression can be y or NULL */
      pcons0_t *c;
      checked_malloc (c, pcons0_t, sizeof (pcons0_t), 1, return NULL;
                      );
      c->type = OTHER_CONS;
      shape_pcons_of_linexpr (pr, rhs, c, intdim, ptrdim);
      r->info.ptr.y = c->info.ptr.x;
      r->info.ptr.offy = OFFSET_NONE;
      r->op = PA_ASSIGN_PTR;
      free (c);
    }
  else if (IS_INTDIM (lhs, intdim, ptrdim)
           && ap_linexpr0_is_integer (rhs, intdim))
    {
      r->info.data.expr = ap_linexpr0_copy (rhs);
      r->info.data.offsets = NULL;
      r->op = PA_ASSIGN_INT;
    }
  else
    arg_assert (0, return NULL;
                );
  // fill hashable fields
  r->lhs = lhs;
  r->lexpr = rhs;
  r->texpr = NULL;

  //TODO:free r ?
  return shape_passign_add (pr, r);
}

passign0_array_t *
shape_passign_of_linexpr_array (shape_internal_t * pr,
                                ap_dim_t * lhs, ap_linexpr0_t ** rhs,
                                size_t size, size_t intdim, size_t ptrdim)
{
  size_t i;
  passign0_array_t *r = NULL;
  checked_malloc (r, passign0_array_t, sizeof (passign0_array_t), 1,
                  return NULL;
                  );
  r->size = size;
  checked_malloc (r->p, passign0_t *, sizeof (passign0_t *), size,
                  return NULL;
                  );
  for (i = 0; i < size; i++)
    r->p[i] = shape_passign_of_linexpr (pr, lhs[i], rhs[i], intdim, ptrdim);

  return r;
}

/* ===================================================================== */
/* Translation of linear assignments.                                    */

/* ===================================================================== */

passign0_t *
shape_passign_of_texpr (shape_internal_t * pr,
                        ap_dim_t lhs, ap_texpr0_t * rhs,
                        size_t intdim, size_t ptrdim)
{

  passign0_t *r;
#ifndef NDEBUG1
  fprintf (stdout, "===== shape_passign_of_texpr: (over intdim=%zu, ptrdim=%zu) of lhs=%d\n",
           intdim, ptrdim, lhs);
  fflush (stdout);
#endif
  arg_assert (rhs && lhs < intdim + ptrdim, return NULL;
              );

  if ((r = shape_passign_search (pr, lhs, intdim, ptrdim, NULL, rhs)) != NULL)
    return r;

  checked_malloc (r, passign0_t, sizeof (passign0_t), 1, return NULL;
                  );

  /* fill lhs part with a ptr dimension and set dereferencing part */
  r->x = lhs;
  r->intdim = intdim;
  r->ptrdim = ptrdim;
  /* special case of free, alloc, and deref NULL */
  if (shape_is_free_stmt (rhs))
    {
      r->op = PA_FREE;
      r->offx = OFFSET_NONE; // WARNING: no dereference in free
    }
  else if (shape_is_alloc_stmt (rhs))
    {
      r->op = PA_ALLOC; // TODO
      r->offx = OFFSET_NONE;
      r->info.alloc.cst = true; // TODO
      r->info.alloc.size.c = 1; // TODO
    }
  else if ((r->offx = shape_is_field_null_stmt (rhs, intdim, ptrdim)) != OFFSET_NONE)
    {
      r->op = PA_ASSIGN_PTR;
      r->info.ptr.y = NULL_DIM;
      r->info.ptr.offy = OFFSET_NONE;
    }
  else
    {
      pcons0_t *c;
      int deref_offset;

      /* other cases than free and alloc */
      checked_malloc (c, pcons0_t, sizeof (pcons0_t), 1, return NULL;
                      );
      c->type = OTHER_CONS;
      c->intdim = intdim;
      c->ptrdim = ptrdim;
      deref_offset = shape_pcons_of_texpr (pr, rhs, c, intdim, ptrdim);
      if (IS_PTRDIM (lhs, intdim, ptrdim))
        {
          /* pointer assignment */
          if (c->type == DATA_CONS
              && (deref_offset <= ((int) intdim + OFFSET_DATA))
              && deref_offset >= 0) /* data assignment */
            {
              r->op = PA_ASSIGN_INT;
              r->offx = deref_offset;
              r->info.data.expr = c->info.data.cons.linexpr0;
              r->info.data.offsets = c->info.data.offsets;
              c->info.data.cons.linexpr0 = NULL; /* to secure free c  */
              c->info.data.offsets = NULL;
            }
          else if (c->type == DATA_CONS)
            { // data constant, thus NULL
              r->op = PA_ASSIGN_PTR;
              r->offx = deref_offset;
              r->info.ptr.y = NULL_DIM;
              r->info.ptr.offy = OFFSET_NONE;
            }
          else
            { // ptr expression
              r->op = PA_ASSIGN_PTR;
              r->offx = deref_offset;
              r->info.ptr.y = c->info.ptr.x;
              r->info.ptr.offy = c->info.ptr.offx;
            }
        }
      else if (IS_INTDIM (lhs, intdim, ptrdim))
        {
          arg_assert (deref_offset == OFFSET_NONE && c->type == DATA_CONS, return NULL;
                      );
          r->op = PA_ASSIGN_INT;
          r->offx = deref_offset;
          r->info.data.expr = c->info.data.cons.linexpr0;
          r->info.data.offsets = c->info.data.offsets;
          c->info.data.cons.linexpr0 = NULL; /* to secure free c  */
          c->info.data.offsets = NULL;
        }
      else
        arg_assert (0, return NULL;
                    );
      free (c);
    }
  r->lhs = lhs;
  r->lexpr = NULL;
  r->texpr = rhs;
  //TODO:free r ?

  return shape_passign_add (pr, r);
}

passign0_array_t *
shape_passign_of_texpr_array (shape_internal_t * pr,
                              ap_dim_t * lhs, ap_texpr0_t ** rhs, size_t size,
                              size_t intdim, size_t ptrdim)
{
  size_t i;
  passign0_array_t *r = NULL;
  checked_malloc (r, passign0_array_t, sizeof (passign0_array_t), 1,
                  return NULL;
                  );
  r->size = size;
  checked_malloc (r->p, passign0_t *, sizeof (passign0_t *), size,
                  return NULL;
                  );
  for (i = 0; i < size; i++)
    r->p[i] = shape_passign_of_texpr (pr, lhs[i], rhs[i], intdim, ptrdim);

  return r;
}


/* ===================================================================== */
/* Apron expressions for Shapes: High level                              */
/* ===================================================================== */

/* build simple tree expressions */

ap_linexpr0_t *
shape_linexpr_null (size_t intdim, size_t ptrdim)
{

  return shape_linexpr_x_y_l (0, 0, 0, 0, 0, 0, intdim, ptrdim);
}

ap_texpr0_t *
shape_texpr_null (size_t intdim, size_t ptrdim)
{

  return ap_texpr0_cst_scalar_double (0.0);
}

ap_linexpr0_t *
shape_linexpr_x (size_t x, size_t intdim, size_t ptrdim)
{

  return shape_linexpr_x_y_l (1, x, 0, 0, 0, 0, intdim, ptrdim);
}

ap_texpr0_t *
shape_texpr_x_next (size_t x, size_t intdim, size_t ptrdim)
{
  // correct only if x >= intdim
  if (IS_PTRDIM (x, intdim, ptrdim))
    return ap_texpr0_unop (AP_TEXPR_CAST,
                           ap_texpr0_binop (AP_TEXPR_MUL,
                                            ap_texpr0_dim (x),
                                            ap_texpr0_cst_scalar_int (OFFSET_NEXT),
                                            AP_RTYPE_REAL, AP_RDIR_ZERO),
                           AP_RTYPE_REAL, AP_RDIR_ZERO);
  else
    return NULL;
}

ap_texpr0_t *
shape_texpr_x_data (size_t x, size_t intdim, size_t ptrdim)
{
  // correct only if x is a ptr var
  if (IS_PTRDIM (x, intdim, ptrdim))
    return ap_texpr0_unop (AP_TEXPR_CAST,
                           ap_texpr0_binop (AP_TEXPR_MUL,
                                            ap_texpr0_dim (x),
                                            ap_texpr0_cst_scalar_int (OFFSET_DATA),
                                            AP_RTYPE_REAL, AP_RDIR_ZERO),
                           AP_RTYPE_INT, AP_RDIR_ZERO);

  else
    return NULL;
}

ap_texpr0_t *
shape_texpr_x_len (size_t x, size_t intdim, size_t ptrdim)
{
  if (IS_PTRDIM (x, intdim, ptrdim))
    return ap_texpr0_unop (AP_TEXPR_CAST,
                           ap_texpr0_dim (x),
                           AP_RTYPE_QUAD, AP_RDIR_UP);
  else
    return NULL;
}

ap_texpr0_t *
shape_texpr_x_ucons (size_t x, size_t intdim, size_t ptrdim)
{
  if (IS_PTRDIM (x, intdim, ptrdim))
    return ap_texpr0_unop (AP_TEXPR_CAST,
                           ap_texpr0_dim (x),
                           AP_RTYPE_QUAD, AP_RDIR_DOWN);
  else
    return NULL;
}

ap_texpr0_t *
shape_texpr_deref_next (ap_texpr0_t * t)
{
  // correct if t != NULL
  // cast_r,0(t *_r,0 next)
  if (t)
    return ap_texpr0_unop (AP_TEXPR_CAST,
                           ap_texpr0_binop (AP_TEXPR_DIV,
                                            t,
                                            ap_texpr0_cst_scalar_int (OFFSET_NEXT),
                                            AP_RTYPE_REAL, AP_RDIR_ZERO),
                           AP_RTYPE_REAL, AP_RDIR_ZERO);
  else
    return NULL;
}

ap_texpr0_t *
shape_texpr_deref_data (ap_texpr0_t * t)
{
  // correct if t != NULL
  // cast_r,0(t *_i,0 data)
  if (t)
    return ap_texpr0_unop (AP_TEXPR_CAST,
                           ap_texpr0_binop (AP_TEXPR_DIV,
                                            t,
                                            ap_texpr0_cst_scalar_int (OFFSET_DATA),
                                            AP_RTYPE_INT, AP_RDIR_ZERO),
                           AP_RTYPE_REAL, AP_RDIR_ZERO);
  else
    return NULL;
}

ap_texpr0_t *
shape_texpr_x_at_i (size_t x, size_t i, size_t intdim, size_t ptrdim)
{
    return ap_texpr0_unop (AP_TEXPR_CAST,
                           ap_texpr0_binop (AP_TEXPR_MUL,
                                            ap_texpr0_dim(x),
                                            ap_texpr0_dim(i),
                                            AP_RTYPE_INT, AP_RDIR_ZERO),
                           AP_RTYPE_REAL, AP_RDIR_NEAREST);
}

ap_lincons0_t
shape_lincons_same_x_y (size_t x, size_t y, size_t intdim, size_t ptrdim)
{
  if (x != NULL_DIM && y != NULL_DIM)
    return shape_lincons_x_y_l (AP_CONS_EQ, 1, x, -1, y, 0, 0, intdim, ptrdim);
  if (x != NULL_DIM && y == NULL_DIM)
    return shape_lincons_x_cst (AP_CONS_EQ, x, 0, intdim, ptrdim);
  if (x == NULL_DIM && y != NULL_DIM)
    return shape_lincons_x_cst (AP_CONS_EQ, y, 0, intdim, ptrdim);

  /* 0 = 0 */
  return ap_lincons0_make (AP_CONS_EQ, ap_linexpr0_alloc (AP_LINEXPR_DENSE, 0), NULL);
}

ap_tcons0_t
shape_tcons_same_x_y (size_t x, bool nx, size_t y, bool ny,
                      size_t intdim, size_t ptrdim)
{

  return shape_tcons_x_y (AP_CONS_EQ, x, nx, y, ny, intdim, ptrdim);

}

bool
shape_tcons_is_true (ap_tcons0_t * t)
{
  // Recognize only 0=0

  return !t || (t->constyp == AP_CONS_EQ &&
                t->texpr0->discr == AP_TEXPR_CST &&
                ap_coeff_zero (&(t->texpr0->val.cst)));
}

bool
shape_tcons_array_is_true (ap_tcons0_array_t * array)
{
  if (!array || !array->size)
    return true;
  size_t i;
  bool r = true;
  for (i = 0; i < array->size && r; i++)
    if (!shape_tcons_is_true (&array->p[i]))
      r = false;

  return r;
}

ap_lincons0_t
shape_lincons_diff_x_y (size_t x, size_t y, size_t intdim, size_t ptrdim)
{
  if (y == NULL_DIM)
    return shape_lincons_x_cst (AP_CONS_DISEQ, x, 0, intdim, ptrdim);

  else
    return shape_lincons_x_y_l (AP_CONS_DISEQ, 1, x, -1, y, 0, 0, intdim, ptrdim);
}

ap_tcons0_t
shape_tcons_diff_x_y (size_t x, bool nx, size_t y, bool ny, size_t intdim, size_t ptrdim)
{

  return shape_tcons_x_y (AP_CONS_DISEQ, x, nx, y, ny, intdim, ptrdim);

}

ap_lincons0_t
shape_lincons_reach_x_y (size_t x, size_t y, size_t intdim, size_t ptrdim)
{
  if (y == NULL_DIM)
    return shape_lincons_x_cst (AP_CONS_SUP, x, 0, intdim, ptrdim);

  else
    return shape_lincons_x_y_l (AP_CONS_SUP, 1, x, -1, y, 0, 0, intdim, ptrdim);
}

ap_tcons0_t
shape_tcons_reach_x_y (size_t x, bool nx, size_t y, bool ny, size_t intdim, size_t ptrdim)
{

  return shape_tcons_x_y (AP_CONS_SUP, x, nx, y, ny, intdim, ptrdim);

}

ap_texpr0_t *
shape_texpr_x_alloc (size_t x, size_t intdim, size_t ptrdim)
{
  // correct if x is a ptr variable
  if (IS_PTRDIM (x, intdim, ptrdim))
    return ap_texpr0_unop (AP_TEXPR_SQRT,
                           ap_texpr0_dim (x),
                           AP_RTYPE_SINGLE, AP_RDIR_NEAREST);

  else
    return NULL;
}

ap_texpr0_t *
shape_texpr_x_free (size_t x, size_t intdim,
                    size_t ptrdim)
{
  // correct if x is a ptr variable
  if (IS_PTRDIM (x, intdim, ptrdim))
    return ap_texpr0_unop (AP_TEXPR_SQRT, ap_texpr0_dim (x),
                           AP_RTYPE_SINGLE, AP_RDIR_ZERO);

  else
    return NULL;
}

/* Random constructors: High level */
const char *exprname[] = {"ptr ", "reach ", "lindata ",
  "next ", "deref ", "data ", "deref_data"};

/* Build randomly linear expressions. */
ap_linexpr0_t *
shape_linexpr_random (shape_internal_t * pr, exprmode_t mode, size_t intdim,
                      size_t ptrdim)
{
  int dim = intdim + ptrdim;
  ap_linexpr0_t *l = ap_linexpr0_alloc (AP_LINEXPR_DENSE, dim);
  long int x, y;
  if (mode <= expr_reach)
    {
      /* x */
      x = RANDOM_PTRDIM (pr, intdim, ptrdim);
      ap_coeff_set_scalar_int (l->p.coeff + x, 1);
      /* - y */
      do
        {
          y = RANDOM_PTRDIM_OR_NULL (pr, intdim, ptrdim);
        }
      while (y == x);
      if (mode >= expr_reach)
        ap_coeff_set_scalar_int (l->p.coeff + y, -1);
    }
  else
    {
      size_t i;
      if (lrand48 () % 20 > 12)
        {
          for (i = 0; i < intdim; i++)
            {
              int n = lrand48 () % 20 - 10;
              ap_coeff_set_scalar_int (l->p.coeff + i, n);
            }
        }
      else
        {
          for (i = 0; i < intdim; i++)
            {

              int n = lrand48 () % 3 - 1;
              ap_coeff_set_scalar_int (l->p.coeff + i, n);
            }
        }
    }
  return l;
}

/* Build randomly tree expressions */
ap_texpr0_t *
shape_texpr_random (shape_internal_t * pr, exprmode_t mode, size_t intdim,
                    size_t ptrdim)
{
  size_t x;

  switch (mode)
    {
    case expr_next: /* generate x->next */
      /* choose a variable x randomly, either if it points to null */
      x = RANDOM_PTRDIM (pr, intdim, ptrdim);
      return shape_texpr_x_next (x, intdim, ptrdim);

    case expr_deref: /* x / next or x->next/next */
      /* choose a variable x randomly, either if it points to null */
      x = RANDOM_PTRDIM_OR_NULL (pr, intdim, ptrdim);
      if (lrand48 () % 10 < 2 && x != NULL_DIM)
        return shape_texpr_deref_next (shape_texpr_x_next (x, intdim, ptrdim));
      else
        return shape_texpr_deref_next (ap_texpr0_dim (x));

    case expr_reach:
      { /* generate x[* next] - y[* next] */
        size_t y;
        ap_texpr0_t *xn, *yn;
        /* choose a variable x randomly, either if it points to null */
        x = RANDOM_PTRDIM (pr, intdim, ptrdim);
        /*
         * choose a variable y randomly but different from x, either if it
         * points to null
         */
        do
          {
            y = RANDOM_PTRDIM (pr, intdim, ptrdim);
          }
        while (x == y);
        if (lrand48 () % 10 > 8)
          {
            xn = shape_texpr_x_next (x, intdim, ptrdim);
            yn = ap_texpr0_dim (y);
          }
        else
          {
            yn = shape_texpr_x_next (y, intdim, ptrdim);
            xn = ap_texpr0_dim (x);
          }
        return ap_texpr0_binop (AP_TEXPR_SUB, xn, yn,
                                AP_RTYPE_SINGLE, AP_RDIR_NEAREST);
      }

    case expr_data: /* v1 + x->data or cst + x->data */
      {
        ap_texpr0_t *e;
        x = RANDOM_PTRDIM (pr, intdim, ptrdim);
        e = shape_texpr_x_data (x, intdim, ptrdim);
        switch (lrand48 () % 4)
          {
          case 0: /* cst * x->data */
            e = ap_texpr0_binop (AP_TEXPR_MUL,
                                 ap_texpr0_cst_scalar_double ((double)
                                                              (lrand48 () %
                                                               5 + 1)), e,
                                 AP_RTYPE_INT, AP_RDIR_NEAREST);
            break;
          case 1: /* cst + cst * x*data */
            e = ap_texpr0_binop (AP_TEXPR_MUL,
                                 ap_texpr0_cst_scalar_double ((double)
                                                              (lrand48 () %
                                                               5 + 1)), e,
                                 AP_RTYPE_INT, AP_RDIR_NEAREST);
            e =
                    ap_texpr0_binop (AP_TEXPR_ADD,
                                     ap_texpr0_cst_scalar_double ((double)
                                                                  (lrand48 () % 5 +
                                                                   1)), e,
                                     AP_RTYPE_INT, AP_RDIR_NEAREST);
            break;
          case 2: /* v + x->data */
            e =
                    ap_texpr0_binop (AP_TEXPR_ADD,
                                     ap_texpr0_dim (lrand48 () % intdim), e,
                                     AP_RTYPE_INT, AP_RDIR_NEAREST);
            break;
          default: /* v + cst*x->data */
            e = ap_texpr0_binop (AP_TEXPR_MUL,
                                 ap_texpr0_cst_scalar_double ((double)
                                                              (lrand48 () %
                                                               5 + 1)), e,
                                 AP_RTYPE_INT, AP_RDIR_NEAREST);
            e =
                    ap_texpr0_binop (AP_TEXPR_ADD,
                                     ap_texpr0_dim (lrand48 () % intdim), e,
                                     AP_RTYPE_INT, AP_RDIR_NEAREST);

            break;
          }
        return e;
      }
    default: /* expr_data / data */
      return shape_texpr_deref_data (shape_texpr_random (pr, expr_data, intdim,
                                                         ptrdim));
    }
  return NULL;
}

/* ===================================================================== */
/* Apron expressions for Shapes: Low level, take care with them!         */
/* ===================================================================== */

/** Build constraint x op cst */
ap_lincons0_t
shape_lincons_x_cst (ap_constyp_t op, size_t x, int cst,
                     size_t intdim, size_t ptrdim)
{
  int dim = intdim + ptrdim;
  ap_linexpr0_t *l = ap_linexpr0_alloc (AP_LINEXPR_DENSE, dim);
  ap_scalar_t *c = ap_scalar_alloc ();
  ap_coeff_set_scalar_int (l->p.coeff + x, 1);
  ap_scalar_set_int (c, cst);

  return ap_lincons0_make (op, l, c);
}

/** Build expression x - y */
ap_linexpr0_t *
shape_linexpr_x_y_l (int cx, size_t x, int cy, size_t y, int cl, size_t l,
                     size_t intdim, size_t ptrdim)
{
  int dim = intdim + ptrdim;
  ap_linexpr0_t *lexpr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, dim);
  ap_coeff_set_scalar_int (lexpr->p.coeff + x, cx);
  ap_coeff_set_scalar_int (lexpr->p.coeff + y, cy);
  ap_coeff_set_scalar_int (lexpr->p.coeff + l, cl);

  return lexpr;
}

ap_lincons0_t
shape_lincons_x_y_l (ap_constyp_t op, int cx, size_t x, int cy, size_t y,
                     int cl, size_t l, size_t intdim, size_t ptrdim)
{

  return ap_lincons0_make (op,
                           shape_linexpr_x_y_l (cx, x, cy, y, cl, l, intdim,
                                                ptrdim), NULL);
}

/* Generate linear constraint:
 * cx*x + cy*y + cv*v + cst op codeop
 */
ap_lincons0_t
shape_lincons_x_y_v_cst (ap_constyp_t op, size_t code, int cx, size_t x,
                         int cy, size_t y, int cv, size_t v, int cst,
                         size_t intdim, size_t ptrdim)
{
  int dim = intdim + ptrdim;
  ap_linexpr0_t *lexpr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, dim);
  ap_coeff_set_scalar_int (lexpr->p.coeff + x, cx);
  ap_coeff_set_scalar_int (lexpr->p.coeff + y, cy);
  ap_coeff_set_scalar_int (lexpr->p.coeff + v, cv);
  ap_linexpr0_set_cst_scalar_int (lexpr, cst);
  ap_scalar_t *codeop = ap_scalar_alloc ();
  ap_scalar_set_int (codeop, code);

  return ap_lincons0_make (op, lexpr, codeop);
}

/** Build expression x(*next)- y(*next) */
ap_texpr0_t *
shape_texpr_x_y (size_t x, bool nx, size_t y, bool ny, size_t intdim,
                 size_t ptrdim)
{
  ap_texpr0_t* tx = NULL;
  ap_texpr0_t* ty = NULL;
  if (x != NULL_DIM)
    tx = shape_texpr_x_n_next (x, (nx) ? 1 : 0, intdim, ptrdim);
  else
    tx = shape_texpr_null (intdim, ptrdim);

  if (y != NULL_DIM)
    ty = shape_texpr_x_n_next (y, (ny) ? 1 : 0, intdim, ptrdim);
  else
    ty = shape_texpr_null (intdim, ptrdim);

  return ap_texpr0_binop (AP_TEXPR_SUB, tx, ty, AP_RTYPE_INT, AP_RDIR_RND);
}

ap_tcons0_t
shape_tcons_x_y (ap_constyp_t op, size_t x, bool offx, size_t y, bool offy,
                 size_t intdim, size_t ptrdim)
{

  return ap_tcons0_make (op,
                         shape_texpr_x_y (x, offx, y, offy, intdim,
                                          ptrdim), NULL);
}

/** Build expression x(->data)+y->data+v+cst */
ap_texpr0_t *
shape_texpr_x_y_v_cst_data (int cx, size_t x, int cy, size_t y,
                            int cv, size_t v, int cst, bool rhs,
                            size_t intdim, size_t ptrdim)
{
  ap_texpr0_t *t;
  if (cx)
    t = ap_texpr0_binop (AP_TEXPR_MUL, ap_texpr0_cst_scalar_int (cx),
                         shape_texpr_x_data (x, intdim, ptrdim), AP_RTYPE_INT, AP_RDIR_RND);
  else
    t = NULL;
  if (cy)
    {
      ap_texpr0_t *ty =
              ap_texpr0_binop (AP_TEXPR_MUL, ap_texpr0_cst_scalar_int (cy),
                               shape_texpr_x_data (y, intdim, ptrdim),
                               AP_RTYPE_INT, AP_RDIR_RND);
      t =
              (t == NULL) ? ty : ap_texpr0_binop (AP_TEXPR_ADD, t, ty, AP_RTYPE_INT,
                                                  AP_RDIR_RND);
    }
  if (cv)
    {
      ap_texpr0_t *tv =
              ap_texpr0_binop (AP_TEXPR_MUL, ap_texpr0_cst_scalar_int (cv),
                               ap_texpr0_dim (v),
                               AP_RTYPE_INT, AP_RDIR_RND);
      t =
              (t == NULL) ? tv : ap_texpr0_binop (AP_TEXPR_ADD, t, tv, AP_RTYPE_INT,
                                                  AP_RDIR_RND);
    }
  if (cst || !t)
    {
      ap_texpr0_t *tcst = ap_texpr0_cst_scalar_int (cst);
      t =
              (t == NULL) ? tcst : ap_texpr0_binop (AP_TEXPR_ADD, t, tcst,
                                                    AP_RTYPE_INT, AP_RDIR_RND);
    }
  if (t && rhs)
    t = shape_texpr_deref_data (t);

  return t;
}

/** Build expression x(->next) */
ap_texpr0_t *
shape_texpr_x_n_next (size_t x, int offx, size_t intdim, size_t ptrdim)
{
  ap_texpr0_t *tx;
  if (offx > 0)
    tx = shape_texpr_x_next (x, intdim, ptrdim);
  else
    {
      tx = ap_texpr0_dim (x);

      if (offx < 0)
        tx = shape_texpr_deref_next (tx);
    }

  return tx;
}

/** Build expression x(->next)- cst */
ap_texpr0_t *
shape_texpr_x_cst (size_t x, bool offx, int cst, size_t intdim,
                   size_t ptrdim)
{
  ap_texpr0_t *tx, *ty;
  tx = shape_texpr_x_n_next (x, (offx) ? 1 : 0, intdim, ptrdim);
  ty = ap_texpr0_cst_scalar_int (cst);

  return ap_texpr0_binop (AP_TEXPR_SUB, tx, ty, AP_RTYPE_INT, AP_RDIR_RND);
}

ap_tcons0_t
shape_tcons_x_cst (ap_constyp_t op, size_t x, bool offx, int cst,
                   size_t intdim, size_t ptrdim)
{

  return ap_tcons0_make (op,
                         shape_texpr_x_cst (x, offx, cst, intdim,
                                            ptrdim), NULL);
}

ap_texpr0_t *
shape_texpr_singleton (int n)
{

  return ap_texpr0_cst_scalar_int (n);
}

/* Create tree expression src/new */
ap_texpr0_t *
shape_texpr_split (size_t src, size_t new)
{
  ap_texpr0_t *expr_src = ap_texpr0_dim ((ap_dim_t) src);
  ap_texpr0_t *expr_new = ap_texpr0_dim ((ap_dim_t) new);

  return ap_texpr0_binop (AP_TEXPR_DIV, expr_src, expr_new, AP_RTYPE_REAL,
                          AP_RDIR_RND);
}

/* Create an array of texpr of size from the array nodes[size][max] given by
 * n_1 / (n2 /  .... / n_max)
 */
ap_texpr0_t **
shape_texpr_merge (size_t ** nodes, size_t size, size_t max)
{
  assert (size > 0 && max > 0);
  size_t i, j;
  ap_texpr0_t **r = (ap_texpr0_t **) malloc (size * sizeof (ap_texpr0_t *));
  for (i = 0; i < size; i++)
    { // create r[i] from nodes[i]
      ap_texpr0_t *tmp = ap_texpr0_dim ((ap_dim_t) nodes[i][max - 1]);
      if (max >= 2)

        for (j = max - 2; j != ((size_t) - 1); j--)
          tmp =
                ap_texpr0_binop (AP_TEXPR_DIV,
                                 ap_texpr0_dim ((ap_dim_t) nodes[i][j]), tmp,
                                 AP_RTYPE_REAL, AP_RDIR_RND);
      r[i] = tmp;
    }
  return r;
}

/* Modify e wrt v2n[ptrdim]
 *
 */
ap_linexpr0_t *
shape_linexpr_of_node (shape_internal_t * pr,
                       ap_linexpr0_t * e,
                       size_t * v2n,
                       size_t size, size_t intdim, size_t ptrdim)
{

  ap_linexpr0_t *l = ap_linexpr0_alloc (AP_LINEXPR_DENSE, intdim + size);
  size_t i, d;
  ap_coeff_t *coeff;

  ap_linexpr0_ForeachLinterm (e, i, d, coeff)
  {
    if (!ap_coeff_zero (coeff))
      {
        arg_assert ((coeff->discr == AP_COEFF_SCALAR), return NULL;
                    );
        if (IS_PTRDIM (d, intdim, ptrdim))
          {
            size_t v = DIM2PTR (d, intdim);
            size_t n = 0;
            arg_assert (!IS_NULLDIM (v), return NULL;
                        );
            n = v2n[v];
            ap_coeff_set (l->p.coeff + intdim + n, coeff);
          }
        else
          ap_coeff_set (l->p.coeff + d, coeff);
      }
  }
  /* constant */
  coeff = ap_coeff_alloc (AP_COEFF_SCALAR);
  ap_linexpr0_get_cst (coeff, e);
  ap_linexpr0_set_cst (l, coeff);
  ap_coeff_free (coeff);
#ifndef NDEBUG2
  fprintf (stdout, "\n====shape_linexpr_of_node: ");
  ap_linexpr0_fprint (stdout, l, NULL);
  fprintf (stdout, "\n");
#endif
  return l;
}

/* Copy the src dimperm in dst from begin */
void
shape_dimperm_copy (ap_dimperm_t * dst, size_t begin, ap_dimperm_t * src)
{
  size_t i;

  for (i = 0; i < src->size && (i + begin) < dst->size; i++)
    dst->dim[begin + i] = begin + src->dim[i];

#ifndef NDEBUG2
  fprintf (stdout, "\n====shape_dimperm_copy: ");
  fprintf (stdout, "\n===src: ");
  ap_dimperm_fprint (stdout, src);
  fprintf (stdout, "\n===dst: ");
  ap_dimperm_fprint (stdout, dst);
  fprintf (stdout, "\n");
#endif

}

bool
shape_dimperm_is_id (ap_dimperm_t * p)
{
  size_t i;
  bool r = true;
  for (i = 0; i < p->size && r; i++)
    if (p->dim[i] != i)
      r = false;

  return r;
}

void
shape_dimchange_sort (ap_dimchange_t * p)
{
  size_t i, j, k, size;
  ap_dim_t d;
  size = p->intdim + p->realdim;
  if (!p || !p->dim || !size)
    return;
  /* insertion sort */
#ifndef NDEBUG2
  fprintf (stdout, "\n====shape_dimchange_sort: [");
  for (i = 0; i < size; i++)
    fprintf (stdout, "%d,", p->dim[i]);
  fprintf (stdout, "]\n");
#endif
  for (i = 1; i < size; i++)
    {
      j = 0;
      while (j != i && p->dim[j] <= p->dim[i])
        j++;
      if (j < i)
        {
          d = p->dim[i];
          for (k = i - 1; k >= j; k--)
            p->dim[k + 1] = p->dim[k];
          p->dim[j] = d;
        }
    }
#ifndef NDEBUG2
  fprintf (stdout, "\n===returns: [");

  for (i = 0; i < size; i++)
    fprintf (stdout, "%d,", p->dim[i]);
  fprintf (stdout, "]\n");
#endif
}


/* ================================================================= */
/* Printing through the Apron interface for Shapes */

/* ================================================================= */

struct name_of_dim_typ
{
  size_t intdim;
  size_t segmdim;
  char** names;
} shape_names;

char*
shape_name_of_dim (size_t dim)
{
  if (dim < shape_names.intdim + shape_names.segmdim
      && shape_names.names != NULL
      && shape_names.names[dim] != NULL)
    return shape_names.names[dim];

  else
    return "dummy_dim";
}

void
shape_init_name_of_dim (size_t intdim, size_t segmdim)
{
  if (shape_names.intdim == intdim
      && shape_names.segmdim == segmdim
      && shape_names.names != NULL)
    /* nothing to do */
    return;

  size_t i, dim;
  if (shape_names.names != NULL)
    {
      dim = shape_names.intdim + shape_names.segmdim;
      for (i = 0; i < dim; i++)
        if (shape_names.names[i])
          free (shape_names.names[i]);
      free (shape_names.names);
      shape_names.names = NULL;
    }
  dim = intdim + segmdim;
  shape_names.intdim = intdim;
  shape_names.segmdim = segmdim;
  shape_names.names =
          (char **) malloc (dim * sizeof (char *));
  for (i = 0; i < intdim; i++)
    {
      shape_names.names[i] = (char *) malloc (5 * sizeof (char));
      sprintf (shape_names.names[i], "x%zu", i);
    }
  for (i = 0; i < segmdim; i++)
    {
      shape_names.names[intdim + i] = (char *) malloc (5 * sizeof (char));
      sprintf (shape_names.names[intdim + i], "n%zu", i);
    }
}

bool shape_print_dot = true;

inline bool
shape_get_print (void)
{
  return shape_print_dot;
}

inline void
shape_set_print_dot ()
{
  shape_print_dot = true;
}

inline void
shape_set_print_smt ()
{
  shape_print_dot = false;
}
