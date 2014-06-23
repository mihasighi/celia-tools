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


#include "shape.h"
#include "shape_internal.h"
#include "apron2shape.h"
#include "shape_macros.h"
#include "ap_generic.h"



/* ============================================================ */
/* Assignement and Substitutions */

/* ============================================================ */

shape_t *
shape_assign_passign_array (shape_internal_t * pr,
                            bool destructive, shape_t * a,
                            passign0_array_t * arr)
{
  ushape_array_t *r = NULL;
  shape_t *rs;
  size_t i;
  for (i = 0; i < a->msize; i++)
    {
      ushape_array_t *rr =
              ushape_assign_passign_array (pr, false, a->m.p[i], arr);
      if (rr)
        {
          r = ushape_array_add_array (pr, true, r, rr);
          ushape_array_init (pr, rr, rr->size);
          free (rr);
        }
    }
  if (!r)
    return NULL;
  /* compute msize */
  for (i = 0; i < r->size; i++)
    if (!r->p[i])
      break;
  checked_malloc (rs, shape_t, sizeof (shape_t), 1, return NULL;
                  );
  rs->m = *r;
  r->p = NULL;
  free (r);
  rs->msize = i;
  rs->set = rs->closed = false;
  rs->intdim = a->intdim;
  rs->realdim = a->realdim;
  if (destructive)
    shape_free_internal (pr, a);
  return rs;
}

shape_t *
shape_assign_linexpr (ap_manager_t * man,
                      bool destructive, shape_t * a,
                      ap_dim_t d, ap_linexpr0_t * expr, shape_t * dest)
{
  if (shape_is_bottom (man, a))
    /* nothing to do */
    return (destructive) ? a : shape_copy (man, a);
  else
    {
      shape_t *b;
      passign0_array_t *op;
      shape_t *r;
      shape_internal_t *pr =
              shape_init_from_manager (man, AP_FUNID_ASSIGN_LINEXPR_ARRAY, 0);
      if (!destructive)
        b = shape_copy_internal (pr, a);
      else
        b = a;
      op =
              shape_passign_of_linexpr_array (pr, &d, &expr, 1, a->intdim,
                                              a->realdim);
#ifndef NDEBUG1
      fprintf (stdout, "****shape_assign_linexpr: with passign=[");
      shape_passign_array_fdump (stdout, op, a->intdim, a->realdim);
      fprintf (stdout, "]\n");
      fflush (stdout);
#endif
      /* go */
      r = shape_assign_passign_array (pr, false, b, op);
      shape_free_internal (pr, b);
      if (dest)
        { /* intersect r with dest */
          shape_t *rr = shape_meet (pr->man, false, r, dest);
          shape_free_internal (pr, r);
          return rr;
        }
#ifndef NDEBUG1
      fprintf (stdout, "****shape_assign_linexpr returns: ");
      shape_fdump (stdout, man, r);
      fflush (stdout);
#endif
      return r;
    }
}

/**
 * Substitute the actual in tdim by the corresponding formals in tnew.
 * @param size length of tnew and tdim
 * @return 
 */
shape_t*
shape_substitute_actuals (shape_internal_t* pr,
                          bool destructive,
                          shape_t* a,
                          ap_dim_t* tdim,
                          ap_dim_t* tnew,
                          size_t size)
{
  if (!a || shape_is_top (pr->man, a))
    return (destructive) ? a : shape_copy (pr->man, a);
  /* go */
  ushape_array_t *r = NULL;
  shape_t *rs;
  size_t i;
  for (i = 0; i < a->msize; i++)
    {
      ushape_array_t *rr = ushape_substitute_actuals (pr, a->m.p[i], tdim, tnew, size);
      if (rr)
        {
          r = ushape_array_add_array (pr, true, r, rr);
          ushape_array_init (pr, rr, rr->size);
          free (rr);
        }
    }
  if (!r)
    return NULL;
  /* compute msize */
  for (i = 0; i < r->size; i++)
    if (!r->p[i])
      break;
  checked_malloc (rs, shape_t, sizeof (shape_t), 1, return NULL;
                  );
  rs->m = *r;
  r->p = NULL;
  free (r);
  rs->msize = i;
  rs->set = rs->closed = false;
  rs->intdim = a->intdim;
  rs->realdim = a->realdim;
  if (destructive)
    shape_free_internal (pr, a);

  // change flag meet_algo if -1    
  if (pr->meet_algo < 0)
    shape_approximate (pr->man, a, 0);

  return rs;
}

shape_t *
shape_substitute_passign_array (shape_internal_t * pr,
                                bool destructive, shape_t * a,
                                passign0_array_t * arr)
{
  ushape_array_t *r = NULL;
  shape_t *rs;
  size_t i;
  for (i = 0; i < a->msize; i++)
    {
      ushape_array_t *rr =
              ushape_substitute_passign_array (pr, a->m.p[i], arr);
      if (rr)
        {
          r = ushape_array_add_array (pr, true, r, rr);
          ushape_array_init (pr, rr, rr->size);
          free (rr);
        }
    }
  if (!r)
    return NULL;
  /* compute msize */
  for (i = 0; i < r->size; i++)
    if (!r->p[i])
      break;
  checked_malloc (rs, shape_t, sizeof (shape_t), 1, return NULL;
                  );
  rs->m = *r;
  r->p = NULL;
  free (r);
  rs->msize = i;
  rs->set = rs->closed = false;
  rs->intdim = a->intdim;
  rs->realdim = a->realdim;
  if (destructive)
    shape_free_internal (pr, a);
  return rs;
}

/* TODO: priority 0 */

/* used for pre-image computation */
shape_t *
shape_substitute_linexpr (ap_manager_t * man,
                          bool destructive, shape_t * a,
                          ap_dim_t d, ap_linexpr0_t * expr, shape_t * dest)
{
  shape_internal_t *pr =
          shape_init_from_manager (man, AP_FUNID_SUBSTITUTE_LINEXPR_ARRAY, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a;
}

shape_t *
shape_assign_linexpr_array (ap_manager_t * man,
                            bool destructive, shape_t * a,
                            ap_dim_t * tdim,
                            ap_linexpr0_t ** texpr,
                            size_t size, shape_t * dest)
{
  if (shape_is_bottom (man, a))
    /* nothing to do */
    return (destructive) ? a : shape_copy (man, a);
  else
    {
      shape_t *b;
      passign0_array_t *op;
      shape_t *r;
      shape_internal_t *pr =
              shape_init_from_manager (man, AP_FUNID_ASSIGN_LINEXPR_ARRAY, 0);
      if (!destructive)
        b = shape_copy_internal (pr, a);
      else
        b = a;
      op =
              shape_passign_of_linexpr_array (pr, tdim, texpr, size, a->intdim,
                                              a->realdim);
#ifndef NDEBUG1
      fprintf (stdout, "\n****shape_assign_linexpr_array: with assign=[");
      size_t i;
      for (i = 0; i < size; i++)
        {
          fprintf (stdout, "x%d := ", tdim[i]);
          ap_linexpr0_fprint (stdout, texpr[i], NULL);
          fprintf (stdout, ",  ");
        }
      fprintf (stdout, "]\n");
      fflush (stdout);
      /*
      fprintf (stdout, " (i.e.  ");
      shape_passign_array_fdump (stdout, op, a->intdim, a->realdim);
      fprintf (stdout, ") \n\t on ");
       */
      fprintf (stdout, " on \n");
      shape_fdump (stdout, man, b);
      fprintf (stdout, "\n");
      fflush (stdout);
#endif
      /* go */
      r = shape_assign_passign_array (pr, false, b, op);
      shape_free_internal (pr, b);
      if (dest)
        { /* intersect r with dest */
          shape_t *rr = shape_meet (pr->man, false, r, dest);
          shape_free_internal (pr, r);
          return rr;
        }
#ifndef NDEBUG1
      fprintf (stdout, "\n****shape_assign_linexpr_array returns: ");
      shape_fdump (stdout, man, r);
      fprintf (stdout, "\n");
      fflush (stdout);
#endif
      return r;
    }
}

/* TODO: priority 0 */

/* used for pre-image computation */
shape_t *
shape_substitute_linexpr_array (ap_manager_t * man,
                                bool destructive, shape_t * a,
                                ap_dim_t * tdim,
                                ap_linexpr0_t ** texpr,
                                size_t size, shape_t * dest)
{
  if (size == 1)
    return shape_substitute_linexpr (man, destructive, a, tdim[0], texpr[0],
                                     dest);

  shape_internal_t *pr =
          shape_init_from_manager (man, AP_FUNID_SUBSTITUTE_TEXPR_ARRAY, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a;
}

shape_t *
shape_assign_texpr_array (ap_manager_t * man,
                          bool destructive, shape_t * a,
                          ap_dim_t * tdim,
                          ap_texpr0_t ** texpr, size_t size, shape_t * dest)
{
  if (shape_is_bottom (man, a))
    /* nothing to do */
    return (destructive) ? a : NULL;
  else
    {
      passign0_array_t *op;
      shape_t *r;
      shape_internal_t *pr =
              shape_init_from_manager (man, AP_FUNID_ASSIGN_TEXPR_ARRAY, 0);
      op = shape_passign_of_texpr_array (pr, tdim, texpr, size, a->intdim,
                                         a->realdim);
#ifndef NDEBUG1
      fprintf (stdout, "\n****shape_assign_texpr_array: with assign=[");
      size_t i;
      for (i = 0; i < size; i++)
        {
          fprintf (stdout, "x%d := ", tdim[i]);
          ap_texpr0_fprint (stdout, texpr[i], NULL);
          fprintf (stdout, ",  ");
        }
      fprintf (stdout, "]\n ");
      fflush (stdout);
      // shape_passign_array_fdump (stdout, op, a->intdim, a->realdim);
      fprintf (stdout, " on ");
      shape_fdump (stdout, man, a);
      fprintf (stdout, "\n");
      fflush (stdout);
#endif
      /* go */
      r = shape_assign_passign_array (pr, false, a, op);
      if (dest)
        { /* intersect r with dest */
          shape_t *rr = shape_meet (pr->man, false, r, dest);
          shape_free_internal (pr, r);
          return rr;
        }
      if (destructive)
        shape_free_internal (pr, a);
#ifndef NDEBUG1
      fprintf (stdout, "\n****shape_assign_texpr_array returns:");
      shape_fdump (stdout, man, r);
      fprintf (stdout, "\n");
      fflush (stdout);
#endif
      return r;
    }
}

/** Used in:
 *  - (for flag meet_algo==-1) implementation of return edge in ICFG.
 *    It is not a simple substitution for ptrdim since edges
 *    pointing to the nodes labeled by tdim are changed to the nodes
 *    in texpr.
 *  - (for flag meet_algo>=0)  implementation of weakest precondition.
 *  @param texpr  is only a variable
 */
shape_t *
shape_substitute_texpr_array (ap_manager_t * man,
                              bool destructive, shape_t * a,
                              ap_dim_t * tdim,
                              ap_texpr0_t ** texpr,
                              size_t size, shape_t * dest)
{
  shape_internal_t *pr =
          shape_init_from_manager (man, AP_FUNID_ASSIGN_TEXPR_ARRAY, 0);
  shape_t *r = NULL;

  if (shape_is_bottom (man, a))
    /* nothing to do */
    return (destructive) ? a : NULL;
  else if (pr->meet_algo < 0)
    {
      /****** substitute actuals with ini formals ******/
      ap_dim_t* tnew;
      size_t i;
      // build array of new dimensions from texpr
      checked_malloc (tnew, ap_dim_t, sizeof (ap_dim_t), size, return NULL;);
      for (i = 0; i < size; i++)
        if (texpr[i]->discr == AP_TEXPR_DIM)
          tnew[i] = texpr[i]->val.dim;
        else
          {
            ERROR ("substitution without dimension",);
            tnew[i] = tdim[i];
          }
#ifndef NDEBUG1
      size_t hedge = 0;
      FILE* fedge = NULL;
      if ((fedge = fopen ("hedge.txt", "r")) != NULL)
        {
          fprintf (stdout, "\n****shape read hedge:\n");
          hedge = fgetc (fedge);
          fprintf (stdout, "\n****shape read hedge: %zu\n", hedge);
          fclose (fedge);
        }
      fprintf (stdout, "\n****shape_substitute_texpr_array: with assign=[");
      for (i = 0; i < size; i++)
        {
          ap_texpr0_fprint (stdout, texpr[i], NULL);
          fprintf (stdout, " (x%d) / x%d, ", tnew[i], tdim[i]);
        }
      fprintf (stdout, "]\n ");
      fflush (stdout);
      // shape_passign_array_fdump (stdout, op, a->intdim, a->realdim);
      fprintf (stdout, " on ");
      shape_fdump (stdout, man, a);
      fprintf (stdout, "\n");
      fflush (stdout);
#endif
      /* go */
      r = shape_substitute_actuals (pr, false, a, tdim, tnew, size);
      // meet_algo flag is set above
    }
  else
    {
      passign0_array_t *op =
              shape_passign_of_texpr_array (pr, tdim, texpr, size, a->intdim,
                                            a->realdim);
#ifndef NDEBUG1
      fprintf (stdout, "\n****shape_substitute_texpr_array: with assign=[");
      size_t i;
      for (i = 0; i < size; i++)
        {
          fprintf (stdout, "x%d := ", tdim[i]);
          ap_texpr0_fprint (stdout, texpr[i], NULL);
          fprintf (stdout, ",  ");
        }
      fprintf (stdout, "]\n ");
      fflush (stdout);
      // shape_passign_array_fdump (stdout, op, a->intdim, a->realdim);
      fprintf (stdout, " on ");
      shape_fdump (stdout, man, a);
      fprintf (stdout, "\n");
      fflush (stdout);
#endif
      /* go */
      r = shape_substitute_passign_array (pr, false, a, op);
    }
  if (dest)
    { /* intersect r with dest */
      shape_t *rr = shape_meet (pr->man, false, r, dest);
      shape_free_internal (pr, r);
      return rr;
    }
  if (destructive)
    shape_free_internal (pr, a);
#ifndef NDEBUG1
  fprintf (stdout, "\n****shape_substitute_texpr_array returns:");
  shape_fdump (stdout, man, r);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  return r;

}
