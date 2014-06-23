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


#include <stdio.h>
#include "shape_options.h"
#include "ushape_internal.h"
#include "shape.h"
#include "shape_internal.h"
#include "shape_macros.h"
#include "ap_generic.h"
#include "box.h"
#include "oct.h"
#include "pk.h"
#include "lsum.h"
#include "lsum_internal.h"
#include "mset.h"
#include "mset_internal.h"
#include "ucons.h"



/* ********************************************************************** */
/* Shapes */
/* ********************************************************************** */

/* for pretty printing */
size_t ushape_number;
/* ============================================================ */
/* Internal management */
/* ============================================================ */

/*
 * generic allocation routine, returns a shape with msize at 0 but size
 * allocated NULL ushapes
 */
inline shape_t *
shape_alloc_internal (shape_internal_t * pr,
                      size_t size)
{

  shape_t * r;

  size_t i;

  checked_malloc (r, shape_t, sizeof (shape_t), 1, return NULL;
                  );

  r->m.p = NULL;

  ushape_array_init (pr, &r->m, size);

  r->msize = 0;

  r->set = r->closed = true;

  r->intdim = pr->intdim;

  r->realdim = pr->realdim;

  return r;

}

/* returns a singleton set with a top ushape */
shape_t *
shape_alloc_top (shape_internal_t * pr, size_t datadim,
                 size_t ptrdim)
{

  shape_t * r = shape_alloc_internal (pr, 1);

  r->m.p[0] = ushape_top (pr->man, datadim, ptrdim);

  r->msize = 1;

  r->intdim = datadim;

  r->realdim = ptrdim;

  if (pr->intdim == 0)
    {

      pr->intdim = datadim;

      pr->realdim = ptrdim;

    }

  return r;

}

shape_t *
shape_make (shape_internal_t * pr, size_t code, size_t datadim,
            size_t ptrdim)
{

  shape_t * r = shape_alloc_internal (pr, 1);

  r->m.p[0] = ushape_make (pr, code, datadim, ptrdim);

  r->msize = 1;

  r->intdim = datadim;

  r->realdim = ptrdim;

  return r;

  /* TODO: generate more than 1 element */
}

shape_t *
shape_random (shape_internal_t * pr, size_t size,

              size_t datadim, size_t ptrdim)
{

  shape_t * r = shape_alloc_internal (pr, 1);

  r->m.p[0] = ushape_random (pr, size, datadim, ptrdim);

  r->msize = 1;

  r->intdim = datadim;

  r->realdim = ptrdim;

  return r;

  /* TODO: more than 1 element */
}

inline void
shape_free_internal (shape_internal_t * pr, shape_t * a)
{

  if (!a)
    return;

  if (a->m.size > 0)

    ushape_array_clear (pr, &a->m, a->m.size);

  free (a);

}

/* TODO */
char
shape_check (shape_internal_t * pr, shape_t * a)
{

  return '.';

}

inline shape_t *
shape_copy_internal (shape_internal_t * pr, shape_t * a)
{

  shape_t * r;

  ushape_array_t * rr;

  if (!a)
    return NULL;

  checked_malloc (r, shape_t, sizeof (shape_t), 1, return NULL;
                  );

  if (a->msize > 0)
    {

      rr = ushape_array_copy (pr, &a->m, a->m.size);

      r->m = *rr;

      rr->p = NULL;

      free (rr);

    }
  else
    {

      r->m.size = 0;

      r->m.p = NULL;

    }

  r->msize = a->msize;

  r->set = a->set;

  r->closed = a->closed;

  r->intdim = a->intdim;

  r->realdim = a->realdim;

  return r;

}




/* ============================================================ */
/* Memory */

/* ============================================================ */

shape_t *
shape_copy (ap_manager_t * man, shape_t * a)
{

  shape_internal_t * pr = shape_init_from_manager (man, AP_FUNID_COPY, 0);

  return shape_copy_internal (pr, a);

}

void
shape_free (ap_manager_t * man, shape_t * a)
{

  if (!a)

    return;

  shape_internal_t * pr = shape_init_from_manager (man, AP_FUNID_FREE, 0);

  shape_free_internal (pr, a);

}

size_t
shape_size (ap_manager_t * man, shape_t * a)
{

  shape_internal_t * pr = shape_init_from_manager (man, AP_FUNID_ASIZE, 0);

  return sizeof (shape_t);

}




/* ============================================================ */
/* Control of internal representation */
/* ============================================================ */

/*
 * Return the set of minimal elements TODO: better complexity if an order of
 * ushapes is fixed
 */
void
shape_minimize (ap_manager_t * man, shape_t * a)
{

  shape_internal_t * pr =
          shape_init_from_manager (man, AP_FUNID_MINIMIZE, 0);

  ushape_array_t r;

  size_t i, j, size;

  r.size = size = 0;

  r.p = NULL;

  ushape_array_init (pr, &r, a->m.size);

  for (i = size = 0; i < a->m.size; i++)

    if (a->m.p[i])
      {

        for (j = 0; j < a->m.size; j++)

          if (i != j && a->m.p[j]
              && ushape_is_leq (pr->man, a->m.p[i], a->m.p[j]))

            break;

        if (j == a->m.size)
          { /* i is in the final set */

            r.p[size++] = ushape_copy (pr->man, a->m.p[i]);

          }

      }

  ushape_array_resize (pr, &r, size);

  ushape_array_clear (pr, &a->m, a->m.size);

  a->m = r;

  a->msize = size;

}

/* Return the set of canonical elements */
void
shape_canonicalize (ap_manager_t * man, shape_t * a)
{

  shape_internal_t * pr =

          shape_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);

  size_t i;

  ushape_array_t * r = ushape_array_make (pr, a->m.size);

  for (i = 0; i < a->m.size; i++)

    if (a->m.p[i])
      {

        ushape_array_t rc = ushape_canonicalize_internal (pr, a->m.p[i]);

        ushape_array_add_array (pr, true, r, &rc);

        ushape_array_init (pr, &rc, rc.size);

        free (rc.p);

      }

  ushape_array_clear (pr, &a->m, a->m.size);

  a->m = *r;

  a->msize = r->size;

  a->closed = a->set = true;

  r->p = NULL;

  free (r);

}

/* TODO: priority 0 */
int
shape_hash (ap_manager_t * man, shape_t * a)
{

  shape_internal_t * pr = shape_init_from_manager (man, AP_FUNID_HASH, 0);

  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,

                              "not implemented");

  return 0;

}

/* Used to signal an approximate meet in hgraphs. */
void
shape_approximate (ap_manager_t * man, shape_t * a, int algorithm)
{

  shape_internal_t * pr =
          shape_init_from_manager (man, AP_FUNID_APPROXIMATE, 0);

#ifndef NDEBUG
  fprintf (stdout, "\n****shape_approximate: algo=%d\n", algorithm);
  fflush (stdout);
#endif
  if (algorithm <= 0)
    // meet parameter
    pr->meet_algo = algorithm;
  else
    {
      // analysis parameters
      // - algorithm is split in 3 parts: (anon, segm, patterns)
      size_t i;
      int anon = (algorithm & 0xF);
      int segm = (algorithm & 0xF0) >> 4;
      pr->max_anon = anon;
      pr->segm_anon = segm;
#ifndef NDEBUG
      fprintf (stdout, "\n****shape_approximate: set max_anon=%zu, segm_anon=%zu\n",
               pr->max_anon, pr->segm_anon);
      fflush (stdout);
#endif
      if (a && pr->man_ucons < pr->size_scons)
        {
          for (i = 0; i < a->m.size; i++)
            if (a->m.p[i])
              // NULL not allowed as parameter of ap_abstract0_approximate
              ushape_approximate (man, a->m.p[i], algorithm);
          /*
          ap_manager_t* man_ucons = pr->man_scons[pr->man_ucons];
          ap_abstract0_t* b = ap_abstract0_bottom(man_ucons, 0,1);
          ap_abstract0_approximate(pr->man_scons[pr->man_ucons], b, algorithm);
          ap_abstract0_free(man_ucons,b);
           */
        }

    }
}

/* TODO: better complexity if sorted set */
bool
shape_is_minimal (ap_manager_t * man, shape_t * a)
{

  shape_internal_t * pr =

          shape_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);

  size_t i, j;

  for (i = 0; i < a->m.size; i++)

    if (a->m.p[i])

      for (j = 0; j < a->m.size; j++)

        if (i != j && a->m.p[j]
            && ushape_is_leq (pr->man, a->m.p[i], a->m.p[j]))

          return false;

  return true;

}

bool
shape_is_canonical (ap_manager_t * man, shape_t * a)
{

  shape_internal_t * pr =

          shape_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);

  size_t i;

  for (i = 0; i < a->m.size; i++)

    if (a->m.p[i] && !ushape_is_canonical (pr->man, a->m.p[i]))

      return false;

  return true;

}



/* ============================================================ */
/* Basic Constructors */

/* ============================================================ */

shape_t *
shape_bottom (ap_manager_t * man, size_t intdim,
              size_t realdim)
{

  shape_internal_t * pr = shape_init_from_manager (man, AP_FUNID_BOTTOM, 0);

  shape_t * r = shape_alloc_internal (pr, 0);

  r->intdim = intdim;

  r->realdim = realdim;

  if (pr->intdim == 0)
    {

      pr->intdim = intdim;

      pr->realdim = realdim;

    }

  return r;

}

shape_t *
shape_top (ap_manager_t * man, size_t intdim, size_t realdim)
{

  shape_internal_t * pr = shape_init_from_manager (man, AP_FUNID_TOP, 0);

  shape_t * r = shape_alloc_top (pr, intdim, realdim);

  return r;

}

/* put constraints on data variables */
shape_t *
shape_of_box (ap_manager_t * man, size_t intdim, size_t realdim,

              ap_interval_t ** t)
{

  shape_internal_t * pr = shape_init_from_manager (man, AP_FUNID_OF_BOX, 0);

  shape_t * r = shape_alloc_top (pr, intdim, realdim);

  if (t)
    {

      ushape_t * r1 = ushape_of_box (pr->man, intdim, realdim, t);

      if (!ushape_is_bottom (pr->man, r1))
        {

          ushape_free_internal (pr, r->m.p[0]);

          r->m.p[0] = r1;

        }
      else
        {

          shape_free_internal (pr, r);

          r = NULL;

        }

    }

  return r;

}

/* ============================================================ */
/* Accessors */

/* ============================================================ */

ap_dimension_t
shape_dimension (ap_manager_t * man, shape_t * a)
{

  shape_internal_t * pr =
          shape_init_from_manager (man, AP_FUNID_DIMENSION, 0);

  ap_dimension_t r;

  r.realdim = 0;

  r.intdim = 0;

  if (a)
    {

      r.intdim = a->intdim;

      r.realdim = a->realdim;

    }
  else
    {

      r.intdim = pr->intdim;

      r.realdim = pr->realdim;

    }

  return r;

}




/* ============================================================ */
/* Managers */

/* ============================================================ */

void
shape_internal_free (shape_internal_t * pr)
{

  /* TODO: free htable */
  pr->hgraphs = NULL;

  pr->pcons = NULL;

  pr->passigns = NULL;

  /* TODO: free segment managers */
  pr->man_scons = NULL;

  free (pr);

}

ap_manager_t *
shape_manager_alloc (void)
{

  size_t i;
  ap_manager_t * man;
  shape_internal_t * pr;
  char *managers;

  pr = (shape_internal_t *) malloc (sizeof (shape_internal_t));
  assert (pr);

  pr->hgraphs = NULL;
  pr->pcons = NULL;
  pr->passigns = NULL;

  /* input parameters */
  pr->intdim = 0;
  pr->realdim = 0;
  pr->meet_algo = 0;
  pr->error_ = 0;
  pr->filenum = 0;
  pr->man_mset = NULL_DIM;
  pr->man_ucons = NULL_DIM;

  /* read parameters from config file */
  FILE * f = fopen ("./cinv.txt", "r");

  managers = (char *) malloc (124 * sizeof (char));
  memset (managers, 0, 124 * sizeof (char));

  if (!f)
    {
      fprintf (stdout, "File cinv.txt not found, set default values for analysis\n");

      pr->max_anon = 0;
      pr->segm_anon = 1;

      pr->size_scons = i = 1;
      pr->man_scons = (ap_manager_t **) malloc (i * sizeof (ap_manager_t *));
      pr->man_scons[0] = lsum_manager_alloc ();

      sprintf (managers, "1.0 with anon (0,1) and abstract domains [-,LSUM,-,-]");
    }
  else
    {
      fprintf (stdout, "File cinv.txt found, read input values for analysis\n");
      /* format %zu,%zu,%s,%s,%s,%s */
      int l = fscanf (f, "%zu", &pr->max_anon);
      fgetc (f);
      l = fscanf (f, "%zu", &pr->segm_anon);
      i = 0;
      bool adomf[4] = {false, false, false, false};
      while (fgetc (f) != EOF)
        {
          char* adom = (char *) malloc (10 * sizeof (char));
          memset (adom, 0, 10 * sizeof (char));
          l = fscanf (f, "%s", adom);
          if (strcmp (adom, "LEN") == 0)
            {
              adomf[0] = true;
              i++;
            }
          else if (strcmp (adom, "LSUM") == 0)
            {
              adomf[1] = true;
              i++;
            }
          else if (strcmp (adom, "MSET") == 0)
            {
              adomf[2] = true;
              i++;
            }
          else if (strcmp (adom, "UCONS") == 0)
            {
              adomf[3] = true;
              i++;
            }
          else if (strlen (adom) != 0)
            fprintf (stdout, "*** unknown domain %s (choose in LEN|LSUM|MSET|UCONS)\n", adom);
          free (adom);
        }
      fclose (f);
      pr->size_scons = i;
      pr->man_scons =
              (i == 0) ? NULL : (ap_manager_t **) malloc (i * sizeof (ap_manager_t *));
      i = 0;
      if (adomf[0])
        pr->man_scons[i++] = oct_manager_alloc ();
      if (adomf[1])
        {
          pr->man_scons[i] = lsum_manager_alloc ();
          i++;
        }
      if (adomf[2])
        {
          pr->man_scons[i] = mset_manager_alloc ();
          mset_internal_t *pr_mset = pr->man_scons[i]->internal;
          pr_mset->segm_anon = pr->segm_anon;
          pr->man_mset = i;
          i++;
        }
      if (adomf[3])
        {
          pr->man_scons[i] = ucons_manager_alloc ();
          ucons_internal_t *pr_ucons = pr->man_scons[i]->internal;
          pr_ucons->segm_anon = pr->segm_anon;
          pr->man_ucons = i;
          i++;
        }

      sprintf (managers, "1.0 with anon (0,1) and abstract domains [%s, %s, %s, %s]",
               (adomf[0]) ? "LEN" : "-",
               (adomf[1]) ? "LSUM" : "-",
               (adomf[2]) ? "MSET" : "-",
               (adomf[3]) ? "UCONS" : "-");

    }

  man = ap_manager_alloc ("shape", managers, pr,
                          (void (*)(void *)) shape_internal_free);
  pr->man = man;


  man->funptr[AP_FUNID_COPY] = &shape_copy;
  man->funptr[AP_FUNID_FREE] = &shape_free;
  man->funptr[AP_FUNID_ASIZE] = &shape_size;
  man->funptr[AP_FUNID_MINIMIZE] = &shape_minimize;
  man->funptr[AP_FUNID_CANONICALIZE] = &shape_canonicalize;
  man->funptr[AP_FUNID_HASH] = &shape_hash;
  man->funptr[AP_FUNID_APPROXIMATE] = &shape_approximate;
  man->funptr[AP_FUNID_FPRINT] = &shape_fprint;
  man->funptr[AP_FUNID_FPRINTDIFF] = &shape_fprintdiff;
  man->funptr[AP_FUNID_FDUMP] = &shape_fdump;
  man->funptr[AP_FUNID_SERIALIZE_RAW] = &shape_serialize_raw;
  man->funptr[AP_FUNID_DESERIALIZE_RAW] = &shape_deserialize_raw;
  man->funptr[AP_FUNID_BOTTOM] = &shape_bottom;
  man->funptr[AP_FUNID_TOP] = &shape_top;
  man->funptr[AP_FUNID_OF_BOX] = &shape_of_box;
  man->funptr[AP_FUNID_DIMENSION] = &shape_dimension;
  man->funptr[AP_FUNID_IS_BOTTOM] = &shape_is_bottom;
  man->funptr[AP_FUNID_IS_TOP] = &shape_is_top;
  man->funptr[AP_FUNID_IS_LEQ] = &shape_is_leq;
  man->funptr[AP_FUNID_IS_EQ] = &shape_is_eq;
  man->funptr[AP_FUNID_IS_DIMENSION_UNCONSTRAINED] =
          &shape_is_dimension_unconstrained;
  man->funptr[AP_FUNID_SAT_INTERVAL] = &shape_sat_interval;
  man->funptr[AP_FUNID_SAT_LINCONS] = &shape_sat_lincons;
  man->funptr[AP_FUNID_SAT_TCONS] = &shape_sat_tcons;
  man->funptr[AP_FUNID_BOUND_DIMENSION] = &shape_bound_dimension;
  man->funptr[AP_FUNID_BOUND_LINEXPR] = &shape_bound_linexpr;
  man->funptr[AP_FUNID_BOUND_TEXPR] = &shape_bound_texpr;
  man->funptr[AP_FUNID_TO_BOX] = &shape_to_box;
  man->funptr[AP_FUNID_TO_LINCONS_ARRAY] = &shape_to_lincons_array;
  man->funptr[AP_FUNID_TO_TCONS_ARRAY] = &shape_to_tcons_array;
  man->funptr[AP_FUNID_TO_GENERATOR_ARRAY] = &shape_to_generator_array;
  man->funptr[AP_FUNID_MEET] = &shape_meet;
  man->funptr[AP_FUNID_MEET_ARRAY] = &shape_meet_array;
  man->funptr[AP_FUNID_MEET_LINCONS_ARRAY] = &shape_meet_lincons_array;
  man->funptr[AP_FUNID_MEET_TCONS_ARRAY] = &shape_meet_tcons_array;
  man->funptr[AP_FUNID_JOIN] = &shape_join;
  man->funptr[AP_FUNID_JOIN_ARRAY] = &shape_join_array;
  man->funptr[AP_FUNID_ADD_RAY_ARRAY] = &shape_add_ray_array;
  man->funptr[AP_FUNID_ASSIGN_LINEXPR_ARRAY] = &shape_assign_linexpr_array;
  man->funptr[AP_FUNID_SUBSTITUTE_LINEXPR_ARRAY] =
          &shape_substitute_linexpr_array;
  man->funptr[AP_FUNID_ASSIGN_TEXPR_ARRAY] = &shape_assign_texpr_array;
  man->funptr[AP_FUNID_SUBSTITUTE_TEXPR_ARRAY] =
          &shape_substitute_texpr_array;
  man->funptr[AP_FUNID_ADD_DIMENSIONS] = &shape_add_dimensions;
  man->funptr[AP_FUNID_REMOVE_DIMENSIONS] = &shape_remove_dimensions;
  man->funptr[AP_FUNID_PERMUTE_DIMENSIONS] = &shape_permute_dimensions;
  man->funptr[AP_FUNID_FORGET_ARRAY] = &shape_forget_array;
  man->funptr[AP_FUNID_EXPAND] = &shape_expand;
  man->funptr[AP_FUNID_FOLD] = &shape_fold;
  man->funptr[AP_FUNID_WIDENING] = &shape_widening;
  man->funptr[AP_FUNID_CLOSURE] = &shape_closure;

  for (i = 0; i < AP_EXC_SIZE; i++)
    ap_manager_set_abort_if_exception (man, i, false);

  return man;
}

shape_t *
shape_of_abstract0 (ap_abstract0_t * a)
{

  return (shape_t *) a->value;

}

ap_abstract0_t *
abstract0_of_shape (ap_manager_t * man, shape_t * a)
{

  ap_abstract0_t * r = malloc (sizeof (ap_abstract0_t));

  assert (r);

  r->value = a;

  r->man = ap_manager_copy (man);

  return r;

}


