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
#include "ushape.h"
#include "ushape_internal.h"
#include "hgraph_internal.h"
#include "shape_manager.h"
#include "shape_options.h"
#include "shape_macros.h"
#include "ap_abstract0.h"
#include "box.h"
#include "oct.h"
#include "pk.h"
//#include "ucons.h"
#include "lsum.h"



/* ********************************************************************** */
/* UShapes */
/* ********************************************************************** */

/* ============================================================ */
/* Internal Management */
/* ============================================================ */

/* returns an empty ushape, i.e., bottom ushape */
inline ushape_t *
ushape_alloc_internal (ushape_internal_t * pr, size_t intdim, size_t realdim)
{
  ushape_t *r;
  checked_malloc (r, ushape_t, sizeof (ushape_t), 1, return NULL;
    );
  r->h = r->closed = NULL;
  r->datadim = intdim;
  r->ptrdim = realdim;

  if (pr->size_scons > 0)
    {
      size_t i, j, datadim_y;
      r->scons = (void **) malloc (sizeof (void *) * pr->size_scons);
      memset (r->scons, 0, pr->size_scons * sizeof (void *));
      /*
         for (i = 0; i < pr->size_scons; i++)
         r->scons[i] = NULL;
       */
    }
  else
    r->scons = NULL;
  return r;
}

ushape_t *
ushape_make (ushape_internal_t * pr, size_t code,
	     size_t datadim, size_t ptrdim)
{
  ushape_t *a;
  hgraph_t *h;
  size_t i;

  /* build the hgraph */
  switch (code)
    {
    case 0:			/* x --> _null and l[x]==_l and l[x]>=1 and data(x) */
      {
	h = hgraph_make (pr, 0, datadim, ptrdim);
	break;
      }
    case 1:			/* x-->y-->null and l[x]+l[y]==_l>=1 and l[x]>=1 and l[y]>=1 and data(xy) */
      {
	h = hgraph_make (pr, 1, datadim, ptrdim);
	break;
      }
    case 2:			/* x -->_null and y -->_null and l[x]==l[y]==_l>=1 and data(x) */
    case 3:			/* x --> _null and y --> _null and l[x]==_l>=1 and l[y]+1<=_l and data(x) and data(y) */
      {
	h = hgraph_make (pr, 2, datadim, ptrdim);
	break;
      }
    case 4:			/* x-->null and y-->null and z-->null and disjoint(x,y,z) and l[x]=l[y]=l[z] and data(x) and data(y) and data(z) */
      {
	h = hgraph_make (pr, 3, datadim, ptrdim);
	break;
      }
    default:			/* TODO: add other cases for ushapes */
      h = hgraph_top (pr->man, datadim, ptrdim);
    }

  a = ushape_alloc_internal (pr, datadim, ptrdim);
  a->h = h;

  size_t raw[3] = { datadim, h->size, code };
  for (i = 0; i < pr->size_scons; i++)
    a->scons[i] = ap_abstract0_deserialize_raw (pr->man_scons[i], raw, NULL);
  ushape_check (pr, a);
#ifndef NDEBUG
  ushape_fdump (stdout, pr->man, a);
#endif
  return a;
}

ushape_t *
ushape_random (ushape_internal_t * pr, size_t size,
	       size_t datadim, size_t ptrdim)
{
  ushape_t *r = ushape_alloc_internal (pr, datadim, ptrdim);
  /* generate a random graph with size nodes */
  r->h = hgraph_random (pr, size, datadim, ptrdim);
  /*
   * generate a random lincons on prg data variables and data fields of ptr
   * variables
   */
  ap_lincons0_t c = ap_lincons0_make ((lrand48 () % 100 >= 80) ? AP_CONS_EQ :
				      (lrand48 () % 100 >= 80) ? AP_CONS_SUP :
				      AP_CONS_SUPEQ,
				      shape_linexpr_random (pr,
							    (lrand48 () %
							     100 >=
							     80) ?
							    expr_lindata :
							    expr_data,
							    datadim, ptrdim),
				      NULL);
  size_t *v2n = hgraph_get_var2node (r->h);
  ap_lincons0_t dc = // set kind (scalar) to data constraint
    shape_lincons_of_node (pr, &c, NULL, v2n, r->h->size, datadim, ptrdim);
  free (v2n);
  ap_lincons0_array_t arr;
  size_t i;
  arr.size = 1;
  arr.p = &dc;
  for (i = 0; i < pr->size_scons; i++)
    r->scons[i] = ap_abstract0_meet_lincons_array (pr->man_scons[i], true,
						   ap_abstract0_top (pr->
								     man_scons
								     [i],
								     datadim,
								     r->h->
								     size),
						   &arr);
  ap_lincons0_clear (&c);
  ap_lincons0_clear (&dc);
  ushape_check (pr, r);
  return r;
}

void
ushape_set_bottom (ushape_internal_t * pr, ushape_t * a)
{
  if (a->h)
    {
      hgraph_free_internal ((ushape_internal_t *) pr, a->h);
      a->h = NULL;
    }
  if (a->closed)
    {
      hgraph_free_internal ((hgraph_internal_t *) pr, a->closed);
      a->closed = NULL;
    }
  if (a->scons)
    {
      size_t i;
      for (i = 0; i < pr->size_scons; i++)
	if (a->scons[i])
	  {
	    ap_abstract0_free (pr->man_scons[i], a->scons[i]);
	    a->scons[i] = NULL;
	  }
      free (a->scons);
      a->scons = NULL;
    }
}

void
ushape_free_internal (ushape_internal_t * pr, ushape_t * a)
{
  if (a)
    {
      ushape_set_bottom (pr, a);
      free (a);
    }
}


char
ushape_check (ushape_internal_t * pr, ushape_t * a)
{
  size_t i, size;
  if (!a)
    {
      ERROR ("ushape null!",;
	);
      return '!';
    }
  size = (a->h) ? a->h->size : 0;
  /* constraints have at least dimension h->size */
  hgraph_check (pr, a->h);
  for (i = 0; i < pr->size_scons; i++)
    {
      ap_dimension_t dc =
	ap_abstract0_dimension (pr->man_scons[i], a->scons[i]);
      if (dc.realdim < size || dc.intdim != a->datadim)
	{
	  ERROR ("ushape segment constraints dimension!",;
	    );
	  return '!';
	}
    }
  /* TODO: other checks? */
  return '.';
}

inline ushape_t *
ushape_copy_internal (ushape_internal_t * pr, ushape_t * a)
{
  ushape_t *r = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  size_t i;
  r->h = hgraph_copy_internal (pr, a->h);
  r->closed = hgraph_copy_internal (pr, a->closed);
  r->datadim = a->datadim;
  r->ptrdim = a->ptrdim;
  ushape_copy_internal_scons(pr, a, r);
  return r;
}

inline void
ushape_copy_internal_scons (ushape_internal_t * pr, ushape_t * a, ushape_t * b)
{
  if (a && b && a->scons && b->scons)
    {
      size_t i;
      for (i = 0; i < pr->size_scons; i++)
	if (a->scons[i])
	  b->scons[i] = ap_abstract0_copy (pr->man_scons[i], a->scons[i]);
    }
}

inline size_t
ushape_get_size (ushape_t * a)
{
  return a->h->size;
}

/*
 * TODO: used? If destructive, returns a with fields h and closed updated and
 * frees the older fields. Otherwise, returns a ushape with the same fields
 * as a and only h and closed changed.
 */
ushape_t *
ushape_set_hgraph (ushape_internal_t * pr, ushape_t * a, hgraph_t * h,
		   hgraph_t * closed, bool destructive)
{
  ushape_t *r;
  if (destructive)
    {
      /* free non-aliased ushapes */
      if (a->h && a->h != h && a->h != closed)
	hgraph_free (pr->man, a->h);
      if (a->closed && a->closed != h && a->closed != closed)
	hgraph_free (pr->man, a->closed);
      r = a;
    }
  else
    {
      /* copy aliased ushapes */
      r = ushape_copy_internal (pr, a);
      if (h && (a->h == h || a->closed == h))
	h = hgraph_copy (pr->man, h);
      if (closed && (a->h == closed || a->closed == closed))
	closed = hgraph_copy (pr->man, closed);
    }
  r->h = h;
  r->closed = closed;
  return r;
}

/* ============================================================ */
/* Memory */
/* ============================================================ */

ushape_t *
ushape_copy (ap_manager_t * man, ushape_t * a)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_COPY, 0);
  return ushape_copy_internal (pr, a);
}

void
ushape_free (ap_manager_t * man, ushape_t * a)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_FREE, 0);
  ushape_free_internal (pr, a);
}

inline size_t
ushape_size (ap_manager_t * man, ushape_t * a)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_SIZE, 0);
  if (!a)
    return 0;
  return sizeof (ushape_t);
}


/* ============================================================ */
/* Control of internal representation */
/* ============================================================ */

/* TODO: priority 3 */
void
ushape_minimize (ap_manager_t * man, ushape_t * a)
{
  ushape_internal_t *pr =
    ushape_init_from_manager (man, AP_FUNID_MINIMIZE, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			      "not implemented");
}

ushape_array_t
ushape_canonicalize_internal (ushape_internal_t * pr, ushape_t * a)
{
  ushape_array_t r;
  ushape_array_init (pr, &r, 1);
  return r;
}

/* Put to NODE_NULL all undefined ptr vars. */
void
ushape_canonicalize (ap_manager_t * man, ushape_t * a)
{
  ushape_internal_t *pr =
    ushape_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);
  hgraph_t* h = hgraph_copy_mem(pr,a->h);
  hgraph_canonicalize(man,h); // SIDE EFFECT!!!!
  //hgraph_free_internal(pr,a->h);
  a->h = hgraph_copy_internal(pr,h);
  hgraph_free_internal(pr,h);
}

/* TODO: priority 0 */
int
ushape_hash (ap_manager_t * man, ushape_t * a)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_HASH, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			      "not implemented");
  return 0;
}

/* Used to change the parameters of the analysis */
void
ushape_approximate (ap_manager_t * man, ushape_t * a, int algorithm)
{
  ushape_internal_t *pr =
    ushape_init_from_manager (man, AP_FUNID_APPROXIMATE, 0);
  // anonymous already changed, change only constraints 
  if (a && a->scons)
    {
      size_t i;
      for (i = 0; i < pr->size_scons; i++)
	if (a->scons[i])
	  // NULL not accepted as input in ap_abstract0_approximate
	  ap_abstract0_approximate (pr->man_scons[i], a->scons[i], algorithm);
    }
}

/* TODO: priority 3 */
bool
ushape_is_minimal (ap_manager_t * man, ushape_t * a)
{
  ushape_internal_t *pr =
    ushape_init_from_manager (man, AP_FUNID_MINIMIZE, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
			      "not implemented");
  return true;
}

bool
ushape_is_canonical (ap_manager_t * man, ushape_t * a)
{
  ushape_internal_t *pr =
    ushape_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);
  if (!a)
    return true;
  return hgraph_is_canonical (man, a->h);
}


/* ============================================================ */
/* Basic Constructors */
/* ============================================================ */

ushape_t *
ushape_bottom (ap_manager_t * man, size_t intdim, size_t realdim)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_BOTTOM, 0);
  ushape_t *r = ushape_alloc_internal (pr, intdim, realdim);
  return r;
}

ushape_t *
ushape_top (ap_manager_t * man, size_t intdim, size_t realdim)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_TOP, 0);
  ushape_t *r = ushape_alloc_internal (pr, intdim, realdim);
  size_t i;
  r->h = hgraph_top (pr->man, intdim, realdim);
  r->closed = NULL;
  r->datadim = intdim;
  r->ptrdim = realdim;
  /* the segment constraints are set to top */
  for (i = 0; i < pr->size_scons; i++)
    r->scons[i] = ap_abstract0_top (pr->man_scons[i], intdim, r->h->size);
  return r;
}

/* Put only constraints on data vars */
ushape_t *
ushape_of_box (ap_manager_t * man, size_t intdim, size_t realdim,
	       ap_interval_t ** t)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_OF_BOX, 0);
  ushape_t *r = ushape_alloc_internal (pr, intdim, realdim);
  size_t i;
  for (i = 0; i < pr->size_scons; i++)
    r->scons[i] = ap_abstract0_of_box (pr->man_scons[i], intdim + 1, 0, t);
  return r;
}


/* ============================================================ */
/* Accessors */
/* ============================================================ */

ap_dimension_t
ushape_dimension_var (ap_manager_t * man, ushape_t * a)
{
  ushape_internal_t *pr =
    ushape_init_from_manager (man, AP_FUNID_DIMENSION, 0);
  ap_dimension_t r;
  r.intdim = 0;
  r.realdim = 0;
  if (a)
    {
      r.intdim = a->datadim;
      r.realdim = a->ptrdim;
    }
  return r;
}

ap_dimension_t
ushape_dimension (ap_manager_t * man, ushape_t * a)
{
  ushape_internal_t *pr =
    ushape_init_from_manager (man, AP_FUNID_DIMENSION, 0);
  ap_dimension_t r;
  r.intdim = 0;
  r.realdim = 0;
  if (a)
    {
      /*
         r.intdim = a->datadim + ((a->h) ? 2*a->h->size : 0);
         r.realdim = a->ptrdim;
       */
      r.intdim = a->datadim;
      r.realdim = (a->h) ? a->h->size : 0;
    }
  return r;
}


/* ============================================================ */
/* Manager */
/* ============================================================ */

void
ushape_internal_free (ushape_internal_t * pr)
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
ushape_manager_alloc (void)
{
  size_t i;
  ap_manager_t *man;
  ushape_internal_t *pr;
  char domain[128];

  pr = (ushape_internal_t *) malloc (sizeof (ushape_internal_t));
  assert (pr);
  pr->hgraphs = NULL;
  pr->pcons = NULL;
  pr->passigns = NULL;

  i = 0;
  if (SHAPE_SCONS_LEN)
    i++;
  if (SHAPE_SCONS_UCONS)
    i++;
  if (SHAPE_SCONS_MSET)
    i++;
  pr->size_scons = i;
  pr->man_scons =
    (i == 0) ? NULL : (ap_manager_t **) malloc (i * sizeof (ap_manager_t *));

  i = 0;
  /*
     if (SHAPE_SCONS_LEN)
     {
     switch (SHAPE_SCONS_LEN_DOMAIN)
     {
     case DOM_BOX:
     pr->man_scons[i] = box_manager_alloc ();
     break;
     case DOM_OCT:
     pr->man_scons[i] = oct_manager_alloc ();
     break;
     default:
     pr->man_scons[i] = pk_manager_alloc (false);
     break;
     }
     i++;
     }
   */
  if (SHAPE_SCONS_UCONS)
    pr->man_scons[i++] = NULL;	/* ucons_manager_alloc (); */
  if (SHAPE_SCONS_MSET)
    pr->man_scons[i++] = NULL;
  if (SHAPE_SCONS_LSUM)
    pr->man_scons[i] = lsum_manager_alloc ();

  pr->meet_algo = 0;
  pr->max_anon = 4;
  pr->segm_anon = 1;
  pr->error_ = 0;

  i = snprintf (domain, 127, "0.1 with dcons=%s and scons = [%s,%s,%s]",
		((SHAPE_DCONS_DOMAIN ==
		  DOM_BOX) ? "Box" : ((SHAPE_DCONS_DOMAIN ==
				       DOM_OCT) ? "Oct" : "Polka")),
		((SHAPE_SCONS_LEN) ? "LEN" : "null"),
		((SHAPE_SCONS_UCONS) ? "UCONS" : "null"),
		((SHAPE_SCONS_MSET) ? "MSET" : "null"));
  domain[i] = '\0';

  man = ap_manager_alloc ("ushape", domain, pr,
			  (void (*)(void *)) ushape_internal_free);

  pr->man = man;

  man->funptr[AP_FUNID_COPY] = &ushape_copy;
  man->funptr[AP_FUNID_FREE] = &ushape_free;
  man->funptr[AP_FUNID_ASIZE] = &ushape_size;
  man->funptr[AP_FUNID_MINIMIZE] = &ushape_minimize;
  man->funptr[AP_FUNID_CANONICALIZE] = &ushape_canonicalize;
  man->funptr[AP_FUNID_HASH] = &ushape_hash;
  man->funptr[AP_FUNID_APPROXIMATE] = &ushape_approximate;
  man->funptr[AP_FUNID_FPRINT] = &ushape_fprint;
  man->funptr[AP_FUNID_FPRINTDIFF] = &ushape_fprintdiff;
  man->funptr[AP_FUNID_FDUMP] = &ushape_fdump;
  man->funptr[AP_FUNID_SERIALIZE_RAW] = &ushape_serialize_raw;
  man->funptr[AP_FUNID_DESERIALIZE_RAW] = &ushape_deserialize_raw;
  man->funptr[AP_FUNID_BOTTOM] = &ushape_bottom;
  man->funptr[AP_FUNID_TOP] = &ushape_top;
  man->funptr[AP_FUNID_OF_BOX] = &ushape_of_box;
  man->funptr[AP_FUNID_DIMENSION] = &ushape_dimension;
  man->funptr[AP_FUNID_IS_BOTTOM] = &ushape_is_bottom;
  man->funptr[AP_FUNID_IS_TOP] = &ushape_is_top;
  man->funptr[AP_FUNID_IS_LEQ] = &ushape_is_leq;
  man->funptr[AP_FUNID_IS_EQ] = &ushape_is_eq;
  man->funptr[AP_FUNID_IS_DIMENSION_UNCONSTRAINED] =
    &ushape_is_dimension_unconstrained;
  man->funptr[AP_FUNID_SAT_INTERVAL] = &ushape_sat_interval;
  man->funptr[AP_FUNID_SAT_LINCONS] = &ushape_sat_lincons;
  man->funptr[AP_FUNID_SAT_TCONS] = &ushape_sat_tcons;
  man->funptr[AP_FUNID_BOUND_DIMENSION] = &ushape_bound_dimension;
  man->funptr[AP_FUNID_BOUND_LINEXPR] = &ushape_bound_linexpr;
  man->funptr[AP_FUNID_BOUND_TEXPR] = &ushape_bound_texpr;
  man->funptr[AP_FUNID_TO_BOX] = &ushape_to_box;
  man->funptr[AP_FUNID_TO_LINCONS_ARRAY] = &ushape_to_lincons_array;
  man->funptr[AP_FUNID_TO_TCONS_ARRAY] = &ushape_to_tcons_array;
  man->funptr[AP_FUNID_TO_GENERATOR_ARRAY] = &ushape_to_generator_array;
  man->funptr[AP_FUNID_MEET] = &ushape_meet;
  man->funptr[AP_FUNID_MEET_ARRAY] = &ushape_meet_array;
  man->funptr[AP_FUNID_MEET_LINCONS_ARRAY] = &ushape_meet_lincons_array;
  man->funptr[AP_FUNID_MEET_TCONS_ARRAY] = &ushape_meet_tcons_array;
  man->funptr[AP_FUNID_JOIN] = &ushape_join;
  man->funptr[AP_FUNID_JOIN_ARRAY] = &ushape_join_array;
  man->funptr[AP_FUNID_ADD_RAY_ARRAY] = &ushape_add_ray_array;
  man->funptr[AP_FUNID_ASSIGN_LINEXPR_ARRAY] = &ushape_assign_linexpr_array;
  man->funptr[AP_FUNID_SUBSTITUTE_LINEXPR_ARRAY] =
    &ushape_substitute_linexpr_array;
  man->funptr[AP_FUNID_ASSIGN_TEXPR_ARRAY] = &ushape_assign_texpr_array;
  man->funptr[AP_FUNID_SUBSTITUTE_TEXPR_ARRAY] =
    &ushape_substitute_texpr_array;
  man->funptr[AP_FUNID_ADD_DIMENSIONS] = &ushape_add_dimensions;
  man->funptr[AP_FUNID_REMOVE_DIMENSIONS] = &ushape_remove_dimensions;
  man->funptr[AP_FUNID_PERMUTE_DIMENSIONS] = &ushape_permute_dimensions;
  man->funptr[AP_FUNID_FORGET_ARRAY] = &ushape_forget_array;
  man->funptr[AP_FUNID_EXPAND] = &ushape_expand;
  man->funptr[AP_FUNID_FOLD] = &ushape_fold;
  man->funptr[AP_FUNID_WIDENING] = &ushape_widening;
  man->funptr[AP_FUNID_CLOSURE] = &ushape_closure;

  for (i = 0; i < AP_EXC_SIZE; i++)
    ap_manager_set_abort_if_exception (man, i, false);

  return man;
}

ushape_t *
ushape_of_abstract0 (ap_abstract0_t * a)
{
  return (ushape_t *) a->value;
}

ap_abstract0_t *
abstract0_of_ushape (ap_manager_t * man, ushape_t * a)
{
  ap_abstract0_t *r = malloc (sizeof (ap_abstract0_t));
  assert (r);
  r->value = a;
  r->man = ap_manager_copy (man);
  return r;
}

/* ********************************************************************** */
/* II. ushape_array_t */
/* ********************************************************************** */

void
ushape_array_init (ushape_internal_t * pr, ushape_array_t * a, size_t size)
{
  size_t i;
  arg_assert (a && (!a->p || (a->p && a->size >= size)), return;
    );
  if (!a->p)
    {
      checked_malloc (a->p, ushape_t *, sizeof (ushape_t *), size, return;
	);
      a->size = size;
    }
  for (i = 0; i < size; i++)
    a->p[i] = NULL;
}

ushape_array_t *
ushape_array_make (ushape_internal_t * pr, size_t size)
{
  ushape_array_t *r;
  size_t i;
  checked_malloc (r, ushape_array_t, sizeof (ushape_array_t), 1, return NULL;
    );
  r->p = NULL;
  r->size = 0;
  ushape_array_init (pr, r, size);
  return r;
}

/* Add a ushape to an array either keeping the set or not */
int
ushape_array_add (ushape_internal_t * pr, bool isset, ushape_array_t * arr,
		  size_t msize, bool docopy, bool destructive, ushape_t * a)
{
  size_t found;
  size_t i, j;
  ushape_t *u = NULL;
  int r;
  if (!a)
    return 0;
  if (!arr)
    {
      if (destructive)
	ushape_free_internal (pr, a);
      return 0;
    }
  found = arr->size;
  j = arr->size;		/* first NULL position */
  /* search if the element is already in the set and also the first position NULL */
  for (i = 0; i < arr->size && found == arr->size; i++)
    {
      if (arr->p[i] == NULL && j == arr->size)
	j = i;
      if (arr->p[i] && isset)
	{
	  u = ushape_join (pr->man, false, a, arr->p[i]);
	  if (u)
	    found = i;
	}
    }
  if (found == arr->size)
    {
      /* now, add the ushape at the position msize, if given */
      if (msize == 0)
	msize = j;
      if (arr->size == msize)
	ushape_array_resize (pr, arr, arr->size + 4);

      arr->p[msize] = (docopy) ? ushape_copy_internal (pr, a) : a;
      r = 1;
    }
  else
    {
      // put u at position found
      arr->p[found] = u;
      r = 0;
    }
  if (destructive)
    ushape_free_internal (pr, a);
  return r;
}

ushape_array_t *
ushape_array_add_array (ushape_internal_t * pr, bool isset,
			ushape_array_t * a, ushape_array_t * b)
{
  ushape_array_t *r;
  if ((!a || !a->size) && (!b || !b->size))
    r = NULL;
  else if (!a && b && b->size)
    r = ushape_array_copy (pr, b, b->size);
  else if (a && !b && a->size)
    r = a;
  else
    {
      size_t i, j;
      for (i = 0; i < a->size; i++)	/* TODO: start with msize? */
	if (!a->p[i])
	  break;
      j = i;
      for (i = 0; i < b->size; i++)
	j += ushape_array_add (pr, isset, a, j, false, false, b->p[i]);
      r = a;
    }
  return r;
}

int
ushape_array_resize (ushape_internal_t * pr, ushape_array_t * a, size_t size)
{
  size_t i;
  if (!a)
    return 0;
  for (i = size; i < a->size; i++)
    {
      ushape_free_internal (pr, a->p[i]);
      a->p[i] = NULL;
    }
  a->p = (ushape_t **) realloc (a->p, size * sizeof (ushape_t *));
  for (i = a->size; i < size; i++)
    {
      a->p[i] = NULL;
    }
  a->size = size;
  return size;
}

void
ushape_array_clear (ushape_internal_t * pr, ushape_array_t * a, size_t size)
{
  size_t i;
  for (i = 0; i < size; i++)
    {
      if (a->p[i])
	ushape_free_internal (pr, a->p[i]);
      a->p[i] = NULL;
    }
  free (a->p);
}

/*
 * Copy the src array into (already allocated) array dst for length size
 */
ushape_array_t *
ushape_array_copy (ushape_internal_t * pr, ushape_array_t * src, size_t size)
{
  size_t i;
  arg_assert (src && 0 < size, return NULL;
    );
  ushape_array_t *r = ushape_array_make (pr, size);
  if (size > 0)
    {
      for (i = 0; i < src->size && i < size; i++)
	if (src->p[i])
	  r->p[i] = ushape_copy_internal (pr, src->p[i]);
	else
	  r->p[i] = NULL;
      for (; i < size; i++)
	r->p[i] = NULL;
    }
  return r;
}
