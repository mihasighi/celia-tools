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


#include "hgraph.h"
#include "hgraph_internal.h"
#include "apron2shape.h"
#include "shape_macros.h"
#include "ap_generic.h"



/* ============================================================ */
/* Tests */
/* ============================================================ */

/*
 * Only needed if hgraphs are used as an independent domain In this case,
 * bottom value is the graph with no nodes
 */
bool
hgraph_is_bottom (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_IS_BOTTOM, 0);
  return (!a || a->size == 0) ? true : false;
}

/*
 * Only needed if hgraphs are used as an independent domain In this case, top
 * value value means that only one node exists (null) and all variables are
 * mapped ont it?
 */
bool
hgraph_is_top (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_IS_TOP, 0);
  size_t i;

  if (!a || (a->size == 0))
    return false;
  for (i = NULL_DIM + 1; i < a->ptrdim; i++)
    if (VAR2NODE (a, i) < a->size)
      return false;
  return true;
}

bool
hgraph_is_equal (hgraph_t * a, hgraph_t * b)
{
  /*
   * same size, dim, succ, n2v and v2n, here we assume the normalized
   * representation of hgraphs
   */
  bool r = (a->size == b->size) && (a->ptrdim == b->ptrdim) &&
          (a->datadim == b->datadim);
  r = r &&
          (memcmp (a->info, b->info, (a->ptrdim + a->size) * sizeof (node_info_t))
           == 0);
  return r;
}

bool
hgraph_is_lt (hgraph_t * a, hgraph_t * b)
{
  node_t i;
  int na, nb;
  bool r;

  /*
   * same nb of ptr vars, same cut nodes, same labeling of nodes, but
   * some anonymous nodes of b are not in a
   */
  /*
   * since the same nb of cut nodes, the field v for cut nodes is the
   * same by canonical representation TODO: explain;
   */
  if (a->ptrdim != b->ptrdim ||
      a->size >= b->size ||
      memcmp (a->info, b->info, a->ptrdim * sizeof (node_info_t)) != 0)
    return false;

  /* verify successor relation transitively by non-cut nodes */
  r = true;
  for (i = 0; i < a->size; i++)
    {
      if (hgraph_node_is_cut (a, i) == false)
        break; /* end of cut nodes */
      if (hgraph_node_get_succ_cut (a, i, &na) !=
          hgraph_node_get_succ_cut (b, i, &nb))
        return false;
      if (na > nb)
        return false;
    }
  return true;
}

/* Comparison (like in C) of two hgraphs */
int
hgraph_cmp (hgraph_t * a, hgraph_t * b)
{
  /* if equality return 0 */
  if (hgraph_is_equal (a, b) == true)
    return 0;

  /* graphs are not comparable (return 2) if ptr vars are not the same */
  if (a->ptrdim != b->ptrdim)
    return 2;

  /* returns -1 and 1 if a included in b resp. b included in a */
  if (hgraph_is_lt (a, b) == true)
    return -1;
  if (hgraph_is_lt (b, a) == true)
    return 1;

  return 2;
}

bool
hgraph_is_leq (ap_manager_t * man, hgraph_t * a1, hgraph_t * a2)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_IS_LEQ, 0);
  if (a1 == NULL || (a2 == NULL && hgraph_is_bottom (man, a1)))
    return true;
  if (a2 == NULL)
    return false;
  arg_assert (a1->ptrdim == a2->ptrdim, return false;
              );
  if (hgraph_is_bottom (man, a1) || hgraph_is_top (man, a2))
    return true;
  if ((hgraph_is_bottom (man, a2) && !hgraph_is_bottom (man, a1)) ||
      (hgraph_is_top (man, a1) && !hgraph_is_top (man, a2)))
    return false;
  return hgraph_is_equal (a1, a2) || hgraph_is_lt (a1, a2);
}

bool
hgraph_is_eq (ap_manager_t * man, hgraph_t * a1, hgraph_t * a2)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_IS_EQ, 0);

  if ((a1 == NULL || hgraph_is_bottom (man, a1)) &&
      (a2 == NULL || hgraph_is_bottom (man, a2)))
    return true;
  if (!a1 || !a2 || (a1->ptrdim != a2->ptrdim))
    return false;
#ifndef NDEBUG
  hgraph_fdump (stdout, man, a1);
  hgraph_fdump (stdout, man, a2);
#endif
  return hgraph_is_equal (a1, a2);
}

/**
 * Return true if all ptr variables non mapped to NODE_T_TOP in a2
 * form a graph isomorphic with a subgraph of a1. 
 * The permutation of nodes in a2 to obtain nodes in a1 is
 * computed in perm2 (of size a2->size).
 */
bool
hgraph_is_spec (ap_manager_t * man,
                hgraph_t * a1, hgraph_t * a2, ap_dimperm_t* perm2)
{
  hgraph_t* r = NULL;
#ifndef NDEBUG1
  fprintf (stdout, "!!!! hgraph_is_spec: a1=(");
  hgraph_fdump (stdout, man, a1);
  fprintf (stdout, ") and a2=(");
  hgraph_fdump (stdout, man, a2);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif
  // Step 1: go from variables not mapped to NODE_T_TOP in a2 to
  // build the mapping to nodes in a1
  bool ok = true;
  size_t v, n1, n2;
  size_t *map12 = (size_t*) malloc (a1->size * sizeof (size_t));
  size_t *map21 = (size_t*) malloc (a2->size * sizeof (size_t));
  memset (map12, 0, a1->size * sizeof (size_t));
  memset (map21, 0, a2->size * sizeof (size_t));
  for (v = 0; (v < a1->ptrdim) && ok; v++)
    {
      if (VAR2NODE (a2, v) == NODE_T_TOP) continue;
      // else
#ifndef NDEBUG1
      fprintf (stdout, "!!!! hgraph_is_spec: vars x%zu\n", v);
      fflush (stdout);
#endif
      bool isEnd = false;
      do
        {
          n1 = VAR2NODE (a1, v);
          n2 = VAR2NODE (a2, v);
          if (n1 == NODE_NULL || n2 == NODE_NULL)
            {
              ok = (n1 == n2);
              isEnd = true;
            }
          else if (map21[n2] != 0)
            {
              ok = (map21[n2] == n1);
              isEnd = true; /* already seen or error */
            }
          else if (map12[n1] != 0)
            {
              ok = (map12[n1] == n2);
              isEnd = true; /* already seen  or error */
            }
          else
            {
              map21[n2] = n1;
              map12[n1] = n2;
#ifndef NDEBUG1
              fprintf (stdout, "!!!! hgraph_is_spec: maps n%zu (a2) to n%zu (a1)\n", n2, n1);
              fflush (stdout);
#endif
              n1 = NODE_NEXT (a1, n1);
              n2 = NODE_NEXT (a2, n2);
            }
        }
      while (ok && !isEnd);

    }
#ifndef NDEBUG1
  fprintf (stdout, "!!!! hgraph_is_spec: maps12=(");
  for (n1 = 0; n1 < a1->size; n1++)
    fprintf (stdout, "%zu -> %zu, ", n1, map12[n1]);
  fprintf (stdout, ")\n and maps21=(");
  for (n2 = 0; n2 < a2->size; n2++)
    fprintf (stdout, "%zu -> %zu, ", n2, map21[n2]);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  if (ok)
    { // there is a mapping from a2 to a subgraph of a1
      // Step 2: build perm2
      ap_dimperm_init (perm2, a1->size);
      ap_dimperm_set_id (perm2);
      for (n2 = 0; n2 < a2->size; n2++)
        perm2->dim[n2] = map21[n2];
      // the added nodes to a2 are mapped (in order) to nodes in a1
      n1 = 1;
      for (; n2 < a1->size; n2++)
        {
          // find next node not mapped in a2
          while ((n1 < a1->size) && (map12[n1] != 0)) n1++;
          perm2->dim[n2] = n1;
	  map12[n1] = n2;
          n1++;
        }
#ifndef NDEBUG1
      fprintf (stdout, "!!!! hgraph_is_spec: generated perm2=(");
      ap_dimperm_fprint (stdout, perm2);
      fprintf (stdout, ")\n");
#endif
    }

  // Step 3: free allocated memory
  free (map12);
  free (map21);
  return ok;
}

bool
hgraph_sat_pcons (hgraph_internal_t * pr, hgraph_t * a, pcons0_t * c)
{
  bool r = false;
  arg_assert (a && c, return false;
              );
  if (c->type == DATA_CONS)
    return false;
  node_t nx = VAR2NODE (a, DIM2PTR (c->info.ptr.x, a->datadim));
  node_t ny = VAR2NODE (a, DIM2PTR (c->info.ptr.y, a->datadim));
  switch (c->type)
    {
    case EQ_CONS:
      r = (nx == ny) ? true : false;
      break;
    case NE_CONS:
      r = (nx != ny) ? true : false;
      break;
    case REACH_CONS:
      r = hgraph_node_is_reachable (a, nx, ny);

      break;
    default:
      r = false;
    }
  return r;
}

/*
 * Returns true iff all anonymous paths have length less equal than pr->max_anon
 * and if yes, the number of anonymous nodes is a multiple of pr->segm_anon
 */
bool
hgraph_is_closed (hgraph_internal_t * pr, hgraph_t * a)
{
  bool r = true;
  size_t i, n;
#ifndef NDEBUG
  fprintf (stdout, "\n!!!!hgraph_is_closed: (%d,%d)\n",
           pr->max_anon, pr->segm_anon);
#endif
  if (a)
    {
      for (i = 1, n = 0; i < a->size; i++)
        if (NODE_VAR (a, i) != NODE_T_TOP && // not a garbage
            NODE_VAR_NEXT (a, i) != NODE_T_TOP &&
            NODE_VAR_NEXT (a, i) > pr->max_anon)
          {
            r = false;
            n++;
          }

      if (!r && pr->segm_anon && (n % pr->segm_anon))
        r = true;
    }
  return r;
}

/*
 * Possibly needed for checking linear constraints representing aliasing
 * between pointer variables, e.g., x = y, in the assert statements.
 */
bool
hgraph_sat_lincons (ap_manager_t * man, hgraph_t * a, ap_lincons0_t * lincons)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_SAT_LINCONS, 0);
  if (hgraph_is_bottom (man, a))
    return false;
  pcons0_t * pcons =
          shape_pcons_of_lincons (pr, lincons, a->datadim, a->ptrdim);

  return hgraph_sat_pcons (pr, a, pcons);
  /* pcons are hashed, so no free of pcons here */
}

/*
 * Possibly needed for checking constraints representing aliasing between
 * pointer variables, e.g., x*next = y, in the assert statements.
 */
bool
hgraph_sat_tcons (ap_manager_t * man, hgraph_t * a, ap_tcons0_t * cons)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_SAT_TCONS, 0);
  if (hgraph_is_bottom (man, a))
    return false;
  pcons0_t * pcons = shape_pcons_of_tcons (pr, cons, a->datadim, a->ptrdim);

  return hgraph_sat_pcons (pr, a, pcons);
  /* pcons are hashed, so no free of pcons here */
}

/* NOT IMPLEMENTED */
bool
hgraph_sat_interval (ap_manager_t * man, hgraph_t * a,
                     ap_dim_t dim, ap_interval_t * i)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_SAT_INTERVAL, 0);
  arg_assert (dim < a->ptrdim, return false;
              );

  return true;
}

bool
hgraph_is_dimension_unconstrained (ap_manager_t * man, hgraph_t * a,
                                   ap_dim_t dim)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_IS_DIMENSION_UNCONSTRAINED, 0);
  size_t pdim = REAL2PTR_DIM (pr, dim);
  arg_assert (pdim < a->ptrdim, return false;
              );
  /* bottom constraints all dimensions */
  if (a->size == 0)
    return false;
  /* top does not constraint dimensions */
  if (a->size == 1)
    return true;
  /*
   * for any other graph, true iff the ptr vars points to a node not
   * defined or null
   */
  if (pdim == 0 || VAR2NODE (a, pdim) == 0
      || VAR2NODE (a, pdim) == NODE_T_TOP)
    return true;

  return false;
}


/* ============================================================ */
/* Extraction of properties */
/* ============================================================ */

/* NOT IMPLEMENTED */
ap_interval_t *
hgraph_bound_linexpr (ap_manager_t * man, hgraph_t * a, ap_linexpr0_t * expr)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_BOUND_LINEXPR, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");

  return NULL;
}

/* NOT IMPLEMENTED */
ap_interval_t *
hgraph_bound_texpr (ap_manager_t * man, hgraph_t * a, ap_texpr0_t * expr)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_BOUND_TEXPR, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");

  return NULL;
}

/* NOT IMPLEMENTED */
ap_interval_t *
hgraph_bound_dimension (ap_manager_t * man, hgraph_t * a, ap_dim_t dim)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_BOUND_DIMENSION, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");

  return NULL;
}

/* NOT IMPLEMENTED */
ap_lincons0_array_t
hgraph_to_lincons_array (ap_manager_t * man, hgraph_t * a)
{
  ap_lincons0_array_t ar;
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_TO_LINCONS_ARRAY, 0);
  ar = ap_lincons0_array_make (1);
  ar.p[0] = ap_lincons0_make_unsat ();

  return ar;
}

/* NOT IMPLEMENTED */
ap_tcons0_array_t
hgraph_to_tcons_array (ap_manager_t * man, hgraph_t * a)
{

  return ap_generic_to_tcons_array (man, a);
}

/* NOT IMPLEMENTED */
ap_interval_t **
hgraph_to_box (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_TO_BOX, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");

  return NULL;
}

/* NOT IMPLEMENTED */
ap_generator0_array_t
hgraph_to_generator_array (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_TO_GENERATOR_ARRAY, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return ap_generator0_array_make (0);
}
