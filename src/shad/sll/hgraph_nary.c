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

/* ============================================================ */
/* Meet and Join */
/* ============================================================ */

/* Used by abstract1_unify */

/* a1 is the callee graph, a2 is the caller graph */

/* variables to be unified point in both graph to node != NODE_T_TOP */
hgraph_t *
hgraph_meet_internal (ap_manager_t * man,
                      hgraph_t * a1,
                      hgraph_t * a2,
                      ap_dimperm_t* perm1,
                      ap_dimperm_t* perm2)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_MEET, 0);
  hgraph_t* r = NULL;
#ifndef NDEBUG1
  fprintf (stdout, "!!!! hgraph_meet_internal: a1=(");
  hgraph_fdump (stdout, man, a1);
  fprintf (stdout, ") and a2=(");
  hgraph_fdump (stdout, man, a2);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif
  // Step 1: check the shared part and compute the number of nodes shared
  //         the shared part is accessible from variables not NULL in both graphs
  size_t v, n1, n2, lshared;
  size_t *shared12 = (size_t*) malloc (a1->size * sizeof (size_t));
  size_t *shared21 = (size_t*) malloc (a2->size * sizeof (size_t));
  memset (shared12, 0, a1->size * sizeof (size_t));
  memset (shared21, 0, a2->size * sizeof (size_t));
  for (v = 0; v < a1->ptrdim; v++)
    {
      n1 = VAR2NODE (a1, v);
      n2 = VAR2NODE (a2, v);
      // dimension is shared iff different from TOP and not only one NULL
      if (n1 != NODE_T_TOP && n2 != NODE_T_TOP)
        {
          if (n1 == NODE_NULL && n2 == NODE_NULL)
            {
              /* nothing to do since NODE_NULL is already shared */
              ;
            }
          else if (n1 != NODE_NULL && n2 != NODE_NULL)
            {
              // shared dimension, check it to be the same and mark nodes
              // check same labeling until a new labeled node
              do
                {
                  shared12[n1] = n2;
                  shared21[n2] = n1;
                  n1 = NODE_NEXT (a1, n1);
                  n2 = NODE_NEXT (a2, n2);
                }
              while (NODE_VAR (a1, n1) == v
                     && NODE_VAR (a2, n2) == v
                     && NODE_VAR_NEXT (a1, n1) == NODE_VAR_NEXT (a2, n2));
              if (NODE_VAR (a1, n1) != NODE_VAR (a2, n2) ||
                  NODE_VAR_NEXT (a1, n1) != NODE_VAR_NEXT (a2, n2))
                {
                  // NOT ISOMORPHIC COMMON PARTS, return
                  goto hgraph_meet_exit;
                }
            }
          else
            {
              // some of nodes is NULL and the other not: NOT ISOMORPHIC
              goto hgraph_meet_exit;
            }
        }
      // else: a non shared variable, continue with a new variable
    }
  lshared = 0;
  for (n1 = 0; n1 < a1->size; n1++)
    if (shared12[n1] != 0)
      lshared++;

#ifndef NDEBUG1
  fprintf (stdout, "!!!! hgraph_meet_internal: (shared=%zu) shared12=(", lshared);
  for (n1 = 0; n1 < a1->size; n1++)
    fprintf (stdout, " %zu -> %zu\n", n1, shared12[n1]);
  fprintf (stdout, ") and shared21=(");
  for (n2 = 0; n2 < a2->size; n2++)
    fprintf (stdout, " %zu -> %zu\n", n2, shared21[n2]);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  // Step 2: build the result by combining graphs
  size_t newsize = a1->size + a2->size - 1 - lshared; // -1 NODE_NULL is shared!
  r = hgraph_alloc_internal (pr, newsize, a1->datadim, a1->ptrdim); // vars are all NODE_T_TOP
  ap_dimperm_t nperm1;
  ap_dimperm_t nperm2;
  ap_dimperm_init (&nperm1, newsize);
  ap_dimperm_set_id (&nperm1);
  ap_dimperm_init (&nperm2, newsize);
  ap_dimperm_set_id (&nperm2);
  // copy all from a2 and fill the added part of perm1
  for (n1 = a1->size, n2 = 0; n2 < a2->size; n2++)
    {
      node_info_copy (&r->info[r->ptrdim + n2], &a2->info[a2->ptrdim + n2]);
      if (n2 != NODE_NULL && shared21[n2] == 0)
        {
          nperm1.dim[n1] = n2;
          n1++;
        }
      else if (shared21[n2] != 0)
        {
          nperm1.dim[shared21[n2]] = n2;
        }
    }
#ifndef NDEBUG1
  fprintf (stdout, "!!!! hgraph_meet_internal: generated r=(");
  hgraph_fdump (stdout, pr->man, r);
  fprintf (stdout, ")\n");
#endif
  // fill sharing information and new nodes ids for a1
  for (n1 = 1, n2 = a2->size; n1 < a1->size; n1++)
    {
      if (shared12[n1] == 0)
        {
          node_info_copy (&r->info[r->ptrdim + n2], &a1->info[a1->ptrdim + n1]);
          nperm1.dim[n1] = n2;
          n2++;
        }
      // else: already done above
    }
#ifndef NDEBUG1
  fprintf (stdout, "!!!! hgraph_meet_internal: generated r=(");
  hgraph_fdump (stdout, pr->man, r);
  fprintf (stdout, ") from a1=(");
  hgraph_fdump (stdout, man, a1);
  fprintf (stdout, ")\n");
#endif
  // then update the successors for the non-shared part of n1, etc.
  for (n2 = a2->size; n2 < newsize; n2++)
    {
      NODE_NEXT (r, n2) = nperm1.dim[NODE_NEXT (r, n2)];
    }
  // init the variable part
  for (v = 0; v < a2->ptrdim; v++)
    {
      size_t n1 = VAR2NODE (a1, v);
      size_t n2 = VAR2NODE (a2, v);
      if (n2 != NODE_T_TOP)
        VAR2NODE (r, v) = n2;
      else if (n1 != NODE_T_TOP)
        VAR2NODE (r, v) = nperm1.dim[n1];
      else
        ERROR ("Bad graphs to unify", return NULL;);
    }

#ifndef NDEBUG1
  fprintf (stdout, "!!!! hgraph_meet_internal: generated nperm1=(");
  ap_dimperm_fprint (stdout, &nperm1);
  fprintf (stdout, "\n) and nperm2=(");
  ap_dimperm_fprint (stdout, &nperm2);
  fprintf (stdout, "\n) and r=(");
  hgraph_fdump (stdout, pr->man, r);
  fprintf (stdout, ")\n");
#endif

  // Step 3: put the result in canonical form and compose permutations
  ap_dimperm_t nodeperm;
  ap_dimperm_init (&nodeperm, newsize);
  ap_dimperm_set_id (&nodeperm);
  hgraph_node_sort (r, 1, &nodeperm);
#ifndef NDEBUG1
  fprintf (stdout, "!!!! hgraph_meet_internal: generated nodeperm=(");
  ap_dimperm_fprint (stdout, &nodeperm);
  fprintf (stdout, ")\n");
#endif
  // alloc perm1, perm2;
  ap_dimperm_init (perm1, newsize);
  ap_dimperm_set_id (perm1);
  ap_dimperm_init (perm2, newsize);
  ap_dimperm_set_id (perm2);
  ap_dimperm_compose (perm1, &nperm1, &nodeperm);
  ap_dimperm_compose (perm2, &nperm2, &nodeperm);
#ifndef NDEBUG1
  fprintf (stdout, "!!!! hgraph_meet_internal: generated perm1=(");
  ap_dimperm_fprint (stdout, perm1);
  fprintf (stdout, ") and perm2=(");
  ap_dimperm_fprint (stdout, perm2);
  fprintf (stdout, ")\n");
#endif
  // free allocated memory
  ap_dimperm_clear (&nperm1);
  ap_dimperm_clear (&nperm2);
  ap_dimperm_clear (&nodeperm);

  /* END */
hgraph_meet_exit:
  free (shared12);
  free (shared21);
  return r;
}

/* TODO: priority 1 */
hgraph_t *
hgraph_meet (ap_manager_t * man, bool destructive, hgraph_t * a1,
             hgraph_t * a2)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_MEET, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a1;
}

/* TODO: priority 1 */
hgraph_t *
hgraph_join (ap_manager_t * man, bool destructive, hgraph_t * a1,
             hgraph_t * a2)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_JOIN, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a1;
}

/* TODO: priority 1 */
hgraph_t *
hgraph_meet_array (ap_manager_t * man, hgraph_t ** tab, size_t size)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_MEET_ARRAY, 0);
  arg_assert (size > 0, return NULL;
              );
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return tab[0];
}

/* TODO: priority 1 */
hgraph_t *
hgraph_join_array (ap_manager_t * man, hgraph_t ** tab, size_t size)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_JOIN_ARRAY, 0);
  arg_assert (size > 0, return NULL;
              );
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return tab[0];
}

/* ============================================================ */
/* Meet constraints and Join generators */
/* ============================================================ */

/* Functions used by meet constraints */
hgraph_t *
hgraph_copy_set_succ (hgraph_internal_t * pr,
                      hgraph_t * a, node_t nsrc, node_t ndst)
{
  hgraph_t *b = hgraph_copy_mem (pr, a);
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->size);
  hgraph_node_set_succ (b, nsrc, ndst, &perm);
  ap_dimperm_clear (&perm);
  return b;
}

hgraph_t *
hgraph_add_succ_fixed (hgraph_internal_t * pr, hgraph_t * a, node_t v,
                       node_t nsucc)
{
  node_t n;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->size);
  hgraph_t *b = hgraph_node_add (pr, a, nsucc, v, 0, &n, &perm);
  ap_dimperm_clear (&perm);
  return b;
}

void
hgraph_add_succ_fixed_all (hgraph_internal_t * pr, hgraph_t * a, size_t v1,
                           size_t v2, hgraph_array_t * r)
{
  node_t n1, n2;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->size);
  hgraph_t *b = hgraph_node_add (pr, a, NODE_T_TOP, v1, 0, &n1, &perm);
  ap_dimperm_clear (&perm);

  hgraph_node_forall (b, n2)
  {
    if (n2 != n1)
      {
        hgraph_t *b1 = hgraph_copy_mem (pr, b);
        hgraph_node_set_var (b1, n2, v2, &perm);
        hgraph_node_set_succ (b1, n1, n2, &perm);
        ap_dimperm_clear (&perm);
        hgraph_array_add (pr, r, true, true, b1);
      }
  }
  hgraph_free_internal (pr, b);
}

hgraph_t *
hgraph_copy_set_var (hgraph_internal_t * pr, hgraph_t * a, node_t n1,
                     size_t v1, node_t n2, size_t v2)
{
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->size);
  hgraph_t *b = hgraph_copy_mem (pr, a);
  hgraph_node_set_var (b, n1, v1, &perm);
  if (v2 && v1 != v2)
    hgraph_node_set_var (b, n2, v2, &perm);
  ap_dimperm_clear (&perm);
  return b;
}

hgraph_t *
hgraph_add_set_var (hgraph_internal_t * pr, hgraph_t * a, node_t nsucc,
                    size_t vmin, size_t vmax)
{
  node_t n;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->size);
  hgraph_t *b = hgraph_node_add (pr, a, NODE_T_TOP, vmin, 0, &n, &perm);
  ap_dimperm_clear (&perm);
  if (vmax && vmin != vmax)
    {
      hgraph_node_set_var (b, n, vmax, &perm);
      ap_dimperm_clear (&perm);
    }
  return b;
}

hgraph_t *
hgraph_copy_set_info (hgraph_internal_t * pr, hgraph_t * a, node_t n,
                      node_t nsucc, size_t v, size_t nn)
{
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->size);
  hgraph_t *b = hgraph_copy_mem (pr, a);
  hgraph_node_set_var (b, n, v, &perm);
  hgraph_node_set_succ (b, n, nsucc, &perm);
  ap_dimperm_clear (&perm);
  return b;
}

hgraph_t *
hgraph_add_set_info (hgraph_internal_t * pr, hgraph_t * a,
                     node_t nsucc, size_t v, size_t nn)
{
  node_t i;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->size);
  hgraph_t *b = hgraph_node_add (pr, a, nsucc, v, nn, &i, &perm);
  ap_dimperm_clear (&perm);
  return b;
}

hgraph_t *
hgraph_set_new_succ (hgraph_internal_t * pr, hgraph_t * a, node_t n, size_t v)
{
  node_t nnew;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->size);
  hgraph_t *b = hgraph_node_add (pr, a, NODE_T_TOP, v, 0, &nnew, &perm);
  ap_dimperm_clear (&perm); /* to clear old size */
  hgraph_node_set_succ (b, n, nnew, &perm);
  ap_dimperm_clear (&perm);
  return b;
}

void
hgraph_set_new_succ_all (hgraph_internal_t * pr, hgraph_t * a, size_t v1,
                         size_t v2, hgraph_array_t * r)
{
  node_t n1, n2;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->size);
  hgraph_t *b = hgraph_node_add (pr, a, NODE_T_TOP, v2, 0, &n2, &perm);
  ap_dimperm_clear (&perm);

  hgraph_node_forall (b, n1)
  {
    if (n1 != n2 && NODE_T_TOP == hgraph_node_get_succ (b, n1))
      {
        hgraph_t *b1 = hgraph_copy_mem (pr, b);
        hgraph_node_set_var (b1, n1, v1, &perm);
        hgraph_node_set_succ (b1, n1, n2, &perm);
        ap_dimperm_clear (&perm);
        hgraph_array_add (pr, r, true, true, b1); /* copy and destroy */
      }
  }
  hgraph_free_internal (pr, b);
}

hgraph_t *
hgraph_add_loop (hgraph_internal_t * pr, hgraph_t * a, size_t vmin,
                 size_t vmax)
{
  node_t n;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->size);
  hgraph_t *b = hgraph_node_add (pr, a, NODE_T_TOP, vmin,
                                 0, &n, &perm);
  ap_dimperm_clear (&perm);
  if (vmax && vmin != vmax)
    hgraph_node_set_var (b, n, vmax, &perm);
  hgraph_node_set_succ (b, n, n, &perm);
  ap_dimperm_clear (&perm);
  return b;
}

hgraph_t *
hgraph_add_edge (hgraph_internal_t * pr, hgraph_t * a, size_t vsrc,
                 size_t vdst)
{
  node_t nsrc, ndst;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->size);
  hgraph_t *b = hgraph_node_add (pr, a, NODE_T_TOP, vdst, 0, &ndst, &perm);
  ap_dimperm_clear (&perm);
  hgraph_t *b1 = b;
  b = hgraph_node_add (pr, b1, ndst, vsrc, 0, &nsrc, &perm);
  hgraph_free_internal (pr, b1);
  ap_dimperm_clear (&perm);
  return b;
}

hgraph_t *
hgraph_add_between (hgraph_internal_t * pr, hgraph_t * a, size_t vsrc,
                    node_t nsrc, node_t ndst)
{
  node_t i;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->size);
  hgraph_t *b = hgraph_node_add (pr, a, ndst, vsrc, 1, &i, &perm);
  ap_dimperm_clear (&perm);
  hgraph_node_set_succ (b, nsrc, i, &perm);
  ap_dimperm_clear (&perm);
  return b;
}

/* Abstract the graph of the disj disjunct of the formula f
 * and return the permutation of nodes generated.
 */
hgraph_t*
hgraph_of_formula (hgraph_internal_t * pr,
                   sh_formula_t* f,
                   size_t disj,
                   ap_dimperm_t* perm)
{
  if (!f || f->form[disj] == NULL)
    return NULL;

#ifndef NDEBUG1
  fprintf (stdout, "!!!! hgraph_of_formula: disjunct %zu\n", disj);
  fflush (stdout);
#endif
  // the environment is given by
  // - the user vars (f->env)
  // - the node vars (f->form[disj].nodes)
  size_t nnodes = (f->form[disj]->nodes == NULL) ? 0 : f->form[disj]->nodes->realdim;
  size_t nvptr = f->env->realdim;
  size_t nvdata = f->env->intdim;

  assert (perm->size == (nnodes + 1));

  hgraph_t* a = hgraph_alloc_internal (pr, nnodes + 1, nvdata, nvptr);

  // Step 1: define graph matrix using edge formulas
  size_t i;
  for (i = 0; i < f->form[disj]->length_eform; i++)
    {
      ap_dim_t nsrc = f->form[disj]->eform[i].src + 1;
      ap_dim_t ndst = (f->form[disj]->eform[i].dst == AP_DIM_MAX) ? NODE_NULL :
              (f->form[disj]->eform[i].dst + 1); // gives 0 for nilNode
      NODE_NEXT (a, nsrc) = ndst;
    }

  // Step 2:
  // - fill info var
  for (i = 0; i < f->form[disj]->length_pform; i++)
    {
      // ptr vars shall start from 0
      ap_dim_t v = DIM2PTR (f->form[disj]->pform[i].var, nvdata);
      ap_dim_t n = (f->form[disj]->pform[i].node == AP_DIM_MAX) ? NODE_NULL :
              (f->form[disj]->pform[i].node + 1); // gives 0 for nilNode
      VAR2NODE (a, v) = n;
      if (n != NODE_NULL)
        {
          NODE_VAR (a, n) = v;
          NODE_VAR_NEXT (a, n) = 0;
        }
    }
  // - label nodes with the minumum ptr var and
  for (i = 0; i < nvptr; i++)
    {
      ap_dim_t n = VAR2NODE (a, i);
      if (n == NODE_NULL || n == NODE_T_TOP) continue;
      ap_dim_t ov = NODE_VAR (a, n);
      if (ov != NULL_DIM && i < ov)
        NODE_VAR (a, n) = i;
    }
  NODE_VAR (a, NODE_NULL) = NULL_DIM;
  NODE_NEXT (a, NODE_NULL) = NODE_NULL;
  // Step 3: signal
  // - error: if nodes with no successors
  for (i = 1; i <= nnodes; i++)
    {
      ap_dim_t n = NODE_NEXT (a, i);
      if (n == NODE_T_TOP)
        NODE_NEXT (a, i) = NODE_NULL;
      // TODO: error
    }
  // - warning: if vars without nodes, set to null
  /* not done to be able to do unify
  for (i = 0; i < nvptr; i++)
    {
      ap_dim_t n = VAR2NODE (a, i);
      if (n == NODE_T_TOP)
        VAR2NODE (a, i) = NODE_NULL;
      // TODO: warning
    }
  */

  // Step 4: sort the nodes from position 1, generate perm (for nodes)

  // sort the graph to obtain the node permutation
  hgraph_node_sort (a, 1, perm);
  return a;
}

/*
 * Compute a meet c and return the set of hgraphs obtained. If not
 * destructive, return the full array of hgraphs built. otherwise, compute an
 * underapproximation.
 */
hgraph_array_t *
hgraph_meet_pcons (hgraph_internal_t * pr, bool destructive, hgraph_t * a,
                   pcons0_t * c)
{
  hgraph_array_t *r = NULL;
  /*
   * TODO: remove ap_dimperm_t perm;
   */
  hgraph_t *b;
  arg_assert (a && c, return NULL;
              );
#ifndef NDEBUG
  fprintf (stdout, "!!!! hgraph_meet_pcons: a = (");
  hgraph_fdump (stdout, pr->man, a);
  fprintf (stdout, "\n\t), pcons = (");
  shape_pcons_fdump (stdout, c);
  fprintf (stdout, "\n\t)\n");
#endif

  if (c->type == DATA_CONS)
    { /* put a in r */
      r = hgraph_array_make (pr, 1);
      r->p[0] = hgraph_copy_internal (pr, a);
      if (destructive)
        hgraph_free_internal (pr, a);
    }
  else if (a && !hgraph_is_bottom (pr->man, a) && c->type < OTHER_CONS)
    {
      node_t vx, nx, nnx, vy, ny, nny;
      /* no more than 1 next dereference */
      if (c->info.ptr.offx != OFFSET_NONE &&
          c->info.ptr.offy != OFFSET_NONE)
        {
          ap_manager_raise_exception (pr->man, AP_EXC_NOT_IMPLEMENTED,
                                      pr->funid, "not implemented");
          return NULL;
        }
      vx = DIM2PTR (c->info.ptr.x, a->datadim);
      nx = VAR2NODE (a, vx);
      nnx = hgraph_node_get_nsucc (a, nx, (c->info.ptr.offx == OFFSET_NONE) ? 0 : 1);
      vy = DIM2PTR (c->info.ptr.y, a->datadim);
      ny = VAR2NODE (a, vy);
      nny = hgraph_node_get_nsucc (a, ny, (c->info.ptr.offy == OFFSET_NONE) ? 0 : 1);
      arg_assert (IS_PTRDIM (vx, a->datadim, a->ptrdim), return NULL;
                  );
      arg_assert (IS_PTRDIM (vy, a->datadim, a->ptrdim), return NULL;
                  );
      if ((nx == 0 && c->info.ptr.offx != OFFSET_NONE)
          || (ny == 0 && c->info.ptr.offy != OFFSET_NONE))
        {
          ERROR ("Null pointer dereference",;
                 );
          return NULL;
        }
      r = hgraph_array_make (pr, 4);
      switch (c->type)
        {
        case EQ_CONS:
          {
            if (nnx != NODE_T_TOP && nny != NODE_T_TOP)
              {
                /* both nodes are defined, check their equality */
                if (nnx == nny)
                  hgraph_array_add (pr, r, true, false, a); /* copy, not distroy */
                /* else bottom result, keep NULL */
              }
            else
              if ((nx != NODE_T_TOP && nnx == NODE_T_TOP && nny != NODE_T_TOP)
                  || (ny != NODE_T_TOP && nny == NODE_T_TOP
                      && nnx != NODE_T_TOP))
              {
                /* fix x->next == y resp. y->next == x */
                size_t n1, n2;
                if (nnx == NODE_T_TOP)
                  {
                    n1 = nx;
                    n2 = ny;
                  }
                else
                  {
                    n1 = ny;
                    n2 = nx;
                  }
                b = hgraph_copy_set_succ (pr, a, n1, n2);
                hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                /*
                 * TODO: remove b = hgraph_copy_mem (pr, a); hgraph_node_set_succ
                 * (b, n1, n2, &perm); ap_dimperm_clear (&perm); hgraph_array_add
                 * (pr, r, hgraph_copy_internal (pr, b)); hgraph_free_internal (pr,
                 * b);
                 */
              }
            else if ((nx != NODE_T_TOP && ny == NODE_T_TOP) ||
                     (ny != NODE_T_TOP && nx == NODE_T_TOP))
              {
                node_t nfixed, nnfixed, nextfixed;
                node_t nany, nextany, vany;
                if (nx != NODE_T_TOP)
                  {
                    nfixed = nx;
                    nnfixed = nnx;
                    nextfixed = (c->info.ptr.offx == OFFSET_NONE) ? 0 : 1;
                    nany = ny;
                    nextany = (c->info.ptr.offy == OFFSET_NONE) ? 0 : 1;
                    vany = vy;
                  }
                else
                  {
                    nfixed = ny;
                    nnfixed = nny;
                    nextfixed = (c->info.ptr.offy == OFFSET_NONE) ? 0 : 1;
                    nany = nx;
                    nextany = (c->info.ptr.offx == OFFSET_NONE) ? 0 : 1;
                    vany = vx;
                  }
                if (nextany == 1)
                  {
                    size_t n;
                    /* case (1) vany->next = vfixed, new node for vany */
                    b = hgraph_add_succ_fixed (pr, a, vany, nfixed);
                    hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                    /*
                     * TODO: remove b = hgraph_node_add (a, nfixed, vany, 0, &n,
                     * &perm); ap_dimperm_clear (&perm); hgraph_array_add (pr, r,
                     * hgraph_copy_internal (pr, b)); hgraph_free_internal (pr, b);
                     */

                    /* case (2) vany->next = vfixed, existing node for vany */
                    hgraph_node_forall (a, n)
                    {
                      if (hgraph_node_get_succ (a, n) == nfixed)
                        {
                          b = hgraph_copy_set_var (pr, a, n, vany, 0, 0);
                          hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                          /*
                           * TODO: remove b = hgraph_copy_mem (pr, a);
                           * hgraph_node_set_var (b, n, vany, &perm); ap_dimperm_clear
                           * (&perm); hgraph_array_add (pr, r, hgraph_copy_internal
                           * (pr, b)); hgraph_free_internal (pr, b);
                           */
                        }
                      else if (hgraph_node_get_succ (a, n) == NODE_T_TOP)
                        {
                          b =
                                  hgraph_copy_set_info (pr, a, n, nfixed, vany, 0);
                          hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                          /*
                           * TODO: remove b = hgraph_copy_mem (pr, a);
                           * hgraph_node_set_var (b, n, vany, &perm);
                           * hgraph_node_set_succ (b, n, nfixed, &perm);
                           * ap_dimperm_clear (&perm); hgraph_array_add (pr, r,
                           * hgraph_copy_internal (pr, b)); hgraph_free_internal (pr,
                           * b);
                           */
                        }
                      else
                        continue;
                    }
                  }
                else if (nnfixed != NODE_T_TOP) /* nextany == 0 */
                  { /* vfixed[->next] = vany */
                    /* set vany to nnfixed */
                    b = hgraph_copy_set_var (pr, a, nnfixed, vany, 0, 0);
                    hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                    /*
                     * TODO: remove b = hgraph_copy_mem (pr, a); hgraph_node_set_var
                     * (b, nnfixed, vany, &perm); ap_dimperm_clear (&perm);
                     * hgraph_array_add (pr, r, hgraph_copy_internal (pr, b));
                     * hgraph_free_internal (pr, b);
                     */
                  }
                else
                  {
                    /* nextany == 0 && nnfixed == NODE_T_TOP */
                    /* case (1): vfixed->next = vany, new node for vany */
                    node_t n;
                    b = hgraph_set_new_succ (pr, a, nfixed, vany);
                    hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                    /*
                     * TODO: remove b = hgraph_node_add (a, NODE_T_TOP, vany, 0, &n,
                     * &perm); ap_dimperm_clear (&perm); hgraph_node_set_succ (b,
                     * nfixed, n, &perm); ap_dimperm_clear (&perm); hgraph_array_add
                     * (pr, r, hgraph_copy_internal (pr, b)); hgraph_free_internal (pr,
                     * b);
                     */

                    /* case (2) vfixed->next = vany, existing node for vany */
                    hgraph_node_forall (a, n)
                    {
                      b = hgraph_copy_set_info (pr, a, n, nfixed, vany, 0);
                      hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                      /*
                       * TODO: remove b = hgraph_copy_mem (pr, a); hgraph_node_set_var
                       * (b, n, vany, &perm); hgraph_node_set_succ (b, nfixed, n,
                       * &perm); ap_dimperm_clear (&perm); hgraph_array_add (pr, r,
                       * hgraph_copy_internal (pr, b)); hgraph_free_internal (pr, b);
                       */
                    }
                  }
              }
            else
              { /* nx == NODE_T_TOP && ny == NODE_T_TOP */
                if (c->info.ptr.offx == OFFSET_NONE &&
                    c->info.ptr.offy == OFFSET_NONE)
                  {
                    /* case (1): set x and y to new node */
                    size_t v1, v2;
                    node_t n;
                    if (vx <= vy)
                      {
                        v1 = vx;
                        v2 = vy;
                      }
                    else
                      {
                        v1 = vy;
                        v2 = vx;
                      }
                    b = hgraph_add_set_var (pr, a, NODE_T_TOP, v1, v2);
                    hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                    /*
                     * TODO: remove b = hgraph_node_add (a, NODE_T_TOP, v1, 0, &n,
                     * &perm); ap_dimperm_clear (&perm); hgraph_node_set_var (b, n, v2,
                     * &perm); ap_dimperm_clear (&perm); hgraph_array_add (pr, r,
                     * hgraph_copy_internal (pr, b)); hgraph_free_internal (pr, b);
                     */

                    /* case (2): set x and y to existing nodes */
                    hgraph_node_forall (a, n)
                    {
                      b = hgraph_copy_set_var (pr, a, n, v1, n, v2);
                      hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                      /*
                       * TODO: remove b = hgraph_copy_mem (pr, a); hgraph_node_set_var
                       * (b, n, v1, &perm); hgraph_node_set_var (b, n, v2, &perm);
                       * hgraph_array_add (pr, r, hgraph_copy_internal (pr, b));
                       * hgraph_free_internal (pr, b); ap_dimperm_clear (&perm);
                       */
                    }
                  }
                else
                  { /* nx == NODE_T_TOP && ny == NODE_T_TOP &&
				 * there is a next */

                    /* v1->next == v2 */
                    size_t v1, v2;
                    node_t n, n2;
                    hgraph_t *b1;
                    if (c->info.ptr.offx != OFFSET_NONE)
                      {
                        v1 = vx;
                        v2 = vy;
                      }
                    else
                      {
                        v1 = vy;
                        v2 = vx;
                      }
                    /* case (1): v1 == v2 and label a new node */
                    b =
                            hgraph_add_loop (pr, a, ((v1 <= v2) ? v1 : v2),
                                             ((v1 < v2) ? v2 : v1));
                    hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                    /*
                     * TODO: remove b = hgraph_node_add (a, NODE_T_TOP, ((v1 <= v2) ? v1
                     * : v2), 0, &n, &perm); ap_dimperm_clear (&perm); if (v1 != v2)
                     * hgraph_node_set_var (b, n, ((v1 < v2) ? v2 : v1), &perm);
                     * hgraph_node_set_succ (b, n, n, &perm); ap_dimperm_clear (&perm);
                     * hgraph_array_add (pr, r, hgraph_copy_internal (pr, b));
                     * hgraph_free_internal (pr, b);
                     */
                    /* case (2): v1 != v2 and both label new nodes */
                    b = hgraph_add_edge (pr, a, v1, v2);
                    hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                    /*
                     * TODO: remove b = hgraph_node_add (a, NODE_T_TOP, v2, 0, &n2,
                     * &perm); ap_dimperm_clear (&perm); b1 = b; b = hgraph_node_add (b1,
                     * n2, v1, 0, &n, &perm); hgraph_free_internal (pr, b1);
                     * ap_dimperm_clear (&perm); hgraph_array_add (pr, r,
                     * hgraph_copy_internal (pr, b)); hgraph_free_internal (pr, b);
                     */

                    /* case (3): v1 and v2 are existing nodes */
                    hgraph_node_forall (a, n)
                    {
                      if (n != hgraph_node_get_null (a))
                        {

                          hgraph_node_forall (a, n2)
                          {
                            if (n2 == hgraph_node_get_succ (a, n))
                              {
                                b =
                                        hgraph_copy_set_var (pr, a, n, v1, n2, v2);
                                hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                                /*
                                 * TODO: remove b = hgraph_copy_mem (pr, a);
                                 * hgraph_node_set_var (b, n, v1, &perm); hgraph_node_set_var
                                 * (b, n2, v2, &perm); ap_dimperm_clear (&perm);
                                 * hgraph_array_add (pr, r, hgraph_copy_internal (pr, b));
                                 * hgraph_free_internal (pr, b);
                                 */
                              }
                          }
                        }
                    }
                    /* case (4): v1 is new, v2 is an existing node */
                    hgraph_add_succ_fixed_all (pr, a, v1, v2, r);
                    /*
                     * TODO: remove b = hgraph_node_add (a, NODE_T_TOP, v1, 0, &n,
                     * &perm); ap_dimperm_clear (&perm); hgraph_node_forall (b, n2) { if
                     * (n2 != n) { b1 = hgraph_copy_mem (pr, b); hgraph_node_set_var (b1,
                     * n2, v2, &perm); hgraph_node_set_succ (b1, n, n2, &perm);
                     * ap_dimperm_clear (&perm); hgraph_array_add (pr, r,
                     * hgraph_copy_internal (pr, b1)); hgraph_free_internal (pr, b1); } }
                     * hgraph_free_internal (pr, b);
                     */
                    /* case (5): v1 is old with succ undefined, v2 is new */
                    hgraph_set_new_succ_all (pr, a, v1, v2, r);
                    /*
                     * TODO: remove b = hgraph_node_add (a, NODE_T_TOP, v2, 0, &n2,
                     * &perm); ap_dimperm_clear (&perm); hgraph_node_forall (b, n) { if
                     * (n != n2 && NODE_T_TOP == hgraph_node_get_succ (b, n)) { b1 =
                     * hgraph_copy_mem (pr, b); hgraph_node_set_var (b1, n, v1, &perm);
                     * hgraph_node_set_succ (b1, n, n2, &perm); ap_dimperm_clear (&perm);
                     * hgraph_array_add (pr, r, hgraph_copy_internal (pr, b1));
                     * hgraph_free_internal (pr, b1); } } hgraph_free_internal (pr, b);
                     */
                  }
              }
            break;
          }
        case NE_CONS:
          {
            node_t i;
            if (nnx != NODE_T_TOP && nny != NODE_T_TOP)
              {
                /* check only equality of nodes */
                if (nnx != nny)
                  hgraph_array_add (pr, r, true, false, a);
                else if ((c->info.ptr.offx != OFFSET_NONE &&
                          c->info.ptr.offy == OFFSET_NONE) ||
                         (c->info.ptr.offx == OFFSET_NONE &&
                          c->info.ptr.offy != OFFSET_NONE))
                  {
                    /* v1->next != v2 : create a new node s.t.
                     * v1-->new node --> v2 */
                    node_t n1, n2;
                    size_t v1, v2;
                    if (c->info.ptr.offx != OFFSET_NONE)
                      {
                        n1 = nx;
                        v1 = vx;
                        n2 = ny;
                        v2 = vy;
                      }
                    else
                      {
                        n1 = ny;
                        v1 = vy;
                        n2 = nx;
                        v2 = vx;
                      }
                    /*
                     * TODO: remove b = hgraph_node_add (a, n2, v1, 1, &i, &perm);
                     * ap_dimperm_clear (&perm); hgraph_node_set_succ (b, n1, i, &perm);
                     * ap_dimperm_clear (&perm); hgraph_array_add (pr, r,
                     * hgraph_copy_internal (pr, b)); hgraph_free_internal (pr, b);
                     */
                    b = hgraph_add_between (pr, a, v1, n1, n2);
                    hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                  }
                /* else bottom result, so keep NULL */
              }
            else
              if ((nx != NODE_T_TOP && nnx == NODE_T_TOP && nny != NODE_T_TOP)
                  || (ny != NODE_T_TOP && nny == NODE_T_TOP
                      && nnx != NODE_T_TOP))
              { /* v1->next != v2, v1->next undefined, v2
				 * fixed */
                node_t n1, n2;
                size_t v1, v2;
                if (nnx == NODE_T_TOP)
                  {
                    n1 = nx;
                    v1 = vx;
                    n2 = ny;
                    v2 = vy;
                  }
                else
                  {
                    n1 = ny;
                    v1 = vy;
                    n2 = nx;
                    v2 = vx;
                  }
                /* case (1): fix v1->next to some new node, so different from n2 */
                /*
                 * TODO: remove b = hgraph_node_add (a, NODE_T_TOP, v1, 1, &i, &perm);
                 * ap_dimperm_clear (&perm); hgraph_array_add (pr, r,
                 * hgraph_copy_internal (pr, b)); hgraph_free_internal (pr, b);
                 */
                b = hgraph_add_set_info (pr, a, NODE_T_TOP, v1, 1);
                hgraph_array_add (pr, r, true, true, b); /* copy and destroy */

                /*
                 * case (2): fix v1->next to some existing node, different from
                 * n2
                 */
                hgraph_node_forall (a, i)
                {
                  if (n2 != hgraph_node_get_succ (a, i))
                    {
                      /*
                       * TODO: remove b = hgraph_copy_mem (pr, a); hgraph_node_set_var
                       * (b, i, v1, &perm); ap_dimperm_clear (&perm); hgraph_array_add
                       * (pr, r, hgraph_copy_internal (pr, b)); hgraph_free_internal (pr,
                       * b);
                       */
                      b = hgraph_copy_set_var (pr, a, i, v1, 0, 0);
                      hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                    }
                }
              }
            else if ((nx == NODE_T_TOP && nny != NODE_T_TOP) ||
                     (ny == NODE_T_TOP && nnx != NODE_T_TOP))
              { /* case v1 != v2[->next] for v1 not yet fixed
				 * and v2[->next] fixed */
                node_t n1, n2;
                size_t v1, v2;
                if (nx == NODE_T_TOP)
                  {
                    n1 = nx;
                    v1 = vx;
                    n2 = nny;
                    v2 = vy;
                  }
                else
                  {
                    n1 = ny;
                    v1 = vy;
                    n2 = nnx;
                    v2 = vx;
                  }
                /* case (1): fix v1 to some new node, so different from n2 */
                b = hgraph_add_set_info (pr, a, NODE_T_TOP, v1, 0);
                hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                /*
                 * TODO: remove b = hgraph_node_add (a, NODE_T_TOP, v1, 0, &i, &perm);
                 * ap_dimperm_clear (&perm); hgraph_array_add (pr, r,
                 * hgraph_copy_internal (pr, b)); hgraph_free_internal (pr, b);
                 */

                /* case (2): fix v1 to some existing node, different from n2 */
                hgraph_node_forall (a, i)
                {
                  if (n2 != i)
                    {
                      b = hgraph_copy_set_var (pr, a, i, v1, 0, 0);
                      hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                      /*
                       * TODO: remove b = hgraph_copy_mem (pr, a); hgraph_node_set_var
                       * (b, i, v1, &perm); ap_dimperm_clear (&perm); hgraph_array_add
                       * (pr, r, hgraph_copy_internal (pr, b)); hgraph_free_internal (pr,
                       * b);
                       */
                    }
                }
              }
            else
              { /* TODO: the same cases like for equality,
				 * but more combinatorial explosion. */
                ap_manager_raise_exception (pr->man, AP_EXC_NOT_IMPLEMENTED,
                                            pr->funid, "not implemented");
                return NULL;
              }
            break;
          }

        case REACH_CONS:
          {
            node_t i;
            /* no next dereference */
            if (c->info.ptr.offx != OFFSET_NONE ||
                c->info.ptr.offy != OFFSET_NONE)
              {
                ERROR
                        ("Not a reachability predicate (no next dereferencing)",;
                         );
                return NULL;
              }
            if (nx != NODE_T_TOP && ny != NODE_T_TOP)
              {
                /* check only reachability of nodes */
                bool found = hgraph_node_is_reachable (a, nx, ny);
                if (found)
                  hgraph_array_add (pr, r, true, false, a);
                /* else bottom result, keep NULL entry */
              }
            else if (nx != NODE_T_TOP && ny == NODE_T_TOP)
              {
                /* set y to some the next successor of nx */
                /* TODO: the test is more complex because of lasso from x */
                int j = 0;
                i = NODE_T_TOP;
                ny = hgraph_node_get_succ (a, nx);
                while (ny != NODE_T_TOP && ny != i && j <= 3)
                  {
                    b = hgraph_copy_set_var (pr, a, ny, vy, 0, 0);
                    hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                    /*
                     * TODO: remove b = hgraph_copy_mem (pr, a); hgraph_node_set_var (b,
                     * ny, vy, &perm); ap_dimperm_clear (&perm); hgraph_array_add (pr, r,
                     * hgraph_copy_internal (pr, b)); hgraph_free_internal (pr, b);
                     */
                    i = ny;
                    j++;
                    ny = hgraph_node_get_succ (a, i);
                  }
              }
            else if (nx == NODE_T_TOP && ny != NODE_T_TOP)
              {
                /* case (1): set x to some new node with the next successor at ny */
                /*
                 * TODO: the case is more complex because several nodes may exist to
                 * ny
                 */
                b = hgraph_add_set_info (pr, a, ny, vx, 0);
                hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                /*
                 * TODO: remove b = hgraph_node_add (a, ny, vx, 0, &nx, &perm);
                 * ap_dimperm_clear (&perm); hgraph_array_add (pr, r,
                 * hgraph_copy_internal (pr, b)); hgraph_free_internal (pr, b);
                 */
                /* case (2): set x to some existing node predecessor of ny */

                /*
                 * TODO: the case is more complex, bellow only immediate predecessor is
                 * chosen
                 */
                hgraph_node_forall (a, i)
                {
                  if (ny == hgraph_node_get_succ (a, i))
                    {
                      b = hgraph_copy_set_var (pr, a, i, vx, 0, 0);
                      hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                      /*
                       * TODO: remove b = hgraph_copy_mem (pr, a); hgraph_node_set_var
                       * (b, i, vx, &perm); ap_dimperm_clear (&perm); hgraph_array_add
                       * (pr, r, hgraph_copy_internal (pr, b)); hgraph_free_internal (pr,
                       * b);
                       */
                    }
                }
              }
            else
              {
                /* nx and ny are undefined, build nx ---> ny */

                size_t v1, v2;
                if (vx <= vy)
                  {
                    v1 = vx;
                    v2 = vy;
                  }
                else
                  {
                    v1 = vy;
                    v2 = vx;
                  }
                /* case (1): x and y are the same new node */
                b = hgraph_add_loop (pr, a, v1, v2); /* v1 <= v2 */
                hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                /*
                 * TODO: remove b = hgraph_node_add (a, NODE_T_TOP, v1, 0, &i, &perm);
                 * ap_dimperm_clear (&perm); hgraph_node_set_succ (b, i, i, &perm);
                 * hgraph_node_set_var (b, i, v2, &perm); ap_dimperm_clear (&perm);
                 * hgraph_array_add (pr, r, hgraph_copy_internal (pr, b));
                 * hgraph_free_internal (pr, b);
                 */
                /* case (2): x and y are different new nodes */
                b = hgraph_add_edge (pr, a, vx, vy);
                hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                /*
                 * TODO: remove b = hgraph_node_add (a, NODE_T_TOP, vy, 0, &ny, &perm);
                 * hgraph_t *b1 = b; b = hgraph_node_add (b1, ny, vx, 0, &nx, &perm);
                 * hgraph_free_internal (pr, b1); ap_dimperm_clear (&perm);
                 * hgraph_array_add (pr, r, hgraph_copy_internal (pr, b));
                 * hgraph_free_internal (pr, b);
                 */
                /* case (3): x new, y existing */
                hgraph_add_succ_fixed_all (pr, a, vx, vy, r);
                /*
                 * TODO: remove b = hgraph_node_add (a, NODE_T_TOP, vx, 0, &nx, &perm);
                 * ap_dimperm_clear (&perm); hgraph_node_forall (b, ny) { if (nx != ny) {
                 * b1 = hgraph_copy_mem (pr, b); hgraph_node_set_var (b1, ny, vy, &perm);
                 * hgraph_node_set_succ (b1, nx, ny, &perm); ap_dimperm_clear (&perm);
                 * hgraph_array_add (pr, r, hgraph_copy_internal (pr, b1));
                 * hgraph_free_internal (pr, b1); } } hgraph_free_internal (pr, b);
                 */
                /* case (4): x existing but no succ, y new */
                /* TODO: simple case here, no intermediate node between x and y */
                hgraph_set_new_succ_all (pr, a, vx, vy, r);
                /*
                 * TODO: remove b = hgraph_node_add (a, NODE_T_TOP, vy, 0, &ny, &perm);
                 * ap_dimperm_clear (&perm); hgraph_node_forall (b, nx) { if (nx != ny &&
                 * NODE_T_TOP == hgraph_node_get_succ (b, nx)) { b1 = hgraph_copy_mem
                 * (pr, b); hgraph_node_set_var (b1, nx, vx, &perm); hgraph_node_set_succ
                 * (b1, nx, ny, &perm); ap_dimperm_clear (&perm); hgraph_array_add (pr,
                 * r, hgraph_copy_internal (pr, b1)); hgraph_free_internal (pr, b1); } }
                 * hgraph_free_internal (pr, b);
                 */
                /* case (5): both x and y are existing nodes */

                /* TODO: more complex cases with more nodes between x and y */
                hgraph_node_forall (a, nx)
                {
                  if (nx != hgraph_node_get_null (a))
                    {

                      hgraph_node_forall (a, ny)
                      {
                        if (ny == hgraph_node_get_succ (a, nx))
                          {
                            b = hgraph_copy_set_var (pr, a, nx, vx, ny, vy);
                            hgraph_array_add (pr, r, true, true, b); /* copy and destroy */
                          }
                      }
                    }
                }
                break;
              }

          case ISO_CONS:
            {
              if (nnx != NODE_T_TOP && nny != NODE_T_TOP)
                {
                  /* both nodes are defined, check their isomorphism */
                  if (hgraph_node_are_iso (pr, a, nnx, nny))
                    hgraph_array_add (pr, r, true, false, a); /* copy, not distroy */
                  /* else bottom result, keep NULL */
                }
              break;
            }
          default:
            break;
          }
          if (destructive)
            {
              hgraph_free_internal (pr, a);
              hgraph_array_resize (pr, r, 4);
            }
        }
    }
#ifndef NDEBUG
  fprintf (stdout, "\n!!!! hgraph_meet_pcons: returns ");
  if (!r)
    fprintf (stdout, "null\n");
  else
    {
      size_t i;
      fprintf (stdout, "array of size %zu and elements \n\t", r->size);
      for (i = 0; i < r->size; i++)
        if (r->p[i])
          hgraph_fdump (stdout, pr->man, r->p[i]);
        else
          fprintf (stdout, "NULL\n");
    }
#endif
  return r;
}

/*
 * Compute the effect of constraints in array (of size elements) on the
 * hgraph a. Better performances are obtained if array is sorted, @see
 * <code>shape_pcons_array_of_lincons_array</code>
 */
hgraph_array_t *
hgraph_meet_pcons_array (hgraph_internal_t * pr, bool destructive,
                         hgraph_t * a, pcons0_array_t * array)
{
  hgraph_array_t *r;
  size_t i, j;
  r = hgraph_array_make (pr, 4);
  r->p[0] = hgraph_copy_mem (pr, a);
  for (i = 0; i < array->size && r; i++)
    {
      //TODO:apply the constraints but not iteratively
      hgraph_array_t *rc = NULL;
      for (j = 0; j < r->size && r->p[j] != NULL; j++)
        {
          hgraph_array_t *lr =
                  hgraph_meet_pcons (pr, false, r->p[j], array->p[i]);
          rc = hgraph_array_merge (pr, true, rc, lr);
        }
      hgraph_array_clear (pr, r);
      if (destructive && rc)
        hgraph_array_resize (pr, rc, 4);
      r = rc;
    }
  if (destructive)
    hgraph_free_internal (pr, a);
  return r;
}

hgraph_t *
hgraph_meet_lincons_array (ap_manager_t * man,
                           bool destructive, hgraph_t * a,
                           ap_lincons0_array_t * array)
{
  if (hgraph_is_bottom (man, a) || !array || array->size == 0)
    /* nothing to do */
    return (destructive) ? a : hgraph_copy (man, a);
  else
    {
      hgraph_t *b;
      pcons0_array_t *arr;
      hgraph_array_t *r;
      hgraph_internal_t *pr =
              hgraph_init_from_manager (man, AP_FUNID_MEET_LINCONS_ARRAY, 0);
      if (!destructive)
        b = hgraph_copy_mem (pr, a);
      else
        b = a;
      /* compute in arr the constraints sorted */
      arr =
              shape_pcons_array_of_lincons_array (pr, array, a->datadim, a->ptrdim);
#ifndef NDEBUG
      printf ("\n========== meet with pcons array:\n");
      shape_pcons_array_fdump (stdout, arr);
#endif
      /* go */
      /* true below means to compute only one element and return it in b */
      r = hgraph_meet_pcons_array (pr, true, b, arr);
      if (r && r->size > 0)
        {
          b = r->p[0];
          r->p[0] = NULL;
          hgraph_array_clear (pr, r);
        }
      else
        b = hgraph_bottom (man, a->datadim, a->ptrdim);
      return b;
    }
}

hgraph_t *
hgraph_meet_tcons_array (ap_manager_t * man,
                         bool destructive, hgraph_t * a,
                         ap_tcons0_array_t * array)
{
  if (hgraph_is_bottom (man, a) || !array || array->size == 0)
    /* nothing to do */
    return (destructive) ? a : hgraph_copy (man, a);
  else
    {
      hgraph_t *b;
      hgraph_array_t *r;
      pcons0_array_t *arr;
      hgraph_internal_t *pr =
              hgraph_init_from_manager (man, AP_FUNID_MEET_TCONS_ARRAY, 0);
      if (!destructive)
        b = hgraph_copy_mem (pr, a);
      else
        b = a;
      /* compute in arr the constraints sorted */
      arr =
              shape_pcons_array_of_tcons_array (pr, array, a->datadim, a->ptrdim);
      /* go */
      /* true below means to compute only one element and return it in b */
      r = hgraph_meet_pcons_array (pr, true, b, arr);
      if (r && r->size > 0)
        {
          b = r->p[0];
          r->p[0] = NULL;
          hgraph_array_clear (pr, r);
        }
      else
        b = hgraph_bottom (pr->man, a->datadim, a->ptrdim);
      hgraph_array_clear (pr, r);
      return b;
    }
}


/* ============================================================ */
/* Meet constraints and Join generators */
/* ============================================================ */

/* NOT IMPLEMENTED */
hgraph_t *
hgraph_add_ray_array (ap_manager_t * man,
                      bool destructive, hgraph_t * a,
                      ap_generator0_array_t * array)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_ADD_RAY_ARRAY, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a;
}



/* ============================================================ */
/* Widening, Narrowing */
/* ============================================================ */

/* TODO: priority 1 */
hgraph_t *
hgraph_widening (ap_manager_t * man, hgraph_t * a1, hgraph_t * a2)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_WIDENING, 0);
  arg_assert (a1->ptrdim == a2->ptrdim, return NULL;
              );
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a2;
}

/* TODO: priority 1 */
hgraph_t *
hgraph_widening_thresholds (ap_manager_t * man,
                            hgraph_t * a1, hgraph_t * a2,
                            ap_scalar_t ** array, size_t nb)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_WIDENING, nb + 1);
  arg_assert (a1->ptrdim == a2->ptrdim, return NULL;
              );
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a2;
}

/* NOT IMPLEMENTED */
hgraph_t *
hgraph_narrowing (ap_manager_t * man, hgraph_t * a1, hgraph_t * a2)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_WIDENING, 0);
  arg_assert (a1->ptrdim == a2->ptrdim, return NULL;
              );
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a2;
}
