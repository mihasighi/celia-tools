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


#include "hgraph_internal.h"
#include "ushape.h"
#include "ushape_internal.h"
#include "shape_macros.h"
#include "ap_generic.h"

/* ============================================================ */
/* Meet and Join */

/* ============================================================ */

ushape_t *
ushape_meet (ap_manager_t * man, bool destructive, ushape_t * a1,
             ushape_t * a2)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_MEET, 0);
  ushape_t *r;

  if (ushape_is_bottom (man, a1))
    return (destructive) ? a1 : ushape_bottom (man, a1->datadim, a1->ptrdim);
  if (ushape_is_bottom (man, a2))
    return (destructive) ? a2 : ushape_bottom (man, a2->datadim, a2->ptrdim);
#ifndef NDEBUG1
  fprintf (stdout, "\n====ushape_meet: with algo=%d and dimensions a1=(%d,%d), a2=(%d,%d)\n",
           pr->meet_algo, a1->datadim, a1->ptrdim, a2->datadim, a2->ptrdim);
  fflush (stdout);
#endif
  size_t i;
  bool isbottom = false;
  if (!hgraph_is_equal (a1->h, a2->h)
      && pr->meet_algo == 0)
    r = NULL; /* bottom */
  else if (pr->meet_algo == -1)
    {
      /****** Unify *******/
      ap_dimperm_t perm1;
      ap_dimperm_t perm2;
      r = ushape_alloc_internal (pr, a1->datadim, a1->ptrdim);
      r->h = hgraph_meet_internal (man, a1->h, a2->h, &perm1, &perm2);
      if (r->h == NULL)
        {
          // if meet cannot be done
          isbottom = true;
          // perms have been deallocated in meet
        }
      else
        {
          // Step1: build dimchange for the added nodes to a1 and a2
          size_t nadded1 = r->h->size - a1->h->size;
          size_t nadded2 = r->h->size - a2->h->size;
          ap_dimchange_t dimchange1;
          ap_dimchange_t dimchange2;
          ap_dimchange_init (&dimchange1, 0, nadded1);
          for (i = 0; i < nadded1; i++) dimchange1.dim[i] = a1->datadim + a1->h->size;
          ap_dimchange_init (&dimchange2, 0, nadded2);
          for (i = 0; i < nadded2; i++) dimchange2.dim[i] = a2->datadim + a2->h->size;
          // apply dimperm_dimchange_n for a1 and a2 to obtain a1p, a2p
          ushape_t *a1p = ushape_alloc_internal (pr, r->datadim, r->ptrdim);
          a1p->h = hgraph_copy_internal (pr, r->h);
          a1p->closed = NULL;
          ushape_apply_dimperm_dimchange_n (pr, a1, a1p, &perm1, &dimchange1);
          ushape_t *a2p = ushape_alloc_internal (pr, r->datadim, r->ptrdim);
          a2p->h = hgraph_copy_internal (pr, r->h);
          a2p->closed = NULL;
          ushape_apply_dimperm_dimchange_n (pr, a2, a2p, &perm2, &dimchange2);
          for (i = 0; i < pr->size_scons && !isbottom; i++)
            {
              r->scons[i] =
                      ap_abstract0_meet (pr->man_scons[i], false,
                                         (a1p->scons) ? a1p->scons[i] : NULL,
                                         (a2p->scons) ? a2p->scons[i] : NULL);
              isbottom = ap_abstract0_is_bottom (pr->man_scons[i], r->scons[i]);
            }
          // free all allocated data
          ap_dimchange_clear (&dimchange1);
          ap_dimchange_clear (&dimchange2);
          ushape_free_internal (pr, a1p);
          ushape_free_internal (pr, a2p);
          // strengthen
          ap_dimperm_clear (&perm1);
          ap_dimperm_clear (&perm2);
        }
    }
  else if (pr->meet_algo == -2 && ushape_is_top (man, a1))
    {
      /***** Meet with a pre-condition, return a2 */
      r = ushape_copy (man, a2);
    }
  else if (pr->meet_algo == -2)
    {
      /***** Meet with a (possibly incomplete) specification a2 ******/
      /* (TODO: beter semantics) test inclusion a1 <= a2 and 
       * return a2 with the graph of a1 
       */
      ap_dimperm_t perm2;
      r = ushape_alloc_internal (pr, a1->datadim, a1->ptrdim);
      r->h = hgraph_copy_internal (pr, a1->h);
      if (!hgraph_is_spec (man, a1->h, a2->h, &perm2))
        // perm2 is allocated above, if needed
        {
          // if meet cannot be done
          isbottom = true;
          // perms have been deallocated in hgraph_is_spec
        }
      else
        {
          // Step1: build dimchange for the added nodes to a2
          size_t nadded2 = r->h->size - a2->h->size;
          ushape_t *a2p = NULL;
          if (nadded2 > 0)
            {
              ap_dimchange_t dimchange2;
              ap_dimchange_init (&dimchange2, 0, nadded2);
              for (i = 0; i < nadded2; i++)
                dimchange2.dim[i] = a2->datadim + a2->h->size;
              // apply dimperm_dimchange_n for a2 to obtain a2p
              a2p = ushape_alloc_internal (pr, r->datadim, r->ptrdim);
              a2p->h = hgraph_copy_internal (pr, r->h);
              a2p->closed = NULL;
              ushape_apply_dimperm_dimchange_n (pr, a2, a2p, &perm2, &dimchange2);
              ap_dimchange_clear (&dimchange2);
              ap_dimperm_clear (&perm2);
            }
          else
            a2p = a2;
          for (i = 0; i < pr->size_scons && !isbottom; i++)
            {
              if (ap_abstract0_is_leq (pr->man_scons[i],
                                       (a1->scons) ? a1->scons[i] : NULL,
                                       (a2p->scons) ? a2p->scons[i] : NULL))
                r->scons[i] =
                      ap_abstract0_copy (pr->man_scons[i], a2p->scons[i]);
              else
                isbottom = true;
            }
          // free all allocated data
          if (nadded2 > 0)
            ushape_free_internal (pr, a2p);
        }
    }
  else
    {
      r = ushape_alloc_internal (pr, a1->datadim, a1->ptrdim);
      r->h = hgraph_copy_internal (pr, a1->h);
      for (i = 0; i < pr->size_scons && !isbottom; i++)
        {
          r->scons[i] =
                  ap_abstract0_meet (pr->man_scons[i], false,
                                     (a1->scons) ? a1->scons[i] : NULL,
                                     (a2->scons) ? a2->scons[i] : NULL);
          isbottom = ap_abstract0_is_bottom (pr->man_scons[i], r->scons[i]);
        }
    }
  if (isbottom)
    {
      ushape_free_internal (pr, r);
      r = ushape_bottom (man, a1->datadim, a1->ptrdim);
    }
  if (destructive)
    {
      ushape_free_internal (pr, a1);
      ushape_free_internal (pr, a2);
    }
#ifndef NDEBUG1
  fprintf (stdout, "\n====ushape_meet: with algo=%d returns r=(\n", pr->meet_algo);
  ushape_fdump (stdout, man, r);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif
  return r;
}

ushape_t *
ushape_join (ap_manager_t * man, bool destructive, ushape_t * a1,
             ushape_t * a2)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_JOIN, 0);

  if (ushape_is_bottom (man, a1) && ushape_is_bottom (man, a2))
    return (destructive) ? a1 : ushape_bottom (man, a1->datadim, a1->ptrdim);
  if (ushape_is_bottom (man, a1))
    return (destructive) ? a2 : ushape_copy (man, a2);
  if (ushape_is_bottom (man, a2))
    return (destructive) ? a1 : ushape_copy (man, a1);

  if (!hgraph_is_equal (a1->h, a2->h))
    return NULL; /* bottom */
  else
    {
      size_t i;
      ushape_t *r = ushape_alloc_internal (pr, a1->datadim, a1->ptrdim);

      r->h = hgraph_copy_internal (pr, a1->h);
      for (i = 0; i < pr->size_scons; i++)
        {
          r->scons[i] =
                  ap_abstract0_join (pr->man_scons[i], false, a1->scons[i],
                                     a2->scons[i]);
        }
      if (destructive)
        {
          ushape_free_internal (pr, a1);
          ushape_free_internal (pr, a2);
        }
      return r;
    }
}

ushape_t *
ushape_meet_array (ap_manager_t * man, ushape_t ** tab, size_t size)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_JOIN, 0);
  arg_assert (size > 0, return NULL;
              );
  ushape_t *r = ushape_copy_internal (pr, tab[0]);
  size_t i;
  bool isbottom = false;

  for (i = 1; i < size && !isbottom; i++)
    {
      ushape_t *rr = ushape_meet (man, false, tab[i], r);
      ushape_free (man, r);
      isbottom = ushape_is_bottom (man, rr);
      r = rr;
    }
  if (isbottom)
    {
      ushape_free (man, r);
      r = NULL;
    }
  return r;
}

ushape_t *
ushape_join_array (ap_manager_t * man, ushape_t ** tab, size_t size)
{
  ushape_internal_t *pr =
          ushape_init_from_manager (man, AP_FUNID_JOIN_ARRAY, 0);
  arg_assert (size > 0, return NULL;
              );

  ushape_t *r = ushape_copy_internal (pr, tab[0]);
  size_t i;
  bool istop = false;

  for (i = 1; i < size && !istop; i++)
    {
      ushape_t *rr = ushape_join (man, false, tab[i], r);
      ushape_free (man, r);
      istop = ushape_is_top (man, rr);
      r = rr;
    }
  if (istop)
    {
      ushape_free (man, r);
      r = ushape_top (man, tab[0]->datadim, tab[0]->ptrdim);
    }
  return r;
}


/* ============================================================ */
/* Meet constraints and Join generators */
/* ============================================================ */

/* Functions used by meet constraints */
ushape_t *
ushape_copy_set_succ (ushape_internal_t * pr,
                      ushape_t * a, node_t nsrc, node_t ndst)
{
  ap_dimperm_t perm;
  ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  hgraph_t *h = hgraph_copy_mem (pr, a->h);
  ap_dimperm_init (&perm, a->h->size);
  hgraph_node_set_succ (h, nsrc, ndst, &perm);
  b->h = hgraph_copy_internal (pr, h);
  hgraph_free_internal (pr, h);
  ushape_apply_dimperm (pr, a, b, &perm);
  return b;
}

ushape_t *
ushape_add_succ_fixed (ushape_internal_t * pr, ushape_t * a, node_t v,
                       node_t nsucc)
{
  node_t n;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->h->size);
  ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  hgraph_t *h = hgraph_node_add (pr, a->h, nsucc, v, 0, &n, &perm);
  b->h = hgraph_copy_internal (pr, h);
  hgraph_free_internal (pr, h);
  ushape_apply_dimperm_dimchange (pr, a, b, &perm, n);
  return b;
}

size_t
ushape_add_succ_fixed_all (ushape_internal_t * pr, ushape_t * a, size_t v1,
                           size_t v2, ushape_array_t * r, size_t rsize)
{
  node_t n1, n2;
  size_t i;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->h->size);
  ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  hgraph_t *h = hgraph_node_add (pr, a->h, NODE_T_TOP, v1, 0, &n1, &perm);
  ushape_apply_dimperm_dimchange (pr, a, b, &perm, n1);
  i = 0;

  hgraph_node_forall (h, n2)
  {
    if (n2 != n1)
      {
        ushape_t *b1 = ushape_alloc_internal (pr, b->datadim, b->ptrdim);
        hgraph_t *h1 = hgraph_copy_mem (pr, h);
        ap_dimperm_init (&perm, h1->size);
        hgraph_node_set_var (h1, n2, v2, &perm);
        hgraph_node_set_succ (h1, n1, n2, &perm);
        b->h = hgraph_copy_internal (pr, h1);
        ushape_apply_dimperm (pr, b, b1, &perm);
        b1->h = hgraph_copy_internal (pr, h1);
        i += ushape_array_add (pr, true, r, (rsize + i), true, true, b1);
        hgraph_free_internal (pr, h1);
        // ushape_free_internal(pr, b1);
      }
  }
  hgraph_free_internal (pr, h);
  ushape_free_internal (pr, b);
  return i;
}

ushape_t *
ushape_copy_set_var (ushape_internal_t * pr, ushape_t * a, node_t n1,
                     size_t v1, node_t n2, size_t v2)
{
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->h->size);
  ap_dimperm_set_id (&perm);
  ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  hgraph_t *h = hgraph_copy_mem (pr, a->h);
  hgraph_node_set_var (h, n1, v1, &perm);
  if (v2 && v1 != v2)
    hgraph_node_set_var (h, n2, v2, &perm);
  b->h = hgraph_copy_internal (pr, h);
  hgraph_free_internal (pr, h);
  ushape_apply_dimperm (pr, a, b, &perm);
  return b;
}

ushape_t *
ushape_add_set_var (ushape_internal_t * pr, ushape_t * a, node_t nsucc,
                    size_t vmin, size_t vmax)
{
  node_t n;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->h->size);
  ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  hgraph_t *h = hgraph_node_add (pr, a->h, NODE_T_TOP, vmin, 0, &n, &perm);
  ushape_apply_dimperm_dimchange (pr, a, b, &perm, n);
  if (vmax && vmin != vmax)
    {
      ap_dimperm_init (&perm, h->size);
      hgraph_node_set_var (h, n, vmax, &perm);
      b->h = hgraph_copy_internal (pr, h);
      ushape_apply_dimperm (pr, b, b, &perm);
      ap_dimperm_clear (&perm);
    }
  else
    b->h = hgraph_copy_internal (pr, h);
  hgraph_free_internal (pr, h);
  return b;
}

ushape_t *
ushape_add_set_info (ushape_internal_t * pr, ushape_t * a,
                     node_t nsucc, size_t v, size_t nn)
{
  node_t i;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->h->size);
  ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  hgraph_t *h = hgraph_node_add (pr, a->h, nsucc, v, nn, &i, &perm);
  b->h = hgraph_copy_internal (pr, h);
  hgraph_free_internal (pr, h);
  ushape_apply_dimperm_dimchange (pr, a, b, &perm, i);
  return b;
}

ushape_t *
ushape_copy_set_info (ushape_internal_t * pr, ushape_t * a, node_t n,
                      node_t nsucc, size_t v, size_t nn)
{
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->h->size);
  ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  hgraph_t *h = hgraph_copy_mem (pr, a->h);
  hgraph_node_set_var (h, n, v, &perm);
  hgraph_node_set_succ (h, n, nsucc, &perm);
  b->h = hgraph_copy_internal (pr, h);
  hgraph_free_internal (pr, h);
  ushape_apply_dimperm (pr, a, b, &perm);
  return b;
}

ushape_t *
ushape_set_new_succ (ushape_internal_t * pr, ushape_t * a, node_t n, size_t v)
{
  node_t nnew;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->h->size);
  ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  hgraph_t *h = hgraph_node_add (pr, a->h, NODE_T_TOP, v, 0, &nnew, &perm);
  ushape_apply_dimperm_dimchange (pr, a, b, &perm, nnew); /* also clear perm */
  hgraph_node_set_succ (h, n, nnew, &perm);
  b->h = hgraph_copy_internal (pr, h);
  hgraph_free_internal (pr, h);
  ushape_apply_dimperm (pr, b, b, &perm);
  return b;
}

size_t
ushape_set_new_succ_all (ushape_internal_t * pr, ushape_t * a, size_t v1,
                         size_t v2, ushape_array_t * r, size_t rsize)
{
  node_t n1, n2;
  size_t i;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->h->size);
  ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  hgraph_t *h = hgraph_node_add (pr, a->h, NODE_T_TOP, v2, 0, &n2, &perm);
  b->h = hgraph_copy_internal (pr, h);
  ushape_apply_dimperm_dimchange (pr, a, b, &perm, n2);
  i = 0;

  hgraph_node_forall (h, n1)
  {
    if (n1 != n2 && NODE_T_TOP == hgraph_node_get_succ (h, n1))
      {
        ushape_t *b1 = ushape_alloc_internal (pr, b->datadim, b->ptrdim);
        hgraph_t *h1 = hgraph_copy_mem (pr, h);
        ap_dimperm_init (&perm, h->size);
        hgraph_node_set_var (h1, n1, v1, &perm);
        hgraph_node_set_succ (h1, n1, n2, &perm);
        b1->h = hgraph_copy_internal (pr, h1);
        ushape_apply_dimperm (pr, b, b1, &perm);
        i += ushape_array_add (pr, true, r, rsize + i, true, true, b1); /* copy and destroy */
        hgraph_free_internal (pr, h1);
      }
  }
  hgraph_free_internal (pr, h);
  ushape_free_internal (pr, b);
  return i;
}

ushape_t *
ushape_add_loop (ushape_internal_t * pr, ushape_t * a, size_t vmin,
                 size_t vmax)
{
  node_t n;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->h->size);
  ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  hgraph_t *h = hgraph_node_add (pr, a->h, NODE_T_TOP, vmin,
                                 0, &n, &perm);
  ushape_apply_dimperm_dimchange (pr, a, b, &perm, n);
  ap_dimperm_init (&perm, h->size);
  if (vmax && vmin != vmax)
    hgraph_node_set_var (h, n, vmax, &perm);
  hgraph_node_set_succ (h, n, n, &perm);
  b->h = hgraph_copy_internal (pr, h);
  hgraph_free_internal (pr, h);
  ushape_apply_dimperm (pr, b, b, &perm);
  return b;
}

ushape_t *
ushape_add_edge (ushape_internal_t * pr, ushape_t * a, size_t vsrc,
                 size_t vdst)
{
  node_t nsrc, ndst;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->h->size);
  ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  hgraph_t *h = hgraph_node_add (pr, a->h, NODE_T_TOP, vdst, 0, &ndst, &perm);
  ushape_apply_dimperm_dimchange (pr, a, b, &perm, ndst);
  hgraph_t *h1 = h;
  ap_dimperm_init (&perm, h->size);
  h = hgraph_node_add (pr, h1, ndst, vsrc, 0, &nsrc, &perm);
  hgraph_free_internal (pr, h1);
  b->h = hgraph_copy_internal (pr, h);
  hgraph_free_internal (pr, h);
  ushape_apply_dimperm_dimchange (pr, b, b, &perm, nsrc);
  return b;
}

ushape_t *
ushape_add_between (ushape_internal_t * pr, ushape_t * a, size_t vsrc,
                    node_t nsrc, node_t ndst)
{
  node_t i;
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, a->h->size);
  ap_dimperm_set_id (&perm);
  ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  hgraph_t *h = hgraph_node_add (pr, a->h, ndst, vsrc, 1, &i, &perm);
  ushape_apply_dimperm_dimchange (pr, a, b, &perm, i);
  ap_dimperm_init (&perm, h->size);
  ap_dimperm_set_id (&perm);
  hgraph_node_set_succ (h, nsrc, i, &perm);
  b->h = hgraph_copy_internal (pr, h);
  hgraph_free_internal (pr, h);
  ushape_apply_dimperm (pr, b, b, &perm);
  return b;
}

/* Meet with a data constraint */
ushape_array_t*
ushape_meet_pcons_data (ushape_internal_t* pr, ushape_t* a, pcons0_t* c, size_t* rsize)
{
  ushape_array_t* r = NULL;
  bool isbot;
  size_t i;
  ap_lincons0_array_t arr = ap_lincons0_array_make (1);
  size_t *v2n = hgraph_get_var2node (a->h);
  arr.p[0] = // set kind (scalar) to data constraint
          shape_lincons_of_node (pr, &c->info.data.cons, c->info.data.offsets, v2n, a->h->size,
                                 a->datadim, a->ptrdim);
  free (v2n);
  r = ushape_array_make (pr, 1);
  r->p[0] = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
  r->p[0]->h = hgraph_copy_internal (pr, a->h);
  r->p[0]->closed = hgraph_copy_internal (pr, a->closed);
  *rsize = 1;
  isbot = false;
  for (i = 0; i < pr->size_scons && !isbot; i++)
    {
      r->p[0]->scons[i] =
              ap_abstract0_meet_lincons_array (pr->man_scons[i], false,
                                               a->scons[i], &arr);
      isbot = ap_abstract0_is_bottom (pr->man_scons[i], a->scons[i]);
    }
  ap_lincons0_array_clear (&arr);
  if (isbot)
    {
      ushape_array_clear (pr, r, r->size);
      free (r);
      r = NULL;
      *rsize = 0;
    }
  return r;
}

/* Test or (if isassume) build an acyclic graph from vx */
ushape_array_t*
ushape_meet_pcons_ptr_acyclic (ushape_internal_t* pr,
                               ushape_t* a, pcons0_t *c,
                               node_t vx, node_t nx, node_t nnx,
                               node_t vy, node_t ny, node_t nny,
                               bool isassume,
                               size_t* rsize)
{
  ushape_array_t* r = ushape_array_make (pr, 4);
  *rsize = 0;
  if (isassume)
    {
      /* build acyclic graph */
      node_t n;
      size_t i;
      ap_dimperm_t perm;
      ap_dimperm_init (&perm, a->h->size);
      // Step 1: add to a->h a node labeled by vx, get the node
      ushape_t *b = ushape_alloc_internal (pr, a->datadim, a->ptrdim);
      hgraph_t *h = hgraph_node_add (pr, a->h, NODE_NULL, vx, 0, &n, &perm);
      b->h = hgraph_copy_internal (pr, h);
      ushape_apply_dimperm_dimchange (pr, a, b, &perm, n);
      hgraph_free_internal (pr, h);
      // Step 2: add the constraint l[n] >= 1
      ap_lincons0_array_t arr = ap_lincons0_array_make (1);
      ap_linexpr0_t *lexpr =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                 a->datadim + a->h->size);
      ap_scalar_t * zero = ap_scalar_alloc ();
      ap_scalar_set_int (zero, OFFSET_LEN);
      ap_linexpr0_set_coeff_scalar_int (lexpr,
                                        b->datadim +
                                        DIM2NODE (b->h, vx),
                                        1);
      ap_linexpr0_set_cst_scalar_int (lexpr, -1);
      arr.p[0] = ap_lincons0_make (AP_CONS_SUPEQ, lexpr, zero);
      for (i = 0; i < pr->size_scons; i++)
        b->scons[i] =
              ap_abstract0_meet_lincons_array (pr->man_scons[i],
                                               true,
                                               b->scons[i],
                                               &arr);
      // Step 3: add to the result and free the constraints allocated
      *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */
      ap_lincons0_array_clear (&arr);
      ap_scalar_free (zero);
    }
  else
    {
      /* do meet ==> bottom or same */
      hgraph_array_t* h_arr = hgraph_meet_pcons (pr, false, a->h, c);
      if (h_arr != NULL && !h_arr->size == 0)
        *rsize += ushape_array_add (pr, true, r, *rsize, true, false, a); /* copy and destroy */
    }
  return r;
}

/* Test or (if isassume) build an acyclic graph from nx */
ushape_array_t*
ushape_meet_pcons_ptr_cyclic (ushape_internal_t* pr,
                              ushape_t* a, pcons0_t *c,
                              node_t vx, node_t nx, node_t nnx,
                              node_t vy, node_t ny, node_t nny,
                              bool isassume,
                              size_t* rsize)
{
  ushape_array_t* r = ushape_array_make (pr, 4);
  *rsize = 0;
  if (isassume)
    {
      /* build cyclic graph */
      node_t n;
      size_t i;
      // Step 1: build the graph
      ushape_t* b = ushape_add_loop (pr, a, vx, vx); // permutation of nodes done inside!
      // Step 2: add the length constraint l[nx]>=1
      ap_lincons0_array_t arr = ap_lincons0_array_make (1);
      ap_linexpr0_t *lexpr =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                 a->datadim + a->h->size);
      ap_scalar_t * zero = ap_scalar_alloc ();
      ap_scalar_set_int (zero, OFFSET_LEN);
      ap_linexpr0_set_coeff_scalar_int (lexpr,
                                        b->datadim +
                                        DIM2NODE (b->h, vx),
                                        1);
      ap_linexpr0_set_cst_scalar_int (lexpr, -1);
      arr.p[0] = ap_lincons0_make (AP_CONS_SUPEQ, lexpr, zero);
      for (i = 0; i < pr->size_scons; i++)
        b->scons[i] =
              ap_abstract0_meet_lincons_array (pr->man_scons[i],
                                               true,
                                               b->scons[i],
                                               &arr);
      *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b);
      ap_lincons0_array_clear (&arr);
      ap_scalar_free (zero);
    }
  else
    {
      /* do meet ==> bottom or same */
      hgraph_array_t* h_arr = hgraph_meet_pcons (pr, false, a->h, c);
      if (h_arr != NULL && h_arr->size != 0)
        *rsize += ushape_array_add (pr, true, r, *rsize, true, false, a); /* copy and destroy */
    }
  return r;
}

/* Meet with iso graphs from @code{vx} and @code{vy}.
 * @code{isassume} is not used.
 */
ushape_array_t*
ushape_meet_pcons_ptr_iso (ushape_internal_t* pr,
                           ushape_t* a, pcons0_t *c,
                           node_t vx, node_t nx, node_t nnx,
                           node_t vy, node_t ny, node_t nny,
                           bool isassume,
                           size_t* rsize)
{
  ushape_array_t* r = ushape_array_make (pr, 4);
  *rsize = 0;
  size_t i;
  // test sub-graphs iso in a->h ==> bottom or the same
  hgraph_array_t* h_arr = hgraph_meet_pcons (pr, false, a->h, c);
  if (h_arr != NULL && h_arr->size != 0)
    {
      // Subgraphs are isomorphic, put constraint on data words
      ushape_t* b = ushape_copy_internal (pr, a);
      // Meet with equality constraints in each data domain
      ap_lincons0_array_t arr = ap_lincons0_array_make (1);
      ap_linexpr0_t *lexpr =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                 a->datadim + a->h->size); // Warning: why 2*a->h->size?
      arr.p[0] = ap_lincons0_make (AP_CONS_EQMOD, lexpr, NULL);
      ap_linexpr0_set_coeff_scalar_int (lexpr,
                                        b->datadim +
                                        nx,
                                        1);
      ap_linexpr0_set_coeff_scalar_int (lexpr,
                                        b->datadim +
                                        ny,
                                        -1);
      // ap_linexpr0_set_cst_scalar_int(lexpr, -1);
      for (i = 0; i < pr->size_scons; i++)
        b->scons[i] =
              ap_abstract0_meet_lincons_array (pr->man_scons[i],
                                               true,
                                               b->scons[i],
                                               &arr);
      *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */
      ap_lincons0_array_clear (&arr);
    }
  return r;
}

/* Meet with or (if isassume) build equality constraint between pointers */
ushape_array_t*
ushape_meet_pcons_ptr_eq (ushape_internal_t* pr,
                          ushape_t* a, pcons0_t *c,
                          node_t vx, node_t nx, node_t nnx,
                          node_t vy, node_t ny, node_t nny,
                          bool isassume,
                          size_t* rsize)
{
  ushape_t* b;
  ushape_array_t* r = ushape_array_make (pr, 4);
  *rsize = 0;
#ifndef NDEBUG1
  fprintf (stdout, "\n++++ushape_meet_pcons_ptr_eq: nnx=%zu, nny=%zu\n", nnx, nny);
#endif
  if (nnx != NODE_T_TOP && nny != NODE_T_TOP)
    {
      /* both nodes are defined, check their equality */
      if (nnx == nny)
        *rsize += ushape_array_add (pr, true, r, *rsize, true, false, a); /* copy, not distroy */
      /* else bottom result, keep NULL */
    }
  else
    if ((nx != NODE_T_TOP && nnx == NODE_T_TOP
         && nny != NODE_T_TOP) || (ny != NODE_T_TOP
                                   && nny == NODE_T_TOP
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
      b = ushape_copy_set_succ (pr, a, n1, n2);
      *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */
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
          b = ushape_add_succ_fixed (pr, a, vany, nfixed);
          *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */

          /* case (2) vany->next = vfixed, existing node for vany */
          hgraph_node_forall (a->h, n)
          {
            if (hgraph_node_get_succ (a->h, n) == nfixed)
              {
                b = ushape_copy_set_var (pr, a, n, vany, 0, 0);
                *rsize += ushape_array_add (pr, *rsize, r, true, true, true, b); /* copy and destroy */
              }
            else if (hgraph_node_get_succ (a->h, n) ==
                     NODE_T_TOP)
              {
                b =
                        ushape_copy_set_info (pr, a, n, nfixed, vany,
                                              0);
                *rsize += ushape_array_add (pr, *rsize, r, true, true, true, b); /* copy and destroy */
              }
            else
              continue;
          }
        }
      else if (nnfixed != NODE_T_TOP) /* nextany == 0 */
        { /* vfixed[->next] = vany */
          /* set vany to nnfixed */
          b = ushape_copy_set_var (pr, a, nnfixed, vany, 0, 0);
          *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */
        }
      else
        { /* nextany == 0 && nnfixed == NODE_T_TOP */
          /* case (1): vfixed->next = vany, new node for vany */
          node_t n;
          b = ushape_set_new_succ (pr, a, nfixed, vany);
          *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */

          /* case (2) vfixed->next = vany, existing node for vany */
          hgraph_node_forall (a->h, n)
          {
            b =
                    ushape_copy_set_info (pr, a, n, nfixed, vany, 0);
            *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */
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
          b = ushape_add_set_var (pr, a, NODE_T_TOP, v1, v2);
          *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */

          /* case (2): set x and y to existing nodes */
          hgraph_node_forall (a->h, n)
          {
            b = ushape_copy_set_var (pr, a, n, v1, n, v2);
            *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */
          }
        }
      else
        {
          /*
           * nx == NODE_T_TOP && ny == NODE_T_TOP && there is a next
           */
          /* v1->next == v2 */
          size_t v1, v2;
          node_t n, n2;
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
          b = ushape_add_loop (pr, a, ((v1 <= v2) ? v1 : v2),
                               ((v1 < v2) ? v2 : v1));
          *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */
          /* case (2): v1 != v2 and both label new nodes */
          b = ushape_add_edge (pr, a, v1, v2);
          *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */

          /* case (3): v1 and v2 are existing nodes */
          hgraph_node_forall (a->h, n)
          {
            if (n != hgraph_node_get_null (a->h))
              {

                hgraph_node_forall (a->h, n2)
                {
                  if (n2 == hgraph_node_get_succ (a->h, n))
                    {
                      b =
                              ushape_copy_set_var (pr, a, n, v1, n2,
                                                   v2);
                      *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */
                    }
                }
              }
          }
          /* case (4): v1 is new, v2 is an existing node */
          *rsize +=
                  ushape_add_succ_fixed_all (pr, a, v1, v2, r, *rsize);
          /* case (5): v1 is old with succ undefined, v2 is new */
          ushape_set_new_succ_all (pr, a, v1, v2, r, *rsize);
        }
    }
  return r;
}

ushape_array_t*
ushape_meet_pcons_ptr_neq (ushape_internal_t* pr,
                           ushape_t* a, pcons0_t *c,
                           node_t vx, node_t nx, node_t nnx,
                           node_t vy, node_t ny, node_t nny,
                           bool isassume,
                           size_t* rsize)
{
  ushape_t* b;
  node_t i;
  ushape_array_t* r = ushape_array_make (pr, 4);
  *rsize = 0;
#ifndef NDEBUG1
  fprintf (stdout, "\n++++ushape_meet_pcons_ptr_neq: nnx=%zu, nny=%zu\n", nnx, nny);
#endif
  if (nnx != NODE_T_TOP && nny != NODE_T_TOP)
    {
      /* check only equality of nodes */
      if (nnx != nny)
        *rsize +=
              ushape_array_add (pr, true, r, *rsize, true, false, a);
      else if ((c->info.ptr.offx != OFFSET_NONE &&
                c->info.ptr.offy == OFFSET_NONE) ||
               (c->info.ptr.offx == OFFSET_NONE &&
                c->info.ptr.offy != OFFSET_NONE))
        { /* v1->next != v2 : create a new node s.t. v1
                   * -->new node --> v2 */
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
          b = ushape_add_between (pr, a, v1, n1, n2);
          *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */
        }
      /* else bottom result, so keep NULL */
    }
  else
    if ((nx != NODE_T_TOP && nnx == NODE_T_TOP
         && nny != NODE_T_TOP) || (ny != NODE_T_TOP
                                   && nny == NODE_T_TOP
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
      b = ushape_add_set_info (pr, a, NODE_T_TOP, v1, 1);
      *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */

      /* case (2): fix v1->next to some existing node, different from n2 */
      hgraph_node_forall (a->h, i)
      {
        if (n2 != hgraph_node_get_succ (a->h, i))
          {
            b = ushape_copy_set_var (pr, a, i, v1, 0, 0);
            *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */
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
      b = ushape_add_set_info (pr, a, NODE_T_TOP, v1, 0);
      *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */

      /* case (2): fix v1 to some existing node, different from n2 */
      hgraph_node_forall (a->h, i)
      {
        if (n2 != i)
          {
            b = ushape_copy_set_var (pr, a, i, v1, 0, 0);
            *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */
          }
      }
    }
  else
    { /* TODO: the same cases like for equality,
                   * but more combinatorial explosion. */
      ap_manager_raise_exception (pr->man,
                                  AP_EXC_NOT_IMPLEMENTED,
                                  pr->funid, "not implemented");
      return NULL;
    }
  return r;
}

ushape_array_t*
ushape_meet_pcons_ptr_reach (ushape_internal_t* pr,
                             ushape_t* a, pcons0_t *c,
                             node_t vx, node_t nx, node_t nnx,
                             node_t vy, node_t ny, node_t nny,
                             bool isassume,
                             size_t* rsize)
{
  ushape_t* b;
  node_t i;
  ushape_array_t* r = ushape_array_make (pr, 4);
  *rsize = 0;
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
      bool found = hgraph_node_is_reachable (a->h, nx, ny);
      if (found && c->type == REACH_CONS)
        *rsize +=
              ushape_array_add (pr, true, r, *rsize, true, false, a);
      /* else bottom result, keep NULL entry */
    }
  else if (nx != NODE_T_TOP && ny == NODE_T_TOP)
    {
      /* set y to some next successor of nx */
      /* TODO: the test is more complex because of lasso from x */
      int j = 0;
      ap_lincons0_array_t arr = ap_lincons0_array_make (1);
      ap_linexpr0_t *lexpr =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                 a->datadim + 2 * a->h->size);
      ap_linexpr0_set_coeff_scalar_int (lexpr,
                                        a->datadim +
                                        a->h->size + nx, 1);
      i = NODE_T_TOP;
      ny = hgraph_node_get_succ (a->h, nx);
      while (ny != NODE_T_TOP && ny != i && j <= 3)
        {
          b = ushape_copy_set_var (pr, a, ny, vy, 0, 0);
          // add constraint l[nx]>=1
          arr.p[0] = ap_lincons0_make (AP_CONS_SUPEQ, lexpr, NULL);
          ap_linexpr0_set_cst_scalar_int (lexpr, -1);
          for (i = 0; i < pr->size_scons; i++)
            b->scons[i] =
                  ap_abstract0_meet_lincons_array (pr->
                                                   man_scons[i],
                                                   true,
                                                   b->scons[i],
                                                   &arr);

          *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); /* copy and destroy */
          i = ny;
          j++;
          ap_linexpr0_set_coeff_scalar_int (lexpr,
                                            a->datadim +
                                            a->h->size + i, 1);
          ny = hgraph_node_get_succ (a->h, i);
        }
    }
  else if (nx == NODE_T_TOP && ny != NODE_T_TOP)
    {
      /* case (1): set x to some new node with the next successor at ny */
      /* TODO: the case is more complex since several nodes may exist to ny */
      b = ushape_add_set_info (pr, a, ny, vx, 0);
      // add constraint l[nx]>=1
      ap_lincons0_array_t arr = ap_lincons0_array_make (1);
      arr.p[0] =
              shape_lincons_x_y_v_cst (AP_CONS_SUPEQ, OFFSET_LEN,
                                       1,
                                       b->datadim +
                                       DIM2NODE (b->h, vx), 0, 0,
                                       0, 0, -1,
                                       b->datadim, b->h->size);
      for (i = 0; i < pr->size_scons; i++)
        b->scons[i] =
              ap_abstract0_meet_lincons_array (pr->man_scons[i],
                                               true,
                                               b->scons[i],
                                               &arr);
      ap_lincons0_array_clear (&arr);
      *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); // copy and destroy

      /* case (2): if ny is not NULL, */
      /* set x to some existing node predecessor of ny */
      /* TODO: the case is more complex, we choose only immediate predecessor of ny */
      if (ny != NODE_NULL)
        {

          hgraph_node_forall (a->h, i)
          {
            if (ny == hgraph_node_get_succ (a->h, i))
              {
                b = ushape_copy_set_var (pr, a, i, vx, 0, 0);
                // add constraint l[nx]>=1
                ap_lincons0_array_t arr =
                        ap_lincons0_array_make (1);
                arr.p[0] =
                        shape_lincons_x_y_v_cst (AP_CONS_SUPEQ,
                                                 OFFSET_LEN, 1,
                                                 b->datadim +
                                                 DIM2NODE (b->h, vx),
                                                 0, 0, 0,
                                                 0, -1,
                                                 b->datadim,
                                                 b->h->size);
                for (i = 0; i < pr->size_scons; i++)
                  b->scons[i] =
                        ap_abstract0_meet_lincons_array (pr->
                                                         man_scons
                                                         [i], true,
                                                         b->
                                                         scons[i],
                                                         &arr);
                ap_lincons0_array_clear (&arr);
                *rsize += ushape_array_add (pr, true, r, *rsize, true, true, b); //copy and destroy

              }
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
      ap_lincons0_array_t arr = ap_lincons0_array_make (1);
      ap_linexpr0_t *lexpr =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                 a->datadim + 2 * a->h->size);
      /* case (1): x and y are the same new node */
      b = ushape_add_loop (pr, a, v1, v2); /* v1 <= v2 */
      // add constraint l[nx]>=1
      arr.p[0] = ap_lincons0_make (AP_CONS_SUPEQ, lexpr, NULL);
      ap_linexpr0_set_coeff_scalar_int (lexpr,
                                        b->datadim +
                                        b->h->size +
                                        DIM2NODE (b->h, vx),
                                        1);
      ap_linexpr0_set_cst_scalar_int (lexpr, -1);
      for (i = 0; i < pr->size_scons; i++)
        b->scons[i] =
              ap_abstract0_meet_lincons_array (pr->man_scons[i],
                                               true,
                                               b->scons[i],
                                               &arr);

      *rsize += ushape_array_add (pr, true, r, *rsize, true, false, b);

      ushape_free_internal (pr, b);
      /* case (2): x and y are different new nodes */
      b = ushape_add_edge (pr, a, vx, vy);
      // add constraint l[nx]>=1
      for (i = 0; i < pr->size_scons; i++)
        b->scons[i] =
              ap_abstract0_meet_lincons_array (pr->man_scons[i],
                                               true,
                                               b->scons[i],
                                               &arr);

      *rsize += ushape_array_add (pr, true, r, *rsize, true, false, b); /* copy and destroy */
      // add constraint l[ny]>=1
      ap_linexpr0_set_coeff_scalar_int (lexpr,
                                        b->datadim +
                                        b->h->size +
                                        DIM2NODE (b->h, vx),
                                        0);
      ap_linexpr0_set_coeff_scalar_int (lexpr,
                                        b->datadim +
                                        b->h->size +
                                        DIM2NODE (b->h, vy),
                                        1);

      for (i = 0; i < pr->size_scons; i++)
        b->scons[i] =
              ap_abstract0_meet_lincons_array (pr->man_scons[i],
                                               true,
                                               b->scons[i],
                                               &arr);

      *rsize += ushape_array_add (pr, true, r, *rsize, true, false, b); /* copy and destroy */

      ap_linexpr0_set_coeff_scalar_int (lexpr,
                                        b->datadim +
                                        b->h->size +
                                        DIM2NODE (b->h, vy),
                                        0);

      ushape_free_internal (pr, b);
      /* case (3): x new, y existing */
      /* TODO: add length constraints */
      *rsize +=
              ushape_add_succ_fixed_all (pr, a, vx, vy, r, *rsize);
      /* case (4): x existing but no succ, y new */
      /* TODO: simple case here, no intermediate node between x and y */
      /* TODO: add length constraints */
      ushape_set_new_succ_all (pr, a, vx, vy, r, *rsize);
      /* case (5): both x and y are existing nodes */

      /* TODO: more complex cases with more nodes between x and y */
      hgraph_node_forall (a->h, nx)
      {
        if (nx != hgraph_node_get_null (a->h))
          {

            hgraph_node_forall (a->h, ny)
            {
              if (ny == hgraph_node_get_succ (a->h, nx))
                {
                  b = ushape_copy_set_var (pr, a, nx, vx, ny, vy);
                  // add constraint l[nx]>=1
                  ap_linexpr0_set_coeff_scalar_int (lexpr,
                                                    b->
                                                    datadim
                                                    +
                                                    b->h->
                                                    size +
                                                    DIM2NODE
                                                    (b->h,
                                                     vx),
                                                    1);
                  for (i = 0; i < pr->size_scons; i++)
                    b->scons[i] =
                          ap_abstract0_meet_lincons_array (pr->
                                                           man_scons
                                                           [i],
                                                           true,
                                                           b->
                                                           scons
                                                           [i],
                                                           &arr);
                  *rsize += ushape_array_add (pr, true, r, *rsize, true, false, b); /* copy and destroy */

                  ap_linexpr0_set_coeff_scalar_int (lexpr,
                                                    b->
                                                    datadim
                                                    +
                                                    b->h->
                                                    size +
                                                    DIM2NODE
                                                    (b->h,
                                                     vx),
                                                    0);
                  ushape_free_internal (pr, b);
                }
            }
          }
      }
    }
  return r;
}

ushape_array_t *
ushape_meet_pcons (ushape_internal_t * pr, bool destructive,
                   ushape_t * a, pcons0_t * c, bool isassume)
{
  ushape_array_t *r = NULL;
  size_t rsize = 0;

  arg_assert (a && c, return NULL;);

#ifndef NDEBUG
  fprintf (stdout, "\n++++ushape_meet_pcons: with pcons=(");
  shape_pcons_fdump (stdout, c);
  fprintf (stdout, ")\n");
  ushape_fdump (stdout, pr->man, a);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif

  if (!ushape_is_bottom (pr->man, a))
    {
      if (c->type == SL3_CONS)
        {
          ERROR ("Ushape_meet_pcons: bad constraint (SL3)!",;
                 );
          return NULL;
        }
      else if (c->type == DATA_CONS)
        r = ushape_meet_pcons_data (pr, a, c, &rsize);
      else if (c->type < OTHER_CONS)
        {
          node_t vx, nx, nnx, vy, ny, nny;
          /* OLD: no more than 1 next dereference */
          // if (c->info.ptr.nextx + c->info.ptr.nexty > 1) {
          /* NEW: no next derefencing in constraints */
          if (c->info.ptr.offx != OFFSET_NONE ||
              c->info.ptr.offy != OFFSET_NONE)
            {
              ERROR ("Bad constraint: next reference is not allowed!",;
                     );
              //ap_manager_raise_exception(pr->man, AP_EXC_NOT_IMPLEMENTED,
              //        pr->funid, "not implemented");
              return NULL;
            }
          vx = DIM2PTR (c->info.ptr.x, a->datadim);
          nx = DIM2NODE (a->h, vx);
          // OLD nnx = hgraph_node_get_nsucc(a->h, nx, c->info.ptr.nextx);
          nnx = nx;
          vy = DIM2PTR (c->info.ptr.y, a->datadim);
          ny = DIM2NODE (a->h, vy);
          // OLD nny = hgraph_node_get_nsucc(a->h, ny, c->info.ptr.nexty);
          nny = ny;
          arg_assert (vx < a->ptrdim || IS_NULLDIM (vx), return NULL;
                      );

          if (vx == 0 && vy == 0)
            {
              // empty constraint, return a
              rsize += ushape_array_add (pr, true, r, rsize, true, false, a); /* copy and destroy */
            }
          else /* if (vx < a->datadim && vy < a->datadim) {
                ap_manager_raise_exception(pr->man, AP_EXC_NOT_IMPLEMENTED,
                        pr->funid, "not implemented");
                return NULL;
            } else */
            {
              arg_assert (vy < a->ptrdim || IS_NULLDIM (vy), return NULL;
                          );

              if ((nx == 0 && c->info.ptr.offx != OFFSET_NONE)
                  || (ny == 0 && c->info.ptr.offy != OFFSET_NONE))
                {
                  ERROR ("Null pointer dereference",;
                         );
                  return NULL;
                }
              switch (c->type)
                {
                case ACYCLIC_CONS:
                  r = ushape_meet_pcons_ptr_acyclic (pr, a, c, vx, nx, nnx, vy, ny, nny,
                                                     isassume, &rsize);
                  break;
                case CYCLIC_CONS:
                  r = ushape_meet_pcons_ptr_cyclic (pr, a, c, vx, nx, nnx, vy, ny, nny,
                                                    isassume, &rsize);
                  break;
                case ISO_CONS:
                  r = ushape_meet_pcons_ptr_iso (pr, a, c, vx, nx, nnx, vy, ny, nny,
                                                 isassume, &rsize);
                  break;
                case EQ_CONS:
                  r = ushape_meet_pcons_ptr_eq (pr, a, c, vx, nx, nnx, vy, ny, nny,
                                                isassume, &rsize);
                  break;
                case SEP_CONS:
                case NE_CONS:
                  r = ushape_meet_pcons_ptr_neq (pr, a, c, vx, nx, nnx, vy, ny, nny,
                                                 isassume, &rsize);
                  break;
                case REACH_CONS:
                  r = ushape_meet_pcons_ptr_reach (pr, a, c, vx, nx, nnx, vy, ny, nny,
                                                   isassume, &rsize);
                  break;

                default:
                  break;
                }
            }
        }
    }
  if (destructive)
    ushape_free_internal (pr, a);
  if (rsize)
    ushape_array_resize (pr, r, rsize);
  else if (r)
    {
      ushape_array_clear (pr, r, r->size);
      free (r);
      r = NULL;
    }

#ifndef NDEBUG
  fprintf (stdout, "\n++++ushape_meet_pcons returns: ");
  if (!r)
    fprintf (stdout, "null\n");
  else
    {
      size_t i;
      fprintf (stdout, "array of size %zu and elements [\n", r->size);
      for (i = 0; i < r->size; i++)
        if (r->p[i])
          ushape_fdump (stdout, pr->man, r->p[i]);
        else
          fprintf (stdout, "NULL\n");
      fprintf (stdout, "]\n");
    }
#endif
  return r;
}

ushape_array_t *
ushape_meet_pcons_array (ushape_internal_t * pr,
                         bool destructive, ushape_t * a,
                         pcons0_array_t * array)
{
  ushape_array_t *r;
  size_t i, j;
  bool isassume = ushape_is_top (pr->man, a);
  r = ushape_array_make (pr, 4);
  r->p[0] = ushape_copy_internal (pr, a);
  for (i = 0; i < array->size && r; i++)
    {
      //TODO:apply the constraints but not iteratively
      ushape_array_t *rc = NULL;
      for (j = 0; j < r->size && r->p[j] != NULL; j++)
        {
          ushape_array_t *lr =
                  ushape_meet_pcons (pr, false, r->p[j], array->p[i], isassume);
          if (lr)
            {
              rc = ushape_array_add_array (pr, true, rc, lr); /* reused elements of lr */
              ushape_array_init (pr, lr, lr->size);
              free (lr);
            }
        }
      ushape_array_clear (pr, r, r->size);
      r = rc;
    }
  if (isassume)
    // put to NULL all unconstrained ptr vars in r
    for (j = 0; j < r->size; j++)
      ushape_canonicalize (pr->man, r->p[j]);

  if (destructive)
    ushape_free_internal (pr, a);
  return r;
}

ushape_t *
ushape_meet_lincons_array (ap_manager_t * man,
                           bool destructive, ushape_t * a,
                           ap_lincons0_array_t * array)
{
  if (ushape_is_bottom (man, a) || !array || array->size == 0)
    /* nothing to do */
    return (destructive) ? a : ushape_copy (man, a);
  else
    {
      ushape_t *b;
      pcons0_array_t *arr;
      ushape_array_t *r;
      ushape_internal_t *pr =
              ushape_init_from_manager (man, AP_FUNID_MEET_LINCONS_ARRAY, 0);
      if (!destructive)
        b = ushape_copy_internal (pr, a);
      else
        b = a;
      /* compute in arr the constraints sorted */
      arr =
              shape_pcons_array_of_lincons_array (pr, array, a->datadim, a->ptrdim);
#ifndef NDEBUG
      fprintf (stdout, "\n++++ushape_meet_lincons_array: with pcons_array=(\n");
      shape_pcons_array_fdump (stdout, arr);
      fprintf (stdout, ")\n");
#endif
      /* go */
      /* true below means to compute only one element and return it in b */
      r = ushape_meet_pcons_array (pr, true, b, arr);
      if (r && r->size > 0)
        {
          b = r->p[0];
          r->p[0] = NULL;
          ushape_array_clear (pr, r, r->size);
        }
      else
        b = ushape_bottom (man, a->datadim, a->ptrdim);
      return b;
    }
}

ushape_t *
ushape_meet_tcons_array (ap_manager_t * man,
                         bool destructive, ushape_t * a,
                         ap_tcons0_array_t * array)
{
  if (ushape_is_bottom (man, a) || !array || array->size == 0)
    /* nothing to do */
    return (destructive) ? a : ushape_copy (man, a);
  else
    {
      ushape_t *b;
      pcons0_array_t *arr;
      ushape_array_t *r;
      ushape_internal_t *pr =
              ushape_init_from_manager (man, AP_FUNID_MEET_LINCONS_ARRAY, 0);
      if (!destructive)
        b = ushape_copy_internal (pr, a);
      else
        b = a;
      /* compute in arr the constraints sorted */
      arr =
              shape_pcons_array_of_tcons_array (pr, array, a->datadim, a->ptrdim);
#ifndef NDEBUG
      fprintf (stdout, "\n++++ushape_meet_tcons_array: with pcons_array=[");
      shape_pcons_array_fdump (stdout, arr);
      fprintf (stdout, "]\n");
#endif
      /* go */
      /* true below means to compute only one element and return it in b */
      r = ushape_meet_pcons_array (pr, true, b, arr);
      if (r && r->size > 0)
        {
          b = r->p[0];
          r->p[0] = NULL;
          ushape_array_clear (pr, r, r->size);
        }
      else
        b = ushape_bottom (man, a->datadim, a->ptrdim);
      return b;
    }
}

/* Abstract a conjunction of constraints. Based on meet, like in ap_abstract0. */
ushape_t*
ushape_of_lincons_array (ap_manager_t* man,
                         size_t intdim, size_t realdim,
                         ap_lincons0_array_t* array)
{
  ushape_t* res = ushape_top (man, intdim, realdim);
  res = ushape_meet_lincons_array (man, true, res, array);
  return res;
}

/** Abstract a disjunct of an SL3 formulas.
 *  The disjunct is given by the param f.
 */
ushape_t*
ushape_of_formula_aux (ushape_internal_t* pr,
                       sh_formula_t* f,
                       size_t disj)
{
#ifndef NDEBUG1
  fprintf (stdout, "==== ushape_of_formula_aux: disjunct %zu\n", disj);
  fflush (stdout);
#endif
  // fix the dimensions of the ushape
  ap_environment_t* env = f->env;
  size_t datadim = env->intdim;
  size_t ptrdim = env->realdim;
  size_t segmdim = f->form[disj]->nodes->realdim + 1; // added nilNode
  ushape_t *b = ushape_alloc_internal (pr, datadim, ptrdim);

  // fix the disjunct in a constraint
  ap_lincons0_array_t arr = ap_lincons0_array_make (1);
  ap_linexpr0_t *lexpr = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                            datadim + segmdim);
  ap_scalar_t * zero = ap_scalar_alloc ();
  ap_scalar_set_int (zero, OFFSET_SL3);
  ap_linexpr0_set_cst_scalar_int (lexpr, disj); // disj in sh_crt
  arr.p[0] = ap_lincons0_make (AP_CONS_SUPEQ, lexpr, zero);

  // Step 1: build the graph part 
  ap_dimperm_t perm;
  ap_dimperm_init (&perm, segmdim); // NULL already added
  ap_dimperm_set_id (&perm);
  hgraph_t *h = hgraph_of_formula (pr, f, disj, &perm);
  b->h = hgraph_copy_internal (pr, h);
  hgraph_free_internal (pr, h);

  // Step 2: add the data constraints taking into account node permutation
  ap_dimperm_t dimperm;
  ap_dimperm_init (&dimperm, datadim + segmdim);
  ap_dimperm_set_id (&dimperm);
  shape_dimperm_copy (&dimperm, datadim, &perm); // perm = nodes -> new nodes
  size_t i;
  for (i = 0; i < pr->size_scons; i++)
    {
      // meet with top, no need to apply the reverse permutation 
      // build the data constraint from the shadform over nodes
      ap_abstract0_t* top_i = ap_abstract0_top (pr->man_scons[i], datadim, segmdim);
      b->scons[i] = ap_abstract0_meet_lincons_array (pr->man_scons[i], true,
                                                     top_i,
                                                     &arr);
      // apply perm extended to data
      b->scons[i] = ap_abstract0_permute_dimensions (pr->man_scons[i],
                                                     true,
                                                     b->scons[i], &dimperm);
    }
  // free allocated data
  ap_dimperm_clear (&perm);
  ap_dimperm_clear (&dimperm);
  ap_lincons0_array_clear (&arr); // free also zero!
#ifndef NDEBUG1
  fprintf (stdout, "==== ushape_of_formula_aux: returns r=(");
  ushape_fdump (stdout, pr->man, b);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif
  return b;
}

ushape_array_t*
ushape_of_formula (ushape_internal_t* pr,
                   sh_formula_t* f,
                   size_t* rsize)
{
  ushape_array_t *r = NULL;
  size_t sz = 0;
  /* The queue of formulas is considered in order.
   * No need to keep track of dealt formulas.
   */
  size_t i, lst;
  lst = f->size;
  *rsize = 0;
  if (lst == 0) return NULL;
  // all formulas shall be in the array, so
  r = ushape_array_make (pr, lst);
  // build the ushape for each shadform
  for (i = 0; i < lst; i++)
    {
      ushape_t *rr = ushape_of_formula_aux (pr, f, i);
      if (rr)
        { // put rr in r (directly, without copy)
          sz += ushape_array_add (pr, true, r, sz, false, false, rr); /* not copy, not distroy */
        }
    }
  if (sz < lst)
    *rsize = ushape_array_resize (pr, r, sz);
  else
    *rsize = sz;
  return r;
}

/* Abstract a conjunction of constraints. Based on meet, like in ap_abstract0. */
ushape_t*
ushape_of_tcons_array (ap_manager_t* man,
                       size_t intdim, size_t realdim,
                       ap_tcons0_array_t* array)
{
  ushape_t* res = ushape_top (man, intdim, realdim);
  res = ushape_meet_tcons_array (man, true, res, array);
  return res;
}

/* NOT IMPLEMENTED */
ushape_t *
ushape_add_ray_array (ap_manager_t * man,
                      bool destructive, ushape_t * a,
                      ap_generator0_array_t * array)
{
  ushape_internal_t *pr =
          ushape_init_from_manager (man, AP_FUNID_ADD_RAY_ARRAY, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a;
}

/* ============================================================ */
/* Widening, Narrowing */

/* ============================================================ */

ushape_t *
ushape_widening (ap_manager_t * man, ushape_t * a1, ushape_t * a2)
{
  ushape_internal_t *pr =
          ushape_init_from_manager (man, AP_FUNID_WIDENING, 0);
  arg_assert (a1->datadim == a2->datadim
              && a1->ptrdim == a2->ptrdim, return NULL;
              );
  // widening is done only between ushape with isomorphic graphs
  if (!hgraph_is_eq (pr->man, a1->h, a2->h))
    return NULL;
  ushape_t *b = ushape_alloc_internal (pr, a1->datadim, a1->ptrdim);
  size_t i;
  bool isbot = false;
  b->h = hgraph_copy_internal (pr, a1->h);
  for (i = 0; i < pr->size_scons && !isbot; i++)
    {
      if (a1->scons[i] != NULL && a2->scons[i] != NULL &&
          ap_abstract0_is_leq (pr->man_scons[i], a1->scons[i], a2->scons[i]))
        /* the test is needed to deal with function calls from specs */
        b->scons[i] =
              ap_abstract0_widening (pr->man_scons[i], a1->scons[i], a2->scons[i]);
      else
        b->scons[i] =
              ap_abstract0_copy (pr->man_scons[i], a2->scons[i]);

      isbot = ap_abstract0_is_bottom (pr->man_scons[i], b->scons[i]);
    }
  if (isbot)
    {
      ushape_free_internal (pr, b);
      b = NULL;
    }
  return b;
}

/* TODO: priority 1 */
ushape_t *
ushape_widening_thresholds (ap_manager_t * man,
                            ushape_t * a1, ushape_t * a2,
                            ap_scalar_t ** array, size_t nb)
{
  ushape_internal_t *pr =
          ushape_init_from_manager (man, AP_FUNID_WIDENING, nb + 1);
  arg_assert (a1->datadim == a2->datadim, return NULL;
              );
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a2;
}

/* NOT IMPLEMENTED */
ushape_t *
ushape_narrowing (ap_manager_t * man, ushape_t * a1, ushape_t * a2)
{
  ushape_internal_t *pr =
          ushape_init_from_manager (man, AP_FUNID_WIDENING, 0);
  arg_assert (a1->datadim == a2->datadim, return NULL;
              );
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return a2;
}
