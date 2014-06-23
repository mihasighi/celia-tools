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
#include "hgraph.h"
#include "hgraph_internal.h"
#include "apron2shape.h"
#include "shape_macros.h"


/* ********************************************************************** */
/* Heap-graphs */
/* ********************************************************************** */

/*
 * We consider adjacency matrices with nodes labeled by variables and 
 * the NULL constant. Node 0 is always NULL node.
 *
 * The matrix is ordered  by the labeling given by variables,
 * with convention that NULL < prg. variables.
 */

/* ============================================================ */
/* Internal management */
/* ============================================================ */

/* ------------------------------------------------------------ */
/* Nodes */
/* ------------------------------------------------------------ */

/*
 * all program ptr variables are mapped to node n,
 * NULL to node 0
 */
void
node_init_v2n (hgraph_internal_t * pr, node_info_t * info, size_t ptrdim,
               node_t n)
{
  size_t i;
  memset (info, 0, (ptrdim) * sizeof (node_info_t));
  for (i = 0; i < ptrdim; i++)
    info[i].s = n;
}

void
node_info_init (hgraph_internal_t * pr, node_info_t * info, size_t ptrdim,
                size_t size)
{
  size_t i;
  for (i = 0; i < size; i++)
    {
      info[ptrdim + i].s = 0;
      info[ptrdim + i].v = 0;
      info[ptrdim + i].nn = 0;
    }
  // Does below really begins at info[ptrdim]?
  // memset (info + ptrdim, 0, size * sizeof (node_info_t));
  info[ptrdim + NODE_NULL].v = NULL_DIM;
}

void
node_info_copy (node_info_t * dst, node_info_t * src)
{
  dst->s = src->s;
  dst->v = src->v;
  dst->nn = src->nn;
}

/*
 * Returns the graph with only null node and ptr variables
 * pointing to null! Not in the Htable.
 */
hgraph_t *
hgraph_alloc_null (hgraph_internal_t * pr, size_t datadim, size_t ptrdim)
{
  hgraph_t *r;
  r = hgraph_alloc_internal (pr, 1, datadim, ptrdim);
  // make ptr vars to point to node null (0)
  memset (r->info, 0, ptrdim * sizeof (node_info_t));
  node_info_init (pr, r->info, ptrdim, 1);
  r->closed = true;
  return r;
}

/*
 * Returns the graph with size nodes and
 * program pointer variables outside the node set.
 * NOT IN THE HTABLE, NOT CANONICAL.
 */
hgraph_t *
hgraph_alloc_internal (hgraph_internal_t * pr, size_t size,
                       size_t datadim, size_t ptrdim)
{
  hgraph_t *r;
  checked_malloc (r, hgraph_t,
                  (sizeof (hgraph_t) +
                   (ptrdim + size) * sizeof (node_info_t)), 1, return NULL;
                  );
  r->closed = false;
  r->size = size;
  r->ptrdim = ptrdim;
  r->datadim = datadim;
  memset (r->info, 0, (ptrdim + size) * sizeof (node_info_t));
  // set ptr vars to point anywhere (NODE_T_TOP)
  node_init_v2n (pr, r->info, ptrdim, NODE_T_TOP);
  // set nodes to undefined variable values except node 0
  // DO NOT COMMENT: node_info_init sets important info for NULL
  node_info_init (pr, r->info, ptrdim, size);
  return r;
}

/* Build graphs used in tests. */
hgraph_t *
hgraph_make (hgraph_internal_t * pr, size_t code, size_t datadim,
             size_t realdim)
{
  size_t ptrdim = REAL2PTR_DIM (pr, realdim);
  hgraph_t *a;
  size_t x = 0;
  size_t y = 2;
  size_t z = 4;
  size_t i;
  switch (code)
    {
    case 0: /* x--> null */
      {
        a = hgraph_alloc_internal (pr, 2, datadim, ptrdim);
        /* set node 1 to x  */
        NODE_VAR (a, 1) = x;
        NODE_NEXT (a, 1) = NODE_NULL;
        NODE_VAR_NEXT (a, 1) = 0;
        VAR2NODE (a, x) = NODE_NULL + 1;
        /* set other vars to null */
        for (i = x + 1; i < ptrdim; i++)
          VAR2NODE (a, i) = NODE_NULL;
        break;
      }
    case 1: /* x-->y-->null and all other ptr vars to
				 * null */
      {
        a = hgraph_alloc_internal (pr, 3, datadim, ptrdim);
        /* set node 1 to x */
        NODE_VAR (a, 1) = x;
        NODE_NEXT (a, 1) = NODE_NULL + 2;
        NODE_VAR_NEXT (a, 1) = 0;
        VAR2NODE (a, x) = NODE_NULL + 1;
        /* set node 2 to y  */
        NODE_VAR (a, 2) = y;
        NODE_NEXT (a, 2) = NODE_NULL;
        NODE_VAR_NEXT (a, 2) = 0;
        VAR2NODE (a, y) = NODE_NULL + 2;
        /* set other vars to null */
        for (i = 0; i < ptrdim; i++)
          if (i != x && i != y)
            VAR2NODE (a, i) = NODE_NULL;
        break;
      }
    case 2: /* x--> null and y --> _null */
      {
        a = hgraph_alloc_internal (pr, 3, datadim, ptrdim);
        /* set node 1 to x  */
        NODE_VAR (a, 1) = x;
        NODE_NEXT (a, 1) = NODE_NULL;
        NODE_VAR_NEXT (a, 1) = 0;
        VAR2NODE (a, x) = NODE_NULL + 1;
        /* set node 2 to y  */
        NODE_VAR (a, 2) = y;
        NODE_NEXT (a, 2) = NODE_NULL;
        NODE_VAR_NEXT (a, 2) = 0;
        VAR2NODE (a, y) = NODE_NULL + 2;
        /* set other vars to null */
        for (i = 0; i < ptrdim; i++)
          if (i != x && i != y)
            VAR2NODE (a, i) = NODE_NULL;
        break;
      }
    case 3: /* x--> null and y --> null and z --> null */
      {
        a = hgraph_alloc_internal (pr, 4, datadim, ptrdim);
        /* set node 1 to x  */
        NODE_VAR (a, 1) = x;
        NODE_NEXT (a, 1) = NODE_NULL;
        NODE_VAR_NEXT (a, 1) = 0;
        VAR2NODE (a, x) = NODE_NULL + 1;
        /* set node 2 to y  */
        NODE_VAR (a, 2) = y;
        NODE_NEXT (a, 2) = NODE_NULL;
        NODE_VAR_NEXT (a, 2) = 0;
        VAR2NODE (a, y) = NODE_NULL + 2;
        /* set node 3 to z  */
        NODE_VAR (a, 3) = z;
        NODE_NEXT (a, 3) = NODE_NULL;
        NODE_VAR_NEXT (a, 3) = 0;
        VAR2NODE (a, z) = NODE_NULL + 3;
        /* set other vars to null */
        for (i = 0; i < ptrdim; i++)
          if (i != x && i != y && i != z)
            VAR2NODE (a, i) = NODE_NULL;
        break;
      }
    case 4: /* x--> null and x==y and all other ptr vars
				 * to null */
      {
        a = hgraph_alloc_internal (pr, 2, datadim, ptrdim);
        /* set node 1 to x  */
        NODE_VAR (a, 1) = x;
        NODE_NEXT (a, 1) = NODE_NULL;
        NODE_VAR_NEXT (a, 1) = 0;
        VAR2NODE (a, x) = NODE_NULL + 1;
        VAR2NODE (a, y) = VAR2NODE (a, x);
        /* set other vars to null */
        for (i = 0; i < ptrdim; i++)
          if (i != x && i != y)
            VAR2NODE (a, i) = NODE_NULL;
        break;
      }
    default:
      a = hgraph_top (pr->man, datadim, ptrdim);
      break;
    }
  if (hgraph_check (pr, a) == '.')
    return hgraph_copy_internal (pr, a);
  else
    return hgraph_top (pr->man, datadim, ptrdim);
}

hgraph_t *
hgraph_random (hgraph_internal_t * pr, size_t size, size_t datadim,
               size_t realdim)
{
  size_t ptrdim = REAL2PTR_DIM (pr, realdim);
  hgraph_t *a = hgraph_alloc_internal (pr, size, datadim, ptrdim);
  size_t i, j, ncut, pi, nn, max_ptr;
  /*
   * select randomly the number of cut nodes, null is
   * already fixed by alloc_internal to 0, program ptr vars are fixed by
   * alloc_internal to size!
   */
  for (j = 0, i = 1, ncut = 1; j < ptrdim;
          j++)
    {
      VAR2NODE (a, j) = NODE_NULL; /* put all prg ptr vars to null */
      if (i < size && lrand48 () % 10 > 1)
        {
          /* map node i to variable j */
          NODE_VAR (a, i) = j;
          NODE_VAR_NEXT (a, i) = 0;
          //a->info[i].s == 0 by alloc
          VAR2NODE (a, j) = i;
          i++;
        }
      ncut = i;
      /* the next nodes are successors of the cut nodes */
      if (i < size)
        {
          size_t ni;
          pi = VAR2NODE (a, j);
          nn = lrand48 () % 4;
          for (ni = 1; ni < nn && i < size; ni++)
            {
              NODE_NEXT (a, pi) = i;
              NODE_VAR (a, i) = j;
              NODE_VAR_NEXT (a, i) = ni;
              pi = i;
              i++;
            }
          NODE_NEXT (a, pi) = 0; /* last node to null */
        }
    }
  /*
   * if still there are nodes not reached, force them to be the successors of
   * the last cut node
   */
  if (i < size)
    {
      int ppi;
      ppi = ncut - 1;
      pi = NODE_NEXT (a, ppi);
      /* last anonymous successor of the last cut node */
      while (NODE_NEXT (a, pi) >= ncut)
        {
          ppi = pi;
          pi = NODE_NEXT (a, pi);
        }
      if (pi < ncut || pi == 0)
        while (i < size)
          {
            NODE_NEXT (a, ppi) = i;
            NODE_VAR (a, i) = NODE_VAR (a, ppi);
            NODE_VAR_NEXT (a, i) = NODE_VAR_NEXT (a, ppi) + 1;
            ppi = i;
            i++;
          }
    }
  /* fill next field of cut nodes not yet filled */
  for (i = 1; i < ncut; i++)
    if (NODE_NEXT (a, i) == 0)
      /* select a cut node radomly */
      NODE_NEXT (a, i) = lrand48 () % ncut;

  hgraph_check (pr, a);
#ifndef NDEBUG2
  fprintf (stdout, "\n@!!!!hgraph_random returns: ");
  hgraph_fdump (stdout, pr->man, a);
  fprintf (stdout, "\n");
#endif
  return a;
}

inline void
hgraph_free_mem (hgraph_internal_t * pr, hgraph_t * a)
{
  if (a)
    free (a);
}

/* Free memory only if not in the htable */
inline void
hgraph_free_internal (hgraph_internal_t * pr, hgraph_t * a)
{
  hgraph_t *r;
  unsigned keylen;
  if (a)
    {
      /* search in the htable */
      keylen =
              offsetof (hgraph_t,
                        info) + (a->ptrdim + a->size) * sizeof (node_info_t) -
              offsetof (hgraph_t, size);
      r = NULL;
      HASH_FIND (hh, pr->hgraphs, &a->size, keylen, r);
      if (!r || r != a)
        hgraph_free_mem (pr, a);
    }
}

/** Check the integrity of the canonical representation
 *  for heap graphs:
 *  - each node (except null) is labeled by the min dimension
 *  - nodes are ordered wrt dimensions and succ relation
 *  - cut nodes are really cut
 *  - if closed is true, then no more than ms->max_anon nodes
 *    between cut nodes
 *
 * @return:
 *  - '!' if not canonical
 *  - 'c' if closed
 *  - '.' if everything is ok
 */
char
hgraph_check (hgraph_internal_t * pr, hgraph_t * a)
{
  size_t i;
  if (!a)
    {
      ERROR ("graph null!",;
             );
      return '!';
    }

  /*
   * the array a->info is sorted wrt lexicographical order on (nn,v)
   */
  if ((i = hgraph_node_is_sorted (a)) < a->size)
    {
      ERROR ("invalid graph: sorted nodes",;
             );
      fprintf (stderr,
               "node %zu (v=%zu,nn=%zu) > node %zu (v=%zu,nn=%zu)\n",
               i, NODE_VAR (a, i), NODE_VAR_NEXT (a, i),
               i + 1, NODE_VAR (a, i + 1), NODE_VAR_NEXT (a, i + 1));
      hgraph_fdump (stderr, pr->man, a);
      return '!';
    }
  /*
   * each node (except null) is labeled by
   * the min dimension or ptrdim * (nb of next) + min dim
   */
  for (i = 0; i < a->size; i++)
    {
      if (NODE_VAR_NEXT (a, i) == 0 && NODE_VAR (a, i) != NULL_DIM)
        {
          /* labeled cut node, verify that v is the min dim for i */
          size_t j;
          for (j = 0; j < NODE_VAR (a, i); j++)
            if (VAR2NODE (a, j) == i)
              {
                ERROR ("invalid graph: cut",;
                       );
                fprintf (stderr,
                         "cut node %zu, labeled by %zu, also labeled by %zu\n",
                         i, NODE_VAR (a, i), j);
                return '!';
              }
        }
      else if (NODE_VAR_NEXT (a, i) > 0 && NODE_VAR (a, i) == NULL_DIM)
        {
          ERROR ("invalid graph: null",;
                 );
          return '!';
        }
      else if (NODE_VAR (a, i) != NULL_DIM)
        {
          size_t ncut = VAR2NODE (a, NODE_VAR (a, i));
          size_t j;
          /* verify that i is nsucc of ncut */
          for (j = 0; j < NODE_VAR_NEXT (a, i); j++)
            ncut = NODE_NEXT (a, ncut);
          if (ncut != i)
            {
              ERROR ("invalid graph: scut",;
                     );
              return '!';
            }
        }
    }
  /* verify closed flag */
  if (a->closed)
    {
      size_t n;
      for (i = 1, n = 0; i < a->size; i++)
        if (NODE_VAR_NEXT (a, i) > pr->max_anon)
          n++;
      if (n && !(n % pr->segm_anon))
        {
          ERROR ("invalid closure",;
                 );
          return 'c';
        }
    }
  /* regular graph */
  return '.';
}

bool
hgraph_isequal_internal (hgraph_internal_t * pr, hgraph_t * a, hgraph_t * b)
{
  hgraph_t *ra, *rb;
  unsigned keylen_a, keylen_b;
  /* search in the htable */
  keylen_a =
          offsetof (hgraph_t,
                    info) + (a->ptrdim + a->size) * sizeof (node_info_t) -
          offsetof (hgraph_t, size);
  keylen_b =
          offsetof (hgraph_t,
                    info) + (b->ptrdim + b->size) * sizeof (node_info_t) -
          offsetof (hgraph_t, size);
  ra = rb = NULL;
  HASH_FIND (hh, pr->hgraphs, &a->size, keylen_a, ra);
  HASH_FIND (hh, pr->hgraphs, &b->size, keylen_b, rb);
  if (ra == rb)
    return true;
  return false;
}

/* TODO: needed? */

/* Comparison using htable */
int
hgraph_cmp_internal (hgraph_internal_t * pr, hgraph_t * a, hgraph_t * b)
{
  /* if equality in htable return 0 */
  if (hgraph_isequal_internal (pr, a, b) == true)
    return 0;

  return hgraph_cmp (a, b);
}

hgraph_t *
hgraph_copy_mem (hgraph_internal_t * pr, hgraph_t * a)
{
  hgraph_t *r;
  size_t i;

  if (!a)
    return NULL;
  r = hgraph_alloc_internal (pr, a->size, a->datadim, a->ptrdim);
  for (i = 0; i < (a->ptrdim + a->size); i++)
    node_info_copy (&r->info[i], &a->info[i]);
  r->closed = a->closed;
  return r;
}

/* Return a if in htable, or put a in htable and return its representant */
hgraph_t *
hgraph_copy_internal (hgraph_internal_t * pr, hgraph_t * a)
{
  hgraph_t *r;
  unsigned keylen;
  if (!a)
    return NULL;

  /* search in the htable */
  keylen =
          offsetof (hgraph_t,
                    info) + (a->ptrdim + a->size) * sizeof (node_info_t) -
          offsetof (hgraph_t, size);
  r = NULL;
  HASH_FIND (hh, pr->hgraphs, &a->size, keylen, r);
  if (!r)
    {
      /* put a in htable */
      HASH_ADD (hh, pr->hgraphs, size, keylen, a);
      HASH_FIND (hh, pr->hgraphs, &a->size, keylen, r);
    }
  return r;
}

hgraph_array_t *
hgraph_array_make (hgraph_internal_t * pr, size_t size)
{
  hgraph_array_t *r;
  size_t i;
  checked_malloc (r, hgraph_array_t, sizeof (hgraph_array_t), 1, return NULL;
                  );
  r->size = size;
  r->p =
          (size == 0) ? NULL : (hgraph_t **) malloc (size * sizeof (hgraph_t *));
  for (i = 0; i < size; i++)
    r->p[i] = NULL;
  return r;
}

void
hgraph_array_clear (hgraph_internal_t * pr, hgraph_array_t * a)
{
  size_t i;
  if (!a)
    return;
  if (a->p)
    {
      for (i = 0; i < a->size; i++)
        if (a->p[i])
          {
            hgraph_free_internal (pr, a->p[i]);
            a->p[i] = NULL;
          }
      free (a->p);
    }
  free (a);
}

/* Get the same content with new size @code{size}; NULL positions may be added. */
int
hgraph_array_resize (hgraph_internal_t * pr, hgraph_array_t * array,
                     size_t size)
{
  size_t i;
  if (!array)
    return 0;
  for (i = size; i < array->size; i++)
    {
      hgraph_free_internal (pr, array->p[i]);
    }
  array->p = (hgraph_t **) realloc (array->p, size * sizeof (hgraph_t *));
  for (i = array->size; i < size; i++)
    {
      array->p[i] = NULL;
    }
  array->size = size;
  return size;
}

hgraph_array_t *
hgraph_array_copy (hgraph_internal_t * pr, hgraph_array_t * a)
{
  hgraph_array_t *r;
  if (!a)
    r = NULL;
  else
    {
      size_t i;
      r = hgraph_array_make (pr, a->size);
      for (i = 0; i < a->size; i++)
        r->p[i] = hgraph_copy (pr->man, a->p[i]);
    }
  return r;
}

/* Add array @code{a} to @code{dst) by avoiding multiples copies.
 * @code{dst->p} is changed. */
int
hgraph_array_add (hgraph_internal_t * pr, hgraph_array_t * arr, bool docopy,
                  bool destructive, hgraph_t * a)
{
  bool found;
  size_t i, j;
  int r;
  if (!a)
    return 0;
  if (!arr)
    {
      if (destructive)
        hgraph_free_internal (pr, a);
      return 0;
    }
  found = false;
  j = arr->size;
  for (i = 0; i < arr->size && !found && j == arr->size; i++)
    if (arr->p[i] == NULL)
      j = i;
    else if (hgraph_is_equal (arr->p[i], a))
      found = true;
  if (!found)
    {
      if (j == arr->size)
        hgraph_array_resize (pr, arr, arr->size + 4);
      arr->p[j] = (docopy) ? hgraph_copy_internal (pr, a) : a;
      r = 1;
    }
  else
    r = 0;
  if (destructive)
    hgraph_free_internal (pr, a);
  return r;
}

/*
 * Merge both arrays (which may be NULL).
 *
 * @return an array built from both or NULL
 */
hgraph_array_t *
hgraph_array_merge (hgraph_internal_t * pr, bool destructive,
                    hgraph_array_t * a, hgraph_array_t * b)
{
  hgraph_array_t *r;
  if (!a && !b)
    r = NULL;
  else if (!a && b)
    r = (destructive) ? b : hgraph_array_copy (pr, b);
  else if (a && !b)
    r = (destructive) ? a : hgraph_array_copy (pr, a);
  else
    {
      size_t i, j;
      r = hgraph_array_make (pr, a->size + 4);
      j = 0;
      for (i = 0; i < a->size; i++)
        j += hgraph_array_add (pr, r, false, false, a->p[i]);
      for (i = 0; i < b->size; i++)
        j += hgraph_array_add (pr, r, false, false, b->p[i]);
      if (destructive)
        {
          hgraph_array_clear (pr, a);
          hgraph_array_clear (pr, b);
        }
      hgraph_array_resize (pr, r, j);
    }
  return r;
}

/* ============================================================ */
/* Direct Management: Use with care! */
/* ============================================================ */

/* Return the successor in a of node n */
inline node_t
hgraph_node_get_succ (hgraph_t * a, node_t n)
{
  if (n < a->size)
    return NODE_NEXT (a, n);
  return NODE_T_TOP;
}

/*
 * The cut node successor of the node n and returns in anon the number of
 * anonymous nodes reached.
 */
node_t
hgraph_node_get_succ_cut (hgraph_t * a, node_t n, int *anon)
{
  node_t s = NODE_NEXT (a, n);
  *anon = 0;
  while (s != n && hgraph_node_is_cut (a, s) == false)
    {
      s = NODE_NEXT (a, s);
      *anon = *anon + 1;
    }
  return s;
}

/* Return the @code{nsucc} successor on n in a.
 * Cycles are traversed.
 */
size_t
hgraph_node_get_nsucc (hgraph_t * a, node_t n, size_t nsucc)
{
  size_t r = n;
  while (nsucc > 0 && r < a->size)
    {
      r = NODE_NEXT (a, r);
      nsucc--;
    }
  return r;
}

/* Compute the number of predecessors for each node */
size_t *
hgraph_node_get_prev (hgraph_t * a)
{
  size_t *r = NULL;
  size_t i;
  r = (size_t *) malloc (a->size * sizeof (size_t));
  memset (r, 0, a->size * sizeof (size_t));
  for (i = 0; i < a->size; i++)
    if (NODE_NEXT (a, i) < a->size)
      r[NODE_NEXT (a, i)] += 1;
  return r;
}

/* Return the position of a node labeled by (nn,v) in a */
size_t
hgraph_node_get_position (hgraph_t * a, size_t nn, size_t v)
{
  size_t i;
  node_info_t n;
  assert (v < a->ptrdim || v == NULL_DIM);
  n.nn = nn;
  n.v = v;
  i = 0;
  // the array a->info[ptrdim..ptrdim+size] is sorted
  while (i < a->size)
    {
      if (hgraph_node_is_lt (&a->info[a->ptrdim + i], &n))
        i++;
      else
        break;
    }
  return i;
}

/* TODO: priority 3 define correctly cut node */
bool
hgraph_node_is_cut (hgraph_t * a, node_t n)
{
  if (n < a->size
      && (NODE_VAR (a, n) < a->ptrdim || n == 0)
      && NODE_VAR_NEXT (a, n) == 0)
    return true;
  return false;
}

/* TODO: priority 3 define correctly used node */
inline bool
hgraph_node_is_used (hgraph_t * a, node_t n)
{
  return true;
}

bool
hgraph_node_is_reachable (hgraph_t * a, node_t init, node_t target)
{
  size_t i;
  bool *visited;
  if (!a || (a->size <= 1))
    return false;
  visited = (bool *) malloc (a->size * sizeof (bool));
  for (i = 0; i < a->size; i++)
    visited[i] = false;
  visited[init] = true;
  i = hgraph_node_get_succ (a, init);
  while (!visited[i] && i != target)
    {
      visited[i] = true;
      i = hgraph_node_get_succ (a, i);
    }
  free (visited);
  return (i == target);
}

bool
hgraph_node_is_lt (node_info_t * a, node_info_t * b)
{
  // The lexicographic ordering (nn,v) is not good for folding
  // because nodes folded shall be ordered.
  // if (a->nn < b->nn || (a->nn == b->nn && a->v < b->v))
  // Then use lexicographic ordering on (v,nn).
  if ((a->v == NULL_DIM)
      || (a->v < b->v)
      || (a->v == b->v && a->nn < b->nn))
    return true;
  else
    return false;
}

size_t
hgraph_node_is_sorted (hgraph_t * a)
{
  size_t i;
  for (i = 0; i < a->size && (i + 1) < a->size; i++)
    if (hgraph_node_is_lt
        (&a->info[a->ptrdim + i], &a->info[a->ptrdim + i + 1]) == false)
      return i;
  return a->size;
}

bool
hgraph_node_are_iso (hgraph_internal_t * pr, hgraph_t * a, node_t n1, node_t n2)
{
  bool result = true;
  node_t nn1, nn2;
  // arrays of mapping between nodes visited from n1 and n2
  size_t * visited = (size_t*) malloc (a->size * sizeof (size_t));
  memset (visited, 0, a->size * sizeof (size_t));
  nn1 = n1;
  nn2 = n2;
  do
    {
      visited[nn1] = nn2;
      nn1 = NODE_NEXT (a, nn1);
      nn2 = NODE_NEXT (a, nn2);
    }
  while (nn1 != NODE_NULL && nn2 != NODE_NULL
         && nn1 != nn2 /* for graphs disjoints?? */
         && visited[nn1] != nn2);
  if (visited[nn1] != nn2)
    result = false;
  free (visited);
  return result;
}

/** Add node with info (s,v,nn) to a and return its position.
 * The permutation produced has (old) a->size nodes
 */
hgraph_t *
hgraph_node_add (hgraph_internal_t * pr, hgraph_t * a, node_t s, size_t v,
                 size_t nn, node_t * n, ap_dimperm_t * perm)
{
  node_t i, j;
  hgraph_t *r;
  assert (a && perm && a->size == perm->size);
  /* init permutation */
  ap_dimperm_set_id (perm);
  /* search the position to do the insertion */
  i = hgraph_node_get_position (a, nn, v);
  /* test that the new node can be added */
#ifndef NDEBUG1
  fprintf (stdout, "!!!! hgraph_node_add: position to insert i=%zu for (x%zu,nn=%zu)\n",
           i, v, nn);
  fflush (stdout);
#endif
  assert (!
          (i < a->size && NODE_VAR_NEXT (a, i) == nn
           && NODE_VAR (a, i) == v));
  /* resize a->info to a->size+1 */
  r = hgraph_resize (pr, a, a->size + 1);
  /* shift right of one position until i in a->info and perm */
  j = a->size;
  while (i < j)
    {
      /* copy j-1 into j */
      NODE_NEXT (r, j) =
              (NODE_NEXT (a, j - 1) < i) ? NODE_NEXT (a,
                                                      j - 1) : (NODE_NEXT (a,
                                                                           j - 1) +
                                                                1);
      NODE_VAR (r, j) = NODE_VAR (a, j - 1);
      NODE_VAR_NEXT (r, j) = (NODE_VAR (a, j - 1) == v
                              && NODE_VAR_NEXT (a,
                                                j - 1) >=
                              nn) ? (NODE_VAR_NEXT (a,
                                                    j - 1) +
                                     1) : NODE_VAR_NEXT (a, j - 1);
      perm->dim[j - 1] = j;
      j--;
    }
  /* shift nodes >= i in a->v2n */
  for (j = 0; j < r->ptrdim; j++)
    if (VAR2NODE (r, j) < r->size && VAR2NODE (r, j) >= i)
      VAR2NODE (r, j) += 1;
  /* set infos for position i */
  NODE_NEXT (r, i) = s;
  NODE_VAR (r, i) = v;
  NODE_VAR_NEXT (r, i) = nn;
  if (nn == 0)
    VAR2NODE (r, v) = i;
  /* update successors */
  for (j = 0; j <= i; j++)
    if (NODE_NEXT (r, j) >= i)
      NODE_NEXT (r, j) = NODE_NEXT (r, j) + 1;
  *n = i;
  return r;
}

/*
 * Add a node between n and its successor labeled by (v,nn);
 * produce permutation of size a->size+1, the node added is given at a->size.
 */
hgraph_t *
hgraph_node_expand (hgraph_internal_t * pr, hgraph_t * a, node_t n,
                    size_t v, size_t nn, ap_dimperm_t * perm)
{
  node_t i, n_succ;
  hgraph_t *r;
  assert (a && perm && (a->size + 1) == perm->size);
  /* init permutation */
  ap_dimperm_set_id (perm);
  /* resize a->info to a->size+1 */
  r = hgraph_resize (pr, a, a->size + 1);
  /* set info about the new node */
  i = a->size; // node added
  n_succ = NODE_NEXT (a, n); // succ of n in the old graph
  NODE_NEXT (r, n) = i;
  NODE_NEXT (r, i) = n_succ;
  NODE_VAR (r, i) = v;
  NODE_VAR_NEXT (r, i) = nn;
  if (nn == 0) VAR2NODE (r, v) = i;
  if (NODE_VAR_NEXT (r, n_succ) != 0)
    {
      NODE_VAR (r, n_succ) = v;
      NODE_VAR_NEXT (r, n_succ) = nn + 1;
    }
  // sort the graph to obtain the node permutation
  hgraph_node_sort (r, 1, perm);
  return r;
}

/*
 * Start sorting a->info from position i (inclusive) and
 * return permutation of nodes in perm.
 * Sorting of nodes is done wrt order defined in @see{hgraph_node_is_lt}.
 * Requires: nodes are labeled by the min dimension.
 * TODO: deal correctly with garbage nodes!!!
 */
void
hgraph_node_sort (hgraph_t * a, size_t i, ap_dimperm_t * perm)
{
  size_t j;
#ifndef NDEBUG2
  fprintf (stdout, "\n!!!!hgraph_node_sort: on a=(");
  hgraph_fdump (stdout, NULL, a);
  fprintf (stdout, ") from %zu with perm=(", i);
  ap_dimperm_fprint (stdout, perm);
  fprintf (stdout, ")\n");
#endif
  ap_dimperm_t perm1, perm2, permj;
  ap_dimperm_init (&perm1, a->size);
  ap_dimperm_init (&perm2, a->size);
  ap_dimperm_init (&permj, a->size);
  /* init perm to identity */
  ap_dimperm_set_id (&perm1);
  ap_dimperm_set_id (&perm2);
  /* start insertion sort from position i */
  for (j = i + 1; j < a->size; j++)
    // insert only nodes which are not garbage
    if (perm && perm->dim[j] != 0)
      {
        ap_dimperm_set_id (&permj);
        node_info_t nj = a->info[a->ptrdim + j];
        node_t k = j - 1;
        while (k >= i && hgraph_node_is_lt (&nj, &a->info[a->ptrdim + k]))
          {
            a->info[a->ptrdim + k + 1] = a->info[a->ptrdim + k];
            permj.dim[k] = k + 1;
            k = k - 1;
          }
        a->info[a->ptrdim + k + 1] = nj;
        permj.dim[j] = k + 1;
        ap_dimperm_compose (&perm2, &perm1, &permj);
#ifndef NDEBUG2
        fprintf (stdout, "\n!!!!hgraph_node_sort: (perm1, permj, perm2)=(");
        ap_dimperm_fprint (stdout, &perm1);
        ap_dimperm_fprint (stdout, &permj);
        ap_dimperm_fprint (stdout, &perm2);
        fprintf (stdout, ")\n");
#endif
        memcpy (perm1.dim, perm2.dim, perm1.size * sizeof (ap_dim_t));
      }
  ap_dimperm_clear (&perm2);
  ap_dimperm_clear (&permj);
  /* change a->v2n according to perm1 */
  for (j = 0; j < a->ptrdim; j++)
    if (VAR2NODE (a, j) < a->size)
      VAR2NODE (a, j) = perm1.dim[VAR2NODE (a, j)];
  /* change a->info according to perm1 */
  for (j = 0; j < a->size; j++)
    NODE_NEXT (a, j) = perm1.dim[NODE_NEXT (a, j)];
  /* deal with permutations */
  if (perm && perm->size == perm1.size)
    {
      /* compose permutations */
      ap_dimperm_init (&perm2, perm1.size);
      ap_dimperm_compose (&perm2, perm, &perm1);
      memcpy (perm->dim, perm2.dim, perm->size * sizeof (ap_dim_t));
      ap_dimperm_clear (&perm2);
    }
  else if (perm)
    {
      /* copy perm1 in perm */
      ap_dimperm_clear (perm);
      ap_dimperm_init (perm, perm1.size);
      memcpy (perm->dim, perm1.dim, perm->size * sizeof (ap_dim_t));
    }
#ifndef NDEBUG2
  fprintf (stdout, "\n!!!!hgraph_node_sort: returns a=(");
  hgraph_fdump (stdout, NULL, a);
  fprintf (stdout, ") and perm=(");
  ap_dimperm_fprint (stdout, perm);
  fprintf (stdout, ")\n");
#endif
  ap_dimperm_clear (&perm1);
}

/*
 * Set n to be labeled by v. Permutation perm is composed with the produced
 * permutation of nodes to obtain the result.
 */
void
hgraph_node_set_var (hgraph_t * a, node_t n, size_t v, ap_dimperm_t * perm)
{
  size_t i, nn;
  node_t s;
  bool sorted = false;

  assert (a != NULL && perm != NULL && n < a->size && v < a->ptrdim);
  /* case (1): nothing is changed for nodes */
  if (n == 0 || NODE_VAR (a, n) <= v)
    {
      VAR2NODE (a, v) = n;
      return;
    }
  /* compute in i the position of the node labeled by (0,v) */
  i = hgraph_node_get_position (a, 0, v);
  if (n == i)
    sorted = true;
  /* change the labeling of v */
  VAR2NODE (a, v) = n;
  /* change the labeling of n */
  NODE_VAR (a, n) = v;
  NODE_VAR_NEXT (a, n) = 0;
  /* change the labeling of anonymous successors of n if needed */
  nn = 1;
  s = NODE_NEXT (a, n);
  while (s < a->size &&
         (NODE_VAR_NEXT (a, s) >= nn
          || (NODE_VAR_NEXT (a, s) == nn && NODE_VAR (a, s) > v)))
    {
      NODE_VAR (a, s) = v;
      NODE_VAR_NEXT (a, s) = nn;
      s = NODE_NEXT (a, s);
      nn++;
    }
  /* case (2): info changes but not need to move nodes */
  if (!sorted)
    {
      /* case (3): a->info shall be re-sorted */
      hgraph_node_sort (a, i, perm);
    }
}

/*
 * Set the successor of n to be m.
 * Requires: the current successor of n is not defined
 * Ensures: If the canonical form is broken, re-sort nodes
 * and produce a permutation in perm. If perm is not initially empty, the
 * permutation produced is composed with perm.
 * TODO: since n cannot be old node, this is not detecting garbage or anonymous paths!
 */
void
hgraph_node_set_succ (hgraph_t * a, node_t n, node_t m, ap_dimperm_t * perm)
{
  size_t nn = 0;
  node_t s;
  ap_dimperm_t *perm2;
  // assert (a != NULL && perm != NULL && n < a->size && a->info[n].s == NODE_T_TOP && m < a->size);
  assert (a != NULL && perm != NULL && n < a->size && m < a->size);
  s = NODE_NEXT (a, n) = m;
  /* change the labeling of anonymous successors of n if needed */
  nn = NODE_VAR_NEXT (a, n);
  while (s < a->size && NODE_VAR_NEXT (a, s) > 0
         && NODE_VAR (a, s) > NODE_VAR (a, n))
    {
      NODE_VAR (a, s) = NODE_VAR (a, n);
      NODE_VAR_NEXT (a, s) = ++nn;
      s = NODE_NEXT (a, s);
    }
  /* check sorted */
  s = hgraph_node_is_sorted (a);
  if (s >= a->size)
    return;
  else
    /* do sorting procedure starting with s */
    hgraph_node_sort (a, s, perm);
  return;
}

/* Return the node labeled by v */
node_t
hgraph_get_var_node (hgraph_t * a, size_t v)
{
  if (v == NULL_DIM)
    return 0;
  if (v < a->ptrdim)
    {
      node_t n = VAR2NODE (a, v);
      if (n < a->size)
        return n;
    }
  return hgraph_node_get_null (a);
}

/* Getter for the mapping (allocated) var 2 nodes */
size_t *
hgraph_get_var2node (hgraph_t * a)
{
  size_t *r = (size_t *) malloc (a->ptrdim * sizeof (size_t));
  size_t i;
  for (i = 0; i < a->ptrdim; i++)
    r[i] = VAR2NODE (a, i);
  return r;
}

/* Copy a with size nodes; nodes added are not initialized. */
hgraph_t *
hgraph_resize (hgraph_internal_t * pr, hgraph_t * a, size_t size)
{
  hgraph_t *r;
  size_t i;
  if (a == NULL)
    return NULL;
  r = hgraph_alloc_internal (pr, size, a->datadim, a->ptrdim);
  r->closed = false;
  r->size = size;
  for (i = 0; i < (a->ptrdim + size) && i < (a->ptrdim + a->size); i++)
    node_info_copy (&r->info[i], &a->info[i]);
  return r;
}

/* ============================================================ */
/* Memory */

/* ============================================================ */

hgraph_t *
hgraph_copy (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_COPY, 0);
  return hgraph_copy_internal (pr, a);
}

void
hgraph_free (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_FREE, 0);
  hgraph_free_internal (pr, a);
}

size_t
hgraph_size (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_ASIZE, 0);
  return sizeof (hgraph_t) + (a->ptrdim + a->size) * sizeof (node_info_t);
}


/* ============================================================ */
/* Control of internal representation */
/* ============================================================ */

/* TODO: priority 3 */

/* Eliminate anonymous nodes */
void
hgraph_minimize (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_MINIMIZE, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
}

/* Put to NULL all ptr vars on labeling TOP */
void
hgraph_canonicalize (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);
  size_t v;
  for (v = 0; v < a->ptrdim; v++)
    if (VAR2NODE (a, v) == NODE_T_TOP)
      VAR2NODE (a, v) = NODE_NULL;
  return;
}

/* TODO: priority 0 */
int
hgraph_hash (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_HASH, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return 0;
}

/* NOT IMPLEMENTED: do nothing */
void
hgraph_approximate (ap_manager_t * man, hgraph_t * a, int algorithm)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_APPROXIMATE, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
}

/* TODO: priority 3 */
bool
hgraph_is_minimal (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return true;
}

/* TODO: priority 3 */
bool
hgraph_is_canonical (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_CANONICALIZE, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return true;
}


/* ============================================================ */
/* Basic Constructors */

/* ============================================================ */

hgraph_t *
hgraph_bottom (ap_manager_t * man, size_t intdim, size_t realdim)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_BOTTOM, 0);
  hgraph_t *r =
          hgraph_alloc_internal (pr, 0, intdim, REAL2PTR_DIM (pr, realdim));
  return r;
}

/*
 * Top value is a hgraph with one node (null) and prg ptr vars pointing
 * outside this node.
 */
hgraph_t *
hgraph_top (ap_manager_t * man, size_t intdim, size_t realdim)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_TOP, 0);
  hgraph_t *r = hgraph_copy_internal (pr,
                                      hgraph_alloc_internal (pr, 1, intdim,
                                                             REAL2PTR_DIM (pr,
                                                                           realdim)));
  return r;
}

/* put constraints on data variables */
hgraph_t *
hgraph_of_box (ap_manager_t * man, size_t intdim, size_t realdim,
               ap_interval_t ** t)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_OF_BOX, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return NULL;
}


/* ============================================================ */
/* Getters */

/* ============================================================ */

ap_dimension_t
hgraph_dimension (ap_manager_t * man, hgraph_t * a)
{
  hgraph_internal_t *pr =
          hgraph_init_from_manager (man, AP_FUNID_DIMENSION, 0);
  ap_dimension_t r;
  r.intdim = a->datadim;
  r.realdim = PTR2REAL_DIM (pr, a->ptrdim);
  return r;
}


/* ============================================================ */
/* Managers */

/* ============================================================ */

void
hgraph_internal_free (hgraph_internal_t * pr)
{
  /* TODO: free htable */
  pr->hgraphs = NULL;
  pr->pcons = NULL;
  pr->passigns = NULL;
  free (pr);
}

ap_manager_t *
hgraph_manager_alloc (void)
{
  size_t i;
  ap_manager_t *man;
  hgraph_internal_t *pr;

  pr = (hgraph_internal_t *) malloc (sizeof (hgraph_internal_t));
  assert (pr);
  pr->hgraphs = NULL;
  pr->pcons = NULL;
  pr->passigns = NULL;

  pr->size_scons = 0;
  pr->man_scons = NULL;

  pr->meet_algo = 0;
  pr->max_anon = 4;
  pr->segm_anon = 1;
  pr->error_ = 0;

  man = ap_manager_alloc ("shape", "0.1 with hgraphs", pr,
                          (void (*)(void *)) hgraph_internal_free);
  pr->man = man;

  man->funptr[AP_FUNID_COPY] = &hgraph_copy;
  man->funptr[AP_FUNID_FREE] = &hgraph_free;
  man->funptr[AP_FUNID_ASIZE] = &hgraph_size;
  man->funptr[AP_FUNID_MINIMIZE] = &hgraph_minimize;
  man->funptr[AP_FUNID_CANONICALIZE] = &hgraph_canonicalize;
  man->funptr[AP_FUNID_HASH] = &hgraph_hash;
  man->funptr[AP_FUNID_APPROXIMATE] = &hgraph_approximate;
  man->funptr[AP_FUNID_FPRINT] = &hgraph_fprint;
  man->funptr[AP_FUNID_FPRINTDIFF] = &hgraph_fprintdiff;
  man->funptr[AP_FUNID_FDUMP] = &hgraph_fdump;
  man->funptr[AP_FUNID_SERIALIZE_RAW] = &hgraph_serialize_raw;
  man->funptr[AP_FUNID_DESERIALIZE_RAW] = &hgraph_deserialize_raw;
  man->funptr[AP_FUNID_BOTTOM] = &hgraph_bottom;
  man->funptr[AP_FUNID_TOP] = &hgraph_top;
  man->funptr[AP_FUNID_OF_BOX] = &hgraph_of_box;
  man->funptr[AP_FUNID_DIMENSION] = &hgraph_dimension;
  man->funptr[AP_FUNID_IS_BOTTOM] = &hgraph_is_bottom;
  man->funptr[AP_FUNID_IS_TOP] = &hgraph_is_top;
  man->funptr[AP_FUNID_IS_LEQ] = &hgraph_is_leq;
  man->funptr[AP_FUNID_IS_EQ] = &hgraph_is_eq;
  man->funptr[AP_FUNID_IS_DIMENSION_UNCONSTRAINED] =
          &hgraph_is_dimension_unconstrained;
  man->funptr[AP_FUNID_SAT_INTERVAL] = &hgraph_sat_interval;
  man->funptr[AP_FUNID_SAT_LINCONS] = &hgraph_sat_lincons;
  man->funptr[AP_FUNID_SAT_TCONS] = &hgraph_sat_tcons;
  man->funptr[AP_FUNID_BOUND_DIMENSION] = &hgraph_bound_dimension;
  man->funptr[AP_FUNID_BOUND_LINEXPR] = &hgraph_bound_linexpr;
  man->funptr[AP_FUNID_BOUND_TEXPR] = &hgraph_bound_texpr;
  man->funptr[AP_FUNID_TO_BOX] = &hgraph_to_box;
  man->funptr[AP_FUNID_TO_LINCONS_ARRAY] = &hgraph_to_lincons_array;
  man->funptr[AP_FUNID_TO_TCONS_ARRAY] = &hgraph_to_tcons_array;
  man->funptr[AP_FUNID_TO_GENERATOR_ARRAY] = &hgraph_to_generator_array;
  man->funptr[AP_FUNID_MEET] = &hgraph_meet;
  man->funptr[AP_FUNID_MEET_ARRAY] = &hgraph_meet_array;
  man->funptr[AP_FUNID_MEET_LINCONS_ARRAY] = &hgraph_meet_lincons_array;
  man->funptr[AP_FUNID_MEET_TCONS_ARRAY] = &hgraph_meet_tcons_array;
  man->funptr[AP_FUNID_JOIN] = &hgraph_join;
  man->funptr[AP_FUNID_JOIN_ARRAY] = &hgraph_join_array;
  man->funptr[AP_FUNID_ADD_RAY_ARRAY] = &hgraph_add_ray_array;
  man->funptr[AP_FUNID_ASSIGN_LINEXPR_ARRAY] = &hgraph_assign_linexpr_array;
  man->funptr[AP_FUNID_SUBSTITUTE_LINEXPR_ARRAY] =
          &hgraph_substitute_linexpr_array;
  man->funptr[AP_FUNID_ASSIGN_TEXPR_ARRAY] = &hgraph_assign_texpr_array;
  man->funptr[AP_FUNID_SUBSTITUTE_TEXPR_ARRAY] =
          &hgraph_substitute_texpr_array;
  man->funptr[AP_FUNID_ADD_DIMENSIONS] = &hgraph_add_dimensions;
  man->funptr[AP_FUNID_REMOVE_DIMENSIONS] = &hgraph_remove_dimensions;
  man->funptr[AP_FUNID_PERMUTE_DIMENSIONS] = &hgraph_permute_dimensions;
  man->funptr[AP_FUNID_FORGET_ARRAY] = &hgraph_forget_array;
  man->funptr[AP_FUNID_EXPAND] = &hgraph_expand;
  man->funptr[AP_FUNID_FOLD] = &hgraph_fold;
  man->funptr[AP_FUNID_WIDENING] = &hgraph_widening;
  man->funptr[AP_FUNID_CLOSURE] = &hgraph_closure;

  for (i = 0; i < AP_EXC_SIZE; i++)
    ap_manager_set_abort_if_exception (man, i, false);

  return man;
}

hgraph_t *
hgraph_of_abstract0 (ap_abstract0_t * a)
{
  return (hgraph_t *) a->value;
}

ap_abstract0_t *
abstract0_of_hgraph (ap_manager_t * man, hgraph_t * a)
{
  ap_abstract0_t *r = malloc (sizeof (ap_abstract0_t));
  assert (r);
  r->value = a;
  r->man = ap_manager_copy (man);
  return r;
}
