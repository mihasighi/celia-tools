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
#include "shape_macros.h"


/* All closures are in-place. */


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* Full Closure */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/* Mark to remove garbage nodes from a, compact the representation of a,
 * and return the number of garbage nodes.
 *
 * Requires: perm has size a->size
 *
 * Ensures: a is a new hgraph with garbage nodes at the end of the a->info
 *          perm is mapping old nodes to new nodes or 0 if these nodes are garbage
 *          @result is the number of garbage nodes
 */
size_t
hgraph_close_garbage (hgraph_internal_t * pr, hgraph_t * a,
                      ap_dimperm_t * perm, bool cutpoint) /* TODO */
{
  node_info_t *newinfo;
  size_t i, v, ngarbage, newsize;
  arg_assert (perm && perm->size == a->size, return 0;);

  // compute the newinfo by a hgraph traversal from the cut nodes labeled by variables.
  // init of newinfo
  checked_malloc (newinfo, node_info_t, sizeof (node_info_t), a->size,
                  return 0;);
  // node 0 is null
  newinfo[0].v = NULL_DIM;
  newinfo[0].nn = 0;
  newinfo[0].s = 0;
  for (i = 1; i < a->size; i++)
    {
      newinfo[i].v = NODE_T_TOP;
      newinfo[i].nn = NODE_T_TOP;
      newinfo[i].s = 0;
    }
  // start traversal and mark newinfo for labeled nodes
  for (v = 0; v < a->ptrdim; v++)
    {
      i = VAR2NODE (a, v);
      if (i != NODE_T_TOP &&
          i != NODE_NULL &&
          newinfo[i].v == NODE_T_TOP)
        {
          newinfo[i].v = v;
          newinfo[i].nn = 0;
          newinfo[i].s = NODE_NEXT (a, i);
        }
    }
#ifndef NDEBUG2
  fprintf (stdout, "\n!!!!hgraph_close_garbage: newinfo=(");
  // print new info
  for (i = 0; i < a->size; i++)
    fprintf (stdout, "\nnewinfo[%zu]=(v=%zu,nn=%zu,s=%zu)", i, newinfo[i].v,
             newinfo[i].nn, newinfo[i].s);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif
  // mark now anonymous nodes reachable from labeled nodes
  for (v = 0; v < a->ptrdim; v++)
    {
      size_t nn = 1;
      i = DIM2NODE (a, v);
      i = (i > 0 && i < a->size) ? newinfo[i].s : i;
      while (i > 0 && i < a->size && newinfo[i].v == NODE_T_TOP)
        {
          newinfo[i].v = v;
          newinfo[i].nn = nn;
          newinfo[i].s = NODE_NEXT (a, i);
          i = NODE_NEXT (a, i);
          nn++;
        }
    }
#ifndef NDEBUG2
  fprintf (stdout, "\n!!!!hgraph_close_garbage: after marking newinfo=(");
  // print new info
  for (i = 0; i < a->size; i++)
    fprintf (stdout, "\nnewinfo[%zu]=(v=%zu,nn=%zu,s=%zu)", i, newinfo[i].v,
             newinfo[i].nn, newinfo[i].s);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif
  // copy newinfo in a->info and compute garbage
  ngarbage = 0;
  // null is never garbage
  for (i = NODE_NULL + 1; i < a->size; i++)
    {
      a->info[a->ptrdim + i] = newinfo[i];
      if (newinfo[i].v == NODE_T_TOP)
        {
          ngarbage++;
          perm->dim[i] = 0;
          if (cutpoint
              && newinfo[NODE_NEXT (a, i)].v < a->ptrdim // node accessible from a labeled node
              && newinfo[NODE_NEXT (a, i)].nn > 0) // in at least one edge
            {
              ERROR ("Bad change of environment: not cut-point free graph!",;
                     );
            }
        }
    }
  free (newinfo);
#ifndef NDEBUG2
  fprintf (stdout, "\n!!!!hgraph_close_garbage: before sort=(\n");
  hgraph_fdump (stdout, pr->man, a);
  fprintf (stdout, ") with ngarbage=%zu and perm=(", ngarbage);
  ap_dimperm_fprint (stdout, perm);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif
  // sort the info, garbage shall be last ngarbage nodes, any record may change!
  hgraph_node_sort (a, 1, perm);
  return ngarbage;
}

/*
 * Compute anonymous paths in anon and remove garbage.
 * Version: works only for anonymous paths starting from cut nodes
 *          labeled by ptrvars
 */
hgraph_t *
hgraph_close_anonymous (hgraph_internal_t * pr, hgraph_t * a,
                        ap_dimperm_t * perm, ap_dim_array2_t * anon,
                        size_t ngarbage)
{
  size_t i, j, v, nv, ns, newsize, maxanon, nanon;
  arg_assert (perm && perm->size == a->size && anon, return NULL;);

  nanon = 0;

  if (hgraph_is_closed (pr, a))
    {
      if (ngarbage == 0)
        return hgraph_copy_internal (pr, a);
      else
        goto hgraph_close_remove_garbage;
    }

  // revert perm to put original nodes in anon s
  ap_dimperm_t *nperm = ap_dimperm_alloc (a->size);
  ap_dimperm_invert (nperm, perm);

  // compute cut information
  size_t *prevn = hgraph_node_get_prev (a);

  // allocate a number of anonymous segments equal to number of ptrvars
  anon->size = a->ptrdim;
  anon->p = (ap_dim_array_t *) malloc (anon->size * sizeof (ap_dim_array_t));
  memset (anon->p, 0, anon->size * sizeof (ap_dim_array_t));

  // maximal number of anonymous nodes on each path
  maxanon = (pr->max_anon + 1) * pr->segm_anon;
  // consider all cut nodes labeled by ptrvars
  for (v = 0, i = 0; v < a->ptrdim;
          v++)
    {
      // consider that the i^th anonymous path starts from v
      // init anon for this path
      nanon = maxanon + 1;
      anon->p[i].p = (ap_dim_t *) malloc (nanon * sizeof (node_t));
      anon->p[i].size = nanon;
      memset (anon->p[i].p, 0, nanon * sizeof (node_t));
      nv = DIM2NODE (a, v);
      anon->p[i].p[0] = nv; /* changed below to recall old nodes (ap_dim_t) nperm->dim[nv]; */
      j = 1;
      ns = NODE_NEXT (a, nv);
      while (NODE_VAR (a, ns) == v && prevn[ns] == 1 && nanon >= 1)
        {
          anon->p[i].p[j] = ns; /* changed below to recall old nodes */
          ns = NODE_NEXT (a, ns);
          j++;
          nanon--;
        }
      if (j >= (pr->max_anon + 2))
        { // at least pr->max_anon anonymous on this path!
          // realloc the array at the good size and set the size (used below)
          anon->p[i].p = realloc (anon->p[i].p, j * sizeof (ap_dim_array_t));
          anon->p[i].size = j;
          i++;
        }
      else
        {
          free (anon->p[i].p);
          anon->p[i].p = NULL;
          anon->p[i].size = 0;
        }
    } // end for each ptr var
  // TODO: add management of anonymous paths starting from anonymous cut nodes
  // i == anon_segm gives the real number of anonymous paths, resize anon->p
  // permutation used to to sorting of nodes after anonymous elimination
  ap_dimperm_t *perm1 = ap_dimperm_alloc (a->size);
  ap_dimperm_set_id (perm1);
  nanon = 0;
  if (i == 0)
    {
      free (anon->p);
      anon->p = NULL;
      anon->size = 0;
    }
  else
    {
      anon->p = realloc (anon->p, i * sizeof (ap_dim_array_t));
      anon->size = i;

#ifndef NDEBUG2
      fprintf (stdout, "\n!!!!hgraph_close_anonymous: infos about anons=(");
      for (i = 0; i < anon->size; i++)
        {
          fprintf (stdout, "[%zu]: %d-->", i, anon->p[i].p[0]);
          for (j = 1; j < (pr->max_anon + 2); j++)
            fprintf (stdout, "%d,", anon->p[i].p[j]);
        }
      fprintf (stdout, ")\n");
#endif

      // set a->info to eliminate anonymous nodes of paths found
      for (i = 0; i < anon->size; i++)
        {
          // TODO: this has to be changed in order to deal with anonymous cut nodes
          // the successor of node starting the path (nv) is
          // the successor (ns) of the last anonymous in the path, which may change of label
          nv = anon->p[i].p[0];
          ns = NODE_NEXT (a, anon->p[i].p[anon->p[i].size - 1]);
          NODE_NEXT (a, nv) = ns;
          if (NODE_VAR (a, ns) == NODE_VAR (a, nv))
            NODE_VAR_NEXT (a, ns) = 1;
          anon->p[i].p[0] = nperm->dim[nv]; // changed to speak about old nodes
          // the anonymous nodes of the path are labeled by max
          for (j = 1; j < anon->p[i].size; j++)
            {
              node_t nj = anon->p[i].p[j];
              NODE_VAR (a, nj) = NODE_T_TOP;
              NODE_VAR_NEXT (a, nj) = NODE_T_TOP;
              NODE_NEXT (a, nj) = 0;
              anon->p[i].p[j] = nperm->dim[nj];
              // problem with codee below because sorting does not compose well with perm
              // perm->dim[nperm->dim[nj]] = 0; // perm is updated to eliminate this node
              perm1->dim[nj] = 0; // eliminate this node in sort
            }
          nanon += (anon->p[i].size - 1);
        }
#ifndef NDEBUG2
      fprintf (stdout, "\n!!!!hgraph_close_anonymous: new hgraph (unsorted)=(");
      hgraph_fdump (stdout, pr->man, a);
      fprintf (stdout, ") with nanon=%zu and perm1=(", nanon);
      ap_dimperm_fprint (stdout, perm1);
      fprintf (stdout, ")\n");
#endif
    }
  // sort a, then all anonymous nodes are put at the end with the garbage
  hgraph_node_sort (a, 1, perm1); // not perm because perm stores also the initial permutation
  ap_dimperm_compose (perm, perm, perm1);
  ap_dimperm_clear (perm1);

hgraph_close_remove_garbage:
  // resize a
  newsize = a->size - ngarbage - nanon;
  hgraph_t *newa = (hgraph_t *) malloc ((sizeof (hgraph_t) +
                                         (a->ptrdim +
                                          newsize) * sizeof (node_info_t)));
  memcpy (newa, a,
          (sizeof (hgraph_t) + (a->ptrdim + newsize) * sizeof (node_info_t)));
  newa->closed = true;
  newa->size = newsize;
#ifndef NDEBUG2
  fprintf (stdout, "\n!!!!hgraph_close_anonymous: new hgraph (sorted)=(");
  hgraph_fdump (stdout, pr->man, newa);
  fprintf (stdout, ") and new perm (sorted)=(");
  ap_dimperm_fprint (stdout, perm);
  fprintf (stdout, ")\n");
#endif
  hgraph_t *r = hgraph_copy_internal (pr, newa);
  hgraph_free_internal (pr, newa);
  return r;
}

/*
 * Compute anonymous paths in anon.
 * Version: works as well for anonymous paths starting from cut nodes
 *          labeled by ptrvars than for unlabeled cut nodes.
 * TODO: have bugs.
 */
void
hgraph_close_anonymous1 (hgraph_internal_t * pr, hgraph_t * a,
                         ap_dimperm_t * perm, ap_dim_array2_t * anon,
                         size_t ngarbage)
{
  node_info_t *newinfo;
  size_t i, j, v, nv, newsize, anon_segm;
  arg_assert (perm && perm->size == a->size, return;);

  // revert perm to put original nodes in anon s
  ap_dimperm_t *nperm = ap_dimperm_alloc (a->size);
  ap_dimperm_invert (nperm, perm);

  // compute cut information
  size_t *prevn = hgraph_node_get_prev (a);
  newsize = (a->size - ngarbage);
  anon_segm = 0; // number of anonymous segments

  // permutation storing the anonymous information:
  // an anonymous node in an anonymous path is mapped into the head of the path
  ap_dimperm_t *anons_info = ap_dimperm_alloc (a->size);
  ap_dimperm_set_id (anons_info);
  // anon_orig[i] == true iff i is the beginning of a path of anonymous
  bool *anon_orig = (bool *) malloc (newsize * sizeof (bool));
  // allocate newinfo for non-garbage nodes
  checked_malloc (newinfo, node_info_t, sizeof (node_info_t), newsize,
                  return;);
  for (i = 0; i < newsize; i++)
    {
      anon_orig[i] = false;
      newinfo[i] = a->info[a->ptrdim + i];
      if (NODE_VAR_NEXT (a, i) > pr->max_anon && anons_info->dim[i] == i) // not already set in an anonymous path
        { // i is in a path with more than pr->max_anon nodes, check if all nodes are anonymous
          size_t nanon = 0; // number of anonymous from the last cut node
          v = NODE_VAR (a, i);
          nv = VAR2NODE (a, v); // last cut node
          j = NODE_NEXT (a, nv);
          while (j != i && j < a->size)
            {
              if (prevn[j] == 1 && anons_info->dim[j] == j)
                nanon++;
              else if (prevn[j] > 1)
                {
                  nv = j;
                  nanon = 0;
                }
              else // anons_info->dim[j] != j
                {
                  nv = anons_info->dim[j]; // TODO: here is not really true, find the first after this anonymous path
                  nanon = 0;
                }
              j = NODE_NEXT (a, j);
            }
          // add i to anonymous nodes
          nanon++;
          assert (j == i);
          if (nanon > pr->max_anon)
            { // the path is built only from anonymous nodes
              // set anons_info to map nodes on the anon path to nv
              anon_orig[nv] = true;
              anon_segm++;
              newinfo[nv].s = NODE_NEXT (a, i);
              j = NODE_NEXT (a, nv);
              while (j != i)
                { // eliminate the anonymous nodes of the path from the graph
                  newinfo[j].v = NODE_T_TOP; /* max label */
                  newinfo[j].nn = NODE_T_TOP;
                  newinfo[j].s = 0;
                  // set info about its anon path
                  anons_info->dim[j] = nv;
                  perm->dim[nperm->dim[j]] = 0;
                  j = NODE_NEXT (a, j);
                }
              newinfo[i].v = NODE_T_TOP; /* max label */
              newinfo[i].nn = NODE_T_TOP;
              newinfo[i].s = 0;
              anons_info->dim[i] = nv;
              perm->dim[nperm->dim[i]] = 0; // set to garbage this node in perm!
            }
        }
    }
#ifndef NDEBUG2
  fprintf (stdout, "\nInfos about anons:");
  fprintf (stdout, "\n\t newinfo: ");
  for (i = 0; i < newsize; i++)
    fprintf (stdout, "%zu-->(v=%zu,nn=%zu,s=%zu), ", i, newinfo[i].v,
             newinfo[i].nn, newinfo[i].s);
  fprintf (stdout, "\n\t anon_orig: ");
  for (i = 0; i < newsize; i++)
    fprintf (stdout, "%zu-->%d, ", i, anon_orig[i]);
  fprintf (stdout, "\n\t anon_info: ");
  ap_dimperm_fprint (stdout, anons_info);
  fprintf (stdout, "\n");
#endif
  // we have here in anon_segm the number of anonymous segments, in anons_info the mapping
  // fill now the anon parameter
  anon->size = anon_segm;
  if (anon_segm > 0)
    {
      anon->p =
              (ap_dim_array_t *) malloc (anon_segm * sizeof (ap_dim_array_t));
      nv = 0;
      for (j = 0; j < newsize && nv < anon_segm; j++)
        if (anon_orig[j])
          {
            size_t k;
            anon->p[nv].p =
                    (ap_dim_t *) malloc ((pr->max_anon + 2) * sizeof (node_t));
            anon->p[nv].p[0] = (ap_dim_t) nperm->dim[j];
            // put nodes in order!
            j = NODE_NEXT (a, nv);
            for (k = 0; k <= pr->max_anon; k++)
              {
                anon->p[nv].p[k + 1] = (ap_dim_t) nperm->dim[j];
                j = NODE_NEXT (a, j);
              }
            nv++;
          }

      // remove anonymous nodes and set perm!
      // copy newinfo in a->info
      for (i = 0; i < newsize; i++)
        a->info[a->ptrdim + i] = newinfo[i];
      free (newinfo);
      // sort the info, garbage shall be the last ngarbage+anon_segm*pr->max_anon nodes
      hgraph_node_sort (a, 1, perm);

      // build a new a which eliminate garbage and anonymous
      newsize -= pr->max_anon * anon_segm;
      // since size is decreasing, reallocation shall work!
      hgraph_t *newa = (hgraph_t *) realloc (a,
                                             (sizeof (hgraph_t) +
                                              (a->ptrdim +
                                               newsize) *
                                              sizeof (node_info_t)));
      assert (newa == a);
      a->closed = true;
      // a->v2n is not changed
      a->size = newsize;
      // a->datadim is not changed
      // a->ptrdim is not changed
      // a->info is the first newsize entries!
#ifndef NDEBUG2
      fprintf (stdout, "\nInfos about anons:");
      fprintf (stdout, "\n\t anon: ");
      for (nv = 0; nv < anon_segm; nv++)
        {
          fprintf (stdout, "\n\t\t %zu: %d, ", nv, anon->p[nv].p[0]);
          for (i = 0; i < pr->max_anon; i++)
            fprintf (stdout, "%d, ", anon->p[nv].p[i + 1]);
        }
#endif
    }
}

/** Remove from a the garbage nodes mapped by perm into 0.
 * Then, find paths containing more than pr->max_anon anonymous nodes:
 * these paths are put in anon and then eliminated from the graph.
 *
 * Requires: perm has size a->size and maps to 0 all garbage nodes
 *
 * Ensures: a is a new hgraph without garbage and paths of anonymous nodes
 *          perm is mapping old nodes to new nodes or 0 if these nodes are no more in a
 *          anon maps each OLD node of a to an ordered list of pr->max_anon nodes
 */
hgraph_t *
hgraph_close (hgraph_internal_t * pr, hgraph_t * a, ap_dimperm_t * perm,
              ap_dim_array2_t * anon, bool cutpoint)
{
  arg_assert (perm && perm->size == a->size, return NULL;);

  size_t ngarbage = hgraph_close_garbage (pr, a, perm, cutpoint);

  hgraph_t *r = hgraph_close_anonymous (pr, a, perm, anon, ngarbage);
  return r;
}


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* Incremental Closure */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/*
 * TODO: priority 1 Semantics: h must equal to a closed graph, except for
 * constraints involving variable v
 *
 * ?? time. ?? space.
 */

bool
hgraph_close_incremental (hgraph_t * h, size_t dim, size_t v)
{
  return false;
}



/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* Sanity Check */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/* TODO: priority 1 */
bool
hgraph_check_closed (hgraph_t * h, size_t dim)
{
  return false;
}


/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* Topological closure operation */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/* TODO: priority 3 */

/* Eliminate anonimous nodes */
hgraph_t *
hgraph_closure (ap_manager_t * man, bool destructive, hgraph_t * a)
{
  hgraph_internal_t *pr = hgraph_init_from_manager (man, AP_FUNID_CLOSURE, 0);
  if (destructive)
    return a;
  else
    return hgraph_copy_internal (pr, a);
}
