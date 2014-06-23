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


/*
 * Projection, changes of dimension, variable permutation.
 */

#include "uthash.h"
#include "ucons.h"
#include "ucons_fun.h"
#include "ucons_internal.h"
#include "shape_macros.h"
#include "apron2shape.h"

//#if defined (UCONS_DCONS_OCT_P21) || defined (UCONS_DCONS_POLY_P21)
ap_dim_t * fold_dim = NULL;
size_t size_fold_dim = 0;
size_t first = 1;
//#endif





/* ============================================================ */
//* Projections */
/* ============================================================ */

/* TODO: priority 0 */

/*
 * not used because we suppose that all ptr variables are declared from the
 * beginning
 */
ucons_t *
ucons_forget_array (ap_manager_t * man,
                    bool destructive, ucons_t * a,
                    ap_dim_t * tdim, size_t size, bool project)
{
  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_FORGET_ARRAY, 0);



  arg_assert (a && tdim && 0 < size && size < a->segmentdim, return NULL;);


#ifndef NDEBUG1
  size_t temp;
  fprintf (stdout, "  \n\t FORGET_ARRAY dimensions:\n [ ");
  for (size_t temp = 0; temp < size; temp++)
    {
      fprintf (stdout, "%zu ", tdim[temp]);
    }
  fprintf (stdout, "\n]");
  ucons_fprint (stdout, pr->man, a, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  // the dimension is kept
  ucons_t *r = ucons_copy_internal (pr, a);

  ap_dim_t *ntdim;
  size_t i;
  // GO for data and length dimensions

  ntdim = (ap_dim_t *) malloc (2 * size * sizeof (ap_dim_t));
  for (i = 0; i < size; i++)
    {
      ntdim[i] = tdim[i];
      ntdim[size + i] = r->segmentdim + tdim[i];
    }
  size = 2 * size;
  r->econs =
          ap_abstract0_forget_array (pr->man_dcons, true, r->econs, ntdim,
                                     size, project);
  pattern_t *s = r->udcons;
  while (s != NULL)
    {
      s->dcons =
              ap_abstract0_forget_array (pr->man_dcons, true, s->dcons, ntdim,
                                         size, project);
    }
  if (destructive)
    ucons_free_internal (pr, a);
  return r;

}


/* ============================================================ */
/* Change and permutation of dimensions */
/* ============================================================ */

/* T#ifndef NDEBUG1
        size_t temp;
        fprintf(stdout,"  \n\t add dimensions: \n");
        for(size_t temp = 0; temp < dimchange->intdim; temp++){
                fprintf(stdout,"%zu",dimchange->dim[temp]);
        }
        fprintf(stdout,"\n");
        //ucons_fprint(stdout,pr->man,a,NULL);
        //fprintf(stdout,"\n");
        fflush(stdout);
#endifODO: priority 0 */
ucons_t *
ucons_add_dimensions (ap_manager_t * man,
                      bool destructive, ucons_t * a,
                      ap_dimchange_t * dimchange, bool project)
{
  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_ADD_DIMENSIONS, 0);


  arg_assert (a && dimchange && dimchange->dim != NULL, return NULL;);

#ifndef NDEBUG1
  size_t temp;
  fprintf (stdout, "  \n\t ucons_add dimensions: \n");
  for (size_t temp = 0; temp < dimchange->intdim; temp++)
    {
      fprintf (stdout, "%zu", dimchange->dim[temp]);
    }
  fprintf (stdout, "\n");
  ucons_fprint (stdout, pr->man, a, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  size_t i;
  size_t addi = dimchange->intdim;
  size_t nbseg = dimchange->realdim;
  ucons_t *r = ucons_alloc_internal (pr, a->datadim + addi, a->segmentdim + nbseg);

  ap_dimchange_t dimadd;
  dimadd.realdim = 0;
  dimadd.intdim = addi + 2 * nbseg;
  dimadd.dim = (ap_dim_t *) malloc (((2 * nbseg) + addi) * sizeof (ap_dim_t));
  for (i = 0; i < addi; i++)
    dimadd.dim[i] = dimchange->dim[i];

  for (i = 0; i < nbseg; i++)
    {
      dimadd.dim[addi + i] = a->datadim + a->segmentdim;
    }
  for (i = 0; i < nbseg; i++)
    {
      dimadd.dim[addi + nbseg + i] = a->datadim + 2 * a->segmentdim;
    }

  r->econs =
          ap_abstract0_add_dimensions (pr->man_dcons, false, a->econs, &dimadd,
                                       project);

  pattern_t *s, *rt = NULL;
  size_t u_seg, e_seg;
  size_t intdim, nr_y;
  unsigned keylen;



  for (s = a->udcons; s != NULL; s = s->hh.next)
    {

      nr_y = pr->PI[s->key.type].nr_y;
      u_seg = pr->PI[s->key.type].u_seg;
      e_seg = pr->PI[s->key.type].e_seg;
      unsigned keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);


      HASH_FIND (hh, r->udcons, &s->key, keylen, rt);
      if (rt)
        {
          rt->dcons = ap_abstract0_add_dimensions (pr->man_dcons, destructive, s->dcons, &dimadd,
                                                   project);
        }
      else
        {
          size_t ii;

          checked_malloc (rt, pattern_t, 1, sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t), return NULL;);
          memset (rt, 0, sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t));

          rt->dcons = ap_abstract0_add_dimensions (pr->man_dcons, destructive, s->dcons, &dimadd,
                                                   project);

          rt->key.type = s->key.type;

          for (ii = 0; ii < (u_seg + e_seg); ii++)
            rt->key.segments[ii] = s->key.segments[ii];


          HASH_ADD (hh, r->udcons, key, keylen, rt);
        }

    }
  ap_dimchange_clear (&dimadd);

  r->n2p = pattern_key_set_copy (pr, a->n2p, a->segmentdim);
  checked_realloc (r->n2p, pattern_key_set_t, (a->segmentdim + nbseg),
                   sizeof (pattern_key_set_t), return NULL;);

  for (size_t i = 0; i < nbseg; i++)
    {
      r->n2p[a->segmentdim + i].p = NULL;
      r->n2p[a->segmentdim + i].size = 0;
    }
  //	size_t dim = a->segmentdim+1;

  if (destructive)
    ucons_free_internal (pr, a);
#ifndef NDEBUG1
  fprintf (stdout, "  \n\t ucons_add dimensions returns: \n");
  for (size_t temp = 0; temp < dimchange->intdim; temp++)
    {
      fprintf (stdout, "%zu", dimchange->dim[temp]);
    }
  fprintf (stdout, "\n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif

  return r;
}

ucons_t*
add_pattern_n2p (ucons_internal_t *pr, ucons_t *r, pattern_key_t * key)
{

#ifndef NDEBUG1
  printf ("\n add_pattern_n2p \n");
  //	ucons_fprint(stdout,pr->man,r,NULL);
  printf ("\n");
  printf ("key->type %d \n", key->type);
  fflush (stdout);
  fprintf (stdout, "\n");
#endif

  if (!key->segments)
    {
#ifndef NDEBUG1
      fprintf (stderr, "@@@@ add_pattern_n2p: returns NULL (!key->segments).\n");
#endif
      return NULL;
    }
  if (!r)
    {
#ifndef NDEBUG1
      fprintf (stderr, "@@@@ add_pattern_n2p: returns NULL (!r).\n");
#endif
      return NULL;
    }
  if (!r->n2p)
    {
#ifndef NDEBUG1
      fprintf (stderr, "@@@@ add_pattern_n2p: returns NULL (!r->n2p).\n");
#endif
      return NULL;
    }

  size_t u_seg, e_seg, dim, size;

  size_t i, j;

  u_seg = pr->PI[key->type].u_seg;
  e_seg = pr->PI[key->type].e_seg;

  for (i = 0; i < (u_seg + e_seg); i++)
    {
      if (!(key->segments[i] < r->segmentdim)) return NULL;
      dim = key->segments[i];

      size = r->n2p[dim].size;
      checked_realloc (r->n2p[dim].p, pattern_key_t*,
                       (r->n2p[dim].size + 1),
                       sizeof (pattern_key_t *),
                       fprintf (stdout, "realloc.\n"); return NULL;);
      checked_malloc (r->n2p[dim].p[size], pattern_key_t,
                      1,
                      sizeof (pattern_key_t)+(u_seg + e_seg) * sizeof (size_t),
                      fprintf (stdout, "malloc.\n"); return NULL;);

      r->n2p[dim].size = r->n2p[dim].size + 1;
      r->n2p[dim].p[r->n2p[dim].size - 1]->type = key->type;

      for (j = 0; j < (u_seg + e_seg); j++)
        {
          size_t a = key->segments[j];
          r->n2p[dim].p[r->n2p[dim].size - 1]->segments[j] = key->segments[j];
        }
    }
#ifndef NDEBUG1
  printf ("\n add_pattern_n2p returns \n");
  //ucons_fprint(stdout,pr->man,r,NULL);
  fprintf (stdout, "\n");
#endif

  return r;

}

ucons_t*
remove_pattern_n2p (ucons_internal_t *pr, ucons_t *r, pattern_key_t * key)
{
#ifndef NDEBUG1
  fprintf (stdout, "\n remove_pattern_n2p \n");
  //ucons_fprint(stdout,pr->man,r,NULL);
  fprintf (stdout, "\n");
  fprintf (stdout, "key->type %zu \n", key->type);
  fprintf (stdout, "key->segm[0] %zu \n", key->segments[0]);
  fflush (stdout);
  fprintf (stdout, "\n");

#endif

  if (!r) return NULL;
  if (!r->n2p) return NULL;

  size_t u_seg, e_seg, dim, size, type;

  size_t i, j;

  u_seg = pr->PI[key->type].u_seg;
  e_seg = pr->PI[key->type].e_seg;
  type = key->type;

  for (i = 0; i < (u_seg + e_seg); i++)
    {
      if (!(key->segments[i] < r->segmentdim)) return NULL;
      dim = key->segments[i];

      size = r->n2p[dim].size;
      size_t u_seg_ctr, e_seg_ctr, type_ctr;
      /* parcurg tabela lui dim si sterg key (daca) de unde apare */
      if (size > 0)
        {
          for (size_t curr = 0; curr < size; curr++)
            {
              type_ctr = r->n2p[dim].p[curr]->type;
              u_seg_ctr = pr->PI[type_ctr].u_seg;
              e_seg_ctr = pr->PI[type_ctr].e_seg;
              if (type_ctr == type && u_seg_ctr == u_seg && e_seg_ctr == e_seg)
                {
                  if (size == 1)
                    {
                      free (r->n2p[dim].p);
                      r->n2p[dim].p = NULL;
                      r->n2p[dim].size = 0;
                    }
                  else
                    {
                      free (r->n2p[dim].p[curr]);
                      r->n2p[dim].p[curr] = NULL;
                      //		because of remove and all really remove the pattern from the list
                      for (size_t cpi = curr; cpi < size - 1; cpi++)
                        {
                          r->n2p[dim].p[cpi] = r->n2p[dim].p[cpi + 1];
                        }
                      r->n2p[dim].size -= 1;
                      checked_realloc (r->n2p[dim].p, pattern_key_t*, (r->n2p[dim].size), sizeof (pattern_key_t *), return NULL;);
                      curr = size;

                    }
                }//end pattern found
            }//end checking the patterns of dim
        }

    }
#ifndef NDEBUG1
  fprintf (stdout, "\n remove_pattern_n2p returns \n");
  //ucons_fprint(stdout,pr->man,r,NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  return r;
}

ucons_t *
ucons_remove_dimensions (ap_manager_t * man,
                         bool destructive, ucons_t * a,
                         ap_dimchange_t * dimchange)
{
  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_REMOVE_DIMENSIONS, 0);

  arg_assert (dimchange, return NULL;);
  if (!a)
    return NULL;

#ifndef NDEBUG1
  fprintf (stdout, "  \n\t ucons_remove: \n");
  for (size_t temp = 0; temp < dimchange->intdim; temp++)
    {
      fprintf (stdout, "%zu", dimchange->dim[temp]);
    }
  printf ("\n");
  ucons_fprint (stdout, pr->man, a, NULL);

  fprintf (stdout, "\n");
  fflush (stdout);
#endif

  // for one dimension realdim to remove, remove 2 dimensions
  ucons_t *r =
          ucons_alloc_internal (pr, a->datadim - dimchange->intdim, a->segmentdim - dimchange->realdim);

  size_t i, j;
  ap_dimchange_t dimrm;
  dimrm.realdim = 0;

  dimrm.intdim = dimchange->intdim + 2 * dimchange->realdim;
  dimrm.dim = (ap_dim_t *) malloc (dimrm.intdim * sizeof (ap_dim_t));

  for (i = 0; i < dimchange->intdim; i++)
    dimrm.dim[i] = dimchange->dim[i];

  for (i = 0; i < dimchange->realdim; i++)
    {
      dimrm.dim[dimchange->intdim + 2 * i + 0] = dimchange->dim[dimchange->intdim + i];
      dimrm.dim[dimchange->intdim + 2 * i + 1] = a->segmentdim + dimchange->dim[dimchange->intdim + i];
    }
  // dimrm shall be sorted!
  shape_dimchange_sort (&dimrm);
  r->econs =
          //ap_abstract0_remove_dimensions (pr->man_dcons, destructive, a->econs, &dimrm);
          ap_abstract0_remove_dimensions (pr->man_dcons, false, a->econs, &dimrm);

  size_t nr_y, u_seg, e_seg;
  pattern_t * s, *rt, *ra, *new_udcons;
  bool change;
  unsigned keylen;

  //#if defined (UCONS_DCONS_OCT_P11) || defined (UCONS_DCONS_POLY_P11) || defined (UCONS_DCONS_OCT_P12) || defined (UCONS_DCONS_POLY_P12)
  //	for(s=r->udcons;s!=NULL;s=s->hh.next){
  //
  //			change = false ;
  //
  //			u_seg=pr->PI[s->key.type].u_seg;
  //			e_seg=pr->PI[s->key.type].e_seg;
  //			//todo not really the good definition of removing when e_seg != 0
  //			keylen = (u_seg+e_seg)*sizeof(size_t)+sizeof(pattern_key_t);
  //
  //			checked_malloc(rt,pattern_t,1,sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t),return NULL;);
  //			memset(rt, 0, sizeof(pattern_t)+(u_seg+e_seg)*sizeof(size_t));
  //
  //			rt->key.type=s->key.type;
  //			rt->dcons = NULL;
  //			for(i=0; i<(u_seg+e_seg); i++){
  //				rt->key.segments[i] = s->key.segments[i];
  //				j=0;
  //				while(j<dimchange->realdim){
  //					if(s->key.segments[i] >= (dimchange->dim[j] - a->datadim))
  //					{
  //						rt->key.segments[i]++ ;
  //						change = true;
  //						j++;
  //						while((j<dimchange->realdim) && (dimchange->dim[j-1]+1 == dimchange->dim[j])  ){
  //							rt->key.segments[i]++ ;
  //							j++;
  //						}
  //					}
  //					else j++;
  //				}
  //			}
  //
  //			if(change){
  //
  //				HASH_FIND(hh,a->udcons,&rt->key, keylen,ra);
  //				if(ra)
  //					s->dcons = ap_abstract0_remove_dimensions(pr->man_dcons,destructive,ra->dcons,&dimrm);
  //				else
  //					s->dcons = ap_abstract0_bottom(pr->man_dcons,r->datadim,r->segmentdim);
  //				free(rt);
  //			}
  //			else {
  //				HASH_FIND(hh,a->udcons,&s->key, keylen,ra);
  //				if(ra)
  //					s->dcons = ap_abstract0_remove_dimensions(pr->man_dcons,destructive,ra->dcons,&dimrm);
  //				else
  //					s->dcons = ap_abstract0_bottom(pr->man_dcons,r->datadim,r->segmentdim);
  //				free(rt);
  //			}
  //
  //		}
  //#else


  for (s = a->udcons; s != NULL; s = s->hh.next)
    {

      u_seg = pr->PI[s->key.type].u_seg;
      e_seg = pr->PI[s->key.type].e_seg;

      size_t new_segments_size = 0;

      for (size_t i = 0; i < (u_seg); i++)
        {
          j = 0;
          bool del = false;
          while (j < dimchange->realdim)
            {
              if (s->key.segments[i] == (dimchange->dim[j] - a->datadim))
                {
                  j++;
                  del = true;
                }
              else j++;
            }
          if (del == false) new_segments_size++;
        }
      /*todo if pattern y... y=l[x] and remove x then look for z with l[z] = l[x]
if it doesn't exist remove otherwise modify
       */
      size_t *new_e_seg;
      new_e_seg = (size_t*) malloc (e_seg * sizeof (size_t));
      for (size_t i = 0; i < e_seg; i++)
        new_e_seg[i] = 0;

      for (size_t i = u_seg; i < (u_seg + e_seg); i++)
        {
          bool del = false;
          j = 0;
          while (j < dimchange->realdim)
            {
              if (s->key.segments[i] == (dimchange->dim[j] - a->datadim))
                {
                  j++;
                  del = true;
                }
              else j++;
            }
          if (del == false)
            {
              new_segments_size++;
              new_e_seg[i - u_seg] = s->key.segments[i];
            }
          else
            {
              // look for z with l[x] = l[z]
              bool found_lz = false;
              size_t dim = s->key.segments[i];
              //			for(size_t j=1; j<r->segmentdim && found_lz==false; j++){
              //				if(j!=s->key.segments[i]){
              //					ap_linexpr0_t* linexpr =
              //							ap_linexpr0_alloc (AP_LINEXPR_DENSE,  a->datadim + 2 * a->segmentdim);
              //					//linexpr.scalar = NULL;
              //					ap_linexpr0_set_coeff_scalar_int (linexpr,
              //							a->datadim + a->segmentdim + dim, 1);
              //					ap_linexpr0_set_coeff_scalar_int (linexpr,
              //							a->datadim + a->segmentdim + j, -1);
              //					ap_linexpr0_set_cst_scalar_int (linexpr, 0);
              //
              //					ap_lincons0_t cons = ap_lincons0_make (AP_CONS_EQ, linexpr, NULL);
              //					if (ap_abstract0_sat_lincons(pr->man_dcons, a->econs, &cons )){
              //						found_lz=true;
              //						new_e_seg[i] = j;
              //					}
              //					ap_lincons0_clear (&cons);
              //				}
              //			}
              if (found_lz)
                {
                  new_segments_size++;
                }
            }
        }

      /* if a pattern is build over some node to remove then the pattern is removed */
      if (new_segments_size == (u_seg + e_seg))
        {
          keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);
          checked_malloc (rt, pattern_t, 1, sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t), return NULL;);
          memset (rt, 0, sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t));
          rt->key.type = s->key.type;
          rt->dcons = NULL;

          for (i = 0; i < (u_seg + e_seg); i++)
            {
              if (i < u_seg) rt->key.segments[i] = s->key.segments[i];
              else rt->key.segments[i] = new_e_seg[i - u_seg];
              j = 0;
              while (j < dimchange->realdim)
                {
                  if ((i < u_seg && s->key.segments[i] >= (dimchange->dim[j] - a->datadim)) ||
                      (i >= u_seg && new_e_seg[i - u_seg] >= (dimchange->dim[j] - a->datadim)))
                    {
                      rt->key.segments[i]--;
                      change = true;
                      j++;
                      while ((j < dimchange->realdim) && (dimchange->dim[j - 1] + 1 == dimchange->dim[j]))
                        {
                          rt->key.segments[i]--;
                          j++;
                        }
                    }
                  else j++;
                }
            }

#ifndef NDEBUG1
          printf ("  \n\t remove: \n");
          for (size_t temp = 0; temp < dimchange->intdim; temp++)
            {
              printf ("%zu", dimchange->dim[temp]);
            }
          printf ("  \n\t pattern looked in the list: \n");
          //printf("\n");
          for (size_t temp = 0; temp < (u_seg + e_seg); temp++)
            {
              printf ("rt->key.segments[%zu]=%zu ", temp, rt->key.segments[temp]);
            }

          printf ("\n");
#endif

          pattern_t *r_pattern;
          HASH_FIND (hh, r->udcons, &rt->key, keylen, r_pattern);
          if (r_pattern)
            {
              r_pattern->dcons = ap_abstract0_remove_dimensions (pr->man_dcons, 
                                                                 false, s->dcons, &dimrm);
            }
          else
            {

              checked_malloc (r_pattern, pattern_t, 1, (sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t)), return NULL;);
              memset (r_pattern, 0, (sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t)));

              r_pattern->dcons = ap_abstract0_remove_dimensions (pr->man_dcons, 
                                                                 false, s->dcons, &dimrm);
              r_pattern->key.type = rt->key.type;

              for (size_t i = 0; i < (u_seg + e_seg); i++)
                {
                  r_pattern->key.segments[i] = rt->key.segments[i];
                }

              HASH_ADD (hh, r->udcons, key, keylen, r_pattern);
              r = add_pattern_n2p (pr, r, &r_pattern->key);
            }
          free (rt);

        }
      //todo add the situation corresponding to pattern removal in case of P11 P12 when new_segments_size!= u_seg +e_seg
    }


  //#endif


  //#ifndef NDEBUG1
  //	fprintf(stdout,"\n remove dimensions returns \n");
  //	//ucons_fprint(stdout,man,r,NULL);
  //	fprintf(stdout,"\n");
  //	fflush(stdout);
  //#endif

  ap_dimchange_clear (&dimrm);
  if (destructive)
    ucons_free_internal (pr, a);

#ifndef NDEBUG1
  fprintf (stdout, "\n ucons_remove dimensions returns: \n");
  ucons_fprint (stdout, man, r, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  return r;
}

/* TODO: priority 0 */
ucons_t *
ucons_permute_dimensions (ap_manager_t * man,
                          bool destructive, ucons_t * a,
                          ap_dimperm_t * perm)
{
  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_ADD_DIMENSIONS, 0);
  if (!a)
    return NULL;
  assert (perm && perm->size == (a->datadim + a->segmentdim));

  if (ucons_is_bottom (pr->man, a))
    return a;


#ifndef NDEBUG1
  fprintf (stdout, "  \n\t ucons_permute: \n");
  fprintf (stdout, "\n");
  ucons_fprint (stdout, pr->man, a, NULL);
  printf ("\n");
  fflush (stdout);
#endif

  ucons_t *r;
  /*
   * Remove node dimensions mapped to a->datadim (except a->datadim itself) and build a new permutation
   */
  ap_dimchange_t dimrm;
  ap_dimperm_t newperm;
  ap_dimperm_t new_node_perm;
  size_t i, j, k, rmsize, newsize;
  dimrm.intdim = 0;
  dimrm.realdim = perm->size;
  dimrm.dim = (ap_dim_t *) malloc (perm->size * sizeof (ap_dim_t));
  newperm.size = perm->size;
  newperm.dim = (ap_dim_t *) malloc (perm->size * sizeof (ap_dim_t));
  /* new_node_perm is the new permutation of nodes */
  new_node_perm.size = perm->size - a->datadim;
  new_node_perm.dim = (ap_dim_t *) malloc ((perm->size - a->datadim) * sizeof (ap_dim_t));

  for (i = 0, j = 0, rmsize = 0; i < perm->size; i++)
    if (i > a->datadim && perm->dim[i] == a->datadim)
      {
        dimrm.dim[rmsize] = i;
        rmsize++;
      }
    else
      {
        newperm.dim[j] = perm->dim[i];
        if (j >= a->datadim) new_node_perm.dim[j - a->datadim] = perm->dim[i] - a->datadim;
        j++;
      }
  if (rmsize > 0)
    {
      dimrm.intdim = 0;
      dimrm.realdim = rmsize;
      dimrm.dim =
              (ap_dim_t *) realloc (dimrm.dim, rmsize * sizeof (ap_dim_t));

      r = ucons_remove_dimensions (man, false, a, &dimrm);
      pr = ucons_init_from_manager (man, AP_FUNID_PERMUTE_DIMENSIONS, 0);
      newperm.size = perm->size - rmsize;
      newperm.dim =
              (ap_dim_t *) realloc (newperm.dim, newperm.size * sizeof (ap_dim_t));

      new_node_perm.size = new_node_perm.size - rmsize;
      new_node_perm.dim =
              (ap_dim_t *) realloc (new_node_perm.dim, new_node_perm.size * sizeof (ap_dim_t));

    }
  else
    r = ucons_copy_internal (pr, a);
  free (dimrm.dim);

#ifndef NDEBUG2
  fprintf (stdout, "  \n\t actually ucons_permute: \n");
  fprintf (stdout, "\n");
  ucons_fprint (stdout, pr->man, r, NULL);

  fprintf (stdout, "\n");
  fflush (stdout);
#endif


  /*
   * Build permutations corresponding to each set of constraints.
   */
  ap_dimperm_t consperm;
  size_t newsegmsize = (newperm.size - a->datadim);

  consperm.size = a->datadim + 2 * newsegmsize;
  consperm.dim = (ap_dim_t *) malloc (consperm.size * sizeof (ap_dim_t));
  ap_dimperm_set_id (&consperm);
  for (i = 0; i < newsegmsize; i++)
    {
      consperm.dim[a->datadim + i] = newperm.dim[a->datadim + i];
      consperm.dim[a->datadim + newsegmsize + i] =
              newsegmsize + newperm.dim[a->datadim + i];
    }

  r->econs = ap_abstract0_permute_dimensions (pr->man_dcons, true, r->econs,
                                              &consperm);

  free (consperm.dim);


  size_t nr_y, u_seg, e_seg;
  pattern_t * s, *ra, *rt;
  unsigned keylen;

  pattern_t * perm_r = NULL;

  for (s = r->udcons; s != NULL; s = s->hh.next)
    {

      ra = NULL;
      u_seg = pr->PI[s->key.type].u_seg;
      e_seg = pr->PI[s->key.type].e_seg;

      checked_malloc (ra, pattern_t, 1, sizeof (pattern_t)+ (u_seg + e_seg) * sizeof (size_t), return NULL;);
      memset (ra, 0, sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t));


      ra->dcons = (s->dcons) ? ap_abstract0_copy (pr->man_dcons, s->dcons) : NULL;

      ra->key.type = s->key.type;
      /* ra->key.segments[] = perm ( s->key.segments[] ) */
      for (i = 0; i < (u_seg + e_seg); i++)
        ra->key.segments[i] = new_node_perm.dim[s->key.segments[i]];


      /*TODO general case: sort segments and apply the resulting permutation on universals */

      bool sort = true;
      //#if defined (UCONS_DCONS_OCT_P21) || defined (UCONS_DCONS_POLY_P21)
      for (size_t ii = 1; ii < u_seg; ii++)
        {
          size_t jj = 0;
          while (jj != ii && ra->key.segments[jj] <= ra->key.segments[ii])
            jj++;
          if (jj < ii)
            {
              size_t d = ra->key.segments[ii];
              size_t tmp;
              size_t kk;
              for (tmp = ii; tmp > jj; tmp--)
                {
                  kk = tmp - 1;
                  ra->key.segments[kk + 1] = ra->key.segments[kk];
                }
              //				for (size_t  kk = ii - 1; kk >= jj; kk--){
              //					ra->key.segments[kk + 1] = ra->key.segments[kk];
              //				}
              ra->key.segments[jj] = d;
              sort = false;
            }
        }



      keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);
      nr_y = pr->PI[s->key.type].nr_y;

      newsegmsize = (newperm.size - a->datadim);
      consperm.size = a->datadim + 2 * newsegmsize + 2 * nr_y;
      consperm.dim = (ap_dim_t *) malloc (consperm.size * sizeof (ap_dim_t));
      ap_dimperm_set_id (&consperm);
      for (i = 0; i < newsegmsize; i++)
        {
          consperm.dim[a->datadim + i] = newperm.dim[a->datadim + i];
          consperm.dim[a->datadim + newsegmsize + i] =
                  newsegmsize + newperm.dim[a->datadim + i];
        }

      /*TODO general case: sort segments and apply the resulting permutation on universals */

      //#if defined (UCONS_DCONS_OCT_P21) || defined (UCONS_DCONS_POLY_P21)
      if (sort == false)
        {
          for (i = 0; i < u_seg; i++)
            {
              consperm.dim[a->datadim + 2 * newsegmsize + 2 * i] = a->datadim + 2 * newsegmsize + (2 * i) + 1;
              consperm.dim[a->datadim + 2 * newsegmsize + 2 * i + 1] = a->datadim + 2 * newsegmsize + (2 * i);
            }

          if (pr->PI[ra->key.type].kind == pattern_2_1_lx)
            ra->key.type = get_pattern_type (pr, u_seg, e_seg, nr_y, pattern_2_1_mlx);
          else if (pr->PI[ra->key.type].kind == pattern_2_1_mlx)
            ra->key.type = get_pattern_type (pr, u_seg, e_seg, nr_y, pattern_2_1_lx);

        }

      //#endif
      ra->dcons = ap_abstract0_permute_dimensions (pr->man_dcons, true, ra->dcons,
                                                   &consperm);
      HASH_ADD (hh, perm_r, key, keylen, ra);
      free (consperm.dim);

    }


  while (r->udcons)
    {
      s = r->udcons;
      HASH_DEL (r->udcons, s);
      free (s);
    }
  r->udcons = NULL;


  for (s = perm_r; s != NULL; s = s->hh.next)
    {

      u_seg = pr->PI[s->key.type].u_seg;
      e_seg = pr->PI[s->key.type].e_seg;

      keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);

      checked_malloc (rt, pattern_t, 1, sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t), return NULL;);
      memset (rt, 0, sizeof (pattern_t)+(u_seg + e_seg) * sizeof (size_t));

      rt->key.type = s->key.type;
      for (i = 0; i < (u_seg + e_seg); i++)
        rt->key.segments[i] = s->key.segments[i];
      rt->dcons = ap_abstract0_copy (pr->man_dcons, s->dcons);

      HASH_ADD (hh, r->udcons, key, keylen, rt);
    }

  while (perm_r)
    {
      s = perm_r;
      HASH_DEL (perm_r, s);
      free (s);
    }




#ifndef NDEBUG2
  printf ("  \n\t node to pattern : \n");
  printf ("\n");
  ucons_fprint (stdout, pr->man, r, NULL);
  printf ("\n");
#endif


  pattern_key_set_t * old_n2p;
  old_n2p = pattern_key_set_copy (pr, r->n2p, r->segmentdim);

  /* clean old list of patterns */
  for (size_t i = 0; i < r->segmentdim; i++)
    {
      r->n2p[i].size = 0;
      free (r->n2p[i].p);
      r->n2p[i].p = NULL;
    }
  /*rebuild list with for new entrances*/
  for (size_t i = 0; i < r->segmentdim; i++)
    {
      size_t new_i = new_node_perm.dim[i];

      size_t size_new_i = old_n2p[i].size;
      checked_malloc (r->n2p[new_i].p, pattern_key_t*, size_new_i, sizeof (pattern_key_t), return NULL;);

      r->n2p[new_i].size = size_new_i;

      for (size_t j = 0; j < size_new_i; j++)
        {
          u_seg = pr->PI[old_n2p[i].p[j]->type].u_seg;
          e_seg = pr->PI[old_n2p[i].p[j]->type].e_seg;
          checked_malloc (r->n2p[new_i].p[j], pattern_key_t, 1,
                          sizeof (pattern_key_t)+ (u_seg + e_seg) * sizeof (size_t), return NULL;);

          r->n2p[new_i].p[j]->type = old_n2p[i].p[j]->type;

          for (size_t k = 0; k < (u_seg + e_seg); k++)
            {
              size_t i2 = old_n2p[i].p[j]->segments[k];
              size_t i3 = new_node_perm.dim[i2];
              r->n2p[new_i].p[j]->segments[k] = i3;
            }
          /* sorting universal quantif segments */
          for (size_t ii = 1; ii < u_seg; ii++)
            {
              size_t jj = 0;
              while (jj != ii && r->n2p[new_i].p[j]->segments[jj] <= r->n2p[new_i].p[j]->segments[ii])
                jj++;
              if (jj < ii)
                {
                  size_t d = r->n2p[new_i].p[j]->segments[ii];
                  size_t tmp;
                  size_t kk;
                  for (tmp = ii; tmp > jj; tmp--)
                    {
                      kk = tmp - 1;
                      r->n2p[new_i].p[j]->segments[kk + 1] = r->n2p[new_i].p[j]->segments[kk];
                    }
                  //							for (size_t  kk = ii - 1; kk >= jj; kk--)
                  //								r->n2p[new_i].p[j]->segments[kk + 1] = r->n2p[new_i].p[j]->segments[kk];
                  r->n2p[new_i].p[j]->segments[jj] = d;
                }
            }
          //todo sort e_seg also

        }

    }
  /*clean old_n2p*/
  for (size_t i = 0; i < r->segmentdim; i++)
    {
      old_n2p[i].size = 0;
      free (old_n2p[i].p);
      old_n2p[i].p = NULL;
    }
  free (old_n2p);
  //if(new_node_perm.dim) free(new_node_perm.dim);



  //	size_t i1,i2,i3;
  //	int jj,kk,ii;
  //	size_t d;
  //
  //	for(i=0;i<r->segmentdim;i++){
  //		i1 = perm_n2p[new_node_perm.dim[i]].size;
  //		r->n2p[i].size = perm_n2p[new_node_perm.dim[i]].size;
  //		free (r->n2p[i].p);
  //		r->n2p[i].p = NULL;
  //		checked_malloc(r->n2p[i].p, pattern_key_t*,r->n2p[i].size,sizeof(pattern_key_t),return NULL;);
  //		for(j=0;j<r->n2p[i].size;j++){
  //			u_seg = pr->PI[perm_n2p[i].p[j]->type].u_seg;
  //			//e_seg = pr->PI[perm_n2p[i].p[j]->type].e_seg;
  //
  //			checked_malloc(r->n2p[i].p[j], pattern_key_t,1,sizeof(pattern_key_t)+ (u_seg)*sizeof(size_t),return NULL ;);
  //
  //			r->n2p[i].p[j]->type = perm_n2p[i].p[j]->type;
  //			for(k=0; k < (u_seg); k++){
  //				i2 = r->n2p[i].p[j]->segments[k] ;
  //				i3 = new_node_perm.dim[perm_n2p[new_node_perm.dim[i]].p[j]->segments[k]];
  //				r->n2p[i].p[j]->segments[k]=new_node_perm.dim[perm_n2p[new_node_perm.dim[i]].p[j]->segments[k]];
  //			}
  //
  //			for (ii = 1; ii < (int)u_seg; ii++)
  //			{
  //				jj = 0;
  //				while (jj != ii && r->n2p[i].p[j]->segments[jj] <= r->n2p[i].p[j]->segments[ii])
  //					jj++;
  //				if (jj < ii)
  //				{
  //					d = r->n2p[i].p[j]->segments[ii];
  //					for (kk = ii - 1; kk >= jj; kk--)
  //						r->n2p[i].p[j]->segments[kk + 1] = r->n2p[i].p[j]->segments[kk];
  //					r->n2p[i].p[j]->segments[jj] = d;
  //				}
  //			}
  //
  //		}
  //	}


  if (destructive)
    ucons_free_internal (pr, a);

#ifndef NDEBUG1
  fprintf (stdout, "  \n\t ucons_permute returns : \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  return r;
}


/* ============================================================ */
/* Expansion and folding of dimensions */

/* ============================================================ */


ucons_t *
ucons_singleton (ucons_internal_t * pr, bool destructive, ucons_t * a,
                 ap_dim_t dim)
{

#ifndef NDEBUG1
  fprintf (stdout, "\n singleton \n");
  ucons_fprint (stdout, pr->man, a, NULL);
  fprintf (stdout, "\n");
#endif

  ucons_t *r = ucons_copy_internal (pr, a);
  /*
   * intersect with constraints l[dim] = 1 and r->udcons[]->dcons = 0
   */

  ap_lincons0_array_t arr = ap_lincons0_array_make (1);
  // l[dim]-1 == 0
  arr.p[0].constyp = AP_CONS_EQ;
  arr.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, a->datadim + 2 * a->segmentdim);
  arr.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                    a->datadim + a->segmentdim + dim, 1);
  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);

  r->econs =
          ap_abstract0_meet_lincons_array (pr->man_dcons, false, a->econs, &arr);
  ap_lincons0_array_clear (&arr);


  //	if (ap_abstract0_is_bottom(pr->man_dcons,r->econs))
  //		return ucons_bottom (pr->man, a->datadim, a->segmentdim);

  /* TODO: add modifications on r->udcons. */

  size_t i, nr_y, u_seg, e_seg, total_dim;
  pattern_t * s;
  bool found;


  for (s = r->udcons; s != NULL;)
    {
      u_seg = pr->PI[s->key.type].u_seg;
      e_seg = pr->PI[s->key.type].e_seg;

      found = false;
      for (i = 0; i < (u_seg); i++)
        if (s->key.segments[i] == dim) found = true;

      total_dim = a->datadim + 2 * a->segmentdim + 2 * pr->PI[s->key.type].nr_y;

      arr = ap_lincons0_array_make (1);
      // l[dim]-1 == 0
      arr.p[0].constyp = AP_CONS_EQ;
      arr.p[0].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, total_dim);
      arr.p[0].scalar = NULL;
      ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                        a->datadim + a->segmentdim + dim, 1);
      ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);

      if (found)
        {
          //#if defined (UCONS_DCONS_OCT_P12) || defined (UCONS_DCONS_POLY_P12)	|| defined (UCONS_DCONS_OCT_P11) || defined (UCONS_DCONS_POLY_P11)
          //			s->dcons = ap_abstract0_bottom(pr->man_dcons,total_dim,0);
          //			s->dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, s->dcons, &arr);
          //			s = s->hh.next;
          //#else
          remove_pattern_n2p (pr, r, &s->key);
          pattern_t *t = s;
          s = s->hh.next;
          HASH_DEL (r->udcons, t);
          free (t);
          //#endif
        }
      else
        {
          s->dcons =
                  ap_abstract0_meet_lincons_array (pr->man_dcons, true, s->dcons, &arr);
          if (ap_abstract0_is_bottom (pr->man_dcons, s->dcons))
            {
              remove_pattern_n2p (pr, r, &s->key);
              pattern_t *t = s;
              s = s->hh.next;
              HASH_DEL (r->udcons, t);
              free (t);
            }
          else
            {
              s = s->hh.next;
            }
        }
      ap_lincons0_array_clear (&arr);
    }


  if (destructive)
    ucons_free_internal (pr, a);

#ifndef NDEBUG1
  fprintf (stdout, "\n singleton returns: \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, "\n");
#endif

  return r;
}

bool
test_pattern_sat (ucons_internal_t * pr, ucons_t *r, pattern_key_t *p, int s)
{

  size_t u_seg, e_seg, nr_yi;
  size_t i, li;
  bool flag;
  ap_lincons0_t cons;
  ap_linexpr0_t *expr;

  u_seg = pr->PI[p->type].u_seg;
  e_seg = pr->PI[p->type].e_seg;

  flag = true;
  if (s == 0)
    {
      for (i = 0; i < (u_seg) && flag == true; i++)
        {
          // l(segments[i]) - pr->PI[p->type].uvar[i].size > 0
          li = r->datadim + r->segmentdim + p->segments[i];
          //	TODO nr_yi = pr->PI[p->type].uvar[i].size; sa ii transform in cei obligatorii cum pt <= e numai unul necesar
          nr_yi = 1;
          expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
          ap_linexpr0_set_coeff_scalar_int (expr, li, 1);
          ap_linexpr0_set_cst_scalar_int (expr, -nr_yi - 1);
          cons = ap_lincons0_make (AP_CONS_SUPEQ, expr, NULL);



#ifndef NDEBUG1
          fprintf (stdout, "\t test if : \n");
          if (r->econs) ap_abstract0_fprint (stdout, pr->man_dcons, r->econs, NULL);
          fprintf (stdout, "\n");
          fprintf (stdout, "\t satisfies : \n");
          ap_lincons0_fprint (stdout, &cons, NULL);
          fflush (stdout);
#endif

          flag = ap_abstract0_sat_lincons (pr->man_dcons, r->econs, &cons);
          ap_lincons0_clear (&cons);

        }
      return flag;
    }
  else if (s == 1)
    {
      li = r->datadim + r->segmentdim + p->segments[0];
      expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
      ap_linexpr0_set_coeff_scalar_int (expr, li, 1);
      ap_linexpr0_set_cst_scalar_int (expr, -1);
      cons = ap_lincons0_make (AP_CONS_EQ, expr, NULL);
      flag = ap_abstract0_sat_lincons (pr->man_dcons, r->econs, &cons);
      ap_lincons0_clear (&cons);

      return flag;
    }
  else
    {
      if (s != 3)
        {
          li = r->datadim + r->segmentdim + p->segments[0];
          expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
          ap_linexpr0_set_coeff_scalar_int (expr, li, 1);
          ap_linexpr0_set_cst_scalar_int (expr, -1);
          cons = ap_lincons0_make (AP_CONS_SUP, expr, NULL);
          flag = ap_abstract0_sat_lincons (pr->man_dcons, r->econs, &cons);
          ap_lincons0_clear (&cons);

          return flag;
        }
      else
        { //s == 3
          // -l(segments[i]) + pr->PI[p->type].uvar[i].size >= 0

          for (i = 0; i < u_seg && flag == true; i++)
            {
              li = r->datadim + r->segmentdim + p->segments[i];
              nr_yi = pr->PI[p->type].uvar[i].size;

              expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
              ap_linexpr0_set_coeff_scalar_int (expr, li, -1);
              ap_linexpr0_set_cst_scalar_int (expr, nr_yi);
              cons = ap_lincons0_make (AP_CONS_SUPEQ, expr, NULL);


#ifndef NDEBUG1
              fprintf (stdout, "\t test if : \n");
              if (r->econs) ap_abstract0_fprint (stdout, pr->man_dcons, r->econs, NULL);
              fprintf (stdout, "\n");
              fprintf (stdout, "\t satisfies : \n");
              ap_lincons0_fprint (stdout, &cons, NULL);
              fflush (stdout);
#endif

              flag = ap_abstract0_sat_lincons (pr->man_dcons, r->econs, &cons);
              ap_lincons0_clear (&cons);
            }
          return !flag;
        }
    }
}

/**
 * Split the node @code{dim} in @code{a}.
 * 
 * @param pr
 * @param destructive   true if a shall be modified
 * @param a
 * @param dim
 * @return 
 */
ucons_t *
ucons_split (ucons_internal_t * pr, bool destructive, ucons_t * a, ap_dim_t dim)
{
  size_t i, j, k;
  pattern_t *aux;
  ap_abstract0_t *n2_aux, *econs_aux;

  unsigned keylen;
  size_t u_seg, e_seg;

  ap_dim_t * tdim;
  size_t contor;
  size_t pos_y = 0; /* the dimension in r->udcons of y to be eliminated , data(y) = r->datadim + 2*r->segmentdim + nr_y + pos_y*/
  size_t nr_y; /* number of universally quantified variables */

  ap_dimchange_t* dimchange;
  ap_lincons0_array_t arr;

  bool pattern_sat;
  bool found_pos_y;

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_split: of %zu in a=(\n", dim);
  ucons_fprint (stdout, pr->man, a, NULL);
  fprintf (stdout, ")\n");
#endif
  /*
   * The node added by expansion has dimension a->segmdim
   * Then apply substitutions/assignments:
   * l[n1] == l[n2] + 1;
   * l[n1] := 1;
   * meet l[n2]>=1
   */
  // Add dimension a->segmentdim
  ap_dimchange_t dimadd;
  ap_dimchange_init (&dimadd, 0, 1);
  dimadd.dim = (ap_dim_t *) malloc (1 * sizeof (ap_dim_t));
  dimadd.dim[0] = a->datadim + a->segmentdim;

  ucons_t *r = ucons_add_dimensions (pr->man, false, a, &dimadd, false);
  ap_dimchange_clear (&dimadd);

  // Build statements and apply them
  ap_dim_t n1, n2, dn2, ln1, ln2;
  ap_linexpr0_t *expr;

  n1 = dim;
  n2 = a->segmentdim;
  dn2 = a->datadim + n2;

  /*TODO: add precision on universal constraints, when working with sets of patterns */


  /* adding the constraints on dn2 in r->econs from universal formulas */

  pattern_key_set_t n1_pat_set = a->n2p[n1];

  /* l[n2]=l[n1]-1  \land l[n1]=1 */


  /*
   * meet with l[n1] > 1 on existential plus universals 
   * (this is the precondition under which split is called)
   * it is true in the graph so it has to be stated in the conditions also
   *
   * */
  ///////*********************** r->econs
  ln1 = r->datadim + r->segmentdim + n1;
  expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
  ap_linexpr0_set_coeff_scalar_int (expr, ln1, 1);
  ap_linexpr0_set_cst_scalar_int (expr, -2);
  ap_lincons0_array_t arrr = ap_lincons0_array_make (1);
  arrr.p[0] = ap_lincons0_make (AP_CONS_SUPEQ, expr, NULL);
  // meet l[n1]-2 >= 0
  r->econs =
          ap_abstract0_meet_lincons_array (pr->man_dcons, true, r->econs, &arrr);
  ap_lincons0_array_clear (&arrr);
  ///////***********************r->udcons

  pattern_t * sp;

  for (sp = r->udcons; sp != NULL; sp = sp->hh.next)
    {

      //		u_seg = pr->PI[sp->key.type].u_seg;
      //		e_seg = pr->PI[sp->key.type].e_seg;
      nr_y = pr->PI[sp->key.type].nr_y;

      ln1 = r->datadim + r->segmentdim + n1;
      expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2 * nr_y);
      ap_linexpr0_set_coeff_scalar_int (expr, ln1, 1);
      ap_linexpr0_set_cst_scalar_int (expr, -2);
      arrr = ap_lincons0_array_make (1);
      arrr.p[0] = ap_lincons0_make (AP_CONS_SUPEQ, expr, NULL);

      if (!ap_abstract0_is_bottom (pr->man_dcons, sp->dcons))
        {
          sp->dcons =
                  ap_abstract0_meet_lincons_array (pr->man_dcons, true, sp->dcons, &arrr);
        }
      ap_lincons0_array_clear (&arrr);
    }

  ///////***********************
#ifndef NDEBUG2
  fprintf (stdout, "\n@@@@ ucons_split: actually r = (\n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, ")\n");
#endif

  if (ap_abstract0_is_bottom (pr->man_dcons, r->econs))
    {
      if (destructive)
        ucons_free_internal (pr, a);
      //		ucons_free_internal(pr,r);
      //		r = NULL;
      return r;
    }

  econs_aux = NULL;
  /* r->econs == add_dimension(a->econs) + l[n1]>2 */
  /* econs_aux = join n2_aux, pattern(y=1 y\in n1) ==> n2_aux
   * n2_aux is defined over y and d(y) where y is the pos_y quantified variable */
  /* computing r->econs = add_dimension(a->econs) \meet n2_aux */

  pattern_key_t *pattern_j;

  for (j = 0; j < n1_pat_set.size; j++)
    {


      if (n1_pat_set.p[j] != NULL)
        {
          /* for every pattern on the segment n1 do :
           * 1. compute the data constraint n2_aux to collect information on dn2 from universal formulas
           * 2. transfer the universal constraint on n2
           * 3. set to bottom the universal formula on n1
           * 4. l'[n2] = l[n1] - 1 /\  l'[n1] = 1 /\ define dn2 in r->econs and r->ucons
           */

          pattern_j = n1_pat_set.p[j];

#ifndef NDEBUG1
          fprintf (stdout, "\n@@@@ ucons_split: pattern: \n");
          //pattern_key_fprint (stdout, pr, pattern_j, NULL);
          fprintf (stdout, "\n type %d segment %d ",
                   pattern_j->type, pattern_j->segments[0]);
          fflush (stdout);
#endif

          u_seg = pr->PI[pattern_j->type].u_seg;
          e_seg = pr->PI[pattern_j->type].e_seg;
          keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);

          /* test pattern validity */
          nr_y = pr->PI[pattern_j->type].nr_y;
          pattern_sat = true;

          HASH_FIND (hh, r->udcons, pattern_j, keylen, aux);
          if (aux)
            {
              /*
               * temporary removed test_pattern_sat
               * because splitting DOES NOT consider sub-patterns
               */
              ap_dim_t n = pattern_j->segments[0];
              pattern_sat = is_pattern_inst (pr, r, pattern_j);
              //!test_singleton(pr->man_dcons,r->econs, r->datadim, r->segmentdim , n);
              //test_pattern_sat(pr,r,pattern_j,0);

              if (pattern_sat)
                {
                  pos_y = 0;
                  found_pos_y = false;
                  /* 1.reinforcing r->econs */
                  for (i = 0; i < u_seg && !found_pos_y; i++)
                    {
                      if (pattern_j->segments[i] == n1)
                        {
                          pos_y += 1;
                          found_pos_y = true;
                        }
                      else pos_y += pr->PI[pattern_j->type].uvar[i].size;
                    }

                  if (nr_y > 1)
                    {
                      /* if the pattern is different from \forall y */
                      checked_malloc (tdim, ap_dim_t, 2 * (nr_y - 1), sizeof (ap_dim_t), return NULL;);
                      /*TO DO does not take into consideration the cases when y1!=1 or y1=y2=1*/
                      contor = 0;
                      for (i = 0; i < 2 * nr_y; i++)
                        {
                          if (i != (pos_y - 1) && i != nr_y + (pos_y - 1))
                            {
                              tdim[contor] = r->datadim + 2 * r->segmentdim + i;
                              contor += 1;
                            }
                        }
                      n2_aux = ap_abstract0_copy (pr->man_dcons, aux->dcons);
                      //n2_aux = ap_abstract0_forget_array(pr->man_dcons,true,n2_aux,tdim,2*(nr_y-1),true);

                      dimchange = ap_dimchange_alloc (2 * (nr_y - 1), 0);
                      for (k = 0; k < 2 * (nr_y - 1); k++)
                        {
                          dimchange->dim[k] = tdim[k];
                        }
                      n2_aux = ap_abstract0_remove_dimensions (pr->man_dcons, true, n2_aux, dimchange);

                      free (tdim);
                      ap_dimchange_free (dimchange);
                    }
                  else
                    {
                      n2_aux = ap_abstract0_copy (pr->man_dcons, aux->dcons);
                    }


                  //#ifndef NDEBUG1
                  //	fprintf(stdout,"\t existential data collected for n2_aux: \n");
                  //	ap_abstract0_fprint(stdout,pr->man->dcons,n2_aux, NULL);
                  //	fprintf(stdout,"\n");
                  //#endif

                  /* n2_aux size is r->datadim + 2*r->segmentdim + 2 */
                  pos_y = r->datadim + 2 * r->segmentdim;
                  /* y=1 */
                  arr = ap_lincons0_array_make (2);
                  // y-1 == 0 where y denotes n2 in universal
                  arr.p[0].constyp = AP_CONS_EQ;
                  arr.p[0].linexpr0 =
                          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
                  arr.p[0].scalar = NULL;
                  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, pos_y, 1);
                  ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -1);

                  //  dn2 already added but unconstrained;
                  //  dn2 = pos_y + 1
                  arr.p[1].constyp = AP_CONS_EQ;
                  arr.p[1].linexpr0 =
                          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
                  arr.p[1].scalar = NULL;
                  ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, (pos_y + 1), -1);
                  ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dn2, 1);
                  ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);


                  n2_aux = ap_abstract0_meet_lincons_array (pr->man_dcons, false, n2_aux, &arr);
                  ap_lincons0_array_clear (&arr);

                  /* eliminate y and dy to get constraints on existential */
                  checked_malloc (tdim, ap_dim_t, 2, sizeof (ap_dim_t), return NULL;);
                  tdim[0] = pos_y;
                  tdim[1] = pos_y + 1;
                  n2_aux = ap_abstract0_forget_array (pr->man_dcons, false, n2_aux, tdim, 1, false);
                  //  remove dimension pos_y, pos_y+1
                  dimchange = ap_dimchange_alloc (2, 0);
                  dimchange->dim[0] = pos_y;
                  dimchange->dim[1] = pos_y + 1;
                  n2_aux = ap_abstract0_remove_dimensions (pr->man_dcons, true, n2_aux, dimchange);
                  ap_dimchange_free (dimchange);

                  /* n2_aux size is r->datadim + 2*r->segmentdim */

                  /* econs_aux data constraint deduced from patterns  */
                  if (econs_aux == NULL) econs_aux = ap_abstract0_copy (pr->man_dcons, n2_aux);
                  else econs_aux = ap_abstract0_meet (pr->man_dcons, true, econs_aux, n2_aux);

                  ap_abstract0_free (pr->man_dcons, n2_aux);
                }//end if pattern_sat
              else
                {
#ifndef	NDEBUG
                  fprintf (stdout, "\n@@@@ ucons_split: Warning!!! pattern not sat in split. \n If type different from 6 it's not suppose to be here!! \n");
                  fflush (stdout);
#endif
                }
            }//end aux pattern find in hash table
        }//end if pattern not null
    }//end for every pattern

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_split: econs_aux=(\n");
  if (econs_aux)
    ap_abstract0_fprint (stdout, pr->man_dcons, econs_aux, NULL);
  fprintf (stdout, ")\n");
#endif

  /* 2. transforming the universal constraints */
  /* case one: searching for the same pattern where n1 replaces with n2
   * check if condition */

  /* TODO the cases when multiple segments */
#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_split: patterns of a=(\n");
  ucons_fprint (stdout, pr->man, a, NULL);
  fprintf (stdout, "\n) resulting r=(\n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, "\n)\t n1_pat_set.size is %zu \n", n1_pat_set.size);
  fflush (stdout);
#endif

  n1_pat_set = a->n2p[n1];
  for (j = 0; j < n1_pat_set.size; j++)
    {
      /* for every pattern on the segment n1 do :
       * 2. transfer the universal constraint on n2
       * 3. set to bottom the universal formula on n1
       * 4. l'[n2] = l[n1] - 1 /\  l'[n1] = 1 /\ define dn2 in r->econs and r->ucons
       */

      pattern_j = n1_pat_set.p[j];

#ifndef NDEBUG1
      fprintf (stdout, "\n\t n1 = %zu j= %zu \n", n1, j);
      fprintf (stdout, "\t w.r.t. pattern: \n");
      //pattern_key_fprint (stdout, pr, pattern_j, NULL);
      //fprintf(stdout,"\n type %d segment %d \n ", pattern_j->type, pattern_j->segments[0]);
      if (pattern_j->type == 1) fprintf (stdout, "\n type %d segment %zu %zu \n ", pattern_j->type,
                                         pattern_j->segments[0], pattern_j->segments[1]);
      fflush (stdout);
      //		pattern_key_fprint(stdout, pr, n1_pat_set.p[j], NULL);
      //		fflush(stdout);
      //		pattern_key_fprint(stdout, pr, pattern_j, NULL);
      //		fprintf(stdout," \n");
      //		fflush(stdout);
#endif


      //#if defined (UCONS_DCONS_OCT_P21) || defined (UCONS_DCONS_POLY_P21)

      if (pattern_j->type == 1 ||
          (pr->PI[pattern_j->type].kind == pattern_2_1_lx && pattern_j->segments[1] == n1) ||
          (pr->PI[pattern_j->type].kind == pattern_2_1_mlx && pattern_j->segments[0] == n1))
        {
          /*
           *  \forall y1\in n1, y2\in n2. y1=y2
           *  \forall y2\in n2, y1\in n1. y2= y1 + l[x_1]+..+l[x_e_seg] (n2 < n1)
           *  \forall y1\in n1, y2\in n2. y1+l[x_1]+..+l[x_e_seg]= y2 (n1 < n2)
           */

          r = split_with_pattern_P21 (pr, r, pattern_j, n1, n2);

        }
      else

        if (pr->PI[pattern_j->type].kind == pattern_succ_1_2 ||
            pr->PI[pattern_j->type].kind == pattern_1_l1 ||
            pr->PI[pattern_j->type].kind == pattern_1_lx_1)
        r = split_with_pattern_succ_P12 (pr, r, pattern_j, n1, n2);
      else
        r = split_with_pattern_P12_P11 (pr, r, pattern_j, n1, n2);

    }//end checking all patterns


  /*
   * 4.define dn2 in r->econs and all udcons
   */

  if (econs_aux)
    {
      r->econs = ap_abstract0_meet (pr->man_dcons, true, econs_aux, r->econs);

      ap_dimchange_t dimadd;
      ap_abstract0_t *udcons_aux;

      size_t contorr;
      pattern_t *s;
      for (s = r->udcons; s != NULL; s = s->hh.next)
        {

          nr_y = pr->PI[s->key.type].nr_y;
          dimadd.realdim = 0;
          dimadd.intdim = 2 * nr_y;
          dimadd.dim = (ap_dim_t *) malloc ((2 * nr_y) * sizeof (ap_dim_t));
          for (contorr = 0; contorr < (2 * nr_y); contorr++)
            {
              dimadd.dim[contorr] = r->datadim + 2 * r->segmentdim;
            }

          udcons_aux =
                  ap_abstract0_add_dimensions (pr->man_dcons, false, econs_aux, &dimadd,
                                               false);
          //free(dimadd.dim);
          ap_dimchange_clear (&dimadd);
          s->dcons = ap_abstract0_meet (pr->man_dcons, true, udcons_aux, s->dcons);
          udcons_aux = NULL;
        }
    }


#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_split: before update lenghts r=(\n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, ")\n");
  fflush (stdout);
#endif

  /*
   * 5. l'[n2] = l[n1] - 1 /\  l'[n1] = 1 in r->econs and all udcons
   */

  // update the lengths on the existential constraint
  ln1 = r->datadim + r->segmentdim + n1;
  ln2 = r->datadim + r->segmentdim + n2;
  expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);

  // l[n1] == l[n2] + 1;
  ap_linexpr0_set_coeff_scalar_int (expr, ln2, 1);
  ap_linexpr0_set_cst_scalar_int (expr, 1);
  r->econs =
          ap_abstract0_substitute_linexpr (pr->man_dcons, true, r->econs, ln1, expr,
                                           NULL);

  // l[n1] := 1;
  ap_linexpr0_set_coeff_scalar_int (expr, ln2, 0);

  r->econs =
          ap_abstract0_assign_linexpr (pr->man_dcons, true, r->econs, ln1, expr,
                                       NULL);

  // meet l[n2]-1 >= 0
  ap_linexpr0_set_coeff_scalar_int (expr, ln2, 1);
  ap_linexpr0_set_cst_scalar_int (expr, -1);
  arr = ap_lincons0_array_make (1);
  arr.p[0] = ap_lincons0_make (AP_CONS_SUPEQ, expr, NULL);
  r->econs =
          ap_abstract0_meet_lincons_array (pr->man_dcons, true, r->econs, &arr);
  ap_lincons0_array_clear (&arr);

  /*
   * 4. l'[n2] = l[n1] - 1 /\  l'[n1] = 1 /\ define dn2 in r->udcons
   */

  // update the lengths on the universal constraint

  pattern_t *s;
  for (s = r->udcons; s != NULL; s = s->hh.next)
    {
      nr_y = pr->PI[s->key.type].nr_y;
      expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2 * nr_y);
      // l[n1] == l[n2] + 1;
      ap_linexpr0_set_coeff_scalar_int (expr, ln2, 1);
      ap_linexpr0_set_cst_scalar_int (expr, 1);
      s->dcons =
              ap_abstract0_substitute_linexpr (pr->man_dcons, true, s->dcons, ln1, expr,
                                               NULL);
      // l[n1] := 1;
      ap_linexpr0_set_coeff_scalar_int (expr, ln2, 0);
      s->dcons =
              ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons, ln1, expr,
                                           NULL);

      //		ap_linexpr0_set_coeff_scalar_int (expr, ln2, 1);
      //		ap_linexpr0_set_cst_scalar_int (expr, -1);
      //		ap_lincons0_array_t arr = ap_lincons0_array_make (1);
      //		arr.p[0] = ap_lincons0_make (AP_CONS_SUPEQ, expr, NULL);
      //
      //		s->dcons=
      //				ap_abstract0_meet_lincons_array (pr->man_dcons, true, s->dcons, &arr);
      //		ap_lincons0_array_clear (&arr);
    }

#ifndef NDEBUG1
  fprintf (stdout, "\n@@@@ ucons_split: returns r=(\n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, ")\n");
#endif

  return r;
}

/*
 * Used with n=1 to simulate the singleton operation.
 * Used with n=2 to simulate the split operation.
 */

ucons_t *
ucons_expand (ap_manager_t * man,
              bool destructive, ucons_t * a, ap_dim_t dim, size_t n)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_EXPAND, 0);

  if (!a)
    return NULL;
  arg_assert (n <= 2 && dim > 0 && dim < a->segmentdim, return NULL;);
  ucons_t *r;
  if (n == 1)
    r = ucons_singleton (pr, destructive, a, dim);
  else
    r = ucons_split (pr, destructive, a, dim);

  if (destructive)
    ucons_free_internal (pr, a);

#ifndef NDEBUG1
  fprintf (stdout, "\n expand returns \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  return r;
}

/**
 * Fold in @param{a} dimensions from tdim in tdim[0].
 * The restured value have dimension a->segmentsize - size + 1
 * @param man
 * @param destructive
 * @param a
 * @param tdim
 * @param size
 * @return 
 */
ucons_t *
ucons_fold (ap_manager_t * man,
            bool destructive, ucons_t * a, ap_dim_t * tdim, size_t size)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_FOLD, 0);

  if (!a)
    return NULL;
  arg_assert (tdim && 0 < size && size < a->segmentdim, return NULL;);

#ifndef NDEBUG1
  fprintf (stdout, "\t pr->segm_anon = %zu \n", pr->segm_anon);
  fflush (stdout);
#endif

#ifndef NDEBUG1
  fprintf (stdout, "\t fold: \n");
  fprintf (stdout, "\t tdim =  ");
  for (size_t o = 0; o < size; o++)
    fprintf (stdout, "\t %zu ", tdim[o]);
  fprintf (stdout, "\n on\n");
  ucons_fprint (stdout, pr->man, a, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif


  ucons_t* r = NULL;
  if (pr->active_patterns[2])
    {
#ifndef NDEBUG1
      fprintf (stdout, "\t size_fold_dim=%zu \t fold_dim: \n", size_fold_dim);
      for (size_t o = 0; o < size_fold_dim; o++)
        fprintf (stdout, "\t %zu ", fold_dim[o]);
      fflush (stdout);
#endif
      r = NULL;
      checked_realloc (fold_dim, ap_dim_t, (size_fold_dim + size),
                       sizeof (ap_dim_t), return NULL;);
      for (size_t ii = 0; ii < size; ii++)
        {
          fold_dim[ii + size_fold_dim] = tdim[ii];
        }
      size_fold_dim += size;

      if (first < pr->segm_anon)
        {
          first++;

#ifndef NDEBUG1
          fprintf (stdout, "\t fold P21 returns : \n");
          ucons_fprint (stdout, pr->man, a, NULL);
          fprintf (stdout, "\n");
          fflush (stdout);
#endif
          return a;
        }
      else
        {

          if (first == 1)
            {
              size_t up_len = pr->nr_active;

              r = fold_with_closoure_of_P21 (pr, a, fold_dim, size_fold_dim, 0,
                                             (up_len == 1));
              up_len--;

              if (pr->active_patterns[1])
                {
                  if (pr->active_patterns[4]) up_len--;
                  ucons_t* rr = fold_with_closure_P11_or_P12 (pr, r,
                                                              fold_dim, size_fold_dim, (up_len == 1));
                  ucons_free_internal (pr, r);
                  r = NULL;
                  r = rr;
                  up_len--;
                }
              if (pr->active_patterns[16])
                {
                  ucons_t* rr = fold_with_closure_succ_P12 (pr, r,
                                                            fold_dim, size_fold_dim, (up_len == 1));
                  ucons_free_internal (pr, r);
                  r = NULL;
                  r = rr;
                  up_len--;
                }
            }
          else
            {
              /* segm_anon = 2 therefore complicated fold with/without
               *  closure and (possibly) with the other patterns */
              size_t up_len = pr->nr_active;

              r = fold_without_closoure_of_P21 (pr, a, fold_dim,
                                                size_fold_dim, (up_len == 1));
              up_len--;
              if (pr->active_patterns[1])
                {
                  if (pr->active_patterns[4]) up_len--;
                  //TODO 1. decompose fold_dim and 2, call funtions on each segment
                  ucons_t* rr = fold_with_closure_P11_or_P12 (pr, r, fold_dim,
                                                              size_fold_dim, (up_len == 1));
                  ucons_free_internal (pr, r);
                  r = NULL;
                  r = rr;
                  up_len--;
                }
              if (pr->active_patterns[16])
                {
                  //TODO  1. decompose fold_dim and 2. call funtions on each segment
                  ucons_t* rr = fold_with_closure_succ_P12 (pr, r, fold_dim,
                                                            size_fold_dim, (up_len == 1));
                  ucons_free_internal (pr, r);
                  r = NULL;
                  r = rr;
                  up_len--;
                }

            }//end else first!=1

        }//end size == 1

      //if (fold_dim) free(fold_dim);
      fold_dim = NULL;
      size_fold_dim = 0;
      first = 1;
    }
  else
    {

#ifndef NDEBUG1
      fprintf (stdout, "\t \n size_fold_dim=%zu \t fold_dim:", size);
      for (size_t o = 0; o < size; o++)
        fprintf (stdout, "\t %zu ", tdim[o]);
      fflush (stdout);
#endif

      r = NULL;
      size_t up_len = pr->nr_active;
      if (pr->active_patterns[1])
        {
          if (pr->active_patterns[4]) up_len--;
          r = fold_with_closure_P11_or_P12 (pr, a, tdim, size, (up_len == 1));
          up_len--;
        }
      if (pr->active_patterns[16])
        {
          if (r == NULL)
            r = fold_with_closure_succ_P12 (pr, a, tdim, size, (up_len == 1));
          else
            {
              ucons_t* rr = fold_with_closure_succ_P12 (pr, r, tdim, size, (up_len == 1));
              ucons_free_internal (pr, r);
              r = NULL;
              r = rr;
            }
          up_len--;
        }
    }


  if (destructive)
    ucons_free_internal (pr, a);
#ifndef NDEBUG1
  fprintf (stdout, "\t fold returns : \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif

  return r;

}

/*  generates the constraints corresponding to the pattern y1<y2  */
ucons_t *
merge_succesors_2 (ucons_internal_t *pr, ucons_t * r, ap_dim_t * gdim,
                   size_t size_gdim, bool update_lenght)
{

  size_t i, j;
  ap_abstract0_t * udcons, *g_dcons;
  pattern_t * n_dcons, *n0_dcons;
  size_t contor, u_seg;
  ap_dim_t ly, dy;
  ap_dim_t ly1, ly2, dy1, dy2;

  ap_dimchange_t * dimchange;
  pattern_key_t * look;

  bool chainge;
  ly1 = r->datadim + 2 * r->segmentdim;
  ly2 = r->datadim + 2 * r->segmentdim + 1;
  dy1 = r->datadim + 2 * r->segmentdim + 2;
  dy2 = r->datadim + 2 * r->segmentdim + 3;
  /* generating data constraint with the pattern \forall y for the node gdim[0] from gdim[1..size_gdim] */

  u_seg = pr->PI[0].u_seg;
  checked_malloc (look, pattern_key_t, sizeof (pattern_key_t) + u_seg * sizeof (size_t), 1, return NULL;);
  memset (look, 0, sizeof (pattern_key_t) + u_seg * sizeof (size_t));
  look->type = 0;
  look->segments[0] = gdim[0];
  unsigned keylen = u_seg * sizeof (size_t) + sizeof (pattern_key_t);

  HASH_FIND (hh, r->udcons, look, keylen, n_dcons);

  if (!n_dcons)
    {
      checked_malloc (n_dcons, pattern_t, 1, sizeof (pattern_t) + u_seg * sizeof (size_t), return NULL;);
      memset (n_dcons, 0, sizeof (pattern_t) + u_seg * sizeof (size_t));
      n_dcons->key.type = look->type;
      for (size_t i = 0; i < (u_seg); i++)
        n_dcons->key.segments[i] = look->segments[i];
      n_dcons->dcons = NULL;

      HASH_ADD (hh, r->udcons, key, keylen, n_dcons);
      r = add_pattern_n2p (pr, r, look);
    }


  ap_abstract0_t * aux_p1_ndcons = NULL;
  /* aux var to keep the constraint with p11*/
  for (i = 1; i < size_gdim; i++)
    {
      dimchange = ap_dimchange_alloc (2, 0);
      dimchange->dim[0] = r->datadim + 2 * r->segmentdim;
      dimchange->dim[1] = r->datadim + 2 * r->segmentdim;
      ly = r->datadim + 2 * r->segmentdim;
      dy = r->datadim + 2 * r->segmentdim + 1;
      udcons = NULL;
      udcons = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, dimchange, false);

      ap_dimchange_free (dimchange);

      ap_lincons0_array_t arr = ap_lincons0_array_make (2);
      // l[y]-i == 0
      arr.p[0].constyp = AP_CONS_EQ;
      arr.p[0].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
      arr.p[0].scalar = NULL;
      ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                        r->datadim + 2 * r->segmentdim, 1);
      ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -i);
      // d(y) - d(i) ==0
      arr.p[1].constyp = AP_CONS_EQ;
      arr.p[1].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
      arr.p[1].scalar = NULL;
      ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                        r->datadim + gdim[i], 1);
      ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                        r->datadim + 2 * r->segmentdim + 1, -1);
      ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);

      udcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, udcons, &arr);

      ap_lincons0_array_clear (&arr);

      if (aux_p1_ndcons == NULL)
        {
          aux_p1_ndcons = ap_abstract0_copy (pr->man_dcons, udcons);
          ap_abstract0_free (pr->man_dcons, udcons);
        }
      else
        {
          aux_p1_ndcons = ap_abstract0_join (pr->man_dcons, true, aux_p1_ndcons, udcons);
        }
      udcons = NULL;
    }
  //}
  /*end generation for the first pattern */
  /* TODO function for each pattern type */

  /* generating data constraint with the pattern \forall y1,y2. y1<y2 for the node gdim[0] from gdim[1..size_gdim] */
  /**********************************************************/

#ifndef NDEBUG1
  fprintf (stdout, "\n merge_1");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif


  pattern_t * n2_dcons;
  ap_abstract0_t *udcons2;
  size_t u_seg2;
  pattern_key_t *look2;
  unsigned keylen2;

  u_seg2 = pr->PI[3].u_seg;
  checked_malloc (look2, pattern_key_t, sizeof (pattern_key_t) + u_seg2 * sizeof (size_t), 1, return NULL;);
  memset (look2, 0, sizeof (pattern_key_t) + u_seg2 * sizeof (size_t));
  look2->type = 3;
  look2->segments[0] = gdim[0];
  keylen2 = u_seg2 * sizeof (size_t) + sizeof (pattern_key_t);

  HASH_FIND (hh, r->udcons, look2, keylen2, n2_dcons);

  if (!n2_dcons)
    {
      checked_malloc (n2_dcons, pattern_t, 1, sizeof (pattern_t) + u_seg2 * sizeof (size_t), return NULL;);
      memset (n2_dcons, 0, sizeof (pattern_t) + u_seg2 * sizeof (size_t));
      n2_dcons->key.type = look2->type;
      for (size_t i = 0; i < (u_seg2); i++)
        n2_dcons->key.segments[i] = look2->segments[i];
      n2_dcons->dcons = NULL;
      HASH_ADD (hh, r->udcons, key, keylen2, n2_dcons);
      r = add_pattern_n2p (pr, r, look2);
    }

  HASH_FIND (hh, r->udcons, look2, keylen2, n2_dcons);

  if (n2_dcons)
    {
      n2_dcons->dcons = NULL;
      //chainge = false;
      for (i = 1; i < size_gdim; i++)
        {
          for (j = i; j < size_gdim; j++)
            {
              dimchange = ap_dimchange_alloc (4, 0);
              dimchange->dim[0] = r->datadim + 2 * r->segmentdim;
              dimchange->dim[1] = r->datadim + 2 * r->segmentdim;
              dimchange->dim[2] = r->datadim + 2 * r->segmentdim;
              dimchange->dim[3] = r->datadim + 2 * r->segmentdim;
              ly1 = r->datadim + 2 * r->segmentdim;
              ly2 = r->datadim + 2 * r->segmentdim + 1;
              dy1 = r->datadim + 2 * r->segmentdim + 2;
              dy2 = r->datadim + 2 * r->segmentdim + 3;

              udcons2 = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, dimchange, false);
#ifndef NDEBUG1
              printf ("\n merge_2 udcons ");
              if (udcons2) ap_abstract0_fprint (stdout, pr->man_dcons, udcons2, NULL);
              printf ("\n");
#endif
              ap_dimchange_free (dimchange);

              ap_lincons0_array_t arr = ap_lincons0_array_make (4);

              // l[y1]-i == 0
              arr.p[0].constyp = AP_CONS_EQ;
              arr.p[0].linexpr0 =
                      ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
              arr.p[0].scalar = NULL;
              ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ly1, 1);
              ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -i);
              //l[y2]-j=0
              arr.p[1].constyp = AP_CONS_EQ;
              arr.p[1].linexpr0 =
                      ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
              arr.p[1].scalar = NULL;
              ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, ly2, 1);
              ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, -j);
              // d(y1) - d(i) ==0
              arr.p[2].constyp = AP_CONS_EQ;
              arr.p[2].linexpr0 =
                      ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
              arr.p[2].scalar = NULL;
              ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0,
                                                r->datadim + gdim[i], 1);
              ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0, dy1, -1);
              ap_linexpr0_set_cst_scalar_int (arr.p[2].linexpr0, 0);
              // d(y2) - d(j) ==0
              arr.p[3].constyp = AP_CONS_EQ;
              arr.p[3].linexpr0 =
                      ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
              arr.p[3].scalar = NULL;
              ap_linexpr0_set_coeff_scalar_int (arr.p[3].linexpr0,
                                                r->datadim + gdim[j], 1);
              ap_linexpr0_set_coeff_scalar_int (arr.p[3].linexpr0, dy2, -1);
              ap_linexpr0_set_cst_scalar_int (arr.p[3].linexpr0, 0);

              udcons2 = ap_abstract0_meet_lincons_array (pr->man_dcons, true, udcons2, &arr);

#ifndef NDEBUG1
              printf ("\n merge_2 udcons ");
              if (udcons2) ap_abstract0_fprint (stdout, pr->man_dcons, udcons2, NULL);
              printf ("\n");
#endif

              ap_lincons0_array_clear (&arr);

              if (n2_dcons->dcons == NULL)
                {
                  n2_dcons->dcons = ap_abstract0_copy (pr->man_dcons, udcons2);
                  ap_abstract0_free (pr->man_dcons, udcons2);
                }
              else n2_dcons->dcons = ap_abstract0_join (pr->man_dcons, true, n2_dcons->dcons, udcons2);

#ifndef NDEBUG1
              printf ("\n merge_2");
              ucons_fprint (stdout, pr->man, r, NULL);
              printf ("\n");
#endif

            }//end for
        }//end for

      if (n2_dcons->dcons == NULL)
        n2_dcons->dcons = ap_abstract0_bottom (pr->man_dcons, r->datadim + 2 * r->segmentdim + 4, 0);

    }/* end generation with pattern pattern_2_1 y1<y2 for gdim[0] */

  HASH_FIND (hh, r->udcons, look, keylen, n_dcons);
  if (n_dcons)
    {
      n_dcons->dcons = aux_p1_ndcons;
    }
#ifndef NDEBUG1
  printf ("\n merge_2 before update lengths");
  ucons_fprint (stdout, pr->man, r, NULL);
  printf ("\n");
#endif

  if (update_lenght)
    {
      /* l(gdim[0]) == 1 ; */
      ap_linexpr0_t *expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
      ap_linexpr0_set_cst_scalar_int (expr, 1);

      ap_linexpr0_t *expr_y = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
      //ap_linexpr0_set_cst_scalar_int (expr_y, 1);

      ap_linexpr0_t *expr_2y = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
      //ap_linexpr0_set_cst_scalar_int (expr_2y, 1);

      ap_dim_t n, ln;
      n = gdim[0];
      ln = r->datadim + r->segmentdim + n;
      r->econs = ap_abstract0_substitute_linexpr (pr->man_dcons, true, r->econs, ln, expr, NULL);
      // l(gdim[0]) := size_gdim
      ap_linexpr0_set_cst_scalar_int (expr, size_gdim);
      r->econs = ap_abstract0_assign_linexpr (pr->man_dcons, true, r->econs, ln, expr, NULL);

      ap_linexpr0_free (expr);

#ifndef NDEBUG1
      ap_lincons0_t cons1 = cons_k_m (r->datadim, r->segmentdim, 1, 2);
      if (ap_abstract0_sat_lincons (pr->man_dcons, r->econs, &cons1))
        printf ("\n \t k<m \n");
      else
        printf ("\n \t !!!(k<m) \n");
#endif

      pattern_t *s;
      size_t nr_y;
      for (s = r->udcons; s != NULL; s = s->hh.next)
        {
          nr_y = pr->PI[s->key.type].nr_y;
          if (nr_y == 1)
            {
              ap_linexpr0_set_cst_scalar_int (expr_y, 1);
              s->dcons = ap_abstract0_substitute_linexpr (pr->man_dcons, true, s->dcons, ln, expr_y, NULL);
              // l(gdim[0]) := size_gdim
              ap_linexpr0_set_cst_scalar_int (expr_y, size_gdim);
              s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons, ln, expr_y, NULL);

#ifndef NDEBUG1
              n = s->key.segments[0];
              ap_lincons0_t cons2 = cons_dy0 (r->datadim, r->segmentdim, r->datadim + 2 * r->segmentdim + 1, n + r->datadim);
              if (ap_abstract0_sat_lincons (pr->man_dcons, s->dcons, &cons2))
                printf ("\n \t d(n)<=dy0 \n");
              else
                printf ("\n \t !!!(d(n)<=dy0) \n");
#endif

            }
          else
            {
              ap_linexpr0_set_cst_scalar_int (expr_2y, 1);
              s->dcons = ap_abstract0_substitute_linexpr (pr->man_dcons, true, s->dcons, ln, expr_2y, NULL);
              // l(gdim[0]) := size_gdim
              ap_linexpr0_set_cst_scalar_int (expr_2y, size_gdim);
              s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons, ln, expr_2y, NULL);

#ifndef NDEBUG1
              ap_lincons0_t cons = cons_dy0_dy1 (r->datadim, r->segmentdim, dy1, dy2);
              if (ap_abstract0_sat_lincons (pr->man_dcons, s->dcons, &cons))
                printf ("\n \t dy0<dy1 \n");
              else
                printf ("\n \t !!!(dy0<dy1) \n");
#endif
            }
        }

      ap_linexpr0_free (expr_y);
      ap_linexpr0_free (expr_2y);
    }
#ifndef NDEBUG1
  printf ("\n merge returns");
  ucons_fprint (stdout, pr->man, r, NULL);
  printf ("\n");
#endif
  return r;
}

/* for the first pattern : \forall y generate the constraints */
ucons_t *
merge_succesors_1 (ucons_internal_t *pr, ucons_t * r, ap_dim_t * gdim,
                   size_t size_gdim, bool update_lenght)
{


  arg_assert (r, return NULL;);
  arg_assert (gdim && size_gdim > 1, return r;);

#ifndef NDEBUG1
  printf ("  \n\t merge_succesors_1 : \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  printf ("\n");
#endif

  size_t i, j;
  ap_abstract0_t * udcons, *g_dcons;
  pattern_t * n_dcons;
  size_t contor, u_seg;
  ap_dim_t ly, dy;
  ap_dim_t ly1, ly2, dy1, dy2;

  ap_dimchange_t * dimchange;
  pattern_key_t * look;

  /* generating data constraint witn the pattern \forall y for the node gdim[0] from gdim[1..size_gdim] */

  u_seg = pr->PI[0].u_seg;
  checked_malloc (look, pattern_key_t, sizeof (pattern_key_t) + u_seg * sizeof (size_t), 1, return NULL;);
  memset (look, 0, sizeof (pattern_key_t) + u_seg * sizeof (size_t));
  look->type = 0;
  look->segments[0] = gdim[0];
  unsigned keylen = u_seg * sizeof (size_t) + sizeof (pattern_key_t);

  HASH_FIND (hh, r->udcons, look, keylen, n_dcons);

  if (!n_dcons)
    {
      checked_malloc (n_dcons, pattern_t, 1, sizeof (pattern_t) + u_seg * sizeof (size_t), return NULL;);
      memset (n_dcons, 0, sizeof (pattern_t) + u_seg * sizeof (size_t));
      n_dcons->key.type = look->type;
      for (size_t i = 0; i < (u_seg); i++)
        n_dcons->key.segments[i] = look->segments[i];
      n_dcons->dcons = NULL;
      HASH_ADD (hh, r->udcons, key, keylen, n_dcons);
      r = add_pattern_n2p (pr, r, look);
    }

  HASH_FIND (hh, r->udcons, look, keylen, n_dcons);

  if (n_dcons)
    {

      n_dcons->dcons = NULL;
      for (i = 1; i < size_gdim; i++)
        {
          dimchange = ap_dimchange_alloc (2, 0);
          dimchange->dim[0] = r->datadim + 2 * r->segmentdim;
          dimchange->dim[1] = r->datadim + 2 * r->segmentdim;
          ly = r->datadim + 2 * r->segmentdim;
          dy = r->datadim + 2 * r->segmentdim + 1;

          udcons = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, dimchange, false);

          ap_dimchange_free (dimchange);

          ap_lincons0_array_t arr = ap_lincons0_array_make (2);
          // l[y]-i == 0
          arr.p[0].constyp = AP_CONS_EQ;
          arr.p[0].linexpr0 =
                  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
          arr.p[0].scalar = NULL;
          ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                            r->datadim + 2 * r->segmentdim, 1);
          ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, -i);
          // d(y) - d(i) ==0
          arr.p[1].constyp = AP_CONS_EQ;
          arr.p[1].linexpr0 =
                  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
          arr.p[1].scalar = NULL;
          ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                            r->datadim + gdim[i], 1);
          ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                            r->datadim + 2 * r->segmentdim + 1, -1);
          ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);

          udcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, udcons, &arr);

          ap_lincons0_array_clear (&arr);

          if (n_dcons->dcons == NULL)
            n_dcons->dcons = ap_abstract0_copy (pr->man_dcons, udcons);
          else n_dcons->dcons = ap_abstract0_join (pr->man_dcons, false, n_dcons->dcons, udcons);


          ap_abstract0_free (pr->man_dcons, udcons);
        }
    }
  /*end generation for the first pattern */
  /* TODO function for each pattern type */

  if (update_lenght)
    {
      /* l(gdim[0]) == 1 ; */
      ap_linexpr0_t *expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
      ap_linexpr0_set_cst_scalar_int (expr, 1);

      ap_linexpr0_t *expr_y = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
      ap_linexpr0_set_cst_scalar_int (expr_y, 1);

      ap_dim_t n, ln;
      n = gdim[0];
      ln = r->datadim + r->segmentdim + n;
      r->econs = ap_abstract0_substitute_linexpr (pr->man_dcons, true, r->econs, ln, expr, NULL);
      // l(gdim[0]) := size_gdim
      ap_linexpr0_set_cst_scalar_int (expr, size_gdim);
      r->econs = ap_abstract0_assign_linexpr (pr->man_dcons, true, r->econs, ln, expr, NULL);

      ap_linexpr0_free (expr);


      pattern_t *s;
      size_t nr_y;
      for (s = r->udcons; s != NULL; s = s->hh.next)
        {
          s->dcons = ap_abstract0_substitute_linexpr (pr->man_dcons, true, s->dcons, ln, expr_y, NULL);
          // l(gdim[0]) := size_gdim
          ap_linexpr0_set_cst_scalar_int (expr_y, size_gdim);
          s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons, ln, expr_y, NULL);
          ap_linexpr0_set_cst_scalar_int (expr_y, 1);
        }

      ap_linexpr0_free (expr_y);
    }

#ifndef NDEBUG1
  printf ("  \n\t merge_succesor_P1 returns: \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  printf ("\n");
#endif

  free (look);
  return r;
}

ap_abstract0_t *
concat_P11 (ucons_internal_t * pr, ucons_t * r, ap_dim_t * cdim,
            size_t size_cdim)
{

  ap_abstract0_t * aux_dcons = NULL;
  ap_abstract0_t * fy_dcons;
  size_t u_seg;
  size_t i, j, k;
  ap_abstract0_t * ni_dcons;
  ap_dim_t li, ly, di, dy;

  pattern_t * n_dcons;

  ap_dimchange_t dimchange;
  pattern_key_t * look;

  arg_assert (r, return NULL;);
  arg_assert (cdim && size_cdim > 1, return NULL;);

  u_seg = pr->PI[0].u_seg;
  checked_malloc (look, pattern_key_t, sizeof (pattern_key_t) + u_seg * sizeof (size_t), 1, return NULL;);
  memset (look, 0, sizeof (pattern_key_t) + u_seg * sizeof (size_t));
  look->type = 0;
  unsigned keylen = u_seg * sizeof (size_t) + sizeof (pattern_key_t);

#ifndef NDEBUG1
  printf ("  \n\t concat_P11 : \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  printf ("\n");
#endif

  //aux_dcons=ap_abstract0_copy(pr->man_dcons,n_dcons->dcons);
  for (i = 1; i < size_cdim; i++)
    {
      look->segments[0] = cdim[i];
      HASH_FIND (hh, r->udcons, look, keylen, n_dcons);


      if (!n_dcons)
        {
          if (test_singleton (pr->man_dcons, r->econs, r->datadim, r->segmentdim, cdim[i]) == false)
            {
#ifndef NDEBUG2
              printf ("  \n\t concat_P11 return with pattern info insufficient: \n");
              ucons_fprint (stdout, pr->man, r, NULL);
              printf ("\n");
#endif
              return NULL;
            }
        }

      if (n_dcons)
        {
          /* TODO add element corresponding to the existential */
          /* join with the universal from cdim[i] not needed for this pattern */
          if (n_dcons->dcons)
            {

              /* */
              ap_linexpr0_t *expr_y = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
              ap_linexpr0_set_cst_scalar_int (expr_y, 0); /*TODO check the definition of the sustitution */
              for (j = 0; j < i; j++)
                {
                  li = r->datadim + r->segmentdim + cdim[j];
                  ap_linexpr0_set_coeff_scalar_int (expr_y, li, -1); /*TODO check the definition of the sustitution */
                }
              ly = r->datadim + 2 * r->segmentdim;
              ap_linexpr0_set_coeff_scalar_int (expr_y, ly, 1);
              ap_abstract0_t *aux_ndcons = ap_abstract0_substitute_linexpr (pr->man_dcons,
                                                                            false, n_dcons->dcons, ly, expr_y, NULL);
              if (aux_dcons)
                aux_dcons = ap_abstract0_join (pr->man_dcons, true, aux_dcons, aux_ndcons);
              else
                {
                  aux_dcons = ap_abstract0_copy (pr->man_dcons, aux_ndcons);
                  ap_abstract0_free (pr->man_dcons, aux_ndcons);
                }

              ap_linexpr0_free (expr_y);
            }
          else
            {
#ifndef NDEBUG2
              fprintf (stdout, "  \n\t concat_P11: \n");
              fprintf (stdout, "P(n%zu)==>null ", i);
              fprintf (stdout, "\n");
              fflush (stdout);
#endif
            }
        }
#ifndef NDEBUG1
      fprintf (stdout, "  \n\tconcat_P11 aux_dcons: \n");
      if (aux_dcons) ap_abstract0_fprint (stdout, pr->man_dcons, aux_dcons, NULL);
      fprintf (stdout, "\n");
      fflush (stdout);
#endif

      /* ni_dcons  =  econs + y == l[cdim[0]] + ... l[cdim[i]]  + d(y) = d(cdim[i]) */

      dimchange.intdim = 2;
      dimchange.realdim = 0;
      dimchange.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
      dimchange.dim[0] = r->datadim + 2 * r->segmentdim;
      dimchange.dim[1] = r->datadim + 2 * r->segmentdim;

      ni_dcons = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimchange, false);
      free (dimchange.dim);


      ap_lincons0_array_t arr = ap_lincons0_array_make (2);
      // l[y] == l[cdim[0]] + ... l[cdim[i-1]]
      arr.p[0].constyp = AP_CONS_EQ;
      arr.p[0].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
      arr.p[0].scalar = NULL;

      for (j = 0; j < i; j++)
        {
          li = r->datadim + r->segmentdim + cdim[j];
          ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, li, -1);
        }

      ly = r->datadim + 2 * r->segmentdim;
      ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ly, 1);
      ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

      // d(y) - d(cdim[i]) ==0
      arr.p[1].constyp = AP_CONS_EQ;
      arr.p[1].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
      arr.p[1].scalar = NULL;
      ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                        r->datadim + cdim[i], -1);
      ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                        r->datadim + 2 * r->segmentdim + 1, 1);

      ni_dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, ni_dcons, &arr);

      ap_lincons0_array_clear (&arr);

      if (aux_dcons)
        aux_dcons = ap_abstract0_join (pr->man_dcons, true, aux_dcons, ni_dcons);
      else
        {
          aux_dcons = ap_abstract0_copy (pr->man_dcons, ni_dcons);
          ap_abstract0_free (pr->man_dcons, ni_dcons);
        }


    }

#ifndef NDEBUG1
  fprintf (stdout, "\n CONCAT_P11 retuns: \t ");
  if (aux_dcons) ap_abstract0_fprint (stdout, pr->man_dcons, aux_dcons, NULL);
  fprintf (stdout, "\n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fflush (stdout);
#endif
  return aux_dcons;
}

ap_abstract0_t *
concat_P12 (ucons_internal_t * pr, ucons_t * r, ap_dim_t * cdim,
            size_t size_cdim)
{

  size_t i, j, k, wl;
  ap_abstract0_t * aux_2y_dcons = NULL; // the abstract attached to cdim[0] with P12
  ap_abstract0_t * fy2_dcons;

  pattern_t *found2_dcons, *found_right_dcons;
  size_t contor, u_seg, u_seg2;
  ap_dim_t ly1, dy1, ly2, dy2;

  ap_abstract0_t *dy1_econs;
  ap_dim_t li, ly, di, dy;

  ap_dimchange_t dimchange;
  pattern_key_t * look, *look2;

  arg_assert (r, return NULL;);
  arg_assert (cdim && size_cdim > 1, return NULL;);

  u_seg = pr->PI[0].u_seg;
  checked_malloc (look, pattern_key_t, sizeof (pattern_key_t) + u_seg * sizeof (size_t), 1, return NULL;);
  memset (look, 0, sizeof (pattern_key_t) + u_seg * sizeof (size_t));
  look->type = 0;
  unsigned keylen = u_seg * sizeof (size_t) + sizeof (pattern_key_t);

  /* concat segments w.r.t. the pattern \forall y1,y2. y1<y2 for the node gdim[0] from gdim[1..size_gdim] */
  /* aux_2y_dcons is the resulting data constraint */
  u_seg2 = pr->PI[3].u_seg;
  checked_malloc (look2, pattern_key_t, sizeof (pattern_key_t) + u_seg2 * sizeof (size_t), 1, return NULL;);
  memset (look2, 0, sizeof (pattern_key_t) + u_seg2 * sizeof (size_t));
  look2->type = 3;
  unsigned keylen2 = u_seg2 * sizeof (size_t) + sizeof (pattern_key_t);

#ifndef NDEBUG1
  fprintf (stdout, "  \n\t concat_P12 : \n");
  for (size_t dd = 0; dd < size_cdim; dd++)
    fprintf (stdout, " %zu ", cdim[dd]);
  fprintf (stdout, "on: \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  ap_abstract0_t* border, *p11_datacons = NULL;
  border = NULL;
  pattern_t *p11_dcons = NULL;
  ;
  ly1 = r->datadim + 2 * r->segmentdim;
  ly2 = r->datadim + 2 * r->segmentdim + 1;
  dy1 = r->datadim + 2 * r->segmentdim + 2;
  dy2 = r->datadim + 2 * r->segmentdim + 3;

  for (i = 1; i < size_cdim; i++)
    {

      /* 1. join universal P12 formulas */
      look2->segments[0] = cdim[i];
      HASH_FIND (hh, r->udcons, look2, keylen2, found2_dcons);

      look->segments[0] = cdim[i];
      HASH_FIND (hh, r->udcons, look, keylen, p11_dcons);
      if (!p11_dcons)
        {
          if (!test_singleton (pr->man_dcons, r->econs, r->datadim, r->segmentdim, cdim[i]))
            {
#ifndef NDEBUG2
              fprintf (stdout, "\t insufficient information \n ");
              fprintf (stdout, "forall y in cdim[%zu] not defined \n ", cdim[i]);
              fprintf (stdout, "fold_with_P12 returns: \n");
              ucons_fprint (stdout, pr->man, r, NULL);
              fflush (stdout);
#endif
              return NULL;
            }
        }
      else
        {

          ap_dimchange_t dimadd;

          ap_dimchange_init (&dimadd, 2, 0);
          dimadd.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
          dimadd.dim[0] = r->datadim + 2 * r->segmentdim + 0;
          dimadd.dim[1] = r->datadim + 2 * r->segmentdim + 1;


          p11_datacons = ap_abstract0_copy (pr->man_dcons, p11_dcons->dcons);
#ifndef NDEBUG1
          fprintf (stdout, "\n forall y in cdim[%zu] p11 without without extra dims ", i);
          if (p11_datacons) ap_abstract0_fprint (stdout, pr->man_dcons, p11_datacons, NULL);
          fprintf (stdout, "\n");
#endif
          p11_datacons = ap_abstract0_add_dimensions (pr->man_dcons, true, p11_datacons, &dimadd, false);
#ifndef NDEBUG1
          fprintf (stdout, "\n forall y in cdim[%zu] p11 with without extra dims ", i);
          if (p11_datacons) ap_abstract0_fprint (stdout, pr->man_dcons, p11_datacons, NULL);
          fprintf (stdout, "\n");
#endif
          ap_dimchange_clear (&dimadd);
          /*
           * added y -- > y2 add y1
           */
          //			ap_linexpr0_t *expr_2y =  ap_linexpr0_alloc (AP_LINEXPR_DENSE,
          //					r->datadim + 2 * r->segmentdim + 4);
          //			ap_linexpr0_set_cst_scalar_int (expr_2y, 0);//!!!
          //			for(j = 0; j < i; j++){
          //				li = r->datadim + r->segmentdim + cdim[j] ;
          //				ap_linexpr0_set_coeff_scalar_int (expr_2y, li, -1);//!!!
          //			}
          //			ap_linexpr0_set_coeff_scalar_int (expr_2y, ly2, 1);
          //
          //			p11_datacons = ap_abstract0_substitute_linexpr (pr->man_dcons, true, p11_datacons,ly2 , expr_2y, NULL);
          //
          //			ap_linexpr0_free(expr_2y);


#ifndef NDEBUG1
          fprintf (stdout, "\n forall y in cdim[%zu] p11 with extra dims ", i);
          if (p11_datacons) ap_abstract0_fprint (stdout, pr->man_dcons, p11_datacons, NULL);
          fprintf (stdout, "\n");
#endif
        }

      if (!found2_dcons)
        {
          if (!test_singleton (pr->man_dcons, r->econs, r->datadim, r->segmentdim, cdim[i]))
            {
#ifndef NDEBUG2
              fprintf (stdout, "\t insufficient information \n ");
              fprintf (stdout, "forall y1<y2 in cdim[%zu] not defined \n ", cdim[i]);
              fprintf (stdout, "fold_with_P12 returns: \n");
              ucons_fprint (stdout, pr->man, r, NULL);
              fflush (stdout);
#endif
              return NULL;
            }
        }
      else
        {
          if (found2_dcons->dcons)
            {
              /* y1 = y1 + l(cdim[0]) + ... +l(cdim[i-1]) + 1 */
              /* y2 = y2 + l(cdim[0]) + ... +l(cdim[i-1]) + 1*/
              ap_linexpr0_t *expr_2y = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                          r->datadim + 2 * r->segmentdim + 4);
              ap_linexpr0_set_cst_scalar_int (expr_2y, 0); //!!!
              for (j = 0; j < i; j++)
                {
                  li = r->datadim + r->segmentdim + cdim[j];
                  ap_linexpr0_set_coeff_scalar_int (expr_2y, li, -1); //!!!
                }

              ap_linexpr0_set_coeff_scalar_int (expr_2y, ly1, 1);

              fy2_dcons = ap_abstract0_copy (pr->man_dcons, found2_dcons->dcons);
              fy2_dcons = ap_abstract0_substitute_linexpr (pr->man_dcons, false, fy2_dcons, ly1, expr_2y, NULL);

              ap_linexpr0_set_coeff_scalar_int (expr_2y, ly1, 0);
              ap_linexpr0_set_coeff_scalar_int (expr_2y, ly2, 1);

              fy2_dcons = ap_abstract0_substitute_linexpr (pr->man_dcons, true, fy2_dcons, ly2, expr_2y, NULL);

              ap_linexpr0_free (expr_2y);

              if (aux_2y_dcons)
                aux_2y_dcons = ap_abstract0_join (pr->man_dcons, true, aux_2y_dcons, fy2_dcons);
              else
                {
                  aux_2y_dcons = ap_abstract0_copy (pr->man_dcons, fy2_dcons);
                  ap_abstract0_free (pr->man_dcons, fy2_dcons);
                }
            }
          else
            {
              if (!test_singleton (pr->man_dcons, r->econs, r->datadim, r->segmentdim, cdim[i]))
                {
#ifndef NDEBUG2
                  fprintf (stdout, "\t insufficient information \n ");
                  fprintf (stdout, "forall y1<y2 in cdim[%zu] is null \n ", cdim[i]);
                  fprintf (stdout, "fold_with_P12 returns: \n");
                  ucons_fprint (stdout, pr->man, r, NULL);
                  fflush (stdout);
#endif
                  return NULL;
                }
            }
        }


#ifndef NDEBUG1
      fprintf (stdout, "\n concat 1 val currenta cu patternul P12. %zu", i);
      if (aux_2y_dcons) ap_abstract0_fprint (stdout, pr->man_dcons, aux_2y_dcons, NULL);
      fprintf (stdout, "\n");
      fflush (stdout);
#endif



      /*  join with y1 = cdim[i] y2 = \forall y \in cdim[i]*/

      ap_abstract0_t *vi = NULL;
      if (p11_datacons)
        {



          ap_abstract0_t *econs_i = ap_abstract0_copy (pr->man_dcons, r->econs);

          ap_dimchange_t dimadd2;

          ap_dimchange_init (&dimadd2, 4, 0);
          dimadd2.intdim = 4;
          dimadd2.realdim = 0;
          dimadd2.dim = (ap_dim_t*) malloc (4 * sizeof (ap_dim_t));
          dimadd2.dim[0] = r->datadim + 2 * r->segmentdim;
          dimadd2.dim[1] = r->datadim + 2 * r->segmentdim;
          dimadd2.dim[2] = r->datadim + 2 * r->segmentdim;
          dimadd2.dim[3] = r->datadim + 2 * r->segmentdim;


          econs_i = ap_abstract0_add_dimensions (pr->man_dcons, true, econs_i, &dimadd2, false);

          ap_dimchange_clear (&dimadd2);

          ap_lincons0_array_t arr = ap_lincons0_array_make (2);

          // l[y1] - l(cdim[0]) + ... + l(cdim[i-1]) ==0
          arr.p[0].constyp = AP_CONS_EQ;
          arr.p[0].linexpr0 =
                  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
          arr.p[0].scalar = NULL;
          for (j = 0; j < i; j++)
            {
              li = r->datadim + r->segmentdim + cdim[j];
              ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, li, 1);
            }
          ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ly1, -1);
          //ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

          // d(y1) - d(cdim[i]) ==0
          arr.p[1].constyp = AP_CONS_EQ;
          arr.p[1].linexpr0 =
                  ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
          arr.p[1].scalar = NULL;
          ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                            r->datadim + cdim[i], 1);
          ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy1, -1);
          //ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);


          econs_i = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                                     true, econs_i, &arr);
          ap_lincons0_array_clear (&arr);
#ifndef NDEBUG1
          fprintf (stdout, "\n cdim[%zu] and forall y in cdim[%zu] ", i, i);
          if (econs_i) ap_abstract0_fprint (stdout, pr->man_dcons, econs_i, NULL);
          if (p11_datacons) ap_abstract0_fprint (stdout, pr->man_dcons, p11_datacons, NULL);
          fprintf (stdout, "\n");
          fflush (stdout);
#endif

          ap_abstract0_t *aux_i_ydconsi = ap_abstract0_meet (pr->man_dcons,
                                                             true, p11_datacons, econs_i);
#ifndef NDEBUG1
          fprintf (stdout, "\n cdim[%zu] and forall y in cdim[%zu] result: ", i, i);
          if (aux_i_ydconsi) ap_abstract0_fprint (stdout, pr->man_dcons, aux_i_ydconsi, NULL);
          fprintf (stdout, "\n");
          fflush (stdout);
#endif

          if (aux_2y_dcons)
            aux_2y_dcons = ap_abstract0_join (pr->man_dcons, true, aux_2y_dcons, aux_i_ydconsi);
          else
            {
              aux_2y_dcons = ap_abstract0_copy (pr->man_dcons, aux_i_ydconsi);
              ap_abstract0_free (pr->man_dcons, aux_i_ydconsi);
            }
        }
      /* 2. join all  (y,cdim[wl]) /\ y\in cdim[i] for all wl in 0..i-1*/
      pattern_t * found_left_dcons = NULL;
      ap_abstract0_t* dy2_econs = NULL;



      for (wl = 0; wl < i; wl++)
        {

          look->segments[0] = cdim[wl];
          HASH_FIND (hh, r->udcons, look, keylen, found_left_dcons);
          if (!found_left_dcons)
            {
              if (!test_singleton (pr->man_dcons, r->econs, r->datadim, r->segmentdim, cdim[wl]))
                {
#ifndef NDEBUG2
                  fprintf (stdout, "\t insufficient information \n ");
                  fprintf (stdout, "forall y in cdim[%zu] is null \n ", cdim[wl]);
                  fprintf (stdout, "fold_without_closure returns: \n");
                  ucons_fprint (stdout, pr->man, r, NULL);
                  fflush (stdout);
#endif
                  return NULL;
                }
            }
          if (found_left_dcons)
            {
              if (found_left_dcons->dcons)
                {

                  /* border = add dimensions to  \forall y in cdim[i]*/
                  dimchange.intdim = 2;
                  dimchange.realdim = 0;
                  dimchange.dim = (ap_dim_t*) malloc (2 * sizeof (ap_dim_t));
                  dimchange.dim[0] = r->datadim + 2 * r->segmentdim + 1;
                  dimchange.dim[1] = r->datadim + 2 * r->segmentdim + 2;

                  border = ap_abstract0_copy (pr->man_dcons, found_left_dcons->dcons);
                  border = ap_abstract0_add_dimensions (pr->man_dcons, true, border, &dimchange, false);
                  free (dimchange.dim);

#ifndef NDEBUG2
                  printf ("\n border : %zu", i);
                  if (border) ap_abstract0_fprint (stdout, pr->man_dcons, border, NULL);
                  printf ("\n");
#endif


                  /**************************
                   * border /\p11_datacons < = >
                   * forall y1 in cdim[wl] /\ forall y2 in cdim[i]
                   */
                  ap_abstract0_t *aux_p11_datacons = NULL;
                  if (p11_datacons)
                    {

                      aux_p11_datacons = ap_abstract0_meet (pr->man_dcons, false, p11_datacons, border);
                      ap_abstract0_free (pr->man_dcons, p11_datacons);
                      p11_datacons = aux_p11_datacons;
#ifndef NDEBUG2
                      printf ("\n border left P11 : %zu", i);
                      if (border) ap_abstract0_fprint (stdout, pr->man_dcons, aux_p11_datacons, NULL);
                      printf ("\n");
#endif
                      if (aux_2y_dcons)
                        aux_2y_dcons = ap_abstract0_join (pr->man_dcons, true, aux_2y_dcons, p11_datacons);
                      else
                        {
                          aux_2y_dcons = ap_abstract0_copy (pr->man_dcons, p11_datacons);
                          ap_abstract0_free (pr->man_dcons, p11_datacons);
                        }
                    }
                  /*
                   *  \forall y \in cdim[wl] /\ cdim[i]
                   */

                  /* dy2_econs = add dimensions to econs and constrain dy2 = data(cdim[i])*/

                  dimchange.intdim = 4;
                  dimchange.realdim = 0;
                  dimchange.dim = (ap_dim_t*) malloc (4 * sizeof (ap_dim_t));
                  dimchange.dim[0] = r->datadim + 2 * r->segmentdim;
                  dimchange.dim[1] = r->datadim + 2 * r->segmentdim;
                  dimchange.dim[2] = r->datadim + 2 * r->segmentdim;
                  dimchange.dim[3] = r->datadim + 2 * r->segmentdim;

                  dy2_econs = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimchange, false);
                  free (dimchange.dim);
                  ap_lincons0_array_t arr = ap_lincons0_array_make (2);

                  // l[y2] - l(cdim[0]) + ... + l(cdim[i-1]) ==0
                  arr.p[0].constyp = AP_CONS_EQ;
                  arr.p[0].linexpr0 =
                          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
                  arr.p[0].scalar = NULL;
                  for (j = 0; j < i; j++)
                    {
                      li = r->datadim + r->segmentdim + cdim[j];
                      ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, li, 1);
                    }
                  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ly2, -1);
                  //ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

                  // d(y2) - d(cdim[i]) ==0
                  arr.p[1].constyp = AP_CONS_EQ;
                  arr.p[1].linexpr0 =
                          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
                  arr.p[1].scalar = NULL;
                  ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                                    r->datadim + cdim[i], 1);
                  ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy2, -1);
                  //ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);

                  dy2_econs = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                                               true, dy2_econs, &arr);

                  ap_abstract0_t *data_wl_i_econs = ap_abstract0_copy (pr->man_dcons, dy2_econs);

                  border = ap_abstract0_meet (pr->man_dcons, true, dy2_econs, border);
                  ap_lincons0_array_clear (&arr);
                  /**/
#ifndef NDEBUG2
                  fprintf (stdout, "\n border_left_i: %zu", i);
                  if (border) ap_abstract0_fprint (stdout, pr->man_dcons, border, NULL);
                  fprintf (stdout, "\n");
                  fflush (stdout);
#endif


                  if (aux_2y_dcons)
                    aux_2y_dcons = ap_abstract0_join (pr->man_dcons, true, border, aux_2y_dcons);
                  else
                    {
                      aux_2y_dcons = ap_abstract0_copy (pr->man_dcons, border);
                      ap_abstract0_free (pr->man_dcons, border);
                    }
                  border = NULL;

                  /*
                   *  y1= cdim[wl] y2 = cdim[i] when wl != head of the segment
                   *
                   * */
                  if (wl != 0)
                    {
                      arr = ap_lincons0_array_make (2);

                      // l[y1] - l(cdim[0]) + ... + l(cdim[i-1]) ==0
                      arr.p[0].constyp = AP_CONS_EQ;
                      arr.p[0].linexpr0 =
                              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
                      arr.p[0].scalar = NULL;
                      for (j = 0; j < wl; j++)
                        {
                          li = r->datadim + r->segmentdim + cdim[j];
                          ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, li, 1);
                        }
                      ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ly1, -1);
                      //ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

                      // d(y1) - d(cdim[i]) ==0
                      arr.p[1].constyp = AP_CONS_EQ;
                      arr.p[1].linexpr0 =
                              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
                      arr.p[1].scalar = NULL;
                      ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                                        r->datadim + cdim[wl], 1);
                      ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, dy1, -1);
                      //ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);

                      data_wl_i_econs = ap_abstract0_meet_lincons_array (pr->man_dcons,
                                                                         true, data_wl_i_econs, &arr);

                      ap_lincons0_array_clear (&arr);

#ifndef NDEBUG2
                      fprintf (stdout, "\n data_wl_i_econs: %zu", i);
                      if (aux_2y_dcons) ap_abstract0_fprint (stdout, pr->man_dcons, aux_2y_dcons, NULL);
                      fprintf (stdout, "\n");
                      fflush (stdout);
#endif
                      if (aux_2y_dcons)
                        aux_2y_dcons = ap_abstract0_join (pr->man_dcons, true, data_wl_i_econs, aux_2y_dcons);
                      else
                        {
                          aux_2y_dcons = ap_abstract0_copy (pr->man_dcons, data_wl_i_econs);
                          ap_abstract0_free (pr->man_dcons, data_wl_i_econs);
                        }
                      data_wl_i_econs = NULL;
                    }
                  else
                    {
                      ap_abstract0_free (pr->man_dcons, data_wl_i_econs);
                      data_wl_i_econs = NULL;
                    }
                } // end found_left_dcons->dcons!=NULL
            }// end found_left_dcons!=NULL, i.e. end join with (y,cdim[i]) with y\in cdim[i-1]
        }// end for wl<=i-1


#ifndef NDEBUG1
      printf ("\n concat 3: %zu", i);
      if (aux_2y_dcons) ap_abstract0_fprint (stdout, pr->man_dcons, aux_2y_dcons, NULL);
      printf ("\n");
#endif
      /* 4. join with (cdim[i],cdim[i]) */
      /* 5. join with y1=i dy1=data(cdim[i])*/
      ap_abstract0_t* dii_econs;

      dimchange.intdim = 4;
      dimchange.realdim = 0;
      dimchange.dim = (ap_dim_t*) malloc (4 * sizeof (ap_dim_t));
      dimchange.dim[0] = r->datadim + 2 * r->segmentdim;
      dimchange.dim[1] = r->datadim + 2 * r->segmentdim;
      dimchange.dim[2] = r->datadim + 2 * r->segmentdim;
      dimchange.dim[3] = r->datadim + 2 * r->segmentdim;

      dii_econs = ap_abstract0_copy (pr->man_dcons, r->econs);
      dii_econs = ap_abstract0_add_dimensions (pr->man_dcons, true,
                                               dii_econs, &dimchange, false);
      free (dimchange.dim);
      ap_lincons0_array_t arr = ap_lincons0_array_make (4);

      // l[y1] - l(cdim[0]) + ... + l(cdim[i-1]) ==0
      // l[y2] - l(cdim[0]) + ... + l(cdim[i-1]) ==0

      arr.p[0].constyp = AP_CONS_EQ;
      arr.p[0].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
      arr.p[0].scalar = NULL;
      arr.p[1].constyp = AP_CONS_EQ;
      arr.p[1].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
      arr.p[1].scalar = NULL;
      for (j = 0; j < i; j++)
        {
          li = r->datadim + r->segmentdim + cdim[j];
          ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, li, 1);
          ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, li, 1);
        }
      ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ly1, -1);
      ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0, ly2, -1);
      ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);
      ap_linexpr0_set_cst_scalar_int (arr.p[1].linexpr0, 0);


      // d(y1) - d(cdim[i]) ==0
      // d(y2) - d(cdim[i]) ==0

      arr.p[2].constyp = AP_CONS_EQ;
      arr.p[2].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
      arr.p[2].scalar = NULL;

      ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0,
                                        r->datadim + cdim[i], 1);
      ap_linexpr0_set_coeff_scalar_int (arr.p[2].linexpr0, dy1, -1);
      ap_linexpr0_set_cst_scalar_int (arr.p[2].linexpr0, 0);

      arr.p[3].constyp = AP_CONS_EQ;
      arr.p[3].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
      arr.p[3].scalar = NULL;

      ap_linexpr0_set_coeff_scalar_int (arr.p[3].linexpr0,
                                        r->datadim + cdim[i], 1);
      ap_linexpr0_set_coeff_scalar_int (arr.p[3].linexpr0, dy2, -1);
      ap_linexpr0_set_cst_scalar_int (arr.p[3].linexpr0, 0);


      dii_econs = ap_abstract0_meet_lincons_array (pr->man_dcons, true, dii_econs, &arr);


      ap_lincons0_array_clear (&arr);

      if (aux_2y_dcons)
        aux_2y_dcons = ap_abstract0_join (pr->man_dcons, true,
                                          dii_econs, aux_2y_dcons);
      else
        {
          aux_2y_dcons = ap_abstract0_copy (pr->man_dcons, dii_econs);
          ap_abstract0_free (pr->man_dcons, dii_econs);
        }
      dii_econs = NULL;
      // end 4. join with (cdim[i],cdim[i])


#ifndef NDEBUG1
      printf ("\n concat 4.%zu", i);
      if (aux_2y_dcons) ap_abstract0_fprint (stdout, pr->man_dcons, aux_2y_dcons, NULL);
      printf ("\n");
#endif
    }//end for
  // end concat with pattern succ_1
#ifndef NDEBUG1
  fprintf (stdout, "\n \t concat P12 returns \n ");
  if (aux_2y_dcons) ap_abstract0_fprint (stdout, pr->man_dcons, aux_2y_dcons, NULL);
  fprintf (stdout, "\n");
  ucons_fprint (stdout, pr->man, r, NULL);
  fflush (stdout);
#endif

  /*y2 - y1 >= 0 */
  ap_lincons0_array_t arr = ap_lincons0_array_make (1);
  arr.p[0].constyp = AP_CONS_SUPEQ;
  arr.p[0].linexpr0 =
          ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
  arr.p[0].scalar = NULL;
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                    ly1, -1);
  ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0,
                                    ly2, 1);

  aux_2y_dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, aux_2y_dcons, &arr);
  ap_lincons0_array_clear (&arr);

  return aux_2y_dcons;
}

/*
 * the function concatenates the segments from cdim preserving information
 * w.r.t the patterns \forall y and \forall y1,y2. y1 < y2
 */

ucons_t*
concat_nodes_2 (ucons_internal_t *pr, ucons_t * r, ap_dim_t* cdim,
                size_t size_cdim, bool update_lenght)
{

  size_t i, j, k;
  ap_abstract0_t * aux_dcons = NULL; // the abstract attached to cdim[0] with P1
  ap_abstract0_t * aux_2y_dcons = NULL; // the abstract attached to cdim[0] with P12
  ap_abstract0_t * dy1_econs;
  ap_abstract0_t * fy_dcons;
  ap_abstract0_t * fy2_dcons;

  pattern_t * n_dcons;
  pattern_t *found2_dcons, *found_dcons;
  size_t contor, u_seg, u_seg2;
  ap_dim_t ly1, dy1, ly2, dy2;

  ap_abstract0_t * ni_dcons;
  ap_dim_t li, ly, di, dy;

  u_seg = 1;
  u_seg2 = 1;
  ap_dimchange_t dimchange;
  pattern_key_t * look, *look2;

  arg_assert (r, return NULL;);
  if (!(cdim && size_cdim > 1)) return r;

  dy1 = r->datadim + 2 * r->segmentdim + 2;
  dy2 = r->datadim + 2 * r->segmentdim + 3;

#ifndef NDEBUG1
  printf ("\n CONCAT: ");
  ucons_fprint (stdout, pr->man, r, NULL);
  printf ("\n");
#endif
  /* concat segments w.r.t. the pattern \forall y for the node cdim[0] from cdim[1..size_cdim]
   * aux_dcons is the resulting data constraint */
  aux_dcons = concat_P11 (pr, r, cdim, size_cdim);
  if (aux_dcons == NULL)
    {
      pattern_key_t *looka;
      pattern_t *aux;
      checked_malloc (looka, pattern_key_t,
                      sizeof (pattern_key_t) + u_seg2 * sizeof (size_t), 1, return NULL;);
      memset (looka, 0, sizeof (pattern_key_t) + u_seg2 * sizeof (size_t));
      looka->type = 0;
      looka->segments[0] = cdim[0];
      unsigned keylen = u_seg * sizeof (size_t) + sizeof (pattern_key_t);
      HASH_FIND (hh, r->udcons, looka, keylen, aux);
      if (aux)
        {
          HASH_DEL (r->udcons, aux);
          free (aux);
          remove_pattern_n2p (pr, r, looka);
        }
      looka->type = 3;
      HASH_FIND (hh, r->udcons, looka, keylen, aux);
      if (aux)
        {
          HASH_DEL (r->udcons, aux);
          free (aux);
          remove_pattern_n2p (pr, r, looka);
        }
      /* stop folding; missing information */
#ifndef NDEBUG1
      printf ("\n CONCAT returns: ");
      ucons_fprint (stdout, pr->man, r, NULL);
      printf ("\n");
#endif
      //goto update_lenghts
      if (update_lenght) update_lenghts (pr, r, cdim, size_cdim);
      return r;
    }
  aux_2y_dcons = concat_P12 (pr, r, cdim, size_cdim);
  if (aux_2y_dcons == NULL)
    {
      /* stop folding; missing information */
      pattern_key_t *looka;
      pattern_t *aux;
      checked_malloc (looka, pattern_key_t,
                      sizeof (pattern_key_t) + u_seg2 * sizeof (size_t), 1, return NULL;);
      memset (looka, 0, sizeof (pattern_key_t) + u_seg2 * sizeof (size_t));
      unsigned keylen = u_seg * sizeof (size_t) + sizeof (pattern_key_t);
      looka->type = 3;
      looka->segments[0] = cdim[0];
      HASH_FIND (hh, r->udcons, looka, keylen, aux);
      if (aux)
        {
          HASH_DEL (r->udcons, aux);
          free (aux);
          remove_pattern_n2p (pr, r, looka);
        }
#ifndef NDEBUG1
      printf ("\n CONCAT returns: ");
      ucons_fprint (stdout, pr->man, r, NULL);
      printf ("\n");
#endif
      if (update_lenght) update_lenghts (pr, r, cdim, size_cdim);
      return r;
    }

  u_seg = pr->PI[0].u_seg;
  checked_malloc (look, pattern_key_t,
                  sizeof (pattern_key_t) + u_seg * sizeof (size_t), 1, return NULL;);
  memset (look, 0, sizeof (pattern_key_t) + u_seg * sizeof (size_t));
  look->type = 0;
  unsigned keylen = u_seg * sizeof (size_t) + sizeof (pattern_key_t);

  /* concat segments w.r.t. the pattern \forall y1,y2. y1<y2 for the node gdim[0] from gdim[1..size_gdim] */
  /* aux_2y_dcons is the resulting data constraint */
  u_seg2 = pr->PI[3].u_seg;
  checked_malloc (look2, pattern_key_t,
                  sizeof (pattern_key_t) + u_seg2 * sizeof (size_t), 1, return NULL;);
  memset (look2, 0, sizeof (pattern_key_t) + u_seg2 * sizeof (size_t));
  look2->type = 3;
  unsigned keylen2 = u_seg2 * sizeof (size_t) + sizeof (pattern_key_t);

  /* pt \forall y inlocuirea vechii constrangeri n_dcons->dcons cu cea nou calculata aux_dcons */
  look->segments[0] = cdim[0];
  HASH_FIND (hh, r->udcons, look, keylen, n_dcons);

  if (n_dcons)
    {
      if (n_dcons->dcons)
        {
          if (aux_dcons)
            n_dcons->dcons = ap_abstract0_join (pr->man_dcons, false,
                                                aux_dcons, n_dcons->dcons);
        }
      else
        {
          if (!test_singleton (pr->man_dcons, r->econs, r->datadim, r->segmentdim, cdim[0]))
            {

              remove_pattern_n2p (pr, r, look);
              HASH_DEL (r->udcons, n_dcons);
              free (n_dcons);
              if (update_lenght) update_lenghts (pr, r, cdim, size_cdim);
              return r;
            }
          else
            {
              if (aux_dcons == NULL)
                {
                  remove_pattern_n2p (pr, r, look);
                  HASH_DEL (r->udcons, n_dcons);
                  free (n_dcons);
                  if (update_lenght) update_lenghts (pr, r, cdim, size_cdim);
                  return r;
                }
              else n_dcons->dcons = ap_abstract0_copy (pr->man_dcons, aux_dcons);
            }
        }
    }
  else
    {
      if (!test_singleton (pr->man_dcons, r->econs, r->datadim, r->segmentdim, cdim[0]))
        {
          if (update_lenght) update_lenghts (pr, r, cdim, size_cdim);
          return r;
        }
      else
        {
          checked_malloc (n_dcons, pattern_t, 1, sizeof (pattern_t) + u_seg * sizeof (size_t), return NULL;);
          memset (n_dcons, 0, sizeof (pattern_t) + u_seg * sizeof (size_t));
          n_dcons->key.type = look->type;
          for (size_t i = 0; i < (u_seg); i++)
            n_dcons->key.segments[i] = look->segments[i];
          n_dcons->dcons = ap_abstract0_copy (pr->man_dcons, aux_dcons);
          HASH_ADD (hh, r->udcons, key, keylen, n_dcons);
          r = add_pattern_n2p (pr, r, look);
        }
    }
  /*************************************/

  /* pt \forall y1,y2. y1<y2 inlocuirea vechii constrangeri found2_dcons->dcons cu cea nou calculata aux_2y_dcons */
  look2->segments[0] = cdim[0];
  HASH_FIND (hh, r->udcons, look2, keylen2, found2_dcons);

  if (found2_dcons)
    {
      if (found2_dcons->dcons)
        {
          if (aux_2y_dcons)
            found2_dcons->dcons = ap_abstract0_join (pr->man_dcons, false, aux_2y_dcons, found2_dcons->dcons);
        }
      else
        {
          if (!test_singleton (pr->man_dcons, r->econs, r->datadim, r->segmentdim, cdim[0]))
            {

              remove_pattern_n2p (pr, r, look2);
              HASH_DEL (r->udcons, found2_dcons);
              free (found2_dcons);
              if (update_lenght) update_lenghts (pr, r, cdim, size_cdim);
              return r;
            }
          else
            {
              if (aux_2y_dcons == NULL)
                {
                  remove_pattern_n2p (pr, r, look2);
                  HASH_DEL (r->udcons, found2_dcons);
                  free (found2_dcons);
                  if (update_lenght) update_lenghts (pr, r, cdim, size_cdim);
                  return r;
                }
              else found2_dcons->dcons = ap_abstract0_copy (pr->man_dcons, aux_2y_dcons);
            }
        }
    }
  else
    {
      if (!test_singleton (pr->man_dcons, r->econs, r->datadim, r->segmentdim, cdim[0]))
        {
          if (update_lenght) update_lenghts (pr, r, cdim, size_cdim);
          return r;
        }
      else
        {
          checked_malloc (found2_dcons, pattern_t, 1, sizeof (pattern_t) + u_seg * sizeof (size_t), return NULL;);
          memset (found2_dcons, 0, sizeof (pattern_t) + u_seg * sizeof (size_t));
          found2_dcons->key.type = look2->type;
          for (size_t i = 0; i < (u_seg); i++)
            found2_dcons->key.segments[i] = look2->segments[i];
          found2_dcons->dcons = ap_abstract0_copy (pr->man_dcons, aux_2y_dcons);
          ;
          HASH_ADD (hh, r->udcons, key, keylen2, found2_dcons);
          r = add_pattern_n2p (pr, r, look2);
        }
    }
  /*************************************/

  if (aux_dcons) ap_abstract0_free (pr->man_dcons, aux_dcons);
  if (aux_2y_dcons) ap_abstract0_free (pr->man_dcons, aux_2y_dcons);
  //end concat for pattern with succ_1

#ifndef NDEBUG1
  printf ("\n concat_nodes_2 before recalc lengths \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  printf ("\n");
#endif

  /* recalculate length of tdim[0] */
  if (update_lenght)
    {
      ap_dim_t l0;
      /* update length for the merged segments */
      /* l(cdim[0]) = l(cdim[0]) + l(cdim[1]) + ... + l(cdim[last])*/

      /* for econs */
      ap_linexpr0_t *expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
      ap_linexpr0_set_cst_scalar_int (expr, 0);
      for (i = 0; i < size_cdim; i++)
        {
          li = r->datadim + r->segmentdim + cdim[i];
          ap_linexpr0_set_coeff_scalar_int (expr, li, 1);
        }
      /* for udcons */
      ap_linexpr0_t *expr_1y = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                  r->datadim + 2 * r->segmentdim + 2);
      ap_linexpr0_set_cst_scalar_int (expr_1y, 0);
      for (i = 0; i < size_cdim; i++)
        {
          li = r->datadim + r->segmentdim + cdim[i];
          ap_linexpr0_set_coeff_scalar_int (expr_1y, li, 1);
        }
      ap_linexpr0_t *expr_2y = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                  r->datadim + 2 * r->segmentdim + 4);
      ap_linexpr0_set_cst_scalar_int (expr_2y, 0);
      for (i = 0; i < size_cdim; i++)
        {
          li = r->datadim + r->segmentdim + cdim[i];
          ap_linexpr0_set_coeff_scalar_int (expr_2y, li, 1);
        }
      /* TODO allocate linexpr for each pattern type */

      l0 = r->datadim + r->segmentdim + cdim[0];
      r->econs = ap_abstract0_assign_linexpr (pr->man_dcons, true,
                                              r->econs, l0, expr, NULL);
      ap_linexpr0_free (expr);

      pattern_t *s;
      size_t nr_y;
      for (s = r->udcons; s != NULL; s = s->hh.next)
        {
#ifndef NDEBUG1
          fprintf (stdout, "\n pattern_type := %zu", s->key.type);
          fflush (stdout);
#endif
          nr_y = pr->PI[s->key.type].nr_y;
          if (nr_y == 1)
            {
#ifndef NDEBUG1
              ap_linexpr0_fprint (stdout, expr_1y, NULL);
              ap_abstract0_fprint (stdout, pr->man_dcons, s->dcons, NULL);
              fprintf (stdout, "\n 1 nr_y = %zu,l0=%zu\n", nr_y, l0);
              fflush (stdout);
#endif
              s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true,
                                                      s->dcons, l0, expr_1y, NULL);
            }
          if (nr_y == 2)
            {
#ifndef NDEBUG1
              ap_linexpr0_fprint (stdout, expr_2y, NULL);
              ap_abstract0_fprint (stdout, pr->man_dcons, s->dcons, NULL);
              fprintf (stdout, "\n 2 nr_y = %zu, l0=%zu\n", nr_y, l0);
              fflush (stdout);
#endif
              s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons, l0, expr_2y, NULL);

            }
        }
      ap_linexpr0_free (expr_1y);
      ap_linexpr0_free (expr_2y);
    }
#ifndef NDEBUG2
  printf ("\n concat_nodes_2 returs \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  printf ("\n");
#endif
  return r;

}

void
update_lenghts (ucons_internal_t *pr, ucons_t *r, ap_dim_t* cdim, size_t size_cdim)
{
  ap_dim_t l0, li;
  size_t i;
  /* update length for the merged segments */
  /* l(cdim[0]) = l(cdim[0]) + l(cdim[1]) + ... + l(cdim[last])*/

  /* for econs */
  ap_linexpr0_t *expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
  ap_linexpr0_set_cst_scalar_int (expr, 0);
  for (i = 0; i < size_cdim; i++)
    {
      li = r->datadim + r->segmentdim + cdim[i];
      ap_linexpr0_set_coeff_scalar_int (expr, li, 1);
    }
  /* for udcons */
  ap_linexpr0_t *expr_1y = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                              r->datadim + 2 * r->segmentdim + 2);
  ap_linexpr0_set_cst_scalar_int (expr_1y, 0);
  for (i = 0; i < size_cdim; i++)
    {
      li = r->datadim + r->segmentdim + cdim[i];
      ap_linexpr0_set_coeff_scalar_int (expr_1y, li, 1);
    }
  ap_linexpr0_t *expr_2y = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                              r->datadim + 2 * r->segmentdim + 4);
  ap_linexpr0_set_cst_scalar_int (expr_2y, 0);
  for (i = 0; i < size_cdim; i++)
    {
      li = r->datadim + r->segmentdim + cdim[i];
      ap_linexpr0_set_coeff_scalar_int (expr_2y, li, 1);
    }
  /* TODO allocate linexpr for each pattern type */

  l0 = r->datadim + r->segmentdim + cdim[0];

  r->econs = ap_abstract0_assign_linexpr (pr->man_dcons, true,
                                          r->econs, l0, expr, NULL);

  ap_linexpr0_free (expr);

  pattern_t *s;
  size_t nr_y;
  for (s = r->udcons; s != NULL; s = s->hh.next)
    {
#ifndef NDEBUG1
      fprintf (stdout, "\n pattern_type := %zu", s->key.type);
      fflush (stdout);
#endif
      nr_y = pr->PI[s->key.type].nr_y;
      if (nr_y == 1)
        {
#ifndef NDEBUG1
          ap_linexpr0_fprint (stdout, expr_1y, NULL);
          ap_abstract0_fprint (stdout, pr->man_dcons, s->dcons, NULL);
          fprintf (stdout, "\n 1 nr_y = %zu,l0=%zu\n", nr_y, l0);
          fflush (stdout);
#endif
          s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true,
                                                  s->dcons, l0, expr_1y, NULL);
        }
      if (nr_y == 2)
        {
#ifndef NDEBUG1
          ap_linexpr0_fprint (stdout, expr_2y, NULL);
          ap_abstract0_fprint (stdout, pr->man_dcons, s->dcons, NULL);
          fprintf (stdout, "\n 2 nr_y = %zu, l0=%zu\n", nr_y, l0);
          fflush (stdout);
#endif
          s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons, l0, expr_2y, NULL);

        }
    }
  ap_linexpr0_free (expr_1y);
  ap_linexpr0_free (expr_2y);
}

ucons_t*
concat_nodes_1 (ucons_internal_t *pr, ucons_t * r, ap_dim_t* cdim,
                size_t size_cdim, bool update_lenght)
{

  size_t i, j;
  ap_abstract0_t * aux_dcons = NULL;
  ap_abstract0_t * aux_2y_dcons = NULL;
  ap_abstract0_t * dy1_econs;
  pattern_t * n_dcons;
  pattern_t *found2_dcons, *found_dcons;
  size_t contor, u_seg, u_seg2;
  ap_dim_t ly1, dy1;

  ap_abstract0_t * ni_dcons;
  ap_dim_t li, ly, di, dy;

  ap_dimchange_t dimchange;
  pattern_key_t * look, *look2;

  arg_assert (r, return NULL;);
  arg_assert (cdim && size_cdim > 1, return r;);

#ifndef NDEBUG2
  printf ("  \n\t concat_nodes_1 : \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  printf ("\n");
#endif

  /* concat segments w.r.t. the pattern \forall y for the node cdim[0] from cdim[1..size_cdim]
   * aux_dcons is the resulting data constraint */
  u_seg = pr->PI[0].u_seg;
  checked_malloc (look, pattern_key_t, sizeof (pattern_key_t) + u_seg * sizeof (size_t), 1, return NULL;);
  memset (look, 0, sizeof (pattern_key_t) + u_seg * sizeof (size_t));
  look->type = 0;
  unsigned keylen = u_seg * sizeof (size_t) + sizeof (pattern_key_t);
  //new

  //aux_dcons=ap_abstract0_copy(pr->man_dcons,n_dcons->dcons);
  for (i = 1; i < size_cdim; i++)
    {
      look->segments[0] = cdim[i];
      HASH_FIND (hh, r->udcons, look, keylen, n_dcons);


      if (!n_dcons)
        {
          if (!test_singleton (pr->man_dcons, r->econs, r->datadim, r->segmentdim, cdim[i]))
            {
#ifndef NDEBUG2
              printf ("  \n\t concat_nodes_1 return with pattern info insufficient: \n");
              ucons_fprint (stdout, pr->man, r, NULL);
              printf ("\n");
#endif
              if (update_lenght)
                update_lenghts (pr, r, cdim, size_cdim);
              return r;
            }
        }

      if (n_dcons)
        {
          /* TODO add element corresponding to the existential */
          /* join with the universal from cdim[i] not needed for this pattern */
          if (n_dcons->dcons)
            {

              /* */
              ap_linexpr0_t *expr_y = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
              ap_linexpr0_set_cst_scalar_int (expr_y, 0); /*TODO check the definition of the sustitution */
              for (j = 0; j < i; j++)
                {
                  li = r->datadim + r->segmentdim + cdim[j];
                  ap_linexpr0_set_coeff_scalar_int (expr_y, li, -1); /*TODO check the definition of the sustitution */
                }
              ly = r->datadim + 2 * r->segmentdim;
              ap_linexpr0_set_coeff_scalar_int (expr_y, ly, 1);
              n_dcons->dcons = ap_abstract0_substitute_linexpr (pr->man_dcons, true, n_dcons->dcons, ly, expr_y, NULL);
              if (aux_dcons)
                aux_dcons = ap_abstract0_join (pr->man_dcons, true, aux_dcons, n_dcons->dcons);
              else aux_dcons = ap_abstract0_copy (pr->man_dcons, n_dcons->dcons);

              ap_linexpr0_free (expr_y);
            }
          else
            {
#ifndef NDEBUG2
              fprintf (stdout, "  \n\t concat_nodes_1: \n");
              fprintf (stdout, "P(n%zu)==>null ", i);
              fprintf (stdout, "\n");
              fflush (stdout);
#endif
            }
        }
#ifndef NDEBUG1
      fprintf (stdout, "  \n\t concat_nodes_1 aux_dcons: \n");
      if (aux_dcons) ap_abstract0_fprint (stdout, pr->man_dcons, aux_dcons, NULL);
      fprintf (stdout, "\n");
      fflush (stdout);
#endif

      /* ni_dcons  =  econs + y == l[cdim[0]] + ... l[cdim[i]]  + d(y) = d(cdim[i]) */

      dimchange.intdim = 2;
      dimchange.realdim = 0;
      dimchange.dim = (ap_dim_t *) malloc (2 * sizeof (ap_dim_t));
      dimchange.dim[0] = r->datadim + 2 * r->segmentdim;
      dimchange.dim[1] = r->datadim + 2 * r->segmentdim;

      ni_dcons = ap_abstract0_add_dimensions (pr->man_dcons, false, r->econs, &dimchange, false);
      free (dimchange.dim);


      ap_lincons0_array_t arr = ap_lincons0_array_make (2);
      // l[y] == l[cdim[0]] + ... l[cdim[i-1]]
      arr.p[0].constyp = AP_CONS_EQ;
      arr.p[0].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
      arr.p[0].scalar = NULL;

      for (j = 0; j < i; j++)
        {
          li = r->datadim + r->segmentdim + cdim[j];
          ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, li, 1);
        }

      ly = r->datadim + 2 * r->segmentdim;
      ap_linexpr0_set_coeff_scalar_int (arr.p[0].linexpr0, ly, -1);
      ap_linexpr0_set_cst_scalar_int (arr.p[0].linexpr0, 0);

      // d(y) - d(cdim[i]) ==0
      arr.p[1].constyp = AP_CONS_EQ;
      arr.p[1].linexpr0 =
              ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
      arr.p[1].scalar = NULL;
      ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                        r->datadim + cdim[i], 1);
      ap_linexpr0_set_coeff_scalar_int (arr.p[1].linexpr0,
                                        r->datadim + 2 * r->segmentdim + 1, -1);

      ni_dcons = ap_abstract0_meet_lincons_array (pr->man_dcons, true, ni_dcons, &arr);

      ap_lincons0_array_clear (&arr);

      if (aux_dcons)
        aux_dcons = ap_abstract0_join (pr->man_dcons, true, aux_dcons, ni_dcons);
      else aux_dcons = ap_abstract0_copy (pr->man_dcons, ni_dcons);


    }

#ifndef NDEBUG1
  fprintf (stdout, "  \n\t concat_nodes_1 aux_dcons: \n");
  if (aux_dcons) ap_abstract0_fprint (stdout, pr->man_dcons, aux_dcons, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif
  /*end computation for the first pattern */


  /* pt \forall y inlocuirea vechii constrangeri n_dcons->dcons cu cea nou calculata aux_dcons */
  look->segments[0] = cdim[0];
  HASH_FIND (hh, r->udcons, look, keylen, n_dcons);
  if (n_dcons)
    {
      if (aux_dcons == NULL)
        {
          /* eliminated from the table since over-aproximated by true */
          remove_pattern_n2p (pr, r, look);
          HASH_DEL (r->udcons, n_dcons);
          free (n_dcons);
          n_dcons = NULL;
          // do not return before all modifications are done on exitentials
          // return r;
        }
      else
        {
          /* combine with property on initial cdim[0]*/
          if (n_dcons->dcons)
            {
#ifndef NDEBUG1
              fprintf (stdout, "  \n\t concat_nodes_1 n1 not singleton found in hash join with  aux_dcons: \n");
              fflush (stdout);
#endif
              n_dcons->dcons = ap_abstract0_join (pr->man_dcons, false, aux_dcons, n_dcons->dcons);
            }
          else
            {
              if (!test_singleton (pr->man_dcons, r->econs, r->datadim, r->segmentdim, cdim[0]))
                {
                  remove_pattern_n2p (pr, r, &n_dcons->key);
                  HASH_DEL (r->udcons, n_dcons);
                  free (n_dcons);
                  n_dcons = NULL;
#ifndef NDEBUG1
                  fprintf (stdout, "  \n\t concat_nodes_1 return with pattern info insufficient: \n");
                  ucons_fprint (stdout, pr->man, r, NULL);
                  fprintf (stdout, "\n");
                  fflush (stdout);
#endif
                }
              else
                {
#ifndef NDEBUG1
                  fprintf (stdout, "  \n\t concat_nodes_1 n1 singleton found in hash copy aux_dcons: \n");
                  fflush (stdout);
#endif
                  n_dcons->dcons = ap_abstract0_copy (pr->man_dcons, aux_dcons);
                }
            }
        }
    }
  else
    {
      if (!test_singleton (pr->man_dcons, r->econs, r->datadim, r->segmentdim, cdim[0]))
        {
#ifndef NDEBUG1
          printf ("  \n\t concat_nodes_1 return with pattern info insufficient: \n");
          ucons_fprint (stdout, pr->man, r, NULL);
          printf ("\n");
#endif
        }
      else
        {
#ifndef NDEBUG1
          fprintf (stdout, "  \n\t concat_nodes_1 n1 singleton aux_dcons: \n");
          fflush (stdout);
#endif
          checked_malloc (n_dcons, pattern_t, 1, sizeof (pattern_t) + u_seg * sizeof (size_t), return NULL;);
          memset (n_dcons, 0, sizeof (pattern_t) + u_seg * sizeof (size_t));
          n_dcons->key.type = look->type;
          for (size_t i = 0; i < (u_seg); i++)
            n_dcons->key.segments[i] = look->segments[i];
          n_dcons->dcons = ap_abstract0_copy (pr->man_dcons, aux_dcons);
          HASH_ADD (hh, r->udcons, key, keylen, n_dcons);
          r = add_pattern_n2p (pr, r, look);
        }
    }

  ap_abstract0_free (pr->man_dcons, aux_dcons);

  /* recalculate length of tdim[0] */
  if (update_lenght)
    {
      ap_dim_t l0;
      /* update length for the merged segments */
      /* l(cdim[0]) = l(cdim[0]) + l(cdim[1]) + ... + l(cdim[last])*/

      /* for econs */
      ap_linexpr0_t *expr = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim);
      ap_linexpr0_set_cst_scalar_int (expr, 0);
      for (i = 0; i < size_cdim; i++)
        {
          li = r->datadim + r->segmentdim + cdim[i];
          ap_linexpr0_set_coeff_scalar_int (expr, li, 1);
        }
      /* for udcons */
      ap_linexpr0_t *expr_1y = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 2);
      ap_linexpr0_set_cst_scalar_int (expr_1y, 0);
      for (i = 0; i < size_cdim; i++)
        {
          li = r->datadim + r->segmentdim + cdim[i];
          ap_linexpr0_set_coeff_scalar_int (expr_1y, li, 1);
        }
      ap_linexpr0_t *expr_2y = ap_linexpr0_alloc (AP_LINEXPR_DENSE, r->datadim + 2 * r->segmentdim + 4);
      ap_linexpr0_set_cst_scalar_int (expr_2y, 0);
      for (i = 0; i < size_cdim; i++)
        {
          li = r->datadim + r->segmentdim + cdim[i];
          ap_linexpr0_set_coeff_scalar_int (expr_2y, li, 1);
        }
      /* TODO allocate linexpr for each pattern type */

      l0 = r->datadim + r->segmentdim + cdim[0];
      r->econs = ap_abstract0_assign_linexpr (pr->man_dcons, true, r->econs, l0, expr, NULL);

      ap_linexpr0_free (expr);

      pattern_t *s;
      size_t nr_y;
      for (s = r->udcons; s != NULL; s = s->hh.next)
        {
          nr_y = pr->PI[s->key.type].nr_y;
          if (nr_y == 1)
            s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons, l0, expr_1y, NULL);
          if (nr_y == 2)
            s->dcons = ap_abstract0_assign_linexpr (pr->man_dcons, true, s->dcons, l0, expr_2y, NULL);
        }
    }

#ifndef NDEBUG1
  printf ("  \n\t concat_nodes_1 returns: \n");
  ucons_fprint (stdout, pr->man, r, NULL);
  printf ("\n");
#endif
  return r;
}




