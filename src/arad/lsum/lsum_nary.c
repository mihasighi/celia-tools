/**************************************************************************/
/*                                                                        */
/*  CELIA Tools / LSUM Abstract Domain                                    */
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


#include "lsum.h"
#include "lsum_internal.h"
#include "sh_macros.h"
#include "ap_generic.h"
#include "ap_abstract0.h"

/* ============================================================ */
/* Meet and Join */
/* ============================================================ */

lsum_t *
lsum_meet (ap_manager_t * man, bool destructive, lsum_t * a1, lsum_t * a2)
{
  if (lsum_is_bottom (man, a1) || lsum_is_bottom (man, a2))
    return NULL;
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_MEET, 0);
  arg_assert (a1->datadim == a2->datadim
              && a1->segmdim == a2->segmdim, return NULL;
    );
  lsum_t *r = lsum_alloc_internal (pr, a1->datadim, a1->segmdim);
  r->dcons =
    ap_abstract0_meet (pr->man_dcons, destructive, a1->dcons, a2->dcons);
  if (destructive)
    {
      a1->dcons = a2->dcons = NULL;
      free (a1);
      free (a2);
    }
  return r;
}

lsum_t *
lsum_join (ap_manager_t * man, bool destructive, lsum_t * a1, lsum_t * a2)
{
  if (lsum_is_bottom (man, a1) || lsum_is_bottom (man, a2))
    return NULL;
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_JOIN, 0);
  arg_assert (a1->datadim == a2->datadim
              && a1->segmdim == a2->segmdim, return NULL;
    );
  lsum_t *r = lsum_alloc_internal (pr, a1->datadim, a1->segmdim);
  r->dcons =
    ap_abstract0_join (pr->man_dcons, destructive, a1->dcons, a2->dcons);
  if (destructive)
    {
      a1->dcons = a2->dcons = NULL;
      free (a1);
      free (a2);
    }
  return r;
}

lsum_t *
lsum_meet_array (ap_manager_t * man, lsum_t ** tab, size_t size)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_MEET_ARRAY, 0);
  arg_assert (size > 0, return NULL;
    );
  lsum_t *r = lsum_copy_internal (pr, tab[0]);
  size_t i;
  for (i = 1; i < size && !lsum_is_bottom (man, r); i++)
    {
      lsum_t *rr = lsum_meet (man, false, tab[i], r);
      lsum_free (man, r);
      r = rr;
    }
  return r;
}

lsum_t *
lsum_join_array (ap_manager_t * man, lsum_t ** tab, size_t size)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_JOIN_ARRAY, 0);
  arg_assert (size > 0, return NULL;
    );
  lsum_t *r = lsum_copy_internal (pr, tab[0]);
  size_t i;
  for (i = 1; i < size; i++)
    {
      lsum_t *rr = lsum_join (man, false, tab[i], r);
      lsum_free (man, r);
      r = rr;
    }
  return r;
}


/* ============================================================ */
/* Widening, Narrowing */
/* ============================================================ */

lsum_t *
lsum_widening (ap_manager_t * man, lsum_t * a1, lsum_t * a2)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_WIDENING, 0);
  arg_assert (a1 && a2 && a1->datadim == a2->datadim
              && a1->segmdim == a2->segmdim, return NULL;
    );
  lsum_t *r = lsum_alloc_internal (pr, a1->datadim, a1->segmdim);
  r->dcons = ap_abstract0_widening (pr->man_dcons, a1->dcons, a2->dcons);
#ifndef NDEBUG1
  fprintf (stdout, "\n????? lsum_widening: with a1=(\n");
  ap_abstract0_fprint (stdout, pr->man_dcons, a1->dcons, NULL);
  fprintf (stdout, ") widen a2=(\n");
  ap_abstract0_fprint (stdout, pr->man_dcons, a2->dcons, NULL);
  fprintf (stdout, ") returns r=(\n");
  ap_abstract0_fprint (stdout, pr->man_dcons, r->dcons, NULL);
  fprintf (stdout, ") \n");
  fflush (stdout);

#endif
  return r;
}

lsum_t *
lsum_widening_threshold (ap_manager_t * man,
                         lsum_t * a1, lsum_t * a2,
                         ap_lincons0_array_t * array)
{
  lsum_internal_t *pr =
    lsum_init_from_manager (man, AP_FUNID_WIDENING, array->size + 1);
  arg_assert (a1 && a2 && a1->datadim == a2->datadim
              && a1->segmdim == a2->segmdim, return NULL;
    );
  lsum_t *r = lsum_alloc_internal (pr, a1->datadim, a1->segmdim);
  r->dcons =
    ap_abstract0_widening_threshold (pr->man_dcons, a1->dcons, a2->dcons,
                                     array);
#ifndef NDEBUG1
#endif
  return r;
}

/* NOT IMPLEMENTED */
lsum_t *
lsum_narrowing (ap_manager_t * man, lsum_t * a1, lsum_t * a2)
{
  lsum_internal_t *pr = lsum_init_from_manager (man, AP_FUNID_WIDENING, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  if (a1 != a1)
    return NULL;                /* to remove warning on unsed parameter */
  return a2;
}
