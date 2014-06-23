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
 * Predicates on ucons, extraction functions
 */

#include "ucons.h"
#include "ucons_internal.h"
#include "ucons_fun.h"
#include "ap_generic.h"

#include <stdio.h>



/* ============================================================ */
/* Tests */
/* ============================================================ */

/**
 * Tests if the pattern @param{p} can be instantiated in @param{r}.
 * Used in ucons_resize.
 */
bool
is_pattern_inst (ucons_internal_t * pr, ucons_t *r, pattern_key_t *p)
{

  size_t type = p->type;
  size_t useg = pr->PI[type].u_seg;

  ap_dim_t n = p->segments[0];

  bool res = true;

  switch (pr->PI[type].kind)
    {
    case pattern_2_1:
    case pattern_1:
    case pattern_1_l1:
      {
        res = !test_singleton (pr->man_dcons, r->econs,
                               r->datadim, r->segmentdim, n);
        break;
      }
    case pattern_1_lx_1:
      {
        //pattern strengthens econs in split only if it's size is 2
        res = (!test_singleton (pr->man_dcons, r->econs,
                                r->datadim, r->segmentdim, n)
               && test_l_leq_v (pr->man_dcons, r->econs,
                                r->datadim, r->segmentdim, n, 2));
        break;
      }
    case pattern_succ_1_2:
      {
        //res == false is fine also. 
        //not need to strengthen econs with this pattern
        res = test_g_leq_v (pr->man_dcons, r->econs,
                            r->datadim, r->segmentdim, n, 3);
        break;
      }

    default:
      {
        bool flag = true;
        //TODO: test fully the patterns satisfiability when is defined on multiple segment
        for (size_t i = 0; i < (useg) && flag == true; i++)
          {
            n = p->segments[i];
            flag = !test_singleton (pr->man_dcons, r->econs,
                                    r->datadim, r->segmentdim, n);
          }
        res = flag;
      }
    }

  return res;
}

/**
 * Test if bottom value.
 * NB: the bottom value is the empty set. 
 */
bool
ucons_is_bottom (ap_manager_t * man, ucons_t * a)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_IS_BOTTOM, 0);
  return (!a || ap_abstract0_is_bottom (pr->man_dcons, a->econs));
}

/**
 * Test if top value.
 * NB: the top value is top ucons in the first position. 
 */
bool
ucons_is_top (ap_manager_t * man, ucons_t * a)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_IS_TOP, 0);

  // E shall be top
  bool res = a && ap_abstract0_is_top (pr->man_dcons, a->econs);

  // for each P=>U with P not empty, U shall be top
  for (pattern_t *s = a->udcons; s != NULL && res; s = s->hh.next)
    res = (s->dcons != NULL) &&
    ap_abstract0_is_top (pr->man_dcons, s->dcons);

  return res;
}

/**
 * Test inclusion of values.
 */
bool
ucons_is_leq (ap_manager_t * man, ucons_t * a1, ucons_t * a2)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_IS_LEQ, 0);
  bool res = true;
  bool option_sat = false;

#ifndef NDEBUG2
  fprintf (stdout, "\n@@@@ ucons_is_leq: generated problem ");
  ucons_print_smt (man, a1, a2, NULL);
  //      fprintf(stdout,"\t a1 : \n");
  //      ucons_fprint(stdout,pr->man,a1, NULL);
  //      fprintf(stdout,"\n");
  //      fprintf(stdout,"\t implies a2 : \n");
  //      ucons_fprint(stdout,pr->man,a2, NULL);
  fprintf (stdout, "\n");
  fflush (stdout);
#endif

  // Basic cases
  if (ucons_is_bottom (man, a1))
    {
      fprintf (stdout, "\n a1 bottom");
      fflush (stdout);
      return true;
    }
  if (ucons_is_bottom (man, a2))
    {
      fprintf (stdout, "\n a2 bottom");
      fflush (stdout);
      return false;
    }
  if (a1->datadim != a2->datadim || a1->segmentdim != a2->segmentdim)
    return false;

  // SYNTACTICAL TEST STARTS HERE
  // variable used to store saturations of a1
  ucons_t* a1s = ucons_copy_internal (pr, a1);
  // variable used in the saturation to store nodes already saturated
  bool* nsat = NULL;

  if (option_sat)
    {
      // initialize nodes saturated
      nsat = (bool*) malloc (sizeof (bool) * a2->segmentdim);
      for (size_t i = 0; i < a2->segmentdim; i++)
        nsat[i] = false;
      // saturate existentials
      a1s = ucons_saturation (pr, a1s, NULL, 0);
    }

  // Test existentials: a1s.E <= a2.E
  if (!ap_abstract0_is_leq (pr->man_dcons, a1s->econs, a2->econs))
    {
#ifndef NDEBUG
      fprintf (stdout, "\n@@@@ ucons_is_leq: exists not leq, return false.\n");
      fflush (stdout);
#endif
      res = false;
      goto ucons_is_leq_return;
    }

  // Test forall s2 in a2, exists s1 in a1 such that s1(P) = s2(P) and 
  //        (1) s1 != null(top) and (s1(U) <= s2(U) or sat(s1(U)) <= s2(U)) 
  //     or (2) s1 == null(top) and [(2a) s1(E) ==> length(n)<length(P) or
  //                                 (2b) s1(E) ==> length(n)>= length(P) and
  //                                      sat(s1(U)) <= s2(U)]
  for (pattern_t* s2 = a2->udcons; s2 != NULL; s2 = s2->hh.next)
    {
      pattern_t* s1;
      size_t u_seg = pr->PI[s2->key.type].u_seg;
      size_t e_seg = pr->PI[s2->key.type].e_seg;
      unsigned keylen = (u_seg + e_seg) * sizeof (size_t) + sizeof (pattern_key_t);
      HASH_FIND (hh, a1s->udcons, &s2->key.type, keylen, s1);
      // for saturation 
      bool dosat = option_sat;
      if (dosat && s1)
        {
          // check that saturation is not already done for the nodes
          dosat = false;
          for (size_t i = 0; i < u_seg && !dosat; i++)
            if (nsat[s1->key.segments[i]] == false)
              dosat = true;
        }
      // The test
      if (s1 != NULL)
        { // a1s contains some constraint for this pattern pr->PI[s2->key.type]
          if (ap_abstract0_is_leq (pr->man_dcons, s1->dcons, s2->dcons))
            {
              // inclusion satisfied, go to another s2 pattern
              continue;
            }
          // !ap_abstract0_is_leq (pr->man_dcons, s1->dcons, s2->dcons)
          // inclusion not satisfied
          if (!dosat)
            {
              // nothing can be done to obtain entailment
#ifndef NDEBUG2
              fprintf (stdout, "\n@@@@ ucons_is_leq: (1) returns false for pattern type %zu, ",
                       s2->key.type);
              fprintf (stdout, "nodes: ");
              for (size_t i = 0; i < u_seg; i++)
                fprintf (stdout, "n%zu, ", s2->key.segments[i]);
              fprintf (stdout, "\n\t\t because U1 = (");
              ap_abstract0_fprint (stdout, pr->man_dcons, s1->dcons, NULL);
              fprintf (stdout, "\n\t\t) not <= than U2 = (");
              ap_abstract0_fprint (stdout, pr->man_dcons, s2->dcons, NULL);
              fprintf (stdout, "\n\t\t)\n");
              fflush (stdout);
#endif
              res = false;
              goto ucons_is_leq_return;
            }
          // else, saturation at the end
        }
      else // s1 == NULL
        {
          // i.e., a1s does not contain a constraint for this pattern pr->PI[s2->key.type]
          // Then only if the pattern is instantiable 
          // (i.e., in a1s the segments in s2 universal quantified have lengths > 1)
          // return false
          bool sat_pattern = true;
          // TODO: put len correctly for each pattern and each node
          int len = (pr->PI[s2->key.type].kind == pattern_succ_1_2) ? 2 : 1;

          for (size_t i = 0; i < u_seg && sat_pattern; i++)
            {
              ap_linexpr0_t* linexpr = ap_linexpr0_alloc (AP_LINEXPR_DENSE,
                                                          a1->datadim + 2 * a1->segmentdim);
              ap_linexpr0_set_cst_scalar_int (linexpr, len);
              ap_linexpr0_set_coeff_scalar_int (linexpr,
                                                a1->datadim + a1->segmentdim +
                                                s2->key.segments[i], -1);
              ap_lincons0_t cons = ap_lincons0_make (AP_CONS_SUPEQ, linexpr, NULL);
              sat_pattern = !ap_abstract0_sat_lincons (pr->man_dcons, a1s->econs, &cons);
              ap_lincons0_clear (&cons); // free also linexpr
            }

          if (sat_pattern && !ap_abstract0_is_top (pr->man_dcons, s2->dcons))
            {
              // patern can be instantiated, but s1 is top and s2 not!
              // bad, but maybe saturation can do something
              if (!dosat)
                {
#ifndef NDEBUG2
                  fprintf (stdout, "\n@@@@ ucons_is_leq: (2) returns false for pattern type %zu, ",
                           s2->key.type);
                  fprintf (stdout, "nodes: ");
                  for (size_t i = 0; i < u_seg; i++)
                    fprintf (stdout, "n%zu, ", s2->key.segments[i]);
                  fprintf (stdout, "\n\t have length > %d in E1 = (", len);
                  ap_abstract0_fprint (stdout, pr->man_dcons, a1s->econs, NULL);
                  fprintf (stdout, "\n\t\t) and U1 = (top) while U2 = (");
                  ap_abstract0_fprint (stdout, pr->man_dcons, s2->dcons, NULL);
                  fprintf (stdout, "\n\t\t)\n");
                  fflush (stdout);
#endif
                  res = false;
                  goto ucons_is_leq_return;
                }
              // else, see saturation and the end
            }
          else
            { // good, the test succeded for this pattern, continue
              continue;
            }
        }
      // there is a problem but maybe saturation will solve it
      if (dosat)
        {
          // SATURATE a1s for nodes involved in the pattern
          ucons_t* aux = ucons_saturation (pr, a1s, s2->key.segments, u_seg);
          ucons_free_internal (pr, a1s);
          a1s = aux;
          HASH_FIND (hh, a1s->udcons, &s2->key.type, keylen, s1);
          if (s1 != NULL &&
              ap_abstract0_is_leq (pr->man_dcons, s1->dcons, s2->dcons))
            {
              // saturation solved the problem, continue with another pattern
              continue;
            }
          else
            {
#ifndef NDEBUG2
              fprintf (stdout, "\n@@@@ ucons_is_leq: (3) returns false (after saturation) for pattern type %zu, ",
                       s2->key.type);
              fprintf (stdout, "nodes: ");
              for (size_t i = 0; i < u_seg; i++)
                fprintf (stdout, "n%zu, ", s2->key.segments[i]);
              fprintf (stdout, "\n\t for U1=(");
              if (s1)
                ap_abstract0_fprint (stdout, pr->man_dcons, s1->dcons, NULL);
              else
                fprintf (stdout, "top");
              fprintf (stdout, "\n\t\t) not <= than  U2 = (");
              ap_abstract0_fprint (stdout, pr->man_dcons, s2->dcons, NULL);
              fprintf (stdout, "\n\t\t)\n");
              fflush (stdout);
#endif
              res = false;
              goto ucons_is_leq_return;
            }
        } // saturation done
    }
ucons_is_leq_return:
  // free allocated memory
  ucons_free_internal (pr, a1s);
  if (nsat) free (nsat);

#ifndef NDEBUG2
  fprintf (stdout, "\n@@@@ ucons_is_leq: returns %zu\n", res);
  fflush (stdout);
#endif
  return res;

}

/**
 * Test equality of abstract values.
 * Based on ucons_is_leq.
 */
bool
ucons_is_eq (ap_manager_t * man, ucons_t * a1, ucons_t * a2)
{
  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_IS_EQ, 0);

  if ((ucons_is_bottom (man, a1) && ucons_is_bottom (man, a2)) ||
      (ucons_is_top (man, a1) && ucons_is_top (man, a2)))
    return true;
  if (a1->datadim != a2->datadim || a1->segmentdim != a2->segmentdim)
    return false;

  if (!ap_abstract0_is_eq (pr->man_dcons, a1->econs, a2->econs))
    return false;

  return ucons_is_leq (man, a1, a2) && ucons_is_leq (man, a2, a1);
}

/*
 * TODO: Possibly needed for checking linear constraints representing - aliasing
 * between pointer variables, e.g., x = y, - or constraints between the
 * program variables (lengths and data) in the assert statements.
 */

/* NOT IMPLEMENTED */
bool
ucons_sat_lincons (ap_manager_t * man, ucons_t * a, ap_lincons0_t * lincons)
{

  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_SAT_LINCONS, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return true;
}

/* NOT IMPLEMENTED */
bool
ucons_sat_tcons (ap_manager_t * man, ucons_t * a, ap_tcons0_t * cons)
{

  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_SAT_TCONS, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return true;
}


/* TODO: Interval constraints are only between non-pointer variables */

/* NOT IMPLEMENTED */
bool
ucons_sat_interval (ap_manager_t * man, ucons_t * a,
                    ap_dim_t dim, ap_interval_t * i)
{

  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_SAT_INTERVAL, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return true;
}

/* NOT IMPLEMENTED */
bool
ucons_is_dimension_unconstrained (ap_manager_t * man, ucons_t * a,
                                  ap_dim_t dim)
{

  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_IS_DIMENSION_UNCONSTRAINED, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return false;
}

/* ============================================================ */
/* Extraction of properties */
/* ============================================================ */

/* NOT IMPLEMENTED */
ap_interval_t *
ucons_bound_linexpr (ap_manager_t * man, ucons_t * a, ap_linexpr0_t * expr)
{

  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_BOUND_LINEXPR, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return NULL;
}

/* NOT IMPLEMENTED */
ap_interval_t *
ucons_bound_texpr (ap_manager_t * man, ucons_t * a, ap_texpr0_t * expr)
{

  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_BOUND_TEXPR, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return NULL;
}

/* NOT IMPLEMENTED */
ap_interval_t *
ucons_bound_dimension (ap_manager_t * man, ucons_t * a, ap_dim_t dim)
{

  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_BOUND_DIMENSION, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return NULL;
}

/* NOT IMPLEMENTED */
ap_lincons0_array_t
ucons_to_lincons_array (ap_manager_t * man, ucons_t * a)
{

  ap_lincons0_array_t ar;
  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_TO_LINCONS_ARRAY, 0);
  ar = ap_lincons0_array_make (1);
  ar.p[0] = ap_lincons0_make_unsat ();
  return ar;
}

/* NOT IMPLEMENTED */
ap_tcons0_array_t
ucons_to_tcons_array (ap_manager_t * man, ucons_t * a)
{

  return ap_generic_to_tcons_array (man, a);
}

/* NOT IMPLEMENTED */
ap_interval_t **
ucons_to_box (ap_manager_t * man, ucons_t * a)
{

  ucons_internal_t *pr = ucons_init_from_manager (man, AP_FUNID_TO_BOX, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return NULL;
}

/* NOT IMPLEMENTED */
ap_generator0_array_t
ucons_to_generator_array (ap_manager_t * man, ucons_t * a)
{
  ucons_internal_t *pr =
          ucons_init_from_manager (man, AP_FUNID_TO_GENERATOR_ARRAY, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return ap_generator0_array_make (0);
}
