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
#include "ushape.h"
#include "ushape_internal.h"


/* ============================================================ */
/* Printing */

/* ============================================================ */

void ushape_fprint_dot (FILE * stream, shape_internal_t * man,
                        ushape_t * a, char **name_of_dim);
void ushape_fprint_smt (FILE * stream, shape_internal_t * man,
                        ushape_t * a, char **name_of_dim);

void
ushape_fprint (FILE * stream, ap_manager_t * man,
               ushape_t * a, char **name_of_dim)
{
  ushape_internal_t *pr = ushape_init_from_manager (man, AP_FUNID_FPRINT, 0);
  bool isdot = shape_get_print ();
  if (isdot)
    ushape_fprint_dot (stream, pr, a, name_of_dim);
  else
    ushape_fprint_smt (stream, pr, a, name_of_dim);
}

void
ushape_fprintdiff (FILE * stream, ap_manager_t * man,
                   ushape_t * a1, ushape_t * a2, char **name_of_dim)
{
  ushape_internal_t *pr =
          ushape_init_from_manager (man, AP_FUNID_FPRINTDIFF, 0);
  size_t i;
  hgraph_fprintdiff (stream, man, (a1) ? a1->h : NULL, (a2) ? a2->h : NULL,
                     name_of_dim);
  for (i = 0; i < pr->size_scons; i++)
    {
      fprintf (stream, "\t[univ constraint:\n");
      ap_abstract0_fprintdiff (stream, pr->man_scons[i],
                               (a1 && a1->scons) ? a1->scons[i] : NULL, (a2
                                                                         &&
                                                                         a2->
                                                                         scons)
                               ? a2->scons[i] : NULL, name_of_dim);
      fprintf (stream, "\n\t]\n");
    }
}

void
ushape_fdump (FILE * stream, ap_manager_t * man, ushape_t * a)
{
  ushape_fprint (stream, man, a, NULL);
}

void
ushape_fprint_dot (FILE * stream, shape_internal_t *pr,
                   ushape_t * a, char **name_of_dim)
{
  size_t i;
  fprintf (stream, "\tcolor = black ;\n");
  if (!a)
    {
      fprintf (stream, "\tlabel = \"ushape empty\" ;\n");
      return;
    }
  fprintf (stream, "\tlabel = \"ushape of dim (%zu,%zu)\" ;\n", a->datadim, a->ptrdim);
  fflush (stream);
  fprintf (stream, "\tsubgraph cluster_hgraph_%zu {\n ", ushape_number);
  if (a->h)
    hgraph_fprint (stream, pr->man, a->h, name_of_dim);
  fprintf (stream, "\n\t}\n");
  if (a->scons)
    {
      for (i = 0; i < pr->size_scons; i++)
        if (a->scons[i])
          ap_abstract0_fprint (stream, pr->man_scons[i], a->scons[i],
                               name_of_dim);
        else
          fprintf (stream, "\tscons_%zu_%zu [label=\"[NULL]\"] ;\n",
                   ushape_number, i);
    }
  else
    fprintf (stream, "\tsubgraph cluster_scons_%zu { label = \"scons [NULL]\" ; }\n",
             ushape_number);
}

void
ushape_fprint_smt (FILE* stream, hgraph_internal_t* pr,
                   ushape_t* a, char** name_of_dim)
{
  size_t i;

  if (!a)
    {
      fprintf (stream, "(true)\n");
      return;
    }

  if (a->h)
    {
      // print the node variables used
      fprintf (stream, "\t(exists (");
      for (i = 1; i < a->h->size; i++)
        fprintf (stream, " (?n%zu Node) ", i);
      fprintf (stream, ")\n\t(and ");
      // print the graph
      hgraph_fprint (stream, pr->man, a->h, name_of_dim);
    }
  else
    // only data constraints
    fprintf (stream, "\t(and ");
  if (a->scons)
    {
      for (i = 0; i < pr->size_scons; i++)
        if (a->scons[i])
          ap_abstract0_fprint (stream, pr->man_scons[i], a->scons[i],
                               name_of_dim);
        else
          fprintf (stream, "\t(true)\n");
    }
  else
    fprintf (stream, "\t(true)\n");

  // end and
  fprintf (stream, "\n\t) ;; end and\n");
  if (a->h)
    fprintf (stream, "\n\t) ;; end exists\n"); // exists

}


/* ============================================================ */
/* Serialization */
/* ============================================================ */

/* NOT IMPLEMENTED: do nothing */
ap_membuf_t
ushape_serialize_raw (ap_manager_t * man, ushape_t * a)
{
  ushape_internal_t *pr =
          ushape_init_from_manager (man, AP_FUNID_SERIALIZE_RAW, 0);
  ap_membuf_t buf;
  buf.size = 0;
  buf.ptr = NULL;
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return buf;
}

/* NOT IMPLEMENTED: do nothing */
ushape_t *
ushape_deserialize_raw (ap_manager_t * man, void *ptr, size_t * size)
{
  ushape_internal_t *pr =
          ushape_init_from_manager (man, AP_FUNID_DESERIALIZE_RAW, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return NULL;
}
