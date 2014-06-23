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


#include "shape_macros.h"
#include "shape.h"
#include "shape_internal.h"


/* ============================================================ */
/* Printing */

/* ============================================================ */

FILE * shape_fdump_new_file (FILE * stream, shape_internal_t * pr, bool isdot);
void shape_fprint_dot (FILE * shape_stream, shape_internal_t *pr,
                       shape_t * a, char **name_of_dim);
void shape_fprint_smt (FILE* f, shape_internal_t *pr,
                       shape_t* a, char** name_of_dim, bool ispos);

void
shape_fprint (FILE * stream, ap_manager_t * man,
              shape_t * a, char **name_of_dim)
{
  shape_internal_t *pr = shape_init_from_manager (man, AP_FUNID_FPRINT, 0);
  bool isdot = shape_get_print ();
  FILE *shape_stream = shape_fdump_new_file (stream, pr, isdot);
  if (isdot)
    shape_fprint_dot (shape_stream, pr, a, name_of_dim);
  else
    shape_fprint_smt (shape_stream, pr, a, name_of_dim, true);
  fclose (shape_stream);
}

void
shape_fprintdiff (FILE * stream, ap_manager_t * man,
                  shape_t * a1, shape_t * a2, char **name_of_dim)
{
  shape_internal_t *pr =
          shape_init_from_manager (man, AP_FUNID_FPRINTDIFF, 0);
  fprintf (stream, "shape1 (size %zu) and shape2 (%zu)\n",
           a1->msize, a2->msize);
  shape_fprint (stream, man, a1, name_of_dim);
  fprintf (stream, "\n<------>\n");
  shape_fprint (stream, man, a2, name_of_dim);
  fflush (stream);
}

void
shape_fdump (FILE * stream, ap_manager_t * man, shape_t * a)
{
  shape_fprint (stream, man, a, NULL);
}

/* ============================================================ */
/* Printing internals */

/* ============================================================ */

FILE *
shape_fdump_new_file (FILE * stream, shape_internal_t * pr, bool isdot)
{
  FILE *r = NULL;
  char filename[24];
  memset (filename, 0, 24 * sizeof (char));
  snprintf (filename, 24, "pan/%s_%ld.%s",
            (isdot) ? "f" : "sl3",
            pr->filenum++,
            (isdot) ? "shp" : "smt");
  r = fopen (filename, "w");
  if (!r)
    {
      ERROR ("Unable to open pan file! quit.", return NULL;);
    }
#ifndef NDEBUG1
  fprintf (stream, "<file:%s>\n", filename);
#endif
  return r;
}

void
shape_fprint_dot (FILE * shape_stream, shape_internal_t *pr, shape_t * a,
                  char **name_of_dim)
{
  size_t i;
  fprintf (shape_stream, "digraph shape {\n");
  fprintf (shape_stream, "\tgraph [rankdir=\"LR\"];\n");
  if (!a)
    {
      fprintf (shape_stream, "label=\"shape EMPTY\";\n}\n");
      return;
    }
  fprintf (shape_stream, "\tlabel=\"shape size %zu ",
           a->msize);
  if (name_of_dim)
    {
      fprintf (shape_stream, " over variables rev[");
      for (i = a->intdim + a->realdim; i > 0; i--)
        fprintf (shape_stream, "%s ", name_of_dim[i]);
      fprintf (shape_stream, "]");
    }
  fprintf (shape_stream, "\" ;\n");
  for (i = 0; i < a->msize; i++)
    {
      fprintf (shape_stream, "\nsubgraph cluster_%zd {\n", i);
      ushape_number = i;
      ushape_fprint (shape_stream, pr->man, a->m.p[i], name_of_dim);
      fprintf (shape_stream, "}\n");
    }
  ushape_number = 0;
  fprintf (shape_stream, "}\n");
  fflush (shape_stream);
}

void
shape_fprint_smt_prelude (FILE* f, shape_internal_t *pr,
                          shape_t* a, char** name_of_dim)
{
  // logic
  fprintf (f, "(set-logic SL3)\n");
  // sorts
  fprintf (f, "(declare-sort Ptr 0)\n");
  fprintf (f, "(declare-sort Node 0)\n");
  // logic functions
  fprintf (f, "(declare-fun nilNode () Node)\n");
  fprintf (f, "(declare-fun len (Node) Int)\n");
  fprintf (f, "(declare-fun data (Node Int) Int)\n");
  fprintf (f, "(declare-fun sum (Node) Int)\n");
  fprintf (f, "(declare-fun mset (Node) Int)\n");
  fprintf (f, "(declare-fun sep (Bool Bool) Bool)\n");
  fprintf (f, "(declare-fun ls (Node Node) Bool)\n");
  fprintf (f, "(declare-fun label (Ptr Node) Bool)\n");
  fprintf (f, "(declare-fun Gall (Int Node) Bool)\n");
  fprintf (f, "(declare-fun Gle2 (Int Int Node) Bool)\n");
  fprintf (f, "(declare-fun Gsucc2 (Int Int Node) Bool)\n");
  fprintf (f, "(declare-fun Gfst (Int Node) Bool)\n");
  fprintf (f, "(declare-fun Glst (Int Node) Bool)\n");
  fprintf (f, "(declare-fun Geq2 (Int Node Int Node) Bool)\n");
  fprintf (f, "(declare-fun Gsl2 (Int Node Int Node Int) Bool)\n");
  fprintf (f, "(declare-fun Gsr2 (Int Node Int Node Int) Bool)\n");
  // variables
  size_t i;
  for (i = 0; i < (a->intdim + a->realdim); i++)
    {
      fprintf (f, "(declare-fun ");
      if (name_of_dim)
        fprintf (f, "%s", name_of_dim[i]);
      else
        fprintf (f, "x%zu", i);
      fprintf (f, " () ");
      if (i < a->intdim)
        fprintf (f, "Int");
      else
        fprintf (f, "Ptr");
      fprintf (f, ")\n");
    }
}

/**
 * Print the abstract value in SMTLIB2 format.
 */
void
shape_fprint_smt (FILE* f, shape_internal_t *pr,
                  shape_t* a, char** name_of_dim, bool ispos)
{
  if (ispos)
    shape_fprint_smt_prelude (f, pr, a, name_of_dim);
  fprintf (f, "(assert \n");

  if (!ispos)
    fprintf (f, "(not \n");

  // print all ushapes
  if (a->msize > 1)
    fprintf (f, "(or \n");
  size_t i;
  for (i = 0; i < a->msize; i++)
    ushape_fprint (f, pr->man, a->m.p[i], name_of_dim);
  if (a->msize > 1) // end or
    fprintf (f, ") ;; end or\n");

  if (!ispos) // end not
    fprintf (f, ") ;; end not\n");

  // end assert
  fprintf (f, ") ;; end assert\n");
  fflush (f);

}

void
shape_fdump_le (FILE * stream, shape_internal_t * pr,
                shape_t * a1, shape_t * a2)
{
  bool isdot = shape_get_print ();
  if (isdot)
    {
      shape_fdump (stream, pr->man, a1);
      fprintf (stream, " ==> ");
      shape_fdump (stream, pr->man, a2);
      fflush (stream);
    }
  else
    {
      FILE *shape_stream = shape_fdump_new_file (stream, pr, isdot);
      shape_fprint_smt (shape_stream, pr, a1, NULL, true);
      shape_fprint_smt (shape_stream, pr, a2, NULL, false);
      fprintf (shape_stream, "(check-sat)\n");
      fflush (shape_stream);
    }
}

/* ============================================================ */
/* Serialization */
/* ============================================================ */

/* TODO: priority 0 */

/* NOT IMPLEMENTED: do nothing */
ap_membuf_t
shape_serialize_raw (ap_manager_t * man, shape_t * a)
{
  shape_internal_t *pr =
          shape_init_from_manager (man, AP_FUNID_SERIALIZE_RAW, 0);
  ap_membuf_t buf;
  buf.size = 0;
  buf.ptr = NULL;
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return buf;
}

/* TODO: priority 0 */

/* NOT IMPLEMENTED: do nothing */
shape_t *
shape_deserialize_raw (ap_manager_t * man, void *ptr, size_t * size)
{
  shape_internal_t *pr =
          shape_init_from_manager (man, AP_FUNID_DESERIALIZE_RAW, 0);
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return NULL;
}
