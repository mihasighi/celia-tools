/**************************************************************************/
/*                                                                        */
/*  CELIA Tools / MSET Abstract Domain                                    */
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


#include "mset.h"
#include "mset_fun.h"
#include "mset_internal.h"
#include "apron2shape.h"
#include "apron2smtlib.h"
#include "apron2acsl.h"


/* ============================================================ */
/* Printing */

/* ============================================================ */

void
mset_fprint_acsl (FILE * stream, mset_internal_t * pr, mset_t * a,
                  char **name_of_dim);
void
mset_fprint_dot (FILE * stream, mset_internal_t * pr, mset_t * a,
                 char **name_of_dim);
void
mset_fprint_smtlib (FILE * stream, mset_internal_t * pr, mset_t * a,
                    char **name_of_dim);

void
mset_fprint (FILE * stream, ap_manager_t * man, mset_t * a,
             char **name_of_dim)
{
  mset_internal_t *pr = mset_init_from_manager (man, AP_FUNID_FPRINT, 0);
  if (sh_print_is_dot ())
    mset_fprint_dot (stream, pr, a, name_of_dim);
  else if (sh_print_is_smtlib ())
    mset_fprint_smtlib (stream, pr, a, name_of_dim);
  else if (sh_print_is_acsl ())
    mset_fprint_acsl (stream, pr, a, name_of_dim);
  else
    assert (0);
}

void
mset_fprintdiff (FILE * stream, ap_manager_t * man,
                 mset_t * a1, mset_t * a2, char **name_of_dim)
{
  mset_internal_t *pr = mset_init_from_manager (man, AP_FUNID_FPRINTDIFF, 0);
  if (name_of_dim == NULL)
    shape_init_name_of_dim (a1->datadim, a1->segmdim);
  fprintf (stream, "mset1 (size %zu) and mset2 (size %zu)\n",
           a1->segmdim, a2->segmdim);
}

void
mset_fdump (FILE * stream, ap_manager_t * man, mset_t * a)
{
  mset_fprint (stream, man, a, NULL);
}

/* ============================================================ */
/* Printing internals */
/* ============================================================ */

/**
 * @brief Initialize names of data for this domain.
 * 
 * Allocate and initialize an array of strings of size
 * corresponding to the data expressed in @p a.
 *
 * @param[in] a  The abstract value
 * @result       An array of strings
 */
char **
mset_init_name_data (mset_t * a, char **name_of_dim)
{
  char **dname;
  size_t i, size;
  size = (a->datadim + 2 * a->segmdim);
  dname = (char **) malloc (size * sizeof (char *));
  for (i = 0; i < a->datadim; i++)
    dname[i] = (name_of_dim) ? name_of_dim[i] : shape_name_of_dim (i);
  for (i = 0; i < a->segmdim; i++)
    {
      char *name = shape_name_of_dim (a->datadim + i);
      size_t lsize = (7 + strlen (name));
      dname[a->datadim + i] = (char *) malloc (lsize * sizeof (char));
      snprintf (dname[a->datadim + i], lsize, "%s[0]", name);
      dname[a->datadim + a->segmdim + i] =
        (char *) malloc (lsize * sizeof (char));
      snprintf (dname[a->datadim + a->segmdim + i], lsize, "len(%s)", name);
    }

  return dname;
}


/**
 * @brief Initialize names of mset constraints for this domain.
 * 
 * Allocate and initialize an array of strings of size
 * corresponding to the mset variables expressed in @p a.
 *
 * @param[in] a  The abstract value
 * @result       An array of strings
 */
char **
mset_init_name_mset (mset_t * a, char **name_of_dim)
{
  char **dname;
  size_t i, size;
  size = a->datadim + 2 * a->segmdim;
  dname = (char **) malloc (size * sizeof (char *));
  for (i = 0; i < a->datadim; i++)
    dname[i] = (name_of_dim) ? name_of_dim[i] : shape_name_of_dim (i);
  for (i = 0; i < a->segmdim; i++)
    {
      char *name = shape_name_of_dim (a->datadim + i);
      size_t lsize = (7 + strlen (name));
      dname[a->datadim + i] = (char *) malloc (lsize * sizeof (char));
      snprintf (dname[a->datadim + i], lsize, "mshd(%s)", name);
      dname[a->datadim + a->segmdim + i] =
        (char *) malloc (lsize * sizeof (char));
      snprintf (dname[a->datadim + a->segmdim + i], lsize, "mstl(%s)", name);
    }
  return dname;
}

/**
 * @brief Free memory allocated in @p t.
 * 
 * @param[in] a  The abstract value gives the size of @p dname
 * @param[in] t  An array of strings to be deallocated
 */
void
mset_free_name (mset_t * a, char **t)
{
  if (t)
    {
      size_t i, size;
      size = a->datadim + 2 * a->segmdim;
      for (i = a->datadim; i < size; i++)
        if (t[i])
          free (t[i]);
      memset (t, 0, size * sizeof (char *));
      free (t);
    }
}


/**
 * @brief Print the abstract value using the ACSL format.
 * 
 * Use the ACSL format to print the abstract value in the
 * file @p stream.
 * 
 * @param[in] stream   Output stream
 * @param[in] pr       Internal manager
 * @param[in] a        Abstract value to be print
 * @param[in] name_of_dim Names for the array represented
 */
void
mset_fprint_acsl (FILE * stream, mset_internal_t * pr, mset_t * a,
                  char **name_of_dim)
{
  if (!a)
    fprintf (stream, "EMPTY");
  else
    {
      if (!name_of_dim)
        shape_init_name_of_dim (a->datadim, a->segmdim);
      char **dname = mset_init_name_data (a, name_of_dim);
      ap_lincons0_array_t arr =
        ap_abstract0_to_lincons_array (pr->man_dcons, a->dcons);
      ap_lincons0_array_fprint_acsl (stream, &arr, dname);
      ap_lincons0_array_clear (&arr);
      mset_free_name (a, dname);

      dname = mset_init_name_mset (a, name_of_dim);
      arr = ap_abstract0_to_lincons_array (pr->man_mscons, a->mscons);
      ap_lincons0_array_fprint_acsl (stream, &arr, dname);
      ap_lincons0_array_clear (&arr);
      mset_free_name (a, dname);
    }
}


/**
 * @brief Print the abstract value using the DOT format.
 * 
 * Use the DOT format to print the abstract value in the
 * file @p stream.
 *  
 * @param[in] stream   Output stream
 * @param[in] pr       Internal manager
 * @param[in] a        Abstract value to be print
 * @param[in] name_of_dim Names for the array represented
 */
void
mset_fprint_dot (FILE * stream, mset_internal_t * pr, mset_t * a,
                 char **name_of_dim)
{
  fprintf (stream,
           "\n\tsubgraph cluster_mset_%zu {\n\tnode [shape=record] ;\n",
           ushape_number);
  if (!a)
    fprintf (stream, "\tlabel=\"mset %zu EMPTY\" ;\n}\n", ushape_number);
  else
    {
      fprintf (stream,
               "\tlabel=\"mset %zu of (datadim=%zu, segmdim=%zu)\" ;\n",
               ushape_number, a->datadim, a->segmdim);
      if (!name_of_dim)
        shape_init_name_of_dim (a->datadim, a->segmdim);
      if (!a->dcons)
        fprintf (stream, "\tmset_dcons_%zu [label=\"bot\"] ;\n",
                 ushape_number);
      else
        {
          char **dname = mset_init_name_data (a, name_of_dim);
          fprintf (stream, "\tmset_dcons_%zu [label=<<table><tr><td>dcons: ",
                   ushape_number);
          ap_abstract0_fprint (stream, pr->man_dcons, a->dcons, dname);
          fprintf (stream, "</td></tr></table>> ] ;\n");
          mset_free_name (a, dname);
        }
      if (!a->mscons)
        fprintf (stream, "\tmset_mscons_%zu [label=\"bot\"] ;\n",
                 ushape_number);
      else
        {
          char **dname = mset_init_name_mset (a, name_of_dim);
          fprintf (stream,
                   "\tmset_mscons_%zu [label=<<table><tr><td>mscons: ",
                   ushape_number);
          ap_abstract0_fprint (stream, pr->man_mscons, a->mscons, dname);
          fprintf (stream, "</td></tr></table>> ] ;\n");
          mset_free_name (a, dname);
        }
      fprintf (stream, "\t}\n");
    }
}

void
mset_fprint_smtlib (FILE * stream, mset_internal_t * pr, mset_t * a,
                    char **name_of_dim)
{
  if (!a)
    {
      fprintf (stream, "\t\t;; mset constraints empty\n");
      return;
    }

  size_t i, size;
  size = (a->datadim + 2 * a->segmdim);
  // print data constraints
  // initialize data and len node names
  char **dname =
    ap_dimension_to_smt (size, a->datadim, a->segmdim, 0, name_of_dim);
  // print the list of constraints
  ap_lincons0_array_t arr =
    ap_abstract0_to_lincons_array (pr->man_dcons, a->dcons);
  ap_lincons0_array_fprint_smtlib (stream, &arr, dname);
  ap_lincons0_array_clear (&arr);

  // print the mset constraints
  // initialize mset node names
  for (i = 0; i < a->segmdim; i++)
    {
      char *name = shape_name_of_dim (a->datadim + i);
      size_t lsize = (8 + strlen (name));
      if (dname[a->datadim + a->segmdim + i])
        free (dname[a->datadim + a->segmdim + i]);
      dname[a->datadim + a->segmdim + i] =
        (char *) malloc (lsize * sizeof (char));
      snprintf (dname[a->datadim + a->segmdim + i], lsize, "(mset %s)", name);
    }

  // print the list of constraints
  arr = ap_abstract0_to_lincons_array (pr->man_mscons, a->mscons);
  ap_lincons0_array_fprint_smtlib (stream, &arr, dname);
  ap_lincons0_array_clear (&arr);

  mset_free_name (a, dname);
}

/* ============================================================ */
/* Serialization */
/* ============================================================ */

/* NOT IMPLEMENTED: do nothing */
ap_membuf_t
mset_serialize_raw (ap_manager_t * man, mset_t * a)
{
  mset_internal_t *pr =
    mset_init_from_manager (man, AP_FUNID_SERIALIZE_RAW, 0);
  ap_membuf_t buf;
  buf.size = 0;
  buf.ptr = (size_t *) a;       /* to remove warning */
  ap_manager_raise_exception (man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
                              "not implemented");
  return buf;
}

/**
 * @brief Build an abstract value using the code in @p ptr.
 * 
 * Read at @p ptr tree size_t numbers representing
 * the code (datasim, segmentdim, code) then build the
 * abstract value for the code.
 *
 * @param[in] man   Global manager
 * @param[in] ptr   Address storing the code
 * @param[in] size  Number of size_t to be read
 * @return          The abstract value represented by the code
 */
mset_t *
mset_deserialize_raw (ap_manager_t * man, void *ptr, size_t * size)
{
  mset_internal_t *pr =
    mset_init_from_manager (man, AP_FUNID_DESERIALIZE_RAW, 0);

  if (size != size)
    return NULL;                /* to remove warning on unused parameter */

  size_t *mset_raw = (size_t *) ptr;
  size_t datadim = mset_raw[0];
  size_t segmdim = mset_raw[1];
  size_t code = mset_raw[2];
  ap_dim_t l = 0;
  ap_dim_t S = 2;               // TODO: take care of collapse with the sum domain
  ap_dim_t T = 3;
  size_t nodex = 1;
  size_t nodey = 2;
  size_t nodez = 3;
  mset_t *r = mset_top (man, datadim, segmdim);

  /* common constraints l[x]==_l>=1 and M[x]==S */
  if (code != 1)
    {
      assert (segmdim >= 2);
      ap_lincons0_array_t arr = ap_lincons0_array_make (3);
      arr.p[0] =
        shape_lincons_x_y_v_cst (AP_CONS_EQ, OFFSET_LEN, 1,
                                 datadim + nodex, 0, 0, -1, l, 0, datadim,
                                 segmdim);
      arr.p[1] =
        shape_lincons_x_y_v_cst (AP_CONS_SUPEQ, OFFSET_LEN, 1,
                                 datadim + nodex, 0, 0, 0, 0, -1, datadim,
                                 segmdim);
      arr.p[2] =
        shape_lincons_x_y_v_cst (AP_CONS_EQ, OFFSET_MSET, 1, datadim + nodex,
                                 0, 0, -1, S, 0, datadim, segmdim);
      r = mset_meet_lincons_array (man, true, r, &arr);
      ap_lincons0_array_clear (&arr);
    }
  /* common constraint l[x]==l[y] and l[y]>=1 */
  if (code == 2 || code == 4)
    {
      assert (segmdim >= 3);
      ap_lincons0_array_t arr = ap_lincons0_array_make (2);
      arr.p[0] =
        shape_lincons_x_y_v_cst (AP_CONS_EQ, OFFSET_LEN, 1,
                                 datadim + nodex, -1, datadim + nodey, 0, 0,
                                 0, datadim, segmdim);
      arr.p[1] =
        shape_lincons_x_y_v_cst (AP_CONS_SUPEQ, OFFSET_LEN, 1,
                                 datadim + nodey, 0, 0, 0, 0, -1, datadim,
                                 segmdim);
      r = mset_meet_lincons_array (man, true, r, &arr);
      ap_lincons0_array_clear (&arr);
    }
  /* specific constraints */
  switch (code)
    {
    case 1:                    /* l[x]+l[y]==_l and l[x]>=1 and l[y]>=1 and M[x]+M[y]==S */
      {
        assert (segmdim >= 2);
        ap_lincons0_array_t arr = ap_lincons0_array_make (4);
        arr.p[0] =
          shape_lincons_x_y_v_cst (AP_CONS_EQ, OFFSET_LEN, 1,
                                   datadim + nodex, 1, datadim + nodey, -1, l,
                                   0, datadim, segmdim);
        arr.p[1] =
          shape_lincons_x_y_v_cst (AP_CONS_SUPEQ, OFFSET_LEN, 1,
                                   datadim + nodex, 0, 0, 0, 0, -1, datadim,
                                   segmdim);
        arr.p[2] =
          shape_lincons_x_y_v_cst (AP_CONS_SUPEQ, OFFSET_LEN, 1,
                                   datadim + nodey, 0, 0, 0, 0, -1, datadim,
                                   segmdim);
        arr.p[3] =
          shape_lincons_x_y_v_cst (AP_CONS_EQ, OFFSET_MSET, 1,
                                   datadim + nodex, 1, datadim + nodey, -1, S,
                                   0, datadim, segmdim);
        r = mset_meet_lincons_array (man, true, r, &arr);
        ap_lincons0_array_clear (&arr);
        break;
      }
    case 3:                    /* non-equal length lists: l[y]+1<=_l and l[y]>=1 */
      {
        assert (segmdim >= 3);
        ap_lincons0_array_t arr = ap_lincons0_array_make (2);
        arr.p[0] =
          shape_lincons_x_y_v_cst (AP_CONS_SUPEQ, OFFSET_LEN, -1,
                                   datadim + nodey, 0, 0, 1, l, -1, datadim,
                                   segmdim);
        arr.p[1] =
          shape_lincons_x_y_v_cst (AP_CONS_SUPEQ, OFFSET_LEN, 1,
                                   datadim + nodey, 0, 0, 0, 0, -1, datadim,
                                   segmdim);
        r = mset_meet_lincons_array (man, true, r, &arr);
        ap_lincons0_array_clear (&arr);
        break;
      }
    case 4:                    /* equal length lists: l[z]=l[x] and M[y]=T */
      {
        assert (segmdim >= 4 && datadim >= 4);
        ap_lincons0_array_t arr = ap_lincons0_array_make (2);
        arr.p[0] =
          shape_lincons_x_y_v_cst (AP_CONS_EQ, OFFSET_LEN, -1,
                                   datadim + nodex, 1, datadim + nodez, 0, 0,
                                   0, datadim, segmdim);
        arr.p[1] =
          shape_lincons_x_y_v_cst (AP_CONS_EQ, OFFSET_MSET, 1,
                                   datadim + nodey, 0, 0, -1, T, 0, datadim,
                                   segmdim);
        r = mset_meet_lincons_array (man, true, r, &arr);
        ap_lincons0_array_clear (&arr);
        break;
      }
    default:
      break;
    }
  return r;
}
