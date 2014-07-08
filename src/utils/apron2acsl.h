/**************************************************************************/
/*                                                                        */
/*  CELIA Tools / Utilities for Abstract Domains                          */
/*                                                                        */
/*  Copyright (C) 2014                                                    */
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


#ifndef __APRON2ACSL_H_
#define __APRON2ACSL_H_

#include "apron2shape.h"
#include "ap_linexpr0.h"
#include "ap_lincons0.h"
#include "ap_texpr0.h"
#include "ap_tcons0.h"

/* *INDENT-OFF* */
#ifdef __cplusplus
extern "C"
{
#endif
  /* *INDENT-ON* */

  /* ================================================================== */
  /* Expressions and constraints */
  /* ================================================================== */

void ap_linexpr0_fprint_acsl (FILE * stream, ap_linexpr0_t * expr,
                              char **name_of_dim);
void ap_texpr0_fprint_acsl (FILE * stream, ap_texpr0_t * expr,
                            char **name_of_dim);
  /* Print an APRON expression with the ACSL syntax */

void ap_lincons0_fprint_acsl (FILE * stream, ap_lincons0_t * cons,
                              char **name_of_dim);
void ap_lincons0_array_fprint_acsl (FILE * stream, ap_lincons0_array_t * arr,
                                    char **name_of_dim);
void ap_tcons0_fprint_acsl (FILE * stream, ap_tcons0_t * cons,
                            char **name_of_dim);
void ap_tcons0_array_fprint_acsl (FILE * stream, ap_tcons0_array_t * arr,
                                  char **name_of_dim);
  /* Print an APRON constraint with the ACSL syntax */

char **ap_dimension_to_acsl (size_t size, size_t datadim, size_t segmdim,
                             size_t offset, char **name_of_dim);
  /* Initializes the names of dimensions for printing in the ACSL format */

  /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */

#endif /* __APRON2ACSL_H_ */
