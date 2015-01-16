/* SHAD - Library of shape abstract domains
 * Copyright (C) 2012-2013 LIAFA
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * Please read the COPYING file packaged in the distribution.
 */

/* 
 * sh_manager.h: global manager passed to all functions
 */

#ifndef _SH_MANAGER_H_
#define _SH_MANAGER_H_

//#include "ap_coeff.h"
#include "sh_loc.h"
#include "sh_typ.h"
#include "sh_env.h"

#ifdef __cplusplus
extern "C" {
#endif


/* ********************************************************************** */
/* I. Types */
/* ********************************************************************** */

/* ====================================================================== */
/* I.1 Identifying functions */
/* ====================================================================== */

typedef enum sh_funid_t {
  SH_FUNID_UNKNOWN,
  SH_FUNID_COPY,
  SH_FUNID_FREE,
  SH_FUNID_ASIZE, /* For avoiding name conflict with SH_FUNID_SIZE */
  SH_FUNID_MINIMIZE,
  SH_FUNID_CANONICALIZE,
  SH_FUNID_HASH,
  SH_FUNID_APPROXIMATE,
  SH_FUNID_FPRINT,
  SH_FUNID_FPRINTDIFF,
  SH_FUNID_FDUMP,
  SH_FUNID_SERIALIZE_RAW,
  SH_FUNID_DESERIALIZE_RAW,
  SH_FUNID_BOTTOM,
  SH_FUNID_TOP,
  SH_FUNID_OF_BOX,
  SH_FUNID_DIMENSION,
  SH_FUNID_IS_BOTTOM,
  SH_FUNID_IS_TOP,
  SH_FUNID_IS_LEQ,
  SH_FUNID_IS_EQ,
  SH_FUNID_IS_DIMENSION_UNCONSTRAINED,
  SH_FUNID_SAT_TCONS,
  SH_FUNID_SAT_TFORM,
  SH_FUNID_TO_TCONS_ARRAY,
  SH_FUNID_TO_TFORM_ARRAY,
  SH_FUNID_MEET,
  SH_FUNID_MEET_ARRAY,
  SH_FUNID_MEET_TCONS_ARRAY,
  SH_FUNID_JOIN,
  SH_FUNID_JOIN_ARRAY,
  SH_FUNID_ASSIGN_TEXPR_ARRAY,
  SH_FUNID_SUBSTITUTE_TEXPR_ARRAY,
  SH_FUNID_ADD_DIMENSIONS,
  SH_FUNID_REMOVE_DIMENSIONS,
  SH_FUNID_PERMUTE_DIMENSIONS,
  SH_FUNID_FORGET_ARRAY,
  SH_FUNID_EXPAND,
  SH_FUNID_FOLD,
  SH_FUNID_WIDENING,
  SH_FUNID_CLOSURE,
  SH_FUNID_SIZE,
  SH_FUNID_CHANGE_ENVIRONMENT,
  SH_FUNID_RENAME_ARRAY,
  SH_FUNID_SIZE2
} sh_funid_t;

extern const char* sh_name_of_funid[SH_FUNID_SIZE2];
/* give the name of a function identifier */


/* ====================================================================== */
/* I.2 Exceptions */
/* ====================================================================== */

/* Exceptions (public type) */
typedef enum sh_exc_t {
  SH_EXC_NONE,             /* no exception detected */
  SH_EXC_TIMEOUT,          /* timeout detected */
  SH_EXC_OUT_OF_SPACE,     /* out of space detected */
  SH_EXC_OVERFLOW,         /* magnitude overflow detected */
  SH_EXC_INVALID_ARGUMENT, /* invalid arguments */
  SH_EXC_NOT_IMPLEMENTED,  /* not implemented */
  SH_EXC_SIZE
} sh_exc_t;

extern const char* sh_name_of_exception[SH_EXC_SIZE];

/* Exception log */
typedef struct sh_exclog_t {
  sh_exc_t exn;
  sh_funid_t funid;
  char* msg;                   /* dynamically allocated */
  struct sh_exclog_t* tail;
} sh_exclog_t;

/* Exceptions and other indications (out) (opaque type) */
typedef struct sh_result_t {
  sh_exclog_t* exclog; /* history of exceptions */
  sh_exc_t exn;        /* exception for the last called function */
  gboolean flag_exact;  /* result is mathematically exact or don't know */
  gboolean flag_best;   /* result is best correct approximation or don't know */
} sh_result_t;


/* ====================================================================== */
/* I.2 Options */
/* ====================================================================== */

/* Option associated to each function (public type) */
typedef struct sh_funopt_t {
  int algorithm;
  /* Algorithm selection:
     - 0 is default algorithm;
     - MAX_INT is most accurate available;
     - MIN_INT is most efficient available;
     - otherwise, no accuracy or speed meaning
  */
  size_t timeout; /* unit !? */
  /* Above the given computation time, the function may abort with the
     exception flag flag_time_out on.
  */
  size_t max_object_size; /* in abstract object size unit. */
  /* If during the computation, the size of some object reach this limit, the
     function may abort with the exception flag flag_out_of_space on.
  */
  gboolean flag_exact_wanted;
  /* return information about exactitude if possible
  */
  gboolean flag_best_wanted;
  /* return information about best correct approximation if possible
  */
} sh_funopt_t;

/* Options (in) (opaque type) */
typedef struct sh_option_t {
  sh_funopt_t funopt[SH_FUNID_SIZE];
  gboolean abort_if_exception[SH_EXC_SIZE];
  // ap_scalar_discr_t scalar_discr; /* Preferred type for scalars */
} sh_option_t;


/* ====================================================================== */
/* I.3 Manager */
/* ====================================================================== */

/* Manager (opaque type) */
struct sh_manager_s {
  char* library;                 /* name of the effective library */
  char* version;                 /* version of the effective library */
  void* internal;                /* library dependent,
				    should be different for each thread
				    (working space) */
  void* funptr[SH_FUNID_SIZE];   /* Array of function pointers,
				    initialized by the effective library */
  sh_option_t option;            /* Options (in) */
  sh_result_t result;            /* Exceptions and other indications (out) */
  void (*internal_free)(void*);  /* deallocation function for internal */
  size_t count;                  /* reference counter */

  /* Informations about the analyzed program. */
  sh_filenv_t* filenv;            /* Files managed */
  sh_typenv_t* typenv;            /* Types managed and their specification */
  sh_stack_t*  varenv;            /* Sets of variables by frame */

};

/* ********************************************************************** */
/* II. User Functions */
/* ********************************************************************** */

/*
 * Constructor and destructor for result.
 */
  sh_exclog_t* sh_exc_cons(sh_exc_t exn,
			   sh_funid_t funid, const char* msg,
			   sh_exclog_t* tail);
  void sh_exclog_free(sh_exclog_t* head);
  
  void sh_manager_clear_exclog(sh_manager_t* man);
  /* erase the current log of exception */


  /* 
   * Constructor and destructor for manager 
   */
  sh_manager_t* sh_manager_alloc(char* library, char* version,
				 void* internal,
				 void (*internal_free)(void*));
  void sh_manager_free(sh_manager_t* man);
  /* dereference the counter,
     and possibly free internal field if it is not yet put to NULL */

  /*
   * Other manager functions
   */
  const char* sh_manager_get_library(sh_manager_t* man);
  const char* sh_manager_get_version(sh_manager_t* man);
  sh_filenv_t* sh_manager_get_filenv(sh_manager_t* man);
  sh_typenv_t* sh_manager_get_typenv(sh_manager_t* man);
  sh_stack_t* sh_manager_get_varenv(sh_manager_t* man);
  /* Reading fields */
  void sh_manager_fdump(FILE* stream, sh_manager_t* man);
  /* Dump informations */

  int sh_manager_get_frame_of_loc(sh_manager_t* man, sh_loc_t* loc);
  /* Computes the frame of the location */

  sh_funopt_t sh_manager_get_funopt(sh_manager_t* man, sh_funid_t funid);
  gboolean sh_manager_get_abort_if_exception(sh_manager_t* man, sh_exc_t exn);

  sh_exc_t sh_manager_get_exception(sh_manager_t* man);
  /* Get the last exception raised */
  sh_exclog_t* sh_manager_get_exclog(sh_manager_t* man);
  /* Get the full log of exception */
  gboolean sh_manager_get_flag_exact(sh_manager_t* man);
  gboolean sh_manager_get_flag_best(sh_manager_t* man);

  /* Settings fields */
  void sh_funopt_init(sh_funopt_t* fopt);
  void sh_manager_set_funopt(sh_manager_t* man, sh_funid_t funid, sh_funopt_t* funopt);
  void sh_manager_set_abort_if_exception(sh_manager_t* man, sh_exc_t exn, gboolean flag);

  static inline
  sh_manager_t* sh_manager_copy(sh_manager_t* man);
  /* Increment the reference counter and return its argument */
  void sh_manager_raise_exception(sh_manager_t* man,
				  sh_exc_t exn, sh_funid_t funid, const char* msg);
  /* raise an exception and put fields
     man->result.flag_exact and man->result.flag_best to
     false
  */

  static inline
  sh_manager_t* sh_manager_copy(sh_manager_t* man)
  { man->count++; return man; }

  
/* ********************************************************************** */
/* III. FPU init */
/* ********************************************************************** */

  gboolean sh_fpu_init(void);
  /* tries to set the FPU rounding-mode towards +oo, returns true if successful */

#ifdef __cplusplus
}
#endif

#endif /* _SH_MANAGER_H */
