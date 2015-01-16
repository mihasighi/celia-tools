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
 * sh_manager.c: global manager passed to all functions
 */

#include <math.h>
#include <assert.h>
#include "sh_manager.h"

/* ********************************************************************** */
/* I. Types */
/* ********************************************************************** */

/* ====================================================================== */
/* I.1 Identifying functions */
/* ====================================================================== */

const char* sh_name_of_funid[SH_FUNID_SIZE2] = {
  "unknown",
  "copy",
  "free",
  "size",
  "minimize",
  "canonicalize",
  "hash",
  "approximate",
  "print",
  "printdiff",
  "dump",
  "serialize_raw",
  "deserialize_raw",
  "bottom",
  "top",
  "of_box",
  "dimension",
  "is_bottom",
  "is_top",
  "is_leq",
  "is_eq",
  "is_dimension_unconstrained",
  "sat_tcons",
  "sat_tform",
  "to_tcons_array",
  "to_tform_array",
  "meet",
  "meet_array",
  "meet_tcons_array",
  "join",
  "join_array",
  "assign_texpr_array",
  "substitute_texpr_array",
  "add_dimensions",
  "remove_dimensions",
  "permute_dimensions",
  "forget_array",
  "expand",
  "fold",
  "widening",
  "closure",
  "unknown",
  "change_environment",
  "rename"
};


/* ====================================================================== */
/* I.2 Exceptions */
/* ====================================================================== */

const char* sh_name_of_exception[SH_EXC_SIZE] = {
  "NONE",
  "TIMEOUT",
  "OUT_OF_SPACE",
  "OVERFLOW",
  "INVALID_ARGUMENT",
  "NOT_IMPLEMENTED"
};

/* ********************************************************************** */
/* I.2 Options: constructor and destructor */
/* ********************************************************************** */

void sh_funopt_init(sh_funopt_t* opt)
{
  opt->algorithm = 0;
  opt->timeout = 0;
  opt->max_object_size = 0;
  opt->flag_exact_wanted = FALSE;
  opt->flag_best_wanted = FALSE;
}

void sh_option_init(sh_option_t* opt)
{
  sh_funid_t funid;
  sh_exc_t exn;

  for (funid=0; funid<SH_FUNID_SIZE; funid++){
    sh_funopt_init(&opt->funopt[funid]);
  }
  for (exn=0; exn<SH_EXC_SIZE; exn++){
    opt->abort_if_exception[exn] = TRUE;
  }
  // opt->scalar_discr = AP_SCALAR_DOUBLE;
}

/* ********************************************************************** */
/* II. User Functions */
/* ********************************************************************** */

/*
 * Constructor and destructor for result.
 */
sh_exclog_t* sh_exc_cons(sh_exc_t exn, 
			 sh_funid_t funid,
			 const char* msg, 
			 sh_exclog_t* tail)
{
  sh_exclog_t* head = (sh_exclog_t*)malloc(sizeof(sh_exclog_t));
  head->exn = exn;
  head->funid = funid;
  head->msg = g_strdup(msg ? msg : "");
  head->tail = tail;
  return head;
}

void sh_exclog_free(sh_exclog_t* head)
{
  sh_exclog_t* p = head;
  while (p!=NULL) {
    sh_exclog_t* tail = p->tail;
    free(p->msg);
    free(p);
    p = tail;
  }
}

void sh_result_add_exception(sh_result_t* result, sh_exc_t exn, sh_funid_t funid, const char* msg)
{
  result->exclog = sh_exc_cons(exn,funid,msg,result->exclog);
  result->exn = exn;
}

void sh_result_init(sh_result_t* result)
{
  result->exclog = NULL;
  result->exn = SH_EXC_NONE;
  result->flag_exact = FALSE;
  result->flag_best = FALSE;
}
void sh_result_clear(sh_result_t* result)
{
  sh_exclog_free(result->exclog);
  sh_result_init(result);
}

/* 
 * Constructor and destructor for manager 
 */
sh_manager_t* 
sh_manager_alloc(char* library, char* version, 
		 void* internal, 
		 void (*internal_free)(void*))
{
  sh_manager_t* man;

  // assert(sizeof(gboolean)==1);

  man = (sh_manager_t*)malloc(sizeof(sh_manager_t));
  man->library = library;
  man->version = version;
  man->internal = internal;
  man->internal_free = internal_free;
  man->count = 1;
  sh_option_init(&man->option);
  sh_result_init(&man->result);

  /* init local informations about the program */
  man->filenv = sh_filenv_alloc_empty(man);
  man->typenv = sh_typenv_alloc_empty(man);
  man->varenv = sh_stack_alloc_empty(man);

  return man;
}

void 
sh_manager_free(sh_manager_t* man)
{
  if (man->count>1){
    man->count--;
  }
  else {
    if (man->internal != NULL){
      man->internal_free(man->internal);
      man->internal = NULL;
    }
    sh_result_clear(&man->result);
    /* free allocated data for program */
    sh_filenv_free(man->filenv);
    sh_typenv_free(man->typenv);
    sh_stack_free(man->varenv);
    free(man);
  }
}

/*
 * Other manager functions
 */
const char* sh_manager_get_library(sh_manager_t* man)
{ return man->library; }
const char* sh_manager_get_version(sh_manager_t* man)
{ return man->version; }
sh_filenv_t* sh_manager_get_filenv(sh_manager_t* man)
{ return  man->filenv; }
sh_typenv_t* sh_manager_get_typenv(sh_manager_t* man)
{ return man->typenv; }
sh_stack_t* sh_manager_get_varenv(sh_manager_t* man)
{ return man->varenv; }

int 
sh_manager_get_frame_of_loc(sh_manager_t* man, sh_loc_t* loc)
{
  assert (NULL != man);
  if (loc == NULL || 
      man->varenv->frames == NULL)
    return -1;
  for (guint i = 0; i < man->varenv->frames->len; i++)
    {
      sh_frame_t* f = g_ptr_array_index (man->varenv->frames,i);
      if (sh_loc_cmp(f->start,loc) <= 0 &&
	  sh_loc_cmp(loc, f->end))
	return i;
    }
  return -1;
}

void 
sh_manager_fdump(FILE* stream, sh_manager_t* man)
{
  fprintf (stream, "Manager %s (version %s):\n", man->library, man->version);
  fprintf (stream, "\tFiles:\n");
  sh_filenv_fdump(stream, man->filenv);
  fprintf (stream, "\tTypes:\n");
  sh_typenv_fdump(stream, man->typenv);
  fprintf (stream, "\tStack:\n");
  sh_stack_fdump(stream, man->varenv);
}

sh_funopt_t sh_manager_get_funopt(sh_manager_t* man, sh_funid_t funid)
{
  if (funid<SH_FUNID_SIZE) 
    return man->option.funopt[funid];
  else {
    fprintf(stderr,"sh_manager.c: sh_manager_get_funopt: funid should be less than SH_FUNID_SIZE\n");
    abort();
  }
}
gboolean sh_manager_get_abort_if_exception(sh_manager_t* man, sh_exc_t exn)
{ return man->option.abort_if_exception[exn]; }
gboolean sh_manager_get_flag_exact(sh_manager_t* man)
{ return man->result.flag_exact; }
gboolean sh_manager_get_flag_best(sh_manager_t* man)
{ return man->result.flag_best; }


void sh_manager_set_funopt(sh_manager_t* man, sh_funid_t funid, sh_funopt_t* funopt)
{ if (funid<SH_FUNID_SIZE) man->option.funopt[funid] = *funopt; }

void sh_manager_set_abort_if_exception(sh_manager_t* man, sh_exc_t exn, gboolean flag)
{ man->option.abort_if_exception[exn] = flag; }
void sh_manager_clear_exclog(sh_manager_t* man)
{
  sh_exclog_free(man->result.exclog);
  man->result.exclog = NULL;
}  

void sh_manager_raise_exception(sh_manager_t* man, 
				sh_exc_t exn, 
				sh_funid_t funid, 
				const char* msg)
{
  gboolean pabort;

  if (exn!=SH_EXC_NONE){
    pabort = man->option.abort_if_exception[exn];
    if (pabort){
      fprintf(stderr,"Apron: Abort because of following exception:\nexception %s in function %s:\n%s\n",
	      sh_name_of_exception[exn], sh_name_of_funid[funid],
	      msg);
      abort();
    }
    else {
      sh_result_add_exception(&man->result,exn,funid,msg);
      man->result.flag_exact = man->result.flag_best = FALSE;
    }
  }
  return;
}


/* ********************************************************************** */
/* III. FPU init */
/* ********************************************************************** */

/* simple run-time test that fpu behaves correctly */
static gboolean test_fpu(void)
{
  int i;
  long double d = 1., dd;
  /* find the minimal long double, as the fixpoint of x -> x/2 with rounding
     towards +oo;
     the max iteration value should be enough for 128-bit floating-point */
  for (i=0;i<5000000;i++) {
    dd = d;
    d /= 2;
    if (d==dd || d==0.) break;
  }
  /* fails if flush to 0 */
  if (d!=dd) { fprintf(stderr,"test_fpu failed test #1 after %i iterations\n",i); return FALSE; }
  /* fails if long double rounding is not towards +oo */
  if (d*0.25!=dd) { fprintf(stderr,"test_fpu failed test #2\n"); return FALSE; }
  /* fails if double rounding is not towards +oo */
  if ((double)d<dd) { fprintf(stderr,"test_fpu failed test #3\n"); return FALSE; }
  /* fails if float rounding is not towards +oo */
  if ((float)d<dd) { fprintf(stderr,"test_fpu failed test #4\n"); return FALSE; }
  return TRUE;
}

#if defined(__ppc__)
gboolean sh_fpu_init(void) 
{ 
  __asm volatile ("mtfsfi 7,2");
  return test_fpu();
}

#elif defined(__linux) || defined (__APPLE__)
#include <fenv.h>
gboolean sh_fpu_init(void) 
{ 
  if (!fesetround(FE_UPWARD)) return test_fpu();
  fprintf(stderr,"could not set fpu rounding mode: fesetround failed\n");
  return FALSE;
}

#elif defined(__FreeBSD__) || defined(sun)
#include <ieeefp.h>
gboolean sh_fpu_init(void)
{ 
  fpsetround(FP_RP); 
  return test_fpu();
}

#else
gboolean sh_fpu_init(void)
{
  fprintf(stderr,"could not set fpu rounding mode: platform not supported\n");
  return FALSE;
}

#endif
