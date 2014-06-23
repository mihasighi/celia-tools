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


#ifndef __SHAPE_MACROS_H_
#define __SHAPE_MACROS_H_

/* *INDENT-OFF* */
#ifdef __cplusplus
extern          "C"
{
#endif
  /* *INDENT-ON* */

  /* invalid argument exception */
#define arg_assert(cond,action)                                         \
do { if (!(cond)) {                                                   \
    char buf_[1024];                                                  \
    snprintf(buf_,sizeof(buf_),                                       \
             "assertion (%s) failed in %s at %s:%i",                  \
             #cond, __func__, __FILE__, __LINE__);            \
    ap_manager_raise_exception(pr->man,AP_EXC_INVALID_ARGUMENT,       \
                               pr->funid,buf_);                       \
    fprintf(stdout,"%s\n",buf_); fflush(stdout); \
    action }                                                          \
} while(0)

  /* feedback */
#define FEEDBACK(msg,action)                                       \
  {                                                          \
    fprintf(stdout,"[shad:] warning: %s\n", msg); fflush(stdout);  \
    action                                        \
  } 

  /* internal errors */
#define ERROR(msg,action)                                       \
  do {                                                          \
    char buf_[1024];                                            \
    snprintf(buf_,sizeof(buf_),                                 \
             "exception (%s) raised in %s at %s:%i",            \
              msg, __func__, __FILE__, __LINE__);               \
    ap_manager_raise_exception(pr->man,AP_EXC_NOT_IMPLEMENTED,  \
                               pr->funid,buf_);                 \
    fprintf(stdout,"%s\n",buf_); fflush(stdout);  \
    action                                        \
    pr->error_++;                                 \
  } while (0)

  /* malloc with safe-guard */
#define checked_malloc(ptr,t,size,nb,action)                               \
do {                                                                  \
  (ptr) = (t*)malloc(size*(nb));                                 \
  if (!(ptr)) {                                                       \
    char buf_[1024];                                                  \
    snprintf(buf_,sizeof(buf_),                                       \
             "cannot allocate %s[%lu] for %s in %s at %s:%i",         \
             #t, (long unsigned)(nb), #ptr,                           \
             __func__, __FILE__, __LINE__);                           \
    ap_manager_raise_exception(pr->man,AP_EXC_OUT_OF_SPACE,           \
                               pr->funid,buf_);                       \
    action }                                                          \
} while(0)


  /* calloc with safe-guard, needed by hash table */
#define checked_calloc(ptr,t,size,nb,action)                               \
do {                                                                  \
  (ptr) = (t*)calloc(size,(nb));                                 \
  if (!(ptr)) {                                                       \
    char buf_[1024];                                                  \
    snprintf(buf_,sizeof(buf_),                                       \
             "cannot allocate %s[%lu] for %s in %s at %s:%i",         \
             #t, (long unsigned)(nb), #ptr,                           \
             __func__, __FILE__, __LINE__);                           \
    ap_manager_raise_exception(pr->man,AP_EXC_OUT_OF_SPACE,           \
                               pr->funid,buf_);                       \
    action }                                                          \
} while(0)

#define checked_realloc(ptr,t,size,nb,action)                               \
do {                                                                  \
  (ptr) = (t*)realloc(ptr,size*(nb));                                 \
  if (!(ptr)) {                                                       \
    char buf_[1024];                                                  \
    snprintf(buf_,sizeof(buf_),                                       \
             "cannot allocate %s[%lu] for %s in %s at %s:%i",         \
             #t, (long unsigned)(nb), #ptr,                           \
             __func__, __FILE__, __LINE__);                           \
    ap_manager_raise_exception(pr->man,AP_EXC_OUT_OF_SPACE,           \
                               pr->funid,buf_);                       \
    action }                                                          \
} while(0)


  /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */


#endif /* SHAPE_MACROS_H_ */
