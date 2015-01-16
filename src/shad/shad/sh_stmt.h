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
 * sh_stmt.h: management of the program stack, binding of vars to dimensions
 */

#ifndef _SH_STMT_H_
#define _SH_STMT_H_

#include "sh_manager.h"

#ifdef __cplusplus
extern "C" {
#endif


/* ********************************************************************** */
/* I. Types */
/* ********************************************************************** */

/* ====================================================================== */
/* I.1 Statements */
/* ====================================================================== */

  typedef enum {
    SH_STMT_SKIP = 0,
    SH_STMT_ASSUME_EQ,
    SH_STMT_ASSUME_NE,
    SH_STMT_NEW,
    SH_STMT_FREE,
    SH_STMT_ASSIGN,
    SH_STMT_PCALL,
    SH_STMT_PRETURN,
    SH_STMT_UNKNOWN /* NOT USED */
  } sh_stmt_e;

  typedef struct sh_stmt_s {
    sh_loc_t* loc;  /* location for this statement */
    size_t* frame;  /* frame of this statement */
    sh_stmt_e kind; /* kind of statement */
    size_t tid;     /* type-id of the value assigned/tested */
    union {
      struct {
	size_t left;     /* left hand side variable */
	GArray* offset_l; /* left hand side dereferencing, array of field-ids */
	size_t right;    /* right hand side variable */
	GArray* offset_r;
      } binary;
      struct {
	char* pname;     /* procedure name */
	size_t pframe;   /* frame of the called procedure */
	GArray args;     /* actuals for this procedure call, only var-ids */
        GArray argk;     /* kind of actuals TRUE: value, FALSE: reference */
	size_t* ret;     /* return var, if != null */
      } call;
    } info;
  } sh_stmt_t;

/* ********************************************************************** */
/* II. User Functions */
/* ********************************************************************** */

/* ====================================================================== */
/* II.1 Statements */
/* ====================================================================== */

  sh_stmt_t* sh_stmt_alloc(sh_manager_t* man, sh_stmt_e kind, sh_loc_t* loc);
  /* Constructor a new variable */
  void sh_stmt_free(sh_stmt_t* a);
  /* Destructor */

  /* Building different kinds of statements */
  sh_stmt_t* sh_stmt_new_skip(sh_manager_t* man, sh_loc_t* loc);
  sh_stmt_t* sh_stmt_new_assume(sh_manager_t* man, sh_loc_t* loc,
				sh_stmt_e kind,
				size_t x, size_t y);
  sh_stmt_t* sh_stmt_new_new(sh_manager_t* man, sh_loc_t* loc, 
			     size_t x);
  sh_stmt_t* sh_stmt_new_free(sh_manager_t* man, sh_loc_t* loc, 
			      size_t x);
  sh_stmt_t* sh_stmt_new_assign_x_y(sh_manager_t* man, sh_loc_t* loc, 
				    size_t x, size_t y);
  /* Assignment x:=y */
  sh_stmt_t* sh_stmt_new_assign_x_f_y(sh_manager_t* man, sh_loc_t* loc, 
				      size_t x, size_t f, size_t y);
  /* Assignment x.f:=y */
  sh_stmt_t* sh_stmt_new_assign_x_y_f(sh_manager_t* man, sh_loc_t* loc, 
				      size_t x, size_t y, size_t f);
  /* Assignment x:=y.f */
  sh_stmt_t* sh_stmt_new_pcall(sh_manager_t* man, sh_loc_t* loc,
			       char* pname, size_t pframe,
			       GArray args, size_t ret);
  /* Statement ret := pname(args) */

  void sh_stmt_fdump (FILE* stream, sh_manager_t* man, sh_stmt_t* a);
  /* Dumping information about a variable */

 
#ifdef __cplusplus
}
#endif

#endif
