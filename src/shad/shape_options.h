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


#ifndef __SHAPE_OPTIONS_H_
#define __SHAPE_OPTIONS_H_

/* Options for the selection of domains */

/* constant names for domains */
typedef enum
{
  DOM_BOX = 0, DOM_OCT, DOM_OCT_P11, DOM_OCT_P12, DOM_OCT_P21, DOM_POLY, DOM_POLY_P11, DOM_POLY_P12, DOM_POLY_P21
} dom_t;


/* the domain for existentials (prg variables and heap nodes) */
#if defined(SHAPE_DCONS_BOX)
#define SHAPE_DCONS_DOMAIN  DOM_BOX
#elif defined (SHAPE_DCONS_OCT)
#define SHAPE_DCONS_DOMAIN  DOM_OCT
#elif defined (SHAPE_DCONS_OCT_P11)
#define SHAPE_DCONS_DOMAIN  DOM_OCT_P11
#elif defined (SHAPE_DCONS_OCT_P12)
#define SHAPE_DCONS_DOMAIN  DOM_OCT_P12
#elif defined (SHAPE_DCONS_OCT_P21)
#define SHAPE_DCONS_DOMAIN  DOM_OCT_P21
#elif defined (SHAPE_DCONS_POLY_P11)
#define SHAPE_DCONS_DOMAIN  DOM_POLY_P11
#elif defined (SHAPE_DCONS_POLY_P12)
#define SHAPE_DCONS_DOMAIN  DOM_POLY_P12
#elif defined (SHAPE_DCONS_POLY_P21)
#define SHAPE_DCONS_DOMAIN  DOM_POLY_P21
#else
#define SHAPE_DCONS_DOMAIN  DOM_POLY
#endif

/* the domains for segments (working in parallel) */
/* for lengths */
#if defined(SHAPE_SCONS_LEN)
#define SHAPE_SCONS_LEN 1
#if defined(SHAPE_SCONS_LEN_BOX)
#define SHAPE_SCONS_LEN_DOMAIN DOM_BOX
#elif defined(SHAPE_SCONS_LEN_OCT)
#define SHAPE_SCONS_LEN_DOMAIN DOM_OCT
#else
#define SHAPE_SCONS_LEN_DOMAIN DOM_POLY
#endif
#else
#define SHAPE_SCONS_LEN 0
#endif

/*
 * for universal constraints on data (options are set inside the domain, see
 * below)
 */
#if defined(SHAPE_SCONS_UCONS)
#define SHAPE_SCONS_UCONS 1
#else
#define SHAPE_SCONS_UCONS 0
#endif

/* TODO: for multiset constraints on data */
#if defined(SHAPE_SCONS_MSET)
#define SHAPE_SCONS_MSET 1
#else
#define SHAPE_SCONS_MSET 0
#endif

/*
 * for sum constraints
 */
#if defined(SHAPE_SCONS_LSUM)
#define SHAPE_SCONS_LSUM 1
#else
#define SHAPE_SCONS_LSUM 0
#endif

/* for pretty printer */
extern size_t ushape_number;

#endif /* __SHAPE_OPTIONS_H_ */
