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
 * noll.h: open declaration for the noll abstract domain
 */


#ifndef __NOLL_H
#define __NOLL_H

#include "sh_manager.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* ============================================================ */
  /* I.1 Manager */
  /* ============================================================ */

  void noll_manager_init(sh_manager_t* man);
  /* Local setting of the manager for noll */


#ifdef __cplusplus
}
#endif

#endif /* __NOLL_H */
