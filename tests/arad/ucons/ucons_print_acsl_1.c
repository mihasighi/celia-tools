/***********************************************************************
 *  CINV Library / Shape Domain 
 *  Copyright (C) 2014 
 *  
 * LIAFA (University of Paris Diderot and CNRS)
 *  
 * you can redistribute it and/or modify it under the terms of the GNU
 * Lesser General Public License as published by the Free Software
 * Foundation, version 3.
 * 
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 * 
 * See the GNU Lesser General Public License version 3.
 * for more details (enclosed in the file LICENSE).   
 * 
 **********************************************************************/
 
#include "ucons.h"
#include "ucons_internal.h"
#include "apron2shape.h"

/*
 * Unit test of the printing procedure for ucons domain.
 * 
 * Test of ucons_fprint_acsl 
 * with build_const_2.
 *  
 */

int
main (void) {

  /* Initilisation of test values */
  ap_manager_t* man_ucons = ucons_manager_alloc();
  ucons_internal_t* pr = ucons_init_from_manager (man_ucons, AP_FUNID_FPRINT, 0);
  
  ucons_t* a = ucons_top (man_ucons, 2,1);
  ap_linexpr0_t* lexpr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, 4);
  ap_linexpr0_set_coeff_scalar_int(lexpr, 2, 3);
  a = build_const_2(pr, a, lexpr);//Type of build_const 'assertion.c'
  
  /* Call the tested method */
  char* nameofdim[] = { "i","j","T1","T2"};
  
  sh_print_set_acsl();
  ucons_fprint(stdout,man_ucons, a, nameofdim);
  
  return 0;
}
