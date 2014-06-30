

#include "ucons.h"
#include "ucons_internal.h"

/*
 * Unit test of the printing procedure for ucons domain.
 */

int
main (void) {

  /* Initilisation des valeurs du test */
  ap_manager_t* man_ucons = ucons_manager_alloc();
  ucons_internal_t* pr = ucons_init_from_manager (man_ucons, AP_FUNID_FPRINT, 0);
  
  ucons_t* a = ucons_top (man_ucons, 2, 2);
  ap_linexpr0_t* lexpr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, 4);
  ap_linexpr0_set_coeff_scalar_int(lexpr, 2, 1);
  a = build_const_1(pr, a, lexpr);
  
  /* Call the tested method */
  char* nameofdim[] = { "i", "j", "T1", "T2" };
  ucons_fprint(stdout, man_ucons, a, nameofdim);
  
  return 0;
}
