#include "ehgraph_util.h"

#include <time.h>
#include <stdlib.h>

bool bfield_get(bit_field a, size_t pos) {
  return (a>>pos)&1;
}

void bfield_set(bit_field* a, size_t pos) {
  *a |= (1<<pos);
}

void bfield_unset(bit_field* a, size_t pos) {
  *a &= ~((bit_field)(1<<pos));
}

bool bvect_exists(bool* v, size_t size) {
  for (size_t i = 0; i < size; ++i)
    if (v[i])
      return true;
  return false;
}

bool random_permutation(int* v, int size) {
  bool bRet = false;

  do {
    CheckNotNULL(v);
    
    srand(13);
    for (int i = 1; i < size; ++i) {
      int p = rand() % (i + 1);
      if (p != i) {
	int aux = v[p];
	v[p] = v[i];
	v[i] = aux;
      }
    }

    bRet = true;
  } while(0);

  return bRet;
}

					 
