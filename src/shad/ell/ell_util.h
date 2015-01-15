#ifndef __EHGRAPH_UTIL_H__
#define __EHGRAPH_UTIL_H__

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

/* bitfields */
typedef unsigned int bit_field;

#define EMPTY_SET ((bit_field) 0)

/* WARNING: these functions don't check if pos small enough */
bool bfield_get(bit_field a, size_t pos);
void bfield_set(bit_field *a, size_t pos);
void bfield_unset(bit_field *a, size_t pos);


/* logging. inspiration from wikipedia */
extern FILE* flog;
extern FILE* fdbg;

#define WHERESTR  "[file %s, line %d]: "
#define WHEREARG  __FILE__, __LINE__


/* Debugging macros */
#define DB_ERROR(...) fprintf(flog, "[ERROR] " WHERESTR, WHEREARG); fprintf(flog, __VA_ARGS__); fprintf(flog, "\n");


#define Man_Error(man,exn,funid,message,...)  {char s[100];		\
  sprintf(s, "[ERROR][file %s, line %d]:" message, __FILE__, __LINE__, __VA_ARGS__); \
  hg_manager_raise_exception(man, exn, HG_FUNID_ ## funid, s);  break; }

#define Man_CheckNotNULL(p, man, funid) if (!p) {Man_Error(man, DH_EXC_NOT_IMPLEMENTED, funid, \
							   "null(%s)", #p); break;}

#define Man_CheckArgNotNULL(p, man, funid) if (!p) {Man_Error(man, DH_EXC_INVALID_ARGUMENT, funid, \
							      "null argument(%s)", #p); break;}

#define Man_CheckAllocation(p, man, funid) if (!p) {Man_Error(man, DH_EXC_OUT_OF_SPACE , funid, \
							      "allocation error(%s)", #p);break;}
  
#define Man_CheckedMalloc(v, type, nElem, man, funid) v = (type*) malloc(sizeof(*v) * nElem); Man_CheckAllocation(v, man, funid);

#define Man_CheckedCall(call, man, funid) if (!call) {/*Man_Error(man, DH_EXC_NOT_IMPLEMENTED, funid, "call failed");*/ break;}

  

#define EHGRAPH_DEBUG

#ifdef EHGRAPH_DEBUG
  #define d_printf(...)  fprintf(fdbg,__VA_ARGS__)
#else
  #define d_printf(...)
#endif

/* miscellaneous */
#define MERGE(a,b) a##b
#define UNIQUE_ID() MERGE(unique,__COUNTER__)

#define Check(x) if (!x) {DB_ERROR("check failed");break;}
#define CheckedCall(call) if (!call) {DB_ERROR("call failed");break;}

#define CheckLessThan(a,b)    if (a >= b) {DB_ERROR("out of bounds %s(%d) >= %s(%d)", #a ,a, #b, b); break;}
#define CheckEqual(a,b)       if (a != b) {DB_ERROR("%s(%d) != %s(%d)", #a, a, #b, b); break;}
#define CheckNotEqual(a,b)    if (a == b) {DB_ERROR("%s(%d) == %s(%d)", #a, a, #b, b); break;}
#define CheckNotNULL(a)       if (NULL == a) {DB_ERROR("%s is NULL", #a); break;}
#define CheckedIf(var,call)   bool var; CheckedCall(call); if (var)

#define CheckedMalloc(v, type, nElem) v = (type*) malloc(sizeof(*v) * nElem);  if (!v) {DB_ERROR("allocation error: " #v); break;}
#define SafeFree(p) if (p) {free(p); p = NULL;}

#define MAKE_STRING(var,...) char var[30]; sprintf(var,__VA_ARGS__);

/* boolean vector */

/* true if there is a true value anywhere in the vector */
bool bvect_exists(bool *v, size_t size);

bool random_permutation(int* v, int size);


#endif //__EHGRAPH_UTIL_H__
