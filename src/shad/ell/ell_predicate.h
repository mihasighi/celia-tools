#ifndef __EHGRAPH_PREDICATE__
#define __EHGRAPH_PREDICATE__

#include "ehgraph_internal.h"
#include "ehgraph_util.h"

typedef size_t predicate_t;

extern size_t g_nPredicates;

/* if the predicate does not exist it will be inserted in the collection */
size_t getPredicateIndex(size_t linkIndex, enode_t target);

void getPredicate(int index, size_t* linkIndex, enode_t* target);

#endif // __EHGRAPH_PREDICATE__
