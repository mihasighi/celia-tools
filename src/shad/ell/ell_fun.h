#ifndef __EHGRAPH_FUN_H__
#define __EHGRAPH_FUN_H__

#include "hg_manager.h"
#include "ehgraph_internal.h"

/* intdim and realdim are not needed, their are used just to have a definition compliant 
   with the manager */
ehgraph_t* ehgraph_top(hg_manager_t* manager, size_t intdim, size_t realdim); 

ehgraph_t* ehgraph_bottom(hg_manager_t* manager, size_t intdim, size_t realdim);

bool ehgraph_is_top(hg_manager_t* manager, ehgraph_t* graph);

bool ehgraph_is_bottom(hg_manager_t* manager, ehgraph_t* graph);


/* WARNING: will modify the graph provided as argument */
/* If this is not desired there are not many changes that are needed */
ehgraph_t* ehgraph_close(hg_manager_t* man, ehgraph_t* graph);

bool ehgraph_is_eq(hg_manager_t* man, ehgraph_t* graph1, ehgraph_t* graph2);

bool ehgraph_is_leq(ehgraph_t* graph1, ehgraph_t* graph2, bool* result);

#endif //__EHGRAPH_FUN_H__
