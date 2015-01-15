#include "ehgraph_predicate.h"
#include "ehgraph_fun.h"
#include "ehgraph_util.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


/* TO-DO: implement hash tables for predicates and dynamic arrays in enode_info */
#define MAX_PREDICATES (31) 

size_t g_linkIndex[MAX_PREDICATES];
size_t g_target[MAX_PREDICATES];

size_t g_nPredicates = 0;


size_t getPredicateIndex(size_t linkIndex, enode_t target) {
  for (size_t i = 0; i < g_nPredicates; ++i)
    if (g_linkIndex[i] == linkIndex && g_target[i] == target)
      return i;
  assert(g_nPredicates < MAX_PREDICATES);
  int ret = g_nPredicates++;
  g_linkIndex[ret] = linkIndex;
  g_target[ret] = target;
  return ret;
}

void getPredicate(int index, size_t* linkIndex, enode_t* target) {
  *linkIndex = g_linkIndex[index];
  *target = g_target[index];
}

bool ehgraph_iso_dfs(ehgraph_t* graph1, ehgraph_t* graph2,
		     int node1, int node2, int* as, int* invas, bool* result) {
		    
  bool bRet = false;

  do {
    CheckNotNULL(graph1);
    CheckNotNULL(graph2);
    CheckNotNULL(result);
    CheckNotNULL(as);
    CheckNotNULL(invas);

    *result = false;
    bRet = true;

    if ((-1 == as[node1]) ^ (-1 == invas[node2]))
      break;
    
    if (-1 == as[node1]) {
      as[node1] = node2;
      invas[node2] = node1;
    } else {
      *result = (as[node1] == node2 && invas[node2] == node1);
      break;
    }

    *result = true;
    *result &= (GET_NODE(graph1, node1).outDoubleLink == GET_NODE(graph2, node2).outDoubleLink);
    *result &= (GET_NODE(graph1, node1).inDoubleLink == GET_NODE(graph2, node2).inDoubleLink);
    
    for (size_t l = 0; l < graph1->nLinks && bRet && (*result); ++l) {
      bool pR = false;
      bRet &= ehgraph_iso_dfs(graph1, graph2, GET_LINK(graph1, node1, l), GET_LINK(graph2, node2, l), as, invas, &pR);
      (*result) &= pR;

      if ((*result) && GET_NODE(graph1, node1).superEdge[l] && GET_NODE(graph2, node2).superEdge[l]) {
	for (size_t lp = 0; lp < graph1->nLinks && bRet && (*result); ++lp) {
	  (*result) &= (GET_NODE(graph1, node1).isPredicate[l][lp] == GET_NODE(graph2, node2).isPredicate[l][lp]);

	  if ((*result) && GET_NODE(graph1, node1).isPredicate[l][lp]) {
	    bool pR = false;
	    bRet &= ehgraph_iso_dfs(graph1, graph2, GET_NODE(graph1, node1).predicateTarget[l][lp],
				    GET_NODE(graph2, node2).predicateTarget[l][lp], as, invas, &pR);
	    (*result) &= pR;
	    
	  }
	}
      }
    }

  } while(0);

  return bRet;
}

bool ehgraph_is_eq(hg_manager_t* man, ehgraph_t* graph1, ehgraph_t* graph2) {
  int *as = NULL, *invas = NULL;
  bool bRet = false, closed1 = false, closed2 = false;

  do {
    CheckNotNULL(man);
    man->result.flag_exact = false;
    Man_CheckArgNotNULL(graph1, man, IS_EQ);
    Man_CheckArgNotNULL(graph2, man, IS_EQ);

    man->result.flag_exact = true;

    if (graph1->isTop ^ graph2->isTop) 
      break;

    if (graph1->isBottom ^ graph2->isBottom) 
      break;

    if (graph1->isTop || graph1->isBottom) {
      bRet = true;
      break;
    }

    man->result.flag_exact = false;
    if (!graph1->closed) {
      graph1 = ehgraph_close(man, graph1);
      Man_CheckNotNULL(graph1, man, IS_EQ);
      closed1 = true;
    }

    if (!graph2->closed) {
      graph2 = ehgraph_close(man, graph2);
      Man_CheckNotNULL(graph2, man, IS_EQ);
      closed2 = true;
    }
 
    if (graph1->size != graph2->size || graph1->ptrdim != graph2->ptrdim || graph1->nLinks != graph2->nLinks) {
      man->result.flag_exact = true;
      break;
    }

    man->result.flag_exact = false;
    size_t size = graph1->size, nVars = graph1->ptrdim;
    Man_CheckedMalloc(as, int, size, man, IS_EQ);
    Man_CheckedMalloc(invas, int, size, man, IS_EQ);
    memset(as, -1, sizeof(*as) * size);
    memset(invas, -1, sizeof(*invas) * size);

    man->result.flag_exact = true;
    bRet = true;
    for(enode_t i = 0; i < nVars && bRet && man->result.flag_exact; ++i) {
      bool pR = false;
      man->result.flag_exact &= ehgraph_iso_dfs(graph1, graph2, VAR2NODE(graph1, i), VAR2NODE(graph2, i), as, invas, &pR);
      bRet &= pR;
    }

  } while(0);

  if (closed1)
    SafeFree(graph1);
  if (closed2)
    SafeFree(graph2);

  SafeFree(as);
  SafeFree(invas);

  return bRet;
}


bool ehgraph_is_leq(ehgraph_t* graph1, ehgraph_t* graph2, bool* result) {
  return false;
  /* bool bRet = false; */
  /* ehgraph_t* copyG2 = NULL; */

  /* do { */
  /*   CheckNotNULL(graph1); */
  /*   CheckNotNULL(graph2); */
  /*   CheckNotNULL(result); */
    
  /*   bRet = true; */
  /*   *result = false; */

  /*   if (graph1->isBottom ^ graph2->isBottom)  */
  /*     break; */

  /*   if (graph1->isTop ^ graph2->isTop)  */
  /*     break; */

  /*   if (graph1->isTop || graph1->isBottom) { */
  /*     *result = true; */
  /*     break; */
  /*   } */

  /*   if (graph1->nLinks != graph2->nLinks)  */
  /*     break; */

  /*   for (size_t i = 0; i < graph1->ptrdim && bRet; ++i)  */
  /*     if (NODE_NULL == VAR2NODE(graph1, i) && i < graph2->ptrdim && NODE_NULL !=  */
  /* 	  VAR2NODE(graph2, i)) { */
	
  /* 	if (NULL == copyG2) { */
  /* 	  copyG2 = ehgraph_copy(graph2); */
  /* 	  if (NULL == copyG2) { */
  /* 	    DB_ERROR("copy failed"); */
  /* 	    bRet = false; */
  /* 	    break; */
  /* 	  } */
  /* 	} */

  /* 	VAR2NODE(copyG2, i) = NODE_NULL; */
  /* 	copyG2->closed = false; */
  /*     } */

  /*   if (!bRet) */
  /*     break; */

  /*   bRet = false; */
  /*   CheckedCall(ehgraph_is_eq(graph1, NULL == copyG2 ? graph2 : copyG2, result)); */
  /*   bRet = true; */
  /* } while(0); */

  /* SafeFree(copyG2); */
  /* return bRet; */
}
