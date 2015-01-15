#include "ehgraph_fun.h"
#include "ehgraph_util.h"

#include <assert.h>
#include <string.h>

void printQueue(char* message, enode_t* q, size_t start, size_t limit) {
  assert(q);
  
  if (message)
    d_printf("%s\n", message);

  if (start + 1 >= limit)
    return;

  d_printf("%d", q[start]);

  for (size_t i = start + 1; i < limit; ++i)
    d_printf(" %d", q[i]);
  d_printf("\n");
}

bool compare_node_predicates(ehgraph_t* graph, enode_t n1, enode_t n2, size_t sl) {
  if (!graph) {
    DB_ERROR("null argument");
    return false;
  }

  for (size_t l = 0; l < graph->nLinks; ++l) {
    if (l == sl)
      continue;

    if (GET_LINK(graph, n1, l) != GET_LINK(graph, n2, l))
      return false;
  }
  
  return true;
}

bool compare_node_predicates_dl(ehgraph_t* graph, enode_t n1, enode_t n2) {
  if (!graph) {
    DB_ERROR("null argument");
    return false;
  }

  for (size_t l = 2; l < graph->nLinks; ++l) {
    if (GET_LINK(graph, n1, l) != GET_LINK(graph, n2, l))
      return false;
  }
  
  return true;
}



int go_simple(ehgraph_t* graph, enode_t start, size_t link, bool* simple, bool* simpleDL) {
  int nRet = -1;

  do {
    CheckNotNULL(graph);
    CheckNotNULL(simple);
    CheckNotNULL(simpleDL);

    bool s = simple[start], sD = simpleDL[start];

    for (; s == simple[start] && sD == simpleDL[start]; start = GET_LINK(graph, start, link));

    nRet = start;
  } while(0);

  return nRet;
}

void makeSuperEdges(ehgraph_t* source, enode_t start, enode_t firstSimple, 
		    enode_t nodeDest, ehgraph_t* gDest, int* mapping, bit_field* visited) {

  if (!source || !gDest || !mapping || !visited) {
    DB_ERROR("null argument %p %p %p %p", source, gDest, mapping, visited);
    return;
  }

  int mStart = mapping[start];
  int mDest = mapping[nodeDest];
  

  for (size_t l = 0; l < source->nLinks; ++l) 
    if (GET_LINK(source, start, l) == firstSimple) {

      GET_LINK(gDest, mStart, l) = mDest;

      bfield_set(visited, l);

      GET_NODE(gDest, mStart).superEdge[l] = true;

      for (size_t lP = 0; lP < source->nLinks; ++lP)
	if (GET_LINK(source, start, lP) != firstSimple) {
	  GET_NODE(gDest, mStart).isPredicate[l][lP] = true;
	  enode_t target = GET_LINK(source, firstSimple, lP);
	  GET_NODE(gDest, mStart).predicateTarget[l][lP] = mapping[ target ];
	}
    }
}

bool copy_edge(ehgraph_t* dest, ehgraph_t *src, enode_t node, size_t l, int* mapping) {
  bool bRet = false;

  do {
    CheckNotNULL(dest);
    CheckNotNULL(src);
    CheckNotNULL(mapping);

    int mNode = mapping[node];
    enode_t target = GET_LINK(src, node, l);
    int mTarget = mapping[ target ];

    GET_LINK(dest, mNode, l) = mTarget;
    
    bool sDL = (0 == l && GET_NODE(src, node).outDoubleLink);
    bool sE = GET_NODE(src, node).superEdge[ l ];

    if (!sDL && !sE) {
      bRet = true;
      break;
    }

    if (sE)
      GET_NODE(dest, mNode).superEdge[l] = true;
    if (sDL) {
      GET_NODE(dest, mNode).outDoubleLink = true;
      GET_NODE(dest, mTarget).inDoubleLink = true;
    }

    for (size_t lt = 0; lt < src->nLinks; ++lt) {
      GET_NODE(dest, mNode).isPredicate[ l ][ lt ] = GET_NODE(src, node).isPredicate[ l ][ lt ];
      GET_NODE(dest, mNode).predicateTarget[ l ][ lt ] = mapping[ GET_NODE(src, node).predicateTarget[ l ][ lt ] ];
    }
  
    bRet = true;
  } while(0);

  return bRet;
}

bool make_super_edge(ehgraph_t* src, enode_t node, size_t l, enode_t firstSimple, enode_t target
		     , ehgraph_t* dest, int* mapping)  {

  bool bRet = false;

  do {
    CheckNotNULL(src);
    CheckNotNULL(dest);
    CheckNotNULL(mapping);

    int mNode = mapping[node], mTarget = mapping[target];
    GET_LINK(dest, mNode, l) = mTarget;
    GET_NODE(dest, mNode).superEdge[l] = true;

    for (size_t lp = 0; lp < src->nLinks; ++lp) 
      if (l != lp) {
	GET_NODE(dest, mNode).isPredicate[l][lp] = true;
	GET_NODE(dest, mNode).predicateTarget[l][lp] = mapping[ GET_LINK(src, firstSimple, lp) ];
      }
    
    bRet = true;
  } while(0);

  return bRet;
}

bool make_super_edge_dl(ehgraph_t* src, enode_t node, enode_t firstSimple, enode_t target
		     , ehgraph_t* dest, int* mapping)  {

  bool bRet = false;

  do {
    CheckNotNULL(src);
    CheckNotNULL(dest);
    CheckNotNULL(mapping);

    int mNode = mapping[node], mTarget = mapping[target];
    GET_LINK(dest, mNode, 0) = mTarget;
    GET_LINK(dest, mTarget, 1) = mNode;
    GET_NODE(dest, mNode).outDoubleLink = true;
    GET_NODE(dest, mTarget).inDoubleLink = true;

    for (size_t lp = 2; lp < src->nLinks; ++lp) {
	GET_NODE(dest, mNode).isPredicate[0][lp] = true;
	GET_NODE(dest, mTarget).isPredicate[1][lp] = true;

	GET_NODE(dest, mNode).predicateTarget[0][lp] = mapping[ GET_LINK(src, firstSimple, lp) ];
	GET_NODE(dest, mTarget).predicateTarget[1][lp] = mapping[ GET_LINK(src, firstSimple, lp) ];
      }
    
    bRet = true;
  } while(0);

  return bRet;
}


bool create_structure(ehgraph_t* dest, ehgraph_t* src, int* mapping, bool* garbage, bool* simple, bool* simpleDL) {
  bool bRet = false;

  do {
    CheckNotNULL(dest);
    CheckNotNULL(src);
    CheckNotNULL(mapping);
    CheckNotNULL(garbage);
    CheckNotNULL(simple);
    CheckNotNULL(simpleDL);

    
    bool allOk = true;
    for (enode_t i = 0; i < src->size && allOk; ++i) {
      if (garbage[i] || simple[i] || simpleDL[i])
	continue;

      for (size_t l = 0; l < src->nLinks; ++l) {
	enode_t t = GET_LINK(src, i, l);
	
	bool concat = simple[t] || (simpleDL[t] && 0 == l);
	
	if (!concat) {
	  if (!copy_edge(dest, src, i, l, mapping))
	    allOk = false;

	  continue;
	}

	enode_t cutPoint = go_simple(src, t, l, simple, simpleDL);

	if (simple[t]) {
	  if (!make_super_edge(src, i, l, t, cutPoint, dest, mapping))
	    allOk = false;
	} else {
	  if (!make_super_edge_dl(src, i, t, cutPoint, dest, mapping))
	    allOk = false;
	}
      }
    }
    
    bRet = true;
  } while(0);

  return bRet;
}

/* transform any go-return edge in a double link edge without predicates on it */
bool make_dl_edges(ehgraph_t* graph) {
  bool bRet = false;

  do {
    CheckNotNULL(graph);

    for (enode_t i = 0; i < graph->size; ++i) {
      enode_t t = GET_LINK(graph, i, 0);
      
      if (i == GET_LINK(graph, t, 1) && i != t) {
	GET_NODE(graph, i).outDoubleLink = true;
	GET_NODE(graph, t).inDoubleLink = true;
      }
    }
    
    bRet = true;
  } while(0);

  return bRet;
}

/* determines the accesible nodes (in the argument garbage_
   the nodes labeled by variables will become cutpoints (in both simple and dl).  */
bool bfs(ehgraph_t* graph, bool* garbage, bool* simple, bool* simpleDL) {
  bool bRet = false;
  enode_t* q = NULL;
  

  do {
    CheckNotNULL(graph);
    CheckNotNULL(garbage);
    CheckNotNULL(simple);
    CheckNotNULL(simpleDL);

    size_t size = graph->size, nLinks = graph->nLinks, nVar = graph->ptrdim;    
    q = (enode_t*) malloc(sizeof(enode_t) * size);
    
    if (NULL == q) {
      DB_ERROR("allocation error");
      break;
    }

    assert(sizeof(bool) == 1);
    memset(garbage, true, sizeof(bool) * size);

    size_t l = 0, r = 0;
    garbage[NODE_NULL] = false;
    q[r++] = NODE_NULL;
    simple[NODE_NULL] = false;
    simpleDL[NODE_NULL] = false;

    garbage[NODE_DANGLING] = false;
    q[r++] = NODE_DANGLING;
    simple[NODE_DANGLING]= false;
    simpleDL[NODE_DANGLING] = false;

    for (enode_t i = 0; i < nVar; ++i) {
      enode_t t = VAR2NODE(graph, i);
      if (garbage[t]) {
	garbage[t] = false;
	simple[t] = false;
	simpleDL[t] = false;
	q[r++] = t;
      }
    }
    printQueue("Nodes labeled with variables", q, 0, r);

    int nrVars = r;

    for (l = 0; l < r; ++l) {
      enode_info_t* node = &GET_NODE(graph, q[l]);
      for (size_t link = 0; link < nLinks; ++link) {
	enode_t target = node->links[link];
	if (garbage[target]) {
	  garbage[target] = false;
	  q[r++] = target;
	}
      }
    }

    printQueue("Nodes reached but not labeled by variables", q, nrVars, r);

    bRet = true;
  } while(0);

  SafeFree(q);
  return bRet;
}


bool make_cutpoints_pred_target(ehgraph_t* graph, bool* garbage, bool* simple, bool* simpleDL) {
  bool bRet = false;

  do {
    CheckNotNULL(graph);
    CheckNotNULL(garbage);
    CheckNotNULL(simple);
    CheckNotNULL(simpleDL);

    size_t size = graph->size, nLinks = graph->nLinks;

    for (enode_t i = 0; i < size; ++i) {
      if (garbage[i])
	continue;

      enode_info_t* node = &GET_NODE(graph, i);
      
      for (size_t l = 0; l < nLinks; ++l)
	if (node->superEdge[l] || (l == 0 && node->outDoubleLink)) {
	  for (size_t p = 0; p < nLinks; ++p) 
	    if (node->isPredicate[l][p]) {
	      enode_t t = node->predicateTarget[l][p];
	      simple[t] = false;
	      simpleDL[t] = false;
	    }
	}
    }

    bRet = true;
  } while(0);

  return bRet;
}


bool determine_in_nodes(ehgraph_t* graph, bool* garbage, bool* simple, int *inNode, int *inEdge) {
  bool bRet = false;

  do {
    CheckNotNULL(graph);
    CheckNotNULL(garbage);
    CheckNotNULL(simple);
    CheckNotNULL(inNode);
    CheckNotNULL(inEdge);

    size_t size = graph->size, nLinks = graph->nLinks;
    memset(inNode, -1, sizeof(int) * size);
    memset(inEdge, -1, sizeof(int) * size);

    for (enode_t i = 0; i < size; ++i) {
      if (garbage[i]) /* garbage */
	continue;

      for (size_t l = 0; l < nLinks; ++l) {
	enode_t target = GET_LINK(graph, i, l);
 
	if (-1 == inNode[target]) 
	  inNode[target] = i;
	
	if (-1 == inEdge[target]) 
	  inEdge[target] = l;

	if (inNode[target] != (int) i || inEdge[target] != (int) l)
	  simple[target] = false;
      }
    }

    bRet = true;
  } while(0);

  return bRet;
}


bool determine_dl_simple_nodes(ehgraph_t* graph, bool* garbage, bool* simpleDL) {
  bool bRet = false;
  int *inNext = NULL, *inPrev = NULL;

  do {
    CheckNotNULL(graph);
    CheckNotNULL(garbage);
    CheckNotNULL(simpleDL);

    size_t nSize = graph->size, nLinks = graph->nLinks;

    CheckedMalloc(inNext, int, nSize);
    CheckedMalloc(inPrev, int, nSize);

    memset(inNext, -1, sizeof(inNext) * nSize);
    memset(inPrev, -1, sizeof(inPrev) * nSize);

    for (enode_t i = 0; i < nSize; ++i) {
      if (garbage[i])
	continue;

      int dest = GET_LINK(graph, i, 0);
      if (-1 != inNext[dest])
	simpleDL[dest] = false;
      
      inNext[dest] = i;

      dest = GET_LINK(graph, i, 1);
      if (-1 != inPrev[dest]) 
	simpleDL[dest] = false;

      inPrev[dest] = i;
    }
    
    for (enode_t i = 0; i < nSize; ++i) 
      if (simpleDL[i]) {
	simpleDL[i] = false;
	do {
	  if (-1 == inNext[i])
	    break;

	  if (GET_LINK(graph, i, 1) != (enode_t) inNext[i])
	    break;

	  if (-1 == inPrev[i]) {
	    if (GET_LINK(graph, i, 0) != NODE_NULL)
	      break;
	  } else {
	    if (GET_LINK(graph, i, 0) != (enode_t) inPrev[i])
	      break;
	  }

	  if (GET_NODE(graph, i).inDoubleLink) {
	    bool ok = true;
	    enode_t n = GET_LINK(graph, i, 1);
	    for (size_t l = 0; l < nLinks && ok; ++l)
	      if (GET_NODE(graph, n).isPredicate[0][l] && GET_NODE(graph, n).predicateTarget[0][l] != GET_LINK(graph, i, l))
		ok = false;

	    if (!ok)
	      break; /* cutpoint */
	  }

	  if (GET_NODE(graph, i).outDoubleLink) {
	    bool ok = true;

	    for (size_t l = 0; l < nLinks && ok; ++l)
	      if (GET_NODE(graph, i).isPredicate[0][l] && GET_NODE(graph, i).predicateTarget[0][l] != GET_LINK(graph, i, l))
		ok = false;

	    if (!ok)
	      break; /* cutpoint */
	  }
	  
	  simpleDL[i] = true;
	} while(0);
      }

    bRet = true;
  } while(0);
  
  SafeFree(inNext);
  SafeFree(inPrev);

  return bRet;
}

bool self_loops(ehgraph_t* graph, bool* garbage, bool* simple, bool* simpleDL) {
  bool bRet = false; 
  
  do {
    CheckNotNULL(graph);
    CheckNotNULL(garbage);
    CheckNotNULL(simple);
    CheckNotNULL(simpleDL);

    for (enode_t i = 0; i < graph->size; ++i) {
      for (size_t l = 0; l < graph->nLinks; ++l)
	if (GET_LINK(graph, i, l) == i) {
	  simple[i] = false;
	  simpleDL[i] = false;
	}
    }
    bRet = true;
  } while(0);

  return bRet;
}


bool path_simple(ehgraph_t* graph, bool* garbage, bool* simple, bool* simple2) {
  bool bRet = false;

  do {
    CheckNotNULL(graph);
    CheckNotNULL(garbage);
    CheckNotNULL(simple);
    CheckNotNULL(simple2);

    size_t size = graph->size, nLinks = graph->nLinks;

    for (enode_t i = 0; i < size; ++i) {
      if (garbage[i] || simple[i])
	continue;
      
      for (size_t l = 0; l < nLinks; ++l) {
	enode_t target = GET_LINK(graph, i, l);
	
	if (!simple[ target ])
	  continue;

	size_t nSimple = 0;
	enode_t node;

	for (node = target; simple[node]; node = GET_LINK(graph, node, l), nSimple++) {
	  if (!compare_node_predicates(graph, target, node, l)) {
    	    if (nSimple <= MAX_SIMPLE_NODES) {
    	      for (enode_t nIt = target; nIt != node; nIt = GET_LINK(graph, nIt, l))
    		simple2[nIt] = false;
    	    }

    	    simple2[node] = false;
    	    target = GET_LINK(graph, node, l);
    	    nSimple = 0;
	  }
	}

    	if (nSimple <= MAX_SIMPLE_NODES) {
    	  for (enode_t nIt = target; nIt != node; nIt = GET_LINK(graph, nIt, l))
    	    simple2[nIt] = false;
    	}
      }
    }
    
    bRet = true;
  } while(0);

  return bRet;
}

bool path_simple_dl(ehgraph_t* graph, bool* garbage, bool* simpleDL, bool* simple2) {
  bool bRet = false;

  do {
    CheckNotNULL(graph);
    CheckNotNULL(garbage);
    CheckNotNULL(simpleDL);
    CheckNotNULL(simple2);

    size_t size = graph->size, nLinks = graph->nLinks;

    for (enode_t i = 0; i < size; ++i) {
      if (garbage[i] || simpleDL[i])
	continue;

      enode_t target = GET_LINK(graph, i, 0);

      if (!simpleDL[target])
	continue;

      simple2[i] = false;

      size_t nSimple = 0;
      enode_t node = target;

      for (enode_t node = target; simpleDL[node]; node = GET_LINK(graph, node, 0), nSimple++) {
	if (!compare_node_predicates_dl(graph, target, node)) {
	  if (nSimple <= MAX_SIMPLE_NODES) {
	    for (enode_t nIt = target; nIt != node; nIt = GET_LINK(graph, nIt, 0))
	      simple2[nIt] = false;
	  }

	  simple2[node] = false;
	  target = GET_LINK(graph, node, 0);
	  nSimple = 0;
	}
      }

      if (nSimple <= MAX_SIMPLE_NODES) {
	for (enode_t nIt = target; nIt != node; nIt = GET_LINK(graph, nIt, 0))
	  simple2[nIt] = false;
      }
      
    }
    
    bRet = true;
  } while(0);

  return bRet;
}

bool merge_simple_arrays(size_t size, bool* simple, bool* simpleDL, bool* simple2) {
  bool bRet = false;
  do {
    CheckNotNULL(simple);
    CheckNotNULL(simpleDL);
    CheckNotNULL(simple2);

    for (size_t i = 0; i < size; ++i) {
      simple[i] &= simple2[i];
      simpleDL[i] &= simple2[i];
    }

    bRet = true;
  } while(0);

  return bRet;
}

bool print_node_statistics(ehgraph_t* graph, bool* garbage, bool* simple, bool* simpleDL) {
  bool bRet = false; 

  do {
    CheckNotNULL(graph);
    CheckNotNULL(garbage);
    CheckNotNULL(simple);
    CheckNotNULL(simpleDL);

    d_printf("Category of nodes:\n");
    for (enode_t i = 0; i < graph->size; ++i) {
      d_printf("%d -> ", i);
      if (garbage[i]) {
	d_printf("garbage\n");
	continue;
      }

      if (simple[i]) {
	d_printf("simple\n");
	continue;
      }

      if (simpleDL[i]) {
	d_printf("simple DL\n");
	continue;
      }

      d_printf("cutpoint\n");
    }

    bRet = true;
  } while(0);

  return bRet;
}

/* returns the size of the graph, -1 if error */
int create_mapping(ehgraph_t* graph, int* mapping, bool* garbage, bool* simple, bool* simpleDL) {
  int nRet = -1;

  do {
    CheckNotNULL(graph);
    CheckNotNULL(mapping);
    CheckNotNULL(garbage);
    CheckNotNULL(simple);
    CheckNotNULL(simpleDL);

    size_t size = graph->size, nLinks = graph->nLinks;
    memset(mapping, -1, sizeof(*mapping) * size);
    
    int index = 0; 
    for (enode_t i = 0; i < size; ++i) 
      if (!garbage[i] && !simple[i] && !simpleDL[i]) 
	mapping[i] = index++;

    nRet = index;
  } while(0);

  return nRet;
}

bool set_vars(ehgraph_t* dest, ehgraph_t* src, int * mapping) {
  bool bRet = false;

  do {
    CheckNotNULL(dest);
    CheckNotNULL(src);
    CheckNotNULL(mapping);

    size_t nVars = src->ptrdim;

    for (size_t v = 0; v < nVars; ++v) 
      VAR2NODE(dest, v) = mapping[ VAR2NODE(src, v) ];

    bRet = true;
  } while(0);

  return bRet;
}


/* TO-DO:  */
/* Double-link edges */
ehgraph_t* ehgraph_close(hg_manager_t* man, ehgraph_t* graph) {
  
  bool bRet = false;
  size_t size = 0, nLinks = 0; /* graph parameters */

  bool* garbage = NULL;
  bool *simple = NULL, *simpleDL = NULL, *simple2 = NULL;
  int *inNode = NULL, *inEdge = NULL;
  int* mapping = NULL;  
  ehgraph_t* ret = NULL;

  do {
    CheckNotNULL(man);
    Man_CheckArgNotNULL(graph, man, CLOSURE);

    size = graph->size;
    nLinks = graph->nLinks;
    
    Man_CheckedMalloc(garbage, bool, size, man, CLOSURE);
    Man_CheckedMalloc(simple, bool, size, man, CLOSURE);
    Man_CheckedMalloc(simpleDL, bool, size, man, CLOSURE);

    assert(1 == sizeof(bool));
    memset(simple, true, sizeof(bool) * size);
    memset(simpleDL, true, sizeof(bool) * size);

    Man_CheckedCall(bfs(graph, garbage, simple, simpleDL), man, CLOSURE);

    Man_CheckedCall(self_loops(graph, garbage, simple, simpleDL), man, CLOSURE);
    Man_CheckedCall(make_cutpoints_pred_target(graph, garbage, simple, simpleDL), man, CLOSURE);
    Man_CheckedCall(determine_dl_simple_nodes(graph, garbage, simpleDL), man, CLOSURE);

    Man_CheckedMalloc(inNode, int, size, man, CLOSURE);
    Man_CheckedMalloc(inEdge, int, size, man, CLOSURE);
    Man_CheckedCall(determine_in_nodes(graph, garbage, simple, inNode, inEdge), man, CLOSURE);

    Man_CheckedMalloc(simple2, bool, size, man, CLOSURE);
    assert(1 == sizeof(bool));
    memset(simple2, true, sizeof(bool) * size);
    
    Man_CheckedCall(path_simple(graph, garbage, simple, simple2), man, CLOSURE);
    Man_CheckedCall(path_simple_dl(graph, garbage, simpleDL, simple2), man, CLOSURE);
    Man_CheckedCall(merge_simple_arrays(size, simple, simpleDL, simple2), man, CLOSURE);
    
    Man_CheckedCall(print_node_statistics(graph, garbage, simple, simpleDL), man, CLOSURE);
   
    Man_CheckedMalloc(mapping, int, size, man, CLOSURE);
    int newSize = create_mapping(graph, mapping, garbage, simple, simpleDL);
    
    if (-1 == newSize) {
      //Man_Error(man, DH_EXC_NOT_IMPLEMENTED, CLOSURE, "call failed");
      break;
    }

    /* TO-DO return top if many nodes */
    if (newSize > MAX_GRAPH_SIZE) {
      ret = ehgraph_top(man, 0, 0);
      bRet = true;
      break;
    }

    ret = ehgraph_alloc(man, newSize, graph->ptrdim, nLinks);    
    if (NULL == ret) {
      DB_ERROR("allocation error");
      break;
    }
    
    /* TO-DO: chains of simple nodes */
    CheckedCall(create_structure(ret, graph, mapping, garbage, simple, simpleDL));
    CheckedCall(set_vars(ret, graph, mapping));
    CheckedCall(make_dl_edges(ret));
    
    bRet = true;
  } while(0);

  SafeFree(garbage);
  SafeFree(inNode);
  SafeFree(inEdge);
  SafeFree(mapping);
  SafeFree(simple);
  SafeFree(simple2) ;
  SafeFree(simpleDL);

  if (!bRet && NULL != ret) {
    //TO-DO clean ret graph
  }

  return bRet ? ret : NULL;
}
