#include "ehgraph_internal.h"
#include "ehgraph_util.h"
#include "ehgraph_fun.h"

#include <assert.h>

bool simplify_link(ehgraph_t* graph, enode_t node, size_t l);

/* vx = vy. vy must no label NULL */
ehgraph_t* ehgraph_assign_x_y(ehgraph_t* graph, size_t vx, size_t vy) {
  ehgraph_t* ret = NULL;
  bool bRet = false;

  do {
    CheckNotNULL(graph);
    CheckLessThan(vx, graph->ptrdim);
    CheckLessThan(vy, graph->ptrdim);

    CheckNotEqual(VAR2NODE(graph, vy), NODE_NULL);

    ret = ehgraph_copy(graph);
    CheckNotNULL(ret);

    VAR2NODE(ret, vx) = VAR2NODE(ret, vy);
    ret->closed = false;

    bRet = true;
  } while(0);

  if (!bRet)
    SafeFree(ret);

  return ret;
}


/* vx must point to NULL */
ehgraph_t *ehgraph_assign_x_new(ehgraph_t * graph, size_t vx) {
  ehgraph_t* ret = NULL;
  bool bRet = false;

  do {
    CheckNotNULL(graph);
    CheckLessThan(vx, graph->ptrdim);
    CheckEqual(VAR2NODE(graph, vx), NODE_NULL);

    ret = ehgraph_alloc(graph->size + 1, graph->ptrdim, graph->nLinks);
    CheckNotNULL(ret);

    ret->closed = graph->closed;
    CheckedCall(ehgraph_transfer_info(ret, graph, graph->size + graph->ptrdim));
    /* by default, in a new allocated graph all the links will point to NULL */

    VAR2NODE(ret, vx) = graph->size;

    bRet = true;
  } while(0);

  if (!bRet)
    SafeFree(ret);

  return ret;
}

ehgraph_t *ehgraph_assign_x_free(ehgraph_t * graph, size_t vx) {
  ehgraph_t* ret = NULL;
  bool bRet = false;

  do {
    CheckNotNULL(graph);
    CheckLessThan(vx, graph->ptrdim);

    CheckNotEqual(VAR2NODE(graph, vx), NODE_NULL);

    ret = ehgraph_alloc(graph->size - 1, graph->ptrdim, graph->nLinks);
    CheckNotNULL(ret);
    
    size_t nVars = graph->ptrdim, size = graph->size, nLinks = graph->ptrdim;
    enode_t node = VAR2NODE(graph, vx);

    bRet = true;

    for (enode_t i = 0; i < nVars + node && bRet; ++i)
      bRet = ehgraph_transfer_node_info(ret, i, graph, i);

    for (enode_t i = nVars + node + 1; i < nVars + size && bRet; ++i) 
      bRet = ehgraph_transfer_node_info(ret, i - 1, graph, i);

    if (!bRet)
      break;

    for (enode_t i = 0; i < size; ++i) {
      if (i == node)
	continue;

      int map = i;
      if (i > node)
	--map;

      for (size_t l = 0; l < nLinks; ++l) {
	if (GET_LINK(graph, i, l) == node) 
	  GET_LINK(ret, map, l) = NODE_DANGLING;

	if (!GET_NODE(graph, i).superEdge[l])
	  continue;

	for (size_t ld = 0; ld < nLinks; ++ld) 
	  if (GET_NODE(graph, i).isPredicate[l][ld] && GET_NODE(graph, i).predicateTarget[l][ld] == node) 
	    GET_NODE(ret, map).predicateTarget[l][ld] = NODE_DANGLING;
      }
    }

    ret->closed = false;
  } while(0);

  if (!bRet)
    SafeFree(ret);

  return ret;
}


ehgraph_t*  ehgraph_assign_x_null(ehgraph_t* graph, size_t vx) {
  bool bRet = false;
  ehgraph_t* ret = NULL;

  do {
    CheckNotNULL(graph);
    CheckLessThan(vx, graph->ptrdim);

    ret = ehgraph_copy(graph);
    CheckNotNULL(ret);

    VAR2NODE(ret, vx) = NODE_NULL;
    ret->closed = false;

    bRet = true;
  } while(0);

  if (!bRet)
    SafeFree(ret);

  return ret;
}

/* will return false if the function fails, the result of the call will be left in the last argument */
bool split_needed(ehgraph_t* graph, size_t vx, size_t link, bool* result) {
  bool bRet = false;

  do {
    CheckNotNULL(graph);
    CheckNotNULL(result);
    CheckLessThan(vx, graph->ptrdim);
    CheckLessThan(link, graph->nLinks);

    *result = false;
    if (VAR2NODE_REF(graph, vx).superEdge[link])
      *result = true;

    if (0 == link && VAR2NODE_REF(graph, vx).outDoubleLink)
      *result = true;

    if (1 == link && VAR2NODE_REF(graph, vx).inDoubleLink)
      *result = true;

    bRet = true;
  } while(0);

  return bRet;
}

bool clear_edge_predicates(ehgraph_t* graph, enode_t node, size_t link) {
  bool bRet = false;

  do {
    CheckNotNULL(graph);
    CheckLessThan(node, graph->size);
    CheckLessThan(link, graph->nLinks);

    bRet = true;
  } while(0);

  return bRet;
}

bool clear_super_edge(ehgraph_t* graph, enode_t node, size_t link) {
  bool bRet = false;

  do {
    CheckNotNULL(graph);
    CheckLessThan(node, graph->size);
    CheckLessThan(link, graph->nLinks);

    if (!GET_NODE(graph, node).superEdge[link]) {
      bRet = true;
      break;
    }

    GET_NODE(graph, node).superEdge[link] = false;
    for (size_t ld = 0; ld < graph->nLinks; ++ld) {
      GET_NODE(graph, node).isPredicate[link][ld] = false;
      GET_NODE(graph, node).predicateTarget[link][ld] = false;
    }
      
    bRet = true;
  } while(0);

  return bRet;
}

/* will split the graph from a variable with a link only if needed */
ehgraph_t* ehgraph_split(ehgraph_t* graph, size_t vx, size_t link, enode_t* outNode, bool oo) {
  ehgraph_t* ret = NULL;
  bool bRet = false; 
  
  do {
    CheckNotNULL(graph);
    CheckNotNULL(outNode);
    CheckLessThan(vx, graph->ptrdim);
    CheckLessThan(link, graph->nLinks);
    
    bool splitNeeded = false;
    CheckedCall(split_needed(graph, vx, link, &splitNeeded));

    if (splitNeeded) {

      /* split is needed */
      ret = ehgraph_alloc(graph->size + 1, graph->ptrdim, graph->nLinks);
      CheckNotNULL(ret);

      /* TO-DO: query on the data constraints */
      bool splitOneOne = oo;

      CheckedCall(ehgraph_transfer_info(ret, graph, graph->size + graph->ptrdim));
      
      enode_t node = VAR2NODE(graph, vx), newNode = graph->size;
      enode_t target = GET_LINK(graph, node, link);

      ehgraph_transfer_edge(ret, newNode, graph, node, link);
      
      for (size_t ld = 0; ld < graph->nLinks; ++ld)
	if (GET_NODE(graph, node).isPredicate[link][ld]) {
	  GET_LINK(ret, newNode, ld) = GET_NODE(graph, node).predicateTarget[link][ld];
	}
      

      do {
	if (0 == link && VAR2NODE_REF(graph, vx).outDoubleLink) {
	  GET_LINK(ret, newNode, 1) = node;
	  GET_LINK(ret, target, 1) = newNode;
	  bRet = true;
	  break;
	}

	if (1 == link && VAR2NODE_REF(graph, vx).inDoubleLink) {
	  GET_LINK(ret, newNode, 0) = node;
	  GET_LINK(ret, target, 0) = newNode;
	  bRet = true;
	  break;
	}
	
	bRet = true;
      } while(0);

      Check(bRet);
      bRet = false;

      CheckedCall(clear_super_edge(ret, node, link));
      if (splitOneOne) 
	CheckedCall(clear_super_edge(ret, newNode, link));

      GET_LINK(ret, node, link) = newNode;
      ret->closed = false;

      *outNode = newNode;
      bRet = true;
      break; /* end of split needed */
    } 

    /* split not needed */
    ret = ehgraph_copy(graph);
    CheckNotNULL(ret);

    ret->closed = graph->closed;
    *outNode = GET_LINK(ret, VAR2NODE(ret, vx), link);
    bRet = true;
  } while(0);

  if (!bRet)
    SafeFree(ret);

  return ret;
}

ehgraph_t* ehgraph_assign_x_link_null(ehgraph_t* graph, size_t vx, size_t link) {
  ehgraph_t* ret = NULL;
  bool bRet = false;

  do {
    CheckNotNULL(graph);
    CheckLessThan(vx, graph->ptrdim);
    CheckLessThan(link, graph->nLinks);

    enode_t newNode = 0;
    ret = ehgraph_split(graph, vx, link, &newNode, false);
    CheckNotNULL(ret);

    GET_LINK(ret, VAR2NODE(graph, vx), link) = NODE_NULL;
    ret->closed = false;

    bRet = true;
  } while(0);

  if (!bRet)
    SafeFree(ret);

  return ret;
}

/* x is not NULL but x.link points to null */
ehgraph_t *ehgraph_assign_x_link_y(ehgraph_t * graph, size_t vx, size_t link, size_t vy) {
  ehgraph_t* ret = NULL;
  bool bRet = false;

  do {
    CheckNotNULL(graph);
    CheckLessThan(vx, graph->ptrdim);
    CheckLessThan(vy, graph->ptrdim);
    CheckLessThan(link, graph->nLinks);

    CheckNotEqual(NODE_NULL, VAR2NODE(graph, vx));
    CheckEqual(NODE_NULL, GET_LINK(graph, VAR2NODE(graph, vx), link));

    ret = ehgraph_copy(graph);
    CheckNotNULL(ret);

    GET_LINK(ret, VAR2NODE(ret, vx), link) = VAR2NODE(ret, vy);
    ret->closed = false;

    bRet = true;
  } while(0);

  if (!bRet)
    SafeFree(ret);

  return ret;
}


bool simplify_link(ehgraph_t* graph, enode_t node, size_t l) {
  bool bRet = false;
  
  do {
    CheckNotNULL(graph);
    CheckLessThan(node, graph->size);
    CheckLessThan(l, graph->nLinks);

    if (0 == l) 
      GET_NODE(graph, node).outDoubleLink = false;

    if (1 == l) 
      GET_NODE(graph, node).inDoubleLink = false;

    GET_NODE(graph, node).superEdge[l] = false;

    for (size_t ld = 0; ld < graph->nLinks; ++ld) {
      GET_NODE(graph, node).isPredicate[l][ld] = false;
      GET_NODE(graph, node).predicateTarget[l][ld] = NODE_NULL; //the value doesn't matter that much
    }
    
    bRet = true;
  } while(0);
  
  return bRet;
}



ehgraph_t *ehgraph_assign_x_y_link(ehgraph_t * graph, size_t vx, size_t vy, size_t link) {
  ehgraph_t* ret = NULL;
  bool bRet = false;
  
  do {
    CheckNotNULL(graph);
    CheckLessThan(vx, graph->ptrdim);
    CheckLessThan(vy, graph->ptrdim);
    CheckLessThan(link, graph->nLinks);
    CheckNotEqual(VAR2NODE(graph, vy), NODE_NULL);

    if (NODE_DANGLING == GET_VAR_LINK(graph, vy, link)) {
      ret = ehgraph_bottom();
      CheckNotNULL(ret);
      bRet = true;
      break;
    }

    enode_t newNode = 0;
    ret = ehgraph_split(graph, vy, link, &newNode, false);
    CheckNotNULL(ret);

    VAR2NODE(ret, vx) = newNode;
    ret->closed = false;

    bRet = true;
  } while(0);

  if (!bRet)
    SafeFree(ret);

  return ret;
}
