#include "ehgraph_internal.h"
#include "ehgraph_util.h"
#include "ehgraph_fun.h"

#include <string.h>
#include <stdlib.h>

bool ehgraph_register_to_hg_manager(hg_manager_t* manager) {
  if (NULL == manager) 
    return false;

  manager->funptr[HG_FUNID_TOP] = ehgraph_top;
  manager->funptr[HG_FUNID_IS_TOP] = ehgraph_is_top;
  manager->funptr[HG_FUNID_BOTTOM] = ehgraph_bottom;
  manager->funptr[HG_FUNID_IS_BOTTOM] = ehgraph_is_bottom;
  manager->funptr[HG_FUNID_SIZE] = ehgraph_size;
  manager->funptr[HG_FUNID_FPRINT] = ehgraph_fprint;
  manager->funptr[HG_FUNID_IS_EQ] = ehgraph_is_eq;
  return true;
}

void enode_init(enode_info_t* node) {
  if (!node)
    return;

  memset(node->links, 0, sizeof(enode_t) * MAX_LINKS);
  node->inDoubleLink = false;
  node->outDoubleLink = false;
  memset(node->superEdge, 0, sizeof(node->superEdge));
    
  for (int i = 0; i < MAX_LINKS; ++i)
    for (int j = 0; j < MAX_LINKS; ++j) {
      node->isPredicate[i][j] = false;
      node->predicateTarget[i][j] = NODE_NULL; /* the value doesn't matter */
    }     
}


ehgraph_t* ehgraph_alloc(hg_manager_t* manager, size_t size, size_t ptrdim, size_t nlinks) {
  ehgraph_t* ret = NULL;

  do {
    CheckNotNULL(manager);

    ret = (ehgraph_t*) malloc(sizeof(ehgraph_t) + sizeof(enode_info_t) * (size + ptrdim));
    
    Man_CheckAllocation(ret, manager, UNKNOWN);

    ret->closed = false;
    ret->ptrdim = ptrdim;
    ret->size = size;
    ret->nLinks = nlinks;
    ret->linkNames = NULL;
    
    for (size_t i = 0; i < ptrdim + size; ++i) 
      enode_init((ret->info) + i);
      
  } while(0);

  return ret;
}

bool ehgraph_set_pred(ehgraph_t* graph, enode_t node, size_t edge, size_t link, size_t target) {
  if (NULL == graph) {
    DB_ERROR("null argument");
    return false;
  }

  GET_NODE(graph, node).superEdge[edge] = true;
  GET_NODE(graph, node).isPredicate[edge][link] = true;
  GET_NODE(graph, node).predicateTarget[edge][link] = target;
  return true;
}

bool ehgraph_del_pred(ehgraph_t* graph, enode_t node, size_t edge, size_t link) {
  if (NULL == graph) {
    DB_ERROR("null argument");
    return false;
  }

  GET_NODE(graph, node).isPredicate[edge][link] = false;
  return true;
}

ehgraph_t* permute_graph(ehgraph_t* graph, int* perm) {
  ehgraph_t* ret = NULL, *b = NULL;

  do {
    CheckNotNULL(graph);
    CheckNotNULL(perm);

    if (0 != perm[0]) {
      DB_ERROR("the null element should remain unchanged in the permutation");
      break;
    }
    
    b = ehgraph_alloc(NULL, graph->size, graph->ptrdim, graph->nLinks);
    if (!b) {
      DB_ERROR("allocation failed");
      break;
    }

    b->closed = graph->closed;

    for (enode_t i = 0; i < graph->ptrdim; ++i) 
      VAR2NODE(b, i) = perm[ VAR2NODE(graph, i) ];

    for (enode_t i = 0; i < graph->size; ++i) {
      GET_NODE(b, perm[i]).inDoubleLink = GET_NODE(graph, i).inDoubleLink;
      GET_NODE(b, perm[i]).outDoubleLink = GET_NODE(graph, i).outDoubleLink;

      for (size_t l = 0; l < graph->nLinks; ++l) {
	GET_LINK(b, perm[i], l) = perm[ GET_LINK(graph, i, l) ];
	GET_NODE(b, perm[i]).superEdge[l] = GET_NODE(graph, i).superEdge[l];

	for (size_t lp = 0; lp < graph->nLinks; ++lp) {
	  GET_NODE(b, perm[i]).isPredicate[l][lp] = GET_NODE(graph, i).isPredicate[l][lp];
	  GET_NODE(b, perm[i]).predicateTarget[l][lp] = perm[ GET_NODE(graph, i).predicateTarget[l][lp] ];
	}
      }
    }
	   
    ret = b;
  } while(0);
  
  if (!ret)
    SafeFree(b);
  
  return ret;
}

bool ehgraph_transfer_edge(ehgraph_t* dest, enode_t n1, ehgraph_t* src, enode_t n2, size_t link) {
  bool bRet = false; 
  do {
    CheckNotNULL(dest);
    CheckNotNULL(src);
    CheckLessThan(n1, dest->size);
    CheckLessThan(n2, src->size);
    CheckEqual(src->nLinks, dest->nLinks);
    CheckLessThan(link, src->nLinks);

    size_t nL = src->nLinks;
    
    enode_info_t* nDst = &GET_NODE(dest, n1);
    enode_info_t* nSrc = &GET_NODE(src, n2);
    
    if (link == 0) 
      nDst->outDoubleLink = nSrc->outDoubleLink;
    
    if (1 == link)
      nDst->inDoubleLink = nSrc->inDoubleLink;
    

    nDst->links[link] = nSrc->links[link];
    nDst->superEdge[link] = nSrc->superEdge[link];
      
    for (size_t ld = 0; ld < nL; ++ld) {
      nDst->isPredicate[link][ld] = nSrc->isPredicate[link][ld];
      nDst->predicateTarget[link][ld] = nSrc->predicateTarget[link][ld];
    }

    bRet = true;
  } while(0);

  return bRet;
}

ehgraph_t* ehgraph_top(hg_manager_t* manager, size_t intdim, size_t realdim) {
  ehgraph_t* ret = NULL;
  do {
    CheckNotNULL(manager);

    ret = ehgraph_alloc(manager, 0, 0, 0);
    CheckNotNULL(ret);

    ret->isTop = true;
  } while(0);

  return ret;
}

ehgraph_t* ehgraph_bottom(hg_manager_t* manager, size_t intdim, size_t realdim) {
  ehgraph_t* ret = NULL;
  do {
    CheckNotNULL(manager);

    ret = ehgraph_alloc(manager, 0, 0, 0);
    CheckNotNULL(ret);

    ret->isBottom = true;
  } while(0);

  return ret;
}

bool ehgraph_is_top(hg_manager_t* manager, ehgraph_t* graph) {
  bool bRet = false;

  do {
    CheckNotNULL(manager);
    manager->result.flag_exact = false;

    Man_CheckArgNotNULL(graph, manager, IS_TOP);

    bRet = graph->isTop;
    manager->result.flag_exact = true;
  } while(0);

  return bRet;
}


bool ehgraph_is_bottom(hg_manager_t* manager, ehgraph_t* graph) {
  bool bRet = false;

  do {
    CheckNotNULL(manager);
    manager->result.flag_exact = false;

    Man_CheckArgNotNULL(graph, manager, IS_BOTTOM);

    bRet = graph->isBottom;
    manager->result.flag_exact = true;
  } while(0);

  return bRet;
}

size_t ehgraph_size(hg_manager_t* manager, ehgraph_t* graph) {
  size_t ret = 0;

  do {
    CheckNotNULL(manager);
    manager->result.flag_exact = false;

    Man_CheckArgNotNULL(graph, manager, SIZE);

    ret = graph->size;
    manager->result.flag_exact = true;
  } while(0);

  return ret;
}


bool ehgraph_transfer_node_info(ehgraph_t* dest, enode_t n1, ehgraph_t* src, enode_t n2) {
  bool bRet = false;

  do {
    CheckNotNULL(dest);
    CheckLessThan(n1, dest->size + dest->ptrdim);

    CheckNotNULL(src);
    CheckLessThan(n2, src->size + src->ptrdim);
    
    CheckEqual(src->nLinks, dest->nLinks);
    
    size_t nL = src->nLinks;
    enode_info_t* nDst = dest->info + n1;
    enode_info_t* nSrc = src->info + n2;
    
    nDst->outDoubleLink = nSrc->outDoubleLink;
    nDst->inDoubleLink = nSrc->inDoubleLink;
    
    for (size_t l = 0; l < nL; ++l) {
      nDst->links[l] = nSrc->links[l];
      nDst->superEdge[l] = nSrc->superEdge[l];
      
      for (size_t ld = 0; ld < nL; ++ld) {
	  nDst->isPredicate[l][ld] = nSrc->isPredicate[l][ld];
	  nDst->predicateTarget[l][ld] = nSrc->predicateTarget[l][ld];
      }
    }
    
    bRet = true;
  } while(0);

  return bRet;
}

bool ehgraph_transfer_info(ehgraph_t* dest, ehgraph_t* src, size_t lim) {
  bool bRet = false;

  do {

    CheckNotNULL(src);
    CheckNotNULL(dest);
    if (dest->ptrdim != src->ptrdim || dest->nLinks != src->nLinks) {
      DB_ERROR("graphs are not compatible");
      break;
    }

    if (lim > (src->size + src->ptrdim) || lim > (dest->size + dest->ptrdim)) {
      DB_ERROR("invalid arg");
      break;
    }

    size_t nL = src->nLinks;

    bRet = true;
    for (size_t i = 0; i < lim && bRet; ++i) 
      bRet = ehgraph_transfer_node_info(dest, i, src, i);
    
  } while(0);
  
  return bRet;
}

ehgraph_t* ehgraph_copy(ehgraph_t* graph) {
  bool bRet = false;
  ehgraph_t* ret = NULL;

  do {
    CheckNotNULL(graph);

    ret = ehgraph_alloc(NULL, graph->size, graph->ptrdim, graph->nLinks);
    if (NULL == ret) {
      DB_ERROR("allocation error");
      break;
    }

    ret->closed = graph->closed;
    ret->isTop = false;
    ret->isBottom = false;
    CheckedCall(ehgraph_transfer_info(ret, graph, graph->ptrdim + graph->size));
    bRet = true;
  } while(0);

  if (!bRet) 
    SafeFree(ret);
   
  return ret;
}

