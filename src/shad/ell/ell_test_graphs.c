#include "ehgraph_test_graphs.h"
#include "ehgraph_internal.h"


ehgraph_t* create_graph_1(hg_manager_t* man) {
  ehgraph_t* g = NULL;

  do {
    assert(man);
    size_t ptrdim = 0;
    size_t size = 5;
    size_t links = 3;

    g = ehgraph_alloc(man, size, ptrdim, links);
    if (!g) {
      DB_ERROR("allocation error");
      break;
    }

    for (size_t i = 1; i < size; ++i) {
      for(size_t l = 0; l < links; ++l)
  	GET_NODE(g, i).links[l] = i;
    }

    for (size_t l = 0; l < links; ++l)
      GET_NODE(g, l + 1).links[l] = l + 2;
      
  } while(0);

  return g;  
}


ehgraph_t* create_graph_2(hg_manager_t* man) {
  ehgraph_t* g = NULL;

  do {
    size_t ptrdim = 4;
    size_t size = 6;
    size_t links = 3;

    g = ehgraph_alloc(man, size, ptrdim, links);
    if (!g) {
      DB_ERROR("allocation error");
      break;
    }

    VAR2NODE(g, 0) = 2;
    VAR2NODE(g, 1) = 2;
    VAR2NODE(g, 2) = 3;
    VAR2NODE(g, 3) = 1;

    for (size_t i = 1; i < size; ++i) {
      for(size_t l = 0; l < links; ++l)
  	GET_NODE(g, i).links[l] = i;
    }

    GET_NODE(g,2).links[0] = 3;
    GET_NODE(g,2).superEdge[0] = true;
    ehgraph_set_pred(g, 2, 0, 1, 1);
    ehgraph_set_pred(g, 2, 0, 2, NODE_NULL);

    GET_NODE(g,1).links[0] = NODE_NULL;
    GET_NODE(g,1).outDoubleLink = true;
    ehgraph_set_pred(g, 2, 0, 1, 1);

    GET_NODE(g, 3).links[2] = 4;
    GET_NODE(g, 3).superEdge[2] = true;
    ehgraph_set_pred(g, 3, 2, 1, 1);

    GET_NODE(g, 5).links[1] = NODE_NULL;
    GET_NODE(g, 5).superEdge[1] = true;
    ehgraph_set_pred(g, 5, 1, 0, 0);  

  } while(0);

  return g;  
}

ehgraph_t* create_graph_3(hg_manager_t* man) {
  ehgraph_t* g = NULL;

  do {
    size_t ptrdim = 1;
    size_t size = 4;
    size_t links = 3;
       
    g = ehgraph_alloc(man, size, ptrdim, links);
    if (NULL == g) {
      DB_ERROR("allocation error");
      break;
    }

    VAR2NODE(g, 0) = 2;

    for (size_t i = 1; i < size; ++i) {
      for(size_t l = 0; l < links; ++l)
  	GET_NODE(g, i).links[l] = i;
    }

    GET_NODE(g, 1).links[0] = NODE_NULL;

    GET_NODE(g, 2).links[0] = 3;
    
    GET_NODE(g, 3).links[0] = 0;
    GET_NODE(g, 3).superEdge[0] = true;

    ehgraph_set_pred(g, 3, 0, 1, NODE_NULL);
    ehgraph_set_pred(g, 3, 0, 2, NODE_NULL);

    GET_NODE(g, 3).links[1] = NODE_NULL;
    GET_NODE(g, 3).links[2] = NODE_NULL;
  } while(0);

  return g;
}

ehgraph_t* create_graph_4(hg_manager_t* man) {
  ehgraph_t* g = NULL;

  do {
    size_t ptrdim = 1;
    size_t size = 7;
    size_t links = 3;
       
    g = ehgraph_alloc(man, size, ptrdim, links);
    if (NULL == g) {
      DB_ERROR("allocation error");
      break;
    }

    g->closed = false;

    VAR2NODE(g, 0) = 1;
    /* VAR2NODE(g, 1) = 6; */

    for (size_t i = 1; i < size; ++i) {
      for(size_t l = 0; l < links; ++l)
  	GET_NODE(g, i).links[l] = i;
    }

    GET_LINK(g, 1, 0) = 2;
    //GET_LINK(g, 1, 1) = 6;
    
    GET_LINK(g, 2, 0) = 3;
    GET_LINK(g, 2, 2) = GET_LINK(g, 2, 1) = NODE_NULL;

    GET_LINK(g, 1, 2) = 4;
    
    GET_LINK(g, 4, 2) = 5;
    GET_LINK(g, 5, 2) = 0;

    GET_LINK(g, 4, 0) = GET_LINK(g, 4, 1) = 1;
    GET_LINK(g, 5, 0) = GET_LINK(g, 5, 1) = 1;

    GET_LINK(g, 6, 0) = 5;
  } while(0);

  return g;
}


ehgraph_t* create_graph_5(hg_manager_t* man) {
  ehgraph_t* g = NULL;

  do {
    size_t ptrdim = 1;
    size_t size = 4;
    size_t links = 3;
       
    g = ehgraph_alloc(man, size, ptrdim, links);
    if (NULL == g) {
      DB_ERROR("allocation error");
      break;
    }

    VAR2NODE(g, 0) = 3;
    /* VAR2NODE(g, 1) = 6; */

    for (size_t i = 1; i < size; ++i) {
      for(size_t l = 0; l < links; ++l)
  	GET_NODE(g, i).links[l] = i;
    }

    GET_LINK(g, 1, 0) = 2;
    GET_LINK(g, 1, 1) = 3;
    GET_LINK(g, 1, 2) = 3;

    GET_LINK(g, 2, 0) = 3;
    GET_LINK(g, 3, 0) = NODE_NULL;

    GET_LINK(g, 2, 1) = 1;
    GET_LINK(g, 3, 1) = 2;
    
    GET_LINK(g, 2, 2) = 3;
    GET_LINK(g, 3, 2) = 3;
  } while(0);

  return g;
}



ehgraph_t* create_graph_6(hg_manager_t* man) {
  ehgraph_t* g = NULL;

  do {
    size_t ptrdim = 1;
    size_t size = 4;
    size_t links = 2;
       
    g = ehgraph_alloc(man, size, ptrdim, links);
    if (NULL == g) {
      DB_ERROR("allocation error");
      break;
    }

    VAR2NODE(g, 0) = 3;

    for (size_t i = 1; i < size; ++i) {
      for(size_t l = 0; l < links; ++l)
  	GET_NODE(g, i).links[l] = i;
    }

    GET_LINK(g, 1, 0) = NODE_NULL;
    GET_LINK(g, 2, 0) = 1;
    GET_LINK(g, 3, 0) = 2;

    GET_LINK(g, 1, 1) = NODE_NULL;
    GET_LINK(g, 2, 1) = NODE_NULL;
    GET_LINK(g, 3, 1) = NODE_NULL;
    
    ehgraph_set_pred(g, 3, 0, 1, 0);
    ehgraph_set_pred(g, 1, 0, 1, 0);
  } while(0);

  return g;
}

ehgraph_t* create_graph_7(hg_manager_t* man) {
 ehgraph_t* ret = NULL;

  do {
    size_t ptrdim = 1;
    size_t size = 3;
    size_t links = 2;
       
    ret = ehgraph_alloc(man, size, ptrdim, links);
    CheckNotNULL(ret);
    ret->linkNames = calloc(2, sizeof(char*));
    ret->linkNames[0] = strdup("next");
    ret->linkNames[1] = strdup("prev");
    

    VAR2NODE(ret, 0) = 1;
    ehgraph_set_pred(ret, 1, 0, 1, 2);
  } while(0);
  
  return ret;
}


ehgraph_t* create_graph_8(hg_manager_t* man) {
 ehgraph_t* ret = NULL;

  do {
    size_t ptrdim = 1;
    size_t size = 3;
    size_t links = 2;
       
    ret = ehgraph_alloc(man, size, ptrdim, links);
    CheckNotNULL(ret);

    VAR2NODE(ret, 0) = 1;
    GET_NODE(ret, 1).outDoubleLink = true;
    ehgraph_set_pred(ret, 1, 0, 1, 2);
  } while(0);
  
  return ret;
}

ehgraph_t* create_graph_9(hg_manager_t* man) {
  ehgraph_t* ret = NULL;

  do {
    size_t ptrdim = 1;
    size_t size = 2;
    size_t links = 2;
       
    ret = ehgraph_alloc(man, size, ptrdim, links);
    CheckNotNULL(ret);

    VAR2NODE(ret, 0) = 1;
    GET_NODE(ret, 1).outDoubleLink = true;
    GET_NODE(ret, 1).superEdge[0] = true;
    ehgraph_set_pred(ret, 1, 0, 1, 0);

    ret->closed = true;
  } while(0);

  return ret;
}

ehgraph_t* create_graph_10(hg_manager_t* man) {
  ehgraph_t* ret = NULL;

  do {
    size_t ptrdim = 1;
    size_t size = 2;
    size_t links = 2;
       
    ret = ehgraph_alloc(man, size, ptrdim, links);
    CheckNotNULL(ret);

    VAR2NODE(ret, 0) = 1;
    GET_NODE(ret, 1).outDoubleLink = true;
    
    ret->closed = true;
  } while(0);

  return ret;
}


ehgraph_t* create_graph_11(hg_manager_t* man) {
  ehgraph_t* ret = NULL;

  do {
    size_t ptrdim = 2;
    size_t size = 2;
    size_t links = 2;
       
    ret = ehgraph_alloc(man, size, ptrdim, links);
    CheckNotNULL(ret);

    VAR2NODE(ret, 0) = 1;
    VAR2NODE(ret, 1) = NODE_NULL;

    ret->closed = true;
  } while(0);

  return ret;
}




ehgraph_t* create_graph_12(hg_manager_t* man) {
  ehgraph_t* ret = NULL;

  do {
    size_t ptrdim = 2;
    size_t size = 3;
    size_t links = 2;
       
    ret = ehgraph_alloc(man, size, ptrdim, links);
    CheckNotNULL(ret);

    VAR2NODE(ret, 0) = 1;
    VAR2NODE(ret, 1) = 2;
    
    ret->closed = true;
  } while(0);

  return ret;
}


int g_nGraphs = 12;

ehgraph_t* create_graph(hg_manager_t* man, int index) {
  assert(man);
  if (index < 1 || index > g_nGraphs)
    return NULL;

  if (1 == index)
    return create_graph_1(man);
  if (2 == index)
    return create_graph_2(man);
  if (3 == index)
    return create_graph_3(man);
  if (4 == index) 
    return create_graph_4(man);
  if (5 == index)
    return create_graph_5(man);
  if (6 == index)
    return create_graph_6(man);
  if (7 == index)
    return create_graph_7(man);
  if (8 == index)
    return create_graph_8(man);
  if (9 == index)
    return create_graph_9(man);
  if (10 == index)
    return create_graph_10(man);
  if (11 == index)
    return create_graph_11(man);
  if (12 == index)
    return create_graph_12(man);

  return NULL;
}
