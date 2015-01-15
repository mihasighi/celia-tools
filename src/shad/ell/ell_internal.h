#ifndef __EHGRAPH_INTERNAL_H_
#define __EHGRAPH_INTERNAL_H_

#include "hg_manager.h"
#include "ehgraph_util.h"
#include <stdlib.h>

#define MAX_LINKS (5)
#define MAX_SIMPLE_NODES (0)
#define MAX_GRAPH_SIZE (10)

typedef size_t enode_t;
typedef size_t var_t;


struct _enode_info_t {
  /* links[0] will hold the link that will induce the list structure for nodes and
     and the value for the variables (first part of info) */
  enode_t links[MAX_LINKS]; 
  bool outDoubleLink;
  bool inDoubleLink; 
  bool superEdge[MAX_LINKS];

  bool isPredicate[MAX_LINKS][MAX_LINKS];
  enode_t predicateTarget[MAX_LINKS][MAX_LINKS];
};

typedef struct _enode_info_t enode_info_t;

struct _ehgraph_t {
  bool closed; 

  size_t ptrdim; /* number of variables */
  size_t size; /* number of nodes */

  char** linkNames;
  size_t nLinks; /* the number of links that each node has. nLinks = 1 <=> single linked list */
  bool isBottom;
  bool isTop;
  enode_info_t info[]; 
};

typedef struct _ehgraph_t ehgraph_t;

#define NODE_NULL (0)
#define NODE_DANGLING (1)

#define VAR2NODE(graph, var) ( (graph)->info[(var)].links[0] )
#define GET_NODE(graph,node) ( (graph)->info[ (graph)->ptrdim + (node) ] )
#define GET_LINK(graph,node,link) ( GET_NODE(graph,node).links[(link)] )
#define GET_VAR_LINK(graph,var,link) ( GET_NODE(graph, VAR2NODE(graph, var)).links[(link)] )

#define VAR2NODE_REF(graph, var) ( GET_NODE( (graph), VAR2NODE( (graph), (var) )) )

/* false only if graph is NULL, otherwise in the resulting graph all it means that the edge 'edge' (ex next)
   is a super edge with all the nodes from it pointing with the link 'link' to the node target */
bool ehgraph_set_pred(ehgraph_t* graph, enode_t node, size_t edge, size_t link, size_t target);
bool ehgraph_del_pred(ehgraph_t* graph, enode_t node, size_t edge, size_t link);


enode_t var_get_node(ehgraph_t* graph, var_t var);
void var_set_node(ehgraph_t* graph, var_t var, enode_t value);

/* TO-DO: fix bug. Doesn't work well when multiple super edges from node */
void ehgraph_fprint(FILE * stream, hg_manager_t * man, ehgraph_t * graph, char **name_of_dim);

void ehgraph_fprint_dot(FILE* stream, ehgraph_t* a, char* name);

ehgraph_t* ehgraph_alloc(hg_manager_t* manager, size_t size, size_t ptrdim, size_t nlinks);

size_t ehgraph_size(hg_manager_t* man, ehgraph_t* a);

ehgraph_t* permute_graph(ehgraph_t* graph, int* perm);

bool ehgraph_transfer_edge(ehgraph_t* dest, enode_t n1, ehgraph_t* src, enode_t n2, size_t link);

bool ehgraph_transfer_node_info(ehgraph_t* dest, enode_t n1, ehgraph_t* src, enode_t n2);

/* copy the from the info array of the source graph to a destination graph. will apply only to the first
   'limit' positions. the first ptrdim positions are for variables */
bool ehgraph_transfer_info(ehgraph_t* dest, ehgraph_t* src, size_t limit);

ehgraph_t* ehgraph_copy(ehgraph_t* graph);

ehgraph_t* ehgraph_split(ehgraph_t* graph, size_t vx, size_t link, enode_t* outNode, bool oo);

ehgraph_t *ehgraph_assign_x_y(ehgraph_t * a, size_t vx,
			      size_t vy); /* DONE */

/* if y->link is a dangling pointer bottom will be returned */
ehgraph_t *ehgraph_assign_x_y_link(ehgraph_t * a, size_t vx, size_t vy, size_t link); /* DONE */

ehgraph_t *ehgraph_assign_x_null(ehgraph_t * a, size_t vx); /* DONE */

ehgraph_t *ehgraph_assign_x_new(ehgraph_t * a, size_t vx); /* DONE */

ehgraph_t *ehgraph_assign_x_free(ehgraph_t * a, size_t vx); /* DONE */

ehgraph_t *ehgraph_assign_x_link_null(ehgraph_t * a, size_t vx, size_t link); /* DONE */

ehgraph_t *ehgraph_assign_x_link_y(ehgraph_t * a, size_t vx, size_t link, size_t vy); /* DONE */

bool ehgraph_register_to_hg_manager(hg_manager_t* manager);

#endif // __EHGRAPH_INTERNAL_H_
