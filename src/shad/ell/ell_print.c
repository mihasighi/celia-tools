#include "ehgraph_internal.h"
#include "ehgraph_util.h"
#include "ehgraph_predicate.h"

#include <assert.h>

const char* get_edge_dot_style(int index) {
  if (0 == index)
    return "solid";

  if (1 == index) 
    return "dashed";
  else 
    return "dotted";  
}

/* print edge to stream */ 
void printEdgeArg(FILE* stream, size_t edgeIndex, bool superEdge, char** names) {
  fprintf(stream, "[");
  if (superEdge)
    fprintf(stream, "color=red, ");
  fprintf(stream, "style=%s", get_edge_dot_style(edgeIndex));

  if (edgeIndex > 0) {
    if (NULL != names)
      fprintf(stream, ", label=\"%s\"", names[edgeIndex]);
    else
      fprintf(stream, ", label=\"l%d\"", edgeIndex);
  }
  fprintf(stream, "]\n");
}

void ehgraph_fprint_dot(FILE* stream, ehgraph_t* graph, char* name) {
  if (NULL == graph || NULL == stream || NULL == name) {
    DB_ERROR("null argument");
    return; 
  }

  /* assure that we have enough bits */
  assert(graph->ptrdim < sizeof(bit_field) * 8 - 1); 

  bit_field* node_vars = (bit_field*) calloc(graph->size, sizeof(bit_field));
  if (NULL == node_vars) {
    DB_ERROR("allocation error");
    return;
  }

  for(int i = 0; i < (int) graph->ptrdim; ++i) {
    enode_t target = VAR2NODE(graph, i);
    bfield_set(node_vars + target, i);
  }

  fprintf(stream, "digraph %s{\n", name);
  fprintf(stream, "  node[shape=box]\n");
  for (int i = 0; i < (int) graph->size; ++i) {
    bit_field f = node_vars[i];
    bool first = true;
    fprintf(stream, "  n%d[label=\"", i);  
    if (NODE_NULL == i) /* null */
      fprintf(stream, "#\\n");
    else if (NODE_DANGLING == i) 
      fprintf(stream, "dangling\\n");
#ifdef EHGRAPH_DEBUG
    else 
      fprintf(stream,"n%d\\n", i);
#endif //EHGRAPH_DEBUG
      

    for (int j = 0; j < (int) graph->ptrdim; ++j)
      if (bfield_get(node_vars[i], j)) {
	if (!first) 
	  fprintf(stream, ", ");
	first = false;

	fprintf(stream, "x%d", j);
      }

    fprintf(stream, "\"]\n");
  }

  /* print the edges */
  for (size_t i = 0; i < graph->size; ++i) {
    enode_info_t* node = &GET_NODE(graph, i);
    bit_field inhibit = 0;
    char label[100];
    label[0] = 0;

     
    if (node->inDoubleLink)
      bfield_set(&inhibit, 1);

    for (size_t j = 0; j < graph->nLinks; ++j)
      if (!bfield_get(inhibit, j) && (node->superEdge[j] || (0 == j && node->outDoubleLink))) {
	bfield_set(&inhibit, j);
	enode_t dest = node->links[j];

	bool pred = false; 
	for (size_t k = 0; k < graph->nLinks && !pred; ++k)
	  pred |= node->isPredicate[j][k];

	if (pred)
	  fprintf(stream, "n%d -> edgen%dn%d -> n%d", i
		  , i, dest, dest);
	else
	  fprintf(stream, "n%d -> n%d", i, dest);
	
	if (0 == j && node->outDoubleLink) {
	  fprintf(stream, "[dir=both, color=\"blue:red\"]\n");
	} else {
	  printEdgeArg(stream, j, true, graph->linkNames);
	}

	bool firstPredicate = true, firstNull = true;
	/* draw the edges for the predicates */
	for (size_t linkIndex = 0; linkIndex < graph->nLinks; ++linkIndex) {
	  if (node->isPredicate[j][linkIndex]) {
	    enode_t target = node->predicateTarget[j][linkIndex];

	    if (firstPredicate) {
	      fprintf(stream, "edgen%dn%d[shape=ellipse,label=\"super edge\"]\n", i, dest);
	      firstPredicate = false;
	    }
	    
	    //bfield_set(&inhibit, linkIndex);
	    if (NODE_NULL == target) {
	      if (firstNull) {
		 /* TO-DO put the variables in the label of null */
		fprintf(stream, "nulln%dn%d[label=\"#\"]\n", i, dest);
		firstNull = false;
	      }
	      fprintf(stream, "edgen%dn%d -> nulln%dn%d", i, dest, i, dest);		      
	    } else {
	      fprintf(stream, "edgen%dn%d -> n%d", i, dest, target);
	    }	    
	    printEdgeArg(stream, linkIndex, false, graph->linkNames);
	  }
	  
	}
      }

    for (size_t j = 0; j < graph->nLinks; ++j)
      if (!bfield_get(inhibit, j)) {
	enode_t dest = node->links[j];
	fprintf(stream, "n%d -> n%d", i, dest);
	printEdgeArg(stream, j, false, graph->linkNames);
      }
  }      
  fprintf(stream, "  {rank=top; n0}\n}\n");
}

void print_with_datadim(FILE* stream, char* prefix, size_t index, char** datadim) {
  assert(stream);
  if (datadim) 
    fprintf(stream, "%s", datadim[index]);
  else
    fprintf(stream, "%s%d", prefix, index);

}

void print_var_name(FILE* stream, size_t index, char** datadim) {
  print_with_datadim(stream, "x", index, datadim);
}


void ehgraph_fprint(FILE * stream, hg_manager_t * man,
        ehgraph_t * graph, char **name_of_dim) {

  do {
    CheckNotNULL(man);
    Man_CheckArgNotNULL(stream, man, FPRINT);
    Man_CheckArgNotNULL(graph, man, FPRINT);

    
    fprintf(stream, "size: %d, ptrdim: %d, nlinks: %d ", graph->size, graph->ptrdim, graph->nLinks);
    if (graph->linkNames) {
      fprintf(stream, "(%s", graph->linkNames[0]);
      
      for (size_t i = 1; i < graph->nLinks; ++i)
	fprintf(stream, ", %s", graph->linkNames[i]);
      fprintf(stream, ")");
    }
    fprintf(stream, "\n");

    if (graph->isTop) {
      fprintf(stream, "Graph is top");
      break;
    }
    if (graph->isBottom) {
      fprintf(stream, "Graph is bottom");
      break;
    }

    for (size_t i = 0; i < graph->size; ++i) {
      fprintf(stream, "Node %d (n%d)\n", i, i);
      bool firstLabel = true;
      if (NODE_NULL == i) {
	fprintf(stream, "\tLabel: NULL, ");
	firstLabel = false;
      }
      if (NODE_DANGLING == i) {
	fprintf(stream, "\tLabel: Dangling");
	firstLabel = false;
      }

      for (size_t j = 0; j < graph->ptrdim; ++j)
	if (VAR2NODE(graph, j) == i) {
	  if (!firstLabel) {
	    fprintf(stream, ", ");
	    firstLabel = false;
	  }
	  print_var_name(stream, j, name_of_dim);
	}

      if (!firstLabel) fprintf(stream, "\n");
      /* print neighbours */
      for (size_t j = 0; j < graph->nLinks; ++j) {
	if (graph->linkNames) 
	  fprintf(stream, "\t.%s", graph->linkNames[j]);
	else 
	  fprintf(stream, "\t.l%d", j);
	fprintf(stream, " -> n%d, type: %s edge", GET_LINK(graph, i, j), GET_NODE(graph, i).superEdge[j] ? "super" : "simple");

	if ( (0 == j && GET_NODE(graph, i).outDoubleLink) || (1 == j && GET_NODE(graph, i).inDoubleLink) )
	  fprintf(stream, ", double linked");
	fprintf(stream, "\n");
	
	if (GET_NODE(graph, i).superEdge[j] && graph->nLinks > 1) {
	  fprintf(stream, "\t\tPredicates: ");
	  bool firstPredicate = true;
	  for (size_t k = 0; k < graph->nLinks; ++k)
	    if (GET_NODE(graph, i).isPredicate[j][k]) {
	      enode_t t = GET_NODE(graph, i).predicateTarget[j][k];
	      if (!firstPredicate) 
		fprintf(stream, ", ");
	      firstPredicate = false;
	      print_with_datadim(stream, "l", k, graph->linkNames);
	      fprintf(stream, "->n%d", t);
	    }
	  fprintf(stream, "\n");
	}
      }
    }

  } while(0);
}

