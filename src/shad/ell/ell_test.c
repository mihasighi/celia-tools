#include "ehgraph_util.h"
#include "ehgraph_internal.h"
#include "ehgraph_predicate.h"
#include "ehgraph_fun.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

FILE* out;
FILE* flog;
FILE* fdbg;

void gPrintDot(ehgraph_t* graph, const char* fileName) {
  assert(NULL != graph && NULL != fileName);
  if (NULL == graph || NULL == fileName) {
    DB_ERROR("null argument");
    return;
  }

  FILE* fout = fopen(fileName, "w");
  ehgraph_fprint_dot(fout, graph, "test");
  fclose(fout);
}

void print_bfield(FILE* fout, bit_field field,size_t limit) {
  fprintf(fout, "{");
  bool first = true;
  for (size_t i = 0; i < limit; ++i)
    if (bfield_get(field, i)) {
      if (!first)
	fprintf(fout, ", ");
      first = false;
      fprintf(fout, "%d", i);
    }
  fprintf(fout, "}");
}

void test_bfield(void) {
  bit_field b = 0;
  fprintf(out, "Empty set:\n");
  print_bfield(out, b, 10);
  fprintf(out, "\n");
  bfield_set(&b, 3);
  bfield_set(&b, 5);
  bfield_set(&b, 7);
  bfield_set(&b, 9);
  fprintf(out, "Second set:\n");
  print_bfield(out, b, 10);
  fprintf(out, "\n");
  bfield_unset(&b, 3);
  bfield_unset(&b, 7);
  fprintf(out, "Third set:\n");
  print_bfield(out, b, 10);
  fprintf(out, "\n");
}



void test_util(void) {
  bool mat[2][2] = {{true, false}, {false, false}};
  assert(bvect_exists(mat[0], 2));
  assert(!bvect_exists(mat[1], 2));
  printf("test_util finished execution");
}

void test_print_i(int index) {
  MAKE_STRING(file, "test_print_%d.dot", index);
  gPrintDot(create_graph(index), file);
}

void test_print(void) {
  for (int i = 1; i <= g_nGraphs; ++i) 
    test_print_i(i);
}



void test_eq_random_i(int index) {
  ehgraph_t* graph = NULL, *graph2 = NULL;
  int* perm = NULL;

  do {
    graph = create_graph(index);
    assert(NULL != graph);

    size_t size = graph->size, nVars = graph->ptrdim;

    CheckedMalloc(perm, int, size);
    for (enode_t i = 0; i < size; ++i)
      perm[i] = i;

    random_permutation(perm + 1, size - 1);
    d_printf("Testing the closure of graph %d\n", index);
    d_printf("The random permutation is:\n");
    for (size_t i = 0; i < size; ++i)
      d_printf("%d ", perm[i]);
    d_printf("\n");
    
    graph2 = permute_graph(graph, perm);
    if (NULL == graph2) {
      DB_ERROR("failed to permute graph");
      break;
    }
    
    GET_LINK(graph2, 1, 2) = 6;
    MAKE_STRING(file, "permuted_%d.dot", index);
    gPrintDot(graph2, file);

    bool result;
    d_printf("\nResult ... %s\n", ehgraph_is_eq(graph, graph2, &result) ? "Ok" : "Failed");

  } while(0);
}

void test_eq(void) {
  ehgraph_t *top1 = ehgraph_top(), *top2 = ehgraph_top();
  ehgraph_t *bot1 = ehgraph_bottom(), *bot2 = ehgraph_bottom();

  bool result;
  assert(ehgraph_is_eq(top1, top2, &result));
  assert(result);
  assert(ehgraph_is_eq(top1, bot1, &result));
  assert(!result);
  assert(ehgraph_is_eq(bot1, bot2, &result));
  assert(result);

  ehgraph_t* g9 = create_graph(9), *g10 = create_graph(10);
  assert(ehgraph_is_eq(g9, g10, &result));
  assert(result);
}


void test_eq_random(void) {
  test_eq_random_i(4);
}


void test_closure_i(int index) {
  ehgraph_t* closed = ehgraph_close(create_graph(index));
  if (NULL == closed)
    return;

  MAKE_STRING(file, "closed_%d.dot", index);
  gPrintDot(closed, file);
}

void test_closure(void) {
  test_closure_i(4);
  test_closure_i(5);
  test_closure_i(6);
}

void test_split_i(int index, bool splitOneOne) {
  MAKE_STRING(fileIn, "split_%d_in.dot", index);
  ehgraph_t* graph = create_graph(index);
  gPrintDot(graph, fileIn);



  enode_t node;
  ehgraph_t* spl = ehgraph_split(graph, 0, 0, &node, splitOneOne);
  MAKE_STRING(fileOut, "split_%d_out.dot", index);
  printf("Split node: %d\n", node);
  gPrintDot(spl, fileOut);
}

void test_split(void) {
  test_split_i(7, true);
  test_split_i(8, true);
}

void test_leq(void) {
  bool result = false;
  ehgraph_t* top1 = ehgraph_top(), *top2 = ehgraph_top();
  ehgraph_t* bot1 = ehgraph_bottom(), *bot2 = ehgraph_bottom();
  assert(top1 && top2 && bot1 && bot2);
  assert(!ehgraph_is_leq(NULL, top1, &result));
  assert(!ehgraph_is_leq(top1, NULL, &result));
  assert(!ehgraph_is_leq(top1, bot1, NULL));

  assert(ehgraph_is_leq(top1, bot1, &result));
  assert(!result);
  assert(ehgraph_is_leq(top1, top2, &result));
  assert(result);
  assert(ehgraph_is_leq(bot1, top2, &result));
  assert(!result);
  assert(ehgraph_is_leq(bot1, bot2, &result));
  assert(result);
  
  ehgraph_t* g11 = create_graph(11), *g12 = create_graph(12);
  assert(g11 && g12);
  assert(ehgraph_is_eq(g11, g12, &result));
  assert(!result);
  assert(ehgraph_is_leq(g11, g12, &result));
  assert(result);
}

void init(void) {
  out = stdout;
  fdbg = fopen("debug.txt", "w");
  flog = fopen("log.txt", "w");  
}

void finalize(void) {
  fclose(fdbg);
  fclose(flog);
}


int main(void) {
  init();
  out = stdout;

  //test_bfield();
  //test_util();
  //test_print();
  //test_closure();
  //test_eq_random();
  //test_print_i(6);
  //test_split();
  //test_eq();
  test_leq();

  finalize();
  return 0;
}
