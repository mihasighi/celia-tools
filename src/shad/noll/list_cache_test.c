/*
 * list_cache_test.c:
 * Test file encoding the cache data structure from Lee & al. CAV'11.
 */

#include <assert.h>
#include "sh_manager.h"
#include "sh_utils.h"
#include "noll.h"
#include "noll_fun.h"

/*
 * The cache list is a data structure of type struct_C built from
 * - a doubly linked list of type struct_L
 *   which is kept circular (and sorted after an integer field)
 * - a cell of type struct_L*
 *   which maintains the last added cell
 *
 * The functions tested on this data structure are:
 * - add_L(X,Y):
 * 		adds cell X before Y
 * - cut_L(X)
 *      remove X from its list
 * - find_L(Q,K,T)
 *      find the first cell T having the key > k
 * - add_C(K,C)
 *      adds k in the list and keep in cache the pointer to the cell added
 * - find_C(K,C)
 *      look inside C to find the key and start by the cache
 * - create_C()
 *      returns a cached list by successively calling add_C
 * - remove_C(K,C)
 *      remove the key from the list and the cache
 * - main
 *      calls create and then remove some key
 */

/* ********************************************************************** */
/* I. Globals */
/* ********************************************************************** */

/* Manager of the program + abstract domain */
static sh_manager_t* man;

/* source file analyzed */
static char* filename = "07-cache.c";

/* quick recall of ids */
typedef enum {
	ID_FILE = 0,
	ID_INT,
	ID_STRUCT_L,
	ID_STRUCT_L_REF,
	ID_STRUCT_C,
	ID_STRUCT_C_REF,
	ID_FLD_NEXT,
	ID_FLD_PREV,
	ID_FLD_KEY,
	ID_FLD_CACHED,
	ID_FLD_Q,
	ID_UNKNOWN
} id_e;
static size_t ids[ID_UNKNOWN];

/* ICFG informations  and computed annotations */
GPtrArray* icfg = NULL; /* array of sh_icfg_annot_t* indexed by vertices */

/* ********************************************************************** */
/* II. Types */
/* ********************************************************************** */

void register_struct_L(void) {

	/* register
	 * type struct L, struct L*, f
	 * fields next, prev, key
	 * predicate dll_L
	 */
	sh_type_register_dll(man, filename, 3, "L", "next", "prev", "key");

	/* recover ids */
	sh_typenv_t* tenv = sh_manager_get_typenv(man);
	ids[ID_INT] = sh_typenv_lookup_type(tenv, "int");
	ids[ID_STRUCT_L] = sh_typenv_lookup_type_record(tenv, "L");
	int fld = sh_typenv_lookup_field(tenv, "next", "struct_L");
	assert (fld >= 0);
	ids[ID_FLD_NEXT] = (size_t) fld;
	sh_fld_t* fld_t = sh_typenv_get_field(tenv, ids[ID_FLD_NEXT]);
	assert (NULL != fld_t);
	ids[ID_STRUCT_L_REF] = fld_t->tid;
	fld = sh_typenv_lookup_field(tenv, "prev", "struct_L");
	ids[ID_FLD_PREV] = (size_t) fld;
	fld = sh_typenv_lookup_field(tenv, "key", "struct_L");
	ids[ID_FLD_KEY] = (size_t) fld;
}

void register_struct_C(void) {

	sh_filenv_t* fenv = sh_manager_get_filenv(man);
	guint file = sh_filenv_push_fname(fenv, filename);
	ids[ID_FILE] = file;

	sh_typenv_t* tenv = sh_manager_get_typenv(man);
	/* type for record = the struct C */
	gulong line = 11;
	sh_loc_t* loc_str = sh_loc_alloc(file, line, 0);
	size_t tid_str = sh_typenv_add_type_record(tenv, loc_str, "C");
	ids[ID_STRUCT_C] = tid_str;

	/* type for recursive field = struct C* */
	sh_loc_t* loc_ref = sh_loc_alloc(file, line + 4, 0);
	size_t tid_ref = sh_typenv_add_type_ptr0(tenv, loc_str, tid_str, 1);
	ids[ID_STRUCT_C_REF] = tid_ref;

	/* register field cache */
	sh_loc_t* loc_fld_cached = sh_loc_alloc(file, line + 1, 0);
	size_t fid_cached = sh_typenv_add_field0(tenv, loc_fld_cached, "cached",
			tid_str, ids[ID_STRUCT_L_REF]);
	sh_typenv_set_field_recursive(tenv, fid_cached, FALSE);
	ids[ID_FLD_CACHED] = fid_cached;

	/* register field q */
	sh_loc_t* loc_fld_q = sh_loc_alloc(file, line + 2, 0);
	size_t fid_q = sh_typenv_add_field0(tenv, loc_fld_q, "q", tid_str,
			ids[ID_STRUCT_L_REF]);
	sh_typenv_set_field_recursive(tenv, fid_q, FALSE);
	ids[ID_FLD_Q] = fid_q;
}

void register_typ(void) {
	register_struct_L();
	register_struct_C();
}

/* ********************************************************************** */
/* III. Functions */
/* ********************************************************************** */

typedef enum {
	FUN_ADD_C = 0,
	FUN_ADD_L,
	FUN_CREATE_C,
	FUN_CUT_L,
	FUN_FIND_C,
	FUN_FIND_L,
	FUN_MAIN,
	FUN_REMOVE_C,
	FUN_UNKNOWN
} proc_e;

static const char* proc_name[FUN_UNKNOWN] = { "add_c", "add_l", "create_C",
		"cut_L", "find_C", "find_L", "main", "remove_C" };

static guint proc_start[FUN_UNKNOWN];
static guint proc_exit[FUN_UNKNOWN];

void register_add_L(void) {
	/* add(X,Y): add cell X before Y
	 * tmp = (Y)->P;
	 * tmp->N = (X);
	 * (X)->N=Y;
	 * (X)->P=tmp; // (Y)->P
	 * (Y)->P=(X)
	 */
	return; // TODO
}
void register_cut_L(void) {
	return; // TODO
}
void register_find_L(void) {
	return; // TODO
}

void register_add_C(void) {
	return; // TODO
}
void register_find_C(void) {
	return; // TODO
}
void register_remove_C(void) {
	return; // TODO
}
void register_create_C(void) {
	/*
	 * locals: cache, q, tmp
	 * 59: cache=new C();
	 * 60: cache->cached=NULL;
	 * 61: q=new L();
	 * 62: q->key=-1;
	 * 63: q->prev=q
	 * 64: q->next=q;
	 * 65: cache->q=q;
	 * 66: while (random() != cst)
	 * 67:   { tmp=new L();
	 * 68:     tmp->key=-1;
	 * 69:     add_L(tmp,q);
	 * 70:     cache->cached=tmp;
	 * 71:   }
	 * 72: return cache
	 * 73: (project all locals except return)
	 */
	// locations
	gulong start = 59;
	gulong exit = 80;
	sh_loc_t* loc_start = sh_loc_alloc(ids[ID_FILE], start, 0);
	sh_loc_t* loc_exit = sh_loc_alloc(ids[ID_FILE], exit, 0);
	// frame
	sh_frame_t* frame = sh_frame_new(loc_start, loc_exit);
	size_t fid = sh_stack_push_frame(sh_manager_get_varenv(man), frame);

	// variables
	// cache : struct C*
	sh_var_t* v_cache = sh_var_new("cache", fid, ids[ID_STRUCT_C_REF]);
	sh_frame_push(frame,v_cache);
	// q : struct L*
	sh_var_t* v_q = sh_var_new("q", fid, ids[ID_STRUCT_L_REF]);
	sh_frame_push(frame,v_q);
	// tmp : struct L*
	sh_var_t* v_tmp = sh_var_new("tmp", fid, ids[ID_STRUCT_L_REF]);
	sh_frame_push(frame,v_tmp);
	// return : struct C*
	sh_var_t* v_ret = sh_var_new("return", fid, ids[ID_STRUCT_C_REF]);
	sh_frame_push(frame,v_ret);
	// only after get the identifiers of variables
	size_t vid_cache = sh_frame_get_var_id(frame, "cache");
	assert (SH_VAR_NULL != vid_cache);
	size_t vid_q = sh_frame_get_var_id(frame, "q");
	assert (SH_VAR_NULL != vid_q);
	size_t vid_tmp = sh_frame_get_var_id(frame, "tmp");
	assert (SH_VAR_NULL != vid_tmp);
	size_t vid_ret = sh_frame_get_var_id(frame, "return");
	assert (SH_VAR_NULL != vid_ret);

	// for computation of infos
	noll_val_t *pre;
	sh_stmt_t* in_stmt;
	// icfg for this procedure
	// line 59: cache=new C();
	guint vrtx = icfg->len; // start vertex
	sh_icfg_info_t* info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = loc_start;
	info->val = pre = noll_empty(man, fid);
	info->edges = g_ptr_array_new();
	sh_icfg_edge_t* edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx + 1;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_new(man, loc_start, vid_cache);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);
	proc_start[FUN_CREATE_C] = vrtx;

	// line 60: cache->cached=NULL;
	start++;
	vrtx = icfg->len;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	// compute value
	info->val = pre = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
	// otherwise info->val = noll_bottom(man, fid);
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx + 1;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_assign_x_f_y(man, info->loc, vid_cache,
			ids[ID_FLD_CACHED], SH_VAR_NULL);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);

	// line 61: q=new L();
	start++;
	vrtx = icfg->len;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	// compute value
	info->val = pre = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx + 1;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_new(man, info->loc, vid_q);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);

	// line 62: q->key=-1;
	// ignored

	// line 63: q->prev=q
	start++;
	vrtx = icfg->len;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	// compute value
	info->val = pre = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx + 1;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_assign_x_f_y(man, info->loc, vid_q,
			ids[ID_FLD_PREV], vid_q);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);

	// line 64: q->next=q;
	start++;
	vrtx = icfg->len;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	// compute value
	info->val = pre = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx + 1;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_assign_x_f_y(man, info->loc, vid_q,
			ids[ID_FLD_NEXT], vid_q);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);

	// line 65: cache->q=q;
	start++;
	vrtx = icfg->len;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	// compute value
	info->val = pre = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx + 1;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_assign_x_f_y(man, info->loc, vid_cache,
			ids[ID_FLD_Q], vid_q);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);

	// line 66: while (random() != cst)
	start++;
	size_t vrtx_loop = vrtx = icfg->len;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	// compute value
	info->val = pre = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx + 1;
	edge->isprj = FALSE;
	edge->stmt = sh_stmt_new_skip(man, info->loc);
	g_ptr_array_add(info->edges, edge);
	// the exit from loop added afterwards
	g_ptr_array_add(icfg, info);

	// line 67: tmp=new L();
	start++;
	vrtx = icfg->len;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	info->val = pre; // in_stmt = skip
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx + 1;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_new(man, info->loc, vid_tmp);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);

	// line 68:     tmp->key=-1;
	// ignored

	// line 69:     add_L(tmp,q): q->prev->next = tmp
	start++;
	vrtx = icfg->len;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	// compute value
	info->val = pre = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx + 1;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_assign_x_f_y(man, info->loc, vid_q,
			ids[ID_FLD_PREV], vid_tmp);
	g_array_append_val(edge->stmt->info.binary.offset_l,ids[ID_FLD_NEXT]);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);

	// line 69:     add_L(tmp,q): tmp->next = q
	start++;
	vrtx = icfg->len;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	// compute value
	info->val = pre = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx + 1;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_assign_x_f_y(man, info->loc, vid_tmp,
			ids[ID_FLD_NEXT], vid_q);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);

	// line 69:     add_L(tmp,q): tmp->prev = q->prev
	start++;
	vrtx = icfg->len;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	// compute value
	info->val = pre = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx + 1;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_assign_x_f_y(man, info->loc, vid_tmp,
			ids[ID_FLD_PREV], vid_q);
	edge->stmt->info.binary.offset_r
			= g_array_new(FALSE, FALSE, sizeof(size_t));
	g_array_append_val(edge->stmt->info.binary.offset_r,ids[ID_FLD_PREV]);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);

	// line 69:     add_L(tmp,q): q->prev = tmp
	start++;
	vrtx = icfg->len;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	// compute value
	info->val = pre = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx + 1;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_assign_x_f_y(man, info->loc, vid_q,
			ids[ID_FLD_PREV], vid_tmp);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);

	// line 70:     cache->cached=tmp;
	start++;
	vrtx = icfg->len;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	// compute value
	info->val = pre = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx_loop;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_assign_x_f_y(man, info->loc, vid_cache,
			ids[ID_FLD_CACHED], vid_tmp);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);

	// analysis of the loop
	guint it = 0;
	// value at the end of the loop
	pre = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
	info = g_ptr_array_index(icfg, vrtx_loop);
	// join with the value computed at the start of the loop
	noll_val_t* post = noll_join(man, FALSE, info->val, pre); // post1(loop) U init(loop)
	// widen old value (info->val) with the new (join) value
	pre = noll_widening(man, info->val, post);
	// compare to see if fixpoint reached
	while (noll_is_le(man, pre, info->val) == FALSE) {
		it++;
		printf("Iteration %d: start with ", it);
		noll_fprint(stdout, man, pre);
		// redo the loop
		do {
			// update value at the vrtx
			info->val = pre;
			edge = g_ptr_array_index(info->edges,0);
			in_stmt = edge->stmt;
			vrtx = edge->to;
			post = noll_assign_post(man, FALSE, pre, in_stmt, NULL);
			info = g_ptr_array_index(icfg, vrtx);
			pre = noll_join(man, TRUE, info->val, post);
		} while (vrtx != vrtx_loop);
		// new computed value for vrtx_loop
		post = pre;
		pre = noll_widening(man, info->val, post);
	}
	// now we are ready to add the edge from vertex_loop to end of loop
	size_t vrtx_loop_end = icfg->len;
	// add the edge in vrtx_loop
	info = g_ptr_array_index(icfg,vrtx_loop);
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx_loop_end;
	edge->isprj = FALSE;
	edge->stmt = sh_stmt_new_skip(man, info->loc);
	g_ptr_array_add(info->edges, edge);

	// 71: return cached;
	// the edge from the vrtx_loop_end
	start++;
	info = g_new (sh_icfg_info_t,1);
	info->frame = fid;
	info->loc = sh_loc_alloc(ids[ID_FILE], start, 0);
	// compute value
	info->val = pre; // the last widen value in the loop
	info->edges = g_ptr_array_new();
	edge = g_new(sh_icfg_edge_t,1);
	edge->from_call = 0;
	edge->to = vrtx_loop_end + 1;
	edge->isprj = FALSE;
	edge->stmt = in_stmt = sh_stmt_new_assign_x_y(man, info->loc, vid_ret,
			vid_cache);
	g_ptr_array_add(info->edges, edge);
	g_ptr_array_add(icfg, info);

	post = noll_assign_post(man, FALSE, pre, in_stmt, NULL);

	printf("Function create_C ends with ");
	noll_fprint(stdout, man, post);

	// project on local variables
	GArray* locals = g_array_new(FALSE, FALSE, sizeof(size_t));
	g_array_append_val(locals, vid_cache);
	g_array_append_val(locals, vid_q);
	g_array_append_val(locals, vid_tmp);
	post = noll_project_dimensions(man, TRUE, post, *locals);
	printf("Function create_C returns ");
	noll_fprint(stdout, man, post);

	return;
}

void register_main(void) {
	return; // TODO
}

void register_fun(void) {

	icfg = g_ptr_array_new();

	register_add_L();
	register_cut_L();
	register_find_L();

	register_add_C();
	register_find_C();
	register_remove_C();
	register_create_C();

	register_main();
}

/* ********************************************************************** */
/* IV. Analysis */
/* ********************************************************************** */

void compute_fixpoint(char* fname) {
	// the computation is done in the function selected
	return; // TODO
}

/* ********************************************************************** */
/* IV. Test function */
/* ********************************************************************** */

int main(int argv, char** args) {

	/* Test correct initialization of static data */
	assert (NULL != filename);

	man = sh_manager_alloc("shad", "0.1", NULL, NULL);

	/* register file */
	sh_file_register(man, filename);

	/* register types and their predicates */
	register_typ();

	/* print situation for types and preds */
	sh_manager_fdump(stdout, man);

	/* register icfg */
	register_fun();

	/* apply the fixpoint algorithm on icfg */
	if (argv >= 2)
		compute_fixpoint(args[1]); // function name
	else
		compute_fixpoint(NULL); // default

}
