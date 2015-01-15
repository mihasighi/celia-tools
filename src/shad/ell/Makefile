include ../Makefile.config

PREFIX = ../..

DEBUG = _debug


# Use ICFLAGS to specify machine-independent compilation flags.
ICFLAGS = \
-I$(MPFR_PREFIX)/include \
-I$(GMP_PREFIX)/include \
-I$(CAML_PREFIX)/lib/ocaml \
-I$(CAMLIDL_PREFIX)/lib/ocaml \
-I$(MLGMPIDL_PREFIX)/include \
-I$(APRON_PREFIX)/include \
-I$(APRON_PREFIX)/lib


MANAGER_SOURCES=hg_manager.c hg_manager_test.c sh_result.c hg_abstract0.c

MANAGER_OBJECTS=$(subst .c,$(DEBUG).o,$(MANAGER_SOURCES))

MANAGER_DEPS=$(subst .o,.d,$(MANAGER_OBJECTS))

debug_makefile:
	echo $(CFLAGS_DEBUG)

-include $(MANAGER_DEPS)


EHGRAPH_SOURCES=ehgraph_representation.c ehgraph_test_graphs.c ehgraph_print.c ehgraph_util.c \
		ehgraph_predicate.c ehgraph_closure.c

EHGRAPH_OBJECTS=$(subst .c,$(DEBUG).o,$(EHGRAPH_SOURCES))

EHGRAPH_DEPS=$(subst .o,.d,$(EHGRAPH_OBJECTS))

-include $(EHGRAPH_DEPS)


APRON_LIBS = \
-lapron$(DEBUG) \
-lboxD \
-loctD \
-lpolkaMPQ$(DEBUG) 

manager_test: $(MANAGER_OBJECTS) $(EHGRAPH_OBJECTS)
	$(CC) $(CFLAGS_DEBUG) $(ICFLAGS) -o $@ $(MANAGER_OBJECTS) $(EHGRAPH_OBJECTS) \
	-L$(APRON_PREFIX)/lib $(APRON_LIBS) \
	-L$(MPFR_PREFIX)/lib -lmpfr \
	-L$(GMP_PREFIX)/lib -lgmp \
	-lm

%_debug.o: %.c
	$(CC) $(CFLAGS_DEBUG) $(ICFLAGS) -c -o $@ $<
	$(CC) $(CFLAGS_DEBUG) $(ICFLAGS) -MM -MT $@ $*.c > $*_debug.d

%.o: %.c 
	$(CC) $(CFLAGS) $(ICFLAGS) -c -o $@ $<
	$(CC) $(CFLAGS) $(ICFLAGS) -MM $*.c > $*.d
