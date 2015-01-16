#TODO
include $(CELIA)/Makefile.config
#HAS_SHARED=1

#---------------------------------------
# Directories
#---------------------------------------

SRCDIR = $(shell pwd)
#
SHAD_PREFIX = ../shad
#
# C include and lib directories
INCDIR = $(PREFIX)/include
LIBDIR = $(PREFIX)/lib
BINDIR = $(PREFIX)/bin
#
GLIB_IFLAGS=$(shell pkg-config --cflags glib-2.0)
GLIB_LDFLAGS=$(shell pkg-config --libs glib-2.0)

#---------------------------------------
# C part
#---------------------------------------
ICFLAGS = $(GLIB_IFLAGS) \
-I$(SHAD_PREFIX)

LDFLAGS = $(GLIB_LDFLAGS) \
-L$(SHAD_PREFIX) -lshad_debug \
-lm

#---------------------------------------
# Files
#---------------------------------------

H_FILES = \
noll.h \
noll_fun.h \
noll_internal.h \

C_FILES = \
noll_build.c \
noll_closure.c \
noll_nary.c \
noll_query.c \
noll_print.c \
noll_resize.c \
noll_transfer.c \
noll_transform.c \
noll_widen.c  


#---------------------------------------
# Rules
#---------------------------------------

all: libnoll.a libnoll_debug.a
ifneq ($(HAS_SHARED),)
all: libnoll.so libnoll_debug.so
endif

#---------------------------------------
# Misc rules
#---------------------------------------

clean:
	/bin/rm -f *~ *.aux *.bbl *.blg *.dvi *.log *.toc *.ps *.pdf 
	/bin/rm -f *.pg *.fn *.ky *.cp *.tp *.vr *.kys *.info* 
	/bin/rm -f *.o *.a *.cmi *.cmo *.cmx *.cmxa *.cma *.so
	/bin/rm -f *_test

install: all
	mkdir -p $(SHAD_PREFIX)/include
	cp $(H_FILES) $(SHAD_PREFIX)/include
	mkdir -p $(SHAD_PREFIX)/lib
	for i in libnoll.a libnoll_debug.a libnoll.so libnoll_debug.so; do if test -f $$i; then $(INSTALL) $$i $(SHAD_PREFIX)/lib; fi; done

distclean:
	/bin/rm -f $(H_FILES:%=$(SHAD_PREFIX)/include/%)
	/bin/rm -f $(SHAD_PREFIX)/lib/libnoll.a
	/bin/rm -f $(SHAD_PREFIX)/lib/libnoll_debug.a
	/bin/rm -f $(SHAD_PREFIX)/lib/libnoll.so
	/bin/rm -f $(SHAD_PREFIX)/lib/libnoll_debug.so
	/bin/rm -f noll.pdf

#---------------------------------------
# Latex rules
#---------------------------------------

.PHONY : html

noll.pdf: noll.texi 
	$(TEXI2DVI) --pdf -o $@ $<

noll.info: noll.texi 
	$(MAKEINFO) -o $@ $<

html: noll.texi 
	$(TEXI2HTML) -split=section -nonumber -menu -subdir=html $<

#---------------------------------------
# C rules
#---------------------------------------

libnoll.a: $(C_FILES:%.c=%.o) 
	$(AR) rs $@ $(C_FILES:%.c=%.o) 
	$(RANLIB) $@

libnoll_debug.a: $(C_FILES:%.c=%_debug.o) 
	$(AR) rs $@ $(C_FILES:%.c=%_debug.o) 
	$(RANLIB) $@

libnoll.so: $(C_FILES:%.c=%.o) 
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -o $@ $(C_FILES:%.c=%.o) 

libnoll_debug.so: $(C_FILES:%.c=%_debug.o) 
	$(CC) $(CFLAGS_DEBUG) $(LDFLAGS) -shared -o $@ $(C_FILES:%.c=%.o) 

%_test: %_test.c
	$(CC) $(CFLAGS_DEBUG) $(ICFLAGS) -o $@ $< \
	-L../noll -lnoll_debug $(LDFLAGS)

.PHONY : dep

dep: $(H_FILES)
	for i in $(C_FILES); do		\
		cpp $(ICFLAGS) -MM $$i;	\
	done

#--------------------------------------------------------------
# IMPLICIT RULES AND DEPENDENCIES
#--------------------------------------------------------------

.SUFFIXES: .c .h .o

#---------------------------------------
# C generic rules
#---------------------------------------

%.o: %.c 
	$(CC) $(CFLAGS) $(ICFLAGS)  -c $<

%_debug.o: %.c 
	$(CC) $(CFLAGS_DEBUG) $(ICFLAGS)  -c -o $@ $<


#---------------------------------------
# dependencies (generated with make dep)
#---------------------------------------

