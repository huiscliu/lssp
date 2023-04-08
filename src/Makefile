
.SUFFIXES: .cxx

default: lib
all: lib

include ../Makefile.inc
include Makefile.dep

SRC_CXX   = $(wildcard *.cxx)
DIR_CXX   = $(notdir $(SRC_CXX))
OBJS_CXX  = $(patsubst %.cxx,%.o,$(DIR_CXX))

OBJS     := $(OBJS_CXX)

lib: liblssp.a

liblssp.a: ${OBJS}
	${AR} -ru $@ $(OBJS)
	${RANLIB} $@

.PHONY: clean

clean:
	@rm -fv liblssp.a
	@rm -fv *.o 
