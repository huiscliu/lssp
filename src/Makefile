
.SUFFIXES: .cxx

default: lib
all: lib

include ../Makefile.inc
include Makefile.dep

CPPFLAGS     := ${CPPFLAGS} ${BUILD_CPPFLAGS} ${USER_CPPFLAGS}
CFLAGS       := ${CFLAGS} ${BUILD_CFLAGS} ${USER_CFLAGS}
CXXFLAGS     := ${CXXFLAGS} ${BUILD_CXXFLAGS} ${USER_CXXFLAGS}

SRC_CXX   = $(wildcard *.cxx)
DIR_CXX   = $(notdir $(SRC_CXX))
OBJS_CXX  = $(patsubst %.cxx,%.o,$(DIR_CXX))

OBJS     := $(OBJS_CXX)

lib: liblssp$(LIB_SUFFIX)

liblssp.a: ${OBJS}
	${AR} -ru $@ $(OBJS)
	${RANLIB} $@

liblssp.so: ${OBJS}
	${BUILD_SHARED} ${OBJS} ${BUILD_SHARED_LIBS} -o liblssp${LIB_SUFFIX}
.PHONY: clean

clean:
	@rm -fv liblssp.a liblssp.so
	@rm -fv *.o 
