all: clean lib exam

lib:
	@(cd ../src; make)

include ../Makefile.inc

exam.o: exam.cxx ../src/lib${LSSPLIB}${LIB_SUFFIX}

clean:
	@rm -f *.o *.log exam
