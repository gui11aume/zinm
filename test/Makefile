vpath %.c ..
vpath %.c lib
vpath %.h ..
vpath %.h lib

all: testRpackage test

P= runtests

CC= gcc
SOURCES= runtests.c
IMPLICIT_SOURCES= zinm.c
OBJECTS= libunittest.so
INCLUDES= -I.. -Ilib
COVERAGE= -fprofile-arcs -ftest-coverage
CFLAGS= -std=gnu99 -g -Wall -O0 $(INCLUDES) $(COVERAGE)
LDLIBS= -L`pwd` -Wl,-rpath=`pwd` -lm -lunittest

$(P): $(OBJECTS) $(SOURCES) $(IMPLICIT_SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) $(OBJECTS) $(LDLIBS) -o $@

clean:
	rm -f $(P) $(OBJECTS) *.gcda *.gcno *.gcov gmon.out .inspect.gdb

libunittest.so: unittest.c
	$(CC) -fPIC -shared $(CFLAGS) -o libunittest.so lib/unittest.c

../Rpackage/src/Rzinm.so:
	cd ../Rpackage/src && $(MAKE)

test: $(P)
	./$(P)

testRpackage: ../Rpackage/src/Rzinm.so
	R -f testset.R

inspect: $(P)
	gdb --command=.inspect.gdb --args $(P)

valgrind: $(P)
	#valgrind --vgdb=yes --vgdb-error=0 ./$(P) --debug
	valgrind --leak-check=full ./$(P) --debug
