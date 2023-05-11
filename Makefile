ifeq ($(debug),1)
	CFLAGS += -pedantic -O0 -g
else ifeq ($(prof),1)
	CFLAGS += -pedantic -O2 -g
else
	CFLAGS += -pedantic -O3 -DNDEBUG
endif

CFLAGS += -fPIC

GOALS = test print_crs print_array libmatrix.io.a

MPICC ?= mpicc
AR ?= ar

all : $(GOALS)

INCLUDES += -I$(PWD)

libmatrix.io.a : matrixio_crs.o utils.o matrixio_array.o array_dtof.o array_ftod.o
	$(AR) r $@ $^ ; \

test : drivers/test.c matrixio_crs.o utils.o matrixio_array.o
	$(MPICC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) ; \

print_crs : drivers/print_crs.c matrixio_crs.o utils.o matrixio_array.o
	$(MPICC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) ; \

print_array : drivers/print_array.c matrixio_crs.o utils.o matrixio_array.o
	$(MPICC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) ; \

run_test: test
	rm  data/test/dump.*.raw
	mpiexec -np 8 ./test
	diff data/test/dump.raw data/test/rhs.raw
	diff data/test/dump.colindex.raw data/test/lhs.colindex.raw
	diff data/test/dump.value.raw data/test/lhs.value.raw
	diff data/test/dump.rowindex.raw data/test/lhs.rowindex.raw

print_test: run_test
	od -L   data/test/lhs.rowindex.raw 
	od -L   data/test/dump.rowindex.raw 
	od -D   data/test/lhs.colindex.raw 
	od -D   data/test/dump.colindex.raw 
	od -f   data/test/lhs.value.raw 
	od -f   data/test/dump.value.raw 
	od -f   data/test/rhs.raw 
	od -f   data/test/dump.raw 

%.o : %.c
	$(MPICC) $(CFLAGS) $(INCLUDES) -c $<

.SUFFIXES :
.PRECIOUS :

clean:
	rm *.o *.a $(GOALS) 

.SUFFIXES:

.PHONY: clean all
