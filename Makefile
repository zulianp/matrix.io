ifeq ($(debug),1)
	CFLAGS += -pedantic -O0 -g
else ifeq ($(prof),1)
	CFLAGS += -pedantic -O2 -g
else
	CFLAGS += -pedantic -O3 -DNDEBUG
endif

GOALS = test print_crs print_array

CC=mpicc

all : $(GOALS)

test : test.o matrixio_crs.o utils.o matrixio_array.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) ; \

print_crs : print_crs.o matrixio_crs.o utils.o matrixio_array.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) ; \

print_array : print_array.o matrixio_crs.o utils.o matrixio_array.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) ; \

run_test: test
	mpiexec -np 1 ./test
	diff data/test/dump.raw data/test/rhs.raw
	diff data/test/dump.colindex.raw data/test/lhs.colindex.raw
	diff data/test/dump.value.raw data/test/lhs.value.raw
	diff data/test/dump.rowindex.raw data/test/lhs.rowindex.raw

%.o : %.c
	$(CC) $(CFLAGS) -c $<

.SUFFIXES :
.PRECIOUS :

clean:
	rm *.o $(GOALS) 

.SUFFIXES:

.PHONY: clean all
