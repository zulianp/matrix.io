ifeq ($(debug),1)
	CFLAGS += -pedantic -O0 -g
	CXXFLAGS += -O0 -g
	CUFLAGS += -O0 -g
else ifeq ($(prof),1)
	CFLAGS += -pedantic -O2 -g -DNDEBUG
	CXXFLAGS += -O2 -g -DNDEBUG
	CUFLAGS += -O2 -g -DNDEBUG 
else ifeq ($(asan), 1)
	ASAN_FLAGS += -fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g 
	ASAN_FLAGS += -O0
# 	ASAN_FLAGS += -O1
	CXXFLAGS += $(ASAN_FLAGS)
	CFLAGS += -pedantic $(ASAN_FLAGS)
else
	CFLAGS += -Ofast -DNDEBUG
	CXXFLAGS += -Ofast -DNDEBUG
	CUFLAGS += -O3 -DNDEBUG 
endif
ifeq ($(parmetis), 1)
	metis = 1
	CFLAGS += -I$(PARMETIS_DIR)/include -DMATRIX_IO_ENABLE_PARMETIS
	CXXFLAGS += -I$(PARMETIS_DIR)/include -DMATRIX_IO_ENABLE_PARMETIS
	DEPS += -L$(PARMETIS_DIR)/lib -lparmetis

	OBJS += matrixio_parmetis.o
	GOALS += partition_crs
endif

ifeq ($(metis), 1)
	METIS_DIR ?= $(PARMETIS_DIR)/../metis
	GKLIB_DIR ?= $(PARMETIS_DIR)/../gklib

	CFLAGS += -I$(METIS_DIR)/include -DMATRIX_IO_ENABLE_METIS
	CXXFLAGS += -I$(METIS_DIR)/include -DMATRIX_IO_ENABLE_METIS
	
	DEPS += -L$(METIS_DIR)/lib -lmetis
	DEPS += -L$(GKLIB_DIR)/lib -lGKlib
endif

VPATH = graphs:checks
INCLUDES += -Igraphs -Ichecks
DEPS += -lm

CFLAGS += -fPIC
OBJS += matrixio_crs.o utils.o matrixio_array.o array_dtof.o array_ftod.o matrixio_checks.o

GOALS += test print_crs print_array libmatrix.io.a 

MPICXX ?= mpicxx
MPICC ?= mpicc
AR ?= ar

all : $(GOALS)

INCLUDES += -I$(PWD)

libmatrix.io.a : $(OBJS)
	$(AR) r $@ $^ ; \

test : drivers/test.c matrixio_crs.o utils.o matrixio_array.o
	$(MPICC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) ; \

print_crs : drivers/print_crs.c matrixio_crs.o utils.o matrixio_array.o
	$(MPICC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) ; \

partition_crs: drivers/partition_crs.c libmatrix.io.a
	$(MPICC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(DEPS) ; \

print_array : drivers/print_array.c matrixio_crs.o utils.o matrixio_array.o
	$(MPICC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) ; \

matrixio_checks.o : matrixio_checks.cpp
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) $(INTERNAL_CXXFLAGS) -c $<

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
