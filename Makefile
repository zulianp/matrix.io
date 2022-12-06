ifeq ($(debug),1)
	CFLAGS += -pedantic -O0 -g
else ifeq ($(prof),1)
	CFLAGS += -pedantic -O2 -g
else
	CFLAGS += -pedantic -O3 -DNDEBUG
endif

GOALS = test print_crs

CC=mpicc

all : $(GOALS)

test : test.o matrixio_crs.o utils.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) ; \

print_crs : print_crs.o matrixio_crs.o utils.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) ; \

%.o : %.c
	$(CC) $(CFLAGS) -c $<

.SUFFIXES :
.PRECIOUS :

clean:
	rm *.o $(GOALS) 

.SUFFIXES:

.PHONY: clean all
