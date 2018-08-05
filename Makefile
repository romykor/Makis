CC 		= g++

OPT3 	= -O3
CFLAGS 	= -Wall

LDFLAGS = -lm -lquadmath

all: Eotvos_prepdata-quad

direct_v7-opt: Eotvos_prepdata-quad.cpp
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $(OPT3) -fopenmp $+ $(LDFLAGS) -o $@

clean:
	rm -rf Eotvos_prepdata-quad *.o *.out *.err *.prv *.pcf *.row *.sym

