CC 		= g++

OPT3 	= -O3
CFLAGS 	= -Wall

LDFLAGS = -lm -lquadmath

all: Def-est-quad

direct_v7-opt: Def-est-quad.cpp
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $(OPT3) -fopenmp $+ $(LDFLAGS) -o $@

clean:
	rm -rf Def-est-quad *.o *.out *.err *.prv *.pcf *.row *.sym

