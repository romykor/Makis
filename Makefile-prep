CC 		= g++

OPT3 	= -O3
CFLAGS 	= -Wall

LDFLAGS = -lm -lquadmath

all: Def_prep-quad

direct_v7-opt: Def_prep-quad.cpp
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $(OPT3) -fopenmp $+ $(LDFLAGS) -o $@

clean:
	rm -rf Def_prep-quad *.o *.out *.err *.prv *.pcf *.row *.sym

