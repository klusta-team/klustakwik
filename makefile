# MAKE FILE

PROGRAM = KlustaKwik
OBJS = io.o linalg.o log.o parameters.o precomputations.o util.o memorytracking.o klustakwik.o
CC = g++
DEBUG = -g
PROFILE = -pg
OPTIMISATIONS = -O3 -ffast-math -march=native
CFLAGS = -Wall -c -Wno-write-strings $(OPTIMISATIONS) -fopenmp
LFLAGS = -Wall -fopenmp

all: executable

# Adds debug flags
debug: CFLAGS += $(DEBUG)
debug: LFLAGS += $(DEBUG)
debug: executable

# Adds profiling flags
profile: CFLAGS += $(PROFILE)
profile: LFLAGS += $(PROFILE)
profile: executable

# Remove -march=native from optimisations for older gcc versions
oldgcc: OPTIMISATIONS = -O3 -ffast-math
oldgcc: executable

.PHONY: all debug native executable clean

executable: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(PROGRAM)

clean:
	\rm *.o $(PROGRAM)

######## DEPENDENCIES FOR FILES IN PROJECT ###################################

linalg.h: numerics.h
klustakwik.h: parameters.h log.h util.h numerics.h linalg.h memorytracking.h
numerics.h: globalswitches.h
parameters.h: numerics.h
util.h: numerics.h
memorytracking.h: numerics.h

io.o: io.cpp klustakwik.h
	$(CC) $(CFLAGS) $<

linalg.o: linalg.cpp linalg.h
	$(CC) $(CFLAGS) $<

log.o: log.cpp log.h parameters.h
	$(CC) $(CFLAGS) $<

klustakwik.o: klustakwik.cpp klustakwik.h
	$(CC) $(CFLAGS) $<
	
parameters.o: parameters.cpp parameters.h log.h util.h
	$(CC) $(CFLAGS) $<

precomputations.o: precomputations.cpp klustakwik.h
	$(CC) $(CFLAGS) $<

util.o: util.cpp util.h
	$(CC) $(CFLAGS) $<

memorytracking.o: memorytracking.cpp memorytracking.h
	$(CC) $(CFLAGS) $<
