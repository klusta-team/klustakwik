# MAKE FILE, use 'make native' if you are running on a reasonably recent version
# of g++ for best optimisation. When ready for release, change PROGRAM to
# KlustaKwik for backwards compatibility

PROGRAM = KlustaKwik
OBJS = io.o linalg.o log.o parameters.o precomputations.o util.o klustakwik.o
CC = g++
DEBUG = -g
OPTIMISATIONS = -O3 -ffast-math
CFLAGS = -Wall -c -Wno-write-strings $(OPTIMISATIONS)
LFLAGS = -Wall

all: executable

# Adds debug flags
debug: CFLAGS += $(DEBUG)
debug: LFLAGS += $(DEBUG)
debug: executable

# Adds -march=native to optimisations, which only works for recent gcc versions
native: OPTIMISATIONS += -march=native
native: executable

.PHONY: all debug native executable clean

executable: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(PROGRAM)

clean:
	\rm *.o $(PROGRAM)

######## DEPENDENCIES FOR FILES IN PROJECT ###################################

linalg.h: numerics.h
klustakwik.h: parameters.h log.h util.h numerics.h linalg.h
numerics.h: globalswitches.h
parameters.h: numerics.h
util.h: numerics.h

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
