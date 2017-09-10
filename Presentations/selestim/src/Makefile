CFLAGS = -std=c99 -Wall -Wextra -fopenmp -I../dcmt0.6.1b/include
LIBS = -L./dcmt0.6.1b/lib/ -ldcmt -lm
CC = gcc

init.o: init.c init.h
	$(CC) -o init.o -c init.c $(LIBS) $(CFLAGS) $(OPTS)

memory.o: memory.c memory.h
	$(CC) -o memory.o -c memory.c $(LIBS) $(CFLAGS) $(OPTS)

moves.o: moves.c moves.h
	$(CC) -o moves.o -c moves.c $(LIBS) $(CFLAGS) $(OPTS)

mymath.o: mymath.c mymath.h
	$(CC) -o mymath.o -c mymath.c $(LIBS) $(CFLAGS) $(OPTS)

rand.o: rand.c rand.h
	$(CC) -o rand.o -c rand.c $(LIBS) $(CFLAGS) $(OPTS)

read.o: read.c read.h
	$(CC) -o read.o -c read.c $(LIBS) $(CFLAGS) $(OPTS)

write.o: write.c write.h
	$(CC) -o write.o -c write.c $(LIBS) $(CFLAGS) $(OPTS)

all: export OPTS = -O3

debug: export OPTS = -g

fastdebug: export OPTS = -g -O3

profiling: export OPTS = -g -pg

xcode: export OPTS = -ggdb

all debug fastdebug profiling xcode: init.o memory.o moves.o mymath.o rand.o read.o write.o
	cd dcmt0.6.1b/lib/; $(MAKE) $(MFLAGS) 
	$(CC) -o selestim selestim.c init.o memory.o moves.o mymath.o rand.o read.o write.o $(LIBS) $(CFLAGS) $(OPTS)

clean:
	rm -rf *.o selestim
	cd dcmt0.6.1b/lib/; $(MAKE) $(MFLAGS) clean
