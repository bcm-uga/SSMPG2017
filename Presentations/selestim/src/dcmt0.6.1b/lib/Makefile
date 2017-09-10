CC = gcc
CFLAGS = -Wall -Wmissing-prototypes -O3 -std=c99 -I../include
OBJS = check32.o prescr.o mt19937.o eqdeg.o seive.o genmtrand.o init.o

lib : $(OBJS)
	ar -crus libdcmt.a $(OBJS)

init.o : dci.h init.c
	$(CC) $(CFLAGS) -c init.c

check32.o : dci.h check32.c
	$(CC) $(CFLAGS) -c check32.c

prescr.o : dci.h prescr.c
	$(CC) $(CFLAGS) -c prescr.c

mt19937.o : mt19937.c
	$(CC) $(CFLAGS) -c mt19937.c

eqdeg.o : dci.h eqdeg.c
	$(CC) $(CFLAGS) -c eqdeg.c

seive.o : dci.h seive.c
	$(CC) $(CFLAGS) -c seive.c

genmtrand.o : dci.h genmtrand.c
	$(CC) $(CFLAGS) -c genmtrand.c

clean :
	/bin/rm -f *.o libdcmt.a

oclean :
	/bin/rm -f *.o

