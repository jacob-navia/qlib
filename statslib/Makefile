#CC=~/lcc/lccarm64/lcc
CC=gcc
CFLAGS=-O2 -g -c 
all:	test libstats.a 
OBJS=bdtrl.o btdtrl.o chdtrl.o fdtrl.o igamil.o incbetl.o incbil.o \
	kolmogor.o nbdtrl.o ndtril.o pdtrl.o psi.o simplestats.o \
	stdtrl.o betainv.o exponentialDistribution.o
libstats.a:	$(OBJS)
	rm -f libstats.a;ar -r libstats.a $(OBJS)
test:	test.o libstats.a stats.h
	gcc -o test test.o libstats.a -lm
test.o:	test.c stats.h
	$(CC) $(CFLAGS) test.c
clean:
	rm -f *.o libstats.a
	
