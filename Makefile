CC          := mpicc
CFLAGS      := -Wall -Werror -O3 -c -DUSE_RANDOM_GRAPH=1 -DUSE_RANDOM_SEED=123
LFLAGS      := -Wall -O3
ALL         := floyd-warshall-par.exe

all : $(ALL)


%.exe : %.o graph tests
	$(CC) $(LFLAGS) -o $@ $< graph.o tests.o


%.o : %.c
	$(CC) $(CFLAGS) $<

graph:
	$(CC) $(CFLAGS) graph.c

tests:
	$(CC) $(CFLAGS) tests.c

clean :
	rm -f *.o $(ALL) graph.o tests.o

