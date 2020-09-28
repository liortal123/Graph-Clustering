FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: cluster

cluster:  cluster.o matA_sparse.o Utility.o matB.o
	gcc cluster.o matA_sparse.o Utility.o matB.o -o cluster $(LIBS) 

clean:
	rm -rf *.o cluster

cluster.o: cluster.c matA_sparse.h Utility.h matB.h
	gcc $(FLAGS) -c cluster.c 
matA_sparse.o: matA_sparse.c matA_sparse.h Utility.h
	gcc $(FLAGS) -c matA_sparse.c 
Utility.o: Utility.c Utility.h matA_sparse.h matB.h
	gcc $(FLAGS) -c Utility.c 
matB.o: matB.c matA_sparse.h Utility.h
	gcc $(FLAGS) -c matB.c