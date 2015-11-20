PCC= pgcc #change this to pgcc for cplab
GCC = gcc
LIBG= -lm -fopenmp -O3
LIBP= -lm -mp -O3
OBJ= course.c
ORIG= given_loops.c

.PHONY : new
new : new_guided.c 
	gcc $^ -lm -fopenmp  -o new


.PHONY : struct 
struct : struct_arrays.c
	gcc $^ -lm -fopenmp  -o struct

.PHONY : loops_g
loops_g : $(OBJ)
	$(GCC) $(OBJ) $(LIBG) -o $@

.PHONY : loops_p
loops_p : $(OBJ)
	$(PCC) $(OBJ) $(LIBP) -o $@


.PHONY : original
original : $(ORIG)
	$(GCC) $(ORIG) $(LIBG) -o $@
 
 

	
.PHONY : clean
clean :  
	rm loops
	rm *.o
	rm *~
