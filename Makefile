PCC= pgcc #change this to pgcc for cplab
GCC = gcc
LIBG= -lm -fopenmp -O3
LIBP= -lm -mp -O3
OBJ= loops2.c
ORIG= given_loops.c


.PHONY : loops_p
loops_p : $(OBJ)
	$(PCC) $(OBJ) $(LIBP) -o $@

.PHONY : loops_g
loops_g : $(OBJ)
	$(GCC) $(OBJ) $(LIBG) -o $@


.PHONY : original
original : $(ORIG)
	$(GCC) $(ORIG) $(LIBG) -o $@
 

	
.PHONY : clean
clean :  
	rm loops
	rm *.o
	rm *~
