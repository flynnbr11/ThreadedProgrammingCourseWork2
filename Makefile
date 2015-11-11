PCC= pgcc #change this to pgcc for cplab
GCC = gcc
LIBG= -lm -fopenmp
LIBP= -lm -mp
OBJ= loops2.c
ORIG= given_loops.c


#.PHONY : loops_p
loops_p : $(OBJ)
	$(PCC) $(OBJ) $(LIBP) -o $@

.PHONY : loops_g
loops_g : $(OBJ)
	$(GCC) $(OBJ) $(LIBG) -o $@


.PHONY : original
original : $(ORIG)
	$(CC) $(ORIG) $(LIB) -o $@
 
.PHONY : bullshit 
bullshit : given_loops.c
	pgcc given_loops.c -mp -o given
	
	
.PHONY : clean
clean :  
	rm loops
	rm *.o
	rm *~
