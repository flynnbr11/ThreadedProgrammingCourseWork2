CC= gcc #change this to pgcc for cplab
LIB= -lm -fopenmp
OBJ= loops2.o
ORIG= given_loops.c

.c.o: 
	$(CC) -c $<
 
loops : $(OBJ)
	$(CC) $(OBJ) $(LIB) -o $@

.PHONY : original
original : $(ORIG)
	$(CC) $(ORIG) $(LIB) -o $@
 
.PHONY : clean
clean :  
	rm loops
	rm *.o
	rm *~
