#!/usr/bin/env bash 
#$ -cwd -V

rm output.dat 

for proc in 1 2 4 
	do
		export OMP_NUM_THREADS=${proc}
		echo -e "\n For num proc = " ${proc} >> output.dat 
		echo -e "\n" >> output.dat 
	
		./new >> output.dat 
	done
