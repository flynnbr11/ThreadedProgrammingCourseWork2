#include <stdio.h>
#include <math.h>

#include<stdlib.h>

#define N 729
#define reps 100
#include <omp.h> 

double a[N][N], b[N][N], c[N];
int jmax[N];  
void init1(void);
void init2(void);
void runloop(int); 
void loop1chunk(int, int);
void loop2chunk(int, int);
void valid1(void);
void valid2(void);
/*
	We declare functions to work out what iterations a thread should work on next
*/
void my_work(int myid, int* lo, int* hi, int* thread_maxima, int* remainders, int* work_remaining, int nthreads);
void steal_work(int nthreads, int* lo, int* hi,  int* most_loaded_thread, int* work_remaining, int* remainders, int* thread_maxima);

int main(int argc, char *argv[]) { 

  double start1,start2,end1,end2;
  int r;

  init1(); 

  start1 = omp_get_wtime(); 

  for (r=0; r<reps; r++){ 
    runloop(1);
  } 

  end1  = omp_get_wtime();  

  valid1(); 

  printf("Total time for %d reps of loop 1 = %f\n",reps, (float)(end1-start1)); 


  init2(); 

  start2 = omp_get_wtime(); 

  for (r=0; r<reps; r++){ 
    runloop(2);
  } 

  end2  = omp_get_wtime(); 

  valid2(); 

  printf("Total time for %d reps of loop 2 = %f\n",reps, (float)(end2-start2)); 

} 

void init1(void){
  int i,j; 

  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      a[i][j] = 0.0; 
      b[i][j] = 3.142*(i+j); 
    }
  }

}

void init2(void){ 
  int i,j, expr; 

  for (i=0; i<N; i++){ 
    expr =  i%( 3*(i/30) + 1); 
    if ( expr == 0) { 
      jmax[i] = N;
    }
    else {
      jmax[i] = 1; 
    }
    c[i] = 0.0;
  }

  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      b[i][j] = (double) (i*j+1) / (double) (N*N); 
    }
  }
 
} 

void runloop(int loopid)  {
	/*
	 dynamically allocate arrays to hold the current number of remaining 
	 iterations belonging to a thread, as well as the highest numbered iteration 
	 that thread should do 
 */
	int* remainders;
	int* thread_maxima; 
	#pragma omp parallel default(none) shared(loopid, remainders, thread_maxima) 
  {
    int myid  = omp_get_thread_num();
    int nthreads = omp_get_num_threads(); 
    int ipt = (int) ceil((double)N/(double)nthreads); 
    int lo = myid*ipt;
    int hi = (myid+1)*ipt;
    if (hi > N) hi = N; 
    int thread_upper_lim = hi;
    int total_iters= hi-lo;
		#pragma omp single //allocate memory to shared arrays
    {
      remainders = (int*) malloc(sizeof(int)*nthreads); 
      thread_maxima = (int*) malloc(sizeof(int)*nthreads);
    }
 		int i, work_remaining, most_loaded_thread;
	  /*
	  	Initialise the shared arrays
	  */
	  
	  thread_maxima[myid] = hi;
		remainders[myid] = hi - lo;    
		
		/*
			We need all arrays to have filled in the shared arrays before 
			others can try to access them, or else we get segmentation faults. 
			We use a barrier here since it is an inexpensive place to do so.
		*/
		#pragma omp barrier
		int finished_work= 0; // When this "flag" is changed, that means all \
															iterations are complete and the program terminates

		while ( finished_work == 0 ) {  
			work_remaining=0;
				#pragma omp critical //Critical region to determine what chunk of iterations to work on next
				{
					most_loaded_thread = myid; // so that rem[thread] > rem[most_loaded] starts by comparing to 0
					work_remaining = 0;
					/* 
						Work out values for lo and hi in a function
						Which function depends on whether or not a thread has finished 
						its local set. 
					*/

					if(remainders[myid] > 0) { 
						my_work(myid, &lo, &hi, thread_maxima, remainders, &work_remaining, nthreads); 
					} 
					
					else if( remainders[myid] <= 0 ) {
						steal_work(nthreads, &lo, &hi, &most_loaded_thread, &work_remaining, remainders, thread_maxima);
					}
				
				} //end critical region
				/* 
					If this thread has found some iterations to work on,
					it raises a flag (work_remaining)
					Then the loop functions are computed outside the critical region
				*/
			if(work_remaining == 1 ) {	
				switch (loopid) { 
						case 1: loop1chunk(lo,hi); break;
						case 2: loop2chunk(lo,hi); break;
				}
			}
			else {
				finished_work = 1; // when nothing has switched the work_remaining flag 
			}
		} // end while loop
 #pragma omp single
	  {
		  free(remainders);
		  free(thread_maxima);
	  }
	} // end parallel region

} // end runloop function


void my_work(int myid, int* lo, int* hi, int* thread_maxima, int* remainders, int* work_remaining, int nthreads) {
	*lo = thread_maxima[myid] - remainders[myid]; //Reassign the value of lo to the previous hi value
	int dist = (int) ceil( remainders[myid] / nthreads ); //Guided like calculation of what to do next
	if(dist == 0 ) dist= 1;
	*hi = *lo + dist;
	if(*hi > thread_maxima[myid]) *hi = thread_maxima[myid];  //ensure we don't surpass 
	remainders[myid] = thread_maxima[myid] - *hi;
			*work_remaining = 1;
}

void steal_work(int nthreads, int* lo, int* hi,  int* most_loaded_thread, int* work_remaining, int* remainders, int* thread_maxima) {
int thread;
	for( thread = 0 ; thread < nthreads; thread++ ) { 
		if(remainders[thread] > remainders[*most_loaded_thread] ) {
			*most_loaded_thread = thread;
			*work_remaining = 1;
		}
	}
	if(*work_remaining == 1) {
		*lo = thread_maxima[*most_loaded_thread] - remainders[*most_loaded_thread];
		int dist = (int) ceil ((double) remainders[*most_loaded_thread] / (double) nthreads);
		if(dist == 0) dist = 1;
		*hi = *lo + dist;
		if( *hi > thread_maxima[*most_loaded_thread] ) *hi = thread_maxima[*most_loaded_thread] ;
		remainders[*most_loaded_thread] = thread_maxima[*most_loaded_thread] - *hi ;
	}		
}

void loop1chunk(int lo, int hi) { 
  int i,j; 
  for (i=lo; i<hi; i++){ 
     for (j=N-1; j>i; j--){
       a[i][j] += cos(b[i][j]);
     } 
  }
} 



void loop2chunk(int lo, int hi) {
  int i,j,k; 
  double rN2; 

  rN2 = 1.0 / (double) (N*N);  

  for (i=lo; i<hi; i++){ 
    for (j=0; j < jmax[i]; j++){
      for (k=0; k<j; k++){ 
      	c[i] += (k+1) * log (b[i][j]) * rN2;
      } 
    }
  }

}

void valid1(void) { 
  int i,j; 
  double suma; 
  
  suma= 0.0; 
  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      suma += a[i][j];
    }
  }
  printf("Loop 1 check: Sum of a is %lf\n", suma);

} 


void valid2(void) { 
  int i; 
  double sumc; 
  
  sumc= 0.0; 
  for (i=0; i<N; i++){ 
    sumc += c[i];
  }
  printf("Loop 2 check: Sum of c is %f\n", sumc);
} 
 

