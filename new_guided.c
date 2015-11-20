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

void mywork(int);

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


struct all_info {
	int remaining;
	int current_high;	
	int thread_max;
};

void runloop(int loopid)  {
	int counter = 0; // keep track of how many iters done between all threads
	struct all_info *all_threads;
	int finished_thread_count =0;
	int total_remaining_iters = N;
	int* remainders;
	int* thread_maxima;
	//int work_remaining = 1;
#pragma omp parallel default(none) shared(loopid, all_threads, counter, finished_thread_count, total_remaining_iters, remainders, thread_maxima ) 
  {
    int myid  = omp_get_thread_num();
    int nthreads = omp_get_num_threads(); 
    int size_of_arrays = sizeof(struct all_info)*nthreads;
    (all_threads) = (struct all_info*) malloc(size_of_arrays);

    int ipt = (int) ceil((double)N/(double)nthreads); 
    int lo = myid*ipt;
    int hi = (myid+1)*ipt;
    if (hi > N) hi = N; 
    int thread_upper_lim = hi;
    int total_iters= hi-lo;
		#pragma omp single
    {
      remainders = (int*) malloc(sizeof(int)*nthreads); 
      thread_maxima = (int*) malloc(sizeof(int)*nthreads);
    }
 	  int remaining_iters, dist ;
    //int dist = ceil(remaining_iters/nthreads);
    //if(dist==0) dist = 1;
		//hi = lo;
		int work_remaining, thread, biggest_remaining, end_point;
		int value = 1;
		int most_loaded_thread=0;
		int i;
	  thread_maxima[myid] = hi;
		remainders[myid] = hi - lo;    
		int finished_threads = 0;
		
		#pragma omp barrier
		int finished_work= 0;
		while ( finished_work == 0 ) {  
			work_remaining=0;
			
			if(finished_threads == 0) {
				if (remainders[myid] > 0 ) {
						lo = thread_maxima[myid] - remainders[myid];
						dist = (int) ceil( remainders[myid] / nthreads );
						if(dist == 0 ) dist= 1;
						hi = lo + dist;
						if(hi > thread_upper_lim) hi = thread_upper_lim; 
						counter += dist;
						remainders[myid] -= dist;
						if(remainders[myid] <= 0 ) finished_threads ++;
						work_remaining=1;
					}
				} // will all work in parallel before anyone finishes their initial chunk
				
			else {
				#pragma omp critical 
				{
					most_loaded_thread = myid; // so that rem[thread] > rem[most_loaded] starts by comparing to 0
					work_remaining = 0;
					if(remainders[myid] > 0) { 
						lo = thread_maxima[myid] - remainders[myid];
						dist = (int) ceil( remainders[myid] / nthreads );
						if(dist == 0 ) dist= 1;
					
						hi = lo + dist;
						if(hi > thread_upper_lim) hi = thread_upper_lim; 
						counter += dist;
						remainders[myid] -= dist;
						work_remaining = 1;
					} 
					
					else if( remainders[myid] <= 0 ) {
						for( thread = 0 ; thread < nthreads; thread++ ) { 
							if(remainders[thread] > remainders[most_loaded_thread] ) {
								most_loaded_thread = thread;
								work_remaining = 1;
							}
						}
						if(work_remaining == 1) {
							lo = thread_maxima[most_loaded_thread] - remainders[most_loaded_thread];
							dist = (int) ceil ((double) remainders[most_loaded_thread] / (double) nthreads);
							if(dist == 0) dist = 1;
							hi = lo + dist;
							if( hi > thread_maxima[most_loaded_thread] ) hi = thread_maxima[most_loaded_thread] ;
							remainders[most_loaded_thread] -= dist;
						}						
					
					}
				
					else printf("Thread %d lost \n", myid);
				} //end critical region
			}	
			
			if(work_remaining == 1 ) {	
	//			printf("Thread %d going to do %d to %d\n", myid, lo, hi);		
				switch (loopid) { 
						case 1: loop1chunk(lo,hi); break;
						case 2: loop2chunk(lo,hi); break;
				}
			}
			else {
				finished_work = 1;
			}
		} // end while 1
 #pragma omp single
	  {
		  free(remainders);
		  free(thread_maxima);
	  }
	} // end parallel region

} // end runloop function


void mywork(int myid) {


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
 

