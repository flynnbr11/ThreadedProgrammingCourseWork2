#include <stdio.h>
#include <math.h>
#include<stdlib.h>

#define N 729
#define reps 1 //100
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
	int local_high;	
	int thread_num;
};

void runloop(int loopid)  {

struct all_info *all_threads; 
omp_lock_t *locks;
int size_of_arrays;
int count_finished_threads=0;
int count_all_threads=0;
int tracker=1;
 //struct all_info *all_threads = (struct all_info)malloc(sizeof(all_threads)*nthreads);

#pragma omp parallel default(none) private(size_of_arrays) shared(loopid, all_threads, locks, count_finished_threads, count_all_threads, tracker) 
  {
    int myid  = omp_get_thread_num();
    int nthreads = omp_get_num_threads(); 
    size_of_arrays = sizeof(struct all_info)*nthreads;
    (all_threads) = (struct all_info*) malloc(size_of_arrays);
    all_threads[myid].thread_num=myid;
    locks = malloc(size_of_arrays);   
    int ipt = (int) ceil((double)N/(double)nthreads); 
    int lo = myid*ipt;
    int hi = (myid+1)*ipt;
    if (hi > N) hi = N; 
//	printf("thread %d has lo = %d and hi = %d \n", myid, lo, hi);
    int total_iters = hi-lo;
		int remaining_iters = hi-lo;
		int remaining_iterations = hi;
    int dist = ceil(remaining_iters/nthreads);
    int counter=0;
		int l;
    for(l=0; l<nthreads; l++) {
      omp_init_lock(&(locks[l]));
    }
    
    while(remaining_iterations>0) {
			      //if another thread has replaced remaining_iters since this thread finished its last chunk        
      // does this need to be locked since not updating the value in the shared array? don't think it does   
      omp_set_lock(&(locks[myid]));     
      remaining_iters = all_threads[myid].remaining; 
      lo = all_threads[myid].local_high;
   // omp_unset_lock(&(locks[myid])); //sync in case another thread has stolen from me

      dist = floor( remaining_iters / nthreads ) + 1;        
      hi = lo + dist; 
			remaining_iters = total_iters - hi;

    // omp_set_lock(&(locks[myid]));     
      all_threads[myid].local_high = hi;
      all_threads[myid].remaining = remaining_iters;
      omp_unset_lock(&(locks[myid])); 

// 	printf("thread : %d lo = %d hi = %d \n", myid, lo, hi);
      switch (loopid) { 
          case 1: loop1chunk(lo,hi); break;
          case 2: loop2chunk(lo,hi); break;
      } 
      counter += hi-lo;
      remaining_iterations = total_iters - counter;
//      lo = hi;
    }
	 count_all_threads +=counter;  
   count_finished_threads +=1;
  printf("Thread %d finished chunk with %d iters remaining. Now %d threads out of %d are finished. \n", myid, remaining_iterations, count_finished_threads, nthreads);
	
	while(tracker==1) {
		  if(remaining_iterations <= 0 && count_finished_threads < nthreads) {
				#pragma omp critical 
				{
					printf("thread %d in critical region \n", myid);
					int max=0;
					int thread_to_steal_from = nthreads+1;
					int current_thread, steal_low, steal_high;
					omp_set_lock(&locks[thread_to_steal_from]);	
					int stolen_chunk;
					int old_max_thread = thread_to_steal_from;
			
					for(current_thread =0 ; current_thread < nthreads; current_thread++) { 
						omp_set_lock(&locks[current_thread]);
						if(all_threads[current_thread].remaining > max) {
							old_max_thread = thread_to_steal_from;
							max = all_threads[current_thread].remaining; // reassign maximum remainder from any thread
							thread_to_steal_from = current_thread; //update which thread is most loaded
							printf("thread %d finds that thread %d is loaded  with %d remaining iters \n", myid, current_thread, max);
						}
						omp_unset_lock(&locks[old_max_thread]); //unlock old max thread
						omp_unset_lock(&locks[current_thread]);
						omp_set_lock(&locks[thread_to_steal_from]);
					}
					printf("outside for loop in critical region \n");
					
					if(thread_to_steal_from < nthreads) {
							printf("thread %d finds thread %d to be most loaded \n", myid, thread_to_steal_from);
						omp_set_lock(&locks[thread_to_steal_from]);
							printf("after locking thread %d \n", thread_to_steal_from);
							steal_low = all_threads[thread_to_steal_from].local_high;
							steal_high = steal_low + ((all_threads[thread_to_steal_from].remaining)/nthreads);
							stolen_chunk = steal_high - steal_low;
							all_threads[thread_to_steal_from].remaining -= stolen_chunk; //change how many iterations are left now
							all_threads[thread_to_steal_from].local_high = steal_high; //make sure owner thread starts from new upper limit
						omp_unset_lock(&locks[thread_to_steal_from]);
								
						//We have now found the most loaded thread, and assigned it to thread_to_steal_from
						 //Stop when we reach 0
						if(stolen_chunk >0) {
							switch (loopid) { 
						    case 1: loop1chunk(steal_low,steal_high); break;
						    case 2: loop2chunk(steal_low,steal_high); break;
							} 
						}
					printf("Thread %d stole %d iterations from thread %d \n ", myid, stolen_chunk, thread_to_steal_from);

					}
					else {
						printf( "thread %d did not steal any iterations \n", myid);
						tracker=0;
					}
					printf("thread %d leaving critical region \n", myid);
				}
			}
		}
  }
  free( all_threads );
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
 

