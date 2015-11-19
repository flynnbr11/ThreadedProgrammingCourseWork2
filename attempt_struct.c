#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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
	int local_low;
	int local_high;	
	int thread_upper_lim;
	int change_flag;
};

void runloop(int loopid)  {

  struct all_info *all_threads; // array of structs of the form defined above
  omp_lock_t *locks; // array of locks
  int size_of_arrays;
  int finished_threads_counter=0;

  #pragma omp parallel default(none) private(size_of_arrays) shared(loopid, all_threads, locks, finished_threads_counter) 
  {
      int myid  = omp_get_thread_num();
      int nthreads = omp_get_num_threads(); 
      int ipt = (int) ceil((double)N/(double)nthreads); 
      int lo = myid*ipt;
      int hi = (myid+1)*ipt;
      if (hi > N) hi = N; 
      int thread_max =hi;
  //dynamically allocate memory for array of structs
      size_of_arrays = sizeof(struct all_info)*nthreads;
      (all_threads) = (struct all_info*) malloc(size_of_arrays);
      locks = malloc(size_of_arrays);   
      int l;
      for(l=0; l<nthreads; l++) {
        omp_init_lock(&(locks[l]));
      }
      all_threads[myid].thread_upper_lim = hi;
      int steal_low=0; int steal_high=0; int steal_dist=0;
      int total_iters = hi-lo;
      int remaining_iters = hi-lo;
      int dist = ceil(remaining_iters/nthreads);
      int counter=0;
      all_threads[myid].remaining = hi;
      hi = lo;

      all_threads[myid].local_high =lo;
      while(all_threads[myid].remaining > 0) {
        //should this be here? will this wait until the next code has exclusive access to locks[myid]?
 /*
 				omp_set_lock(&locks[myid]);
       if(all_threads[myid].change_flag == 1) { 
          printf("Thread %d found flag raised so now re-assigning values \n", myid);
          printf("remaining in array = %d \n", all_threads[myid].remaining);
          lo = all_threads[myid].local_high;
          dist = floor(all_threads[myid].remaining/nthreads) + 1;
          hi = lo + dist;
          all_threads[myid].change_flag = 0;
          all_threads[myid].local_low = lo;
          all_threads[myid].local_high = hi;
          all_threads[myid].remaining = remaining_iters;        

          omp_unset_lock(&locks[myid]);
          printf( "after steal, thread %d Reassigned values lo = %d hi = %d, dist = %d \n", myid, lo, hi, dist);
        }
        else {
					lo = hi;
					dist = floor( all_threads[myid].remaining / nthreads ) + 1;
					hi = lo + dist;  
          all_threads[myid].local_low = lo;
          all_threads[myid].local_high = hi;
          all_threads[myid].remaining = remaining_iters;        
          omp_unset_lock(&locks[myid]);
        }
	*/
				omp_set_lock(&locks[myid]);
				lo = all_threads[myid].local_high;
				dist = floor( all_threads[myid].remaining / nthreads ) + 1;
				hi = lo + dist;  
				all_threads[myid].local_low = lo; 
				all_threads[myid].local_high = hi;
				omp_unset_lock(&locks[myid]);
		
        printf( "Thread %d count %d lo = %d hi = %d dist = %d remaining = %d \
        thread limit = %d \n", myid, counter, lo, hi, dist, remaining_iters, thread_max);
     // Sanity check print statement
				if(hi > N) hi = N;
				if(all_threads[myid].remaining <= 0 ) {
					printf("Breaking \n");
					break; 
				}
				else {
		      switch (loopid) { 
		          case 1: loop1chunk(lo,hi); break;
		          case 2: loop2chunk(lo,hi); break;
		      } 
		      
				omp_set_lock(&locks[myid]);     
    		all_threads[myid].remaining -= dist;
				omp_unset_lock(&locks[myid]);
		
				}
				counter ++;   
      } // end while loop

    omp_set_lock(&locks[myid]);
    all_threads[myid].remaining=0;  //to be sure, shouldn't make a difference
    omp_unset_lock(&locks[myid]);

    finished_threads_counter ++;  
    printf("Thread %d has finished allocated work. %d threads finished out of %d \n",myid, finished_threads_counter, nthreads);
    //Now to set up critical region to steal work
  /*  if(finished_threads_counter < nthreads ) { //sanity check
      #pragma omp critical
      {
        printf("Thread %d has entered the critical region\n", myid);
        int current_thread;
        int thread_to_steal_from=0;
        int max_remaining=0;
        int old_max_thread;
        int found_new_max=0;
        steal_high =0;
        for(current_thread=0; current_thread < nthreads; current_thread++) { 
          printf("Inside for loop with current thread = %d \n", current_thread);
          omp_set_lock(&locks[current_thread]);
          if(all_threads[current_thread].remaining > max_remaining) { 
            printf("Inside if loop \n");
            max_remaining = all_threads[current_thread].remaining;
   //         old_max_thread = thread_to_steal_from;
            thread_to_steal_from = current_thread;
            printf(" if-loop: max = %d on thread %d \n", max_remaining, current_thread);
       //    omp_unset_lock(&locks[current_thread]);
       //    omp_unset_lock(&locks[old_max_thread]);
            printf("After unset\n");
       //      omp_set_lock(&locks[thread_to_steal_from]);
           
            printf("After lock reset\n");
            found_new_max = 1;
            printf("New maximum found on thread %d, remaining = %d \n", thread_to_steal_from, max_remaining);
          } //close if loop for new max found
          else {
     //       omp_unset_lock(&locks[current_thread]);
            printf("Did not enter if-loop for thread %d\n", current_thread);
          }
          
        }
        
        if(found_new_max == 0) {
		      //omp_unset_lock(&locks[0]);
		      printf("Thread %d found nothing to steal from\n", myid);
        }
        else { //i.e if found_new_max = 1 => have raised the flag and need to update variables in thread to steal from
          printf("thread %d should steal from thread %d \n", myid, thread_to_steal_from);        
          all_threads[thread_to_steal_from].change_flag=1; //include this during steal
          steal_low = all_threads[thread_to_steal_from].local_high;
          steal_dist = floor(all_threads[thread_to_steal_from].remaining / nthreads) +1;
          steal_high = steal_low + steal_dist;
          all_threads[thread_to_steal_from].local_low = steal_low;
          all_threads[thread_to_steal_from].local_high = steal_high;
          all_threads[thread_to_steal_from].remaining -= steal_dist;
          if(all_threads[thread_to_steal_from].remaining <=0) all_threads[thread_to_steal_from].remaining =0;
          printf("Thread %d stealing %d from thread %d; from lo = %d to hi =%d \n", myid, steal_dist, thread_to_steal_from, steal_low, steal_high);
          //omp_unset_lock(&locks[thread_to_steal_from]);
          
          //run loop for steal_low to steal_high
        }
        
			for(current_thread=0; current_thread<nthreads; current_thread ++) {
 				omp_unset_lock(&locks[current_thread]);					
			}
      
        
      printf("thread %d leaving critical region \n", myid);
      }// end critical region
    }    //close if(finished counter)
    else printf("Thread %d not entering critical region because %d threads have already finished \n", myid, finished_threads_counter);
    
    if(steal_high > 0) { 
      printf("Thread %d doing iterations %d to %d\n", myid, steal_low, steal_high);
      switch (loopid) { 
        case 1: loop1chunk(steal_low, steal_high); break;
        case 2: loop2chunk(steal_low, steal_high); break;
      } 
    }

*/
  } // close parallel region
  free(all_threads);
  free(locks);



} //end runloop function

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
 

