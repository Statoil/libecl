#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <thread_pool.h>
#include <util.h>

struct thread_pool_struct {
  int        pool_size;
  int        jobs_running;
  pthread_t *thread_list;
};





static void thread_pool_resize(thread_pool_type * pool, int new_size) {
  pool->pool_size   = new_size;
  pool->thread_list = realloc(pool->thread_list , new_size * sizeof * pool->thread_list);
}



void thread_pool_join(thread_pool_type * pool) {
  int i;
  if (pool->pool_size == 0)
    return;
  else {
    for (i=0; i < pool->jobs_running; i++)
      pthread_join(pool->thread_list[i] , NULL);  /* Second argument: void **value_ptr */
    pool->jobs_running = 0;
  }
}


thread_pool_type * thread_pool_alloc(int pool_size) {
  thread_pool_type * pool = malloc(sizeof *pool);
  pool->thread_list = NULL;
  thread_pool_resize(pool , pool_size);
  pool->jobs_running = 0;
  return pool;
}


void thread_pool_add_job(thread_pool_type * pool , 
			 void *(start_func) (void *) , void *arg) {

  if (pool->pool_size == 0) 
    start_func(arg);
  else {

    if (pool->jobs_running == pool->pool_size) 
      thread_pool_join(pool);
    
    {
      int pthread_return = pthread_create(&pool->thread_list[pool->jobs_running] , NULL , start_func , arg);
      if (pthread_return != 0) 
	util_abort("%s: failed to add new job pthread_create return value: %d.\n",__func__ , pthread_return);
    }

    pool->jobs_running++;
  }
  
}
  
void thread_pool_free(thread_pool_type * pool) {
  if (pool->thread_list != NULL) free(pool->thread_list);
  free(pool);
}


