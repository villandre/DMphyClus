/*
 This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
* Original coder: mbrossard
* Original code: https://github.com/mbrossard/threadpool
* This library was heavily modified by Calvin Ference
*/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include "threadpool.h"

typedef enum {
  immediate_shutdown = 1,
  graceful_shutdown  = 2
} THREADPOOL_SHUTDOWN_t;


typedef struct {
  void (*function)(void *);
  void *argument;
} threadpool_task_t;


struct threadpool_t {
  pthread_mutex_t lock;
  pthread_cond_t notify;
  pthread_t *threads;
  threadpool_task_t *queue ;
  int thread_count;
  int queue_size;
  int head;
  int tail;
  int count;
  int shutdown;
  int started;
  int threadID;
  bool readQueue;
  pthread_spinlock_t ID_spin;
};


/*
* Generic threadpool thread. This is just a simple worker in an infinite loop.
*/
static void *threadpool_thread(void *threadpool) {
  threadpool_t *pool = (threadpool_t *)threadpool;
  threadpool_task_t task;
  int threadID;
  int * arg;
  
  pthread_spin_lock(&pool->ID_spin);
  threadID = pool->threadID;
  pool->threadID++;
  pthread_spin_unlock(&pool->ID_spin);
  
  fprintf(stdout, "[INFO] (%s) Worker Started - ThreadID %i\n", __func__, threadID);
  
  while(1) {
    // Lock the threadpool
    pthread_mutex_lock(&(pool->lock));
    
    // Wait on condition variable if no tasks are available.
    while((pool->count == 0) && (!pool->shutdown))
      pthread_cond_wait(&(pool->notify), &(pool->lock));
    
    // Check if the threadpool is shutting down. if so, time to kill this thread.
    if((pool->shutdown == immediate_shutdown) || ((pool->shutdown == graceful_shutdown) && (pool->count == 0)))
      break;
    
    // Grab next task Modify this to determine which queue takes the task.
    // [CALVIN] TODO: Refactor this, Kinda ugly and I feel his list structure is slightly broken.
    task.function = pool->queue[pool->head].function;
    task.argument = pool->queue[pool->head].argument;
    pool->head += 1;
    pool->head = (pool->head == pool->queue_size) ? 0 : pool->head;
    pool->count -= 1;
    
    // Add threadID to argument list
    arg = (int *) task.argument;
    arg[0] = threadID;
    task.argument = (void *) arg;
    
    // Unlock 
    pthread_mutex_unlock(&(pool->lock));
    
    // Get to work maggot!! *snaps whip*
    (*(task.function))(task.argument);
    
  }
  
  pool->started--;
  fprintf(stdout, "[INFO] (%s) Worker Stopped - ThreadID %i\n", __func__, threadID);
  pthread_mutex_unlock(&(pool->lock));
  pthread_exit(NULL);
  return(NULL);
}

/*
* Following functions are a couple of getters, nothing fancy.
*/ 
int threadpool_get_thread_count(threadpool_t * pool){
  return pool->thread_count;
}

int threadpool_get_queue_size(threadpool_t * pool){
  return pool->queue_size;
}

int threadpool_is_started(threadpool_t * pool){
  return (pool->started > 0) ? 1 : 0;
}

int threadpool_is_shutdown(threadpool_t * pool){
  return (pool->shutdown > 0) ? 1 : 0;
}

/*
* This creates an empty threadpool and starts the working threads.
*/
threadpool_t *threadpool_create(int thread_count, int queue_size, int flags) {
  threadpool_t *pool = NULL;
  int i;
  
  // Check for stupid values... or maybe just convert this to unsigned... 
  if(thread_count <= 0 || thread_count > 4096 || queue_size <= 0 || queue_size > 400000)
    goto err;
  
  pool = (threadpool_t *) malloc(sizeof(threadpool_t));
  if(pool == NULL)
    goto err;
  
  // Initialize the pool
  pool->thread_count = 0;
  pool->queue_size = queue_size;
  pool->head = pool->tail = pool->count = 0;
  pool->shutdown = pool->started = 0;
  pool->threadID = 0;
  
  // Allocate thread and task queue 
  pool->threads = (pthread_t *) malloc(thread_count * sizeof(pthread_t));
  pool->queue = (threadpool_task_t *)malloc(sizeof(threadpool_task_t) * queue_size);
  
  // Initialize spinlock
  if(pthread_spin_init(&(pool->ID_spin), 0) != 0)
    goto err;
  // Initialize mutex and conditional variable
  if(pthread_mutex_init(&(pool->lock), NULL) != 0)
    goto err;
  if(pthread_cond_init(&(pool->notify), NULL) != 0)
    goto err;
  if(pool->threads == NULL)
    goto err;
  if(pool->queue == NULL)
    goto err;
  
  // worker thread creation
  for(i = 0; i < thread_count; i++) {
    if(pthread_create(&(pool->threads[i]), NULL, threadpool_thread, (void*)pool) != 0) {
      threadpool_destroy(pool, 0);
      goto err;
    }
    pool->thread_count++;
    // TODO refactor started
    pool->started++;
  }
  
  // TODO add scheduler thread
  
  return pool;
  
  err:
    if(pool) {
      threadpool_free(pool);
    }
    return NULL;
}

int threadpool_new_worker(threadpool_t *pool){
  int err = 0;
  pthread_t *new_threads;
  
  
  if(pool == NULL)
    return THREADPOOL_INVALID;
  
  if(pthread_mutex_lock(&(pool->lock)) != 0)
    return THREADPOOL_LOCK_FAILURE;
  
  
  new_threads = (pthread_t *) realloc(pool->threads, (pool->thread_count+1) * sizeof(pthread_t));
  if(new_threads == NULL)
    return THREADPOOL_MEMORY_ERROR;
  
  pool->threads = new_threads;
  
  if(pthread_create(&(pool->threads[pool->thread_count+1]), NULL, threadpool_thread, (void*)pool) != 0)
    err = THREADPOOL_THREAD_FAILURE;
  
  if(!err){
    pool->thread_count++;
    pool->started++;
  }
  
  if(pthread_cond_signal(&(pool->notify)) != 0)
    err = THREADPOOL_LOCK_FAILURE;
  if(pthread_mutex_unlock(&pool->lock) != 0)
    err = THREADPOOL_LOCK_FAILURE;
  
  return err;
}

// UNTESTED
int threadpool_del_worker(threadpool_t *pool){
  int err = 0;
  pthread_t *new_threads;
  int newSize;
  
  newSize = pool->thread_count -1;
  if(newSize < 0)
    return THREADPOOL_INVALID_SIZE;
  
  if(pool == NULL)
    return THREADPOOL_INVALID;
  
  if(pthread_mutex_lock(&(pool->lock)) != 0)
    return THREADPOOL_LOCK_FAILURE;
  
  
  // Shrink the buffer. Since realloc has different behaviours for shrinking
  // depending on the implementation, a simple malloc+memcpy will suffice.
  new_threads = (pthread_t *) malloc(sizeof(pthread_t) * newSize);
  if(new_threads == NULL)
    return THREADPOOL_MEMORY_ERROR;
  
  memcpy(new_threads, pool->threads, newSize);
  pool->threads = new_threads;
  pool->started--;
  pool->thread_count--;
  
  if(pthread_cond_signal(&(pool->notify)) != 0)
    err = THREADPOOL_LOCK_FAILURE;
  if(pthread_mutex_unlock(&pool->lock) != 0)
    err = THREADPOOL_LOCK_FAILURE;         
  return err;
}


int threadpool_add(threadpool_t *pool, void (*function)(void *), void *argument, int flags) {
  int err = 0;
  int next;
  
  if(pool == NULL || function == NULL)
    return THREADPOOL_INVALID;
  
  if(pthread_mutex_lock(&(pool->lock)) != 0)
    return THREADPOOL_LOCK_FAILURE;
  
  next = pool->tail + 1;
  next = (next == pool->queue_size) ? 0 : next;
  
  // Check status of pool
  if(pool->count == pool->queue_size)
    err = THREADPOOL_QUEUE_FULL;
  if(pool->shutdown) 
    err = THREADPOOL_SHUTDOWN;
  
  // Add task to queue
  pool->queue[pool->tail].function = function;
  pool->queue[pool->tail].argument = argument;
  pool->tail = next;
  pool->count += 1;
  
  // wake up any sleeping worker *cracks whip*
  if(pthread_cond_signal(&(pool->notify)) != 0)
    err = THREADPOOL_LOCK_FAILURE;
  if(pthread_mutex_unlock(&pool->lock) != 0)
    err = THREADPOOL_LOCK_FAILURE;
  
  return err;
}

int threadpool_free(threadpool_t *pool) {
  if(pool == NULL || pool->started > 0)
    return -1;
  
  // You are free my minions! Freeeee~~~
  if(pool->threads) {
    free(pool->threads);
    free(pool->queue);
    pool->threads = NULL;
    pool->queue = NULL;
    
    pthread_mutex_lock(&(pool->lock));
    pthread_mutex_destroy(&(pool->lock));
    pthread_cond_destroy(&(pool->notify));
  }
  
  free(pool); 
  pool = NULL;
  return 0;
}

int threadpool_destroy(threadpool_t *pool, int flags) {
  int i, err = 0;
  
  if(pool == NULL)
    return THREADPOOL_INVALID;
  
  // TODO Fix this
  // This was a bug about a mutex unlock not being done... ???
  pthread_mutex_unlock(&(pool->lock));
  if(pthread_mutex_lock(&(pool->lock)) != 0)
    return THREADPOOL_LOCK_FAILURE;
  
  // Already shutting down mate, don't make me slap you...
  if(pool->shutdown)
    err = THREADPOOL_SHUTDOWN;
  
  pool->shutdown = (flags & THREADPOOL_GRACEFUL) ? graceful_shutdown : immediate_shutdown;
  
  if(pool->shutdown == 1)
    fprintf(stdout, "[WARN] (%s) Immediate threadpool shutdown\n", __func__);
  else
    fprintf(stdout, "[INFO] (%s) Gracefully shutting down threadpool\n", __func__);
  
  // Wake up all worker threads to finish any remaining work.
  if((pthread_cond_broadcast(&(pool->notify)) != 0) || (pthread_mutex_unlock(&(pool->lock)) != 0))
    err = THREADPOOL_LOCK_FAILURE;
  
  // Join all worker thread
  for(i = 0; i < pool->thread_count; i++) {
    if(pthread_join(pool->threads[i], NULL) != 0)
      err = THREADPOOL_THREAD_FAILURE;
  }
  
  
  // TODO Fix this
  // Only if everything went well do we deallocate the pool 
  if(!err)
    threadpool_free(pool);
  
  fprintf(stdout, "[INFO] (%s) Threadpool completely shutdown\n", __func__);
  return err;
}