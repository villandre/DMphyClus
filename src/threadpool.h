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


#ifndef _THREADPOOL_H_
#define _THREADPOOL_H_

//TODO Refactor to better defines
// Errors
#define THREADPOOL_INVALID          -1
#define THREADPOOL_LOCK_FAILURE     -2
#define THREADPOOL_QUEUE_FULL       -3
#define THREADPOOL_SHUTDOWN         -4
#define THREADPOOL_THREAD_FAILURE   -5
#define THREADPOOL_MEMORY_ERROR     -6
#define THREADPOOL_INVALID_SIZE     -7

// Flags
#define THREADPOOL_GRACEFUL         1

typedef struct threadpool_t threadpool_t;
threadpool_t *threadpool_create(int thread_count, int queue_size, int flags);
int threadpool_add(threadpool_t *pool, void (*routine)(void *), void *arg, int flags);
int threadpool_destroy(threadpool_t *pool, int flags);
int threadpool_free(threadpool_t *pool);

int threadpool_get_thread_count(threadpool_t * pool);
int threadpool_get_queue_size(threadpool_t * pool);
int threadpool_is_started(threadpool_t * pool);
int threadpool_is_shutdown(threadpool_t * pool);
int threadpool_new_worker(threadpool_t *pool);
int threadpool_del_worker(threadpool_t *pool); // This is untested


#endif 

