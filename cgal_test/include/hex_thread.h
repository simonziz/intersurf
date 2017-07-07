/*-----------------------------------------------------------------------------
**  File:       hex_thread.h
**
**  Author:     Dave Ritchie, 12/09/09
**              Vishwesh Venkatraman
**
**  Purpose:    #include file for the "multi-threading" functions.
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2009 D.W. Ritchie, University of Aberdeen.
**  Copyright (C) 2010 D.W. Ritchie, INRIA.
**
**  This software (or modified copies thereof) is protected by copyright and
**  may not be redistributed in any way without the express permission of the
**  author. Any questions about this copyright notice should be addressed to
**  dave.ritchie@inria.fr.
**
**  If this software was not obtained directly from Dave Ritchie, then it is
**  an unauthorised copy and it should be erased from your computer system,
**  and any associated media should be returned to the author.
**
**---------------------------------------------------------------------------*/

#ifndef hex_thread_h
#define hex_thread_h

#include "hex_sys.h"

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------------*/

#if defined(hex_win32)
typedef HANDLE           hex_thread;
typedef HANDLE           hex_mutex;
typedef HANDLE           hex_cond;
#else
typedef pthread_t        hex_thread;
typedef pthread_mutex_t  hex_mutex;
typedef pthread_cond_t   hex_cond;
#endif

/*---------------------------------------------------------------------------*/

typedef struct _HexTaskSlot {

   int slot_num;                  /* counts from zero for all devices */
   int dev_num;                   /* counts from first of each device type */
   int dev_type;                  /* either HEX_CPU_DEV or HEX_GPU_DEV */
   int dev_memory;                /* GPU device memory */

} HexTaskSlot;

/*---------------------------------------------------------------------------*/

typedef struct {

   int          stack_size;
   int          stack_lock;
   int          stack_top;
   HexTaskSlot *task_stack;

} HexTaskStruct;

/*---------------------------------------------------------------------------*/

/*  Function prototypes... */
/*  ---------------------- */

void        hex_thread_setup     (int max_threads);
void        hex_cond_signal      (int the_cond);
void        hex_cond_wait        (int the_cond, int the_mutex);

void        hex_mutex_acquire    (int locknum);
void        hex_mutex_release    (int locknum);

void        hex_mutex_acquire_w  (int locknum);
void        hex_mutex_release_w  (int locknum);
void        hex_mutex_acquire_r  (int locknum);
void        hex_mutex_release_r  (int locknum);
void        hex_mutex_release_w_ready (int locknum);
void        hex_mutex_release_r_ready (int locknum);

int         hex_mutex_status     (int locknum);

int         hex_thread_slots     ();
int         hex_thread_cpus      ();

// not thread-safe

int         hex_thread_setup_task(int ncpu, int ngpu);
void        hex_thread_push_task (HexTaskSlot ts);
HexTaskSlot hex_thread_pop_task  ();

// thread-safe

HexTaskStruct *hex_makeTaskStruct(int ncpu, int ngpu, int the_lock);
void           hex_freeTaskStruct(HexTaskStruct *t_struct);
void           hex_threadPushTask(HexTaskSlot ts, HexTaskStruct *t_struct);
HexTaskSlot    hex_threadPopTask (HexTaskStruct *t_struct);


void        hex_thread_load_range(int l_min, int l_max);
void        hex_thread_load_list (int nl, int *slots);
int         hex_thread_pop_list  ();

void        hex_thread_start     (int slot, int joinable, 
                                  void *fn(void *), void *data);
void        hex_thread_join      (int slot);
void        hex_thread_exit      ();

#ifdef __cplusplus
}
#endif

#endif  /* hex_thread_h */

/*---------------------------------------------------------------------------*/
