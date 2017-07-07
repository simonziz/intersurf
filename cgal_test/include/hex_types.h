/*-----------------------------------------------------------------------------
**  File:       hex_types.h
**
**  Author:     Dave Ritchie, 10/05/03
**
**  Purpose:    Header file for the main Hex data types and parameters
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2003  D.W. Ritchie, University of Aberdeen.
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

#ifndef hex_types_h
#define hex_types_h

#include <stdint.h>

// ifdef __cplusplus
// extern "C" {
//  endif

/*---------------------------------------------------------------------------*/

#if defined(__GCC__) || defined(hex_icc)
#define HEX_INLINE static __inline__
#else
#define HEX_INLINE static
#endif

#if defined(hex_win32) 
#define HEX_HANDLE HANDLE   
#define HEX_MODULE HMODULE
#define HEX_SYMBOL FARPROC
#else
#define HEX_HANDLE int
#define HEX_MODULE void *
#define HEX_SYMBOL void *
#endif

/*---------------------------------------------------------------------------*/

#define MAX_MOL              4        /* maximum no. of molecules */
#define TWO_MOL              2        /* receptor and ligand only */
#define MAX_PATHNAME       257        /* file name length in bytes */
#define  MAX_OUTPUT_WIDTH  120        /* widest text output line width */
#define  TEXT_BUFFER_DIM  4000        /* message text buffer/growth size */

#define PI            3.14159265358979323846
#define BIG_INT       2147483647  /* largest positive 32-bit integer */

#define FALSE         0
#define TRUE          1

#define MAX_OPLS       100        /* max no. OPLS types (for LJ potentials) */
#define MAX_HBOND       10        /* max no. h-donor/acceptor (12-10) types */

#define TM_MUTEX       0          /* mutexes are used to serialise access to */
#define THREAD_MUTEX   1          /* non-thread-safe critical sections */
#define DIST_MUTEX     2
#define CUDA_MUTEX     3
#define MSG_MUTEX      4
#define YLM_MUTEX      5  
#define MATH_MUTEX     6
#define TEST_MUTEX     7
#define VM_MUTEX       8
#define RM_MUTEX       9
#define ALIGN_MUTEX   10
#define CACHE_MUTEX   11
#define TILE_MUTEX    12
#define TRACE_MUTEX   13
#define SPARE_MUTEX1  14
#define SPARE_MUTEX2  15

#define RESULTS_LOCK   0          /* locks provide a way to synchronise access */
#define TASK_LOCK      1          /* (read/write/wait) to shared data structures */
#define OTHER_LOCK     2
#define FINAL_LOCK     3
#define THREAD_LOCK    4
#define CACHE_LOCK     5
#define SKIN_LOCK      6
#define PROGRESS_LOCK  7
#define REFINE_LOCK    8
#define ROT_LOCK       9
#define ALIGN_LOCK    10
#define TILE_LOCK     11
#define TRACE_LOCK    12
#define MSA_LOCK      13
#define SPARE_LOCK1   14
#define SPARE_LOCK2   15
#define SPARE_LOCK3   16
#define SPARE_LOCK4   17
#define SPARE_LOCK5   18
#define SPARE_LOCK6   19
#define SPARE_LOCK7   20
#define ALIGN_LOCK2   21
#define TRACE_LOCK2   22
#define TILE_LOCK2    23
//  up to max 31 currently


#define CONS_SLOT0   (MAX_THREAD_SLOTS + 0)
#define CONS_SLOT1   (MAX_THREAD_SLOTS + 1)
#define CONS_SLOT2   (MAX_THREAD_SLOTS + 2)
#define CONS_SLOT3   (MAX_THREAD_SLOTS + 3)
#define CONS_SLOT4   (MAX_THREAD_SLOTS + 4)
#define CONS_SLOT5   (MAX_THREAD_SLOTS + 5)
#define CONS_SLOT6   (MAX_THREAD_SLOTS + 6)
#define CONS_SLOT7   (MAX_THREAD_SLOTS + 7)

#define TASK_STRUCT    1          /* readability constants for shared memory */
#define SCAN_STRUCT    2
#define SCAN_BLOCK     3
#define SEARCH_BLOCK   4

#define SHAPE_SCAN     0          /* first task type - not used */

//define SCAN_BLOCK_SIZE  1000    /* size of docking search producer buffer */
#define SCAN_BLOCK_SIZE   512     /* size of docking search producer buffer */

#define MM_BLOCK_SIZE    1024     /* size of docking refinement buffer */

#define SEARCH_BLOCK_SIZE (1638*8*8)  /* gives < 8*32K SearchBlock size */
                                                                                
#define MOL_MC_BIT     0x01       /* main-chain atom */
#define MOL_ACC_BIT    0x02       /* accessible atom */
#define MOL_HD_BIT     0x04       /* h-bond donor atom */
#define MOL_HA_BIT     0x08       /* h-bond acceptor atom */

#define MOL_PDB_TYPE   0
#define MOL_SDF_TYPE   1

#define MAX_FAC     170           /* highest factorial with IEEE doubles */
#define MAX_SPOW    512           /* highest order symbolic power */
#define MAX_SER     400           /* limit for infinite series expansions */
#define MAX_S        32           /* highest sum to n - see sn1...sn4 */
                                                                                            
#define MAX_L        35           /* highest order L for Ylm(theta,phi) */
#define MAX_N   (MAX_L+1)         /* highest order N for Rn(r) */
#define N_DIM   (MAX_N+1)         /* i.e. N dimensions for 0:31 */
#define L_DIM   (MAX_L+1)         /* i.e. L dimensions for 0:30 */
#define K_DIM   (MAX_K+1)         /* i.e. K dimensions for 0:30 */
#define M_DIM (2*MAX_L+1)         /* i.e. M dimensions for -30:30 */
#define MAX_B   (2*N_DIM)         /* largest binomial coefficient */
                                                                                            
#define N_PLM            200      /* accuracy of P[lm] integration */
#define NRK_STEPS        500      /* no. Runge-Kutta (4th-order) steps to use*/
                                  /* when calculating hypergeometric fn.s    */

#define MAX_LAMBDA 16             /* size of sin(theta)^m/2 Taylor expansion */

#define HEX_PROGRESS_START    -1  /* for reporting docking "progress" */
#define HEX_PROGRESS_END       0 
#define HEX_PROGRESS_TOTAL     1
#define HEX_PROGRESS_REFINE    2
#define HEX_PROGRESS_SEARCH    3
#define HEX_PROGRESS_SCAN      4
#define HEX_PROGRESS_SETUP     5

typedef const char          rchar;   /* readonly character (damned C++) */
typedef unsigned char       byte;    /* for convenience */
typedef unsigned int        uint;    /* for convenience */
typedef short int           sint;    /* i.e. 2-byte signed integers */
typedef unsigned short int  usint;   /* i.e. 2-byte unsigned integers */
typedef unsigned long       quadra;  /* i.e. 4-byte integers, same as GLuint */
typedef long double         hexa;    /* i.e. hexadecimal (16) bytes, we hope */

#if defined(hex_win32)
//typedef long int      octa;  /* i.e. 8 byte integers, we hope */
typedef long long      octa;
//typedef int64_t       octa;  
#else
typedef int64_t        octa; 
#endif

/* We *need* to have O_BINARY when opening binary files on Windows. 
   This is a lazy way of stopping errors when compiling on Unix 
*/

#if defined(hex_win32)
#define HEX_O_BINARY O_BINARY
#else
#define HEX_O_BINARY 0
#endif

/*---------------------------------------------------------------------------*/

// ifdef __cplusplus
// }
// endif

#endif  /* hex_types_h */

/*---------------------------------------------------------------------------*/

