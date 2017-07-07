/*-----------------------------------------------------------------------------
**  File:       hex_sys.h
**
**  Author:     Dave Ritchie, 10/05/03
**
**  Purpose:    Header file for the "sys-level" functions.
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2003 D.W. Ritchie, University of Aberdeen.
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

#ifndef hex_sys_h
#define hex_sys_h

// using namespace std;

#include "hex_types.h"

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <cerrno>
#include <cfloat>

//include <stdlib.h>
//include <stdio.h>
//include <stdarg.h>
//include <malloc.h>

#if defined (hex_darwin) || defined(hex_darwin64)
#include <unistd.h>
#include <sys/sysctl.h>
#include <sys/types.h>
#endif

#if defined(hex_win32)

#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
#include <ctype.h>
#include <direct.h>
#include <process.h>
#include <winsock2.h>
#include <windows.h>

#define strncasecmp(x,y,n) strnicmp(x,y,n)
#define strcasecmp(x,y) stricmp(x,y)

// these were removed for C++ version - src6a

// # define  __S_ISTYPE(mode, mask)  (((mode) & S_IFMT) == (mask))
// # define S_ISDIR(mode)   __S_ISTYPE((mode), S_IFDIR)
#ifndef S_ISDIR
# define S_ISDIR(mode)   (((mode) & _S_IFMT) == _S_IFDIR)
#endif

#else  /* unix-style */

#include <unistd.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <sys/times.h>
#include <sys/wait.h> 
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>
/* include <wait.h> */  /* Darwin */
#include <limits.h>
#include <pwd.h>
#include <netdb.h>
#include <libgen.h>       /* don't need ? */
#include <pthread.h>
#include <dirent.h>
#include <fnmatch.h>
#include <dlfcn.h>
#endif /* hex_win32 */

#ifdef hex_irix
#include <invent.h>
#endif

#include <signal.h>
#include <time.h>
#include <ctype.h>
#include <fcntl.h>
//include <errno.h>
//include <math.h>

#include "zlib.h"
#include "bzlib.h"

#define HEX_FILE_READ         0x00
#define HEX_FILE_WRITE        0x01
#define HEX_FILE_BINARY       0x02
#define HEX_FILE_RWMASK       0x03
#define HEX_FILE_PIPE         0x04
                                                                                
#define HEX_FILE_NONE         0x00
#define HEX_FILE_GZ           0x08
#define HEX_FILE_BZ           0x10
#define HEX_FILE_GZMASK       0x18

#define HEX_FILE_EOF           EOF
#define HEX_FILE_EOD            -9

#define MAX_GPU_SLOTS            8
#define MAX_CPU_SLOTS          128 
#define MAX_THREAD_SLOTS       (MAX_CPU_SLOTS+MAX_GPU_SLOTS)

#define HEX_CPU_DEV              0
#define HEX_GPU_DEV              1

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------------*/

/* patch up missing bits in Linux */
                                                                                
#ifdef hex_linux
typedef unsigned short int ushort_t;
#endif

typedef FILE *fFile;
typedef void *bFile;

typedef struct _hex_file {

   byte mode;                /* bitwise-OR of HEX_FILE_? parameters */
   fFile  f;
   gzFile g;
   bFile  b;
   HEX_HANDLE p;
   HEX_HANDLE pid;
   char    *fname;
   int      lineno;

} HexFile;

/*---------------------------------------------------------------------------*/

typedef void (*hex_exit_cb)     (void);
typedef void (*hex_pipe_cb)     (void);
typedef void (*hex_abort_cb)    (char *text);
typedef void (*hex_logger_cb)   (char *text);
typedef void (*hex_wakeup_cb)   (void);
typedef int  (*hex_progress_cb) (int phase, double progress);
typedef void (*hex_update_cb)   (int thing, va_list ap);

/*---------------------------------------------------------------------------*/

int  hex_setup_tasks    (int ntask_requested, int ncpus_requested);
int  hex_get_ntask      ();
int  hex_start_tasks    (int ntask);
void hex_end_task       (int child_no);
void hex_kill_tasks     ();
void hex_abort_tasks    (char *reason);
void hex_task_handlers  (hex_exit_cb exit_cb, hex_abort_cb abort_cb,
                         hex_pipe_cb pipe_cb);

void hex_check_root     (rchar *prog_name);

void hex_setup_signals  (int batch_mode);
void hex_signal         (int signum, void (*handler) (int));

void hex_ctrlc       (int); /* common signals we may want to handle/ignore */
void hex_sighup      (int);
void hex_sigquit     (int);
void hex_sigterm     (int);
void hex_sigint      (int);
void hex_sigstop     (int);
void hex_sigchld     (int);
void hex_sigalrm     (int);
void hex_sigpipe     (int);
void hex_sigcont     (int);

void hex_sigbus      (int);  /* signals usually caused by programming errors */
void hex_sigseg      (int);
void hex_sigsys      (int);
void hex_sigill      (int);
void hex_sigfpe      (int);

void   hex_zero           (void *addr, long size);            /* equivalent to "bzero" */
void   hex_copy           (void *src, void *dest, long size); /* equivalent to "bcopy" */
char  *hex_index          (rchar *s, int c);                 /* equivalent to "index" */

int    hex_check_dmy      (int d, int m, int y);

HEX_MODULE hex_dlopen     (const char *filename);
HEX_SYMBOL hex_dlsym      (HEX_MODULE module, const char *name);
void       hex_dlclose    (HEX_MODULE module);

void   hex_get_date       (char *date);
void   hex_get_user       (char *name);
void   hex_get_host       (char *host);
void   hex_get_cwd        (char *cwd);

int    hex_umask          (int mask);
int    hex_hostid         ();  /* supply unique(?) host ID */
char  *hex_hostaddr       (char *name); /* supply host IP address string */
void   hex_putenv         (rchar *name, rchar *value); /* set env variable */
char  *hex_getenv         (rchar *name); /* environment variable value */
char  *hex_null_env       (rchar *name); /* supply value or "NULL" */
int    hex_bool_env       (rchar *name, int def); 
int    hex_int_env        (rchar *name, int def); 
rchar *hex_char_env       (rchar *name, rchar *def); 
void   hex_nice           (int nice_val);
int    hex_endian         ();  /* 1 for big-endian, 0 for little-endian */
int    hex_cpus_online    ();  /* processor status */
int    hex_cpus_configured();
int    hex_mb_memory      ();

int    find_val           (double dval, double dlist[], int  length);
int    find_min           (double dlist[], int  length);
int    find_max           (double dlist[], int  length);
int    find_int           (int ival, int ilist[], int  length);
int    find_int           (int  lval, int  llist[], int  length);
int    find_string        (char *string, char *strings, int  stride, int  n);

char  *hex_ran_string     ();                     /* random string */
double hex_ran_range      (double, double);       /* pseudo random no. in given range */

void   hex_init_time      ();  /* recommended */
double hex_get_time       ();  /* thread-safe */
void hex_elapsed_string   (double seconds, char *string);


int    note_time     ();                     /* DEPRECATED allows ovelapping timings */
double delta_time    (int index);

#define hex_vm_wipe(x) { if (x != NULL) { hex_vm_free(x); x = NULL; } }

void   hex_vm_debug   (int flag);
void  *hex_vm_get     (long size);           /* virtual memory (malloc) */
void  *hex_vm_get_0   (long size);
void   hex_vm_free    (void *ptr);
void  *hex_vm_mod     (void *ptr, long size);
char  *hex_vm_strdup  (rchar *ptr); 
int    hex_vm_check   (void *ptr, char *label); 
long   hex_vm_total   (); 
long   hex_vm_max     (); 
long   hex_vm_get_size(void *ptr); 
long   hex_vm_get_tag (void *user_addr);
void   hex_vm_set_tag (void *user_addr, long value);


void  *hex_sm_get     (long size);  /* shared memory (shm) access */
void   hex_sm_free    (void *ptr);
void   hex_sm_exit    ();           /* exit handler */
int    hex_sm_id      (void *ptr);
long   hex_sm_total   ();
long   hex_sm_max     (); 
void  *hex_sm_attach  (int sid);
void   hex_sm_detach  (void *ptr);
void   hex_log        (rchar *fmt, ...);
void   hex_msg        (rchar *fmt, ...);
void   hex_wrn        (rchar *fmt, ...);
void   hex_adv        (rchar *fmt, ...);
void   hex_err        (rchar *fmt, ...);
void   hex_oops       (rchar *fmt, ...);

void  hex_logger_exit     ();
int   hex_logger_status   ();
void  hex_logger_state    (int on);
void  hex_logger_buffer   (int on);
void  hex_logger_stderr   (int on);
char *hex_logger_file     ();
void  hex_logger_name     (rchar *filename);
void  hex_logger_err      (rchar *prefix, rchar *buffer);
void  hex_logger_msg      (rchar *buffer);
void  hex_logger_log      (rchar *buffer);
void  hex_msg_brk         ();
void  hex_tick            ();
hex_logger_cb hex_logger_callback(hex_logger_cb);

HexFile *hex_fopen   (rchar *path, rchar *mode);
HexFile *hex_fdopen  (int fd, rchar *mode, int compression);
HexFile *hex_popen_write (rchar *reader_program,  rchar *arg1, rchar *arg2);
HexFile *hex_popen_read  (rchar *writer_program,  rchar *arg1, rchar *arg2);
int      hex_spawn_file(rchar  *executable_file,  rchar *arg1, rchar *arg2);
int      hex_spawn_list(rchar *executables[], rchar *arg1, rchar *arg2);

int      hex_bbread  (rchar *filename, void *buffer, int nb);
int      hex_bbwrite (rchar *filename, void *buffer, int nb);

int      hex_fclose  (HexFile *file);
int      hex_fflush  (HexFile *file);
void     hex_rewind  (HexFile *file);
int      hex_fgetc   (HexFile *file);
int      hex_fputc   (HexFile *file, int c);
int      hex_fputs   (HexFile *file, char *s);
int      hex_fprintf (HexFile *file, rchar *format, ...);
size_t   hex_fread   (HexFile *file, void *buf, size_t count);
size_t   hex_fwrite  (HexFile *file, void *buf, size_t count);
size_t   hex_fwriteln(HexFile *file, rchar *line);

char    *hex_findfile (rchar *filename, char *result, int npaths, ...);
int      hex_ffindln  (HexFile *file, rchar *keyword, int length, char *line);
int      hex_freadln  (HexFile *file, int length, char *line);
int      hex_freadtok (HexFile *file, char *delim, int length, char *token);
int      hex_ftokenise(HexFile *file, char *delim, int length, char *buf,
                       char ***toks);
int      hex_fline    (HexFile *file);
char    *hex_fname    (HexFile *file);

long     hex_filesize       (rchar *filename);
char    *hex_fileslash      ();
char    *hex_pathslash      (rchar *pathname, char *result);
char    *hex_dirname        (rchar *filename, char *dname);
char    *hex_rootname       (rchar *filename, char *rname);
char    *hex_basename       (rchar *filename, char *bname);
char    *hex_extname        (rchar *filename, char *ename);
char    *hex_endname        (rchar *filename, char *ename);
char    *hex_filename       (rchar *filename, rchar *ext, char *fname);
char    *hex_cwdname        (char *cwd);
char    *hex_catpath        (rchar *dirname, rchar *filename, char *result);
char    *hex_make_path      (rchar *env_variable, rchar *hex_path, char *result);
int      hex_mkdir          (rchar *path, int mode);
int      hex_isdir          (rchar *path);
int      hex_rmdir          (rchar *path, int really_do_it);
char   **hex_scandir        (rchar *path, rchar *pattern, int *n_match);

void     hex_socket_startup ();                          /* once on startup */
int      hex_socket_create  (char *hostname);            /* on the "server" */
int      hex_socket_port    (int socket);                /* on the "server" */
int      hex_socket_accept  (int socket);                /* on the "server" */
int      hex_socket_connect (char *hostname, int port);  /* on the "client" */
int      hex_socket_send    (int socket, const char *buf, int num_bytes);
int      hex_socket_recv    (int socket,       char *buf, int max_bytes);

int    hex_get_cpus         ();  /* supply no. of CPUs to use */
int    hex_get_pid          ();  /* current process ID */

void   hex_alloc_locks      ();  /* create read/write locks using semaphores */
void   hex_free_locks       ();

void   hex_write_acquire    (int locknum); /* unix semaphores - deprecated */
void   hex_write_release    (int locknum);
void   hex_read_acquire     (int locknum);
void   hex_read_release     (int locknum);
void   hex_read_drop        (int locknum);

void   hex_flag_set         (int flagnum);
void   hex_flag_clear       (int flagnum);

void   str_field            (char *line, int start, int end, int base, char *field);
void   str_ljust            (char *dest, char *src, int dsize);

void   str_copy             (char *dest, char *src, int dsize);
void   str_rclip            (char *str);
char  *str_lclip            (char *str);
void   str_alpha            (char *in, char *out);
void   str_upper            (char *s);
void   str_lower            (char *s);
void   str_sanitise         (char *s);
void   str_unsanitise       (char *s);

void   hex_sortc  (char  **c, int  n, int  *order);
void   hex_sorto  (octa   *o, int  n, int  *order);
void   hex_sortl  (int    *l, int  n, int  *order);
void   hex_sortu  (int    *l, int  n, int  *order);
void   hex_sortb  (byte   *b, int  n, int  m, int  *order);
void   hex_sortj  (int    *l, int  n, int  m, int  *order);
void   hex_sortf  (float  *f, int  n, int  *order);
void   hex_sortd  (double *d, int  n, int  *order);
   
int    hex_spanf  (float *f, int  n, int  ndim);
int    hex_spanf_s(float *f, int  n, int  ndim, int  *sort);
int    hex_spanl  (int   *l, int  n, int  ndim);
int    hex_spanl_s(int   *l, int  n, int  ndim, int  *sort);
int    hex_spano  (octa  *o, int  n, int  ndim);
int    hex_spano_s(octa  *o, int  n, int  ndim, int  *sort);
int    hex_spani  (int   *i, int  n, int  ndim);
int    hex_spanc  (char  *c, int  n, int  ndim);

void hex_find_ispans  (int ni, int *is, int *nspan, int *spans);
void hex_find_ispans_s(int ni, int *is, int *sort, int *nspan, int *spans);

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif  /* hex_sys_h */

/*---------------------------------------------------------------------------*/
