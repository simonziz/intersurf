/*-----------------------------------------------------------------------------
**  File:       hex_cuda.h
**
**  Author:     Dave Ritchie, 07/08/08
**
**  Purpose:    Header file the C-style Hex-CUDA interface functions.
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2008 D.W. Ritchie, University of Aberdeen.
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

#ifndef hex_cuda_h
#define hex_cuda_h

#if defined(hex_cuda)
#   include "cuda.h"
#   include "cuda_runtime_api.h"
#   include "cufft.h"
#endif

#include "hex_math.h"

#define  HEX_HOST_TO_HOST   0
#define  HEX_HOST_TO_CUDA   1
#define  HEX_CUDA_TO_HOST   2
#define  HEX_CUDA_TO_CUDA   3

#define  HEX_CUDA_MINI      4
#define  HEX_CUDA_HALF      8
#define  HEX_CUDA_BLOCK    16
#define  HEX_CUDA_BLOCK2   32

#ifdef __cplusplus
  extern "C" {
#endif

/*---------------------------------------------------------------------------*/

#if defined(hex_cuda)
typedef struct __align__(8) _cfloat {
#else
typedef struct              _cfloat {
#endif
   float re;
   float im;
} cfloat;

#if defined(hex_cuda)
typedef struct __align__(8) _fpair {
#else
typedef struct              _fpair {
#endif
   float phi;
   float rho;
} fpair;


/*---------------------------------------------------------------------------*/
/* lookup table for converting between "nlm" and "jp" coefficient formats */

typedef struct _Ptable {

    int *p_order;
    int *p_start;
    int *n_table;
    int *l_table;

} Ptable;

Ptable *hex_ptable_make   (int n_max, int p_dim);
void    hex_ptable_free   (Ptable *ptable);

void  hex_cuda_ptable     (int n_max, int p_dim, Ptable *ptable);

// void  hex_cuda_ptable     (int n_max, int p_dim, int *p_order, int *p_start, 
//                         int *l_table);

/*---------------------------------------------------------------------------*/

int   hex_cuda_loadapi    ();
int   hex_cuda_loaddrv    ();

void  hex_cuda_setup      (int first_gpu, int num_gpus);
int   hex_cuda_online     ();
int   hex_cuda_enabled    ();
int   hex_cuda_active_devices();
int   hex_cuda_first_device();
void  hex_cuda_set_device (int device);
void  hex_cuda_sync       ();
void  hex_cuda_info       ();

int   hex_cuda_driver     ();
int   hex_cuda_runtime    ();

void  hex_cuda_device  (int the_device, char *device);
int   hex_cuda_cores   (int the_device);
int   hex_cuda_global  (int the_device);
int   hex_cuda_memory  ();

void *hex_cuda_get     (long size);
void *hex_cuda_get_0   (long size);
void  hex_cuda_free    (void *device_addr, int size);
int   hex_cuda_total   ();
void  hex_cuda_copy    (void *src, void *dst, long size, int direction);
void  hex_cuda_zero    (void *device_addr, long size);

void  hex_cuda_open_1d (int cd, int a_dim);
void  hex_cuda_close_1d(int cd);
void  hex_cuda_open_3d (int cd, int *uvm_dim);
void  hex_cuda_close_3d(int cd);

void hex_cuda_znlm2jp  (int n_vec, int n_max, int j_dim, int p_dim, int k_dim,
                        cfloat *zknlm, cfloat *zkjp);

void hex_cuda_zjp2nlm  (int n_vec, int n_max, int j_dim, int p_dim, int k_dim,
                        cfloat *zkjp, cfloat *zknlm);

void  hex_cuda_fkuvm   (int n_max, int l_dim, int k_dim,
                        int m_min, int m_max, int m_dim,
                        int v_min, int v_max, int v_dim,
                        int u_min, int u_max, int u_dim,
                        cfloat *lambda_mvul, cfloat *s_mvkl, cfloat *f_kuvm);

void  hex_cuda_smvkl   (int n_vec, int n_max, int l_dim, 
                        int k_dim, int j_dim, int p_dim, 
                        cfloat *zakjp, cfloat *zbjp, cfloat *s_mvkl);

void  hex_cuda_skim    (int n_vec, int n_max, int k_dim, 
                        int i_dim, int j_dim, int p_dim, int a_dim,
                        fpair *ajkp, fpair *bjip, cfloat *cs);

void hex_cuda_clmp     (int n_vec, int n_max, int k_dim,
                        cfloat *aknlm, cfloat *bknlm, cfloat *clmp);

void  hex_cuda_ckpm    (int n_vec, int n_max, int k_dim, int p_dim,
                        int j_dim, int q_dim, int q_max, float *omega,
                        cfloat *akjp, cfloat *bkjp, cfloat *cbuf);

void  hex_cuda_qkpm    (int n_vec, int n_max, int k_dim, int p_dim,
                        int j_dim, int q_dim, int q_max, float *omega,
                        cfloat *akjp, cfloat *bkjp, cfloat *qbuf);

void hex_cuda_qkt      (int n_vec, int n_max, int k_dim, int q_dim, int q_max,
                        float *omega, cfloat *aknlm, cfloat *bknlm, cfloat *qbuf);

void  hex_cuda_skim_list(int n_vec, int n_max, int k_dim,
                        int i_dim, int j_dim, int p_dim, int a_dim,
                        int nki_dim, int *nki_list, 
                        fpair *ajkp, fpair *bjip, cfloat *cs);

void  hex_cuda_translate(int n_vec, int *vec_types, int n_max, 
                         int k_pos, int ko_dim, int ki_dim, 
                         int j_dim, int p_dim, float *tm0, float *tm1,
                         fpair *akjp, fpair *tkjp);

void  hex_cuda_setup_fft ();

void  hex_cuda_1d_fft_batch(int cd, int dir, int k_dim, int i_dim, int a_dim, cfloat *cs);

void  hex_cuda_3d_fft_batch(int cd, int dir, int k_dim, cfloat *cs, cfloat *cs_work);

void  hex_cuda_scanner (int k_dim, int i_dim, int a_dim, cfloat *cs,
                        float *scores, int *alphas);

void hex_cuda_1d_fft_setup(int a_dim);

void mat_x_col_fpair   (int n, float mat[], fpair a[], fpair b[]);

void plm_init_cuda     (int l_max, float *pll, float *clmj);

void ylm_all_cuda      (int l_max, float theta, float phi, float *clmj, float *pll, 
                        float *ylm_norm, float *yfn);

void hnl_all_cuda      (int n_max, float r_bar, float *hnl_norm, float *hfn);

void hex_cuda_integrate(int n_max, int n_samp, int n_vec, int a_dim,
                        float *polar, cfloat *weights,
                        float *clmj, float *pll,
                        float *rnl_norm, float *ylm_norm, cfloat *anlm);

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif  /* hex_cuda_h */

/*---------------------------------------------------------------------------*/
