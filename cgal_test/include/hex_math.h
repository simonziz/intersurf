/*-----------------------------------------------------------------------------
**  File:       hex_math.h
**
**  Author:     Dave Ritchie, 10/05/03
**
**  Purpose:    Header file for the "math-level" functions.
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) : D.W. Ritchie, University of Aberdeen.
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

#ifndef hex_math_h
#define hex_math_h

#if defined (hex_win32)
#define _USE_MATH_DEFINES
#endif

// using namespace std;

#include "hex_sys.h"
#include "hex_cmplx.h"
#include "hex_gmp.h"
#include "hex_transform.h"

#include <cmath>
// include <math.h>

#include <algorithm>   // this provides min, max, etc 
#include <complex>

using namespace std;

#ifdef __cplusplus
  extern "C" {
#endif

typedef complex<float>       zfloat;
typedef complex<double>      zdouble;
typedef complex<long double> zzdouble;

typedef struct _cell2 {
   int  i;
   int  j;
} Cell2;
                                                                                
typedef struct _cell3 {
   int  i;
   int  j;
   int  k;
} Cell3;

typedef struct {               /* low-level structure for factorials */

   int     n_top;              /* the highest factor so far */
   int    *factors;            /* ordinary precision translation */

} sp_t_struct;

typedef sp_t_struct sp_t[1];


typedef struct {               /* low-level translation matrix structure */

   int     n_max;              /* order of translation matrix */
   double  r;                  /* requested distance "R" of the translation */
   double  r_cache;            /* actual distance "R" used in the cache */
   double *tm;                 /* ordinary precision translation */
   mpf_t  *mpf;                /* extended precision translation */

} tm_t_struct;

typedef tm_t_struct tm_t[1];

typedef struct {               /* more convenient to work two-together */

    tm_t gto;
    tm_t eto;

} HexTm_struct;

typedef HexTm_struct HexTm[1]; /* high-level translation matrix structure */

typedef struct {               /* similarly, basic rotation matrix structure */

    int     n_max;
    Euler   er;
    double *rlmj;

} HexRm_struct;

typedef HexRm_struct HexRm[1]; /* application-level access to the structure */


typedef struct {             /* canonical attributes and orientation of a molecule */

   int        n_max;         /* GTO expansion order at which calculated */
   Matrix3D   t_align;       /* the transformation (rotation) that aligns Q with axes */
   double     mpole;         /* the mass monopole moment (i.e. the "mass") */
   Vector3D   qpole;         /* principal quadrupole moments: X<Y<Z */
   double     moct[8];       /* octant masses in order of (+x,+y,+z), (-x,+y,+z), ... */
   Vector3D   radii;         /* approximate ellipsoidal radii */

} HexCanon;


#define HEX_GTO_TM   1
#define HEX_ETO_TM   0

#define HEX_RDEP     0
#define HEX_ZINV     1
#define HEX_RINV     2

/*---------------------------------------------------------------------------*/

#define sign(x)      (((x) >= 0) ? 1 : -1)
#define neg1(m)      (((m)/2)*2 == (m) ? 1 : -1)     /* (-1)^m in-lined      */
#define oddm(m)      (((m)/2)*2 == (m) ? 0 :  1)     /* odd(m) in-lined      */
#define evenm(m)     (((m)/2)*2 == (m) ? 1 :  0)     /* even(m) in-lined     */
#define pow2(m)      (1<<(m))                        /* 2^m, m positive      */
//define min(a,b)     (((a) < (b)) ? (a) : (b))       /* min(a,b)             */
//define max(a,b)     (((a) > (b)) ? (a) : (b))       /* max(a,b)             */
//define square(x)    ((x)*(x))                       /* x^2, any x           */
//define abs_(a)      (((a) > 0) ? (a) : -(a))        /* abs(a)               */
/*---------------------------------------------------------------------------*/

/*  Array indexing and dimensioning Macros ... */
/*  ------------------------------------------ */

/* circular functions: M, -L <= +L , trivial, but convenient */

#define ZM(M)       (MAX_L+(M))
#define ZM_DIM      (2*MAX_L+1) 

/* triangular array [l,m] for l>=m>=0, for Legendre polynomials, etc. */

#define PLM(L,M)    ((L)*((L)+1)/2+(M))
#define PLM_DIM(L)  (((L)+1)*((L)+2)/2)

/* Harmonic expansion coefficients, etc.: 2*L+1 triangle: -L<=M<=L */

#define YLM(L,M)    ((L)*(L)+(L)+(M))
#define YLM_DIM(L)  (((L)+1)*((L)+1))

/* (2*L+1)^2 square rotation matrix for a given L where -L<=J,M<=L */

#define DLJM(L,J,M) (((J)+(L))*(2*(L)+1)+(L)+(M))
#define DLJM_DIM(L) ((2*(L)+1)*(2*(L)+1))

/* total coefficient rotation matrices: note M<->J transposition from above */

#define ELMJ(L,M,J) ((L)*(4*(L)*(L)-1)/3 + ((M)+(L))*((L)+(L)+1) + ((J)+(L)))
#define ELMJ_DIM(L) (((L)+1)*(4*((L)+1)*((L)+1)-1)/3)

/* a[n,l,m] coefficient indexing: N > L and -L <= M <= +L, etc... */

#define NLM_DIM(N)  ((N)*((N)+1)*(2*(N)+1)/6) 
#define NLM(N,L,M)  (((N)-1)*(N)*(2*(N)-1)/6 + YLM(L,M))

/* triangular array [n,l] for n>l>=0 */

#define NL(N,L)     (((N)-1)*(N)/2+(L))
#define NL_DIM(N)   ((N)*((N)+1)/2)

/* unpacked translation matrix T[pq,kj,m] with positive m to order N */

#define TM_DIM(N)       (NL_DIM(N)*NL_DIM(N)*(N))
#define TM(N,P,Q,K,J,M) (NL_DIM(N)*(NL_DIM(N)*NL(P,Q)+NL(K,J))+(M))


/* square array lambda[m,n] for m=0...MAX_L, n=0...MAX_LAMBDA */

#define LAMBDA_MN(M,N) ((M)*MAX_LAMBDA+(N))
#define LAMBDA_DIM     (MAX_LAMBDA*L_DIM)


/* obsolete */

#define RNJ(R,J)    (R*N_DIM+J)

/* index/dimension the Legendre coefficients: mu(l,k,m) */
/* define MU(L,M,K) (L*(M*(M+1)/2+K)) */

/* its less memory-efficient to treat mu as a regular volume */

#define MU(L,M,K) (K + M*L_DIM + L*L_DIM*L_DIM)
#define MU_DIM    (L_DIM*L_DIM*L_DIM)

/* better index/dimension for the Legendre coefficients: D(l,m) */

#define CLMJ_DIM(L) (((L)+1)*((L)+2)*((L)+3)/6)
#define CLMJ(L,M,J) ((L)*((L)+1)*((L)+2)/6 + ((L)+1)*(M)-((M)*((M)-1))/2+(J))

/* index/dimension for the angular array A[l,l',|m|] where l>=l'>=m>=0 */

#define ALKM_DIM(L) (((L)+1)*((L)+2)*((L)+3)/6)
#define ALKM(L,K,M) ((L)*((L)+1)*((L)+2)/6 + (K)*((K)+1)/2 + (M))

/* index/dimension for the angular array A[l,l',|m|] (N should be n_top) */

#define ALJM_DIM(N)   ((N)*(N)*(N))
#define ALJM(N,L,J,M) ((N)*(N)*(L) + (N)*(J) + (M))

/* index/dimension for the Laguerre coefficients: c(n,l,k) */

/*
define CNLK_DIM(N) ((N)*((N)+1)*(2*(N)+10)/12)
define CNLK(N,L,K) (((N)-1)*(N)*(2*(N)+8)/12 + ((L)*(2*(N)-(L)+3))/2 + (K))
*/

#define CNLK_DIM(N) ((N)*((N)+1)*((N)+2)/6)
#define CNLK(N,L,K) (((N)-1)*(N)*((N)+1)/6 + ((L)*(2*(N)-(L)+1))/2 + (K))

/* generalised Laguerre polynomials have a different form */

/*
#define CNLM_DIM(N) ((N)*((N)+1)*(2*(N)+1)/6)
#define CNLM(N,L,M) (((N)-1)*(N)*(2*(N)-1)/6 + (N)*(L) + (M))

#define CPLM_DIM(P) ((((P)+1)*((P)+2)*(2*(P)+3))/6)
#define CPLM(P,L,M) (((P)*((P)+1)*(2*(P)+1))/6 + ((P)+1)*(L) + (M))
*/

/*
define CPQK_DIM(P) ((((P)+1)*((P)*(P)+5*(P)+6))/6)
define CPQK(P,Q,K) (((P)*((P)*(P)+3*(P)+2))/6+(Q)*((P)+1)-((Q)*((Q)-1))/2+(K))
*/

/* binomial coefficients b[n,m] where m <= n, etc. */

#define BNM(N,M)    ((N)*((N)+1)/2 + (M))
#define BNM_DIM(N)  (((N)+1)*((N)+2)/2)

/* in-lined version of idx_ijk() - 3d array indexing in hex_math.c */

#define IDX_IJK(I,J,K,NI,NJ)  ((I) +(J)*(NI) + (K)*(NI*NJ))

/*---------------------------------------------------------------------------*/

void    hex_cache_gto_scale (int q);
void    hex_cache_eto_scale (int l);
void    hex_cache_probe     ();
void    hex_cache_mode      (int mode);
int     hex_cache_type      ();
void    hex_write_cache     (int  n_max, int  lr12, int jr12, int type_flag, 
                             int  nc, double *coeffs);
int     hex_read_cache      (int  n_max, int  lr12, int jr12, int type_flag, 
                             int  nc, double *coeffs);
char   *hex_cache_path      ();
int     hex_cache_create    (char *filename);
int     hex_cache_open      (char *filename);
void    hex_cache_close     (int fildes);

HexFile *hex_cache_fcreate  (char *filename);
HexFile *hex_cache_fopen    (char *filename);
void     hex_cache_fclose   (HexFile *file);

         /* K(R) integral stuff - mostly OBSOLETE */

void    hex_kr_increments   (int  n_steps);

void    hex_print_kr        (int  n_max, int  n_low, int  n_top, double *kr);

void    hex_make_kr         (int  n_max, int  ir12, int  ir0);

double *hex_get_kr          (int  n_max, double r12, int  ir, int *status);

void    hex_free_kr         (double **kr);

void    hex_unit_kr         (int  n_max, double *kr);
void    hex_calc_kr         (int  n_min, int  n_max, double r12, int  ir, 
                             double *kr);

void    hex_calc_jr1        (int  n_min, int  n_max, double r12, double *jr);
void    hex_calc_jr2        (int  n_min, int  n_max, double r12, double *jr);

void hex_kr_anlm            (int  n_max, double *kr,
                             double *anlm, double *bnlm);

int    quad_roots    (double a, double b, double c, double *x1, double *x2);

int    kronecker     (int, int);             /* Kronecker delta function */

double gam           (double a);             /* Euler's Gamma fn for real arg*/
double igam          (double a, double x);   /* incomplete gamma function    */
double jgam          (double a, double x);   /* incomplete Gamma function    */

double bet           (double a, double b);             /* Beta function */
double ibet          (double a, double b, double x);   /* incomplete Beta */
double jbet          (double a, double b, double x);   /* incomplete Beta */

double ak_beta       (int k, double beta);   /* special case gamma functions */
double bk_beta       (int k, double beta);
double ck_beta       (int k, double beta);

double fac           (int  n);                       /* factorial(n)         */
double ffac          (int  n, double x);             /* falling factorial    */
double rfac          (int  n, double x);             /* rising factorial     */
double rfac2         (int  n);                       /* i.e. rfac(n,0.5)     */
double gbcf          (double r, int  k);             /* binomial coeff (r,k) */
double gam2          (int  n);                       /* Gamma(n+1/2)         */
double g2mg2         (int  n, int  m);               /* gam2(n)*gam2(m)      */
double g2dg2         (int  n, int  m);               /* gam2(n)/gam2(m)      */
double sqrtn         (int  n);                       /* cached sqrt(n)       */

hexa   hfac          (int  n);                       /* high precision       */
hexa   hrfac         (int  n, hexa x);               /* versions of the above*/
hexa   hffac         (int  n, hexa x);
hexa   hsqrt         (hexa x);
hexa   hsqrtn        (int  n);

void   sp_init       (sp_t sp);
void   sp_fac        (sp_t sp, int op, int  n);               /* symbolic prime factor*/
void   sp_fall       (sp_t sp, int op, int  n, int  k);       /* method of evaluating */
void   sp_rise       (sp_t sp, int op, int  n, int  k);       /* products of primes   */
void   sp_rise2      (sp_t sp, int op, int  k);
void   sp_fall2      (sp_t sp, int op, int  n, int  k);
void   sp_pow        (sp_t sp, int op, int  n, int  k);
hexa   sp_ans        (sp_t sp);
void   sp_free       (sp_t sp);

int    neg1_m        (int m);                         /* (-1)^m               */
int    odd_m         (int m);                         /* true if m is odd     */
int    even_m        (int m);                         /* true if m is even    */
int    roundup_m     (int  m, int modulus);           /* round up to multiple */
int    roundup_bank  (int  m, int modulus, int bank); /* round up with memory bank */

double  h1f1         (double a, double b, double x); /* hypergeometric fn 1F1*/
double  h1f1_ser     (double a, double b, double x); /* by series expansion  */
double  h1f1_rk      (double a, double b, double x); /* by Runge-kutta */
Complex h1f1z        (double a, double b, Complex z, /* complex series expn  */
                      Complex *dz);

int    sn1           (int  n);             /* sum from 0 to n of n^1 */
int    sn2           (int  n);             /* sum from 0 to n of n^2 */
int    sn3           (int  n);             /* sum from 0 to n of n^3 */
int    sn4           (int  n);             /* sum from 0 to n of n^4 */

// void   hex_fft_alloc (int npow);           /* 1D complex FFT stuff */
// void   hex_fft_mode  (int mode);
// void   hex_fft_fs    (int nmax, double *am, double *bm, double *data);
// void   hex_fft       (int direction, double *fr, double *fi);
// void   hex_fft_order (int np, int *sort);
// void   hex_fft_free  ();
// void   hex_fht_fs_thr(int n_max, double *am, double *bm, double *real, double* lwkr);


void   hex_fft_1d    (int np,                 int isign, double *data);
void   hex_fft_2d    (int px, int py,         int isign, double *data);
void   hex_fft_3d    (int px, int py, int pz, int isign, double *data);
void   hex_fft_nd    (int ndim, int *npow, int isign, double *data);
void   hex_fft_ref   (int ndim, int *npow, int isign, double *data);
void   hex_fft_test  ();

double row_x_col     (int, double [], double []);           /* matrix math */
void   col_x_row     (int, double [], double [], double []);
void   mat_x_mat     (int, double [], double [], double []);
void   mat_x_col     (int, double [], double [], double []);
void   val_x_vec     (int n, double val, double a[], double b[]);

void   mat_ident     (int, double []);
double mat_trace     (int, double []);
void   mat_transpose (int, double [], double []);
void   mat_inverse   (int, double [], double []);
void   mat_lsq       (int n, int m, double *y, double *f, double *wt,
                      double *a);
void   mat_eigen_jr  (int n, double *a, double *e, double *mu);
void   mat_eigen_ql  (int n, double *a, double *e, double *mu);
void   mat_eigen_qlip(int n, double *a, double *mu);
void   hex_eigen_tol (double tol);
void   mat_eigen_sort(int n, int with_sign, int ascending, 
                      double *e,      double *mu, 
                      double *e_sort, double *mu_sort);

void   mat_cycle_3d  (int *dims_in, int *dims_out, float *cba, float *acb);

hexa   hrow_x_col    (int, hexa [], hexa []);        /* quad precision */
void   hcol_x_row    (int, hexa [], hexa [], hexa []);
void   hmat_x_mat    (int, hexa [], hexa [], hexa []);
void   hmat_x_col    (int, hexa [], hexa [], hexa []);

void   hmat_ident    (int, hexa []);
hexa   hmat_trace    (int, hexa []);
void   hmat_transpose(int, hexa [], hexa []);
void   hmat_inverse  (int, hexa [], hexa []);
void   hmat_lsq      (int n, int m, double *y, double *f, double *wt,
                      hexa *a);

double rot_to_axis   (double [], double []);

void   vec_scale     (int n, double a[], double scale, double out[]);
void   vec_accum     (int n, double a[], double out[]);
void   vec_add       (int n, double a[], double scale, double b[], double out[]);
double vec_distance  (int n, double a[], double b[]);
double vec_mag       (int n, double a[]);
void   vec_copy      (int, double [], double []);

void   vec_cross3    (double a[], double b[], double c[]);

void   vec_p_vec     (int, double [], double [], double []);
void   vec_m_vec     (int, double [], double [], double []);

void   vec_print     (int, double [], rchar []);
void   mat_print     (int, double [], rchar []);

double mod_angle     (double, double);

int    degenerate_rotn(double, double, double,
                       double, double, double, double);

int    degenerate_axes(Axis, Axis, double);

double hex_bfgs   (double (*fn)(int nvar, double *vars, void *fn_data, int it),
                   int    nvar,         /* no. variables in x[] */
                   int    ntries,       /* no. iterations */
                   void   *fn_data,  /* address of other user data */
                   double xtol,      /* convergence tolerance for |x| */
                   double xstart[],  /* starting guess for variable vector x */
                   double xmin[],    /* vector of lower bounds for x */
                   double xmax[],    /* vector of upper bounds for x */
                   double xstep[],   /* initial guess at step for x */
                   double xans[]);   /* solution vector of x values */


HexCanon hex_gto_canonical(int n_max, double *anlm, int use_quadrupole);

Matrix3D tf_fitAlm  (int  lmax, double *moving_alm, double *fixed_blm);
Matrix3D tf_fitAnlm (int  nmax, double *moving_anlm, double *fixed_bnlm);
Matrix3D tf_fitRot  (int  npts, Point3D *moving_pts, Point3D *fixed_pts);
Matrix3D tf_fitPts  (int  npts, Point3D *moving_pts, Point3D *fixed_pts);
Matrix3D tf_fitPtsWts(int  npts, Point3D *moving_pts, Point3D *fixed_pts,
                      double *weights);

Matrix3D tf_fitQCP   (int npts, Point3D *moving_pts, Point3D *fixed_pts);
Matrix3D tf_fitQCPwts(int npts, Point3D *moving_pts, Point3D *fixed_pts, double *wts);


Matrix3D tf_fitpts_old(int  npts, Point3D *moving_pts, Point3D *fixed_pts);

Plane pl_pts_pt    (int npts, Point3D *pts, Point3D x0);


double rms_d2pts   (int npts,             double  *d2_pts);
double rms_pts     (int npts,             Point3D *moving_pts, Point3D *fixed_pts);
double rms_tf_pts  (int npts, Matrix3D t, Point3D *moving_pts, Point3D *fixed_pts);
void   tf_pts      (int npts, Matrix3D t, Point3D *moving_pts, Point3D *fixed_pts);

void polar_ligand_tf(Matrix3D t,
                     double *r12,    double *beta1, double *gamma1,
                     double *alpha2, double *beta2, double *gamma2);

Cell2 c2_make      (int  i, int  j);
Cell2 c2_add       (Cell2 c1, Cell2 c2);
Cell2 c2_idx       (int  idx, int  ni);
Cell3 c3_make      (int  i, int  j, int  k);
Cell3 c3_byte      (byte b);
Cell3 c3_add       (Cell3 c1, Cell3 c2);
Cell3 c3_idx       (int  idx, int  ni, int  nj);

int   idx_c2       (Cell2 c, int  ni);
int   idx_c3       (Cell3 c, int  ni, int  nj);

double hex_r12_did   (int did);
int    hex_did_r12   (double r12);
double hex_r12_did_max();

int     idx_ijk    (int  i, int  j, int  k, int  ni, int  nj);
void    ijk_idx    (int  idx, int  ni, int  nj, int  *i, int  *j, int  *k);

void    ylm_mode     (int);                          /* direct/recursive     */
Complex cylm         (int, int, double, double);     /* Y[l,m](theta,phi)    */
Complex eim          (int m, double phi);            /* e^{im.phi}/sqrt(2PI) */

void    dljm_all     (int l_max, double beta, double *dbuf); /* d(l)[m',m](theta) */

double  dljm         (int l, int j, int m, double);  /* d(l)[m',m](theta)    */
double  dljm0        (int l, int j, int m, double);  /* d(l)[m',m](theta)    */
double  dljm1        (int l, int j, int m, double);  /* d(l)[m',m](theta)    */
double  delta_ljm    (int l, int j, int m);          /* d(l)[m',m](PI/2)     */

double  ylm          (int, int, double, double);     /* y[l,m](theta,phi)    */
double  eljm         (int l,int j,int m,       /* E(l)[m',m](alpha,beta,gam) */
                      double alpha, double beta,double gamma);
double eljm_trace    (int l,                          /* trace of E(l)(abg)  */
                      double alpha, double beta, double gamma);
double eljm_subtrace (int l, int j,                   /* trace of E(l)(abg)  */
                      double alpha, double beta, double gamma);
Axis   eljm_axis     (int l,
                      double alpha, double beta, double gamma);

double  zm           (int m, double phi);             /* Z[m](phi)           */
double  zm_bar       (int m, double phi);             /* unnormalised zm()   */
double  zm_norm      (int m);                         /* normalisation factor*/
void    zm_all       (int l_max, double phi, double *zfn);
double  dzm          (int m, double);                /* d/d(phi) Z[m](phi)   */
double  izm          (int m, double, double);        /* Intg Z[m].d(phi)     */

double  qlm          (int, int, double);             /* (-1)^m.P[l,m](theta) */
double  plm          (int, int, double);             /* P[l,m](theta)        */
double  plm_norm     (int l, int m);
double  ylm_norm     (int l, int m);
double  iplm         (int l, int m, double th1, double th2);        /* Intg p[l,m].d(theta) */
double  ipklm        (int k, int l, int m, double th1, double th2); /* Intg p[l,m].d(theta) */

double  lna          (int n, double alpha, double r); /* normalised Laguerre */
double  lna_ser      (int n, double alpha, double r); /* normalised Laguerre */
double  lna_rec      (int n, double alpha, double r); /* normalised Laguerre */

void    lna_all_bar  (int n, int mult, double beta,   /* unnormalised array  */
                      double r, double *lfn);

void    hlna_all_bar (int n, int mult, hexa beta,     /* quad precision */
                      hexa r, hexa *lfn);

void    lna_all_vec  (int n, double alpha, double r,  /* unnormalised vector */
                      double *v);                     /* (preferred method)  */
void    hlna_all_vec (int n, hexa alpha, hexa r,      /* extra precision     */
                      hexa *v);

double  lna_2l_two_norm (int n, int l);                     /* for 2l+2 form */

double  lna_1l_half     (int n, int l, double r);          /* for l+1/2 form */
double  lna_1l_half_bar (int n, int l, double r);
double  lna_1l_half_norm(int n, int l);
hexa    hlna_1l_half_norm(int n, int l);

double  hnl          (int n, int l, double r);      /* harmonic radial basis */
double  hnl_norm     (int n, int l); 
hexa    hhnl_norm    (int n, int l); 

double  hnl_monopole  (int n_max, double *anlm);
void    hnl_dipole    (int n_max, double *anlm, double *vec);
void    hnl_quadrupole(int n_max, double *anlm, double *mat);
void    hnl_inertia   (int n_max, double *anlm, double *mat);
void    hnl_hemisphere(int n_max, double *anlm, double *vec);
void    hnl_octants   (int n_max, double *anlm, double *vec);
void    hnl_radial    (int n_max, double *anlm, double *vec);


double  pnl          (int n, int l, double r);       /* Jacobi radial basis */
double  pnl_norm     (int n, int l); 

double  tnl          (int n, int l, double r);       /* Coulomb radial basis */
double  tnl_norm     (int n, int l); 

double  dtnl         (int n, int l, double r);       /* d/dr*T[n,l](r)       */
double  itnl         (int n, int l, double r1, double r2); /* Intg T[n,l](r) */

double  jlr_test     (int l, double r);              /* j[l](r) (sph Bessel) */
double  jlr          (int l, double r);              /* j[l](r) (sph Bessel) */
double  jlr_ser      (int l, double r);              /* by series expansion  */
void    jlr_vec      (int l, double r, double *vec); /* by recursion         */

double  knr          (int n, double r);              /* K[n+1/2](r) Bessel fn*/
double  rnr          (int n, double r);              /* "reduced" K[n+1/2] fn*/
void    rnr_all      (int n, double r, double *rfn);
void    hrnr_all     (int n, double r, hexa   *rfn);
void    hrnr_all_bar (int n, double r, hexa   *rfn);

hexa   *clj_jlr      (int l_max);                    /* supply j[l](r) coeffs*/

hexa    hjlr_ser     (int l_max, hexa r);            /* quad precision  ... */
void    hjlr_vec     (int l_max, hexa r, hexa *jfn); 

void    pnl_all_bar  (int n_max, double r,     double *pfn);
void    tnl_all_bar  (int n_max, double r,     double *tfn); 
void    hnl_all_bar  (int n_max, double r,     double *hfn); 
void    plm_all_bar  (int l_max, double theta, double *pfn);
void    zm_all_bar   (int l_max, double phi,   double *zfn);
void    qnl_all_bar  (int n_max, double r,     double *qfn); 
void    wnl_all_bar  (int n_max, double r,     double *wfn); 
void    jlr_all      (int l_max, double r,     double *jfn);

void    hhnl_all_bar (int n_max, hexa r,       hexa   *hfn); 
void    hjlr_all     (int l_max, hexa r,       hexa   *jfn);

void    ylm_all_bar  (int l_max, double theta, double phi, double *yfn);
void    ylm_all_norm (int l_max,                           double *yfn);

void    hn_all_bar   (int n_max, double r,  double s, double *hfn); 
void    hn_all_bar2  (int n_max, double r2, double s, double *hfn); 
double  hn_norm      (int n, double s);

double pab_reference(int n, double alpha, double beta, double x);
double pab_shift(int n, double alpha, double beta, double x);

double  pab          (int n, double alpha, double beta, /* unnormalised */
                      double x);
double  pab_norm     (int n, double alpha, double beta);
double  jab          (int n, double alpha, double beta, /* normalised */
                      double x);
double  pkl_ser      (int k, int l, double x);
void    pkl_all_vec  (int k_max, int l, double x, double *vec);

double  tlm          (int  l, int m,  double dz);    /* translation operator */

double  wcc_mp       (int  j1, int  j2, int  j,      /* Wigner's coupling    */
                      int  m1, int  m2, int  m);     /* coefficient          */

hexa    wcc          (int  j1, int  j2, int  j,      /* Wigner's coupling    */
                      int  m1, int  m2, int  m);     /* coefficient          */
void    wcc_all      (int  j1, int  j2, int  j, double *cfn);
void    wcc_up       (int  j1, int  j2, int  j,
                      int  ma, int  mb, double *cfn);
void    wcc_down     (int  j1, int  j2, int  j,
                      int  ma, int  mb, double *cfn);

hexa    w3j          (int  j1, int  j2, int  j,      /* Wigner 3-j coeffcient*/
                      int  m1, int  m2, int  m);
double  t3j          (int  j1, int  j2, int  j,      /* Wigner 3-j coeffcient*/
                      int  m1, int  m2, int  m);

double *dkj_integrals(int  n_max, double s, int  *j_stride);
double *dkj_mp       (int  n_max, double s, int  *j_stride);

double *dkjn_mp      (int  n_max, double s, int  *n_stride, int  *jn_stride);


void    rnl_factor   (double factor);                /* set scale for rnl()  */
void    pnl_factor   (double factor);                /* set scale for pnl()  */
void    tnl_factor   (double factor);                /* set scale for tnl()  */
void    hnl_factor   (double factor);                /* set scale for hnl()  */
void    wnl_factor   (double factor);                /* set scale for wnl()  */

double  hnl_scale    ();
double  tnl_scale    ();
double  pnl_scale    ();
double  wnl_scale    ();

double  lpq          (int n, int l, double r);       /* L[p,q](r)            */

double  rnl          (int n, int l, double r);       /* R[n,l](r)            */
double  drnl         (int n, int l, double r);       /* d/dr*R[n,l](r)       */
double  irnl         (int n, int l, double r1, double r2); /* Intg R[n,l](r) */

double  snl          (int n, int l, double beta);    /* S[n,l](beta) integral*/

double  qnl          (int n, int l, double r);       /* Q[n,l](r) integral   */

double  tnjl         (int n, int nn, int l);         /* T[n,n'l] integral    */
double  hnjl         (int n, int nn, int l);         /* H[n,n'l] integral    */

double  glk          (int l, int k);

double  jr_int       (int n1, int n2,                /* electrostatic overlap*/
                      int l1, int l2, int m, double r);

double  rnl_rho      (int n, int l, double r);       /* R[n,l](r,r0)         */
double  rhonl        (int n, int l, double rho);     /* R[n,l](rho)          */

double  tn           (int n, double);                /* T[n](theta)          */

double  jlx          (int l, double x);              /* J[l](x) Bessel fn    */

void    alm_zero     (int order,double[]);           /* zero-out coeffs      */
double  alm_value    (int,double,double,double[]);   /* evaluate expansions  */

void    anlm_zero    (int order,double[]);           /* 3D zero vector */
void    anlm_copy    (int order,double[], double[]); /* copy a vector */
void    anlm_unit    (int order,double[]);           /* 3D unit vector */
double  anlm_value   (int,double,double,double,double[]);
double  anlm_value_t (int,double,double,double,double[]);
double  anlm_value_g (int,double,double,double,double[]);
double  anlm_value_h (int,double,double,double,double[]);

void    anlm_print   (int, double[]);
void    anlm_check   (int n_max, double *anlm, char *where);

void    anlm_n00     (int n_max, hexa *an00, double *anlm);

void    rp_nlm       (int, double r,                /* R[n,l](r)P[l,m](theta)*/
                      double theta,double fnlm[]);

int     tm_mode      (int method);
void    tm_bits      (int n_bits);                  /* set bits of precision */
int     tm_top       ();
void    tm_init      (tm_t t);
void    tm_make      (int  n_max, double r_max);
void    tm_get       (int  n_max, int matrix_type, double r, tm_t t);
void    tm_free      (tm_t t);

void    tm_gto_ser   (int  n_max, double dz, double *tm);
void    tm_eto_ser   (int  n_max, double dz, double *tm);
void    tm_gto_num   (int  n_max, double dz, double *tm);
void    tm_eto_num   (int  n_max, double dz, double *tm);

void    ztm_mult_gto (int  n_max, tm_t t, zdouble *zanlm, zdouble *zanlm_t);
void    ztm_mult_eto (int  n_max, tm_t t, zdouble *zanlm, zdouble *zanlm_t);
void    ztm_mult     (int  n_max, int matrix_type, tm_t t, zdouble *zanlm, 
                      zdouble *zanlm_t);


void    tm_copy      (int n_max, tm_t t, double *tm);

int     rcc_dim      (int  n_max); /* supply radial coupling coeff dimension */
int     cc_dim       (int  n_max); /* alternative coupling coeff dimension */
int     ac_dim       (int  n_max); /* angular coefficient dimension */

int     tv_dim       (int  n_max);                  /* supply V(R) dimension */
int     tm_dim       (int  n_max);                  /* supply T(R) dimension */
int     tm_top       ();                            /* max N for this session*/
double *tm_eto       (int  n_max, double dz);       /* calc Coulomb T(R)     */
double *tm_gto       (int  n_max, double dz);       /* calc T(R) in HO basis */
hexa   *tz_gto       (int  n_max, double dz);       /* quad precision HO-Z(R)*/
hexa   *tu_gto       (int  n_max, double dz);       /* quad precision HO-U(R)*/
double *tm_gto_power (int  n_max, double dz);       /* quad precision HO-T(R)*/
double *tm_eto_power (int  n_max, double dz);       /* quad precision W-T(R)*/
hexa   *tz_gto_power (int  n_max, double dz);       /* quad precision HO-Z(R)*/
hexa   *tu_gto_power (int  n_max, double dz);       /* quad precision HO-U(R)*/
double *tm_gto_hexa  (int  n_max, double dz);       /* quad precision HO-T(R)*/

mpf_t  *tv_gto       (int  n_max, double dz);       /* steric  V(R) overlaps */
mpf_t  *tw_eto       (int  n_max, double dz);       /* electro W(R) overlaps */
mpf_t  *ts_gto       (int  n_max, double dz);       /* generic V(R) overlaps */
mpf_t  *ts_eto       (int  n_max, double dz);       /* generic W(R) overlaps */


void    tm_neg       (int  n_max, double *tm);      /* reverse a translation */
void    tv_neg       (int  n_max, double *tm);      /* reverse a translation */
void    tm_neg_mpf   (int  n_max, mpf_t  *tm);      /* reverse a translation */

void tm_ortho_test   (int  n_max, double *tm);
void tm_resol_test   (int  n_max, double *tm_lo, double *tm_hi);
void tm_unpack_m     (int  n_max, int  m, double *tm, double *mat);

void    tm2mqp       (int n_max, int m_dim, int p_dim, /* unpack a tm array */
                      int *p_order,                 /* to sorted mqp format */
                      double *tm, float *tmqp);        

void    tm_ident     (int  n_max, double *t);
void    tm_product   (int  n_max, double *t1, double *t2, double *t);
int     tm_sign      (int  n_max, int  n1, int  l1, int  n2, int  l2);
int    *tm_table     (int  n_max);
int     tm_index     (int  n_max, int  n1, int  l1, int  n2, int  l2);
double  tm_value     (int  n_max, int  n1, int  l1, int  n2, int  l2, 
                      int  m, double *tm);

void    tm_apply     (int  n_max, double *tm,       /* translate coeff vector*/
                      double *anlm, double *tnlm);  /* by applying T(R) to it*/

void    ztm_apply    (int  n_max, double *tm,       /* translate coeff vector*/
                      zdouble *zanlm, zdouble *ztnlm); 

void ztm_apply_vec   (int n_max, int n_vec, int *vec_types, HexTm htm,
                      zdouble *zanlm, zdouble *zanlm_t);


void    tm_apply_mpf (int  n_max, mpf_t *tm,        /* same thing ... */
                      mpf_t *anlm, mpf_t *tnlm);

void    tv_apply     (int  n_max, mpf_t *tm,        /* translate coeff vector*/
                      double *anlm, double *tnlm);  /* by applying T(R) to it*/

void    tw_apply     (int  n_max, mpf_t *tm,         /* translate coeff vector*/
                      double *anlm, double *tnlm);  /* by applying T(R) to it*/

void    tz_apply     (int  n_max, hexa *tz,         /* apply reduced T(R) to */
                      hexa *an00, double *anlm);    /* a reduced an00 vector */

void   tm_mult       (int  n_max, int matrix_type, tm_t t, double *anlm,
                      double *anlm_t);

void    tm_print     (int  n_min, int  n_max, double *tm);
void    tm_print_all (int  n_min, int  n_max, double *tm);

void    anlm_to_ajlm_gto (int n_max, mpf_t *anlm, mpf_t *ajlm);
void    ajlm_to_anlm_gto (int n_max, mpf_t *ajlm, mpf_t *anlm);

void    anlm_to_ajlm_eto (int n_max, mpf_t *anlm, mpf_t *ajlm);
void    ajlm_to_anlm_eto (int n_max, mpf_t *ajlm, mpf_t *anlm);

void hex_tm_init     (HexTm htm);
void hex_ztm_apply   (int n_max, int n_vec, int *vec_types, 
                      double r, HexTm htm, zdouble *zanlm, zdouble *zbnlm);
void hex_tm_apply    (int n_max, int n_vec, int *vec_types, 
                      double r, HexTm htm, double *anlm, double *bnlm);
void hex_trm_apply   (int n_max, int n_vec, int *vec_types, double r, Euler er, 
                      HexTm htm, HexRm hrm, double *anlm, double *bnlm);
double *hex_tm_get   (int n_max, int vec_type, double r, HexTm htm);
void hex_tm_free     (HexTm htm);

void hex_rm_init     (HexRm hrm);
void hex_rm_apply    (int n_max, int n_vec, Euler er, HexRm hrm, 
                      double *anlm, double *bnlm);
double *hex_rm_get   (int n_max, Euler er, HexRm hrm);
void hex_rm_free     (HexRm hrm);


double  ylm_potential(int,double,                   /* potential expansion   */
                      double,double,double,double[]);

double  ylm_euler    (int,int,double,double,    /* rotate y[l,m](theta,phi)  */
                      double,double,double);    /* by (alpha,beta,gamma)     */

void    alm_rotate   (int,                           /* rotate coefficients  */
                      double alpha, double beta, double gamma,
                      double *alm, double *blm);

void    almkj_rotate (int l_max,                    /* rotate density matrix */
                      double alpha, double beta, double gamma,
                      double *almkj, double *blmkj);

double  anlm_stupid  (int n_vec, int n_max, double *anlm, double *bnlm);

void    anlm_rinv    (int n_max, int n_vec, int rinv, 
                      double *anlm, double *anlm_rinv);

double  anlm_carbo_rinv(int n_max, int n_vec, int rinv, 
                      double *anlm, double *anlm_rinv);

void    anlm_flip    (int n_vec, int n_max, double *anlm, double *anlm_flip);

void    anlm_rotate  (int n_max, int n_vec,          /* rotate coefficients */
                      double alpha, double beta, double gamma,
                      double *anlm, double *bnlm);

void    anlm_rm_apply(int n_max, int n_vec,          /* rotate coefficients */
                      Euler er, double *anlm, double *bnlm);

void    anlm_rot_phi (int  n_max,                    /* specialised rotation */
                      double *cos_mphi, double *sin_mphi,
                      double *anlm, double *bnlm);

double  anlm_sub_phi (int  n_max, double *anlm, double *bnlm,
                      double *cos_mphi, double *sin_mphi);

double  alm_dot      (int  l_max,                    /* sum of products      */
                      double *alm, double *blm);

double  alm_mag      (int  l_max, double *alm);      /* magnitude of vector  */

double  anlm_cross_sum(int  n_max, double *arho, double *atau,
                                   double *brho, double *btau);

void    anlm_accum   (int n_max, double *anlm, double *bnlm);
void    anlm_accum_n (int n_max, int n_nec, double *anlm, double *bnlm);

double  anlm_mag     (int  n_max, double *anlm);

double  anlm_dot     (int  n_max, double *anlm, double *bnlm);

double  anlm_distance(int  n_max, double *anlm, double *bnlm);

double  anlm_carbo   (int  n_max, double *anlm, double *bnlm);

void    anlm_add     (int  n_max, double *anlm, double *bnlm, double factor,
                      double *snlm);

void    anlm_scale   (int  n_max, double *anlm, double factor, double *snlm);
void    anlm_scale_n (int n_vec, int  n_max, double *anlm, double factor, double *snlm);

void    anlm_mult    (int  n_max, double *anlm, double *bnlm, double *cnlm);

void    anlm_mult_bar(int  n_max, double *anlm, double *bnlm, double *cnlm);

void    anlm_lambda  (int  n_max, double *anlm, double *lamdba);

double  lambda_sum   (int  m_max, double *lambda1, double *lambda2);

void    lambda_fs    (int  m_max, double *lambda1, double *lambda2,
                      double *fs1, double *fs2);

void anlm_cross_e    (double e);
void anlm_cross_f    (double f);
void anlm_cross_g    (double g);
void anlm_cross_p    (double p);
void anlm_cross_q    (double q);
void anlm_cross_s    (double s);
void anlm_cross_z    (double z);

void cylm_addn       (int l_max, Point3D x1, Point3D x2);
void  ylm_addn       (int l_max, Point3D x1, Point3D x2);
void  jrl_addn       (int l_max, Point3D x1, Point3D x2);

double plm_test      (int, int);                     /* test integration     */
double plm_test_fac  (int, int);                     /* same with factorials */
double plm_r         (int, int, double);             /* recursive P[l,m]     */

void   set_phi_cache (int max_m, double phi, double **cm, double **sm);

void   set_trig_cache(int);                          /* set cache size       */
double cos_theta_n   (double, int);                  /* cos(theta)^n         */
double sin_theta_n   (double, int);                  /* sin(theta)^n         */
double cos_m_phi     (int, double);                  /* cos(m*phi)           */
double sin_m_phi     (int, double);                  /* sin(m*phi)           */
double cos_m_alpha   (int, double);                  /* cos(m*alpha)         */
double sin_m_alpha   (int, double);                  /* sin(m*alpha)         */
double cos_m_beta    (int, double);                  /* cos(m*beta)          */
double sin_m_beta    (int, double);                  /* sin(m*beta)          */
double cos_m_gamma   (int, double);                  /* cos(m*gamma)         */
double sin_m_gamma   (int, double);                  /* sin(m*gamma)         */
double cos_beta2_n   (double, int);                  /* cos(beta/2)^n        */
double sin_beta2_n   (double, int);                  /* sin(beta/2)^n        */
double cos_delta_n   (double, int);                  /* cos(delta)^n         */
double sin_delta_n   (double, int);                  /* sin(delta)^n         */

void   hex_d2f       (double *d, int n, float *f);
void   hex_f2d       (float  *f, int n, double *d);


/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif  /* hex_math_h */

/*---------------------------------------------------------------------------*/
