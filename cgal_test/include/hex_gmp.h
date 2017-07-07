/*-----------------------------------------------------------------------------
**  File:       hex_gmp.h
**
**  Author:     Dave Ritchie, 11/07/00
**
**  Purpose:    header file for various GNU MPF high precision functions...
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) : D.W. Ritchie, University of Aberdeen.
**  Copyright (C) 2010 D.W. Ritchie, INRIA.
**
**  This software was developed at the University of Aberdeen 1996-2000, by
**  Dave Ritchie, as part of a project funded by the BBSRC.
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

#ifndef hex_gmp_h
#define hex_gmp_h

//  include "mkl_gmp.h"  // Intel MKL only proivides MPZ, not MPF. No good.

#if defined(hex_mpir)
#   include "mpir.h"
#else
#   include "gmp.h"
#endif

typedef struct {               /* low-level structure for factorials */

   int     n_top;              /* the highest factor so far */
   int    *factors;            /* ordinary precision translation */

} mp_t_struct;

typedef mp_t_struct mp_t[1];

// use MPFR for multiple precision floating point arithmetic
// and make the MPFR constants and special functions look like MPF functions

#if defined(hex_mpfr)
#include "mpfr.h"
#include "mpf2mpfr.h"

#define mpf_pi    (ROP)           mpfr_const_pi   (ROP, MPFR_DEFAULT_RND)
#define mpf_euler (ROP)           mpfr_const_euler(ROP, MPFR_DEFAULT_RND)

#define mpf_gamma (ROP, OP)       mpfr_gamma  (ROP, OP, MPFR_DEFAULT_RND)
#define mpf_zeta  (ROP, OP)       mpfr_zeta   (ROP, OP, MPFR_DEFAULT_RND)

#define mpf_log   (ROP, OP)       mpfr_log    (ROP, OP, MPFR_DEFAULT_RND)
#define mpf_log2  (ROP, OP)       mpfr_log2   (ROP, OP, MPFR_DEFAULT_RND)
#define mpf_log10 (ROP, OP)       mpfr_log10  (ROP, OP, MPFR_DEFAULT_RND)

#define mpf_exp   (ROP, OP)       mpfr_exp    (ROP, OP, MPFR_DEFAULT_RND)
#define mpf_exp2  (ROP, OP)       mpfr_exp2   (ROP, OP, MPFR_DEFAULT_RND)
#define mpf_exp10 (ROP, OP)       mpfr_exp10  (ROP, OP, MPFR_DEFAULT_RND)

#define mpf_cos   (ROP, OP)       mpfr_cos    (ROP, OP, MPFR_DEFAULT_RND)
#define mpf_sin   (ROP, OP)       mpfr_sin    (ROP, OP, MPFR_DEFAULT_RND)
#define mpf_tan   (ROP, OP)       mpfr_tan    (ROP, OP, MPFR_DEFAULT_RND)

#define mpf_acos  (ROP, OP)       mpfr_acos   (ROP, OP, MPFR_DEFAULT_RND)
#define mpf_asin  (ROP, OP)       mpfr_asin   (ROP, OP, MPFR_DEFAULT_RND)
#define mpf_atan  (ROP, OP)       mpfr_atan   (ROP, OP, MPFR_DEFAULT_RND)

#define mpf_acosh (ROP, OP)       mpfr_acosh  (ROP, OP, MPFR_DEFAULT_RND)
#define mpf_asinh (ROP, OP)       mpfr_asinh  (ROP, OP, MPFR_DEFAULT_RND)
#define mpf_atanh (ROP, OP)       mpfr_atanh  (ROP, OP, MPFR_DEFAULT_RND)

#endif

/*---------------------------------------------------------------------------*/

int    mp_bits       (int  bits);          /* set no. bits precision */
char  *mp_version    ();

mpf_t *hex_mp_get  (int  n);               /* allocate an array of mpf_t's */
void   hex_mp_zero (mpf_t *array, int  n); /* zero out an array of mpf_t's */
void   hex_mp_free (mpf_t *array, int  n); /* free an array of mpf_t's */

mpq_t *hex_mq_get  (int  n);               /* allocate an array of mpq_t's */
void   hex_mq_zero (mpq_t *array, int  n); /* zero out an array of mpq_t's */
void   hex_mq_free (mpq_t *array, int  n); /* free an array of mpq_t's */

#define hex_mp_wipe(x,n) { if (x != NULL) { hex_mp_free(x,n); x = NULL; }}
#define hex_mq_wipe(x,n) { if (x != NULL) { hex_mq_free(x,n); x = NULL; }}

void   mp_init     (mp_t mp);
void   mp_fac      (mp_t mp, int op, int  n);
void   mp_fall     (mp_t mp, int op, int  n, int  k);
void   mp_rise     (mp_t mp, int op, int  n, int  k);
void   mp_rise2    (mp_t mp, int op, int  n);
void   mp_pow      (mp_t mp, int op, int  n, int  k);
void   mpf_ans     (mp_t mp, mpf_t a);
void   mpq_ans     (mp_t mp, mpq_t a);
void   mp_free     (mp_t mp);

mpf_t *mp_dkj      (int  n_max, double s, int  *j_stride, int  *d_dim);
mpf_t *mp_bkj      (int  n_max, double s, int  *j_stride, int  *b_dim);

mpf_t *mp_dkjn     (int  n_max, double s, int  *n_stride, int  *jn_stride, 
                    int  *kjn_dim);

void mp_wcc_all    (int  j1, int  j2, int  j, mpf_t *cfn);

#endif /* hex_gmp_h */
/*---------------------------------------------------------------------------*/
