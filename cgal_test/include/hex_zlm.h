/*-----------------------------------------------------------------------------
**  File:       hex_zlm.h
**
**  Author:     Dave Ritchie, 20/07/07
**
**  Purpose:    Header file for complex spherical harmonics. 
**
**  This assumes "gcc -std=c99" or equivalent is available to provide built-in
**  complex data types. Seems not to be the case on SunOS "cc". Seems not to
**  be compatible with the g++ complex type, so try to keep separate for now.
**
**-----------------------------------------------------------------------------
**  
**  Copyright (c) 2007 Dave Ritchie, University of Aberdeen.
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

#ifndef hex_zlm_h
#define hex_zlm_h

#include <complex>
using namespace std;

#include "hex_math.h"

// typedef complex<float>       zfloat;
// typedef complex<double>      zdouble;
// typedef complex<long double> zzdouble;

#define I  zdouble(0,1)     // deprecated

zdouble zneg1_m     (int m);
zdouble zi_m        (int m);

zdouble zeim        (int m, double phi);
zdouble zlm         (int l, int m, double theta, double phi);
void    zlm_all_bar (int l_max, double theta, double phi, 
                             zdouble *zfn);
void    zjoin       (int n, double  *a,  double *b, zdouble *z);
void    zsplit      (int n, zdouble *z,  double *a, double *b);
void    zreal       (int n, zdouble *z,  double *a);
void    zimag       (int n, zdouble *z,  double *a);
void    zswap       (int n, zdouble *z,  zdouble *s);
void    zconj       (int n, zdouble *z,  zdouble *c);
void    zconjf      (int n, zfloat  *z,  zfloat  *c);
void    zd2zf       (int n, zdouble *zd, zfloat *zf);
void    zf2zd       (int n, zfloat  *zf, zdouble *zd);

void    zzreal      (int n, zzdouble *z, double *a);
void    zzimag      (int n, zzdouble *z, double *a);

void    ylm_zlm     (int l_max, zdouble *zlm, double         *ylm);
void    clm_zlm     (int l_max, zdouble *zlm, zdouble *clm);

void    zlm_ylm     (int l_max, double         *ylm, zdouble *zlm);
void    zlm_clm     (int l_max, zdouble *clm, zdouble *zlm);

void    alm_qlm     (int l_max, zdouble *qlm, double *alm);
void    clm_qlm     (int l_max, zdouble *qlm, zdouble *clm);
void    qlm_alm     (int l_max, double *alm, zdouble *qlm);
void    qlm_clm     (int l_max, zdouble *clm, zdouble *qlm);
void    qlm_ablm    (int l_max, double *alm, double *blm, zdouble *qlm);
void    anlm_qnlm   (int n_max, zdouble *qnlm, double *anlm);
void    cnlm_qnlm   (int n_max, zdouble *qnlm, zdouble *cnlm);
void    qnlm_anlm   (int n_max, double *anlm, zdouble *qnlm);
void    qnlm_cnlm   (int n_max, zdouble *cnlm, zdouble *qnlm);
void    qnlm_conj   (int n_max, zdouble *qnlm, zdouble *cnlm);
void    qnlm_abnlm  (int n_max, double *anlm, double *bnlm, zdouble *qnlm);

zdouble zanlm_dot   (int n_max, zdouble *anlm, zdouble *bnlm);
zdouble zanlm_stupid(int  n_vec, int  n_max, zdouble *anlm, zdouble *bnlm);

void zanlm_rotate   (int n_max, int n_vec, double alpha, double beta, double gamma,
                     zdouble *zanlm, zdouble *zbnlm);

void zanlm_rot_phi  (int  n_max, double phi, 
                     zdouble *zanlm, zdouble *zbnlm);

void znlm_phi       (int n_max, double phi, zdouble *znlm_phi);
void znlm_mult      (int n_max, zdouble *zanlm, zdouble *zbnlm, zdouble *zcnlm);


#endif  /* hex_zlm_h */

/*---------------------------------------------------------------------------*/
