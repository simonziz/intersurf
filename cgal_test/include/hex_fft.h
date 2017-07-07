/*-----------------------------------------------------------------------------
**  File:       hex_fft.h
**
**  Author:     Dave Ritchie, 19/09/09
**
**  Purpose:    Header file for the "fft" functions.
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2007 D.W. Ritchie, University of Aberdeen.
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

#ifndef hex_fft_h
#define hex_fft_h

#include "hex_math.h"

/*---------------------------------------------------------------------------*/

//ifdef __cplusplus
//extern "C" {
//endif

/*---------------------------------------------------------------------------*/

typedef struct _HexFftPlan {

   int   use_fht;    /* 1 = FFT by FHT, 0 = FFT by FFT */
   int   np;         /* order of current transform */
   int   nn;         /* length of current transform (2^np) */

   int    *order;    /* bit-reversed data order table */
   double *cosx;     /* lookup table for cos(m*alpha) */
   double *sinm;     /* lookup table for sin(m*alpha) */
   double *sinp;     /* lookup table for sin(m*alpha) */
   double *wrkr;     /* scratch array for hex_fft_fs() */
   double *wrki;     /* scratch array for hex_fft_fs() */

} HexFftPlan;

/*---------------------------------------------------------------------------*/

HexFftPlan *hex_fft_open       (int n_dim);
void     hex_fft_mode       (HexFftPlan *fp, int mode);
void     hex_fft_execute_fs (HexFftPlan *fp, int n_max, double *am, double *bm,
                             double *result);
void     hex_fft_execute_1d (HexFftPlan *fp, int direction, double *data);
void     hex_fft_close      (HexFftPlan *fp);

void     hex_fft_order      (int n_pow, int *sort);

/*---------------------------------------------------------------------------*/

//ifdef __cplusplus
//}
//endif

#endif  /* hex_fft_h */

/*---------------------------------------------------------------------------*/
