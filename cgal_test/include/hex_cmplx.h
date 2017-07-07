/*-----------------------------------------------------------------------------
**  File:       hex_cmplx.h
**
**  Author:     Dave Ritchie, 20/07/07
**
**  Purpose:    Header file for complex arithmetic
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

#ifndef hex_cmplx_h
#define hex_cmplx_h

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------------*/

/* old home-made stuff */

typedef struct _my_complex {
   double r;
   double i;
} Complex;


Complex cmake        (double  xr, double xi);
void    cunpack      (Complex z, double  *zr, double *zi);
Complex cneg         (Complex x);
Complex cadd         (Complex x, Complex y);
Complex csub         (Complex x, Complex y);
Complex cscale       (double  s, Complex x);
Complex cmul         (Complex x, Complex y);
Complex cdiv         (Complex x, Complex y);

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif  /* hex_cmplx_h */

/*---------------------------------------------------------------------------*/
