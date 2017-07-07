/*-----------------------------------------------------------------------------
**  File:       hex_dft.h
**
**  Author:     Dave Ritchie, 25/08/07
**
**  Purpose:    Header file for the "dft" functions.
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

#ifndef hex_dft_h
#define hex_dft_h

//ifdef __cplusplus
//extern "C" {
//endif

#define HEX_DFT_SINGLE        1
#define HEX_DFT_DOUBLE        2
#define HEX_DFT_REAL          3 
#define HEX_DFT_COMPLEX       4
#define HEX_DFT_FORWARD       5
#define HEX_DFT_BACKWARD      6
#define HEX_DFT_PLAN          7
#define HEX_DFT_NOPLAN        0

/*---------------------------------------------------------------------------*/

typedef struct _DftPlan {

   int    direction;
   int    precision;
   int    datatype;
   int    planning;
   int    n_dim;
   int    n_elements;
   long  *ldims;   /* needed for MKL  compatibility (32/64 bit) */
   int   *idims;   /* neeedd for FFTW compatibility (32/64 bit) */
   int   *kind;
   void  *data;
   void  *handle;

} DftPlan;

/*---------------------------------------------------------------------------*/

DftPlan *hex_dft_open (int domain, int precision, int datatype, int planning,
                       int n_dim, int *dims, void *buffer);
void hex_dft_execute  (DftPlan *dft);
void hex_dft_close    (DftPlan *dft);

void hex_kft_plan     (DftPlan *dft);
void hex_kft_execute  (DftPlan *dft);
void hex_kft_destroy  (DftPlan *dft);


/*---------------------------------------------------------------------------*/

//ifdef __cplusplus
//}
//endif

#endif  /* hex_dft_h */

/*---------------------------------------------------------------------------*/
