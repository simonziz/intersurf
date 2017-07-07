/*-----------------------------------------------------------------------------
**  File:       hex_roc.h
**
**  Author:     Dave Ritchie, 2012.
**  Update:     
**
**  Purpose:    Declare the types and structures for hex_roc.cpp
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2012 Dave Ritchie, INRIA
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

#ifndef hex_roc_h
#define hex_roc_h

#include "hex_math.h"

/*---------------------------------------------------------------------------*/
// data structure to hold the TPR/FPR data for a plot

typedef struct _hex_roc {

   int      n_roc;        // the number of (TPR,FPR) points in the ROC plot

   int      n_samp;       // the number of data samples (length of cls)
   int      n_pos;        // the number of positive samples
   int      n_neg;        // the number of negative samples (or 1)

   float   *tpr;          // true positive rate
   float   *fpr;          // false positive rate
   byte    *cls;          // list of true or false sample classes

} HexRoc;

/*---------------------------------------------------------------------------*/

HexRoc *hex_makeROC    (int n_samples, int *matches);
float   hex_aucROC     (HexRoc *roc);
void    hex_sampleROC  (HexRoc *roc, int n_samples, float *x, float *y, float *a);
void    hex_freeROC    (HexRoc *roc);

#endif  /* hex_roc_h */
/*---------------------------------------------------------------------------*/

