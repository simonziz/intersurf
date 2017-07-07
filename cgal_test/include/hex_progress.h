/*-----------------------------------------------------------------------------
**  File:       hex_progress.h
**
**  Author:     Dave Ritchie, 05/04/2012
**
**  Purpose:    header for hex_progress.cpp
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2012 D.W. Ritchie, INRIA.
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
---------------------------------------------------------------------------*/

#ifndef hex_progress_h
#define hex_progress_h

#include "hex_sys.h"

#ifdef __cplusplus
  extern "C" {
#endif

/*---------------------------------------------------------------------------*/

void   hex_progress_callback  (hex_progress_cb);

void   hex_progress_filename  (char *filename);
void   hex_progress_weights   (double *wts);
void   hex_progress_start     ();
void   hex_progress_end       ();
int    hex_file_progress      (int phase, double fraction);
int    hex_progress1          (int phase, double fraction);
int    hex_progress2          (int phase, int part, int all);

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
  }
#endif

#endif /* hex_progress_h */
/*---------------------------------------------------------------------------*/
