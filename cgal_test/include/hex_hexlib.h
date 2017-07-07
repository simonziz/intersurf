/*-----------------------------------------------------------------------------
**  File:       hex_hexlib.h
**
**  Author:     Dave Ritchie, 25/08/07
**
**  Purpose:    Header file for the "hexlib" functions.
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

#ifndef hex_hexlib_h
#define hex_hexlib_h

#include "hex_units.h"
#include "hex_types.h"
#include "hex_sys.h"
#include "hex_version.h"
#include "hex_transform.h"
#include "hex_geom.h"
#include "hex_dft.h"
#include "hex_fft.h"
#include "hex_dx.h"
#if defined(hex_pft)
#include "hex_harmonics.h"
#else
#include "hex_math.h"
#include "hex_zlm.h"
#endif

#endif  /* hex_hexlib_h */

/*---------------------------------------------------------------------------*/
