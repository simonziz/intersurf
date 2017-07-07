/*-----------------------------------------------------------------------------
**  File:       hex_dx.h
**
**  Author:     Dave Ritchie, 25/08/07
**
**  Purpose:    Header file for the "dx" functions.
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

#ifndef hex_dx_h
#define hex_dx_h

#include "hex_transform.h"

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------------*/

typedef struct _dxg {  /* DX format grid in ZYX data order */

   Point3D origin;  /* coordinates of corner of cell [0,0,0] */
   float   dg;      /* grid cell size (delta) */
   int     nx;      /* dimensions */
   int     ny;
   int     nz;
   float   *values; /* grid values (at centre of grid cells) */

} DxGrid;


/*---------------------------------------------------------------------------*/

DxGrid *dx_read         (char *filename);
void    dx_write        (char *filename, DxGrid *dxg);
void    dx_free         (DxGrid *dxg);
DxGrid *dx_create_xyz   (int nx, int ny, int nz, Point3D origin, float *values);
DxGrid *dx_create_zyx   (int nx, int ny, int nz, Point3D origin, float *values);

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif  /* hex_dx_h */

/*---------------------------------------------------------------------------*/
