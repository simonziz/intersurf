/*-----------------------------------------------------------------------------
**  File:       hex_multipole.h
**
**  Author:     Dave Ritchie, 24/02/06
**
**  Purpose:    Define data structures for Buckingham Cartesian moments etc.
**
**  NB. The order of elements is as defined in the Vamp SDF files
**      (see mopac/src/sdwriter.f and /Parasurf/src/multipole.f90)
**
**-----------------------------------------------------------------------------
**
**  Copyright(c) 2006 Dave Ritchie, University of Aberdeen.
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

#ifndef hex_multipole_h
#define hex_multipole_h

#include "hex_geom.h"

typedef struct _multipole {

   double q;

   double x;
   double y;
   double z;

   double xx;
   double xy;
   double yy;
   double xz;
   double yz;
   double zz;

   double xxx;
   double xxy;
   double xxz;
   double yyx;
   double yyy;
   double yyz;
   double zzx;
   double zzy;
   double zzz;
   double xyz;

} HexMpole;

void hex_mpole_print(char *msg, HexMpole *mpole);
void hex_mpole2qlm(int l_max, HexMpole *mpole, double *qlm);
void hex_qlm2mpole(int l_max, double *qlm, HexMpole *mpole);
double hex_mpole_eval(int l_max, Point3D pt, HexMpole *mpole);
double hex_qlm_eval(int l_max, Point3D pt, double *qlm);

#endif  /* hex_multipole_h */

/*---------------------------------------------------------------------------*/

