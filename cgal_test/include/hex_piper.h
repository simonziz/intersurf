/*-----------------------------------------------------------------------------
**  File:       hex_piper.h
**
**  Author:     Dave Ritchie, 23/07/07
**
**  Purpose:    #include file for Dima Kozakov's Hex-PIPER interface code
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

#ifndef hex_piper_h
#define hex_piper_h

#include "hex_dx.h"

//ifdef __cplusplus
//extern "C" {
//endif


#ifdef hex_piper

#include "mol.0.0.4.h"

void hex_piper_calc(char *dxname, int n_max, double *anlm, 
                    int nx, int ny, int nz, double dg, Point3D origin);

void piper_freegrids(struct grid **grids, int num);

void piper_set_repulsion_grid (struct grid* g, struct atomgrp* p, struct prms* prms);

struct grid_size *piper_set_gridsize(struct atomgrp* proteinA, float cell_size);

void hex_pe_pipergrids(Molecule *molecule, int mdl, int  n_max,
                  int total, struct grid**  regs);

void hex_pe_piper(Molecule *molecule, int mdl, int  n_max,
                  int which, char *dxname);

void hex_piper_calc(char *dxname, int n_max, double *anlm,
                    int nx, int ny, int nz, double dg, Point3D origin);

struct grid **piper_make_grids(Molecule *mol1, int mdl1,
                               struct prms* prms, int nvalues );

struct atomgrp *hex_mol2piper(Molecule *molecule, int mdl, struct prms *prms);


#endif

//ifdef __cplusplus
//}
//endif

#endif /* hex_piper_h */

/*---------------------------------------------------------------------------*/
