/*-----------------------------------------------------------------------------
**  File:       hex_version.h
**
**  Author:     Dave Ritchie, 16/04/96
**
**  Purpose:    #include file for the various "hex" functions.
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 1996-2007 D.W. Ritchie, University of Aberdeen.
**  Copyright (C) 2010-2011 D.W. Ritchie, INRIA.
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
**---------------------------------------------------------------------------*/

#ifndef hex_version_h
#define hex_version_h

#define HEX_PROGRAM     0
#define PEPSI_PROGRAM   1
#define PFT_PROGRAM     2
#define TDB_PROGRAM     3
#define TDS_PROGRAM     4
#define KPAX_PROGRAM    5
#define SAM_PROGRAM     6

#define HEX_VERSION     "8.1.1"
#define KPAX_VERSION    "5.0.2"

#define TDB_VERSION     "2.11.0"
#define PFT_VERSION     "11.0.3"
#define TDS_VERSION     "2.1.0"
#define SAM_VERSION     "2.0.2"

#define PEPSI_VERSION   "1.0.1"

// ifdef __cplusplus
// extern "C" {
// endif

char *hex_version       ();  /* supply program name + version string */
char *hex_version_only  ();  /* supply version string only */
void  set_version       (int program);

//ifdef __cplusplus
//}
//endif

#endif /* hex_version_h */

/*---------------------------------------------------------------------------*/
