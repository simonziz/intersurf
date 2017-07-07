/*-----------------------------------------------------------------------------
**  File:       pdbs.h
**
**  Author:     Dave Ritchie, 15/08/05
**
**  Purpose:    Main Pepsi include file...
**
**-----------------------------------------------------------------------------
**
**  Copyright(c) 2005 Dave Ritchie, University of Aberdeen.
**
**  This software (or modified copies thereof) is protected by copyright and
**  may not be redistributed in any way without the express permission of the
**  author. Any questions about this copyright notice should be addressed to 
**  dritchie@csd.abdn.ac.uk.
**
**---------------------------------------------------------------------------*/

#ifndef pdbs_h
#define pdbs_h

/*
include <stdio.h>
include <stdlib.h>
include <malloc.h>
include <sys/types.h>
include <sys/stat.h>
include <unistd.h>
include <errno.h>
*/

#include "hex_pdb.h"
#include "hex_sys.h"
#include "hex_util.h"
#include "hex_transform.h"
#include "hex_version.h"
#include "kids.h"

#if defined(__GCC__) || defined(hex_icc)
#define DDI_INLINE static __inline__
#else
#define DDI_INLINE static
#endif

#define DDI_MAX_FILENAME    256
#define DDI_LINE_LENGTH     132

#define DDI_FIT             0


/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------------------------*/

void   pdbs_process     (char *input_file, char *write_file, char *sym_label,
                         int n_sym, int select_chain, int do_fasta, int do_pdb, 
                         int do_sum);

char  *pdbs_findfile    (char *filename, char *def_ext);

#ifdef __cplusplus
}
#endif

#endif  /* pdbs_h */
/*---------------------------------------------------------------------------*/

