/*-----------------------------------------------------------------------------
**  File:       hex_main.h
**
**  Author:     Dave Ritchie, 28/05/03
**
**  Purpose:    #include file for the "main" functions.
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2003 D.W. Ritchie, University of Aberdeen.
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

#ifndef hex_main_h
#define hex_main_h

/*---------------------------------------------------------------------------*/

//ifdef __cplusplus
//extern "C" {
//endif

int hex_main   (int argc, char **argv, char **envp);
int hex_parser (int argc, char **argv, char **envp);
int hex_test   (int argc, char **argv, char **envp);
int hex_plot   (int argc, char **argv, char **envp);

void hex_first_actions ();
void hex_last_actions  ();

/*---------------------------------------------------------------------------*/

//ifdef __cplusplus
//}
//endif

#endif  /* hex_main_h */

/*---------------------------------------------------------------------------*/
