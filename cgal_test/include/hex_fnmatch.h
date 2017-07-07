/*-----------------------------------------------------------------------------
**  File:       hex_fnmatch.h
**
**  Author:     Dave Ritchie, 17/11/2013
**
**              Extracted from binutils "fnmatch.h" 
**
**              The original copyright is included below.
*---------------------------------------------------------------------------*/

/* Copyright 1991, 1992, 1993, 1996 Free Software Foundation, Inc.

NOTE: The canonical source of this file is maintained with the GNU C Library.
Bugs can be reported to bug-glibc@prep.ai.mit.edu.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, 51 Franklin Street - Fifth Floor,
Boston, MA 02110-1301, USA.  
*/

/* -------------------------------------------------------------------------*/

#ifndef	_HEX_FNMATCH_H
#define	_HEX_FNMATCH_H	1

#include <dirent.h>
#include <string.h>
#include <ctype.h>

#ifdef	__cplusplus
extern "C" {
#endif

/* Bits set in the FLAGS argument to `fnmatch'.  */

#define	FNM_PATHNAME	(1 << 0) /* No wildcard can ever match `/'.  */
#define	FNM_NOESCAPE	(1 << 1) /* Backslashes don't quote special chars.  */
#define	FNM_PERIOD	(1 << 2) /* Leading `.' is matched only explicitly.  */
#define	FNM_LEADING_DIR	(1 << 3) /* Ignore `/...' after a match.  */
#define	FNM_CASEFOLD	(1 << 4) /* Compare without regard to case.  */

#define	FNM_FILE_NAME	FNM_PATHNAME /* Preferred GNU name.  */

/* Value returned by "hex_fnmatch" if STRING does not match PATTERN.  */
#define	FNM_NOMATCH	1

int hex_fnmatch(const char *pattern, const char *string, int flags);

int hex_alphasort(const struct dirent **a, const struct dirent **b);

#ifdef	__cplusplus
}
#endif

#endif /* _HEX_FNMATCH_H */
