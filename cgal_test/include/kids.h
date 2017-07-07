/*-----------------------------------------------------------------------------
**  File:       kids.h
**
**  Author:     Dave Ritchie, 27/08/11
**
**  Purpose:    Declarations for "keep it dead simple" array memory management.
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2011 D.W. Ritchie, INRIA.
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
**-----------------------------------------------------------------------------
** Here are some example typedefs from http://publications.gbdirect.co.uk/c_book/
**
**     typedef int aaa, bbb, ccc;
**     typedef int ar[15], arr[9][6];
**     typedef char c, *cp, carr[100];
**     
**     // now declare some objects 
**     
**     aaa     int1;   // int1 is an int
**     bbb     myint;  // so is myint
**     
**     ar      yyy;    // array of 15 ints
**     arr     xxx;    // 9*6 array of int
**     
**     c       ch;     // ch is a char 
**     cp      pnt;    // pnt is a pointer to char
**     carr    chry;   // chry is an array of 100 char
**
**---------------------------------------------------------------------------*/

#ifndef kids_h
#define kids_h

#include "hex_sys.h"

#ifdef __cplusplus
extern "C" {
#endif

// ----------------------------------------------------------------------------
// these unions allow mixed-type access to memory, but are not currently needed 
// ----------------------------------------------------------------------------

typedef union { double d;    int i; } UnionDI;
typedef union { float  f;    int i; } UnionFI;
typedef union { char   c[4]; int i; } UnionCI;

typedef UnionDI *                     KlistDI;
typedef UnionFI *                     KlistFI;

// ----------------------------------------------------------------------------
// these types allow un-decorated use of pointers, pointers-to-pointers, etc.
// they are named systematically, and some equivalent aliases have been 
// provided for the types - e.g. KarrayD/KmatrixD means 2-dimensional array
// or matrix of type "double".
// ----------------------------------------------------------------------------

typedef void            KwordV;
typedef void          * KlistV,     * KvectorV,   * K1V;
typedef void         ** KarrayV,   ** KmatrixV,  ** K2V;

typedef double          KwordD;
typedef double        * KlistD,     * KvectorD,   * K1D;
typedef double       ** KarrayD,   ** KmatrixD,  ** K2D;
typedef double      *** KcubeD,                 *** K3D;

typedef float           KwordF;
typedef float         * KlistF,     * KvectorF,   * K1F;
typedef float        ** KarrayF,   ** KmatrixF,  ** K2F;
typedef float       *** KcubeF,                 *** K3F;

typedef int             KwordI;
typedef int           * KlistI,     * KvectorI,   * K1I;
typedef int          ** KarrayI,   ** KmatrixI,  ** K2I;
typedef int         *** KcubeI,                 *** K3I;

typedef short int       KwordJ;
typedef short int     * KlistJ,     * KvectorJ,   * K1J;
typedef short int    ** KarrayJ,   ** KmatrixJ,  ** K2J;
typedef short int   *** KcubeJ,                 *** K3J;

typedef char            KwordC;          
typedef char          * KlistC,     * KwordS,     * K1C,    * Kstring;
typedef char         ** KarrayC,   ** KlistS,    ** K2C,   ** K1S;
typedef char        *** KcubeC,   *** KarrayS,  *** K3C,  *** K2S;

typedef const char      KwordCC;
typedef const char    * KlistCC,    * KwordSC,    * K1CC,   * KstringC;
typedef const char   ** KarrayCC,  ** KlistSC,   ** K2CC,  ** K1SC;
typedef const char  *** KcubeCC,  *** KarraySC, *** K3CC, *** K2SC;

// ----------------------------------------------------------------------------
// for internal use - this assumes that a pointer can be cast to a "long"
// ----------------------------------------------------------------------------

typedef long          * KlongP1,    * KP1;    
typedef long         ** KlongP2,   ** KP2; 
typedef long        *** KlongP3,  *** KP3; 

// ----------------------------------------------------------------------------
// un-decorated function declarations 
// ----------------------------------------------------------------------------

KlistV   kget1B   (long n_elements);      // these calls need to type-cast
KlistV   kmod1B   (KlistV addr, long n_elements);

KlistD   kget1D   (long n_elements);      // get a list/vector of double words
KlistF   kget1F   (long n_elements);      // floats
KlistI   kget1I   (long n_elements);      // ints
KlistJ   kget1J   (long n_elements);      // short ints

KlistV   kget1B_0 (long n_elements);
KlistD   kget1D_0 (long n_elements);      // same thing, but zero-fill
KlistF   kget1F_0 (long n_elements);
KlistI   kget1I_0 (long n_elements);
KlistJ   kget1J_0 (long n_elements);

KlistC   kget1C   (long n_elements);      // fixed size list of characters 
KlistC   kget1C_0 (long n_elements);

KlistS   kget1S   (long n_elements);      // list of pointers to strings

KarrayD  kget2D   (long n_rows, long n_cols);
KarrayF  kget2F   (long n_rows, long n_cols);
KarrayI  kget2I   (long n_rows, long n_cols);
KarrayJ  kget2J   (long n_rows, long n_cols);
KarrayC  kget2C   (long n_rows, long n_cols);

KarrayD  kget2D_0 (long n_rows, long n_cols);
KarrayF  kget2F_0 (long n_rows, long n_cols);
KarrayI  kget2I_0 (long n_rows, long n_cols);
KarrayJ  kget2J_0 (long n_rows, long n_cols);
KarrayC  kget2C_0 (long n_rows, long n_cols);

KcubeD   kget3D   (long n_planes, long n_rows, long n_cols);
KcubeD   kget3D_0 (long n_planes, long n_rows, long n_cols);
KcubeF   kget3F   (long n_planes, long n_rows, long n_cols);
KcubeF   kget3F_0 (long n_planes, long n_rows, long n_cols);
KcubeI   kget3I   (long n_planes, long n_rows, long n_cols);
KcubeI   kget3I_0 (long n_planes, long n_rows, long n_cols);
KcubeJ   kget3J   (long n_planes, long n_rows, long n_cols);
KcubeJ   kget3J_0 (long n_planes, long n_rows, long n_cols);
KcubeC   kget3C   (long n_planes, long n_rows, long n_cols);
KcubeC   kget3C_0 (long n_planes, long n_rows, long n_cols);

Kstring  kdupS    (KstringC string);       // duplicate a string 
Kstring  kcatS    (KstringC string, ...);  // args must be NULL-terminated

int      kcopy1S  (K1S to_list, int n_to, K1S from_list, int n_from);
int      kscan1S  (Kstring query, K1S list, int n_elements); // scan a list 

KlistS   kscanDir (KstringC path, KstringC pattern, int *n_elements); // scan a directory

void     kfree    (KlistV kthing);     // free memory for all of the above

// ----------------------------------------------------------------------------
// and finally...
// ----------------------------------------------------------------------------

#define  kwipe(x)      { if (x != NULL) { kfree(x); x = NULL; } }

#define  kgetD(x)      kget1D(x)    // these are just convenient aliases
#define  kgetF(x)      kget1F(x)
#define  kgetI(x)      kget1I(x)
#define  kgetJ(x)      kget1J(x)
#define  kgetB(x)      kget1B(x)

#define  kgetD_0(x)    kget1D_0(x)  // more aliases
#define  kgetF_0(x)    kget1F_0(x)
#define  kgetI_0(x)    kget1I_0(x)
#define  kgetJ_0(x)    kget1J_0(x)
#define  kgetB_0(x)    kget1B_0(x)

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif  /* kids_h */
// ----------------------------------------------------------------------------
