/*-----------------------------------------------------------------------------
**  File:       hex_mol.h
**
**  Author:     Dave Ritchie, 30/06/06
**
**  Purpose:    Header file for file utilities ...
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2006 D.W. Ritchie, University of Aberdeen.
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

#ifndef hex_mol_h
#define hex_mol_h

#include "hex_sys.h"
#include "hex_math.h"
#include "hex_util.h"

//ifdef __cplusplus
// extern "C" {
//endif

/*---------------------------------------------------------------------------*/

#define MAX_LINE_LENGTH   133            /* 132-char card image + 1 NULL */

#define DEF_RADIUS        1.5
#define DEF_CHARGE        0.0
#define DEF_OCCUPANCY     1.0
#define DEF_TEMPERATURE  99.99

#define MOL_SYM         240              /* no. symmetry transformations */
#define MOL_STRUCTURES    3              /* no. simultaneous structures */
#define MOL_PROBES        3              /* max no. GRID probe types */

/* these parameters are array increments, not limits :-) */

#define MOL_MODELS     2000              /* no. model structures */
#define MOL_CHAINS     2000              /* no. chains in PDB file */
#define MOL_RESIDUES   2000              /* no. residues in PDB file */
#define MOL_ATOMS     10000              /* no. atoms in PDB file */
#define MOL_BONDS     10000              /* no. bonds between atoms */
#define MOL_PAIRS      1000              /* no. CONECT pairs PDB file */
#define MOL_ENERGIES   5000              /* no. ENERGY records from GRID */

/*---------------------------------------------------------------------------*/

typedef struct _HexAtom {

   int      n_atoms;
   int      n_bonds;
   char    *iname;
   char    *ename;
   float   *radius;
   float   *charge;
   Point3D *pt;
   float   *occ;
   float   *tfac;
   char    *type;
   char    *access;
   char    *hetatm;
   char    *altatm;
   sint    *defcol;
   sint    *colour;
   sint    *opls;
   sint    *ace;
   byte    *pr;
   byte    *lp;
   int     *bond;
   byte    *bflag;     /* 0=auto, 1=calculated, 2=from PDB */
   int     *index;     /* atom neighbour lookup table */
   int     *table;

} HexAtom;

typedef struct _HexResidue {

   int      n_residues;
   char    *name;
   char    *sstype;
   char    *id;
   int     *atoms;

} HexResidue;

typedef struct _HexChain {

   int      n_chains;
   char    *id;
   int     *residues;
   int     *bonds;       
   int     *atoms;       
   sint    *colour;

} HexChain;

typedef struct _HexModel {

   int      n_models;
   int     *atoms;
   Point3D *cog;      /* centre of gravity */
   Point3D *coh;      /* centre of harmonics */
   HexShd  **shd;     /* spherical harmonic data block */
   sint    *id;
   char    *label;

} HexModel;

typedef struct _HexPot {

   int      n_potential;
   sint    *code;
   Point3D *pt;
   float   *value;
   int      probe[MOL_PROBES];
   int      colour[MOL_PROBES];

} HexPot;

typedef struct _HexSym {

   int             n_crystal;
   int             n_biol;
   Matrix3D *crystal;
   Matrix3D *biol;

} HexSym;
   

typedef struct _HexMol {

   int         open;
   int         is3DBlast;
   int         ss_done;
   int         ss_method;
   int         mol_type;
   int         col_type;
   int         sym_type;
   int         sym_no;
   char       *filename;
   char       *label;

   HexAtom    *matm;
   HexResidue *mres;
   HexChain   *mchn;
   HexModel   *modl;
   HexPot     *mpot;
   HexSym     *msym;

   Point3D     com;             /* CoM - centre of mass */
   double      radius;          /* bounding radius wrt CoM */

   Matrix3D t_user;
   Matrix3D t_base;

} HexMol;

/*---------------------------------------------------------------------------*/

int   mol_open          (int mol_index, char *filename);
void  mol_close         (int);
void  mol_select        (int);
void  mol_sse_method    (int method);  
char *mol_filename      ();
int   mol_type          ();
char *mol_code          ();
void  mol_chain_sequence(int, char []);
int   mol_models        (); 
int   mol_model_id      (int model_no);
char *mol_model_label   (int model_no); 
int   mol_lookup_model  (char *model_label);
void  mol_model_atoms   (int model_no, int  *atom1, int  *atom2);
void  mol_structure     (int  *, int  *, int  *, int  *);
void  mol_chain         (int,char[],int  *, int  *, int  *);
char  mol_chain_id      (int);
int   mol_chain_colour  (int);
int   mol_residue_chain (int rid);
void  mol_residue_name  (int rid, char *name, char *seq_id);
int   mol_residue_sstype(int rid);
void  mol_residue_atoms (int rid, int *a1, int *a2);
int   mol_lookup_pattern(char *pattern, int show);

int   mol_resnum_residue(int chain, int resnum);
void  mol_resnum_info  (int chain, int resnum, char *name, char *seq_id,
                         int *a1, int *a2);
int   mol_lookup_chain  (char chain_id);
int   mol_lookup_resnum(int  chain, char *seq_id);

void  mol_atom_set_access(int, int);
void  mol_atom_set_colour(int, int);
int   mol_atom_residue  (int atom);
int   mol_atom_resnum   (int atom);
char *mol_atom_ename    (int);
char *mol_atom_iname    (int);
int   mol_atom_model    (int  atom);
char *mol_atom_name     (int);
Point3D mol_atom_pt     (int);
int   mol_pt_atom       (int heavy, Point3D pt);
void  mol_pt_bond       (int heavy, Point3D pt, int  *chain, int  *bond);
Point3D mol_atom_basept (int);
float mol_atom_charge   (int);
float mol_atom_radius   (int);
char  mol_atom_class    (int);
char  mol_atom_hetatm   (int);
char  mol_atom_altatm   (int);
int   mol_atom_opls     (int);
int   mol_atom_ace      (int);
int   mol_atom_pr       (int);
int   mol_atom_lp       (int);
int   mol_atom_access   (int);
int   mol_atom_colour   (int);
void  mol_colour_mode   (int);
float mol_atom_occupancy(int);
float mol_atom_temperature(int);
int   mol_atom_chain    (int);
int   mol_atom_id       (int, char [], char [], char []);
int   mol_atom_atoms    (int, int  **);
void  mol_atom_set_basept(int, Point3D coords);
void  mol_bond_atoms    (int, int, int  *, int  *);
int   mol_bond_flag     (int, int);
void  mol_header        (char [],char []);

Point3D mol_centroid    ();
Point3D mol_model_cog   (int model);
Point3D mol_model_coh   (int model);
HexShd *mol_model_shd   (int model);
int     mol_3dblast     ();
double  mol_radius      ();
double  mol_radius_sym  ();
void  mol_transform     (Matrix3D);
Matrix3D mol_get_transform ();
Matrix3D mol_get_symmetry  (int sym);
void  mol_symmetry_type (int type);
void  mol_symmetry_num  (int sym);
int   mol_nsymmetry     ();
int   mol_match_atom    (char [], int, int);
char  mol_short_code    (char []);
char *mol_long_code     (char);
int   mol_is_polar      (int  iatom);
int   mol_is_apolar     (int  iatom);
int   mol_is_h          (int  iatom);
int   mol_is_dna        (char *residue_name);
int   mol_is_aa         (char *residue_name);
int   mol_is_aaa        (char *residue_name);
int   mol_energies      ();
int   mol_energy_colour (int ienergy);
int   mol_energy_code   (int ienergy);
float mol_energy_value  (int ienergy);
Point3D mol_energy_pt   (int ienergy);

void  mol_stride_pipe   ();
void  mol_dssp_pipe     ();
void  mol_kpax_pipe     ();

char  mol_aa_code       (char *residue_name);
char *mol_aa_res        (char code);


int   hex_load_tdb      (char *filename, int n_max, double **anlm_out);

int   hex_load_pdb      (char *filename, HexMol *mo);
int   hex_load_sdf      (char *filename, HexMol *mo);



void mol_mem_alloc      (HexMol *mo, int  na, int  nb, int  nr, int  nc,
                                     int  nm, int  ne);
void mol_mem_realloc    (HexMol *mo);
void mol_mem_free       (HexMol *mo);


/*---------------------------------------------------------------------------*/

//ifdef __cplusplus
//}
//endif

#endif  /* hex_mol_h */

/*---------------------------------------------------------------------------*/
