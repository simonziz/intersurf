/*-----------------------------------------------------------------------------
**  File:       hex_pdb.h
**
**  Author:     Dave Ritchie, 2012.
**  Update:     
**
**  Purpose:    Declare the types and structures for hex_pdb.cpp
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2012 Dave Ritchie, INRIA
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

#ifndef hex_pdb_h
#define hex_pdb_h

#include "hex_math.h"

#define  HEX_PDB_DIM     5   // only for the first 5 atom types below

#define  PDB_A           1   // alpha carbon
#define  PDB_N           2   // amino nitrogen
#define  PDB_C           3   // carbonyl carbon
#define  PDB_O           4   // carbonyl oxygen
#define  PDB_B           5   // beta carbon
#define  PDB_OXT         6  // main chain terminal OXT
#define  PDB_OT1         7  // main chain terminal OT1
#define  PDB_OT2         8  // main chain terminal OT2

/*---------------------------------------------------------------------------*/
// this is more-or-less an image of one ATOM line of a PDB file.

typedef struct _pdb_atom {

  char    record[7];
  char    atomName[5];
  char    eName[5];
  char    alt[2];
  char    residue[4];
  char    sequence[5];
  char    insert[2];
  char    chain[2];
  float   occ;
  float   tfac;
  Point3D pt;          // atom [X,Y,Z] coordinates, might be moved

  int     mid;         // internal model number, counting from zero
  int     cid;         // internal chain number, counting from zero
  int     rid;         // internal residue number, counting from zero

  byte    is_ca;       // internal flag = true if atom is CA
  byte    is_mc;       // internal main-chain flag (CA,N,O,C,CB)
  byte    is_bb;       // internal backbone flag (CA,N,O,C)
  byte    is_sc;       // internal side-chain flag

  byte    is_acc;      // becomes 1 if solvent accessible (from hex_accessible)

  double  radius;      // = atom radius from Hex
  char    seqid[6];    // = sequence + alt

} PdbAtom;

/*---------------------------------------------------------------------------*/
// similarly, define just enough information about one PDB residue

typedef struct _pdb_residue {

   int    mask;        // 0 = no special actions, 1 = do not write on output

   int    mid;        // internal model number (0 => no "model")
   int    cid;        // internal chain number, counting from zero

   int    atom1;      // index no. of first and last atom
   int    atom2;
   int    n_sc;       // number of side chain atoms 

// Point3D pt_scm;      // side chain centre of mass
// Point3D pt_scm_orig; // side chain centre of mass

   char   chain[2];    // null-terminated chain letter
   char   code[2];     // null-terminated amino-acid code letter
   char   name[4];     // null-terminated amino-acid name
   char   sse[2];      // null-terminated secondary structure code latter
   char   label[6];    // alpha-numeric residue sequence label

} PdbResidue;

/*---------------------------------------------------------------------------*/
// Similarly, define a structure containing the backbone atom coordinates

typedef struct _pdb_peptide {

   Point3D  *ca;   // alpha carbon
   Point3D  *n;    // amino nitorgen
   Point3D  *c;    // carbonyl carbon
   Point3D  *o;    // carbonyl oxygen
   Point3D  *cb;   // beta carbon
// Point3D  *sm;   // side chain centre of mass
   double   *wt;   // vector of weights
   char     *aa_code;  // amino acid code
   char     *aa_flag;  // amino acid flag (dead or SSE type)

} PdbPeptide;

/*---------------------------------------------------------------------------*/
// load a PDB file into memory (but don't assume n_atoms=n_all > n_bb >n_ca)

typedef struct _pdb_image {

   int      n_models;     // the apparent number of "models" in the PDB file
   int      n_chains;     // the apparent number of chains in the PDB file
   int      n_atoms;      // the total number of all non-hetero atoms in the image
   int      n_hetero;     // the total number of all hetero atoms in the image
   int      n_ca;         // the number of CA atoms in the PDB file
   int      n_bb;         // the number of backbone atoms in the PDB file
   int      n_mc;         // the number of main chain (backbone + CB)
   int      n_sc;         // the number of sidechain atoms in the PDB file

   int      n_cryst;      // the number of SMTRY matrices in the PDB file
   int      n_biomt;      // the number of BIOMT matrices in the PDB file

   int      c_delta;      // chain letter delta for symmetry renumbering

   int      is_fake;      // 1=>artificial structure, 0=>real PDB structure

   char    *name;         // just like TdbThing
   char    *nickname;     // just like TdbThing

   Matrix3D mat_import;   // transformations applied when loading the PDB

   Matrix3D *mat_biomt;   // BIOMT symmetry matrices
   Matrix3D *mat_cryst;   // SMTRY symmetry matrices

   PdbResidue *residue;   // list of amino acid residues (length = n_ca).

   PdbAtom *atom;         // list of "ATOM" records
   PdbAtom *hetero;       // list of "HETATM" records

   Point3D *atom_backup;   // backup copy of atom coordinates, length=n_atoms
   Point3D *hetero_backup; // backup copy of hetero atom coords, length=nhetero_

   float   energy_isc;    // RosettaDock side chain interface energy
   float   energy_hex;    // Hex docking energy (not used)

} PdbImage;

/*---------------------------------------------------------------------------*/

PdbImage  *hex_readPdb       (char *filename, char *nickname);
void       hex_freePdb       (PdbImage *image);
PdbImage  *hex_copyPdb       (PdbImage *image);
void       hex_maskPdb       (PdbImage *image, int mask);
void       hex_movePdb       (PdbImage *image, Matrix3D tf_pdb);
void       hex_resetPdb      (PdbImage *image);
void       hex_freezePdb     (PdbImage *image);
void       hex_movePdbResidue(PdbImage *image, int rid, Matrix3D tf_pdb);
void       hex_labelPdbResidue(PdbImage *image, int rid, int delta);
Point3D    hex_ptPdbResidue  (PdbImage *image, int rid);
void       hex_writePdb      (char *filename, PdbImage *image, int chain_id, 
                              char chain_label, int atom_mode, 
                              int with_sse, int with_ter, Matrix3D t_export);
double     hex_importPdb     (PdbImage *image, int use_random);
void       hex_exportPdb     (char *from_file, char *to_file, int atom_mode,
                              Matrix3D t_export);
int        hex_chainsPdb     (PdbImage *pdb, char *chains, int *the_chains);

void       hex_splitPdb      (PdbImage *image, int n_chains, int *the_chains,
                              char **the_names, PdbImage **the_images);
PdbImage  *hex_extractPdb    (PdbImage *pdb, int n_chains, int *the_chains);

HexFile   *hex_openPdbFile   (char *filename);
int        hex_writePdbFile  (HexFile *file, PdbImage *image, int chain_id, 
                              char chain_label, int atom_mode, int with_sse,  
                              int with_ter, int atom_start, Matrix3D t_export);
//int      hex_writePdbResidue(HexFile *file, PdbImage *image, int rid, 
//                            char chain_label, int atom_mode, int with_ter, 
//                            int atom_start, Matrix3D t_export);
void       hex_closePdbFile  (HexFile *file);

int        hex_writeBinaryPdb(char *filename, PdbImage *image);
PdbImage  *hex_readBinaryPdb (char *filename, char *nickname);

int        hex_getPdbAlphaNos(PdbImage *image, int chain_id, int atom_mode, 
                              int atom_num, int *a_num);

PdbPeptide *hex_allocPdbPep  (int n_peptides);
PdbPeptide *hex_makePdbPep   (PdbImage *image);
PdbPeptide *hex_dupPdbPep    (int n_peptides, PdbPeptide*src);
void        hex_copyPdbPep   (PdbPeptide *src, int start, PdbPeptide *dst, int np);
void        hex_loadPdbPep   (PdbImage *image, PdbPeptide *pep);
void        hex_freePdbPep   (PdbPeptide *pep);

Point3D   *hex_ptsPdb        (PdbImage *image, int atom_mode);
float     *hex_radPdb        (PdbImage *image, int atom_mode);

void       hex_setPdbDelta   (PdbImage *image, int delta);
void       hex_setPdbFake    (PdbImage *image, int fake_flag);

Matrix3D   hex_pdbBiomt      (PdbImage *image, int biomt);
Matrix3D   hex_pdbCryst      (PdbImage *image, int cryst);

void       hex_setPdbChainLabel(PdbImage *image, char new_chain);

void       hex_setPdbChainType(int chain_type);
void       hex_setPdbHydrogens(int hydrogens);
void       hex_setPdbHetRes   (int hetres);
void       hex_setPdbMissing  (int discard);
void       hex_setPdbWarnings (int warnings);

void       hex_setPdbOneModel (int only_one);  // read up to next ENDMDL
void       hex_setPdbOneChain (int only_one);  // read up to next TER

int        hex_isAlphaCarbon (char *atom_name);
int        hex_isMainChain   (char *atom_name);
int        hex_isBackbone    (char *atom_name);
int        hex_isHydrogen    (char *atom_name);
int        hex_isSideChain   (char *atom_name);
char       hex_aaCode        (char *resdiue_name);
char      *hex_aaName        (char residue_code);

int        hex_isAA          (char *resdiue_name);
int        hex_isDNA         (char *resdiue_name);

char       hex_dnaCode       (char *resdiue_name);
int        hex_isBase        (char *atomName);
int        hex_isSugar       (char *atomName);

// not sure if these really belong here ?

typedef struct _buried_blobs {
           Point3D   centre;
           double    radius;
           int       n_ids;
           int      *atom_ids;
} BuriedBlob;

void       hex_getImageAtomRadii(PdbImage *image);

void       hex_accessible      (PdbImage *image, double probe_radius);
int        hex_buried_blobs    (PdbImage *image, BuriedBlob **blobs);
void       hex_free_blobs      (int n_blobs, BuriedBlob *blobs);
void       hex_accessible_write(PdbImage *image, double probe_radius);
void       hex_accessible_test (PdbImage *image, double probe_radius);

#endif  /* hex_pdb_h */
/*---------------------------------------------------------------------------*/
