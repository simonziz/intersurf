/*-----------------------------------------------------------------------------
**  File:       hex_pdb.cpp
**
**  Author:     Dave Ritchie, 2012 
**
**  Update:     Dave Ritchie, 2013, now loads SYMTRY and BIOMT.
**               Dave Ritchie,2014, removed side chain COM; added backup 
**                                  coords for reset after flexible fitting.
**
**  Purpose:    To read and write PDB files. This should be easy enough, but
**              PDB files can contain all kinds of problems including partial
**              residues, missing atoms, duplicate atoms, etc.
**
**              NB. currently, it is up to the caller to handle multiple 
**              MODELs by checking the value of image->residue[r].model, etc.
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2012-2014 Dave Ritchie, INRIA
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

#include "kids.h"
#include "hex_util.h"
#include "hex_version.h"
#include "hex_pdb.h"

#define  MAX_SYMMETRY       242

static int use_hydrogens   = 0;
static int use_hetres      = 0;  
static int use_warnings    = 0;  
static int strip_alt       = 0;
static int save_space      = 1;
static int discard_missing = 0;
static int one_model       = 1;
static int one_chain       = 0;
static int select_chain    = 1;
static int chunk_size      = 5000;

extern int hex_debug;

/* ----------------------------------------------------------------------------- */

void hex_setPdbChainType(int chain_type)
{
   if (chain_type == 0) select_chain = 0;  // protein + DNA/RNA
   if (chain_type == 1) select_chain = 1;  // protein only
   if (chain_type == 2) select_chain = 2;  // DNA/RNA only
}

/* ----------------------------------------------------------------------------- */

void hex_setPdbHydrogens(int hydrogens)
{
   use_hydrogens = hydrogens;
}

/* ----------------------------------------------------------------------------- */
// read only up to the next ENDMDL keyword and then close or return

void hex_setPdbOneModel(int only_one)
{
   one_model = only_one;
}

/* ----------------------------------------------------------------------------- */
// read only up to the next TER keyword and then close or return

void hex_setPdbOneChain(int only_one)
{
   one_chain = only_one;
}

/* ----------------------------------------------------------------------------- */

void hex_setPdbHetRes(int hetres)
{
   use_hetres = hetres;
}

/* ----------------------------------------------------------------------------- */

void hex_setPdbWarnings(int warnings)
{
   use_warnings = warnings;
}

/* ----------------------------------------------------------------------------- */

void hex_setPdbMissing(int discard)
{
   discard_missing = discard;
}

/* ----------------------------------------------------------------------------- */

void hex_setPdbDelta(PdbImage *image, int delta)
{
   image->c_delta = delta;
}

/* ----------------------------------------------------------------------------- */

void hex_setPdbFake(PdbImage *image, int fake)
{
   image->is_fake = fake;
}

/* ----------------------------------------------------------------------------- */

Matrix3D hex_pdbBiomt(PdbImage *image, int biomt)
{
   if (biomt >= 1 && biomt <= image->n_biomt) return(image->mat_biomt[biomt-1]);
   else                                       return(tf_identity());
}

/* ----------------------------------------------------------------------------- */

Matrix3D hex_pdbCryst(PdbImage *image, int cryst)
{
   if (cryst >= 1 && cryst <= image->n_cryst) return(image->mat_cryst[cryst-1]);
   else                                       return(tf_identity());
}

/* -------------------------------------------------------------------- */
// set the all residue masks in a PDB image (mainly just for reference)

void hex_maskPdb(PdbImage *pdb, int mask)
{
   for (int r=0; r<pdb->n_ca; r++) {

      pdb->residue[r].mask = mask;
   }
}

/* ----------------------------------------------------------------------------- */
// helpher function to cycle a chain label or other character

char deltaChar(char the_char, int delta)
{
   if (delta > 0) {

      if (the_char >= 'A' && the_char <= 'Z') {
   
         int x = the_char - 'A' + delta;
   
         if (x >= 26) x = x - (x/26)*26;
   
         the_char = 'A' + (char) x;
   
      } else if (the_char >= 'a' && the_char <= 'z') {
   
         int x = the_char - 'a' + delta;
   
         if (x >= 26) x = x - (x/26)*26;
   
         the_char = 'a' + (char) x;
   
      } else if (the_char >= '0' && the_char <= '9') {
    
         int x = the_char - '0' + delta;
   
         if (x >= 10) x = x - (x/10)*10;
   
         the_char = '0' + (char) x;
      }

   } else if (delta < 0) {  // change case

      if (the_char >= 'A' && the_char  <= 'Z') {

         int x = (the_char - 'A');

         the_char = 'a' + (char) x;

      } else if (the_char >= 'a' && the_char  <= 'z') {

         int x = (the_char - 'a');

         the_char = 'A' + (char) x;
      }
   }

   return(the_char);
}

/* ----------------------------------------------------------------------------- */
// check for valid atom coordinate values, assuming the usual PDB format (f8.3)

static int isBadPdbAtom(Point3D pt)
{
   int is_bad = isnan_pt(pt);

   if (is_bad == 0) {

      if      (pt.x < -1000.0 || pt.x > 10000.0) is_bad = 1;
      else if (pt.y < -1000.0 || pt.y > 10000.0) is_bad = 1;
      else if (pt.z < -1000.0 || pt.z > 10000.0) is_bad = 1;
   }
   
   return(is_bad);
}

/* ----------------------------------------------------------------------------- */
// helper function to re-label a chain according to a given delta

#if 1
static void set_chain(PdbImage *image, int a, char user_chain, char *the_chain)
{
   if (user_chain) the_chain[0] = user_chain;
   else            the_chain[0] = image->atom[a].chain[0];

   if (image->c_delta > 0) the_chain[0] = deltaChar(the_chain[0],
                                                    image->c_delta * image->n_chains);
   else                    the_chain[0] = deltaChar(the_chain[0],
                                                    image->c_delta * 1);
   the_chain[1] = '\0';
/*
   if (image->c_delta > 0) the_chain[0] = deltaChar(image->atom[a].chain[0], 
                                                    image->c_delta * image->n_chains);
   else                    the_chain[0] = deltaChar(image->atom[a].chain[0], 
                                                    image->c_delta * 1);
*/
}
#else
static void set_chain(PdbImage *image, int a, char *the_chain)
{
   strcpy(the_chain, image->atom[a].chain);

   if (image->c_delta > 0) {

      if (the_chain[0] >= 'A' && the_chain[0]  <= 'Z') {

         int x = (the_chain[0] - 'A') + (image->c_delta * image->n_chains);

         if (x >= 26) x = x - (x/26)*26;

         the_chain[0] = 'A' + (char) x;

      } else if (the_chain[0] >= 'a' && the_chain[0]  <= 'z') {

         int x = (the_chain[0] - 'a') + (image->c_delta * image->n_chains);

         if (x >= 26) x = x - (x/26)*26;

         the_chain[0] = 'a' + (char) x;

      } else if (the_chain[0] >= '0' && the_chain[0]  <= '9') {

         int x = (the_chain[0] - '0') + (image->c_delta * image->n_chains);

         if (x >= 10) x = x - (x/10)*10;

         the_chain[0] = '0' + (char) x;
      }

   } else if (image->c_delta < 0) {  // change case

      if (the_chain[0] >= 'A' && the_chain[0]  <= 'Z') {

         int x = (the_chain[0] - 'A');

         the_chain[0] = 'a' + (char) x;

      } else if (the_chain[0] >= 'a' && the_chain[0]  <= 'z') {

         int x = (the_chain[0] - 'a');

         the_chain[0] = 'A' + (char) x;
      }
   }
}
#endif

/* ----------------------------------------------------------------------------- */
// hack to handle ATOM records that overflow the usual 5-digit (I5) field
// 
// the atom name should normally start column 14, but it might contain a 
// leading digit in column 13. If these two columns are both blank, it almost
// certainly means that the atom number field has overflowed to the right.

static inline int bad_atom_hack(char *cline)
{
   int  shift = 0;
   char bad_atom[16];

   if ((strncmp(cline,"ATOM ",5) == 0) || (strncmp(cline,"HETATM ",7) == 0)) {

      if (cline[12] == ' ' && cline[13] == ' ') {

         hex_copy(cline, &bad_atom[0], 14);

         bad_atom[15] = '\0';

         hex_adv("<%s> -- overflows PDB format. Trying to continue\n", bad_atom);
 
         shift = 1;
      }
   }

   return(shift);
}

/* ----------------------------------------------------------------------------- */

static void hex_parsePdbAtom(char *line, int shift, PdbAtom *atom)
{
   char tmp[20];

/* we don't know these values yet: they will be set later */

   atom->mid      = 0;
   atom->cid      = 0;
   atom->rid      = 0;

   atom->is_ca    = 0;
   atom->is_mc    = 0;
   atom->is_bb    = 0;
   atom->is_sc    = 0;

   atom->radius = 0.0;

   str_field(line, 1, 6, 1, atom->record);

/* atom name without leading blanks */

   str_field(line, shift+13, shift+16, 1, atom->atomName); 

/* external name may have leading blanks */

   str_copy(atom->eName, &line[shift+13-1], 5);  

   atom->alt[0] = line[shift+17-1];
   atom->alt[1] = '\0';

   str_field(line, shift+18, shift+20, 1, atom->residue);

   atom->chain[0] = line[shift+22-1];
   atom->chain[1] = '\0';

// fix blank chain labels

   if (atom->chain[0] == ' ') atom->chain[0] = '_';

   str_field(line, shift+23, shift+26, 1, atom->sequence);

   atom->insert[0] = line[shift+27-1];
   atom->insert[1] = '\0';

   strcpy(atom->seqid, atom->sequence);

   if (atom->insert[0] != ' ') {

      strcat(atom->seqid, atom->insert);
   }

   str_field(line, shift+31, shift+38, 1, tmp);
   atom->pt.x = atof(tmp);

   str_field(line, shift+39, shift+46, 1, tmp);
   atom->pt.y = atof(tmp);

   str_field(line, shift+47, shift+54, 1, tmp);
   atom->pt.z = atof(tmp);

   str_field(line, shift+55,shift+60,1, tmp);
   atom->occ = (float) atof(tmp);

   if (atom->occ == 0.0) atom->occ = 1.0;       // default occupancy

   str_field(line, shift+61,shift+66,1, tmp);

   atom->tfac = (float) atof(tmp);
   if (atom->tfac == 0.0) atom->tfac = 99.99;  // default temperature factor
}

/* ----------------------------------------------------------------------------- */
// do the real work of reading a PDB file into an in-memory image 
//
// The strange counting here ensures that we keep the counts of residues and atoms 
// consistent, we do not allow more CA atoms than residues, and we do not allow 
// more than one alternate atom for each of the 5 "main chain" atoms. Multiple 
// side chain atoms can slip past these tests, but compared to the main chain atoms, 
// this is usually not a big problem.

PdbImage *hex_readPdb(char *filename, char *nickname)
{
   int       a, h, c, r, m, a_base, b_base, m_base, s_base, na, nr, nh, is_mc, n_sc;
   int       got_ca, got_cb, got_n, got_c, got_o, got_oxt, got_ot1, got_ot2;
   int       ir, im, n_items, n_vr;
   int       is_aa, is_dna, is_sugar, got_p, got_c3, got_c4, got_c5, shift;
   char      rec_type, sse_type;
   float     t1, t2, t3, t4, occ, tfac;
   Point3D   pt_ca;
   PdbImage *image;
   HexFile *file;
   char      atom_rec[] = "ATOM  ";
   char      pdb_file[MAX_PATHNAME], ctmp[MAX_PATHNAME], rname[MAX_PATHNAME];
   char      line[MAX_PATHNAME];
   char      this_chain[10], this_seqid[10], atomName[10];
   char      cnum[12], smtry[6], label[16], details[20];
   char      hname[6], this_hname[6], aname[6];
   

   image = (PdbImage *) kgetB_0(sizeof(PdbImage));

   image->mat_import = tf_identity();

   if (hex_filesize(filename) >= 0) strcpy(pdb_file, filename);
   else                             strcpy(pdb_file, hex_filename(filename,
                                                                  "pdb",ctmp));

   if ((file = hex_fopen(pdb_file, "r"))) {

      if (hex_debug > 2) hex_msg("Reading PDB: %s\n", pdb_file);

      image->name = kdupS(file->fname);

      if (nickname) image->nickname = kdupS(nickname);
      else          image->nickname = kdupS(hex_rootname(
                                            hex_basename(file->fname, ctmp), rname));

      image->residue = (PdbResidue *) kgetB_0(chunk_size*sizeof(PdbResidue));
      image->atom    = (PdbAtom *)    kgetB_0(chunk_size*sizeof(PdbAtom));
      image->hetero  = (PdbAtom *)    kgetB_0(chunk_size*sizeof(PdbAtom));

      na = chunk_size;
      nr = chunk_size;
      nh = chunk_size;

      a_base = b_base = m_base = s_base = 0;

      got_ca = got_n = got_c = got_o = got_oxt = got_ot1 = got_ot2 = got_cb = 0;

      got_p = got_c3 = got_c4 = got_c5 = 0;

      n_sc = 0; // no. side chain atoms

      n_vr = 0; // no. valid residues in current chain

      pt_ca  = pt_zero();

      strcpy(this_seqid, "none");
      strcpy(this_chain, "none");

      strcpy(this_hname, "none");

      image->n_cryst = 0;
      image->n_biomt = 0;

      image->mat_biomt = (Matrix3D *) kgetB_0(MAX_SYMMETRY*sizeof(Matrix3D));
      image->mat_cryst = (Matrix3D *) kgetB_0(MAX_SYMMETRY*sizeof(Matrix3D));

      m = c = r = a = h = 0;  // model = chain = residue = atom = hetero = 0

      int status = hex_freadln(file, MAX_PATHNAME, line);

      while (status != HEX_FILE_EOF) {

         if      (strncmp(line,"ATOM  ",6) == 0) rec_type = 'A';
         else if (strncmp(line,"HETATM",6) == 0) rec_type = 'H';
         else if (strncmp(line,"MODEL ",6) == 0) rec_type = 'M';
         else if (strncmp(line,"ENDMDL",6) == 0) rec_type = 'E';
         else if (strncmp(line,"CONECT",6) == 0) rec_type = 'C';
         else if (strncmp(line,"REMARK",6) == 0) rec_type = 'R';
         else if (strncmp(line,"ENERGY",6) == 0) rec_type = 'G';
         else if (strncmp(line,"TER   ",6) == 0) rec_type = 'T';
         else                                    rec_type = 'X';
/*
         if (strncmp(&line[13], "OXT", 3) == 0) {

            hex_msg("OXT here: %s\n", &line[0]);
         }
*/

         shift = bad_atom_hack(line);

// try to add non-standard residues that appear as HETATMs

         if (rec_type == 'H' && use_hetres) {

            str_field(line, shift+17, shift+19, 0, hname);

            char hcode = hex_aaCode(hname);

            if (hcode != 'X') { // this means we have a non-standard residue

               strcpy(aname, hex_aaName(hcode));

// so change the line to "ATOM" and carry on

               hex_msg("Changing <%s> to <ATOM ...>\n", &line[0]);

               hex_copy(atom_rec, &line[0], 6);

               rec_type = 'A';
/*
               if (strcmp(hname, this_hname) != 0) {

                  hex_msg("\n");
                  hex_wrn("Adding HETATM residue %s to alignment as %s\n", hname, aname);

                  strcpy(this_hname, hname);
               } 
*/

            } else {  // check for CA nonetheless

               str_field(line,shift+12,shift+15,0,aname);        

               if ((strcmp(aname, "CA") == 0) && (strcmp(hname, "CA") != 0)) {

                  if (use_warnings)  hex_wrn("Ignoring HETATM residue %s with %s atom\n",
                                             hname, aname);
               }
            }
         }

// now handle the different record types...

         if (rec_type == 'R') {  // /* pick out special stuff here */

            str_field(line,  14,18, 1, smtry);
   
            if (strcmp(smtry, "SMTRY") == 0) {
   
               str_field(line, 19,19, 1, cnum);
   
               ir = atoi(cnum) - 1;  /* matrix row number */
    
               str_field(line, 20,23, 1, cnum);
   
               im = atoi(cnum) - 1;  /* matrix number (0 should be identity) */
   
               n_items = sscanf(&line[23], "%f%f%f%f", &t1,&t2,&t3,&t4);
   
               if (ir >= 0 && ir <= 2 && im >= 0 && n_items == 4) {
   
                  if (im < MAX_SYMMETRY) {
   
                     image->mat_cryst[im].element[ir][0] = t1;
                     image->mat_cryst[im].element[ir][1] = t2;
                     image->mat_cryst[im].element[ir][2] = t3;
                     image->mat_cryst[im].element[ir][3] = t4;
   
                     image->n_cryst = max_(im+1, image->n_cryst);
   
                  } else {
   
                     if (use_warnings) hex_wrn("Too many SMTRY matrices in PDB file: %d\n", im);
                  }
   
               } else {
   
                  if (use_warnings) hex_wrn("Ignoring bad SMTRY record in PDB file:\n%s\n", line);
               }
   
            } else if (strcmp(smtry, "BIOMT") == 0) {
   
               str_field(line, 19,19, 1, cnum);
   
               ir = atoi(cnum) - 1;  /* matrix row number */
   
               str_field(line, 20,23, 1, cnum);
   
               im = atoi(cnum) - 1;     /* matrix number (0 should be identity) */
   
               n_items = sscanf(&line[23], "%f%f%f%f", &t1,&t2,&t3,&t4);
   
               if (ir >= 0 && ir <= 2 && im >= 0 && n_items == 4) {
   
                  if (im < MAX_SYMMETRY) {
   
                     image->mat_biomt[im].element[ir][0] = t1;
                     image->mat_biomt[im].element[ir][1] = t2;
                     image->mat_biomt[im].element[ir][2] = t3;
                     image->mat_biomt[im].element[ir][3] = t4;
   
                     image->n_biomt = max_(im+1, image->n_biomt);
   
                  } else {
   
                     if (use_warnings) hex_wrn("Too many BIOMT matrices in PDB file: %d\n", im);
                  }
   
               } else {
   
                  if (use_warnings) hex_wrn("Ignoring bad BIOMT record in PDB file:\n%s\n", line);
               }
            }

         } else if (rec_type == 'X') {   // non-standard or unexpected line

            if (strncmp(line,"I_sc",4) == 0) {

               image->energy_isc = atof(&line[5]);

//             hex_msg("I_sc = %.3f\n", image->energy_isc);
            }

         } else if (rec_type == 'M') {

            image->n_models += 1;

            str_field(line,  6,16, 1, cnum);

            m = atoi(cnum);

            c = 0;
            n_vr = 0;

            strcpy(this_seqid, "none");
            strcpy(this_chain, "none");

         } else if (rec_type == 'E') {

            if (one_model) break;

            strcpy(this_seqid, "none");
            strcpy(this_chain, "none");

         } else if (rec_type == 'T') {

            if (one_chain) break;

            strcpy(this_seqid, "none");
            strcpy(this_chain, "none");

         } else if (rec_type == 'A') {  

// normal case - handle ATOM data here

            if (!use_hydrogens) {
   
// perform a minimal parse to see the atom type
   
                str_field(line,shift+12,shift+15,0,atomName);        
   
                if (hex_isHydrogen(atomName)) goto pdb_foot;
            }
   
            if (a >= na) {
   
               na += chunk_size;
   
               image->atom = (PdbAtom *) kmod1B(image->atom, na*sizeof(PdbAtom));
            }
   
            hex_parsePdbAtom(line, shift, &image->atom[a]); 

            if (strip_alt) {

               if (image->atom[a].alt[0] != ' ' && image->atom[a].alt[0] != 'A') {

                  if (use_warnings) hex_wrn("Ignoring ALT atom %s:%s-%s:%s%s in PDB file %s\n", 
                                            image->atom[a].chain, image->atom[a].residue, 
                                            image->atom[a].seqid, image->atom[a].eName, 
                                            image->atom[a].alt, image->name);

                  goto pdb_foot;
               }
            }


            if (isBadPdbAtom(image->atom[a].pt)) {

               if (use_warnings) {

                  hex_wrn("Bad ATOM coordinates at line %d of file %s\n", 
                                         file->lineno, file->fname);
                  hex_msg("[%s]\n", line);
               }

               image->atom[a].pt = pt_zero();
            }

// NB. by default, all DNA and RNA structures are ignored

            is_aa  = hex_isAA (image->atom[a].residue);
            is_dna = hex_isDNA(image->atom[a].residue);

            if (select_chain == 0 && is_aa  == 0 && is_dna == 0) goto pdb_foot;
            if (select_chain == 1 && is_aa  == 0)                goto pdb_foot;
            if (select_chain == 2 && is_dna == 0)                goto pdb_foot;
/*
            if (is_aa  && select_chain == 2) goto pdb_foot;
            if (is_dna && select_chain == 1) goto pdb_foot;
*/

// is this the first atom of a new residue?
   
            if (strcmp(image->atom[a].seqid, this_seqid) != 0 ||
                strcmp(image->atom[a].chain, this_chain) != 0) {
   
               strcpy(this_seqid, image->atom[a].seqid);

               sse_type = 'D';

               if      (got_ca && got_c  && got_n  && got_o)  sse_type = 'X'; 
               else if (got_p  && got_c3 && got_c4 && got_c5) sse_type = 'X';

// save data for previous residue

               image->residue[r].n_sc = n_sc;
               image->residue[r].mask = 0;

               image->residue[r].mid = m;
               image->residue[r].cid = c;
   
               image->residue[r].atom1 = a_base;
               image->residue[r].atom2 = a - 1;

               strcpy(image->residue[r].chain, image->atom[a_base].chain);
               strcpy(image->residue[r].name,  image->atom[a_base].residue);
               strcpy(image->residue[r].label, image->atom[a_base].seqid);

               image->residue[r].sse[0] = sse_type;  // initially 'X' = UNKNOWN or 'D' = DEAD
               image->residue[r].sse[1] = '\0';
      
               if (hex_isAA(image->atom[a_base].residue)) {

                    image->residue[r].code[0] = hex_aaCode (image->atom[a_base].residue);
                    image->residue[r].code[1] = '\0';

               } else {

                    image->residue[r].code[0] = hex_dnaCode(image->atom[a_base].residue);
                    image->residue[r].code[1] = '\0';
               }

               if (hex_debug > 2) hex_msg(" Residue = %4d, atoms = [%d,%d]\n", r,
                                          image->residue[r].atom1, image->residue[r].atom2);

               if ((((got_ca && got_c  && got_n  && got_o)  || (got_ca && !discard_missing))) ||
                   (((got_p  && got_c3 && got_c4 && got_c5) || (got_p  && !discard_missing)))) {

//             if (got_ca || got_p) {

                  if (got_ca && use_warnings) {

                     if (!got_c) hex_wrn("Missing C atom for residue %s:%s-%s in PDB file %s\n", 
                                         image->residue[r].chain, image->residue[r].label, 
                                         image->residue[r].code, image->name);
                     if (!got_n) hex_wrn("Missing N atom for residue %s:%s-%s in PDB file %s\n", 
                                         image->residue[r].chain, image->residue[r].label, 
                                         image->residue[r].code, image->name);

                     if (!got_o && !(got_oxt || got_ot1 || got_ot2)) {

                        hex_wrn("Missing O atom for residue %s:%s-%s in PDB file %s\n", 
                                         image->residue[r].chain, image->residue[r].label, 
                                         image->residue[r].code, image->name);
                     }

                  } else if (got_p && use_warnings) {

                     if (!got_c3) hex_wrn("Missing C3 atom for base %s:%s-%s in PDB file %s\n", 
                                          image->residue[r].chain, image->residue[r].label, 
                                          image->residue[r].code, image->name);
                     if (!got_c4) hex_wrn("Missing C4 atom for base %s:%s-%s in PDB file %s\n", 
                                          image->residue[r].chain, image->residue[r].label, 
                                          image->residue[r].code, image->name);
                     if (!got_c5) hex_wrn("Missing C5 atom for base %s:%s-%s in PDB file %s\n", 
                                          image->residue[r].chain, image->residue[r].label, 
                                          image->residue[r].code, image->name);
                  }

// make the previous residue "official", and note the first atom of the new residue
   
                 r += 1;      

                 n_vr += 1;
   
                 a_base = a;
   
                 b_base = image->n_bb;
                 m_base = image->n_mc;
                 s_base = image->n_sc;
   
               } else {
   
// discard any accumulated atoms from the previous incomplete residue
   
                 if (hex_isAA(image->atom[a_base].residue)) {

                    if (got_ca || got_c || got_n || got_o || got_cb) {
   
                       if (use_warnings) hex_wrn("Ignoring incomplete residue: %s:%s-%s in file %s\n", 
                                                 image->atom[a_base].chain, image->atom[a_base].residue, 
                                                 image->atom[a_base].seqid, image->nickname);
                    }

                  } else {

                     if (got_p || got_c3 || got_c4 || got_c5) {

                       if (use_warnings) hex_wrn("Ignoring incomplete base: %s:%s-%s in file %s\n", 
                                                 image->atom[a_base].chain, image->atom[a_base].residue, 
                                                 image->atom[a_base].seqid, image->nickname);
                     }
                  }
   
                  a = a_base; 
   
                  image->n_bb = b_base;
                  image->n_mc = m_base;
                  image->n_sc = s_base;
   
// and re-parse the incoming line into the current atom slot
   
                  hex_parsePdbAtom(line, shift, &image->atom[a]); 

                  if (isBadPdbAtom(image->atom[a].pt)) {
      
                     if (use_warnings) {

                        hex_wrn("Bad ATOM coordinates at line %d of file %s\n", file->lineno, file->fname);
                        hex_msg("[%s]\n", line);
                     }
      
                     image->atom[a].pt = pt_zero();
                  }
               }
   
// reset side chain and main chain counters for the next residue

               pt_ca  = pt_zero();
   
               n_sc = 0;
   
               got_ca = got_n = got_c = got_o = got_oxt = got_ot1 = got_ot2 = got_cb = 0;

               got_p = got_c3 = got_c4 = got_c5 = 0;
   
               if (r+1 >= nr) {
   
                  nr += chunk_size;
   
                  image->residue = (PdbResidue *) kmod1B(image->residue,
                                                             nr*sizeof(PdbResidue));
               }
   
               image->residue[r].n_sc   = 0;
            }
   
// when we get a new chain, make a new chain id and reset the count of valid residues
   
            if (strcmp(image->atom[a].chain, this_chain) != 0) {
   
               if (n_vr > 0) {

                  c += 1;

                  n_vr = 0;
               }

               strcpy(this_chain, image->atom[a].chain);
            }

// now save the info about the current atom

            if (is_aa) {

               is_mc = hex_isMainChain(image->atom[a].atomName);
      
               if (is_mc == PDB_N && !got_n) {            // N = main chain + backbone
      
                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
      
                  image->atom[a].is_bb = 1;
                  image->atom[a].is_mc = 1;
      
                  image->n_bb += 1; 
                  image->n_mc += 1; 
      
                  got_n = 1; 
                  a++; 
      
               } else if (is_mc == PDB_A && !got_ca) {   // CA = main chain + backbone
      
                  pt_ca = image->atom[a].pt;
      
                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
   
                  image->atom[a].is_bb = 1;
                  image->atom[a].is_mc = 1;
                  image->atom[a].is_ca = 1;
   
                  image->n_bb += 1; 
                  image->n_mc += 1; 
      
                  got_ca = 1; 
                  a++; 
      
               } else if (is_mc == PDB_C && !got_c) {     // C = main chain + backbone
      
                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
      
                  image->n_bb += 1; 
                  image->n_mc += 1; 
      
                  image->atom[a].is_bb = 1;
                  image->atom[a].is_mc = 1;
   
                  got_c = 1; 
                  a++; 
      
               } else if (is_mc == PDB_O && !got_o) {     // N = main chain + backbone 
      
                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
      
                  image->atom[a].is_bb = 1;
                  image->atom[a].is_mc = 1;
   
                  image->n_bb += 1; 
                  image->n_mc += 1; 
      
                  got_o = 1; 
                  a++; 

               } else if (is_mc == PDB_OXT && !got_oxt) {     // N = main chain + backbone 
      
                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
      
                  image->atom[a].is_bb = 1;
                  image->atom[a].is_mc = 1;
   
                  image->n_bb += 1; 
                  image->n_mc += 1; 
      
                  got_oxt = 1; 
                  a++; 
      
               } else if (is_mc == PDB_OT1 && !got_ot1) {     // N = main chain + backbone 
      
                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
      
                  image->atom[a].is_bb = 1;
                  image->atom[a].is_mc = 1;
   
                  image->n_bb += 1; 
                  image->n_mc += 1; 
      
                  got_ot1 = 1; 
                  a++; 
      
               } else if (is_mc == PDB_OT2 && !got_ot2) {     // N = main chain + backbone 
      
                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
      
                  image->atom[a].is_bb = 1;
                  image->atom[a].is_mc = 1;
   
                  image->n_bb += 1; 
                  image->n_mc += 1; 
      
                  got_ot2 = 1; 
                  a++; 
      
               } else if (is_mc == PDB_B && !got_cb) {   // CB = main chain + sidechain
      
                  n_sc += 1;
      
                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
      
                  image->atom[a].is_mc = 1;
                  image->atom[a].is_sc = 1;
   
                  image->n_mc += 1; 
                  image->n_sc += 1;
      
                  got_cb = 1; 
                  a++; 
      
               } else if (hex_isSideChain(image->atom[a].atomName)) {
      
                  n_sc += 1;
      
                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
   
                  image->atom[a].is_sc = 1;
   
                  image->n_sc += 1;
      
                  a++;
               }

            } else {  // its DNA or RNA

               is_sugar = hex_isSugar(image->atom[a].atomName);

               if (is_sugar == 3 && !got_c3) {
               
                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
      
                  image->atom[a].is_mc = 1;
                  image->atom[a].is_bb = 1;
      
                  image->n_bb += 1; 
                  image->n_mc += 1; 
      
                  got_c3 = 1; 
                  a++; 

               } else if (is_sugar == 4 && !got_c4) {

                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
      
                  image->atom[a].is_mc = 1;
                  image->atom[a].is_bb = 1;
      
                  image->n_bb += 1; 
                  image->n_mc += 1; 
      
                  got_c4 = 1; 
                  a++; 

               } else if (is_sugar == 5 && !got_c5) {

                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
      
                  image->atom[a].is_mc = 1;
                  image->atom[a].is_bb = 1;
      
                  image->n_bb += 1; 
                  image->n_mc += 1; 
      
                  got_c5 = 1; 
                  a++; 

               } else if (is_sugar == 9 && !got_p) {
               
                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
      
                  image->atom[a].is_mc = 1;
                  image->atom[a].is_bb = 1;
      
                  image->n_bb += 1; 
                  image->n_mc += 1; 
      
                  got_p = 1; 
                  a++; 

               } else if ((is_sugar == 1 || is_sugar == 2) || is_sugar > 5) {

                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
      
                  image->atom[a].is_mc = 1;
      
                  image->n_mc += 1; 
      
                  a++; 

               } else {  // assume its a base

                  n_sc += 1;
      
                  image->atom[a].mid = m;
                  image->atom[a].rid = r;
                  image->atom[a].cid = c;
   
                  image->atom[a].is_sc = 1;
   
                  image->n_sc += 1;
      
                  a++;
               }
            }

         } else if (rec_type == 'H') {  

// hetero atoms are stored in a separate list of atoms

            if (!use_hydrogens) {
   
// perform a minimal parse to see the atom type
   
                str_field(line,shift+12,shift+15,0,atomName);        
   
                if (hex_isHydrogen(atomName)) goto pdb_foot;
            }
   
            if (h >= nh) {
   
               nh += chunk_size;
   
               image->hetero = (PdbAtom *) kmod1B(image->hetero, nh*sizeof(PdbAtom));
            }
   
            hex_parsePdbAtom(line, shift, &image->hetero[h]); 

            if (isBadPdbAtom(image->hetero[h].pt)) {

               if (use_warnings) {
                  hex_wrn("Bad HETATM coordinates at line %d of file %s\n", file->lineno, file->fname);
                  hex_msg("[%s]\n", line);
               }

               image->hetero[h].pt = pt_zero();
            }

// NB. just collect up all hetero atoms without testing their type

/*
            is_aa  = hex_isAA (image->hetero[h].residue);
            is_dna = hex_isDNA(image->hetero[h].residue);

            if (select_chain == 0 && is_aa  == 0 && is_dna == 0) goto pdb_foot;
            if (select_chain == 1 && is_aa  == 0)                goto pdb_foot;
            if (select_chain == 2 && is_dna == 0)                goto pdb_foot;
*/

// NB. by default, all DNA and RNA structures are ignored (these tests should fail anyway for HETATM)
/*
            if (is_aa  && select_chain == 2) goto pdb_foot;
            if (is_dna && select_chain == 1) goto pdb_foot;
*/

// is this the first atom of a new residue?
   
            if (strcmp(image->hetero[h].seqid, this_seqid) != 0 ||
                strcmp(image->hetero[h].chain, this_chain) != 0) {
   
               strcpy(this_seqid, image->hetero[h].seqid);

               sse_type = 'D';

               if      (got_ca && got_c  && got_n  && got_o)  sse_type = 'X'; 
               else if (got_p  && got_c3 && got_c4 && got_c5) sse_type = 'X';

// save data for previous residue

               image->residue[r].mask = 0;
               image->residue[r].n_sc = n_sc;

               image->residue[r].mid = m;
               image->residue[r].cid = c;
   
               image->residue[r].atom1 = a_base;
               image->residue[r].atom2 = a - 1;

               strcpy(image->residue[r].chain, image->atom[a_base].chain);
               strcpy(image->residue[r].name,  image->atom[a_base].residue);
               strcpy(image->residue[r].label, image->atom[a_base].seqid);

               image->residue[r].sse[0] = sse_type;  // initially 'X' = UNKNOWN or 'D' = DEAD
               image->residue[r].sse[1] = '\0';
      
               if (hex_isAA(image->atom[a_base].residue)) {

                    image->residue[r].code[0] = hex_aaCode (image->atom[a_base].residue);
                    image->residue[r].code[1] = '\0';

               } else {

                    image->residue[r].code[0] = hex_dnaCode(image->atom[a_base].residue);
                    image->residue[r].code[1] = '\0';
               }


              if (hex_debug > 2) hex_msg(" Residue = %4d, atoms = [%d,%d]\n", r,
                                         image->residue[r].atom1, image->residue[r].atom2);

//             if (got_ca || got_p) {

               if ((((got_ca && got_c  && got_n  && got_o)  || (got_ca && !discard_missing))) ||
                   (((got_p  && got_c3 && got_c4 && got_c5) || (got_p  && !discard_missing)))) {

                  if (got_ca && use_warnings) {

                     if (!got_c) hex_wrn("Missing C atom for residue %s:%s-%s in PDB file %s\n", 
                                         image->residue[r].chain, image->residue[r].label, 
                                         image->residue[r].code, image->name);
                     if (!got_n) hex_wrn("Missing N atom for residue %s:%s-%s in PDB file %s\n", 
                                         image->residue[r].chain, image->residue[r].label, 
                                         image->residue[r].code, image->name);

                     if (!got_o && !(got_oxt || got_ot1 || got_ot2)) {

                        hex_wrn("Missing O atom for residue %s:%s-%s in PDB file %s\n", 
                                         image->residue[r].chain, image->residue[r].label, 
                                         image->residue[r].code, image->name);
                      }

                  } else if (got_p && use_warnings) {

                     if (!got_c3) hex_wrn("Missing C3 atom for base %s:%s-%s in PDB file %s\n", 
                                          image->residue[r].chain, image->residue[r].label, 
                                          image->residue[r].code, image->name);
                     if (!got_c4) hex_wrn("Missing C4 atom for base %s:%s-%s in PDB file %s\n", 
                                          image->residue[r].chain, image->residue[r].label, 
                                          image->residue[r].code, image->name);
                     if (!got_c5) hex_wrn("Missing C5 atom for base %s:%s-%s in PDB file %s\n", 
                                          image->residue[r].chain, image->residue[r].label, 
                                          image->residue[r].code, image->name);
                  }

// make the previous residue "official", and note the first atom of the new residue
   
                 r += 1;      

                 n_vr += 1;
   
                 a_base = a;
   
                 b_base = image->n_bb;
                 m_base = image->n_mc;
                 s_base = image->n_sc;
   
               } else {
   
// discard any accumulated atoms from the previous incomplete residue
   
                 if (hex_isAA(image->atom[a_base].residue)) {

                    if (got_ca || got_c || got_n || got_o || got_oxt || got_cb) {
   
                       if (use_warnings) hex_wrn("Ignoring incomplete residue: %s:%s-%s in file %s\n", 
                                                 image->atom[a_base].chain, image->atom[a_base].residue, 
                                                 image->atom[a_base].seqid, image->nickname);
                    }

                  } else {

                     if (got_p || got_c3 || got_c4 || got_c5) {

                       if (use_warnings) hex_wrn("Ignoring incomplete base: %s:%s-%s in file %s\n", 
                                                 image->atom[a_base].chain, image->atom[a_base].residue, 
                                                 image->atom[a_base].seqid, image->nickname);
                     }
                  }
   
                  a = a_base; 
   
                  image->n_bb = b_base;
                  image->n_mc = m_base;
                  image->n_sc = s_base;
   
// and re-parse the incoming line into the current atom slot
   
                  hex_parsePdbAtom(line, shift, &image->hetero[h]); 

                  if (isBadPdbAtom(image->hetero[h].pt)) {
      
                     if (use_warnings) {
                        hex_wrn("Bad HETATM coordinates at line %d of file %s\n", file->lineno, file->fname);
                        hex_msg("[%s]\n", line);
                     }
      
                     image->hetero[h].pt = pt_zero();
                  }
               }
   
// reset side chain and main chain counters for the next residue

               pt_ca  = pt_zero();
   
               n_sc = 0;
   
               got_ca = got_n = got_c = got_o = got_oxt = got_ot1 = got_ot2 = got_cb = 0;

               got_p = got_c3 = got_c4 = got_c5 = 0;
   
               if (r+1 >= nr) {
   
                  nr += chunk_size;
   
                  image->residue = (PdbResidue *) kmod1B(image->residue,
                                                             nr*sizeof(PdbResidue));
               }
   
               image->residue[r].n_sc   = 0;
            }
   
// when we get a new chain, make a new chain id and reset the count of valid residues
   
#if 1
            if (strcmp(image->hetero[h].chain, this_chain) != 0) {
   
               strcpy(this_chain, image->hetero[h].chain);

               if (n_vr > 0) {

                  c += 1;

                  n_vr = 0;
               }
            }
#endif

// now save the info about the current atom 

            image->hetero[h].mid = m;
            image->hetero[h].rid = r;
            image->hetero[h].cid = c;

            h++;
         }

pdb_foot: 
         status = hex_freadln(file, MAX_PATHNAME, line);
      }  // next line of PDB file

//  finally, either accept the last residue or throw away any incomplete atoms

      sse_type = 'D';

      if      (got_ca && got_c  && got_n  && got_o)  sse_type = 'X'; 
      else if (got_p  && got_c3 && got_c4 && got_c5) sse_type = 'X';

// save the data for the last residue

         image->residue[r].n_sc = n_sc;
         image->residue[r].mask = 0;

         image->residue[r].mid = m;
         image->residue[r].cid = c;

         image->residue[r].atom1 = a_base;
         image->residue[r].atom2 = a - 1;

         strcpy(image->residue[r].chain, image->atom[a_base].chain);
         strcpy(image->residue[r].name,  image->atom[a_base].residue);
         strcpy(image->residue[r].label, image->atom[a_base].seqid);

         image->residue[r].sse[0] = sse_type;
         image->residue[r].sse[1] = '\0';
      
         if (hex_isAA(image->atom[a_base].residue)) {
            image->residue[r].code[0] = hex_aaCode (image->atom[a_base].residue);
            image->residue[r].code[1] = '\0';
         } else {
            image->residue[r].code[0] = hex_dnaCode(image->atom[a_base].residue);
            image->residue[r].code[1] = '\0';
         }

      if ((((got_ca && got_c  && got_n  && got_o)  || (got_ca && !discard_missing))) ||
          (((got_p  && got_c3 && got_c4 && got_c5) || (got_p  && !discard_missing)))) {

//    if (got_ca || got_p) {

         if (got_ca && use_warnings) {
   
            if (!got_c) hex_wrn("Missing C atom for residue %s:%s-%s in PDB file %s\n", 
                                image->residue[r].chain, image->residue[r].label, 
                                image->residue[r].code, image->name);
            if (!got_n) hex_wrn("Missing N atom for residue %s:%s-%s in PDB file %s\n", 
                                image->residue[r].chain, image->residue[r].label, 
                                image->residue[r].code, image->name);

            if (!got_o && !(got_oxt || got_ot1 || got_ot2)) {

               hex_wrn("Missing O atom for residue %s:%s-%s in PDB file %s\n", 
                                image->residue[r].chain, image->residue[r].label, 
                                image->residue[r].code, image->name);
            }

         } else if (got_p && use_warnings) {
   
            if (!got_c3) hex_wrn("Missing C3 atom for base %s:%s-%s in PDB file %s\n", 
                                 image->residue[r].chain, image->residue[r].label, 
                                 image->residue[r].code, image->name);
            if (!got_c4) hex_wrn("Missing C4 atom for base %s:%s-%s in PDB file %s\n", 
                                 image->residue[r].chain, image->residue[r].label, 
                                 image->residue[r].code, image->name);
            if (!got_c5) hex_wrn("Missing C5 atom for base %s:%s-%s in PDB file %s\n", 
                                 image->residue[r].chain, image->residue[r].label, 
                                 image->residue[r].code, image->name);
         }

         if (hex_debug > 2) hex_msg(" Residue = %4d, atoms = [%d,%d]\n", r,
                                      image->residue[r].atom1, image->residue[r].atom2);

         r += 1;      

         n_vr += 1;

      } else {  

// ignore incomplete last residue - i.e. reset to last good atom

         if (hex_isAA(image->atom[a_base].residue)) {

            if (got_ca || got_c || got_n || got_o || got_oxt || got_cb) {
   
               if (use_warnings) hex_wrn("Ignoring incomplete residue: %s:%s-%s in file %s\n", 
                                         image->atom[a_base].chain, image->atom[a_base].residue, 
                                         image->atom[a_base].seqid, image->nickname);
            }

         } else {

            if (got_p || got_c3 || got_c4 || got_c5) {

            if (use_warnings) hex_wrn("Ignoring incomplete base: %s:%s-%s in file %s\n", 
                                      image->atom[a_base].chain, image->atom[a_base].residue, 
                                      image->atom[a_base].seqid, image->nickname);
            }
         }
   
         a = a_base;
      }

// we now have complete backbone atom sets and absolutely consistent atom counts

      image->n_ca     = r;
      image->n_atoms  = a;
      image->n_hetero = h;

      if (n_vr > 0) image->n_chains = c + 1;
      else          image->n_chains = c;

      image->c_delta  = 0;

      if (save_space) {

         image->residue = (PdbResidue *) kmod1B(image->residue,
                                                image->n_ca*sizeof(PdbResidue));

         image->atom    = (PdbAtom *)    kmod1B(image->atom,
                                                image->n_atoms*sizeof(PdbAtom));

         image->hetero  = (PdbAtom *)    kmod1B(image->hetero,
                                                image->n_hetero*sizeof(PdbAtom));
      }

// make a copy of the atom coordinates to support a reset after flexible fitting

      image->atom_backup   = (Point3D *) kgetB(image->n_atoms *sizeof(Point3D));
      image->hetero_backup = (Point3D *) kgetB(image->n_hetero*sizeof(Point3D));

      for (a=0; a<image->n_atoms; a++) {

         image->atom_backup[a] = image->atom[a].pt;
      }

      for (h=0; h<image->n_hetero; h++) {

         image->hetero_backup[h] = image->hetero[h].pt;
      }

// make sure we have at least one well-defined symmetry matrix of each type

      if (image->n_cryst == 0) {

         image->mat_cryst[0] = tf_identity();

         image->n_cryst = 1;

         if (save_space) {

            image->mat_cryst = (Matrix3D *) kmod1B(image->mat_cryst, 
                                                   image->n_cryst*sizeof(Matrix3D));
         }
      }

      if (image->n_biomt == 0) {

         image->mat_biomt[0] = tf_identity();

         image->n_biomt = 1;

         if (save_space) {

            image->mat_biomt = (Matrix3D *) kmod1B(image->mat_biomt, 
                                                   image->n_biomt*sizeof(Matrix3D));
         }
      }

      hex_fclose(file);

   } else {

      hex_err("Failed to open PDB file: %s\n", pdb_file);
   }

   return(image);
}

/* ----------------------------------------------------------------------------- */

PdbImage *hex_copyPdb(PdbImage *image)
{
   if (image) {

      PdbImage *new_image = (PdbImage *) kgetB_0(sizeof(PdbImage));

      hex_copy(image, new_image, sizeof(PdbImage));

      new_image->name     = kdupS(image->name);
      new_image->nickname = kdupS(image->nickname);

      new_image->residue = (PdbResidue *) kgetB(image->n_ca*sizeof(PdbResidue));
      new_image->atom    = (PdbAtom *)    kgetB(image->n_atoms*sizeof(PdbAtom));
      new_image->hetero  = (PdbAtom *)    kgetB(image->n_hetero*sizeof(PdbAtom));

      new_image->atom_backup   = (Point3D *) kgetB(image->n_atoms *sizeof(Point3D));
      new_image->hetero_backup = (Point3D *) kgetB(image->n_hetero*sizeof(Point3D));

      new_image->mat_cryst = (Matrix3D *) kgetB(image->n_cryst*sizeof(Matrix3D));
      new_image->mat_biomt = (Matrix3D *) kgetB(image->n_biomt*sizeof(Matrix3D));

      hex_copy(image->residue, new_image->residue, image->n_ca*sizeof(PdbResidue));
      hex_copy(image->atom,    new_image->atom,    image->n_atoms*sizeof(PdbAtom));
      hex_copy(image->hetero,  new_image->hetero,  image->n_hetero*sizeof(PdbAtom));

      hex_copy(image->atom_backup,   new_image->atom_backup,   image->n_atoms *sizeof(Point3D));
      hex_copy(image->hetero_backup, new_image->hetero_backup, image->n_hetero*sizeof(Point3D));

      hex_copy(image->mat_cryst, new_image->mat_cryst, image->n_cryst*sizeof(Matrix3D));
      hex_copy(image->mat_biomt, new_image->mat_biomt, image->n_biomt*sizeof(Matrix3D));

      return(new_image);
   }

   return(NULL);
}

/* -------------------------------------------------------------------- */

static int lookup_pos(int test, int n_list, int *list)
{
   for (int i=0; i<n_list; i++) if (test == list[i]) return(i);

   return(-1);
}

/* -------------------------------------------------------------------- */
// make a new PdbImage from an existing one by extracting the residues
// and atoms that match a given list of chains (this function never 
// copies hetero atoms, so the new image will not have any hetero atoms).

// NB. hex_pdb works with incremental chain numbers, not chain letters,
// so we can have several "A" chains, but each will have a different 
// internal chain number (and probably a different model number).

PdbImage *hex_extractPdb(PdbImage *pdb, int n_chains, int *the_chains)
{
   int  c, r, a, na, nr, j, pos, r_base, a_base;

   na = chunk_size;
   nr = chunk_size;

   PdbImage *image = (PdbImage *)   kgetB_0(sizeof(PdbImage));

   image->atom        = (PdbAtom *)    kgetB_0(chunk_size*sizeof(PdbAtom));
   image->residue     = (PdbResidue *) kgetB_0(chunk_size*sizeof(PdbResidue));
   image->atom_backup = (Point3D *)    kgetB_0(chunk_size*sizeof(Point3D));

   if (pdb->name)     image->name     = kdupS(pdb->name);
   if (pdb->nickname) image->nickname = kdupS(pdb->nickname);

   image->c_delta  = pdb->c_delta;

   image->n_cryst  = 1;
   image->n_biomt  = 1;

   image->mat_cryst = (Matrix3D *) kgetB_0(1*sizeof(Matrix3D));
   image->mat_biomt = (Matrix3D *) kgetB_0(1*sizeof(Matrix3D));

   image->mat_cryst[0] = tf_identity();
   image->mat_biomt[0] = tf_identity();

   image->mat_import = pdb->mat_import;

// copy out those atoms that match the given chains

   int *found_chains = kget1I(n_chains);

   for (j=0; j<n_chains; j++) found_chains[j] = -1;

   for (j=0,a=0,c=0; j<pdb->n_atoms; j++) {

      pos = lookup_pos(pdb->atom[j].cid, n_chains, the_chains);

      if (pos >= 0) {

         if (found_chains[pos] == -1) {

            image->n_chains++;

            if (image->n_chains > 1) c++;

            found_chains[pos] = c;
         }

         if (a >= na) {
   
            na += chunk_size;
   
            image->atom        = (PdbAtom *) kmod1B(image->atom,        na*sizeof(PdbAtom));
            image->atom_backup = (Point3D *) kmod1B(image->atom_backup, na*sizeof(Point3D));
         }

         image->atom[a]        = pdb->atom[j];
         image->atom_backup[a] = pdb->atom_backup[j];

         if (image->atom[a].is_ca) image->n_ca++;
         if (image->atom[a].is_bb) image->n_bb++;
         if (image->atom[a].is_mc) image->n_mc++;
         if (image->atom[a].is_sc) image->n_sc++;

// the extracted chain numbers count from zero in the new image

         image->atom[a].cid = c;

         a++;
      }
   }

   image->n_atoms = a;

// similarly, copy out the residues that match the given chains

   for (j=0,r=0; j<pdb->n_ca; j++) {

      pos = lookup_pos(pdb->residue[j].cid, n_chains, the_chains);

      if (pos >= 0) {

         if (r >= nr) {
   
            nr += chunk_size;
   
            image->residue = (PdbResidue *) kmod1B(image->residue,
                                                   nr*sizeof(PdbResidue));
         }
   
         image->residue[r] = pdb->residue[j];

         image->residue[r].cid = found_chains[pos];

         r++;
      }
   }

   if (image->n_chains != n_chains && use_warnings) {

      hex_wrn("hex_extractPdb(): no. chains extracted (%d) not equal to no. requested (%d).\n", 
              c, n_chains);
   }

   if (image->n_ca != r && use_warnings) {

      hex_wrn("hex_extractPdb(): no. residues (%d) not equal to No. CA atoms (%d).\n", 
              r, image->n_ca);

      image->n_ca = min_(image->n_ca, r);
   }

// now re-number the internal residue and atom_numbers to count from zero

   r_base = image->atom[0].rid;
   a_base = image->residue[0].atom1;

   for (j=0; j<image->n_atoms; j++) {

      image->atom[j].rid -= r_base;
   }

   for (j=0; j<image->n_ca; j++) {

      image->residue[j].atom1 -= a_base;
      image->residue[j].atom2 -= a_base;
   }

   if (hex_debug > 0) hex_msg("Extracted %d chains, %d residues, %d atoms, "
                              "CA=%d, BB=%d, MC=%d, SC=%d\n", image->n_chains, r, 
                              image->n_atoms, image->n_ca, image->n_bb, image->n_mc, image->n_sc);

   if (save_space) {

      image->atom        = (PdbAtom *)     kmod1B(image->atom,
                                                  image->n_atoms*sizeof(PdbAtom));

      image->residue     = (PdbResidue *)  kmod1B(image->residue,
                                                  image->n_ca*sizeof(PdbResidue));

      image->atom_backup = (Point3D *)     kmod1B(image->atom_backup,
                                                  image->n_atoms*sizeof(Point3D));
   }

   kwipe(found_chains);

   return(image);
}

/* -------------------------------------------------------------------- */
// split a multi-model or multi-chain PdbImage into one or more new 
// images, each one containing just one model or just one chain. 
// It is assumed that the caller already knows (from image->n_models and
// image->n_chains) how many chains are in the image, and which ones he 
// wants to extract (usually all of them).

// NB. hex_pdb works with incremental chain numbers, not chain letters,
// so we can have several "A" chains, but each will have a different 
// internal chain number (and probably a different model number).

void hex_splitPdb(PdbImage *pdb, int n_chains, int *the_chains, 
                  char ** the_names, PdbImage **images)
{
   if (pdb->n_models <= 1) {  // one model, several chains

      for (int i=0; i<n_chains; i++) {
   
         images[i] = hex_extractPdb(pdb, 1, &the_chains[i]);
   
         if (images[i]) {
   
            kwipe(images[i]->nickname);
   
            images[i]->nickname = kdupS(the_names[i]);
         }
      }

   } else {  // several models, one or more chains in each model

      int n_cpm = pdb->n_chains / pdb-> n_models;  // chains per model

      for (int i=0; i<n_chains; i+=n_cpm) {
      
         images[i] = hex_extractPdb(pdb, n_cpm, &the_chains[i]);
      
         if (images[i]) {
      
            kwipe(images[i]->nickname);
      
            images[i]->nickname = kdupS(the_names[i]);
         } 
      }
   }
}

/* -------------------------------------------------------------------- */
// find the chain numbers that first match the given chain letters

int hex_chainsPdb(PdbImage *pdb, char *chains, int *the_chains)
{
   int found = 0;

   for (int n=0; n<strlen(chains); n++) {

      for (int r=0; r<pdb->n_ca; r++) if (pdb->residue[r].mid == 0) {

         if (pdb->residue[r].chain[0] == chains[n]) {

            the_chains[found++] = pdb->residue[r].cid;

            break;
         }
      }
   }

   return(found);
}

/* -------------------------------------------------------------------- */
// shift a PDB image to the origin, and record the inverse transform
// (NB. This does NOT change the SYMTRY or BIOMT matrices)

double hex_importPdb(PdbImage *pdb, int use_random)
{
   double         radius;
   Euler          e;
   Point3D        pt_origin;
   Vector3D       v_move;
   Matrix3D       t_random;
   Point3D       *pts_all;

   radius = 0.0;

   if (pdb->n_atoms > 0) {

      pts_all = hex_ptsPdb(pdb, 0);

      pt_origin = pt_average(pdb->n_atoms, pts_all);

      for (int a=0; a<pdb->n_atoms; a++) {

         radius = fmax(d_ptpt(pts_all[a],pt_origin), radius);
      }

      v_move = v_make(-pt_origin.x, -pt_origin.y, -pt_origin.z);

      hex_movePdb(pdb, tf_v(v_move));

      pdb->mat_import = tf_tftf(pdb->mat_import, tf_v(v_move));

      if (use_random) {

         t_random = tf_random();

         if (hex_debug > 1) {

           e = euler_tf(t_random);

           hex_msg("Applying random rotation to %s = [%.2lf,%.2lf,%.2lf]\n",
                   pdb->nickname, 
                   rad2deg(e.alpha), rad2deg(e.beta), rad2deg(e.gamma));
         }
        
         hex_movePdb(pdb, t_random);

         pdb->mat_import = tf_tftf(pdb->mat_import, t_random);
      }

      kfree(pts_all);
   }

   return(radius);
}

/* ----------------------------------------------------------------------------- */

void hex_freePdb(PdbImage *image)
{
   if (image) {

      kwipe(image->name);
      kwipe(image->nickname);
      kwipe(image->hetero);
      kwipe(image->atom);
      kwipe(image->residue);

      kwipe(image->atom_backup);
      kwipe(image->hetero_backup);

      kwipe(image->mat_biomt);
      kwipe(image->mat_cryst);

      kwipe(image);
   }
}

/* ----------------------------------------------------------------------------- */
// transform all of the atoms in a PDB image

void hex_movePdb(PdbImage *image, Matrix3D tf_pdb)
{
   for (int a=0; a<image->n_atoms; a++) {

      image->atom[a].pt = pt_tfpt(tf_pdb, image->atom[a].pt);
   }

   for (int h=0; h<image->n_hetero; h++) {

      image->hetero[h].pt = pt_tfpt(tf_pdb, image->hetero[h].pt);
   }
}

/* ----------------------------------------------------------------------------- */
// reset all the coordinates in a PDB image back to their original values

void hex_resetPdb(PdbImage *image)
{
   for (int a=0; a<image->n_atoms; a++) {

      image->atom[a].pt = image->atom_backup[a];
   }

   for (int h=0; h<image->n_hetero; h++) {

      image->hetero[h].pt = image->hetero_backup[h];
   }
}

/* ----------------------------------------------------------------------------- */
// freeze the current coordinates, as if they were loaded that way

void hex_freezePdb(PdbImage *image)
{
   image->mat_import = tf_identity();

   for (int a=0; a<image->n_atoms; a++) {

      image->atom_backup[a] = image->atom[a].pt;
   }

   for (int h=0; h<image->n_hetero; h++) {

      image->hetero_backup[h] = image->hetero[h].pt;
   }
}

/* ----------------------------------------------------------------------------- */
// similarly, transform the atoms of a single residue

void hex_movePdbResidue(PdbImage *image, int rid, Matrix3D tf_pdb)
{
   if (rid >= 0 && rid < image->n_ca) {

      int atom1 = image->residue[rid].atom1;
      int atom2 = image->residue[rid].atom2;

      for (int a=atom1; a<=atom2; a++) {

         image->atom[a].pt = pt_tfpt(tf_pdb, image->atom[a].pt);
      }
   }
}

/* ----------------------------------------------------------------------------- */
// similarly, change the chain label of a single residue

void hex_labelPdbResidue(PdbImage *image, int rid, int delta)
{
   if (rid >= 0 && rid < image->n_ca) {

      image->residue[rid].chain[0] = deltaChar(image->residue[rid].chain[0], delta);

      int atom1 = image->residue[rid].atom1;
      int atom2 = image->residue[rid].atom2;

      for (int a=atom1; a<=atom2; a++) {

         image->atom[a].chain[0] = deltaChar(image->atom[a].chain[0], delta);
      }
   }
}

/* ----------------------------------------------------------------------------- */
// change the chain label of all the atoms in a PDB image. Doing this is 
// surprisingly expensive, so it should be avoided if all all possible.

void hex_setPdbChainLabel(PdbImage *image, char new_chain)
{
   for (int rid=0; rid<image->n_ca; rid++) {

      image->residue[rid].chain[0] = new_chain;

      int atom1 = image->residue[rid].atom1;
      int atom2 = image->residue[rid].atom2;

      for (int a=atom1; a<=atom2; a++) {

         image->atom[a].chain[0] = new_chain;
      }
   }
}

/* ----------------------------------------------------------------------------- */
// extract the CA coordinates for a given residue

Point3D hex_ptPdbResidue(PdbImage *image, int rid)
{
   Point3D pt = pt_zero();

   if (rid >= 0 && rid < image->n_ca) {

      int atom1 = image->residue[rid].atom1;
      int atom2 = image->residue[rid].atom2;

      for (int a=atom1; a<=atom2; a++) {

         if (hex_isAlphaCarbon(image->atom[a].atomName)) return(image->atom[a].pt);
      }
   }

   return(pt);
}


/* ----------------------------------------------------------------------------- */

void hex_writePdb(char *filename, PdbImage *image, int chain_id, char chain_label, 
                  int atom_mode, int with_sse, int with_ter, Matrix3D t_export)
{
   HexFile  *file = hex_openPdbFile(filename);

   if (file) {

      int atoms = 0;

      atoms = hex_writePdbFile(file, image, chain_id, chain_label, atom_mode, 
                               with_sse, with_ter, atoms, t_export);

      hex_closePdbFile(file);
   }
}

/* ----------------------------------------------------------------------------- */
// open a new PDB file for writing

HexFile *hex_openPdbFile(char *filename)
{
   HexFile  *file;
   char      pdb_file[MAX_PATHNAME], ctmp[MAX_PATHNAME];

   strcpy(pdb_file, hex_filename(filename, "pdb", ctmp));

   if ((file=hex_fopen(pdb_file, "w"))) {

      if (hex_debug >= 2) hex_msg("Writing PDB: %s\n", pdb_file);

      hex_get_date(ctmp);

      hex_fprintf(file, "REMARK  File written by %s on %s.\n", hex_version(), ctmp);

   } else {

      hex_err("Failed to write PDB file: %s\n", pdb_file);
   }

   return(file);
}

/* ----------------------------------------------------------------------------- */

void hex_closePdbFile(HexFile *file)
{
   if (file) hex_fclose(file);
}

/* ----------------------------------------------------------------------------- */
// write out SSE commands to a PDB file ...
// (NB. even though this concerns residues, we count off atoms in order to be
//  able to apply the same chain re-numbering scheme as for writing atoms)

void hex_writeSSEInfo(HexFile *file, PdbImage *image, 
                      char *sse_type, int chain_id, char chain_label)
{
   int      a1, a2, r, r1, r2, s, s1, s2, s_num, c_num;
   char     chain1[] = " ";
   char     chain2[] = " ";
   char     s_label[8];

   s_num = 1;   // SSE serial number
   c_num = -99; // chain serial number

   r1 = 0;

   while (r1 < image->n_ca) {

      if (strcmp(image->residue[r1].sse, sse_type) == 0) {

         r2 = r1;

         for (s=r1+1; s<image->n_ca; s++) {

            if (strcmp(image->residue[s].sse,   image->residue[r1].sse) != 0) break;
            if (       image->residue[s].cid != image->residue[r1].cid)       break;
      
            r2 = s;
         }

// now have residues r1:r2 of the desired type on same chain -- check chain type and mask

         s1 = s2 = -1;
   
         for (r=r1; r<=r2; r++) { 
   
            if ((image->residue[r].cid == chain_id || chain_id == -1) && 
                 image->residue[r].mask == 0) {
   
               s1 = r;
               break;
            }
         }
   
         for (r=r2; r>=r1; r--) { // check chain type and mask
   
            if ((image->residue[r].cid == chain_id || chain_id == -1) && 
                 image->residue[r].mask == 0) {
   
               s2 = r;
               break;
            }
         }
   
// if we still have a valid range, write the first and last residue of this SSE
   
         if (s1 >= 0 && s2 >= 0) {

            a1 = image->residue[s1].atom1;
            a2 = image->residue[s2].atom1;
   
            PdbAtom *atom1 = &image->atom[a1];
            PdbAtom *atom2 = &image->atom[a2];
   
            set_chain(image, a1, chain_label, chain1);
            set_chain(image, a2, chain_label, chain2);
   
            if (image->residue[s1].cid == c_num) {

// advance SSE serial number on current chain (limit of 99 avoids 2-digit field overflow)

               s_num += 1;

               if (s_num > 99) s_num = 1;  

            } else { // reset serial number for new chain

               c_num = image->residue[s1].cid;
               s_num = 1;
            }

            if (strcmp(sse_type, "A") == 0) {  // HELIX
   
               sprintf(s_label, "H%2.2d", s_num);

               hex_fprintf(file, "%-6s %3d %-3s %3s %1s %4s%1s %3s %1s %4s%1s%2d%-30s %5d\n",
                           "HELIX ",
                           s_num,
                           s_label,
                           atom1->residue,
                           chain1,
                           atom1->sequence,
                           atom1->insert,
                           atom2->residue,
                           chain2,
                           atom2->sequence,
                           atom2->insert,
                           1,
                           " ",
                           s2 - s1 + 1);
            }

            if (strcmp(sse_type, "B") == 0) {  

// write a beta STRAND (don't bother trying to define a "SHEET")

               sprintf(s_label, "S%2.2d", s_num);

               hex_fprintf(file, "%-6s %3d %-3s%2d %3s %1s%4s%1s %3s %1s%4s%s%2d\n",
                           "SHEET ",
                           s_num,
                           s_label,
                           1,
                           atom1->residue,
                           chain1,
                           atom1->sequence,
                           atom1->insert,
                           atom2->residue,
                           chain2,
                           atom2->sequence,
                           atom2->insert,
                           0);
            }
   
            if (strcmp(sse_type, "C") == 0) {  // write a TURN

               sprintf(s_label, "T%2.2d", s_num);

               hex_fprintf(file, "%-6s %3d %-3s %3s %1s%4s%1s %3s %1s%4s%s\n",
                           "TURN  ",
                           s_num,
                           s_label,
                           atom1->residue,
                           chain1,
                           atom1->sequence,
                           atom1->insert,
                           atom2->residue,
                           chain2,
                           atom2->sequence,
                           atom2->insert);
            }
        
         }
   
// advance to next residue
   
         r1 = r2 + 1;

      } else {

         r1 = r1 + 1;
      }
   }
}

/* ----------------------------------------------------------------------------- */
// write out SSE records to a PDB file ...

void hex_writePdbSSEs(HexFile *file, PdbImage *image, int chain_id, char chain_label)
{
   hex_writeSSEInfo(file, image, "A", chain_id, chain_label);
   hex_writeSSEInfo(file, image, "B", chain_id, chain_label);
   hex_writeSSEInfo(file, image, "C", chain_id, chain_label);
}

/* ----------------------------------------------------------------------------- */
// write a PDB image to an already open file stream (initially set atom_start=0)

int hex_writePdbFile(HexFile *file, PdbImage *image, int chain_id, char chain_label,
                     int atom_mode, int with_sse, int with_ter, int atom_num, 
                     Matrix3D t_export)
{
   PdbAtom *atom = NULL;
   char    the_chain[] = " ";

   if (with_sse) { // write SSE records, if there is any SSE information

      hex_writePdbSSEs(file, image, chain_id, chain_label);
   }

// now write the coordinates...

   for (int a=0; a<image->n_atoms; a++) {

      int r = image->atom[a].rid;
/*
      if (strcmp(image->atom[a].atomName, "OXT") == 0) {

         hex_msg("OXT at a = %d\n", a);
      }
*/
      if (image->residue[r].mask != 0) continue;

      if (image->residue[r].cid != chain_id && chain_id != -1) continue;

      if (atom_mode == 4 && !hex_isMainChain  (image->atom[a].atomName)) continue;
      if (atom_mode == 3 && !hex_isSideChain  (image->atom[a].atomName)) continue;
      if (atom_mode == 2 && !hex_isAlphaCarbon(image->atom[a].atomName)) continue;
      if (atom_mode == 1 && !hex_isBackbone   (image->atom[a].atomName)) continue;


      atom = &image->atom[a];

      Point3D pt = pt_tfpt(t_export, atom->pt);

      atom_num++;

      if (atom_num > 99999) atom_num = 1;

      set_chain(image, a, chain_label, the_chain);

      hex_fprintf(file, "%-6s%5d %-4s%1s%-3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                  atom->record,
                  atom_num,
                  atom->eName,
                  atom->alt,
                  atom->residue,
                  the_chain,
                  atom->sequence,
                  atom->insert,
                  pt.x, pt.y, pt.z,
                  atom->occ,
                  atom->tfac);
   }

   if (with_ter && atom) {

      atom_num++;

      hex_fprintf(file, "%-6s%5d %-4s%1s%-3s %1s%4s%1s\n",
                  "TER", atom_num, "    ", " ", atom->residue, the_chain, 
                   atom->sequence, atom->insert);
   }

// append any hetero atoms (using their original chain labels)

   for (int h=0; h<image->n_hetero; h++) {

      atom = &image->hetero[h];

      Point3D pt = pt_tfpt(t_export, atom->pt);

      atom_num++;

      if (atom_num > 99999) atom_num = 1;

      hex_fprintf(file, "%-6s%5d %-4s%1s%-3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                  atom->record,
                  atom_num,
                  atom->eName,
                  atom->alt,
                  atom->residue,
                  atom->chain,
                  atom->sequence,
                  atom->insert,
                  pt.x, pt.y, pt.z,
                  atom->occ,
                  atom->tfac);
   }

   return(atom_num);
}

/* ----------------------------------------------------------------------------- */
// similarly, write one residue of a PDB image to an already open file stream;
// this is mainly for writing a flexibly fitted PDB from Kpax, where each
// residue might have a different transform
//
// NOT USED

int hex_writePdbResidue(HexFile *file, PdbImage *image, int rid, 
                        char chain_label, int atom_mode, 
                        int with_ter, int atom_num, Matrix3D t_export)
{
   PdbAtom *atom = NULL;
   char    the_chain[] = " ";

   if (rid >= 0 && rid < image->n_ca) {

      if (image->residue[rid].mask == 0) {

         int atom1 = image->residue[rid].atom1;
         int atom2 = image->residue[rid].atom2;

         for (int a=atom1; a<=atom2; a++) {

            if (atom_mode == 4 && !hex_isMainChain  (image->atom[a].atomName)) continue;
            if (atom_mode == 3 && !hex_isSideChain  (image->atom[a].atomName)) continue;
            if (atom_mode == 2 && !hex_isAlphaCarbon(image->atom[a].atomName)) continue;
            if (atom_mode == 1 && !hex_isBackbone   (image->atom[a].atomName)) continue;
      
            atom = &image->atom[a];
      
            Point3D pt = pt_tfpt(t_export, atom->pt);
      
            atom_num++;
      
            if (atom_num > 99999) atom_num = 1;
      
            set_chain(image, a, chain_label, the_chain);
      
            hex_fprintf(file, "%-6s%5d %-4s%1s%-3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                        atom->record,
                        atom_num,
                        atom->eName,
                        atom->alt,
                        atom->residue,
                        the_chain,
                        atom->sequence,
                        atom->insert,
                        pt.x, pt.y, pt.z,
                        atom->occ,
                        atom->tfac);
         }
      }
   }

   if (with_ter && atom && rid == image->n_ca-1) {

      atom_num++;

      hex_fprintf(file, "%-6s%5d %-4s%1s%-3s %1s%4s%1s\n",
                  "TER", atom_num, "    ", " ", atom->residue, the_chain, 
                   atom->sequence, atom->insert);
   }

   return(atom_num);
}

/* ----------------------------------------------------------------------------- */
// write a binary PDB file (do not write any hetero atoms here)

int hex_writeBinaryPdb(char *filename, PdbImage *image)
{
   int       n_bytes, status;
   HexFile  *file;
   int       ibuf[12];
   char      pdb_file[MAX_PATHNAME], ctmp[MAX_PATHNAME];

   strcpy(pdb_file, hex_filename(filename, "kpdb", ctmp));

   if ((file = hex_fopen(pdb_file, "wb"))) {

      if (hex_debug > 2) hex_msg("Writing PDB: %s\n", pdb_file);

// first, write a header of 12 integers

      ibuf[0] = image->n_models;
      ibuf[1] = image->n_chains;
      ibuf[2] = image->n_atoms;
      ibuf[3] = image->n_ca;
      ibuf[4] = image->n_bb;
      ibuf[5] = image->n_mc;
      ibuf[6] = image->n_sc;
      ibuf[7] = image->n_cryst;
      ibuf[8] = image->n_biomt;
   
      ibuf[9]  = image->n_hetero;  // New in version 4.0.0

      ibuf[10] = 0;  // spare
      ibuf[11] = 0;  // spare

      status = hex_fwrite(file, ibuf, 12*sizeof(int));

// write the transformation matrix
   
      if (status) {

         status = hex_fwrite(file, &image->mat_import, 1*sizeof(Matrix3D));
      }
   
// write all the atom and residue data
   
      if (status) {

         status = hex_fwrite(file, image->atom, image->n_atoms*sizeof(PdbAtom));
      }

      if (status) {

         status = hex_fwrite(file, image->residue, image->n_ca*sizeof(PdbResidue));
      }

// write the optional crystallographics transforms

      if (status && image->n_cryst > 0) {
   
         status = hex_fwrite(file, image->mat_cryst, image->n_cryst*sizeof(Matrix3D));
      }
     
      if (status && image->n_biomt > 0) {
   
         status = hex_fwrite(file, image->mat_biomt, image->n_biomt*sizeof(Matrix3D));
      }

// write the optional hetero atom coordinates [Added in version 4.0.0]

      if (status && image->n_hetero > 0) {

         status = hex_fwrite(file, image->hetero, image->n_hetero*sizeof(PdbAtom));
      }

      hex_fclose(file);

   } else {  // failed to open

      status = 0;
   }

   return(status);
}

/* ----------------------------------------------------------------------------- */
// read a binary PDB file into an in-memory image

PdbImage *hex_readBinaryPdb(char *filename, char *nickname)
{
   int       n_bytes, n_spare, status;
   HexFile  *file;
   PdbImage *image;
   int       ibuf[12];
   char      pdb_file[MAX_PATHNAME], ctmp[MAX_PATHNAME], rname[MAX_PATHNAME];

   image = (PdbImage *) kgetB_0(sizeof(PdbImage));

   if (hex_filesize(filename) >= 0) strcpy(pdb_file, filename);
   else                             strcpy(pdb_file, hex_filename(filename, 
                                                                  "kpdb",ctmp));
// first, set the name fields
   
   image->name = kdupS(pdb_file);
   
   if (nickname) image->nickname = kdupS(nickname);
   else          image->nickname = kdupS(hex_rootname(
                                         hex_basename(pdb_file, ctmp), rname));
   
   if ((file = hex_fopen(pdb_file, "rb"))) {

      if (hex_debug > 2) hex_msg("Reading PDB: %s\n", pdb_file);

// next, read a header of 12 integers

      status = hex_fread(file, ibuf, 12*sizeof(int));

      if (status) {

         image->n_models = ibuf[0];
         image->n_chains = ibuf[1];
         image->n_atoms  = ibuf[2];
         image->n_ca     = ibuf[3];
         image->n_bb     = ibuf[4];
         image->n_mc     = ibuf[5];
         image->n_sc     = ibuf[6];
         image->n_cryst  = ibuf[7];
         image->n_biomt  = ibuf[8];
   
         image->n_hetero = ibuf[9];  // added in version 4.0.0

         n_spare         = ibuf[10];
         n_spare         = ibuf[11];
      }

// read the transformation matrix
   
      if (status) {

         status = hex_fread(file, &image->mat_import, 1*sizeof(Matrix3D));
      }
   
// read all the atom and residue data
   
      image->atom      = (PdbAtom *)    kgetB(image->n_atoms *sizeof(PdbAtom));
      image->residue   = (PdbResidue *) kgetB(image->n_ca    *sizeof(PdbResidue));
   
      if (status) {

         status = hex_fread(file, image->atom, image->n_atoms*sizeof(PdbAtom));
      }

      if (status) {

         status = hex_fread(file, image->residue, image->n_ca*sizeof(PdbResidue));
      }

      image->atom_backup   = (Point3D *) kgetB(image->n_atoms *sizeof(Point3D));
      image->hetero_backup = (Point3D *) kgetB(image->n_hetero*sizeof(Point3D));

      for (int a=0; a<image->n_atoms; a++) {

         image->atom_backup[a] = image->atom[a].pt;
      }

      for (int h=0; h<image->n_hetero; h++) {

         image->hetero_backup[h] = image->hetero[h].pt;
      }

// read the optional crystallographics transforms

      if (status && image->n_cryst > 0) {
   
         image->mat_cryst = (Matrix3D *) kgetB(image->n_cryst*sizeof(Matrix3D));

         status = hex_fread(file, image->mat_cryst, image->n_cryst*sizeof(Matrix3D));
      }
     
      if (status && image->n_biomt > 0) {
   
         image->mat_biomt = (Matrix3D *) kgetB(image->n_biomt*sizeof(Matrix3D));

         status = hex_fread(file, image->mat_biomt, image->n_biomt*sizeof(Matrix3D));
      }

      if (status && image->n_hetero > 0) {

         image->hetero = (PdbAtom *) kgetB(image->n_hetero*sizeof(PdbAtom));

         status = hex_fread(file, image->hetero, image->n_hetero*sizeof(PdbAtom));
      }

      hex_fclose(file);

      if (status == 0) {  // failed to read

         hex_freePdb(image);

         image = NULL;
      }

   } else {  // failed to open

      hex_freePdb(image);

      image = NULL;
   }

   return(image);
}

/* ----------------------------------------------------------------------------- */
// get the list of atom numbers of the non-masked CA positions by counting off
// atoms in exactly the same was as when actually writing the file (as above).
// This is mainly to allow Kpax to be able to generate CONECT records.

int hex_getPdbAlphaNos(PdbImage *image, int chain_id, int atom_mode, int atom_num, 
                     int *a_num)
{
   int na = 0;

   for (int a=0; a<image->n_atoms; a++) {

      int r = image->atom[a].rid;

      if (image->residue[r].mask != 0) continue;

      if (image->residue[r].cid != chain_id && chain_id != -1) continue;

      if (atom_mode == 4 && !hex_isMainChain  (image->atom[a].atomName)) continue;
      if (atom_mode == 3 && !hex_isSideChain  (image->atom[a].atomName)) continue;
      if (atom_mode == 2 && !hex_isAlphaCarbon(image->atom[a].atomName)) continue;
      if (atom_mode == 1 && !hex_isBackbone   (image->atom[a].atomName)) continue;

      atom_num++;

      if (atom_num > 99999) atom_num = 1;

      if (hex_isAlphaCarbon(image->atom[a].atomName)) {
 
         a_num[na++] = atom_num;
      }
   }

   return(na);
}

/* ----------------------------------------------------------------------------- */
// OBSOLETE - export a PDB from a database back to user space 

// this is obsolete, and should not be used because it is not compatible with
// binary PDB files. Instead, read "from" as binary and write "to" as formatted.

void  hex_exportPdb(char *from_file, char *to_file, int atom_mode, Matrix3D t_export)
{
   PdbAtom  atom;
   HexFile *fin, *fout;
   char     line[MAX_PATHNAME];
   char     src_pdb[MAX_PATHNAME], dst_pdb[MAX_PATHNAME], ctmp[MAX_PATHNAME];

// force ".pdb" onto the given file names

   strcpy(src_pdb, hex_filename(from_file, "pdb", ctmp));  
   strcpy(dst_pdb, hex_filename(to_file,   "pdb", ctmp));

   if ((fin=hex_fopen(src_pdb, "r"))) {

      if ((fout=hex_fopen(dst_pdb, "w"))) {

         if (hex_freadln(fin, 80, line) != HEX_FILE_EOF) {

            hex_get_date(ctmp);

            hex_fprintf(fout, "REMARK  File written by %s on %s.\n", 
                        hex_version(), ctmp);

            for (int a=0;;) {

               if (strncmp(line, "ATOM", 4) == 0) {

                  hex_parsePdbAtom(line, 0, &atom);

                  if (atom_mode == 4 && !hex_isMainChain  (atom.atomName)) continue;
                  if (atom_mode == 3 && !hex_isSideChain  (atom.atomName)) continue;
                  if (atom_mode == 2 && !hex_isAlphaCarbon(atom.atomName)) continue;
                  if (atom_mode == 1 && !hex_isBackbone   (atom.atomName)) continue;

                  a++;

                  atom.pt = pt_tfpt(t_export, atom.pt);

                  hex_fprintf(fout, "%-6s%5d %-4s%1s%-3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                              atom.record,
                              a,
                              atom.eName,
                              atom.alt,
                              atom.residue,
                              atom.chain,
                              atom.sequence,
                              atom.insert,
                              atom.pt.x, 
                              atom.pt.y, 
                              atom.pt.z,
                              atom.occ,
                              atom.tfac);
               }

               if (hex_freadln(fin, 80, line)== HEX_FILE_EOF) break;
            }





         } else {

            hex_err("Empty PDB file: %s\n", src_pdb);
         }

         hex_fclose(fout);

      } else {                       // could happen

         hex_err("Failed to write PDB file: %s\n", dst_pdb);
      }

      hex_fclose(fin);

   } else {                       // should never happen

      hex_err("Failed to read PDB file: %s\n", src_pdb);
   }
}

/* ----------------------------------------------------------------------------- */
// extract atom coordinates

Point3D *hex_ptsPdb(PdbImage *image, int atom_mode)
{
   int       n_pts;
   Point3D  *pts;

   if      (atom_mode == 4) n_pts = image->n_mc;
   else if (atom_mode == 3) n_pts = image->n_sc;
   else if (atom_mode == 2) n_pts = image->n_ca;
   else if (atom_mode == 1) n_pts = image->n_bb;
   else                     n_pts = image->n_atoms;

   pts = (Point3D *) kgetB_0(n_pts*sizeof(Point3D));

   for (int a=0, j=0; a<image->n_atoms; a++) {

      if (atom_mode == 4 && !hex_isMainChain  (image->atom[a].atomName)) continue;
      if (atom_mode == 3 && !hex_isSideChain  (image->atom[a].atomName)) continue;
      if (atom_mode == 2 && !hex_isAlphaCarbon(image->atom[a].atomName)) continue;
      if (atom_mode == 1 && !hex_isBackbone   (image->atom[a].atomName)) continue;

      pts[j++] = image->atom[a].pt;

      if (j >= n_pts) break;
   }

   return(pts);
}

/* ----------------------------------------------------------------------------- */

PdbPeptide *hex_allocPdbPep(int np)
{
   PdbPeptide *pep = (PdbPeptide *) kgetB(sizeof(PdbPeptide));

   pep->ca = (Point3D *) kgetB_0(np*sizeof(Point3D));
   pep->n  = (Point3D *) kgetB_0(np*sizeof(Point3D));
   pep->c  = (Point3D *) kgetB_0(np*sizeof(Point3D));
   pep->o  = (Point3D *) kgetB_0(np*sizeof(Point3D));
   pep->cb = (Point3D *) kgetB_0(np*sizeof(Point3D));
// pep->sm = (Point3D *) kgetB_0(np*sizeof(Point3D));
   pep->aa_code = (char *) kgetB_0(np);
   pep->aa_flag = (char *) kgetB_0(np);
   pep->wt = (double  *) kget1D_0(np);

   return(pep);
}

/* ----------------------------------------------------------------------------- */

void hex_freePdbPep(PdbPeptide *pep)
{
   if (pep) {

      kfree(pep->ca);
      kfree(pep->n);
      kfree(pep->c);
      kfree(pep->o);
      kfree(pep->cb);
//    kfree(pep->sm);

      kfree(pep->aa_code);
      kfree(pep->aa_flag);

      kfree(pep->wt);

      kfree(pep);
   }
}

/*---------------------------------------------------------------------------*/
// copy a given number of peptide coordinates from "src" to "dst"

void hex_copyPdbPep(PdbPeptide *src, int start, PdbPeptide *dst, int np)
{
   hex_copy(&src->ca[start], dst->ca, np*sizeof(Point3D));
   hex_copy(&src->n[start],  dst->n,  np*sizeof(Point3D));
   hex_copy(&src->c[start],  dst->c,  np*sizeof(Point3D));
   hex_copy(&src->o[start],  dst->o,  np*sizeof(Point3D));
   hex_copy(&src->cb[start], dst->cb, np*sizeof(Point3D));
// hex_copy(&src->sm[start], dst->sm, np*sizeof(Point3D));
   hex_copy(&src->wt[start], dst->wt, np*sizeof(double));

   hex_copy(&src->aa_code[start], dst->aa_code, np*sizeof(char));
   hex_copy(&src->aa_flag[start], dst->aa_flag, np*sizeof(char));
}

/*---------------------------------------------------------------------------*/
// duplicate an entire PdpPeptide structure and all its contents

PdbPeptide *hex_dupPdbPep(int np, PdbPeptide*src)
{
   PdbPeptide *pep = (PdbPeptide *) kgetB(sizeof(PdbPeptide));

   pep->ca = (Point3D *) kgetB(np*sizeof(Point3D));
   pep->n  = (Point3D *) kgetB(np*sizeof(Point3D));
   pep->c  = (Point3D *) kgetB(np*sizeof(Point3D));
   pep->o  = (Point3D *) kgetB(np*sizeof(Point3D));
   pep->cb = (Point3D *) kgetB(np*sizeof(Point3D));
// pep->sm = (Point3D *) kgetB(np*sizeof(Point3D));
   pep->aa_code = (char *) kgetB(np*sizeof(char));
   pep->aa_flag = (char *) kgetB(np*sizeof(char));
   pep->wt = (double  *) kget1D(np);

   hex_copy(src->ca, pep->ca, np*sizeof(Point3D));
   hex_copy(src->n,  pep->n,  np*sizeof(Point3D));
   hex_copy(src->c,  pep->c,  np*sizeof(Point3D));
   hex_copy(src->o,  pep->o,  np*sizeof(Point3D));
   hex_copy(src->cb, pep->cb, np*sizeof(Point3D));
// hex_copy(src->sm, pep->sm, np*sizeof(Point3D));

   hex_copy(src->aa_code, pep->aa_code, np*sizeof(char));
   hex_copy(src->aa_flag, pep->aa_flag, np*sizeof(char));

   hex_copy(src->wt, pep->wt, np*sizeof(double));

   return(pep);
}

/* ----------------------------------------------------------------------------- */

void hex_loadPdbPep(PdbImage *image, PdbPeptide *pep)
{
   int         a, r, nn, na, nc, no, np, is_mc;

   np = image->n_ca;

// collect up the side chain COM coordinates

   for (r=0; r<np; r++) {

       pep->aa_code[r] = image->residue[r].code[0];
       pep->aa_flag[r] = image->residue[r].sse[0];
   }

// collect up the remaining main chain and CB atom coordinates

// nn = na = nc = no = 0;

   for (a=0; a<image->n_atoms; a++) {

      is_mc = hex_isMainChain(image->atom[a].atomName);

      r = image->atom[a].rid;

      if (is_mc == PDB_N && r < np) {

         pep->n[r] = image->atom[a].pt;
//       nn++;

      } else if (is_mc == PDB_A && r < np) {

         pep->ca[r] = image->atom[a].pt;
//       na++;

      } else if (is_mc == PDB_C && r < np) {

         pep->c[r] = image->atom[a].pt;
//       nc++;

      } else if (is_mc == PDB_O && r < np) {

         pep->o[r] = image->atom[a].pt;
//       no++;
      
      } else if (is_mc == PDB_B && r < np) {

         pep->cb[r] = image->atom[a].pt;
      }
   }
/*
   if (nn != np) hex_wrn("Some N  atoms are missing in PDB file %s\n", image->name);
   if (na != np) hex_wrn("Some CA atoms are missing in PDB file %s\n", image->name);
   if (nc != np) hex_wrn("Some C  atoms are missing in PDB file %s\n", image->name);
   if (no != np) hex_wrn("Some O  atoms are missing in PDB file %s\n", image->name);
*/
}

/* ----------------------------------------------------------------------------- */
// supply the main chain atoms in a standard format

PdbPeptide *hex_makePdbPep(PdbImage *image)
{
   PdbPeptide *pep = (PdbPeptide *) hex_allocPdbPep(image->n_ca);

   hex_loadPdbPep(image, pep);

   return(pep);
}

/* ----------------------------------------------------------------------------- */
// extract atom radii

float *hex_radPdb(PdbImage *image, int atom_mode)
{
   int      a, j, n_pts;
   char     type, check[2];
   byte     pr, lp;
   sint     opls, ace;
   float    charge;
   float   *rad;
   PdbAtom *atom;

   if      (atom_mode == 4) n_pts = image->n_mc;
   else if (atom_mode == 3) n_pts = image->n_sc;
   else if (atom_mode == 2) n_pts = image->n_ca;
   else if (atom_mode == 1) n_pts = image->n_bb;
   else                     n_pts = image->n_atoms;

   rad = (float *) kgetB_0(n_pts*sizeof(float));

   for (a=0, j=0; a<image->n_atoms; a++) {

      atom = &image->atom[a];

      if (atom_mode == 4 && !hex_isMainChain  (atom->atomName)) continue;
      if (atom_mode == 3 && !hex_isSideChain  (atom->atomName)) continue;
      if (atom_mode == 2 && !hex_isAlphaCarbon(atom->atomName)) continue;
      if (atom_mode == 1 && !hex_isBackbone   (atom->atomName)) continue;

      if (!hex_atom_data(atom->residue, atom->atomName, &rad[j],
                         &charge, &type, &opls, &ace, &pr, &lp)) {

         if (!hex_atom_data((char *) "UNK", atom->atomName, &rad[j],
                            &charge, &type, &opls, &ace, &pr, &lp)) {

            check[0] = atom->atomName[0];
            check[1] = '\0';

            if (!hex_atom_data((char *) "UNK", check, &rad[j],
                               &charge, &type, &opls, &ace, &pr, &lp)) {

               hex_atom_data((char *) "elm", check, &rad[j],
                             &charge, &type, &opls, &ace, &pr, &lp);
            }
         }
      }
      j++;

      if (j >= n_pts) break;
   }

   return(rad);
}


/* ---------------------------------------------------------------------------*/

int hex_isAA(char *res)
{
   if (hex_aaCode(res) != 'X') return(1);  // definitely a standard amino acid
   
   return(0);
}

/* ---------------------------------------------------------------------------*/

int hex_isDNA(char *res)
{
   if (hex_dnaCode(res) != 'X') return(1);  // definitely a standard DNA/RNA base 
   
   return(0);
}

/* ---------------------------------------------------------------------------*/

int hex_isAlphaCarbon(char *atomName)
{
   if (strcmp(atomName, "CA") == 0) return(PDB_A);

   return(0); 
}

/* ----------------------------------------------------------------------------*/

int hex_isMainChain(char *atomName)
{
   if (strcmp(atomName, "CA") == 0) return(PDB_A);
   if (strcmp(atomName, "N")  == 0) return(PDB_N);
   if (strcmp(atomName, "C")  == 0) return(PDB_C);
   if (strcmp(atomName, "O")  == 0) return(PDB_O);

   if (strcmp(atomName, "CB") == 0) return(PDB_B);

   if (strcmp(atomName, "OT")  == 0) return(PDB_OXT);
   if (strcmp(atomName, "OXT") == 0) return(PDB_OXT);

   if (strcmp(atomName, "OT1") == 0) return(PDB_OT1);
   if (strcmp(atomName, "OT2") == 0) return(PDB_OT2);

   return(0);
}

/* ----------------------------------------------------------------------------*/

int hex_isBackbone(char *atomName)
{
   if (strcmp(atomName, "CA") == 0) return(PDB_A);
   if (strcmp(atomName, "N")  == 0) return(PDB_N);
   if (strcmp(atomName, "C")  == 0) return(PDB_C);
   if (strcmp(atomName, "O")  == 0) return(PDB_O);

   if (strcmp(atomName, "OT")  == 0) return(PDB_OXT);
   if (strcmp(atomName, "OXT") == 0) return(PDB_OXT);

   if (strcmp(atomName, "OT1") == 0) return(PDB_OT1);
   if (strcmp(atomName, "OT2") == 0) return(PDB_OT2);

   return(0);
}

/* ---------------------------------------------------------------------------*/

int hex_isHydrogen(char *atomName)
{
   if (strcmp(atomName, "H") == 0) return(1);

   if (strlen(atomName) >= 2) {

      if (isdigit(atomName[0]) && atomName[1] == 'H') return(1);
      if (isdigit(atomName[1]) && atomName[0] == 'H') return(1);

      if (strncmp(atomName, "HA", 2)  == 0) return(1);
      if (strncmp(atomName, "HB", 2)  == 0) return(1);
      if (strncmp(atomName, "HG", 2)  == 0) return(1);
      if (strncmp(atomName, "HD", 2)  == 0) return(1);
      if (strncmp(atomName, "HE", 2)  == 0) return(1);
      if (strncmp(atomName, "HH", 2)  == 0) return(1);
      if (strncmp(atomName, "HZ", 2)  == 0) return(1);
   }

   return(0);
}

/* ---------------------------------------------------------------------------*/

int hex_isSideChain(char *atomName)
{
   if (hex_isBackbone(atomName)) return(0);
   if (hex_isHydrogen(atomName)) return(0);

   return(1);
}

/* ---------------------------------------------------------------------------*/

char *hex_aaName(char letter)
{
   if (letter == 'A') return("ALA"); 
   if (letter == 'R') return("ARG");
   if (letter == 'N') return("ASN");
   if (letter == 'D') return("ASP");
   if (letter == 'C') return("CYS");
   if (letter == 'E') return("GLU");
   if (letter == 'Q') return("GLN");
   if (letter == 'G') return("GLY");
   if (letter == 'H') return("HIS");
   if (letter == 'I') return("ILE");
   if (letter == 'L') return("LEU");
   if (letter == 'K') return("LYS");
   if (letter == 'M') return("MET");
   if (letter == 'F') return("PHE");
   if (letter == 'P') return("PRO");
   if (letter == 'S') return("SER");
   if (letter == 'T') return("THR");
   if (letter == 'W') return("TRP");
   if (letter == 'Y') return("TYR");
   if (letter == 'V') return("VAL");

   if (letter == 'B') return("ASX");  // ASN or ASP
   if (letter == 'Z') return("GLX");  // GLN or GLU
   if (letter == 'J') return("XLE");  // LEU or ILE
   if (letter == 'U') return("SEC");  // selenocysteine
   if (letter == 'O') return("PYL");  // pyrrolysine
   if (letter == 'X') return("UNK");  // Unknown

   return("UNK");
}

/* ---------------------------------------------------------------------------*/

char hex_aaCode(char *res)
{
   if (strcmp(res, "ALA") == 0) return('A');
   if (strcmp(res, "DAL") == 0) return('A');  // D-alanine
   if (strcmp(res, "ABA") == 0) return('A');  // alpha-amino-butyric acid
   if (strcmp(res, "AIB") == 0) return('A');  // 
   if (strcmp(res, "DAL") == 0) return('A');  // 
   if (strcmp(res, "ORN") == 0) return('A');  // ornithine

   if (strcmp(res, "ARG") == 0) return('R');
   if (strcmp(res, "DAR") == 0) return('R');  // D-arginine

   if (strcmp(res, "ASN") == 0) return('N');
   if (strcmp(res, "DSG") == 0) return('N');  // D-asparaginine

   if (strcmp(res, "ASP") == 0) return('D');
   if (strcmp(res, "DSP") == 0) return('D');  // D-aspartate
   if (strcmp(res, "ASX") == 0) return('D');  // 

   if (strcmp(res, "CYS") == 0) return('C');
   if (strcmp(res, "DCY") == 0) return('C');  // D-cysteine
   if (strcmp(res, "NPH") == 0) return('C');  // cysteine methylene-carbmoyl-phenanthroline
   if (strcmp(res, "CEA") == 0) return('C');  // s-hydrocy-cysteine
   if (strcmp(res, "CSD") == 0) return('C');  // sulfino-alanine
   if (strcmp(res, "CSW") == 0) return('C');  // cysteine-s-dioxide
   if (strcmp(res, "CSS") == 0) return('C');  // s-mercapto-cysteine
   if (strcmp(res, "CSO") == 0) return('C');  // s-hydroxy cysteine 
   if (strcmp(res, "CME") == 0) return('C');  // ss-2-hydroxy-thio-cysteine
   if (strcmp(res, "CSE") == 0) return('C');  // seleno-cysteine
   if (strcmp(res, "CSX") == 0) return('C');  // S-oxy cysteine
   if (strcmp(res, "CYG") == 0) return('C');  // 
   if (strcmp(res, "OCS") == 0) return('C');  // cysteine sulphonic acid
   if (strcmp(res, "SMC") == 0) return('C');  // s-methyl-cysteine 
   if (strcmp(res, "CSP") == 0) return('C');  // phosphorylated cysteine 
   if (strcmp(res, "CSR") == 0) return('C');  // 
   if (strcmp(res, "143") == 0) return('C');  // hack for 1jvn (CYS covalently linked to activin)

   if (strcmp(res, "SEC") == 0) return('U');  // seleno-cysteine

   if (strcmp(res, "GLU") == 0) return('E');
   if (strcmp(res, "DGL") == 0) return('E');  // D-gultamate
   if (strcmp(res, "PCA") == 0) return('E');  // pyro-glutamic acid
   if (strcmp(res, "CGU") == 0) return('E');  // gamma-carboxy-glutamic acid
   if (strcmp(res, "5HP") == 0) return('E');  // 5-hydroxyproline
   if (strcmp(res, "GLX") == 0) return('E');  // 

   if (strcmp(res, "GLN") == 0) return('Q');
   if (strcmp(res, "DGN") == 0) return('Q');  // D-gultamine

   if (strcmp(res, "GLY") == 0) return('G');
   if (strcmp(res, "SAR") == 0) return('G');  // sarcosine

   if (strcmp(res, "HIS") == 0) return('H');
   if (strcmp(res, "DHI") == 0) return('H');  // D-histidine
   if (strcmp(res, "HSE") == 0) return('H');
   if (strcmp(res, "HSP") == 0) return('H');
   if (strcmp(res, "HSD") == 0) return('H');
   if (strcmp(res, "HIP") == 0) return('H');
   if (strcmp(res, "HID") == 0) return('H');
   if (strcmp(res, "HIE") == 0) return('H');

   if (strcmp(res, "HSM") == 0) return('H');  // histamine

   if (strcmp(res, "ILE") == 0) return('I');
   if (strcmp(res, "DIL") == 0) return('I');  // D-isoleucine

   if (strcmp(res, "LEU") == 0) return('L');
   if (strcmp(res, "DLE") == 0) return('L');  // D-leucine
   if (strcmp(res, "XLE") == 0) return('L');  // nor-leucine
   if (strcmp(res, "MLE") == 0) return('L');  // n-methyl-leucine
   if (strcmp(res, "NLE") == 0) return('L');  // nor-leucine

   if (strcmp(res, "LYS") == 0) return('K');
   if (strcmp(res, "MLY") == 0) return('K');  // di-methyl Lysine
   if (strcmp(res, "KCX") == 0) return('K');  // 
   if (strcmp(res, "LLP") == 0) return('K');  // 
   if (strcmp(res, "PYL") == 0) return('K');  // pyrrolysine

   if (strcmp(res, "MET") == 0) return('M');
   if (strcmp(res, "CXM") == 0) return('M');  // n-carboxy-methionine
   if (strcmp(res, "MSE") == 0) return('M');  // seleno-methionine
   if (strcmp(res, "FME") == 0) return('M');  // formyl-methionine
   if (strcmp(res, "SAM") == 0) return('M');  // adeno-methyl-methionine
   if (strcmp(res, "OMT") == 0) return('M');  // s-dioxy-methionine


   if (strcmp(res, "PRO") == 0) return('P');
   if (strcmp(res, "DPR") == 0) return('P');  // D-proline
   if (strcmp(res, "HYP") == 0) return('P');

   if (strcmp(res, "PHE") == 0) return('F');
   if (strcmp(res, "DPN") == 0) return('F');  // D-phenylalanine
   if (strcmp(res, "TPQ") == 0) return('F');  // topo-quinine

   if (strcmp(res, "SER") == 0) return('S');
   if (strcmp(res, "DSN") == 0) return('S');  // D-serine
   if (strcmp(res, "SEP") == 0) return('S');  // phospho-serine
   if (strcmp(res, "SEB") == 0) return('S');  // o-benzyl-sulfonyl-serine

   if (strcmp(res, "THR") == 0) return('T');
   if (strcmp(res, "DTH") == 0) return('T');  // D-threonine
   if (strcmp(res, "TPO") == 0) return('T');  // phospho-threonine
   if (strcmp(res, "BMT") == 0) return('T');  // 

   if (strcmp(res, "TRP") == 0) return('W');
   if (strcmp(res, "DTR") == 0) return('W');  // D-tryptophan

   if (strcmp(res, "TYR") == 0) return('Y');
   if (strcmp(res, "DTY") == 0) return('Y');  // D-typrosine
   if (strcmp(res, "PTR") == 0) return('Y');  // o-phospho-tyrosine
   if (strcmp(res, "STY") == 0) return('Y');  // tyrosine-o-sulphonic acid
   if (strcmp(res, "TYS") == 0) return('Y');  // sulphonated-tyrosine

   if (strcmp(res, "VAL") == 0) return('V');
   if (strcmp(res, "DVA") == 0) return('V');  // D-valine
   if (strcmp(res, "DIV") == 0) return('V');  // D-isovaline
   if (strcmp(res, "MVA") == 0) return('V');  // n-methyl-valine

   return('X');
}

/* ---------------------------------------------------------------------------*/

char hex_dnaCode(char *res)
{
   if (strcmp(res, "A") == 0)   return('A');
   if (strcmp(res, "C") == 0)   return('C');
   if (strcmp(res, "I") == 0)   return('I');
   if (strcmp(res, "G") == 0)   return('G');
   if (strcmp(res, "T") == 0)   return('T');
   if (strcmp(res, "U") == 0)   return('U');

   if (strcmp(res, "+A") == 0)   return('A');
   if (strcmp(res, "+C") == 0)   return('C');
   if (strcmp(res, "+I") == 0)   return('I');
   if (strcmp(res, "+G") == 0)   return('G');
   if (strcmp(res, "+T") == 0)   return('T');
   if (strcmp(res, "+U") == 0)   return('U');

   if (strcmp(res, "DA") == 0)   return('A');
   if (strcmp(res, "DC") == 0)   return('C');
   if (strcmp(res, "DI") == 0)   return('I');
   if (strcmp(res, "DG") == 0)   return('G');
   if (strcmp(res, "DT") == 0)   return('T');
   if (strcmp(res, "DU") == 0)   return('U');

   return('X');
}

/* ----------------------------------------------------------------------------*/

int hex_isSugar(char *atomName)
{
   if (strcmp(atomName, "C1\'") == 0) return(1);
   if (strcmp(atomName, "C2\'") == 0) return(2);
   if (strcmp(atomName, "C3\'") == 0) return(3);
   if (strcmp(atomName, "C4\'") == 0) return(4);
   if (strcmp(atomName, "C5\'") == 0) return(5);

   if (strcmp(atomName, "C1*")  == 0) return(1);
   if (strcmp(atomName, "C2*")  == 0) return(2);
   if (strcmp(atomName, "C3*")  == 0) return(3);
   if (strcmp(atomName, "C4*")  == 0) return(4);
   if (strcmp(atomName, "C5*")  == 0) return(5);

   if (strcmp(atomName, "O3\'") == 0) return(6);
   if (strcmp(atomName, "O4\'") == 0) return(7);
   if (strcmp(atomName, "O5\'") == 0) return(8);

   if (strcmp(atomName, "O3*")  == 0) return(6);
   if (strcmp(atomName, "O4*")  == 0) return(7);
   if (strcmp(atomName, "O5*")  == 0) return(8);

   if (strcmp(atomName, "P")    == 0) return(9);
   if (strcmp(atomName, "OP1")  == 0) return(10);
   if (strcmp(atomName, "O1P")  == 0) return(10);
   if (strcmp(atomName, "OP2")  == 0) return(11);
   if (strcmp(atomName, "O2P")  == 0) return(11);

   return(0);
}

/* ----------------------------------------------------------------------------*/

int hex_isBase(char *atomName)
{
   if (strcmp(atomName, "N1")  == 0) return(1);
   if (strcmp(atomName, "N2")  == 0) return(1);
   if (strcmp(atomName, "N3")  == 0) return(1);
   if (strcmp(atomName, "N4")  == 0) return(1);
   if (strcmp(atomName, "N6")  == 0) return(1);
   if (strcmp(atomName, "N7")  == 0) return(1);
   if (strcmp(atomName, "N9")  == 0) return(1);

   if (strcmp(atomName, "O2")  == 0) return(1);
   if (strcmp(atomName, "O4")  == 0) return(1);

   if (strcmp(atomName, "C2")  == 0) return(1);
   if (strcmp(atomName, "C4")  == 0) return(1);
   if (strcmp(atomName, "C5")  == 0) return(1);
   if (strcmp(atomName, "C5M") == 0) return(1);
   if (strcmp(atomName, "C6")  == 0) return(1);
   if (strcmp(atomName, "C8")  == 0) return(1);

   return(0);
}

/* ---------------------------------------------------------------------------*/
