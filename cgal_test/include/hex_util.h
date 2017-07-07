/*-----------------------------------------------------------------------------
**  File:       hex_util.h
**
**  Author:     Dave Ritchie, 10/05/03
**
**  Purpose:    Header file for the "util-level" functions.
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

#ifndef hex_util_h
#define hex_util_h

#include "hex_sys.h"
#include "hex_transform.h"

#define min_(a,b)     (((a) < (b)) ? (a) : (b))
#define max_(a,b)     (((a) > (b)) ? (a) : (b))
#define abs_(a)       (((a) > 0) ? (a) : -(a))
#define square(x)    ((x)*(x))

#ifdef __cplusplus
extern "C" {
#endif

#define TBF_MAX_LINES       5000

#define SDF_MAX_FILENAME    256
#define SDF_LINE_LENGTH     133

#define SDF_SURFACE_TAG     "<SPHERICAL_HARMONIC_SURFACE>"
#define SDF_MEP_TAG         "<SPHERICAL_HARMONIC_MEP>"
#define SDF_IEL_TAG         "<SPHERICAL_HARMONIC_IEL>"
#define SDF_EAL_TAG         "<SPHERICAL_HARMONIC_EAL>"
#define SDF_EA_TAG          "<SPHERICAL_HARMONIC_EA>" /* same as EAL */
#define SDF_ALPHA_TAG       "<SPHERICAL_HARMONIC_ALPHA(L)>"

#define SDF_COG_TAG         "<MOLECULAR_CENTERS>"
#define SDF_DIPOL_TAG       "<DIPOL>"
#define SDF_DIPOLPHS_TAG    "<DIPOL_PHS>"
#define SDF_QUADPOL_TAG     "<QUADPOL>"
#define SDF_OCTUPOL_TAG     "<OCTUPOL>"
#define SDF_MULTIPOLE_TAG   "<ATOMIC MULTIPOLES>"
#define SDF_NAOPC_TAG       "<NAO-PC>"
#define SDF_VAMPBASICS_TAG  "<VAMPBASICS>"
#define SDF_MOPACBASICS_TAG "<MOPACBASICS>"
#define SDF_HAMILTONIAN_TAG "<HAMILTONIAN>"
#define SDF_DMATRIX_TAG     "<DENSITY MATRIX ELEMENTS>"

#define SDF_CASRN_TAG       "<CAS_RN>"
#define SDF_CAS_TAG         "<CAS>"
#define SDF_MOLNAME_TAG     "<MOLNAME>"

#define SDF_NAME_TAG        "<NAME>"

#define SDF_PG_RAD_TAG      "<ParaGraph RAD CRITICAL POINTS>"
#define SDF_PG_MEP_TAG      "<ParaGraph MEP CRITICAL POINTS>"
#define SDF_PG_IEL_TAG      "<ParaGraph IEL CRITICAL POINTS>"
#define SDF_PG_EAL_TAG      "<ParaGraph EAL CRITICAL POINTS>"
#define SDF_PG_POL_TAG      "<ParaGraph POL CRITICAL POINTS>"

#define TDS_DATE_TAG        "<3D_SNAP_DATE>"
#define TDS_ORIGIN_TAG      "<3D_SNAP_ORIGIN>"
#define TDS_MATRIX_TAG      "<3D_SNAP_MATRIX>"
#define TDS_SCALE_TAG       "<3D_SNAP_SCALE>"
#define TDS_ORDER_TAG       "<3D_SNAP_ORDER>"
#define TDS_OFFSET_TAG      "<3D_SNAP_OFFSET>"

#define HEX_YLM_NLOCAL        7

#define HEX_YLM_SURFACE       0  /* these are from ParaSurf */
#define HEX_YLM_MEP           1
#define HEX_YLM_IEL           2
#define HEX_YLM_EAL           3
#define HEX_YLM_ALPHA         4

#define HEX_YLM_R             5  /* these are from Hex */
#define HEX_YLM_Q             6

/*---------------------------------------------------------------------------*/

typedef struct _TextBuf {

   int    n_lines;    /* no. lines currently loaded */
   int    n_buffer;   /* current actual buffer size */
   char **buffer;     /* array of strings from hex_vm_strdup */

} TextBuf;


typedef struct _shd {  /* spherical harmonic data block, one per SDF file */

   int     l_min;
   int     l_max[HEX_YLM_NLOCAL];
   double *alm[HEX_YLM_NLOCAL];

} HexShd;


typedef struct _SdfAtm {

   int      n_atoms;
   int      n_bonds;

   Point3D *pt;
   int     *elem;
   double  *mass;
   double  *rvdw;
   char    *name;
   int     *bond;

} SdfAtm;

typedef struct _SdfFile {

   HexFile *hfile;
   int      status;    /* current file I/O status */
   int      asd;       /* true if its an anonymous SD file (".asd") */

} SdfFile;

typedef struct _org { /* molecular origin data block, one per SDF file */

   Point3D cog;       /* coordinate origin - "CoG" */
   Point3D coh;       /* spherical harmonic origin, usually = CoG */

} SdfCog;


typedef struct _TdsInf {

  char           *date;   /* date of the file creation, also used for authentication */
  Point3D        cog;     /* coordinate + sh origin */
  Transformation t;       /* transformation applied to the 3D coordinates but not the coeffs*/
  double         scale;   /* scaling factor */
  int            order;   /* maximum l */
  int            offset;  /* position of the coeffs in the binary file */
  double         *sig;
  double         *tau;

} TdsInf;

typedef struct _SdfMol {

   int     n_shd;     /* no. spherical harmonic properties actually read */
   int     n_lines;   /* no. lines for this molecule actually read */
   int     seq_no;    /* sequence number in SDF file, counting from 1 */
   int     mw;        /* calculated molecular weight */
   char   *cas_rn;
   char   *name;
   SdfCog *scg;
   SdfAtm *sda;
   HexShd *shd;

   TdsInf *tds;

} SdfMol;

typedef int     Node;
typedef int     Edge;
typedef int     Pole;

typedef struct _image {
   int width;
   int height;
   byte *buf;
} HexImage;


/*---------------------------------------------------------------------------*/

TextBuf *tbf_open       (int allocation);
void     tbf_close      (TextBuf *tb);
void     tbf_writeln    (TextBuf *tb, char *cline);
void     tbf_writev     (TextBuf *tb, int  npl, char *fmt, int  nv, double *v);
int      tbf_size       (TextBuf *tb);
char    *tbf_line       (TextBuf *tb, int n);


HexShd *hex_shd_make    (int  l_min);
int     hex_shd_order   (HexShd *shd, int  l_type);
double *hex_shd_alm     (HexShd *shd, int  l_type);
double  hex_shd_mag     (HexShd *shd, int  l_type, int  l_max);
double  hex_shd_ri      (HexShd *shd, int  l_type, int  l);
void    hex_shd_set     (HexShd *shd, int  l_type, int  l_max, double *alm);
void    hex_shd_reset   (HexShd *shd, int  l_type);
void    hex_shd_free    (HexShd *shd);

int    sdf_read_report  (char *filename);
int    sdf_read_multi   (char *filename, int l_max, SdfMol ***sdm_array);
int    sdf_read_all     (char *filename, int max_molecules,
                         int l_max, SdfMol **sdm_array);
void   sdf_free_all     (SdfMol **sdm_array, int n_molecules);
void   sdf_free         (SdfMol *sdm);

SdfFile *sdf_open       (char *filename, char *rw_mode);
void     sdf_close      (SdfFile *sfile);

SdfMol *sdf_read_mol    (SdfFile *sfile, char *line, int l_max);

SdfMol *sdf_make_consensus(char *molname, int l_max);

int    sdf_readln       (SdfFile *sfile, char *line);
int    sdf_readln_anyway(SdfFile *sfile, char *line);
int    sdf_findln       (SdfFile *sfile, char *line, char *keyword);
int    sdf_findln_copy  (SdfFile *src,   char *line, char *keyword,
                         TextBuf *tb);
int    sdf_find_eod     (SdfFile *sfile);
int    sdf_find_shd     (SdfFile *sfile, char *line, char *keyword,
                         double *alm);
int    sdf_read_shd     (SdfFile *sfile, char *line, double *alm);
int    sdf_read_vec     (SdfFile *sfile, char *line, int npl, int width,
                         int nv, double *vec);

int    sdf_testln       (char *cline, char *ckey);

int    sdf_write        (SdfFile *sfile, char *line);
int    sdf_writeln      (SdfFile *sfile, char *line);

int    sdf_lineno       (SdfFile *sfile);
int    sdf_asd          (SdfFile *sfile);
char  *sdf_filename     (SdfFile *sfile);

SdfAtm *sdf_atm_alloc   (int  n_atoms, int  n_bonds);
void    sdf_atm_free    (SdfAtm *sda);

Point3D sdf_cog         (int n_atoms, Point3D *coords,
                         double *rvdw, double *mass);
Point3D sdf_coh         (int n_atoms, Point3D *coords,
                         double *rvdw, Point3D cog, SdfFile *sfile);

int     sdf_find_element(char *name);
char   *sdf_get_element (int element);
double  sdf_get_rvdw    (int element);
double  sdf_get_atwt    (int element);
int     sdf_get_naos    (int element, char *hamiltonian);

HexImage hex_image_alloc(int w, int h);
void     hex_image_free (HexImage *image);
void     hex_image_load (int x0, int y0, HexImage *image);

void     hex_png_write  (char *filename, HexImage *image);
void     hex_jpeg_write (char *filename, HexImage *image);
HexImage hex_jpeg_read  (char *filename);

int  hex_ps_open        (char *filename, int eps, 
                         double *background, int width, int height);
void hex_ps_image       (double x0, double y0, HexImage *image);
void hex_ps_linewidth   (double width);
void hex_ps_pointsize   (double size);
void hex_ps_point       (double v[3]);
void hex_ps_newline     (double v1[3], double v2[3]);
void hex_ps_line        (double v1[3], double v2[3]);
void hex_ps_poly3       (double v1[3], double v2[3], double v3[3]);
void hex_ps_polyn       (int n, double *v);
void hex_ps_char        (double v[3], char c);
void hex_ps_string      (double v[3], char *string);
void hex_ps_pixel       (double v[3], int  pixel);
void hex_ps_colour      (double c[3]);
void hex_ps_close       ();

double sphere_overlap_volume   (double r1, double r2, double r12);
double sphere_centroid_points  (int npts, Point3D *pts, Point3D *centre);
double sphere_centroid_spheres (int npts, Point3D *pts, double *r_pts,
                                Point3D *centroid);

/* new stuff for GRID handling */

int   hex_kont_open(char *filename);
void  hex_kont_dimension(int *xdim, int *ydim, int *zdim);
void  hex_kont_size(float *ga, float *gx, float *gy, float *gz);
void  hex_kont_data(float **xyz_data);
void  hex_kont_close();

void  hex_open_grid_data();
void  hex_grid_test();
int   hex_grid_data(char *res_name, char *atm_name, char *probe_name);
char *probe_num_to_name (int num);
int   mol_num_probes();
char *mol_probe_name(int p);
char *probe_name_to_opposite (char *name);


void hex_rotlib_open    (char *types_file, char *torsions_file);
void hex_rotlib_print   ();
int  hex_rotlib_count   (char *residue_name);
int  hex_rotlib_rotamer (int  chain_no, int  residue_no, int  rotamer_no,
                         Point3D **rotamer_coords);
void hex_rotlib_close   ();

int    hex_is_aa            (char *residue_name);
int    hex_is_dna           (char *residue_name);
char  hex_short_code       (char *residue_name);

void   hex_open_mesh        (int, double, double);
void   hex_enable_grid      ();
int    hex_mesh_order       ();
int    hex_npoles           ();
int    hex_nnodes           ();
int    hex_nedges           ();
int    hex_ncells           ();
int    hex_nbands           ();
int    hex_nfaces           (int order); /* independent of current mesh */
void   hex_set_radius       (double);
void   hex_set_rotation     (double alpha, double beta, double gamma);
void   hex_clear_rotation   ();
Point3D  hex_pole           (Pole);
Point3D  hex_node           (Node);
void   hex_edge_nodes       (Edge, Node *, Node *);
void   hex_edge_poles       (Edge, Pole *, Pole *);
void   hex_node_poles       (Node, Pole *, Pole *, Pole *);
void   hex_pole_nodes       (Pole, Node []);
double hex_pole_area        (Pole);
double hex_node_area        (Node);
void   hex_mesh_angles      (double *pole_angle, double *node_angle);
int    hex_pole_poles       (Pole p, Pole ps[]);
int    hex_node_nodes       (Node n, Node ns[]);
Pole   hex_nearest_pole     (double theta, double phi);
Pole   hex_nearest_pole_new (double theta, double phi);
void   hex_pole_ring_start  (Pole);
Pole  *hex_pole_ring_next   (int  *);
void   hex_pole_ring_end    ();
double hex_grid_theta       (int iband);
int    hex_grid_cell        (double theta, double phi);
void   hex_grid_region      (int, double *, double *, double *, double *);
void   hex_close_mesh       ();

void   setup_platonic       ();
int    platonic_dual        (int  n_vertices);
float *platonic_vertices    (int  n_vertices);
float *platonic_areas       (int  n_vertices);
float  platonic_angle       (int  n_vertices);
int    platonic_resolution  (double angle);
int   *platonic_dual_mesh   (int  n_vertices, int  *n_dimension);
int    platonic_poles       (int  n_vertices);
int    platonic_nodes       (int  n_vertices);

double platonic_set_radius  (double radius);

int    platonic_samples     (int  nv, Euler **es_grid);
int    box_samples          (double box_size, double cell_size, double r12,
                             int *id_zero, Euler **r1, int **tzi, Euler **r2);
int    polar_samples        (double box_size, double cell_size, double r12,
                             Euler **es_grid, int **id_grid);

int    guess_vertices       (int  n_vertices, float **buffer);
int    near_north_points    (int  n_vertices, float angle, Point3D **points);

float *unit_sphere_points   (int  n_vertices);

double *unit_cube_poles     (double edge_length);

void hex_open_atom_data     (char *filename);
int  hex_atom_data          (char *residue, char *atom,
                             float *radius, float *charge, 
                             char *type, sint *opls, sint *ace,
                             byte *ha, byte *hd);

void hex_open_opls_data     (char *filename);
int  hex_opls_data          (int opls_code, double *sigma, double *epsilon);

void   hex_open_dars_data   (char *filename);
double hex_dars_value       (int code1, int code2);
double hex_dars_reduced     (int code1, int code2, int p_max);
double hex_dars_component   (int code1, int p);
double hex_dars_eigenvalue  (int p);
void   hex_close_dars_data  ();

void hex_open_hbond_data    (char *filename);
int  hex_hbond_data         (int pr_code, int lp_code, double *c, double *d);
void hex_open_atom_colours  (char *filename);
int  hex_atom_polar_hydrogens(char *residue, char *h_atoms);
int  hex_atom_colour         (char *residue, char *atom);
int  hex_chain_colour        (char chain);
char *hex_internal_name (char *rname, char *ename);

void hex_open_bond_template  (char *filename);
int  hex_match_bond_template (char *residue_name, int  *b_start, int  *b_end);
void hex_get_bond_template   (int  t, char **a1, char **a2, byte *order);
void hex_open_template   (char *filename);

void hex_dump_template   ();

int  hex_match_template  (int  tpl, char *residue_name);
void hex_get_template    (int  tpl, char *h, 
                          char *e, char *x, char *x1, char *x2);
int  hex_lookup_template (char *residue_name, char *atom_name,
                          char *e, char *x, char *x1, char*x2);
Point3D hex_apply_template(int  j, Point3D x, Point3D x1,  Point3D x2);

int   hex_rgb_dim  ();                          /* NB an ID is NOT an int */
void  hex_rgb_setn (int id);                    /* set current colour */
void  hex_rgb_n2f  (int id, float *floats);     /* ID to floats */
void  hex_rgb_n2b  (int id, byte *bytes);       /* ID to bytes */
int   hex_rgb_b2n  (byte *bytes);               /* bytes to ID */
int   hex_rgb_s2n  (rchar *name);               /* string to ID */
char *hex_rgb_n2s  (int colnum);                /* ID to string */

int   hex_rgba_b2i (byte *bytes);               /* bytes to packed int */
void  hex_rgba_i2b (int rgba, byte *bytes);     /* packed int to bytes */
int   hex_rgba_d2i (double *doubles);           /* doubles to int */
void  hex_rgba_n2d (int id, double *doubles);   /* ID to doubles */
int   hex_rgba_n2i (int id);                    /* ID to packed int */

int  hex_rgb_classic2i(double theta, double phi); /* classic polar colours */

void hex_rgb_blend_f (float  wt1, float  *col1, float  *col2, float  *col3);
void hex_rgb_blend_d (double wt1, double *col1, double *col2, double *col3);

void probe_ray_poles    (int, double, double *, double *);
void atom_ray_poles     (int, Point3D, double, double, double *);
int  pole_ray_probe     (Pole, Vector3D,Vector3D,double, double *);
int  pole_ray_atom      (Pole, Vector3D, Vector3D, int, double, double *);
int  ray_sphere         (Vector3D, Vector3D, double, Point3D *, Point3D *);

int  hex_clusterAgg     (int method, int n_objects, double threshold,
                         char **names, double *dist, int *clusters);

int  hex_atom_pairs     (int  n_atoms, Point3D *pts,
                         double grid_size, Box box, int  **atom_pairs);

int  hex_interface_pairs(int  n1, Point3D *pts1, int  n2, Point3D *pts2,
                         double grid_size, Box box, int  **atom_pairs);

int  hex_accessible_pairs(int  n_atoms, float *xyz, float *rp, Box box,
                          byte accessible[],
                          int  **atom_pairs, float **squares);


/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif  /* hex_util_h */

/*---------------------------------------------------------------------------*/
