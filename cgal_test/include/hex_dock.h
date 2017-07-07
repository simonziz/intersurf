/*-----------------------------------------------------------------------------
**  File:       hex_dock.h
**
**  Author:     Dave Ritchie, 10/05/03
**
**  Purpose:    Header file for the "dock-level" functions.
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

#ifndef hex_dock_h
#define hex_dock_h

#include "hex_types.h"
#include "hex_math.h"
#include "hex_thread.h"
#include "hex_molecule.h"
#include "hex_zlm.h"

//ifdef __cplusplus
//extern "C" {
//endif

/*---------------------------------------------------------------------------*/

typedef struct _scan_block {      /* just less than 32K results block */

   int    n_values;                /* no. useful values in data structure */
   octa   uid[SCAN_BLOCK_SIZE];    /* orientation ID */
   float  energy[SCAN_BLOCK_SIZE]; /* total energy/score */

} ScanBlock;

/*---------------------------------------------------------------------------*/

/*
typedef struct _search_block {    

   int   n_values;                
   int   oid[SEARCH_BLOCK_SIZE]; 
   int   tid[SEARCH_BLOCK_SIZE];
   float t[SEARCH_BLOCK_SIZE]; 
   float s[SEARCH_BLOCK_SIZE];
   float e[SEARCH_BLOCK_SIZE];

} SearchBlock;
*/

/*---------------------------------------------------------------------------*/

typedef struct {

   int    n_values;                /* no. useful values in data structure */

   int    id[MM_BLOCK_SIZE];
   int    bumps[MM_BLOCK_SIZE];

   float  energy[MM_BLOCK_SIZE];
   float  r12[MM_BLOCK_SIZE];

   Euler  er[MM_BLOCK_SIZE];
   Euler  el[MM_BLOCK_SIZE];
/*
   float  alpha2[MM_BLOCK_SIZE];
   float  beta2[MM_BLOCK_SIZE];
   float  gamma2[MM_BLOCK_SIZE];
   float  beta1[MM_BLOCK_SIZE];
   float  gamma1[MM_BLOCK_SIZE];
*/

} MMBlock;

/*---------------------------------------------------------------------------*/

typedef struct _HexProgress {

   int     stage;
   int     n_done;
   int     n_units;
   int     aborting;

} HexProgress;

/*---------------------------------------------------------------------------*/

typedef struct _HexRefConPlan {

   int n_tasks;

   Solution *solution;

} HexRefConPlan;

/*---------------------------------------------------------------------------*/

typedef struct _HexScanProducer {

   HexTaskSlot  task_slot;
   HexProgress *progress;
 
   int      do_scan;
   int      n_vec;
   int      n_max;

   int      a_dim;
   int      ma;

   int      fft_min[10];
   int      fft_max[10];
   int      fft_dim[10];

   int      job_unit;
   octa    *job_list;

   zfloat  *rbuf_cpu;
   zfloat  *rbuf_gpu;

   zfloat  *lbuf_cpu;
   zfloat  *lbuf_gpu;

   zfloat  *lambda_cpu;

   double  *drbuf;
   double  *dlbuf;

} HexScanProducer;

/*---------------------------------------------------------------------------*/

typedef struct _HexRefineProducer {

   HexTaskSlot  task_slot;
   HexProgress *progress;

   int      n_vec;
   int      n_max;

   int      a_dim;
   int      ma;

   int      fft_min[10];
   int      fft_max[10];
   int      fft_dim[10];

   int      job_unit;
   octa    *job_list;

   zfloat  *lambda_cpu;

} HexRefineProducer;

/*---------------------------------------------------------------------------*/

typedef struct _HexDockConsumer {

   int      n_tasks;
   float    e0;

} HexDockConsumer;

/*---------------------------------------------------------------------------*/

typedef struct _HexSkinConsumer {

   int      n_tasks;
   int      n_done;

} HexSkinConsumer, 
  HexRotConsumer;

/*---------------------------------------------------------------------------*/

typedef struct _HexRotatePlan {

   HexTaskSlot  task_slot;

   HexRotConsumer *consumer;

   int      n_vec;
   int      n_max;

   int      job_unit;
   octa    *job_list;

} HexRotatePlan;

/*---------------------------------------------------------------------------*/

typedef struct _HexSkinWorker {

   HexTaskSlot  task_slot;

   HexProgress     *progress;
   HexSkinConsumer *consumer;

   int      n_max;
   int      n_vec;
   int      n_dars;

   int      ng;
   double   dg;

   Point3D  origin;

   byte    *sbuf;

   double  *abuf;
   double  *bbuf;

   int     *job_list;

   double  *asig;
   double  *atau;
   double  *apsi;

   double  *bsig;
   double  *btau;
   double  *bpsi;

} HexSkinWorker;



/*---------------------------------------------------------------------------*/

void hex_docker                (DockSpec *ds, Docking *docking,
                                Molecule *molecule1, Molecule *molecule2);

void hex_dock_poses            (DockSpec *ds, Docking *docking,
                                Molecule *mol1, Molecule *mol2, int mdl1, int mdl2);

void hex_setDockSpec           (DockSpec *ds);
void hex_setMatchSpec          (DockSpec *ds);

DockData *hex_setDockData      (DockSpec *ds, Molecule *mol1, Molecule *mol2, 
                                int mdl1, int mdl2, int with_message);

void hex_freeDockData          (DockData *dd);

void hex_matcher               (DockSpec *ds, Docking *docking,
                                Molecule *molecule1, Molecule *molecule2);

void hex_matcher_test          (DockSpec *ds, Docking *docking,
                                Molecule *molecule1, Molecule *molecule2);

int  hex_sixd_search           (DockSpec *ds, DockData *dd, Solution *solution);

void hex_alloc_solutions       (Docking *docking, int  max_solutions);
void hex_free_solutions        (Docking *docking);

void hex_coeffs                (DockSpec *ds, 
                                Molecule *molecule1, Molecule *molecule2);

void hex_ligand_move_anlm      (int n_max, int n_vec, int vec_types[],
                                double r12,    Euler  reuler, Euler leuler,
                                double *bnlm,  double *bnlm_new);

void hex_ligand_transform_anlm(int n_max, int n_vec, int *vec_types,
                               Matrix3D t_ligand,
                               double *bnlm,  double *bnlm_new);

double hex_ligand_score        (int n_max, int n_vec, int vec_types[],
                                Matrix3D t_ligand,
                                double *anlm,  double *bnlm);

double hex_point_score         (int n_max, int n_vec, int vec_types[],
                                double r12,    Euler  reuler, Euler leuler,
                                double *anlm,  double *bnlm);

double hex_point_energy        (int n_max, int n_vec, int vec_types[],
                                double r12,    Euler reuler, Euler leuler,
                                double *anlm_phi, double *anlm_rho,
                                double *bnlm_phi, double *bnlm_rho);

int  hex_check_threefold        (double r12,    Euler er, Euler el);

/*
double hex_point_score         (int n_max, int n_vec, int vec_types[],
                                double r12,    double beta1, double gamma1,
                                double alpha2, double beta2, double gamma2,
                                double *anlm,  double *bnlm);

double hex_point_energy        (int n_max, int n_vec, int vec_types[],
                                double r12,    double beta1, double gamma1,
                                double alpha2, double beta2, double gamma2,
                                double *anlm_phi, double *anlm_rho,
                                double *bnlm_phi, double *bnlm_rho);
*/

void tnl_all_norm              (int  n_max, double *tnn);

void hnl_all_norm              (int  n_max, double *hnn);

void hex_function_surface      (Molecule *molecule, int mdl, int with_hydrogens,
                                int culling, int make_msurf, int surface_type, 
                                int  n_max, int  j_max,
                                double grid_size, double softness, 
                                double z_shift, 
                                MeshSurface *msurf);

void hex_setup_electro         (Molecule *molecule, int mdl, int  n_max);
void hex_setup_grid            (Molecule *molecule, int mdl, int  n_max);

void hex_electro_fn            (int  n_pts, int  n_max, 
                                double *phi_nlm, 
                                double *rho_nlm,
                                Point3D *pts, 
                                float *phi, float *rho);

void hex_steric_fn             (int  n_pts, int  n_max, 
                                double *pb1_nlm, 
                                double *pb2_nlm,
                                double *pb3_nlm,
                                double *pb4_nlm,
                                Point3D *pts, 
                                float *pb1, float *pb2, float *pb3, float*pb4);

void hex_grid_fn               (Molecule *molecule, int model_no,
                                int function_type,
                                int  n_max, double *anlm, int z_sym,
                                double grid_size, Point3D grid_origin,
                                int  ng, float *grid);

void hex_steric_compute        (Molecule *molecule, int mdl, int which);
void hex_electro_compute       (Molecule *molecule);
void hex_grid_compute          (Molecule *molecule);

void hex_skin_calc         (Molecule *molecule, int model_no, 
                            int  n_max, int  ng, double dg, Point3D origin, 
                            double skina, double skinp, 
                            byte isig, byte itau, byte *sbuf);

void hex_skin_calc_new     (Molecule *molecule, int model_no, 
                            DockingAtoms *da,
                            int  n_max, int  ng, double dg, Point3D origin, 
                            byte isig, byte itau, byte *sbuf);

void hex_skin_alloc        (Molecule *molecule, int model, int  n_max);

void hex_skin_copy         (Molecule *molecule, int from_model, int to_model, 
                            int  n_max);

void hex_skin_laplacian    (int  n_max, double *anlm,  double *anlml);

void hex_envelope_integrate(int  j_max, 
                            double *r_poles, double *r_nodes,
                            double *eta);

void hex_envelope_integrate2(int  j_max, 
                            double *r_poles1, double *r_nodes1,
                            double *r_poles2, double *r_nodes2,
                            double *eta1, double *eta2);

void   hex_pe_piper        (Molecule *, int model_no, int  n_max,
                            int which, char *dxname);

void   hex_gtau_calc3      (Molecule *molA, int model_no, 
                            Molecule *molB, int  n_max);

void   hex_pe_calc3        (Molecule *, int model_no, int  n_max,
                            double cell_size, double max_radius);

void   hex_pe_copy         (Molecule *, int from_model, int to_model, 
                            int  n_max, int n_dars);

void   hex_pe_alloc        (Molecule *, int model, int  n_max, int n_dars);

void   hex_es_calc         (Molecule *, int model_no, int  n_max,
                            double max_radius);

void   hex_es_alloc        (Molecule *, int model, int  n_max);

void   hex_es_copy         (Molecule *, int from_model, int to_model, 
                            int  n_max);

int    hex_es_test         (Molecule *, int n_max);
void   hex_es_exact_potential(Molecule *);
void   hex_es_exact_energy (Molecule *mol1, Molecule *mol2, double r12);

void hex_docker_reduce      (int n_vec, int n_max, 
                             double *anlm_sig, double *anlm_tau,
                             double *bnlm_que, double *bnlm_tau,
                             double *fs_cos,   double *fs_sin);

void hex_float_reduce      (int n_vec, int n_max, 
                             float *anlm_sig, float *anlm_tau,
                             float *bnlm_que, float *bnlm_tau,
                             float *fs_cos,   float *fs_sin);

void hex_basic_reduce       (int  n_max, double *anlm_sig, double *anlm_tau,
                                         double *bnlm_que, double *bnlm_tau,
                             double *fs_cos, double *fs_sin);

void hex_matcher_reduce     (int  n_max, double *anlm_tau, double *bnlm_tau,
                             double *fs_cos, double *fs_sin);

void hex_min_picker         (int  n, int  m, double *q, byte *peaks);

void hex_min_scanner        (int  n, int  m, double *q, int  *pos, float *val);
void hex_min_scanner_float  (int  n, int  m, float  *q, int  *pos, float *val);
void hex_max_scanner        (int  n, int  m, double *q, int  *pos, float *val);

void hex_delta_density   (int n_max, double radius, hexa *b);
void hex_sphere_density  (int n_max, double radius, hexa *b);
void hex_atom_density    (int n_max, double radius, hexa *b);
void plot_atom_density   (int n_max, double r_min, double r_max, double dr,
                          double *anlm);
void plot_gauss_density  (double radius);
double test_gauss_rms    (int n_max, double radius, hexa *b);
void plot_lj_energy      (double a12, double c12,
                          double r_min, double r_max, double dr);
void plot_hb_energy      (double c12, double d10,
                          double r_min, double r_max, double dr);
void setup_lj_example    (double *a12, double *c12);
void fit_lj_energy       (int n_max, double r_min, double r_max, double dr,
                          double dz);
void fit_lj_delta        (int n_max, double r_min, double r_max, double dr,
                          double dz);

void hex_refine          (Molecule *mol1, Molecule *mol2, 
                          int mdl1, int mdl2, int docking_refine, 
                          int ns, Solution *solution,
                          double e_factor, double es_scale);

// int    hex_refine        (Molecule *mol1, Molecule *mol2, 
//                           int mdl1, int mdl2, int docking_refine, int ns,
//                           float *t_top, float *s_top, float *e_top, 
//                           float *r_top, Euler *er_top, Euler *el_top, int *b_top,
//                           double e_factor, double es_scale);

void ylm_open_solver       (Docking *docking,
                            int docking_mode, Molecule *mol1, Molecule *mol2,
                            int max_solutions, 
                            int propagate_solutions, int docking_bias,
                            int min_harmonic, int max_harmonic, 
                            int inc_harmonic,
                            double contact_angle, 
                            double seg_area, double gyration_angle);

void ylm_close_solver      ();

void ylm_solve_ellipsoid   (Molecule *molecule);

void ylm_guess_similarity  (int order, int n_guess, int n_alpha);

void ylm_refine_similarity (int start, int end, int step);

void   hex_task            (int parent_pid, int child_no); /* task controller */

void save_docking_atoms    (char *filename, Molecule *molecule, int mdl);
void write_docking_atoms   (HexFile *unit, Molecule *molecule, int mdl);


DockingAtoms *hex_docking_atoms_get(Molecule *molecule, int mdl,
                                double skina_polar, double skin_polar);

void hex_docking_atoms_free(DockingAtoms *da);

int  hex_fft_bandlimit     (double angle);
int  hex_fft_pow2          (double angle);

/*---------------------------------------------------------------------------*/

void set_piper_weights(double w0, double w1, double w2, double w3);

/*---------------------------------------------------------------------------*/

// ifdef __cplusplus
// }
// endif

#endif  /* hex_dock_h */

/*---------------------------------------------------------------------------*/
