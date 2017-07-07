/*-----------------------------------------------------------------------------
**  File:       hex_molecule.h
**
**  Author:     Dave Ritchie, 10/05/03
**
**  Purpose:    #include file for the "mol-level" molecule handling functions.
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

#ifndef hex_molecule_h
#define hex_molecule_h

#include "hex_mol.h"

/*---------------------------------------------------------------------------*/

//ifdef __cplusplus
//extern "C" {
//endif

/*---------------------------------------------------------------------------*/

typedef struct _ramp {        /* max_val > top_val > bot_val > min_val */

   char   name[40];

   double min_val;
   double max_val;
   double bot_val;
   double top_val;

   Vector3D min_col;
   Vector3D mid_col;
   Vector3D max_col;

   Vector3D min_col_save;
   Vector3D mid_col_save;
   Vector3D max_col_save;

} Ramp;



typedef struct _gsurface {  /* stores the geometry of a gnomonic surface */

   int       n_poles;
   int       n_nodes;

   Point3D  *pole_points;
   Point3D  *node_points;
   Vector3D *pole_normals;
   Vector3D *node_normals;

   int      *pole_colours;
   int      *node_colours;

   double   *pole_values;
   double   *node_values;

   double    mep_min;
   double    mep_max;
   double    iel_min;
   double    iel_max;
   double    eal_min;
   double    eal_max;
   double    alpha_min;
   double    alpha_max;

   int       colour;          /* avoid-recalculation flags */
   int       model_no;
   int       ramp_no;
   int       l_max;
   int       p_max;
   int       hydrogen;

} Gsurface;


typedef struct _molecule_style { /* controls the display of a molecule */

   int     display;          /* true/false */
   int     colour_scheme;    /* molecule colour scheme */
   int     solvent_colour;   /* solvent colour scheme */
   int     hetero_colour;    /* hetero atom colour scheme */
   int     inaccessible_col; /* special colour for "inaccessible" atoms */
   int     solid_fill;       /* solid fill atoms */
   int     solid_res;        /* of filled spheres */
   int     draw_backbone;    /* true/false */
   int     draw_sidechain;   /* true/false */
   int     draw_solvent;     /* true/false */
   int     draw_hetero;      /* true/false */
   int     draw_hydrogen;    /* true/false */
   int     draw_hbonds;      /* true/false */
   int     symmetry_type;    /* 0=none, 1=crystal, 2=biological */
   int     enable_solvent;   /* true/false */
   int     enable_hetero;    /* true/false */
   int     enable_arglys;    /* true/false */
   int     bond_width;       /* in pixel units, default = 2 */
   int     post_id;          /* draw residue ID */
   int     ellipsoid_arcs;   /* draw molecular ellipsoid */
   int     ellipsoid_axes;   /* draw ellipsoidal axes */
   int     molecular_origin; /* draw molecular origin */
   int     radial_filter;    /* apply radial filtering */
   double  radial_filter_dist; /* the current radial filter distance */
   Point3D   radial_filter_pt; /* the molecule-zero radial filter centroid */
   char    colour_file[MAX_PATHNAME];  /* most recent user colour file */

} MoleculeStyle;

typedef struct _cartoon_style { /* controls protein cartoon display */

   int     display;          /* 0:off, 1:on */
   int     style;            /* 0: edged ribbon, 1: solid ribbon, 2: tube */
   int     sse_method;       /* 0 = stride, 1 = dssp, 2 = kpax */

   int     resolution;       /* 0: low, 1:medium, 2:high, 3:ultra */
   int     colour_scheme;    /* 0: structure, 1: chains */
   int     helix_colour;     /* colour ID for alpha helices */
   int     sheet_colour;     /* colour ID for beta strands */
   int     turn_colour;      /* colour ID for beta turns */
   int     three10_colour;   /* colour ID for 310/pi helices */
   int     loop_colour;      /* colour ID for loops */
   int     helix_edge_colour;/* colour ID for helix edges */
   int     sheet_edge_colour;/* colour ID for sheet edges */
   int     turn_edge_colour; /* colour ID for turn edges */
   int     three10_edge_colour;/* colour ID for 310/pi edges */
   double  radius;           /* loop radius = 1/2 sheet thickness */
   double  width;            /* width of sheets */


} MoleculeCartoon;

typedef struct _poly_style { /* controls surface polygon calculation/display */

   int     display;          /* 0=off, 1=on (only applies to harmonics) */
   int     style;            /* 0=wire-frame, 1=polyframe etc. */

   int     edge_colour;      /* wire frame edge colour */
   int     thickness;        /* line thickness */
   int     colour;           /* colour by: 
                                0 Hex-classic, 1 Hex-chains, 
                                9 PS-classic, 
                                10 PS-MEP, 11, PS-IEL, 
                                12 PS-EAL, 13 PS-ALPHA */
   int     hex_mode;         /* mesh type 0=triangular, 1=hex */
   int     transparency;     /* 0=alpha blending off, 1=on */
   int     cull;             /* OBSOLETE - backface culling: 1=yes, 0=no */
   double  alpha;            /* for alpha blending */

} MoleculePoly;

typedef struct _solid_style { /* controls 3D solid calculation/display */

  int    display;            /* 0: off, 1: on */
  int    style;              /* 1:VDW, 2:ball/stick, 3:licorice */

  double aromatic_thickness;
  double peptide_thickness;

  double sphere_size;
  double cylinder_size;

  int colour_scheme;

  int show_aromatic_rings;
  int show_peptide_planes;

  int aromatic_colour;
  int peptide_colour;
  int aromatic_facets;

} MoleculeSolid;


typedef struct _dot_style { /* controls molecular surface dot display */

   int     display;          /* 0=off, 1=on */
   int     size;             /* dot thickness in pixel units */
   int     colour;           /* dot colour */
   int     style;            /* 3=dots */
   double  probe;            /* probe radius used during surface sampling */
   double  skin;             /* skin thickness: 0=molecular,1=accessible */
   double  grid;             /* contouring grid size */

} MoleculeDots;

typedef struct _mesh_style { /* controls triangulated surface mesh display */

   int     display;          /* 0=off, 1=on */
   int     style;            /* 0=wire, 1=poly, 2=wirepoly, 3=dots */
   int     colour;           /* 0=classic, 1=atoms, 2=potential, 3=charge */
   int     draw;             /* 1=primary, 2=secondary, 0=both */
   int     cull;             /* 1=backface culling, 0=off */
   int     surface_type;     /* requested surface: 0=gaussian, 1=sigma, 2=tau */
   int     transparency;     /* 0=alpha blending off, 1=on */
   double  grid_size;        /* sampling grid size: def = 0.75A */
   double  softness;         /* Gaussian softness factor: 1.0 = default */
   double  viewing_probe;    /* probe radius when viewing */
   double  apolar_probe;     /* probe radius when docking */
   double  polar_probe;      /* probe radius when docking */
   double  sigma_tau_shift;  /* translate sigma tau surfaces in Z */
   double  alpha;            /* for alpha blending */

} MoleculeMesh;

typedef struct _hex_mesh     { /* stored surface mesh */

   int       nvertex;          /* total no. calculated vertices */
   int       ntriangle;        /* total no. surface triangles */
   int       nsurface;         /* total no. distinct surfaces */

   Point3D  *vertices;         /* vertex coordinates as x,y,z triplets */
   Vector3D *normals;          /* triangle normal vectors (normalised) */
   int      *triangles;        /* triangles as triplets of vertex numbers */
   int      *neighbours;       /* neighbours as triplets of triangle numbers */
   sint     *surfaces;         /* surface Id no.s of triangles */

} HexSurface;

typedef struct _mesh_surface { /* stored surface mesh */

   HexSurface   hsurf;            /* the contoured surface mesh */

   int       model_no;         /* model no. when surface calculated */
   int       surface_type;     /* surface type when surface calculated */
   int       with_hydrogens;   /* hydrogen mode when surface calculated */
   int       n_max;            /* Order N when sigma,tau,phi,rho calculated */
   int       culling;          /* cull edges state when surface calculated*/
   double    grid_size;        /* grid size when surface calculated */
   double    softness;         /* softness when surface calculated */
   double    apolar_probe;     /* probe radius when surface calculated */
   double    polar_probe;      /* probe radius when surface calculated */
   double    rmax;             /* radius of sphere that encloses surface */
   double    sigma_tau_shift;  /* current sigma tau Z translation */

   double    z_shift;          /* for receptor shift in protein-ligand mode */

   int      *vcolours;         /* vertex colors packed as 8-bit RGB values */

   float     phi_min;          /* electrostatic potential data range */
   float     phi_max;
   float     rho_min;          /* defines charge density data range */
   float     rho_max;
   float     pb1_min;          /* data ranges for three probes */
   float     pb1_max;
   float     pb2_min;
   float     pb2_max;
   float     pb3_min;
   float     pb3_max;
   float     pb4_min;
   float     pb4_max;

   float    *phi;              /* electrostatic potential at each vertex */
   float    *rho;              /* electrostatic charge at each vertex */
   float    *pb1;              /* probe 1 energy at each vertex, etc. */
   float    *pb2;
   float    *pb3;
   float    *pb4;

} MeshSurface;

typedef struct _docking_atoms {

   int       npts;
   float     rmax;
   Point3D  *apts;
   float    *rvdw;
   float    *rsas;

} DockingAtoms;

typedef struct _rotamer {   /* the data structure of a rotamer */

   struct _rotamer *next;   /* pointer to next rotamer */
   int    chain_no;         /* chain no. in pdb structure */
   int    residue_no;       /* residue no. in pdb structure */
   int    atom_no;          /* first atom no. in pdb structure */
   int    rotamer_no;       /* -1=> saved original pdb coords, else 0,1... */
   int    n_atoms;          /* length of atom coordinate vector */
   Point3D *atom_coords;      /* rotamer atom coords in order of pdb file */

} Rotamer;

typedef struct _model {     /* for model-specific docking */

   int    id;               /* ID no - normally allocated sequentially */
   int    visible;          /* do we display the model - a bit of a hack */

   char   label[MAX_PATHNAME]; /* model name (from PDB model or SDF compound) */

   int    a_min;            /* 1st and last atom no.s in the model */
   int    a_max;

   Matrix3D t_model;        /* places the model in space after any "edits" */

   Matrix3D t_predock;      /* places the model at origin prior to docking */

   Matrix3D t_postdock;     /* transform the complex at origin after docking */

   Matrix3D t_canonical;    /* rotate molecule into canonical pose */

   int  canonical_done;     /* 1 => canonical orientation calculated */

   int  ellipsoid_done;     /* 1 => ellipsoid calculated */

   double ellipsoid_x;      /* ellipsoidal radii */
   double ellipsoid_y;
   double ellipsoid_z;


   Point3D   origin;          /* coords of docking origin - normally zero */

   double *sig;             /* exterior skin expansion coefficients */
   double *tau;             /* interior skin expansion coefficients */

   double *phi;             /* electrostatic expansion coefficients */
   double *rho;             /* electrostatic expansion coefficients */

   double *psi;             /* DARS/GRID potential vectors */
   double *del;             /* DARS/GRID density functions */

/* double *eta; */          /* radial distance coefficients */

   HexShd *shd;             /* molecule-specific spherical harmonic data */

   double  gnome_softness;   /* parameters at which gnome last calculated */
   double  gnome_probe;
   double  gnome_grid;
   
   double *gnome_poles;    /* gnomonic projection distances (poles) */


} Model;

typedef struct _solution {  /* data structure for each docking solution */

   sint    model[2];        /* model no.s of this solution */
   octa    id;              /* solution ID as a 8-byte code number */
   float   e_total;         /* total docking energy in KJ/mol */
   float   e_shape;         /* shape contribution to e_total, KJ/mol */
   float   e_force;         /* force-field or electrostatic contribution to e_total, KJ/mol */
   float   e_air;           /* air boost contribution to e_total, KJ/mol */
   float   v_overlap;       /* skin-interior overlap in A^3 */
   float   v_penalty;       /* interior-interior overlap in A^3 */
   float   rms;             /* RMS of docking orientation */
   float   r12;             /* intermolecular distance */
   Euler   reuler;          /* 0 => receptor & complex, 1 => ligand ... */
   Euler   leuler;          /* 0 => receptor & complex, 1 => ligand ... */
   sint    bumps;           /* no. main-chain bumps in solution */
   sint    flag;            /* post-processing: 0 => all OK */
   sint    omega;           /* Cn symmetry index number from symmetry docking */
   sint    pose;            /* base pose number for Dn symmetry dockig */
   Rotamer *rotamers[2];    /* NOT USED - per-solution linked list of rotamers */

} Solution;

typedef struct _docking_data {   /* base receptor (a) and ligand (b) docking data */

   int       n_pose;             /* number of fused symmetry poses (groups of coefficients) */
   Matrix3D *mat_poseA;          /* symmetry pose set-up matrix for the "receptor" */
   Matrix3D *mat_poseB;          /* symmetry pose set-up matrix for the "ligand" */
   float    *e_pose;             /* pair-wise energy per pose */

   int       n_point;            /* no. of point orientations to evaluate, may be zero */
   Matrix3D *mat_pointA;         /* the transformation to apply to the "A" structure */
   Matrix3D *mat_pointB;         /* the transformation to apply to the "B" structure */

   int       n_vec;              /* number of coefficient vectors */
   int      *vec_types;          /* coeff vector types, 1 => GTO, 0 => ETO, etc. */
   double   *anlm_phi;           /* coeff vectors, size = n_vec * NLM_DIM(n_max) */
   double   *anlm_rho;           /* coeffs are stored in nlm-natural order */
   double   *bnlm_phi;
   double   *bnlm_rho;

} DockData;

/* NB. quite a few flags from hex_render.c should go in here */

typedef struct _docking {   /* docking/superposition-specific results block */

   int     max_solutions;   /* no. solution structures allocated */
   int     num_solutions;   /* current total no. docking solutions */
   int     num_clusters;    /* current no. of clustered docking solutions */

   DockData *dd;            /* the docking data */

   Solution *solution;      /* solution structures allocated later */
   int      *order;         /* sort order for solutions */
   int      *cluster;       /* cluster no. for each solution */
   int      *tops;          /* index of 1st member of each cluster */

} Docking;

typedef struct _docking_spec {   /* controls docking calculations */

   int    max_solutions;         /* max no. solution structures to allocate */
   int    problem_type;          /* 1=similarity, 2=docking, 3=multidocking */
   int    potential_type;        /* 1 = electrostatics, 2 = GRID (exptl) */
   int    threefold_hack;        /* 1 if hacking */
   int    n_min;
   int    n_max;                 /* docking skin expansion order */
   int    j_max;                 /* shape matching expansion order */
   int    fft_type;              /* FFT correlation type 1=1D ... 5=5D */
   int    fft_device;            /* FFT device: 0=CPU, 1=GPU, 2=BOTH */
   int    search_mode;           /* 0: rotational search; 1: translational */
   int    refine;                /* 1 = MM energy, 2 = MM refine */
   int    volume;                /* 1 = calculate final overlap volumes */
   int    receptor_samples;
   int    ligand_samples;
   int    alpha_samples;
   int    receptor_bandlimit;
   int    ligand_bandlimit;
   int    alpha_bandlimit;
   int    r12_substeps;          /* distance substeps during refinement = 2 */
   int    air_mode;              /* 1: apply restraints; 0: no restraints */

   int    sym_omega_method;      /* 0: use supplied value; 1: use supplied range */
   int    sym_omega1;            /* first symmetry order for symmetry docking */
   int    sym_omega2;            /* last  symmetry order for symmetry docking */
   double sym_omega_value;       /* fixed "y" rotation angle for 2D symmetry FFT */

   int    sym_local;             /* 1: local frame coordinates; 0: world frame */
   double sym_LtoW;              /* local-to-world scale factor */
                                 /* = 1/(2xcos((pi+omega)/2)) for Cn */
                                 /* = 1/sqrt(3/8) for Tetrahedral, etc. */
   double sym_phiA;              /* extra z-rotations for high symmetry modes */
   double sym_phiB;              /* (always zero except for T, O, and I cases) */
   

   double score_scale;           /* scoring scale factor (mainly for symmetry) */
   double score_threshold;       /* docking energies worse than this are not saved */

   double receptor_stepsize;
   double ligand_stepsize;
   double alpha_stepsize;

   double air_energy;            /* restraint energy boost */
   double air_distance;          /* restraint threshold distance */


   double r12_step;              /* distance step size = 1.0 */
   double r12_guess;
   double r12_range;             /* +/- 40 Angstrom and start guess */
   double grid_size;
   double potential_grid_size;
   double integration_radius;
   double receptor_range_angle;
   double ligand_range_angle;
   double alpha_range_angle;

   int    sampling_method;       /* 0=range angles, 1=receptor box */
   double box_size;              /* size of receptor box, def=6 Angstrom */
   double z_shift;

   double p_factor;              /* penalty for skin/skin overlap (def=0.0) */
   double q_factor;              /* dimensionless steric penalty term (def=12.0) */
   double e_factor;              /* grid energy scaling factor (def=0.2) */
   double f_factor;              /* scaling factor for 3-probe (def=0.6) */
   double g_factor;              /* desolvation energy density J/A^3 (def=700) */
   double s_factor;              /* shape similarity cross term (def=0.0) */
   double t_factor;              /* electrostatic test scale factor (def=1.0) */
   double z_factor;              /* dimensionless PIPER scale factor (def=5.0) */


} DockSpec;

typedef struct _cluster_spec {   /* controls clustering */
   int    max_clusters;          /* max no. clusters to consider for output */
   int    window;                /* cluster look-ahead window size */
   int    sort_mode;             /* 0: by energy; 1: by cluster */
   int    view_mode;             /* 0: view all;  1: view best */
   int    bumps;                 /* no. bumps to allow in solutions */
   double threshold;             /* RMS clustering threshold */
} ClusterState;


typedef struct _normals {   /* for "macro-molecule" surface sampling */

  int  n_pts;           
  Point3D *top_pts;           /* each "top" and "bot" point defines a vector */
  Point3D *bot_pts;           /* normal to the molecular surface */

} MoleculeNormal;

typedef struct _macro_state {  /* for "macro-molecule" surface sampling */

   int     display_macro_normals; /* display surface normal vectors */
   int     display_macro_spheres; /* display surface spheres  */
   int     macro_auto_coverage;   /* sphere coverage */
   int     macro_min_coverage;    /* sphere coverage */
   double  macro_sphere_radius;   /* display surface spheres  */
   double  macro_weighting;       /* score weighting */
   double  m12_distance;          /* multi-docking starting distance */

} MacroState;

typedef struct _ylm_spec {  /* structure for spherical harmonics state */

   int    Npoles;
   int    Nnodes;
   int    Nedges;

   int    harmonic_display;
   int    harmonic_order;
   int    harmonic_degree;
   int    harmonic_recursion;

   int    harmonic_property;

   double harmonic_angle_alpha;
   double harmonic_angle_beta;
   double harmonic_angle_gamma;

} YlmState;

typedef struct _molecule {  /* the data structure of a molecule */

   int     open;            /* is the rest of the structure legal */
   int     id;              /* molecule id: 0=receptor, 1;ligand */
   int     mol_index;       /* pdb handle to use for current molecule */
   int     mol_type;        /* MOL_PDB_TYPE or MOL_SDF_TYPE */

   int     is3DBlast;        /* 0: normal molecule, 1: special 3D-Blast PDB */

   struct _molecule *clone_molecule; /* ptr to clone parent molecule */
   int               clone_model;    /* model no. of clone parent */

   char    mol_file[MAX_PATHNAME];   /* where molecule came from */
   char    mol_code[MAX_PATHNAME];   /* PDB code or molecule name */

/* molecule-specific coordinate data */

   double    radius;         /* centroid + radius defines bounding box */

   Point3D   centroid;        /* coordinate centroid (initial c.o.m.) */
   Point3D   real_origin;     /* working c.o.m. (after user edits etc.) */
   Point3D   world_origin;    /* position of molecule in graphics world */

   int     origin_atom;     /* atom to use for origin; */
   int     interface_atom;  /* atom to use for z-axis pre-orientation; */
   double  interface_radius; /* distance from origin to interface */

   Euler   sh_rot;           /* for spherical harmonic test rotations (disfavoured) */
/*
   Angle   euler_alpha;    
   Angle   euler_beta;
   Angle   euler_gamma;
*/

/* molecule-specific atom/rotamer quantities */

   byte    *atom_id_text;    /* atom id display flags */

   Rotamer *base_rotamers;  /* for manually added rotamers (not docking) */
   Rotamer *saved_rotamers; /* saved atom coords prior to rotamer substn. */

/* molecule-specific flags & counters */

   int     rnl_order;       /* highest N for steric coefficients */
   int     tnl_order;       /* highest N for electrostatic coefficients */
   int     enl_order;       /* highest N for grid energy coefficients */
   int     num_energies;    /* no. GRID energy values (exptl) */

   int     solvent_culled;  /* 1 if accessible solvent has been culled */

   int    gnome_damage;     /* nasty re-draw flag */
   int    mesh_damage;      /* nasty re-draw flag */

   int    gnome_dlist;      /* display list for harmonic/gnomonic surfaces */
   int    dots_dlist;       /* display list for dots */
   int    atom_dlist;       /* display list for atoms, bonds and labels */
   int    solid_dlist;      /* display list for solid model rendering */
   int    normal_dlist;     /* display list for surface normals */
   int    mesh_dlist;       /* display list for triangulated surface mesh */
   int    cartoon_dlist;    /* display list for protein cartoons */

   int    model_type;       /* 0: none, 1: PDB (e.g. NMR) 2: macrodocking,
                               3: user clone (exptl) */

   int    num_models;       /* no. of "models" in PDB: always at least 1 */
   int    display_model;    /* model no. to use for display: 0=> all */
   int    docking_model;    /* model no. to use for docking: 0=> all */
   int    viewing_model;    /* model no. to display for docking solution */

   int    n_active_airs;    /* no. active restraint residues */
   int    n_passive_airs;   /* no. passive restraint residues */

   long  *active_airs;      /* lists of airs (may be null) */
   long  *passive_airs;

   Model *model;            /* coefficient vectors for each model, etc. */

   MeshSurface    msurf1;    /* surface mesh for docking */

   MeshSurface    msurf2;    /* surface mesh for viewing */

   MeshSurface    msurf3;    /* surface mesh for gnomonic dot surfaces */

   MoleculeStyle  style;     /* style sub-structure */

   MoleculeDots   dot;       /* dot style sub-structure */

   MoleculePoly   poly;      /* polygon style sub-structure */

   MoleculeSolid  solid;     /* solid style sub-structure */

   MoleculeCartoon cartoon;  /* cartoon style sub-structure */

   MoleculeNormal normal;    /* surface normal vectors sub-structure */

   MoleculeMesh   mesh;      /* mesh style sub-structure */

   Matrix3D t_base;

   Gsurface       gsurface;

} Molecule;


/*---------------------------------------------------------------------------*/

typedef struct _graphics_state { /* maintains general graphics state */

   Matrix3D scn_rot_matrix;
   Matrix3D scn_trn_matrix;
   Matrix3D lok_rot_matrix;
   Matrix3D lok_trn_matrix;

   Matrix3D mol_rot_matrix[MAX_MOL];   /* molecule rotation */
   Matrix3D mol_trn_matrix[MAX_MOL];   /* molecule translation */
   Matrix3D mol_org_matrix[MAX_MOL];   /* molecule origin translation */

   double zoom_factor;
   double world_scale_x;            /* size of viewport in world coords */
   double world_scale_y;
   int    world_scale;              /* =1 if size has been calculated */

   int    solid_resolution;          /* 0:low, 1:medium, 2:high */
   int    print_buffer;              /* print buffer size: 16Mb */

   int   Background;
   int   Foreground;                 /* will get values at start-up time */
   int   OffWhite;
   int   Orange;
   int   Green;
   int   Yellow;
   int   Cyan;

   int   pointer_mode;              /* 1 = translate/rotate scene, 
                                       2 = translate/rotate plane, 
                                       3 = pick atom ID, 
                                       4 = edit molecule orientation,
                                       5 = edit molecule origin,
                                       6 = set scene origin */

   int   pointer_focus;             /* 1 or 2 (temorary editing toggle) */
   int   button_motion;
   int   motion_fill;
   int   projection_type;           /* 0=orthographic, 1=perspective */

   int    movie_activated;
   int    movie_type;
   int    movie_cycle;
   double movie_rate;

   int    spin_enabled;
   int    spin_activated;
   int    spin_button;
   double spin_x_motion;
   double spin_y_motion;
   double spin_angle;

   double draw_time;
   double delay_time;
   double running_draw_time;
   int    running_draw_count;

   float left_ortho;
   float right_ortho;
   float bottom_ortho;
   float top_ortho;
   float ortho_width;
   float ortho_height;
   float z_near;
   float z_far;
   float dummy_z_dist;       /* not used */
   float pixel_parallax;
   float world_parallax;
   float surface_line_thickness;

   int   test_dlist;
   int   full_screen;
   int   window_width;
   int   window_height;
   int   use_fog;
   int   use_parallax;
   int   test_hex_mesh;
   int   surface_density;
   int   n_max;                    /* for drawing SPF densities */
   int   j_max;                    /* for drawing SPF densities */
   int   p_max;                    /* for drawing ParaSurf SH properties */
   int   display_background;       /* 1: show background perspective effect */
   int   display_mode;             /* -1=cube,0=test,1=normal,2=logo */
   int   display_axes;
   int   display_score;            /* display docking score in graphics */
   int   moving_molecule;          /* -1=scene, 0=receptor, 1=ligand */
   int   moving_origin;            /* -1=scene, 0=receptor, 1=ligand */
   int   origin_display;           /* =1 when editing origins */
   int   pending_edits;            /* =1 when molecular edits pending */
   int   surface_display;          /* harmonic/gnomonic surfaces */
   int   test_surface_style;       /* wireframe */
   int   molecular_axis;
   int   sld_scope;               /* apply style changes to both */
   int   mol_scope;               /* apply style changes to both */
   int   sol_scope;
   int   dot_scope;
   int   srf_scope;
   int   clp_scope;
   int   edt_scope;               /* 0=edit all, 1=edit model, 2=clone model */
   int   msh_scope;
   int   crt_scope;
   int   atom_density;          /* no. interior atom samples */
   int   user_solution;         /* user solution number */
   int   real_solution;         /* internal solution number */
   int   symmetry_element;      /* symmetry transform for output */
   int   batch_mode;            /* 1 = no graphics or terminal input */
   int   snap_origins;          /* 1 = origins snap to atoms */
   int   auto_position;         /* do we attempt to centre new molecules */
   int   auto_origin;           /* do we write coordinates on original PDB frame */
   int   unify_models;          /* do we save dockings to a single model file */
   int   problem_type;          /* 1=similarity, 2=docking, 3=multidocking */

   int   mramp;                 /* mesh ramp: 0=potential, 1=charge, etc. */
   int   sramp;                 /* surface ramp: 0=MEP 1=IEL 2=EAL 3=ALPHA */

   Ramp  mesh_ramp[6];          /* the actual colour ramp dfefinitions */
   Ramp  surface_ramp[5];


   double r12_distance;         /* 1:1 docking intermolecular distance */
   double r12_shift;            /* extra shift to force separation */
   double default_distance;     /* default separation for 1:1 docking */
   double rms_thres;            /* for sidechain RMS colouring */
   float  ligand_base_rms;      /* wrt fitting to ligand complex */

   int  make_font;
   char font_name[40];
   char font_weight[40];
   char font_size[40];

} GraphicsState;

/*---------------------------------------------------------------------------*/

/*  Function prototypes... */
/*  ---------------------- */

int hex_battery(int argc, char **argv, char **envp);

int   hex_using_stereo();

void hex_accessibility      (int  n_atoms, int  n_samples, 
                             double probe_radius,
                             Box box, float xyz[], float rad[],
                             byte accessible[]);


void   hex_alloc_shared_mem();
void  *hex_shared_mem      (int id);  /* for access to shared memory */

void   hex_start_browser  ();
void   hex_set_browser0   (char *browser_path);
void   hex_set_browserd1  (char *browser_path);
void   hex_set_browserd2  (char *browser_path);
void   hex_set_browseru1  (char *browser_path);
void   hex_set_browseru2  (char *browser_path);
void   hex_set_browserw1  (char *browser_path);
void   hex_set_browserw2  (char *browser_path);

void   hex_gui_activated  (int have_stereo);
void   hex_update_gui     (int thing, ...);
int    hex_gui_enabled    ();
int    hex_gui_stereo     ();
void   hex_run_job        (int the_job);

int   get_receptor_id   ();
int   get_ligand_id     ();

void hex_open_cmd       ();
void hex_write_cmd      (char *fmt, ...);
void hex_close_cmd      ();
int  hex_query_cmd      ();
void hex_save_cmd       (char *filename);
void hex_run_startup    ();
void hex_run_macro      (char *filename);
void hex_run_cmd        (char *cmd, 
                         char *a1, char *a2, char *a3, char*a4, char *a5);
int  hex_macro_running  ();

void hex_mm_open        (double r_max);
int  hex_mm_idx         (double r);
int  hex_mm_idx2        (double r2);
int  hex_mm_lj          (int  idx, int lj_code, double *e);
int  hex_mm_hb          (int  idx, int proton_code, int lp_code, double *e);
int  hex_mm_es          (int  idx, float q1, float q2, double *e);
void hex_mm_close       ();

void open_clipping_planes      (double size);
void close_clipping_planes     ();
void hide_clipping_planes      ();
int  any_active_planes         ();
int  enabled_clipping_planes   ();
void setup_clipping_planes     (int pln, double size);
void commit_clipping_edits     ();
void cancel_clipping_edits     (double size);
void disable_clipping_planes   ();
void scale_clipping_planes     (int sense);
void apply_clipping_planes     ();
void draw_clipping_planes      (Matrix3D scn_rot, double scale,
                                int current_col, int active_col,
                                int marked_col);
void toggle_current_plane      ();
void toggle_unmarked_planes    ();
void toggle_plane_marking      (int pln);
void swap_clipping_sense       ();
void shift_clipping_planes     (double amount);
void activate_clipping_plane   (Matrix3D scn_rot, int pln, double scale);
void rotate_clipping_planes    (Matrix3D scn_rot, Matrix3D tr);
void translate_clipping_world  (double x, double y, double z);
void translate_clipping_scene  (Matrix3D scn_rot,
                                double x, double y, double z);
void translate_clipping_plane  (Matrix3D scn_rot,
                                double x_trans, double y_trans);


void  hex_light_shininess       (double shininess);
void  hex_light_specular        (double specular);
void  hex_light_diffuse         (double r, double g, double b);
void  hex_light_position        (double x, double y, double z);

void  hex_create_background     (double width, double height, 
                                 double near, double far);
void  hex_draw_background       ();
void  hex_delete_background     ();

void  hex_draw_axes             (int mode, int fore, int back, float length);
void  hex_delete_axes           ();

void hex_print_scene            (char *file, int print_type, int print_eps);

void  get_molecule_rotation     (double *alpha, double *beta, double *gamma);


void  setup_gsurface            (Gsurface *s, double alpha, 
                                 int  n_poles, int  n_nodes,
                                 double *r_pole, double *r_node);
void  colour_gsurface           (Gsurface *s, Ramp r, double alpha,
                                 double *v_pole, double *v_node);
void  recolour_gsurface         (Gsurface *s, Ramp r, double alpha);
void  draw_gsurface             (Gsurface *s, int poly_type, int surf_type, 
                                 int edge_col);
void  free_gsurface             (Gsurface *s);

void  hex_canonical_ylm         (Molecule *molecule, int model_no,
                                 int ellipsoid_order,
                                 int  n_poles, int  n_nodes);

void  hex_canonical_gto         (Molecule *molecule, int model_no,
                                 int ellipsoid_order);

/*---------------------------------------------------------------------------*/

/* operations that can be applied to a molecule */

int  open_the_receptor           (char *filename); /* convenience fns */
int  open_the_ligand             (char *filename);
int  open_the_complex            (char *filename);
void open_the_clone              (Molecule *parent, int from_model, int to);
void set_the_origin              (int mol);
void set_the_interface           (int mol);

int  open_molecule               (Molecule *, int id, int mol_id, char []);
void close_molecule              (Molecule *);
void delete_molecule_ylm         (Molecule *, int model_no);
void delete_molecule_gto         (Molecule *, int model_no);
void delete_molecule_eto         (Molecule *, int model_no);
void delete_cartoon_dlist        (Molecule *);
void delete_atom_dlist           (Molecule *);
void delete_molecule_mesh       (Molecule *);
void delete_molecule_samples     (Molecule *);
void delete_molecule_normals     (Molecule *);
void delete_normal_dlist         (Molecule *);
void delete_molecule_airs        (Molecule *, int which);
void delete_molecule_gnome       (Molecule *, int model_no);
void delete_gnome_dlist          (Molecule *);
void delete_molecule_gsurface    (Molecule *);
void delete_dots_dlist           (Molecule *);
void delete_mesh_dlist           (Molecule *);
void delete_molecule_vcolours    (Molecule *);
void assign_molecule_colours     (Molecule *molecule, Molecule *reference,
                                  Matrix3D t_dock, double rms_thresh);
void assign_complex_colours      (Molecule *reference);
void allocate_molecule_model     (Molecule *molecule);
void allocate_molecule_airs      (Molecule *molecule, int active, char *pat);
void deallocate_molecule_model   (Molecule *molecule);
int  blast_molecule_model        (Molecule *);
void load_molecule_colours       (Molecule *molecule);
void set_molecule_colour         (Molecule *molecule, int col_num);
void load_molecule_colour_range  (Molecule *molecule, char chain,
                                  char *res1, char *res2, char *colour);
void reset_molecule_model        (Molecule *molecule, int model_no);

void delete_msurf                (MeshSurface *msurf);

void draw_molecule               (Molecule *, int solid_style);
void hex_save_pdb               (char *, Docking *, int mol, 
                                 Molecule *receptor, Molecule *ligand, 
                                 double r12_shift, int real_soln, int sym_el, 
                                 int auto_origin);
void hex_save_pdbs              (char *, Docking *, 
                                 Molecule *receptor, Molecule *ligand, 
                                 double r12_shift, int real_soln, int auto_origin);
void hex_save_pdbs_hack         (char *, Docking *, 
                                 Molecule *receptor, Molecule *ligand, 
                                 double r12_shift, int real_soln, int auto_origin);
HexFile *save_model_open         (char *filename);
void save_model_close            (HexFile *unit);
void save_model_write            (HexFile *unit, Docking *docking,
                                  Molecule *receptor, Molecule *ligand, 
                                  double r12_shift, int soln, int auto_origin, int model);
void pipe_molecules              (Docking *, 
                                  Molecule *, Molecule *, 
                                  char [], 
                                  double r12_shift, int, int, Matrix3D t_scene);
int  mol_write_pdb               (HexFile *unit, int box_mode, double r12_shift,
                                  Docking *docking, int mol, int  n1,
                                  Molecule *receptor, Molecule *ligand,
                                  int soln, int sym_el, int save_mode,
                                  Matrix3D t_scene, Box *box);
void dump_molecule               (Molecule *, char []);

int  lookup_molecule_id          (Molecule *, int  cnum, int  rnum, int  anum,
                                  char chain_id[], char residue_name[], 
                                  char residue_id[], char atom_name[]);
void mark_molecule_residues      (Molecule *molecule, int mark);
void mark_molecule_label         (Molecule *molecule, int  atom, int mark);
int  molecule_atom_label         (Molecule *molecule, int  atom);
void clear_molecule_labels       (Molecule *molecule);

/*---------------------------------------------------------------------------*/

void   hex_draw_msurface   (Molecule *mol);

void   hex_draw_gsurface   (Molecule *molecule, int l_max, int p_max,
                            int  n_poles, int  n_nodes,
                            double alpha, double beta, double gamma);

void hex_vertex_atom_colour(Molecule *molecule, int model_no, int hydrogens, 
                            double probe_radius, double z_shift, double ablend,
                            int  nv, Point3D *vpoints, int *vcolours);

void hex_depth_sort      (int  nv, Point3D *pts, int  *order);
void hex_triangle_sort   (int  nt, int  *triangles, Point3D *pts, int  *order);

void   clear_mol_transformations(int);

void hex_set_mol_vol       (Molecule *mol);

int  hex_get_mol_pts       (Molecule *mol, Point3D **xptr);

void hex_rho               (Molecule *molecule, int  n_max, 
                            int docking_density, int dot_density,
                            double skin_factor);

void hex_rho_seg           (Molecule *molecule, int  n_max, int dot_density,
                            Matrix3D t, double r12, double r_max,
                            double *anlm, double *bnlm);

void hex_rho_integrate     (int  n_max, int  n_poles, int  nr, 
                            double r_min, double r_step, 
                            byte *rgrid, byte *tgrid,
                            double *anlm, double *bnlm);

void hex_model_integrate    (Molecule *mol1, int mdl1, Molecule *mol2, int mdl2,
                            int n_max, int n_dars, int ng, double dg, 
                            Point3D origin, byte *sbuf, double *abuf, double *bbuf);

void hex_grid_integrate    (Molecule *mol1, int mdl1, Molecule *mol2, int mdl2,
                            int n_max, int n_dars, int ng, double dg, 
                            Point3D origin, byte *sbuf, double *abuf, double *bbuf);

void hex_skin_integrate    (int  n_max, int  ng, double dg,
                               Point3D origin, byte *sbuf,
                               double *anlm, double *bnlm,
                               double *cnlm, double *dnlm);

int  make_rshells          (double rmax, double emax, int  tmax, 
                            int  *n_cells, float **r_shells, int  **t_shells);

int hex_gnome              (Molecule *molecule, int model_no);

void   ylm_rot_xyz         (Molecule *, int, Matrix3D,
                            double [], double [], double []);

void   ylm_sample_surface  (Molecule *, int order, int  Npoles, int  Nnodes);

void   ylm_sample_data     (int order, int  Npoles, int  Nnodes,
                            double *rp, double *rn, double *ylm_c);

double ylm_sample_area     (Molecule *, int order, int  Nnodes);
double ylm_sample_error    (Molecule *, int order, int  Nnodes);

double hex_shs_value       (HexShd *shd, int  l_type, int  l_max,
                            double theta, double phi);

void   hex_shs_rotate      (HexShd *shd, int  l_type, int  l_max,
                            double alpha, double beta, double gamma,
                            double *alm_rot);

void hex_shs_store         (HexShd *shd, int l_type, int l_max,
                            int n_poles, int n_nodes, double *r_poles);

void hex_shs_integrate     (int l_max, int n_poles, int n_nodes, 
                            double *r_poles, double *alm);

void hex_shs_pole_nodes    (int l_max, double *alm, 
                            int  n_poles, int  n_nodes,
                            double *rp, double *rn, 
                            double *r_min, double *r_max);

Vector3D ylm_normal          (int,double,double,double[]);   /* surface normal*/

Vector3D pp_normal           (double,double,double,double,double);  

void   ylm_integral_test   (int l, int m, int  Nnodes);

void   hex_draw_dots       (Molecule *);

void   hex_gaussian_softness (double softness);
double hex_gaussian_threshold();

void   hex_gaussian_surface(Molecule *, int model_no, int with_hydrogens,
                            int culling, 
                            double apolar_probe, double polar_probe, 
                            double grid_size, 
                            double softness, MeshSurface *msurf);

void   hex_msurface_colour  (Molecule *molecule);

void   hex_surface_steric  (Molecule *, int model_no, double probe_radius);

void   hex_surface_main    (Molecule *, int model_no, double probe_radius, 
                            int flag);

void   hex_surface_softness(double softness);

void   hex_dot_surface     (Molecule *, int model_no);

void hex_charge            (Molecule *molecule, int mdl, int  Npoles);

void   hex_cull_solvent    (Molecule *);

void   hex_toggle_picks    (Molecule *);

void   hex_enable_arglys   (Molecule *, int enable);

void   hex_enable_solvent  (Molecule *, int enable);

void   hex_enable_hetero   (Molecule *, int enable);

void   hex_dots            (int  n_atoms, int  *ids,
                            int  n_vertex, double r_probe,
                            Box box, float *xyz, float *rad,
                            unsigned char *accessible, 
                            int    *n_surface_atoms, 
                            int   **surface_atom_ids, 
                            int    *n_surface_pts,
                            int   **atom_pt_ids,
                            float **atom_surface_pts,
                            float **atom_outer_pts);

void hex_setup_8cell       (Box box, double r_max, double grid_size,
                            int  *ni_grid, int  *nj_grid, int  *nk_grid,
                            Point3D *grid_origin);

void hex_sample_8cell      (int  n_atoms, Point3D *pts, float *rad,
                            double grid_size, int  ni, int nj, int nk, 
                            Point3D origin, float *grid_values);

HexSurface hex_contour     (int culling, int  ni, int  nj, int  nk, 
                            double grid_size, Point3D grid_origin, 
                            float contour_level, float *grid_values);

void  hex_surface_free     (HexSurface *hsurf);

void hex_docker_cluster    (Docking *docking, 
                            Molecule *receptor, Molecule *ligand,
                            int cluster_window, int max_clusters,
                            double rms_threshold, int bumps_threshold);

void hex_restraint         (DockSpec *ds, Docking *docking,
                            Molecule *receptor, Molecule *ligand);

void    hex_native_receptor(Molecule *reference, Molecule *receptor);

float   hex_native_ligand  (Molecule *reference, Molecule *ligand, double r12);

void    hex_native_postdock(Molecule *molecule);

float   hex_native_rms     (Matrix3D t_predock_receptor,
                             Matrix3D t_postdock_receptor,
                             Matrix3D t_predock_ligand,
                             int ligand_model, double r12, Matrix3D t_rec, Matrix3D t_lig);
//                           int ligand_model, double r12, Euler er, Euler el);

void    hex_native_close   ();

Matrix3D tf_fit_ligand_complex (Molecule *reference, Molecule *ligand,
                                double r12);

Matrix3D tf_fit_molecules   (Molecule *fixed, Molecule *moving);

void hex_fit_molecules      (Molecule *fixed, Molecule *moving);

float check_transform_rms   (Matrix3D tl,
                             int  nvr, int  nvl, int  na, double dr,
                             double r12);
void run_check_transform     (int  nvr, int  nvl, int  na, double dr, double r12);

Matrix3D get_docking_transform (int mol, int soln, int for_surface);

void    hex_deviations      (Molecule *reference, Molecule *ligand, 
                             Matrix3D t_dock, double rms_tol);


void hex_draw_hbonds     (int the_mol);
int  hex_get_hbonds      (int the_mol, Point3D **pts1, Point3D **pts2);
int  hex_num_hbonds      (int the_mol);
void hex_create_hbonds   (Docking *docking, 
                          Molecule *receptor, Molecule *ligand, 
                          Molecule *reference, int solution_number);
void hex_delete_hbonds   (int the_mol);

int  hex_molecule_hbonds (Molecule *molecule, int use_protons,
                          Point3D **pts1, Point3D **pts2);
int  hex_complex_hbonds  (Docking *docking,
                          Molecule *receptor, int r_protons, 
                          Molecule *ligand, int l_protons,
                          int solution_number, Point3D **pts1, Point3D **pts2);

void hex_insert_rotamer  (Docking *docking, Molecule *molecule, int soln, 
                          char chain_label, char *residue_id, 
                          char *residue_name, int  rotamer_no); /* test fn */

void hex_create_rotamer  (Docking *docking, Molecule *molecule, int soln, 
                          int  rotamer_no,
                         int  chain, int  residue, Point3D *atom_coords);
int  hex_install_rotamers(Docking *docking, Molecule *molecule, int soln);
void hex_copy_rotamers   (Rotamer *source, Rotamer **target);
void hex_delete_rotamers (Rotamer **rotamers);

void hex_model_setup_one           (Molecule *molecule);
void hex_model_setup_two           (Molecule *receptor, Molecule *ligand,
                                    int auto_pos);

Matrix3D hex_mol_to_world    (Molecule *molecule, int model_no);
Matrix3D hex_mol_to_origin   (Molecule *molecule, int model_no);
Matrix3D hex_origin_to_pdb   (Molecule *molecule, int model_no);
Matrix3D hex_model_to_world  (Molecule *molecule, int model_no);
Matrix3D hex_model_to_origin (Molecule *molecule, int model_no);
Matrix3D hex_model_transform (Molecule *molecule, int model_no);
Point3D          hex_model_centroid  (Molecule *molecule, int model_no);

void  hex_model_spawn              (Molecule *receptor, Molecule *ligand,
                                    Point3D pt1, Point3D pt2);

void  hex_transform_model          (Molecule *molecule, int model_no, 
                                    Matrix3D t);

void  hex_model_make_clone         (Molecule *, int from_model, 
                                    int to_model, Molecule *, int);
void  hex_model_commit_clone       (Molecule *clone);

int   hex_n_models                 (Molecule *molecule);
void  hex_model_atoms              (Molecule *molecule, int model_no, 
                                    int  *a1, int  *a2);
void  hex_docking_model_range      (Molecule *molecule, int *mdl1, int *mdl2);
char *hex_model_name               (Molecule *molecule, int model_no);
char *hex_model_label              (Molecule *molecule, int model_no);
int   hex_model_id                 (Molecule *molecule, int model_no);

void hex_clash_enable              (int flag);
void hex_clash_open                (Molecule *receptor, Molecule *ligand,
                                    int mdl1, int mdl2);
int  hex_clash_check               (double r12,    double beta1, double gamma1,
                                    double alpha2, double beta2, double gamma2);
void hex_clash_close               ();

void post_mramp_gui            (Ramp *ramp);
void post_sramp_gui            (Ramp *ramp);

Ramp hex_get_mramp             (int j_ramp);
void hex_rescale_mramp         (int j_ramp);

Ramp hex_get_sramp             (int j_ramp);
void hex_rescale_sramp         (int j_ramp);

int  hex_get_motion_fill       ();

void display_macro_normals     (int mol);
void activate_macro_normals    ();
void set_macro_chains          (char *);

void update_moving_molecule    (int molecule_no, int model_no);
void update_moving_origin      (int);
void update_molecule_rotation  ();

void hex_setup_scene           ();
void hex_set_scene_distance    (double d12);
void hex_set_scene_centre      (int centre_thing);
void hex_set_scene_motion      (int motion_thing);
void hex_set_user_centre       (Point3D user_centre);
void hex_rotate_scene          (Matrix3D rot);
Vector3D hex_get_scene_shift     ();
Point3D  hex_get_scene_centre    ();
int  hex_get_docking_motion    ();
int  hex_get_centre_thing      ();

void setup_main_projection     (int eye);
void setup_space               ();
void setup_molecules           ();
void setup_test                ();
void setup_logo                ();
void setup_world               ();
void setup_mesh                (int new_density);

void hex_open_logo             ();
void hex_draw_logo             (int width, int height);
void hex_close_logo            ();

void select_solution           (int soln);
void print_solution            ();
void print_scores              ();

void hex_ps_feedback            (char *filename, int print_eps, double *back);
void hex_ps_pixmap              (char *filename, int print_eps, double *back);

void render_main_selection     (int button, double x, double y);
void render_main_scene         (int eye);

void rotate_world              ();
void get_world_scale           (double *xs, double *ys);

int  n_molecules               ();  /* counts all molecules */
int  both_molecules            ();  /* check for exactly receptor+ligand */

void hex_delete_display_lists  ();

void remove_all_lists          ();
void remove_test_lists         ();
void remove_molecule_lists     (int);

void hex_init_state            ();

Ramp   hex_ramp_make           (char *name, double u_min, double u_max, 
                                double mid_frac, Vector3D min_colour, 
                                Vector3D mid_colour, Vector3D max_colour);
void   hex_ramp_reset          (Ramp *ramp);
void   hex_ramp_range          (Ramp *ramp, 
                                double u_min, double u_max, double mid_frac);
void   hex_ramp_colour         (Ramp *ramp,
                                Vector3D min_colour, Vector3D mid_colour, 
                                Vector3D max_colour);
int    hex_ramp_value          (Ramp ramp, double alpha, double value);

void shutdown_main_scene       ();

int  hex_next_frame            ();

char *get_solution_text1       ();
char *get_solution_text2       ();
char *get_frame_rate_text      ();

/*---------------------------------------------------------------------------*/

//ifdef __cplusplus
//}
//endif

#endif  /* hex_molecule_h */

/*---------------------------------------------------------------------------*/
