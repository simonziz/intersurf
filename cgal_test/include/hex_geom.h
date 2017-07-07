/*-----------------------------------------------------------------------------
**  File:       hex_geom.h
**
**  Author:     Graham Kemp 1994
**              Dave Ritchie, 16/04/96  - Added Polar stuff, bugfix
**                                      - renamed from Graham's geom.h
**              Dave Ritchie, 1996-2010 - periodically added more bits...
**
**  Purpose:    #include file for the various "hex" functions.
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 1994-2003 G.J.L. Kemp, D.W. Ritchie, University of Aberdeen
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
**--------------------------------------------------------------------------*/

#ifndef hex_geom_h
#define hex_geom_h
                                                                                
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------*/

#define FALSE   0
#define TRUE   1

#define PI               3.14159265358979323846
#define PIf     ((float) 3.14159265358979323846)
#define PIl              3.1415926535897932384626433832795029L

#define degrees_to_radians(X)   ((X) * PI / 180)
#define radians_to_degrees(X)   ((X) * 180 / PI)

/*----------------------------------------------------------------------------
Type definitions
----------------------------------------------------------------------------*/
typedef double   Distance;
typedef double   Angle;

#ifndef Logical
typedef int   Logical;
#endif

typedef struct {
   double   x, y, z;
} Point3D, Vector3D;

typedef struct {
   float x, y, z;
} Point3Df, Vector3Df;


typedef struct {
   double   r, theta, phi;
} Polar;

typedef struct {
   Angle    alpha, beta, gamma;
} Euler;


typedef struct {
   Point3D   line_point;
   Vector3D  line_vector;
} Line3D;

typedef struct {
   Angle      angle;
   Vector3D   normal;
} Axis, Cone3D;

typedef struct {
   double q0, q1, q2, q3;
} Quat3D;


typedef struct {
   /*--------------------------------------------------------------------
   Equation of plane is
      ax + by + cz + d = 0
   where the vector (a,b,c) is perpendicular to the plane.
   --------------------------------------------------------------------*/
   double   a, b, c, d;
} Plane;


typedef struct {
   Point3D   sphere_centre;
   double   sphere_radius;
} Sphere;

typedef struct {
   Point3D   circle_centre;
   double   circle_radius;
   Vector3D   circle_normal;
} Circle;

typedef struct {
   Point3D    box_min;
   Point3D    box_max;
} Box;


/*----------------------------------------------------------------------------
Function Prototypes of all functions in hex_geom.cpp
----------------------------------------------------------------------------*/

int        isnan_pt       (Point3D pt);
Point3D    pt_bspline     (double t, int segment, Point3D *pts);
Vector3D   v_bspline      (double t, int segment, Point3D *pts);
Point3D    pt_bezier      (double u,int npts, Point3D *pts);
Vector3D   v_bezier       (double u,int npts, Point3D *pts);
Vector3D   v_random       ();
Box        box_pt         (Point3D pt);
Box        box_pts        (int npts, Point3D *pts);
Box        box_boxpt      (Box box, Point3D pt);
Box        box_box        (Box box1, Box box2);
Point3D    pt_box         (Box box);
double     d_box          (Box box);
double     d2_ptpt        (Point3D p1, Point3D p2);
Distance   d_ptpt         (Point3D, Point3D);
double     component_vv   (Vector3D, Vector3D);
double     det_vvv        (Vector3D, Vector3D, Vector3D);
double     triangle_area  (Point3D, Point3D, Point3D);
double     mag_vector     (Vector3D);
double     square_vector  (Vector3D);
Vector3D   op_vector      (Vector3D);
Vector3D   v_magv         (double, Vector3D);
Vector3D   v_scale        (double, Vector3D);
Vector3D   v_lambdav      (double, Vector3D);
Vector3D   v_add          (Vector3D, Vector3D);
Vector3D   v_sub          (Vector3D, Vector3D);
Vector3D   v_neg          (Vector3D);
Vector3D   v_pt           (Point3D);
Vector3D   v_ptpt         (Point3D, Point3D);
Vector3D   v_zero         ();
Vector3D   v_unit         (Vector3D);
Vector3D   v_cross        (Vector3D, Vector3D);
Vector3D   unit_vector    (Vector3D);
Vector3D   zero_vector    ();
Vector3D   cross_product  (Vector3D, Vector3D);
double     dot_product    (Vector3D, Vector3D);
Angle      a_vv           (Vector3D, Vector3D);
Point3D    pt_v           (Vector3D);
Point3D    pt_ptv         (Point3D, Vector3D);
Point3D    pt_add         (Point3D, Point3D);
Point3D    pt_sub         (Point3D, Point3D);
Point3D    pt_midpt       (Point3D, Point3D);
Point3D    pt_average     (int npts, Point3D *pts);
Point3D    pt_stddev      (int npts, Point3D *pts, Point3D mu);
double     stddev_pts     (int npts, Point3D *pts, Point3D mu);
Point3D    pt_scale       (double, Point3D);
Point3D    pt_zero        ();
Point3D    pt_nf_pl       (Plane);
Point3D    pt_n_ptpl      (Point3D, Plane);
Point3D    pt_i_lpl       (Line3D, Plane);
Point3D    pt_n_ptl       (Point3D, Line3D);
Distance   d_ptl          (Point3D, Line3D);
Line3D     l_j_ptpt       (Point3D, Point3D);
Line3D     l_pt           (Point3D);
Plane      pl_ptptpt      (Point3D, Point3D, Point3D);
Plane      pl_n_ptv       (Point3D, Vector3D);
Vector3D   v_n_pl         (Plane);
Vector3D   n_ptptpt       (Point3D p1, Point3D p2, Point3D p3);
Line3D     l_i_plpl       (Plane, Plane);
Logical    i_ss           (Sphere, Sphere);
Plane      pl_i_ss        (Sphere, Sphere);
Logical    i_sss          (Sphere, Sphere, Sphere);
Circle     c_i_ss         (Sphere, Sphere);
Plane      pl_c           (Circle);
Logical    c_in_s         (Circle, Sphere);
Logical    i_sc           (Sphere, Circle);
Logical    i_sc_ptpt      (Sphere, Circle, Point3D *, Point3D *);
Angle      a_ptptpt       (Point3D, Point3D, Point3D);
Angle      a_ptptptpt     (Point3D, Point3D, Point3D, Point3D);
Distance   d_pp_pp        (Polar, Polar);
Point3D    pt_pp          (Polar);
Polar      pp_pt          (Point3D);
Point3D    pt_floatv      (float *xyz);
Vector3D   v_floatv       (float *xyz);
Point3D    pt_doubles     (double *xyz);
Vector3D   v_doubles      (double *xyz);
void       floatv_v       (Vector3D v, float *xyz);
void       floatv_pt      (Point3D p, float *xyz);
Polar      pp_make        (double r, double theta, double phi);
Point3D    pt_make        (double x, double y, double z);
Vector3D   v_make         (double x, double y, double z);
double     deg2rad        (double deg);
double     rad2deg        (double rad);
Euler      e_make         (Angle alpha, Angle beta, Angle gamma);

Point3Df   ptf_pt         (Point3D pt);
Point3D    pt_ptf         (Point3Df ptf);

int        vv_cones       (Cone3D cone1, Cone3D cone2, Vector3D *v_soln);

/*---------------------------------------------------------------------------*/
                                                                                
#ifdef __cplusplus
}
#endif

#endif /* hex_geom_h */
                                                                                
/*---------------------------------------------------------------------------*/
