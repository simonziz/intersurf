/*-----------------------------------------------------------------------------
**  File:       hex_transform.h
**
**  Author:     Graham Kemp 1994
**              Dave Ritchie, 16/04/96 - Added Polar stuff, bugfixes
**                                     - renamed from Graham's transform.h
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
** --------------------------------------------------------------------------*/

#ifndef hex_transform_h
#define hex_transform_h

#include "hex_geom.h"

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------------
Paul, R.P. (1981) "Robot manipulators", The MIT Press.

Assuming that we will not be scaling or skewing objects, transformations will
have the form:

	(r11	r12	r13	a)
	(r21	r22	r23	b)
	(r31	r32	r33	c)
	(0	0	0	1)

Since the last row is always [0, 0, 0, 1], we don't need to store it.
----------------------------------------------------------------------------*/
typedef struct {

	double	element[3][4];

} Matrix3D, Transformation;

/*----------------------------------------------------------------------------
Function Prototypes of all functions in transform.c
----------------------------------------------------------------------------*/

void           print_tf        (char *msg, Matrix3D m);
void           print_transformation(Matrix3D);
void           tr_tf           (Matrix3D m, Vector3D *v, Matrix3D *r);
void           rt_tf           (Matrix3D m, Matrix3D *r, Vector3D *v);
Matrix3D       tf_tftf         (Matrix3D, Matrix3D);

Matrix3D       tf2             (Matrix3D, Matrix3D);
Matrix3D       tf3             (Matrix3D, Matrix3D, Matrix3D);
Matrix3D       tf4             (Matrix3D, Matrix3D, Matrix3D, Matrix3D);
Matrix3D       tf5             (Matrix3D, Matrix3D, Matrix3D, Matrix3D,
                                Matrix3D);
Matrix3D       tf6             (Matrix3D, Matrix3D, Matrix3D, Matrix3D,
                                Matrix3D, Matrix3D);
Matrix3D       tf7             (Matrix3D, Matrix3D, Matrix3D, Matrix3D,
                                Matrix3D, Matrix3D, Matrix3D);
Matrix3D       tf8             (Matrix3D, Matrix3D, Matrix3D, Matrix3D,
                                Matrix3D, Matrix3D, Matrix3D, Matrix3D);

Matrix3D       tf_translate    (Vector3D v);
Matrix3D       tf_scale        (double s);
Matrix3D       tf_diagonal     (double a, double b, double c);
Matrix3D       tf_rotx         (Angle);
Matrix3D       tf_roty         (Angle);
Matrix3D       tf_rotz         (Angle);
Matrix3D       tf_rotv         (Angle, Vector3D);
Matrix3D       tf_euler        (Euler);
Matrix3D       tf_v            (Vector3D);
Matrix3D       tf_rtf          (Matrix3D);
Matrix3D       tf_inverse      (Matrix3D);
Matrix3D       tf_transpose    (Matrix3D);
Matrix3D       tf_origin       (Point3D);
Matrix3D       tf_identity     ();
Matrix3D       tf_ptpta        (Point3D, Point3D, Angle);
Matrix3D       tf_origin_negz  (Point3D, Point3D);
Matrix3D       tf_origin_posz  (Point3D, Point3D);
Matrix3D       tf_origin_negz_posx(Point3D, Point3D, Point3D);
Matrix3D       tf_origin_posz_posx(Point3D, Point3D, Point3D);
Matrix3D       tf_origin_poszv(Point3D, Vector3D);
Matrix3D       tf_origin_poszv_posxv(Point3D, Vector3D, Vector3D);
Matrix3D       tf_copy         (Matrix3D);
Matrix3D       tf_3_vectors    (Vector3D, Vector3D, Vector3D);
Matrix3D       tf_3_points     (Point3D, Point3D, Point3D);
Matrix3D       tf_1_point      (Point3D);
Matrix3D       tf_mat44        (double *mat44);
Matrix3D       tf_mat33        (double *mat33);
Matrix3D       tf_random       ();
Vector3D       v_tfv           (Matrix3D, Vector3D);
Point3D        pt_tfpt         (Matrix3D, Point3D);
Euler          euler_tf        (Matrix3D);
double         euler_beta      (double beta);
double         trace_tf        (Matrix3D);
double         gd_tftf         (Matrix3D, Matrix3D);
void           rpzr_tf         (Matrix3D, 
                                Matrix3D *tr1, double *tpz, Matrix3D *tr2);
void           rnzr_tf         (Matrix3D, 
                                Matrix3D *tr1, double *tnz, Matrix3D *tr2);
Vector3D       v_tf            (Matrix3D);
void           mat44_tf        (Matrix3D t, double *mat44);
void           mat33_tf        (Matrix3D t, double *mat33);
Axis           axis_tf         (Matrix3D);
Matrix3D       tf_axis         (Axis);
Matrix3D       tf_axis_trig    (Axis);
double         det_tf          (Matrix3D);

Quat3D         q_tf            (Matrix3D t);
Quat3D         q_make          (Quat3D q, Point3D pt);
Quat3D         q_axis          (Axis ax);
Quat3D         q_euler         (Euler e);
Quat3D         q_inv           (Quat3D q);
Quat3D         q_qq            (Quat3D q1, Quat3D q2);
Quat3D         q_q2            (Quat3D q1, Quat3D q2);
Quat3D         q_q3            (Quat3D q1, Quat3D q2, Quat3D q3);
Quat3D         q_q4            (Quat3D q1, Quat3D q2, Quat3D q3, Quat3D q4);
Point3D        pt_qpt          (Quat3D q, Point3D pt);
Matrix3D       tf_q            (Quat3D q);
Matrix3D       tf_q_axis       (Quat3D q);
Axis           axis_q          (Quat3D q);
double         dot_qq          (Quat3D r, Quat3D s);
double         a_qq            (Quat3D r, Quat3D s);
double         d_qq            (Quat3D r, Quat3D s);
double         d2_qq           (Quat3D r, Quat3D s);
double         gd_qq           (Quat3D r, Quat3D s);
double         ge_qq           (Quat3D r, Quat3D s);
double         gf_qq           (Quat3D r, Quat3D s);


/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif /* hex_transform_h */

/*---------------------------------------------------------------------------*/
