/*-----------------------------------------------------------------------------
**  File:       hex_gl.h
**
**  Author:     Dave Ritchie, 18/03/98
**
**  Purpose:    to provide OpenGL declarations...
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 1998-2003 D.W. Ritchie, University of Aberdeen.
**  Copyright (C) 2010 D.W. Ritchie, INRIA.
**
**  This software was developed at the University of Aberdeen 1996-2000, by
**  Dave Ritchie, as part of a project funded by the BBSRC.
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
---------------------------------------------------------------------------*/

#ifndef hex_gl_h
#define hex_gl_h

#if defined (hex_gui)
#if defined (hex_win32)
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#ifndef GL_RESCALE_NORMAL
#define GL_RESCALE_NORMAL GL_RESCALE_NORMAL_EXT
#endif
#else  /* not hex_win32 */
#if defined (hex_darwin) || defined (hex_darwin64)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <OpenGL/glext.h>
#else  /* unix-style */
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#endif
#endif /* hex_win32 */
#else  /* not hex_gui */
#include "hex_gl_dummy.h"
#endif /* hex_gui */

/*---------------------------------------------------------------------------*/

#define  HEX_GL_LINE_WIDTH  0x1001
#define  HEX_GL_POINT_SIZE  0x1002
#define  HEX_GL_COLOUR      0x1003
#define  HEX_GL_STRING      0x1004

//ifdef __cplusplus
// extern "C" {
//endif

void   hex_print         (char *filename,
                          int print_eps, double *background,
                          int width, int height,
                          GLenum mode,
                          GLint bufsiz, GLfloat *buffer);

int    hex_gl_make_font  (char *name, char *weight, char *size);
int    hex_gl_write_font (char *string);
int    hex_gl_delete_font();

GLint  hex_gl_set_mode   (GLint mode);
GLint  hex_gl_get_mode   ();

void   hex_gl_blend      (int on);

void   hex_gl_linewidth  (GLfloat width);
void   hex_gl_pointsize  (GLfloat size);
void   hex_gl_colour3fv  (GLfloat *colours);
void   hex_gl_colour3f   (GLfloat red, GLfloat green, GLfloat blue);
void   hex_gl_string     (char *string);

void   hex_gl_mult_tf    (Matrix3D t);



void   hex_glVertex      (Point3D v);
void   hex_glNormal      (Vector3D n);
void   hex_glColour      (int rgba);

void   hex_draw_line     (Point3D p1, Point3D p2);
void   hex_draw_cross    (double e, Point3D pt);
void   hex_draw_box      (Box box);
void   hex_draw_circle_xy(int npts, double r);
void   hex_draw_ellipsoid(int npts, int col1, int col2, int col3,
                          double a, double b, double c);


void   hex_create_primitives      (int resolution);
void   hex_delete_primitives      ();

void   prim_create_unit_mesh      (int mesh_resolution);
void   prim_create_unit_sphere    (int resolution);
void   prim_create_unit_cylinder  (int resolution);
void   prim_create_unit_cone      (int resolution);

void   prim_draw_unit_mesh        ();
void   prim_draw_unit_sphere      ();
void   prim_draw_unit_cylinder    ();
void   prim_draw_unit_cone        ();

void   prim_draw_cylinder         (double r, double h);
void   prim_draw_cone             (double r, double h);

void   prim_make_xy_circle        (int npts, double radius, Point3D *pts,
                                                            Vector3D *nrm);
void   prim_draw_xy_circle        (int npts, int cc, Point3D *pts, Vector3D *nrm);
void   prim_draw_xy_cylinder      (int npts, Point3D *pts_bot, Vector3D *nrm_bot,
                                             Point3D *pts_top, Vector3D *nrm_top);
void   prim_draw_xy_strip         (Point3D *pts_bot, Vector3D *nrm_bot,
                                   Point3D *pts_top, Vector3D *nrm_top,
                                   float *face_col, float *side_col);
void   prim_make_xy_rect          (double dx, double dy, Point3D *pts,
                                                         Vector3D *nrm);
void   prim_draw_xy_rect          (int cc, Point3D *pts, Vector3D *nrm);

/*---------------------------------------------------------------------------*/

//ifdef __cplusplus
//}
//#endif

#endif /* hex_gl_h */
/*---------------------------------------------------------------------------*/
