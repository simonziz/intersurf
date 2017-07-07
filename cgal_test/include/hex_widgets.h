/******************************************************************************
**
**  File:    hex_widgets.h
**
**  Purpose: C/C++ interface to the FLTK toolkit classes & functions.
**
**  Ryan Cairns,  2002/03   First implementation.
**  Dave Ritchie, 20/04/03  Added new widget sub-classes & tidied up.
**                          Renamed fltk_* -> gui_* for consistency with
**                          other GUI components, and to avoid possible
**                          namespace clashes with the FLTK toolkit itself.
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 2000-2003 D.W. Ritchie, University of Aberdeen.
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

#ifndef hex_widgets_h
#define hex_widgets_h

#include "hex_types.h"

/* this is what the new fltk 2.0.0 looks like */

/*
include <fltk/run.h>
include <fltk/Window.h>
include <fltk/GlWindow.h>
include <fltk/Button.h>
include <fltk/HorizontalSlider.h>
include <fltk/ToggleButton.h>
include <fltk/HelpDialog.h>
*/

#include <FL/Fl.H>
#include <FL/gl.h>
#include <FL/Fl_Window.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Menu_Item.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Menu_Button.H>
#include <FL/Fl_File_Browser.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Bitmap.H>
#include <FL/Fl_Progress.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Text_Buffer.H>
#include <FL/Fl_Menu_Window.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Pack.H>
#include <FL/Fl_JPEG_Image.H>
#include <FL/Fl_BMP_Image.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Color_Chooser.H>
#include <FL/fl_show_colormap.H>
#include <FL/Fl_Pixmap.H>
#include <FL/Fl_Help_Dialog.H>
#include <FL/Fl_Preferences.H>
#include <FL/Fl_Tile.H>
#include <FL/Fl_File_Browser.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_File_Input.H>
#include <FL/Fl_Return_Button.H>
#include <FL/Fl_Shared_Image.H>
#if defined (hex_linux) || (hex_linux64)
#   include  <FL/x.H>
#endif
#include <FL/fl_ask.H>
#include <FL/filename.H>

//include <FL/glut.H>


/*---------------------------------------------------------------------------*/

/* NB. on Windows, "double-buffered" windows prevent window damage when 
   moving windows around - this doesn't seem to be necessary on Unix !? 
*/

#if defined(hex_win32)
#define gui_window Fl_Double_Window
#else
#define gui_window Fl_Window
#endif

/*---------------------------------------------------------------------------*/

/* define short-hand versions for various FLTK names ... I really don't like
   mixed-case names; this also helps make the main-line code look simpler ...
*/

typedef Fl_Group            gui_group;
typedef Fl_Color            gui_colour;
typedef Fl_Widget           gui_widget;
typedef struct Fl_Menu_Item gui_menu;

typedef void (*gui_cb)  (Fl_Widget *w, void *ptr);
typedef void (*gui_fsb) (rchar *action, rchar *filename);

#define i_toggle(W)  ((int)    (((gui_toggle  *) W)->value()))
#define i_choice(W)  ((int)    (((gui_choice  *) W)->value()))
#define i_slider(W)  ((int)    (((gui_slider  *) W)->value()))
#define d_slider(W)  ((double) (((gui_slider  *) W)->value()))
#define c_input(W)   ((char *) (((gui_input   *) W)->value()))
#define i_browser(W) ((int)    (((gui_browser *) W)->value()))

#define   abs(a)   (((a) > 0) ? (a) : -(a))

#define HEX_GUI_TOGGLE      1
#define HEX_GUI_CHOICE      2
#define HEX_GUI_SLIDER      3
#define MAX_PATHNAME      257  /* file name length in bytes */
#define GUI_MAX_WINDOWS    25  /* max no. of window IDs */
#define GUI_MAX_DRAGS      50  /* max no. of "solid-motion" slider IDs */
#define GUI_MAX_BLOCKS     50  /* max no. of "blockable" widgets */
#define MIN_ZOOM_VALUE    -80
#define MAX_ZOOM_VALUE    120
#define MAX_R12_VALUE     100   /* Angstrom units */
#define MIN_NORMAL_VALUE -100
#define MAX_NORMAL_VALUE  100

#define GUI_MESSAGE_SIZE   14
#define GUI_LABEL_SIZE     13
#define GUI_TEXT_SIZE      13
#define GUI_MENU_SIZE      13

// define GUI_MESSAGE_SIZE   12
// define GUI_LABEL_SIZE     12
// define GUI_TEXT_SIZE      12
// define GUI_MENU_SIZE      12

/*---------------------------------------------------------------------------*/

/* sub-class some of the FLTK classes to give Hex-flavoured widgets */

/*---------------------------------------------------------------------------*/

/* main GUI window constructor (sub-class of FLTK's Fl_Window) */

class gui_win : public gui_window {
public:
   gui_win(rchar *txt, int x, int y, int w, int h);
   void dismiss();
   int b();
};

/*---------------------------------------------------------------------------*/

/* graphics window constructor (sub-class of FLTK's Fl_Gl_Window) */
                                                                                
class gui_ogl : public Fl_Gl_Window {
public:
   gui_ogl(int x, int y, int w, int h);
   void draw();
   int handle(int event);
};

/*---------------------------------------------------------------------------*/

/* main GUI */
                                                                                
class gui_display : public gui_win {
public:
   gui_display();
   void dismiss();
// int handle(int event);
   int w();
   int h();
   int b();
private:
   gui_win *hex_main_window;  /* SCRAP THIS */
};

/*---------------------------------------------------------------------------*/

class gui_panel : public gui_win {
public:
   gui_panel(rchar *txt, int w, int h);
   gui_panel(rchar *txt, int x, int y, int w, int h);
   void dismiss();
   int dx();
   int dy();
   int xpos();
   int ypos();
   void xpos(int x_new);
   void ypos(int y_new);
   void large();
   void big();
   void normal();
   void little();
   void tiny();
};

/*---------------------------------------------------------------------------*/

class gui_toggle : public Fl_Button {
private:
   void init(int val, gui_cb cb);
public:
  gui_toggle(int x, int y, rchar *txt, int val, gui_cb cb);
  gui_toggle(rchar *txt, int val, gui_cb cb);
  int w();
  int h();
  int dx();
  int dy();
};

/*---------------------------------------------------------------------------*/

class gui_button : public Fl_Button {
public:
   gui_button(int x, int y, int w, int h, rchar *name, 
              rchar *tip, Fl_Image *image, gui_cb cb);
};

/*---------------------------------------------------------------------------*/

class gui_action : public Fl_Button {
private:
   void init(gui_cb cb);
public:
  gui_action(int x, int y, int w, int h, rchar *txt, gui_cb cb);
  gui_action(int x, int y, rchar *txt, gui_cb cb);
  gui_action(rchar *txt, gui_cb cb);
//  int w();
//  int h();
  int dx();
  int dy();
};

/*---------------------------------------------------------------------------*/

class gui_slider : public Fl_Value_Slider {
private:
   void init(double val, double lo, double hi, int digits, gui_cb cb);
public:
  gui_slider(int x, int y, rchar *txt, double val, double lo, double hi,
              int digits, gui_cb cb);
  gui_slider(rchar *txt, double val, double lo, double hi,
              int digits, gui_cb cb);
  int w();
  int h();
  int dx();
  int dy();
};

/*---------------------------------------------------------------------------*/

class gui_choice : public Fl_Choice {
private:
   void init(int val, gui_menu *items, gui_cb cb);
public:
  gui_choice(int x, int y, rchar *txt, int val, gui_menu *items,
             gui_cb cb);
  gui_choice(rchar *txt, int val, gui_menu *items, gui_cb cb);
  int w();
  int h();
  int dx();
  int dy();
};

/*---------------------------------------------------------------------------*/

class gui_input : public Fl_Input {
private:
   void init(gui_cb cb);
public:
  gui_input(int x, int y, int w, int h, rchar *txt, gui_cb cb);
  gui_input(int x, int y, rchar *txt, gui_cb cb);
  gui_input(rchar *txt, gui_cb cb);
  int w();
  int h();
  int dx();
  int dy();
};

/*---------------------------------------------------------------------------*/

class gui_box : public Fl_Box {
public:
   gui_box(int x, int y, int w, int h, rchar *txt);
};

/*---------------------------------------------------------------------------*/

class gui_colbox : public Fl_Widget {
  void draw();
public:
  gui_colbox(int x, int y, int w, int h)
  :Fl_Widget(x, y, w, h) { box(FL_ENGRAVED_FRAME);}
  void rgb(unsigned char *);
private:
  uchar rgb_cols[4];
};

/*---------------------------------------------------------------------------*/

class gui_gradient : public gui_box {
  void draw();
public:
  gui_gradient(int x, int y, int w, int h, rchar *l, Ramp r)
    : gui_box(x, y, w, h, l) { set_ramp(r); box(FL_DOWN_FRAME); }
  void set_ramp(Ramp r);
private:
   Ramp ramp;
};

/*---------------------------------------------------------------------------*/

class gui_browser : public Fl_File_Browser {
public:
   gui_browser(int x, int y, int w, int h, gui_cb cb);
};
                                                                                
/*---------------------------------------------------------------------------*/

class gui_vscale : public Fl_Slider {
public:
  gui_vscale(int x, int y, int w, int h, rchar *txt, 
             double val, double lo, double hi, gui_cb cb);
};

/*---------------------------------------------------------------------------*/

/* these GUI functions are compiled with C++, and are NOT for main-line Hex */

int    gui_open_visual      (int graphics);
void   gui_plain_visual     ();
void   gui_process_events_thread();
int    gui_process_events   (double time_to_wait);
void   gui_register_drag    (Fl_Widget *w);
void   gui_register_window  (gui_win * w);
void   gui_register_block   (Fl_Widget *w);
void   gui_post_rotation    (int mode, int button, int x,int y);
void   gui_post_translation (int,int);
void   gui_post_pick        (int button, int x, int y);
void   gui_set_cursor       (int button);
void   gui_close_panels     ();
void   gui_simple_colour_chooser(int target);
void   gui_fancy_colour_chooser(rchar *txt, int target);
int    gui_key_down         (int key);
int    gui_key_up           (int key);
int    gui_full_screen      ();
void   gui_close_windows    ();
char  *gui_icon             ();
char  *gui_message_identifier();
int    gui_ogl_handler      (int event);
void   gui_ogl_resize       ();
void   gui_ogl_redraw       ();
void   gui_toggle_gui       ();
void   gui_menu_block       (int block);
void   gui_create_fsb(rchar *action, char *dir, char *pat, gui_fsb cb);
void   gui_create_save();
void   gui_create_buttons();
void   gui_create_menu();
void   gui_create_molecule();
void   gui_create_hetero();
void   gui_create_similarity();
void   gui_create_docking();
void   gui_create_matching();
void   gui_create_orientation();
void   gui_create_test();
void   gui_create_cluster();
void   gui_create_model();
void   gui_create_cartoon();
void   gui_create_solid();
void   gui_create_surfaces();
void   gui_create_movie();
void   gui_create_dot();
void   gui_create_restraint();
void   gui_create_projection();
void   gui_create_macrosampling();
void   gui_create_harmonics();
void   gui_create_lighting();
void   gui_create_separation();
void   gui_create_normal();
void   gui_create_mramp();
void   gui_create_sramp();
void   gui_create_about();
void   gui_create_hotkeys();
void   gui_create_progress();

/* NB. FLTK colours, Hex colour IDs, and Hex packed ints are ALL DIFFERENT */

gui_colour gui_rgba_b2c (uchar *rgba);               /* RGB bytes to FLTK */ 
void       gui_rgba_c2b (gui_colour c, uchar *rgba); /* FLTK to RGB bytes */
gui_colour gui_rgb_s2c  (rchar *colname);            /* Hex colour to FLTK */
int        gui_rgb_c2n  (gui_colour c);             /* FLTK to Hex colour ID */
int        gui_rgb_c2i  (gui_colour c);             /* FLTK to Hex packed int */
void       gui_rgb_b2d  (uchar *rgb, double *r, double *g, double *b);
void       gui_rgb_d2b  (double r, double g, double b, uchar *rgb);

#endif /* hex_widgets_h */

/*---------------------------------------------------------------------------*/
