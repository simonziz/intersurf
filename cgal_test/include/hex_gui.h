/*-----------------------------------------------------------------------------
**  File:       hex_gui.h
**
**  Author:     Dave Ritchie, 04/07/96
**
**  Purpose:    To define constants & function prototypes for the interface
**              between the Motif GUI and the main application code...
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 1996-2003 D.W. Ritchie, University of Aberdeen.
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

#ifndef hex_gui_h
#define hex_gui_h

/*---------------------------------------------------------------------------*/
//ifdef __cplusplus
//extern "C" {
//endif

/*---------------------------------------------------------------------------*/

/* the main GUI entry point, called from hex_main() when in graphics mode */

void hex_open_gui(char *prog_name, int argc, char **argv);

/* the secondary GUI entry point, called by main() to make a message window */

void gui_create_messages(char *host_name, int message_port);

/* the function that execs the Hex for message mode */

void gui_spawn_messages(char *prog_name, char *host_name, int message_port);

/* siwtch off message redirection if pipe fails (happens on trogon?) */

void gui_close_messages();

/* the main "down-call" from the command processor to Hex */

void hex_update_state(int thing, ...);

/* these called as "up-calls" from Hex to the GUI */

void     post_hex_job              (int);
void     post_animation_mode       (int);
void     post_frame_rate           (double);
void     post_matching_solution    (int);
void     post_docking_solution     (int);
void     post_docking_max          (int);
void     post_matching_max         (int);
void     post_moving_molecule      (int thing);
void     post_centre_thing         (int thing);
void     post_molecule_rotation    (double alpha, double beta, double gamma);
void     post_display_model        (int id, int model);
void     post_docking_model        (int id, int model);
void     post_normal_display       (int);
void     post_normal_coverage      (int);
void     post_normal_auto          (int);
void     post_redraw               ();         /* ask GUI to re-draw scene soon */
void     post_animation            (int mode); /* tell GUI about animation */
void     post_resize_window        (int width, int height);
void     post_distance_value       (double x);
void     post_shift_value          (double x);
void     post_block_mode           (int block_mode);
void     post_drag_mode            (int solid_motion);
void     post_windows_mouse        (int solid_motion);
void     post_progress_enable      (int enable);
void     post_progress_popup       (int enable);
void     post_interface_atom       (int mol, int status);
void     post_origin_atom          (int mol, int status);
void     post_auto_position        (int status);
void     post_auto_origin          (int status);

void     post_receptor_active      (char *ptr);
void     post_receptor_passive     (char *ptr);
void     post_ligand_active        (char *ptr);
void     post_ligand_passive       (char *ptr);

void     post_receptor_origin      (char *ptr);
void     post_receptor_interface   (char *ptr);
void     post_ligand_origin        (char *ptr);
void     post_ligand_interface     (char *ptr);

void     post_mramp_gui            (Ramp *ramp);
void     post_sramp_gui            (Ramp *ramp);

/* these are C "down-calls" from the GUI to Hex */

void     gui_exit();
void     gui_pipe();
void     gui_abort(char *reason);
int      gui_progress(int phase, double progress);
void     gui_update(int thing, va_list ap);
double   gui_frame_rate();
int      gui_write_font(char *s);
void     gui_set_screenmode(int mode);
void     gui_toggle_screenmode();
void     gui_set_bordermode(int mode);
void     gui_toggle_bordermode();
void     gui_set_projection(int mode);
void     gui_set_background(int mode);
void     gui_set_parallax(int mode);
int      gui_get_beta(int nbeta); /* default angular increments */
void     gui_raise_messages();
int      gui_redirect_messages();
void     gui_trigger_expose(void *user_data);
void     gui_send_expose();
void     gui_sticky_fsb(int sticky);
void     gui_schedule_expose(int  micro_seconds);
void     gui_request_stereo(int stereo);
int      gui_have_stereo();
void     gui_set_stereo();
void     gui_activate_normal(int);
void     gui_set_viewport(int width, double near, double far);

void     gui_pointer_cursor();
void     gui_edit_cursor();
void     gui_origin_cursor();
void     gui_reset_cursor();
void     gui_busy_cursor(int);

/* should move and/or rename these ... */

char    *get_model_name(int id, int model_num);
int      get_num_models(int id);

void     toggle_cartoon_display();
void     toggle_dot_display();
void     toggle_axis_display();
void     toggle_harmonic_display();
void     toggle_sidechain_display();
void     toggle_model_display();
void     toggle_surface_display();

void     clear_current_plane();
void     draw_main_scene(int eye);
void     finish_main_scene();

void     hex_final_setup(int batch, char *program_file);
int      hex_get_pointer_mode();
void     hex_register_wakeup(hex_wakeup_cb);
void     hex_update_callback(hex_update_cb);

/*---------------------------------------------------------------------------*/
//ifdef __cplusplus
//}
//endif

#endif /* hex_gui_h */
/*---------------------------------------------------------------------------*/
