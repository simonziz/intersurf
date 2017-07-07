/*-----------------------------------------------------------------------------
**  File:       hex_cmd.h
**
**  Author:     Dave Ritchie, 04/07/96
**
**  Purpose:    define constants for Hex's "command interpreter"
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

#ifndef hex_cmd_h
#define hex_cmd_h

#define HARMONIC_ORDER                  1
#define HARMONIC_ANGLE                  2 
#define HARMONIC_DISPLAY                3
#define PROBE_SIZE                      4
#define ZOOM_FACTOR                     5
#define NO_LOGO                         6
#define SURFACE_DISPLAY                 7
#define SURFACE_STYLE                   8
#define MOTION_FILL                     9
#define HEX_MESH                       10 
#define SELECT_HOME_VIEW               11
#define FLUSH_UPDATES                  12
#define BUNDLE_UPDATES                 13
#define PROJECTION_TYPE                14 
#define ORTHO_WIDTH                    15
#define Z_NEAR_FAR                     16
#define OPEN_RECEPTOR                  17
#define CLOSE_RECEPTOR                 18
#define OPEN_LIGAND                    19
#define CLOSE_LIGAND                   20
#define OPEN_TEST                      21
#define WINDOW_RESIZE                  22
#define CURSOR_ACTIVE                  23
#define CURSOR_ROTATION                24
#define CURSOR_TRANSLATION             25
#define ACTIVATE_DOCKING               26
#define DOCKING_ALPHA_SAMPLES          27
#define MOLECULE_STYLE                 28  /* deprecated */
#define MATCH_LIGAND_SAMPLES           29
#define DOCKING_LIGAND_SAMPLES         30
#define DOCKING_RECEPTOR_SAMPLES       31
#define HARMONIC_INTEGRAL              32
#define HARMONIC_RECURSION             33
#define SAVE_RECEPTOR                  34
#define SAVE_LIGAND                    35
#define MODEL_DISPLAY                  36
#define TOGGLE_FOG                     37
#define FONT_NAME                      38
#define FONT_WEIGHT                    39
#define FONT_SIZE                      40
#define BACKGROUND_COLOUR              41
#define DISPLAY_SIDECHAIN              42
#define MOLECULE_BONDWIDTH             43
#define MOLECULE_RESIDUE_ID            44
#define MOLECULE_ATOMFILL              45
#define MOLECULE_ATOMRES               46
#define SURFACE_DATA                   47
#define VIEW_MATCHING_SOLUTION         48
#define SEARCH_METHOD                  49
#define TOGGLE_STENCIL                 50
#define FOREGROUND_COLOUR              51
#define SOLID_ORDER                    52

#define MOLECULE_SEPARATION            54
#define CENTRE_MOLECULE                55
#define SET_COLOUR_FILE                56
#define SET_ROTATION                   57
#define ELLIPSOID_ARCS                 58  /* obsolete */
#define ELLIPSOID_AXES                 59
#define NORMALISE_MOLECULE             60
#define BUTTON_ACTIVE                  61
#define SET_HOME_VIEW                  62
#define COMMIT_VIEW                    63  /* obsolete */
#define TRIG_CACHE                     64
#define DISPLAY_AXES                   65
#define ELLIPSOID_ORIGIN               66
#define MOLECULAR_ORIGIN               67
#define MOVING_THING                   68
#define MOLECULAR_AXIS                 69
#define BUTTON_PICK                    70
#define OPEN_DATA                      71 
#define DISPLAY_CONTACT                72  /* not used */
#define OPEN_CUBE                      73
#define CLOSE_ALL                      74
#define SAVE_BOTH                      75
#define DISPLAY_SOLVENT                76
#define POINTER_MODE                   77
#define DOT_DISPLAY                    78
#define DOT_DENSITY                    79
#define DOT_STYLE                      80
#define DOT_PROBE                      81
#define DOT_COLOUR_NAME                82
#define DOT_SIZE                       83
#define HAMMER_POINTER_MODE            84
#define DOT_SCOPE                      85
#define POLY_COLOUR_NAME               86
#define SURFACE_DOTS                   87  /* obsolete */
#define TEST_ROTATION                  88
#define DOCKING_MOTION                 89
#define POLY_SCOPE                     90
#define LOCK_HOME_VIEW                 91
#define DOCKING_MAIN_SCAN              92
#define DOCKING_MAIN_SEARCH            93
#define DISPLAY_SYMMETRY               94
#define MOLECULE_SCOPE                 95
#define COMMIT_EDITS                   96
#define VIEW_DOCKING_SOLUTION          97
#define WINDOW_EXPOSE                  98
#define PRINT_SIZE                     99
#define MOLECULE_COLOUR_ID            100
#define MATCHING_MAIN_SEARCH          101
#define MATCHING_MAIN_SCAN            102
#define MATCHING_ELECTROSTATICS       103
#define MATCHING_LIGAND_SAMPLES       104
#define MATCHING_RECEPTOR_SAMPLES     105
#define MATCHING_ALPHA_SAMPLES        106
#define MATCHING_GRID_SIZE            107
#define MATCHING_R12_RANGE            108
#define MATCHING_R12_STEP             109
#define MATCHING_SEARCH_MODE          110

#define DOT_COLOUR_ID                 117
#define POLY_COLOUR_ID                118
#define DOT_SKIN                      119
#define PRINT_INFO                    120
#define SOLVENT_COLOUR_ID             121
#define SOLVENT_SCOPE                 122
#define ENABLE_SOLVENT                123
#define CULL_SOLVENT                  124
#define TOGGLE_PICKS                  125
#define MOVE_NORMAL                   126
#define KEYBOARD_SPACE                127
#define KEYBOARD_BACKSPACE            128
#define KEYBOARD_PLUS                 129
#define KEYBOARD_MINUS                130
#define KEYBOARD_EQUAL                131
#define PAGE_UP                       132
#define PAGE_DOWN                     133
#define MOLECULE_DISPLAY              134
#define SURFACE_CULLING               135
#define SURFACE_DENSITY               136

#define R12_SAMPLES                   138  /* obsolete */
#define KEYBOARD_P                    139 

#define DUMP_RECEPTOR                 141
#define DUMP_LIGAND                   142
#define KEYBOARD_T                    143
#define KEYPAD_HOME                   144
#define KEYPAD_END                    145
#define OPEN_COMPLEX                  146
#define CLOSE_COMPLEX                 147
#define SAVE_COMPLEX                  148
#define DOCKING_SEARCH_MODE           149
#define DOCKING_GRID_SIZE             150
#define RUN_MACRO                     152
#define LEARN_LESSON                  153
#define SAVE_LESSON                   154
#define FORGET_LESSON                 155
#define DISPLAY_HETERO                156
#define R12_CALCULATE                 157
#define EXIT                          159
#define ALPHA_RANGE_ANGLE             160
#define RECEPTOR_RANGE_ANGLE          161
#define LIGAND_RANGE_ANGLE            162
#define MEMORY_STATUS                 163
#define MEMORY_REPORTS                164
#define RANDOMISE_MOLECULE            165
#define FACTOR_Q                      166
#define SMALL_MEMORY                  167
#define KEYBOARD_ESC                  168
#define FIT_LIGAND                    169
#define SAVE_SYMMETRY                 170
#define ACTIVATE_MATCHING             172
#define EVALUATE_MATCHING             173
#define PRINT_PS                      175
#define FACTOR_E                      176
#define PRINT_EPS                     177
#define ENABLE_HETERO                 178
#define RESIZE_WINDOW                 179
#define ROTATE                        180
#define TRANSLATE                     181
#define FACTOR_F                      182
#define FACTOR_G                      183
#define DUMP_TEMPLATE                 184
#define DISPLAY_HYDROGEN              185
#define FACTOR_P                      186
#define SYMMETRY_ELEMENT              187
#define SYMMETRY_TYPE                 188
#define DISC_CACHE                    189
#define LAUNCH_BROWSER                190
#define BROWSER_NAME                  191
#define RESET_NORMALISATION           192  /* not used */
#define DOCKING_ELECTROSTATICS        193
#define DEFAULT_SEPARATION            194
#define DOCKING_REFINE                195
#define DISPLAY_HBONDS                196
#define ENABLE_ARGLYS                 197
#define DISPLAY_RECEPTOR_MODEL        198
#define DISPLAY_LIGAND_MODEL          199
#define DOCK_RECEPTOR_MODEL           200
#define DOCK_LIGAND_MODEL             201
#define SESSION_LOG                   202
#define EDIT_MODEL_SCOPE              203
#define DOCKING_R12_STEP              204
#define DOCKING_R12_RANGE             205
#define MACRO_NORMAL_DISPLAY          206
#define MACRO_DOCKING_SEPARATION      207
#define DELETE_MODELS                 208
#define MACRO_SPHERE_DISPLAY          209
#define ACTIVATE_MACRO_MODEL          210
#define MACRO_SPHERE_RADIUS           211
#define MACRO_WEIGHTING               212
#define MACRO_MIN_COVERAGE            213
#define MACRO_AUTO_COVERAGE           214
#define MACRO_MATCH_CHAINS            215
#define MACRO_CHANGE_CHAINS           216
#define DOCKING_SORT_MODE             217
#define DOCKING_CLUSTER_MODE          218
#define DOCKING_CLUSTER_THRESHOLD     219
#define DOCKING_CLUSTER_BUMPS         220
#define PIPE_BOTH                     221
#define SNAP_ORIGINS                  222
#define RECEPTOR_ORIGIN               223
#define LIGAND_ORIGIN                 224
#define RECEPTOR_CHANGE_ORIGIN        225
#define LIGAND_CHANGE_ORIGIN          226
#define AUTO_POSITIONING              227
#define SAVE_RANGE                    228
#define STICKY_FSB                    229
#define MAKE_TM                       230
#define SAVE_SUMMARY                  231
#define SOLID_SCOPE                   232
#define DISPLAY_SOLID_PEPTIDES        233
#define DISPLAY_SOLID_AROMATICS       234
#define SOLID_PEPTIDE_COLOUR          235
#define SOLID_AROMATIC_COLOUR         236
#define SOLID_AROMATIC_FACETS         237
#define SOLID_AROMATIC_THICKNESS      238
#define SOLID_PEPTIDE_THICKNESS       239
#define SOLID_COLOUR_SCHEME           240
#define SOLID_RESOLUTION              241
#define SOLID_SPHERE_SIZE             242
#define SOLID_CYLINDER_SIZE           243
#define LIGHT_POSITION                244
#define LIGHT_COLOUR                  245
#define LIGHT_EXPONENT                246
#define DISPLAY_BACKGROUND            248 
#define PRINT_JPEG                    249
#define PRINT_PNG                     250
#define MESH_DISPLAY                  251
#define MESH_GRID                     252
#define MESH_STYLE                    253
#define MESH_SCOPE                    254
#define MESH_DRAW                     255
#define MESH_CULLING                  256
#define LIGHT_SPECULAR                257
#define MESH_VIEWING_PROBE            258
#define MESH_DOCKING_PROBE            259
#define MESH_SOFTNESS                 260
#define DOCKING_RADIAL_FILTER         261
#define SAVE_DOCKING                  262
#define RADIAL_CUTOFF                 263
#define RADIAL_DISTANCE               264
#define FULL_SCREEN                   265
#define KEYBOARD_F                    266
#define DISPLAY_SCORE                 267
#define MOVIE_TYPE                    268
#define MOVIE_CYCLE                   269
#define MOVIE_RATE                    270
#define ACTIVATE_MOVIE                271
#define ACTIVATE_SPIN                 272
#define SPIN_ANGLE                    273
#define ENABLE_SPIN                   274
#define MESH_COLOUR                   275
#define MESH_TYPE                     276

#define MESH_RAMP_VALUES              277
#define MESH_RAMP_TYPE                278
#define SIGMA_TAU_SHIFT               279
#define ENABLE_PARALLAX               280
#define PIXEL_PARALLAX                281
#define KEYBOARD_S                    282
#define KEYBOARD_M                    283
#define KEYBOARD_B                    284

#define UNIX_BROWSER1                 285
#define UNIX_BROWSER2                 286
#define WINDOWS_BROWSER1              287
#define WINDOWS_BROWSER2              288
#define WINDOWS_MOUSE                 289
#define ENABLE_PROGRESS               290
#define POPUP_PROGRESS                291

#define MESH_RAMP_MINCOL              292
#define MESH_RAMP_MIDCOL              293
#define MESH_RAMP_MAXCOL              294

#define DISPLAY_CARTOON               295
#define CARTOON_SCOPE                 296
#define CARTOON_STYLE                 297
#define SOLID_STYLE                   298  /* formerly MOLECULE_STYLE */
#define CARTOON_HELIX_COLOUR_ID       299
#define CARTOON_SHEET_COLOUR_ID       300
#define CARTOON_TURN_COLOUR_ID        301
#define CARTOON_LOOP_COLOUR_ID        302
#define CARTOON_310_COLOUR_ID         303
#define CARTOON_HELIX_EDGE_COLOUR_ID  304
#define CARTOON_SHEET_EDGE_COLOUR_ID  305
#define CARTOON_TURN_EDGE_COLOUR_ID   306
#define CARTOON_310_EDGE_COLOUR_ID    307
#define CARTOON_RADIUS                308
#define CARTOON_WIDTH                 309
#define CARTOON_RESOLUTION            310
#define CARTOON_COLOUR_SCHEME         311
#define DOCKING_CORRELATION           312
#define POTENTIAL_GRID_SIZE           313
#define THREEFOLD_HACK                314
#define UNIFY_MODELS                  315
#define DOCKING_R12_SUBSTEPS          316
#define SAVE_MATRIX                   317
#define SAVE_TRANSFORM                318
#define MAX_DOCKING_SOLUTIONS         319
#define MAX_MATCHING_SOLUTIONS        320
#define DOCKING_CLUSTER_WINDOW        321

#define DARWIN_BROWSER1               322
#define DARWIN_BROWSER2               323

#define KEYPAD_LEFT                   324
#define KEYPAD_RIGHT                  325
#define KEYPAD_UP                     326
#define KEYPAD_DOWN                   327

#define SURFACE_RAMP_VALUES           328
#define SURFACE_RAMP_TYPE             329
#define SURFACE_RAMP_MINCOL           330
#define SURFACE_RAMP_MIDCOL           331
#define SURFACE_RAMP_MAXCOL           332

#define SURFACE_TRANSPARENCY          333
#define SURFACE_ALPHA                 334
#define MESH_TRANSPARENCY             335
#define MESH_ALPHA                    336

#define HARMONIC_PROPERTY             337

#define RECEPTOR_INTERFACE            338
#define LIGAND_INTERFACE              339
#define DOT_GRID                      340

#define RESTRAINT_DISTANCE            341
#define RESTRAINT_ENERGY              342
#define RECEPTOR_ACTIVE_RESIDUES      343
#define RECEPTOR_PASSIVE_RESIDUES     344
#define LIGAND_ACTIVE_RESIDUES        345
#define LIGAND_PASSIVE_RESIDUES       346
#define RESTRAINT_ENABLE              347

#define DISPLAY_SKELETON              348

#define DOCKING_POLAR_PROBE           349
#define DOCKING_APOLAR_PROBE          350

#define FACTOR_Z                      351
#define PIPER_WEIGHTS                 352

#define DOCKING_FFT_TYPE              353

#define DOCKING_RECEPTOR_STEPSIZE     354
#define DOCKING_LIGAND_STEPSIZE       355
#define DOCKING_ALPHA_STEPSIZE        356

#define MATCHING_RECEPTOR_STEPSIZE    357
#define MATCHING_LIGAND_STEPSIZE      358
#define MATCHING_ALPHA_STEPSIZE       359

#define MAX_DOCKING_CLUSTERS          360
#define MOLECULE_SHIFT                361
#define DOCKING_FFT_DEVICE            362
#define MATCHING_CORRELATION          363
#define MATCHING_FFT_TYPE             364
#define MATCHING_FFT_DEVICE           365
#define AUTO_ORIGIN                   366
#define MATCHING_BOX_SIZE             367
#define MATCHING_RECEPTOR_RANGE       368
#define MATCHING_LIGAND_RANGE         369
#define MATCHING_ALPHA_RANGE          370
#define MATCHING_R12_SUBSTEPS         371
#define MATCHING_SAMPLING_METHOD      372
#define DOCKING_SAMPLING_METHOD       373
#define DOCKING_BOX_SIZE              374
#define CANONICALISE_MOLECULE         375
#define PROGRESS_FILE                 376

#define GUI_PROGRESS_WEIGHTS_XXXX     377  // not used
#define GUI_NORMAL_DISPLAY            378
#define GUI_NORMAL_COVERAGE           379
#define GUI_NORMAL_AUTO               380

#define GUI_ANIMATION_MODE            381
#define GUI_FRAME_RATE                382
#define GUI_PROGRESS_ENABLE           383
#define GUI_PROGRESS_POPUP            384
#define GUI_PROGRESS_FILE_XXXXX       385   // not used
#define GUI_DOCKING_MAXIMUM           386
#define GUI_DOCKING_SOLUTION          387
#define GUI_MATCHING_MAXIMUM          388
#define GUI_MATCHING_SOLUTION         389
#define GUI_ORIGIN_ATOM               390

#define GUI_INTERFACE_ATOM            391
#define GUI_WINDOWS_MOUSE             392
#define GUI_WINDOW_SIZE               393
#define GUI_RECEPTOR_ORIGIN           394
#define GUI_RECEPTOR_INTERFACE        395
#define GUI_RECEPTOR_ACTIVE           396
#define GUI_RECEPTOR_PASSIVE          397
#define GUI_LIGAND_ORIGIN             398
#define GUI_LIGAND_INTERFACE          399
#define GUI_LIGAND_ACTIVE             400

#define GUI_LIGAND_PASSIVE            401
#define GUI_CENTRE_THING              402
#define GUI_DISPLAY_MODEL             403
#define GUI_DOCKING_MODEL             404
#define GUI_MOLECULE_ROTATION         405
#define GUI_R12_SHIFT                 406
#define GUI_R12_DISTANCE              407
#define GUI_REDRAW                    408
#define GUI_MOVING_MOLECULE           409
#define GUI_AUTO_POSITION             410

#define GUI_AUTO_ORIGIN               411
#define GUI_MRAMP                     412
#define GUI_SRAMP                     413
#define GUI_HEX_JOB                   414
#define GUI_EXIT                      415
#define GUI_BACKGROUND                416
#define GUI_PROJECTION                417
#define GUI_PARALLAX                  418
#define GUI_FULLSCREEN                419

#define GUI_ORIGIN_CURSOR             420
#define GUI_EDIT_CURSOR               421
#define GUI_POINTER_CURSOR            422
#define GUI_RESET_CURSOR              423
#define GUI_BUSY_CURSOR               424
#define GUI_VIEWPORT                  425
#define GUI_STICKY_FSB                426
#define GUI_ACTIVATE_NORMAL           427

#define GUI_LIGHT_SPECULAR            428
#define GUI_LIGHT_EXPONENT            429
#define GUI_LIGHT_POSITION            430
#define GUI_LIGHT_COLOUR              431
#define GUI_REDRAW_PICK               432
#define GUI_PRINT_SCENE               433

#define CARTOON_METHOD                434
#define RESIDUE_COLOUR                435

#define DOCKING_SCORE_THRESHOLD       436
#define MATCHING_SCORE_THRESHOLD      437


#endif /* hex_cmd_h */

/*---------------------------------------------------------------------------*/
