/******************************************************************************
**
**  File:    hex_bitmaps.h
**
**  These were copied/inspired by the X11 bit maps supplied in:
**
**     /usr/include/X11/bitmaps or: /usr/openwin/share/include/X11/bitmaps
**
**  using the X11 "bitmap" program.
**
**  Dave Ritchie,  1996  First implementation.
**  Dave Ritchie,  2003  Moved to own include file.
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

#ifndef hex_bitmaps_h
#define hex_bitmaps_h

#define left_ptr_width 16
#define left_ptr_height 16
static unsigned char left_ptr_bits[] = {
   0x00, 0x00, 0x08, 0x00, 0x18, 0x00, 0x38, 0x00, 0x78, 0x00, 0xf8, 0x00,
   0xf8, 0x01, 0xf8, 0x03, 0xf8, 0x07, 0xf8, 0x00, 0xd8, 0x00, 0x88, 0x01,
   0x80, 0x01, 0x00, 0x03, 0x00, 0x03, 0x00, 0x00};

#define tie_fighter_width 16
#define tie_fighter_height 16
static unsigned char tie_fighter_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x08, 0x10, 0x04, 0x20, 0x02, 0x40, 0x02, 0x40,
   0xe2, 0x47, 0x3e, 0x7c, 0x12, 0x48, 0x3e, 0x7c, 0xe2, 0x47, 0x02, 0x40,
   0x42, 0x42, 0x64, 0x26, 0x28, 0x14, 0x00, 0x00};

#define id_width 16
#define id_height 16
static unsigned char id_bits[] = {
   0x00, 0x00, 0xc6, 0x0f, 0xc6, 0x18, 0xc6, 0x30, 0xc6, 0x30, 0xc6, 0x30,
   0xc6, 0x30, 0xc6, 0x30, 0xc6, 0x30, 0xc6, 0x30, 0xc6, 0x30, 0xc6, 0x30,
   0xc6, 0x30, 0xc6, 0x18, 0xc6, 0x0f, 0x00, 0x00};

#define dag_width 16
#define dag_height 16
static unsigned char dag_bits[] = {
   0x00, 0x1c, 0x00, 0x0e, 0x00, 0x07, 0x80, 0x03, 0xc0, 0x01, 0x60, 0x00,
   0xf0, 0x3f, 0xf8, 0x1f, 0x00, 0x0e, 0x00, 0x07, 0x80, 0x03, 0xd0, 0x01,
   0xf0, 0x00, 0x70, 0x00, 0xf0, 0x00, 0x00, 0x00};

#define rdag_width 16
#define rdag_height 16
static unsigned char rdag_bits[] = {
   0xff, 0xe3, 0xff, 0xf1, 0xff, 0xf8, 0x7f, 0xfc, 0x3f, 0xfe, 0x9f, 0xff,
   0x0f, 0xc0, 0x07, 0xe0, 0xff, 0xf1, 0xff, 0xf8, 0x7f, 0xfc, 0x2f, 0xfe,
   0x0f, 0xff, 0x8f, 0xff, 0x0f, 0xff, 0xff, 0xff};

#define hammer_width 16
#define hammer_height 16
static unsigned char hammer_bits[] = {
   0x00, 0x00, 0x00, 0x0e, 0x80, 0x07, 0xe0, 0x01, 0xf8, 0x00, 0xfe, 0x00,
   0xff, 0x01, 0xbf, 0x03, 0x1e, 0x07, 0x0c, 0x0e, 0x00, 0x1c, 0x00, 0x38,
   0x00, 0x70, 0x00, 0xe0, 0x00, 0x00, 0x00, 0x00};

#define dotsurf_width 16
#define dotsurf_height 16
static unsigned char dotsurf_bits[] = {
   0x00, 0x00, 0x78, 0x00, 0xcc, 0x00, 0x86, 0x01, 0x32, 0x01, 0xca, 0x0f,
   0x26, 0x32, 0x0c, 0x29, 0xf8, 0x4c, 0xac, 0x5e, 0xc6, 0x48, 0x82, 0x61,
   0xc6, 0x33, 0x6c, 0x1c, 0x38, 0x00, 0x00, 0x00};


#define cart_width 16
#define cart_height 16
static unsigned char cart_bits[] = {
   0x00, 0x00, 0x00, 0x60, 0x00, 0x70, 0xe0, 0x3f, 0xc0, 0x1f, 0xc0, 0x18,
   0xe0, 0x18, 0x70, 0x18, 0x18, 0x1e, 0x0c, 0x1f, 0x8c, 0x13, 0x98, 0x01,
   0xfc, 0x00, 0x6e, 0x00, 0x06, 0x00, 0x00, 0x00};

#define vdw_width 16
#define vdw_height 16
static unsigned char vdw_bits[] = {
   0x00, 0x00, 0x78, 0x00, 0xfc, 0x00, 0xfe, 0x01, 0xce, 0x01, 0xee, 0x0d,
   0xfe, 0x3f, 0xfc, 0x3f, 0xf8, 0x7f, 0xfc, 0x73, 0xce, 0x7b, 0xee, 0x7f,
   0xfe, 0x3e, 0x7c, 0x1c, 0x38, 0x00, 0x00, 0x00};

#define surf_width 16
#define surf_height 16
static unsigned char surf_bits[] = {
   0x00, 0x00, 0x78, 0x00, 0xfc, 0x01, 0x86, 0x1f, 0x06, 0x3e, 0x06, 0x30,
   0x06, 0x30, 0x0c, 0x30, 0x0c, 0x60, 0x0c, 0x60, 0x06, 0x60, 0x86, 0x63,
   0xc6, 0x3f, 0x7c, 0x1c, 0x38, 0x00, 0x00, 0x00};

#define rsurf_width 16
#define rsurf_height 16
static unsigned char rsurf_bits[] = {
   0xff, 0xff, 0x87, 0xff, 0x03, 0xfe, 0x79, 0xe0, 0xf9, 0xc1, 0xf9, 0xcf,
   0xf9, 0xcf, 0xf3, 0xcf, 0xf3, 0x9f, 0xf3, 0x9f, 0xf9, 0x9f, 0x79, 0x9c,
   0x39, 0xc0, 0x83, 0xe3, 0xc7, 0xff, 0xff, 0xff};

#define harm_width 16
#define harm_height 16
static unsigned char harm_bits[] = {
   0x00, 0x00, 0x00, 0x1c, 0x00, 0x2e, 0x00, 0x25, 0x30, 0x23, 0x70, 0x19,
   0xe0, 0x0e, 0x40, 0x01, 0xb8, 0x02, 0x7c, 0x07, 0x4a, 0x06, 0x66, 0x00,
   0x32, 0x00, 0x1e, 0x00, 0x00, 0x00, 0x00, 0x00};

#define axes_width 16
#define axes_height 16
static unsigned char axes_bits[] = {
   0x00, 0x00, 0x18, 0x00, 0x3c, 0x00, 0x7e, 0x00, 0x18, 0x00, 0x18, 0x00,
   0x18, 0x00, 0x18, 0x00, 0x18, 0x00, 0x18, 0x08, 0x18, 0x18, 0xf8, 0x3f,
   0xf8, 0x3f, 0x00, 0x18, 0x00, 0x08, 0x00, 0x00};

#define axis_width 16
#define axis_height 16
static unsigned char axis_bits[] = {
   0x00, 0x00, 0x00, 0x78, 0x00, 0x70, 0x00, 0x78, 0x00, 0x5c, 0x00, 0x0e,
   0x00, 0x07, 0x80, 0x03, 0xc0, 0x01, 0xe0, 0x00, 0x70, 0x00, 0x3a, 0x00,
   0x1e, 0x00, 0x0e, 0x00, 0x1e, 0x00, 0x00, 0x00};

#define home_width 16
#define home_height 16
static unsigned char home_bits[] = {
   0x00, 0x00, 0x80, 0x0c, 0xc0, 0x0d, 0x60, 0x0f, 0x30, 0x0e, 0x18, 0x0c,
   0x0c, 0x18, 0x06, 0x30, 0xe7, 0x71, 0x2c, 0x19, 0x2c, 0x19, 0x6c, 0x19,
   0x2c, 0x19, 0x2c, 0x19, 0xfc, 0x1f, 0x00, 0x00};

/*
static unsigned char ulock_bits[] = {
   0x00, 0x00, 0x00, 0x07, 0xc0, 0x01, 0x60, 0x00, 0x20, 0x00, 0x20, 0x00,
   0xf8, 0x1f, 0x08, 0x10, 0x08, 0x10, 0x18, 0x18, 0x10, 0x08, 0x30, 0x0c,
   0x20, 0x04, 0x60, 0x06, 0xc0, 0x03, 0x00, 0x00};
*/

#define ulock_width 16
#define ulock_height 16
static unsigned char ulock_bits[] = {
   0x00, 0x00, 0xe0, 0x0f, 0xe0, 0x00, 0x20, 0x00, 0x20, 0x00, 0x20, 0x00,
   0xf8, 0x1f, 0x08, 0x10, 0x88, 0x11, 0x98, 0x19, 0x10, 0x09, 0x30, 0x0d,
   0x20, 0x04, 0x60, 0x06, 0xc0, 0x03, 0x00, 0x00};

/* empty home */

#define ehome_width 16
#define ehome_height 16
/*
static unsigned char ehome_bits[] = {
   0x00, 0x00, 0x80, 0x00, 0xc0, 0x01, 0x60, 0x03, 0x30, 0x06, 0x18, 0x0c,
   0x0c, 0x18, 0x06, 0x30, 0x07, 0x70, 0x0c, 0x18, 0x0c, 0x18, 0x0c, 0x18,
   0x0c, 0x18, 0x0c, 0x18, 0xfc, 0x1f, 0x00, 0x00};
*/

#define fill_width 16
#define fill_height 16
static unsigned char fill_bits[] = {
   0x00, 0x00, 0x00, 0x78, 0x00, 0x70, 0x00, 0x78, 0x20, 0x5c, 0x70, 0x08,
   0xf8, 0x00, 0xfc, 0x01, 0xfe, 0x03, 0xff, 0x07, 0xfe, 0x03, 0xfc, 0x01,
   0xf8, 0x00, 0x70, 0x00, 0x20, 0x00, 0x00, 0x00};

#define edo_width 16
#define edo_height 16
static unsigned char edo_bits[] = {
   0x00, 0x30, 0x00, 0x30, 0x00, 0xfc, 0x00, 0xfc, 0x00, 0x30, 0x00, 0x30,
   0x00, 0x03, 0x80, 0x03, 0xc0, 0x01, 0xe0, 0x00, 0x4c, 0x00, 0x0c, 0x00,
   0x3f, 0x00, 0x3f, 0x00, 0x0c, 0x00, 0x0c, 0x00};

#define redo_width 16
#define redo_height 16
static unsigned char redo_bits[] = {
   0xff, 0xcf, 0xff, 0xcf, 0xff, 0x03, 0xff, 0x03, 0xff, 0xcf, 0xff, 0xcf,
   0xff, 0xfc, 0x7f, 0xfc, 0x3f, 0xfe, 0x1f, 0xff, 0xb3, 0xff, 0xf3, 0xff,
   0xc0, 0xff, 0xc0, 0xff, 0xf3, 0xff, 0xf3, 0xff};

#define lock_width 16
#define lock_height 16
static unsigned char lock_bits[] = {
   0x00, 0x00, 0xc0, 0x03, 0x60, 0x06, 0x20, 0x04, 0x20, 0x04, 0xf8, 0x1f,
   0x08, 0x10, 0x88, 0x11, 0x98, 0x19, 0x10, 0x09, 0x30, 0x0d, 0x20, 0x04,
   0x60, 0x06, 0xc0, 0x03, 0x00, 0x00, 0x00, 0x00};


#define clip_width 16
#define clip_height 16
static unsigned char clip_bits[] = {
   0x00, 0x00, 0xf8, 0x0f, 0x0c, 0x08, 0x06, 0xcc, 0x02, 0xe6, 0xfe, 0xb3,
   0x00, 0x98, 0xfe, 0x8b, 0x02, 0x8a, 0x02, 0x8a, 0x02, 0xca, 0x02, 0x6a,
   0x02, 0x3a, 0xfe, 0x1b, 0x00, 0x00, 0x00, 0x00};
/*
static unsigned char clip_bits[] = {
   0x00, 0x00, 0xf0, 0x00, 0xf8, 0x00, 0xcc, 0x00, 0xc6, 0x00, 0xc6, 0x10,
   0xc6, 0x30, 0xf6, 0x7f, 0xf6, 0x7f, 0xc6, 0x30, 0xc6, 0x10, 0xc6, 0x00,
   0x66, 0x00, 0x3e, 0x00, 0x1e, 0x00, 0x0e, 0x00};
*/

#define schain_width 16
#define schain_height 16
static unsigned char schain_bits[] = {
   0x00, 0x01, 0x80, 0x03, 0xc0, 0x06, 0x60, 0x0c, 0x20, 0x08, 0x20, 0x08,
   0x60, 0x0c, 0xc0, 0x06, 0x80, 0x03, 0x00, 0xc1, 0x00, 0x61, 0x00, 0x31,
   0xf0, 0x1f, 0x18, 0x00, 0x0c, 0x00, 0x06, 0x00};

#define cc_width 16
#define cc_height 16
static unsigned char cc_bits[] = {
   0x00, 0x00, 0xc6, 0x63, 0x7e, 0x7e, 0x1c, 0x38, 0x3c, 0x3c, 0x74, 0x2e,
   0xe6, 0x67, 0xc2, 0x43, 0xc2, 0x43, 0xe6, 0x67, 0x74, 0x2e, 0x3c, 0x3c,
   0x1c, 0x38, 0x7e, 0x7e, 0xc6, 0x63, 0x00, 0x00};

#define six_width 16
#define six_height 16
static unsigned char six_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x30, 0x0c, 0x18, 0x06, 0x0c,
   0x02, 0x04, 0x1a, 0x34, 0x36, 0x6c, 0x22, 0x44, 0x22, 0x44, 0x22, 0x44,
   0xa2, 0x45, 0x9c, 0x39, 0x00, 0x00, 0x00, 0x00};

/*
static unsigned char cc_bits[] = {
   0x00, 0x00, 0xc0, 0x03, 0x78, 0x1e, 0x0c, 0x30, 0x04, 0x20, 0x64, 0x26,
   0xe6, 0x67, 0x42, 0x42, 0x42, 0x42, 0xe6, 0x67, 0x64, 0x26, 0x04, 0x20,
   0x0c, 0x30, 0x78, 0x1e, 0xc0, 0x03, 0x00, 0x00};
*/

#endif /* hex_bitmaps_h */

/*------------------------------------------------------------------------*/
