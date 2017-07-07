/*-----------------------------------------------------------------------------
**  File:       hex_units.h
**
**  Author:     Dave Ritchie, 20/12/98
**
**  Purpose:    Set up symbolic constants for various physical/derived units.
**
**		With proteins, the conventional units are usually lengths in 
**		Angstroms, charges in electron units, etc.
**
**		This file also defines some conversion factors to set these
**		quantities in SI units - usually KJ/mol, etc...
**
**-----------------------------------------------------------------------------
**
**  Copyright (C) 1998-2003 D.W. Ritchie, University of Aberdeen.
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

#ifndef hex_units_h
#define hex_units_h

/* first, the fundamental units ... */

#define K_C           2.997925e8   /* M/s        - speed of light */
#define K_B           1.38066e-23  /* J/K        - Boltzmann constant */
#define K_L           6.02205e23   /* 1/mol      - Avogadro constant */
#define K_P           6.62618e-34  /* Js         - Planck's constant */
#define K_PB          1.05459e-34  /* Js         - Planck's constant / 2PI */
#define K_G           6.672e-11    /* N.M^2/Kg^2 - gravitational const */
#define K_M0         (4*PI*1.0e-7) /* Js^2/C^2/M - vacuum permeability */
#define K_E0          8.854188e-12 /* C^2/J/M    - vacuum permittivity */
#define K_E           1.60219e-19  /* C          - electronic unit charge */
#define K_AMU         1.66057e-27  /* Kg         - atomic mass unit */
#define K_PM          1.67265e-27  /* Kg         - proton mass */
#define K_NM          1.67495e-27  /* Kg         - neutron mass */
#define K_EM          9.10953e-31  /* Kg         - electron mass */
#define K_ATM         1.01325e5    /* N/M^2      - atmospheric pressure */
#define K_F           9.64846e4    /* C/mol      - Faraday = K_L*K_E */
#define K_CAL         4.184        /* J          - calorie */
#define K_D           3.33564e-30  /* C.M        - Debye */
#define K_eV          1.602189e-19 /* J          - electron-volt */
#define K_N           1.0          /* J/M        - Newton */
#define K_V           1.0          /* J/C        - Volt */
#define K_W           1.0          /* J/s        - Watt */
#define K_R           8.31441      /* J/K/mol    - Gas constant = K_L*K_B */
#define K_RT          2478.9       /* J/mol      - RT at 298.15K */
#define K_ANG         1.0e-10      /* M          - Angstrom */
#define K_J_KJM       (K_L*1.0e-3) /* KJ/mol     - Joules to KJ/mol */
#define K_J_KJ        (1.0e-3)     /* KJ/mol     - Joules to KJ */

/* now, some working constants and conversion factors */

#define K_CAL_J       4.184        /* J          - calorie to J */
#define K_KCAL_KJ     4.184        /* KJ         - Kcalorie to KJ */

/* now, some derived constants and conversion factors (from Parasurf) */

#define BohrToAngstr  0.52917706
#define AngstrToBohr  1.88972667
#define MEPconst      332.07           /* charge/Angstrom -> kcal/mol */
#define DebyeToBohr   0.3934301
#define BohrToDebye   2.541748
#define AOtoAngstr    0.529167
#define AngstrToDebye 4.803207

#define K_EPRP        8    /* effective protein rel permittivity */

/* convert internal electrostatic energy (charges, Angstroms) into KJ/mol */

#define K_ES_KJM      (K_E*K_E/K_ANG/(4*PI*K_E0)*K_J_KJM)

#define K_ES_KJM_EPRP (K_ES_KJM/K_EPRP)    /* including rel permittivity */

/* convert internal electrostatic potential into Volts */

#define K_ES_V        (K_E/K_ANG/(4*PI*K_E0))

#define K_ES_V_EPRP   (K_ES_V/K_EPRP)       /* including rel permittivity */

#endif /* hex_units_h */
/*---------------------------------------------------------------------------*/
