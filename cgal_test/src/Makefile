#------------------------------------------------------------------------------
#
#  Makefile.linux - to build the PDBS program...
#
#  Dave Ritchie, 08/03/2013
#
#------------------------------------------------------------------------------

HEX_VERSION=8.1.1

#------------------------------------------------------------------------------

PDBS_HEX_LIBS=../lib
PDBS_HEX_INCLUDES=../include

#------------------------------------------------------------------------------

PDBS_VAR_LIBS=../var/Linux64/lib
PDBS_VAR_INCLUDES=../var/Linux64/include

PDBS_EXT      = x64
PDBS_OPTS     = -Dhex_linux64 -fPIC -fPIE

PDBS_F77      = g77
PDBS_CC       = gcc
PDBS_CPP      = icc
PDBS_CPP      = g++
PDBS_COPTS    = $(PDBS_OPTS) -O2
PDBS_COPTS    = $(PDBS_OPTS) -g
PDBS_INCLUDES = -I. -I$(PDBS_HEX_INCLUDES) -I$(PDBS_VAR_INCLUDES)



PDBS_LD       = icc
PDBS_LD       = g++

PDBS_LDOPTS   = -no-ipo -static-libstdc++
PDBS_LDOPTS   = -no-ipo -static
PDBS_LDOPTS   = -static

#                -lmol_$(HEX_VERSION)_Linux_x64 \

PDBS_LIBS     =  -L$(PDBS_HEX_LIBS)  -L$(PDBS_VAR_LIBS) -L/usr/lib \
                -lhex_8.1.1_Linux_x64  \
                -Wl,-Bstatic \
                -lgmp \
                -lbz2 -lz \
                -lpthread \
                -ldl \
                -lm\


include makefile.include
