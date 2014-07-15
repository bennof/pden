


# Global Variables
ACLOCAL_AMFLAGS  =

AM_CFLAGS        = -pedantic -Wall -std=c99 -O2 $(OPENMP_CFLAGS)
AM_LDFLAGS       = $(OPENMP_CFLAGS)


# Debug 
if IS_DEBUG
AM_CFLAGS       += -g -ggdb -O0
AM_LDFLAGS      += -g
endif

# Precision
if IS_SINGLE
AM_CFLAGS       += -DFLOAT=float   -Dreal=float
else
AM_CFLAGS       += -DDOUBLE=double -Dreal=double
endif

#Lib
LDADD = libpden.a -lfftw3

# Headers
include_HEADERS     = pden.h

# Build library
lib_LIBRARIES    = libpden.a
libpden_a_SOURCES   = babinet.c fft.c gradient.c pden.c calc.c mrc.c xplor.c powerspectra.c render.c sfrefine.c split2fft.c

# Build binaries
bin_PROGRAMS        = pdcalcps pdapplysf pdnormalize pdrefinesf
pdapplysf_SOURCES   = pdapplysf.c
pdcalcps_SOURCES    = pdcalcps.c
pdnormalize_SOURCES = pdnormalize.c
pdrefinesf_SOURCES  = pdrefinesf.c
