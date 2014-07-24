# Global Variables

AM_CFLAGS           = @DEFAULT_CFLAGS@
AM_CFLAGS          += $(OPENMP_CFLAGS) 
AM_LDFLAGS          = $(OPENMP_CFLAGS)


# Debug 
if DEBUG
AM_CFLAGS          += -g @DEBUG_CFLAGS@ -DDEBUG
AM_LDFLAGS         += -g
else
AM_CFLAGS          += @OPT_CFLAGS@
endif

# Verbose
if VERBOSE
AM_CFLAGS          += -DVERBOSE
endif

# Which functions/groups should be parallel
if PARALLEL
AM_CFLAGS          += -D__BABINET_C_OPENMP
AM_CFLAGS          += -D__GRADIENT_C_OPENMP
AM_CFLAGS          += -D__FFT_C_OPENMP
endif




# Precision
if SINGLE

include_HEADERS     = pdenf.h
EXTRA_DIST          = pden.in.h
CLEANFILES          = $(include_HEADERS)
do_subst = sed -e 's/real/float/g' 

AM_CFLAGS          += -DFLOAT=float   -Dreal=float

# Build library
lib_LIBRARIES       = libpdenf.a
libpdenf_a_SOURCES   = babinet.c fft.c gradient.c pden.c calc.c mrc.c xplor.c powerspectra.c render.c sfrefine.c split2fft.c

#Lib
LDADD 		    = libpdenf.a -lfftw3f

else #Double
include_HEADERS     = pden.h
EXTRA_DIST          = pden.in.h
CLEANFILES          = $(include_HEADERS)
do_subst = sed -e 's/real/double/g' 
AM_CFLAGS          += -DDOUBLE=double -Dreal=double

# Build library
lib_LIBRARIES       = libpden.a
libpden_a_SOURCES   = babinet.c fft.c gradient.c pden.c calc.c mrc.c xplor.c structurefactor.c render.c sfrefine.c split2fft.c io.c

#Lib
LDADD               = libpden.a -lfftw3
endif

# Headers
noinst_HEADERS      = atomtypes.h mathtools.h pden_sfrefine.h split2fft.h tools.h types.h



pden.h: pden.in.h
	$(do_subst) < pden.in.h > pden.h

pdenf.h: pden.in.h
	$(do_subst) < pden.in.h > pdenf.h



# Build binaries
bin_PROGRAMS        = pdcalcsf pdapplysf pdnormalize pdrefinesf pdcorr
pdcorr_SOURCES      = pdcorr.c
pdcorr_CFLAGS       = -DVERBOSE $(AM_CFLAGS)
pdapplysf_SOURCES   = pdapplysf.c
pdapplysf_CFLAGS    = -DVERBOSE $(AM_CFLAGS)
pdcalcsf_SOURCES    = pdcalcsf.c
pdcalcsf_CFLAGS     = -DVERBOSE $(AM_CFLAGS)
pdnormalize_SOURCES = pdnormalize.c
pdnormalize_CFLAGS  = -DVERBOSE $(AM_CFLAGS)
pdrefinesf_SOURCES  = pdrefinesf.c
pdrefinesf_CFLAGS   = -DVERBOSE $(AM_CFLAGS)
