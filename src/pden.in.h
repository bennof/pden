/*
 *
 *  Copyright 2014 Benjamin Falkner 
 *
 *  This file is part of PDen.
 *
 *  PDen is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PDen is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PDen.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __PDEN_H
#define __PDEN_H 0.3
#define __PDEN_VERSION__ "0.3"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file pden.in.h
 * @author Benjamin Falkner
 * @brief File containing all functions that should be directly used from this library
 * This header will be installed as pden.h (single precision) or pdenf.h (double precision). 
 * Depending on selected precision real will be replaced by real or double 
 */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief helps to construct 3 dimensional structures
 */
#define VEC3(value) struct{value x; value y; value z;}

/**
  *  @brief bitflag to signal that a density is in phase space
  *  use as bitflags because the other bits are reserved for future usage
  */
#define PDEN_MODE_PHASE_SPACE 0x1
/**
  *  @brief bitflag to signal that a density is in real space
  *  use as bitflags because the other bits are reserved for future usage
  */
#define PDEN_MODE_REAL_SPACE  0x0

#define PDEN_RENDER_NM_COORDS 1

/**
 * @brief Split2FFT object
 * hidden field - only used as reference - do not try to access this field
 */
typedef struct Split2FFT_struct Split2FFT_t;

/** 
  * @brief protein density map
  * structure containing all informations about an protein density map
  */
typedef struct PDen_struct {
	int           mode;    /**< Mode bitflag only for real space / phase space used. (other bits are reserved for future usage) */
	size_t	      n;       /**< number of records of the density in real space */
	VEC3(size_t)  size;    /**< grid size */
	VEC3(real)    apix;    /**< Angstrom per pixel*/
	VEC3(real)    origin;  /**< origin */
	VEC3(size_t)  cell;    /**< cell size */
	VEC3(ssize_t) start;   /**< cell start */
	Split2FFT_t   *fft;    /**< Split2FFT reference */
	real          *data;   /**< map data buffer */
} PDen_t;

// Constructor & Destructor
/** 
 * @brief empty constructor 
 * @return empty reference of a protein density map
 */
PDen_t * pDenNew();

/** 
 * @brief  constructor from reference 
 * @param reference a reference map
 * @return a reference of a desnsity map derived from another map 
 */
PDen_t * pDenNewFrom( PDen_t * reference );

/** 
 * @brief defined constructor 
 * @param n grid size
 * @param apix Angstrom per pixel
 * @param origin an origine vector 
 * @return a reference of the derived desnsity map 
 */
PDen_t * pDenNewE(size_t n[],real apix[], real origin[]);

/** 
 * @brief simple defined constructor 
 * @param n grid size
 * @param min minimal coordinates
 * @param max maximal coordinates 
 * @param apix Angstrom per pixel
 * @param offset a set of empty space around the given min and max
 * @return a reference of the derived desnsity map 
 */
PDen_t * pDenNewESimple(real min[], real max[],real apix, real offset);

/** 
 * @brief constructor from atomic model 
 * @param coords atomic coords
 * @param apix Angstrom per pixel
 * @param origin an origine vector 
 * @param offset a set of empty space around the protein
 * @return a reference of the derived desnsity map 
 */
PDen_t * pDenNewFromMol( real coords[], size_t natoms, real apix, real offset);

/**
  * @brief allocate memory
  * use this function to allocate memory for the protein. Using default functions (malloc/calloc) could cause segmentation faults using fft 
  * param this a density map
  * return 0 - will be changed to (PDen_t*) 
  */

int      pDenAlloc( PDen_t * this );

/**
  *  @brief Explicit setup of the fft
  *  lib will perform this implicit
  *  @return 0;
  */
int      pDenSetupFFT(PDen_t * this,int optimize);

/**
 *  @brief Delete a density map
 *  performes a full deconstruction of the density map and frees 
 *  @param this object to deconstruct
 *  @return null;
 */
PDen_t * pDenDelete ( PDen_t * this );

//Access
/**
  *  @brief useless this is c not c++
  */
real *   getData( PDen_t * this );

// IO (mode is not used)
/**
 * @brief print info to stream 
 */
void     pDenPrint     (PDen_t * this, FILE * stream);

/**
 * @brief read xplor file (mode=0) 
 */
int      pDenReadXPLOR (PDen_t * this, const char * filename, int mode);

/**
 * @brief read mrc file (mode=0)
 */
int      pDenReadMRC   (PDen_t * this, const char * filename, int mode);

/**
 * @brief write mrcfile (mode=0)
 */
int      pDenWriteMRC  (PDen_t * this, const char * filename, int mode);

//(something like) member functions

/**
 * @brief clear map
 */
int      pDenClear(PDen_t * this);

/**
 * @brief calc mean value and square deviation
 */
int      pDenMeanSd(PDen_t * this,real *mean, real *sd);

/**
 * @brief nomalize helper function - applies 
 */
int      pDenNormalizeTo( PDen_t * this, real m, real s );

/**
 * @brief normalize map (mu=0.0, s=1.0) 
 */
int      pDenNormalize( PDen_t * this );

/**
 * @brief scale map by value (inplace) 
 */
int      pDenScaleBy(PDen_t * this, real value);

/**
 * @brief multiply map by another map (inplace)
 */
int      PDenScaleMap(PDen_t * this, PDen_t * map);


//statistics
/**
 * @brief calculate correlatio
 */
real pDenCorr(PDen_t * map1, PDen_t * map2);

/**
 * @brief compare to maps using correlation and rsr
 */
int  pDenMapCompare(PDen_t * map1, PDen_t * map2, real *corr, real *rsr);

/**
 * @brief calculate correlation error usng bootstrapping (slow)
 */
int  pDenMapCorrError(PDen_t * map1, PDen_t * map2, real *corr, real *e_corr, size_t nboot);


// external functions
/**
 * @brief add to maps
 */
PDen_t * pDenAdd           (PDen_t * result, PDen_t * a, PDen_t * b );

/**
 * @brief substract maps
 */
PDen_t * pDenSub           (PDen_t * result, PDen_t * a, PDen_t * b );

/**
 * @brief multiply maps
 */
PDen_t * pDenMult          (PDen_t * result, PDen_t * a, PDen_t * b );

/**
 * @brief devide maps
 */
PDen_t * pDenDiv           (PDen_t * result, PDen_t * a, PDen_t * b );

/**
 * @brief sqrt of a map
 */
PDen_t * pDenSqrt          (PDen_t * result, PDen_t * a );

/**
 * @brief inverse sqrt of a map (1/sqrt(x))
 */
PDen_t * pDenISqrt         (PDen_t * result, PDen_t * a );


//TODO: threshold / mask

/**
 * @brief scale a map
 */
PDen_t * pDenScale         (PDen_t * result, PDen_t * a, real s );

/**
 * @brief add white/grey noise with sigma s to the density
 */
PDen_t * pDenAddNoise      (PDen_t * result, PDen_t * a, real s );

/**
 * @brief add gaussian noise with sigma s to the density  
 */
PDen_t * pDenAddGaussNoise (PDen_t * result, PDen_t * a, real s);

/**
 * @brief create a mask at a level threshold
 */
PDen_t * pDenMkMask        (PDen_t * result, PDen_t * a, real level );

//fft

/**
 * @brief perform fft with normalization only after backtransform
 * @deprecated this function is deprecated because of a calculation overhead 
 */
int      pDenFFT(PDen_t * this,int direction); // do not combine with aother functions - use pDenFFTNormal

/**
 * @brief perform fft with normalization
 */
int      pDenFFTNormal(PDen_t * this,int direction);

/**
 * @brief get the size of the split sets work/free 
 */
int 	 pDenGetSizeSplitCV     (PDen_t * input, real lb1,real ub1,real lb2,real ub2, size_t *nused, size_t *nfree);

/**
 * @brief split map in work/free set for crossvalidation 
 */
PDen_t * pDenSplitCV       (PDen_t * result, PDen_t * free, PDen_t * input,real lb1,real ub1,real lb2,real ub2);

/**
 * @brief calculate R-factor 
 */
real     pDenRFactor(PDen_t * a, PDen_t * b);


/**
 * @brief add the fourier space representaions of two maps
 */
PDen_t * pDenAddFS (PDen_t * result, PDen_t * a, PDen_t * b );

/**
 * @brief substract the fourier space representaions of two maps
 */
PDen_t * pDenSubFS (PDen_t * result, PDen_t * a, PDen_t * b );


/**
 * @brief shift density inside the grid (fourier based)
 */
PDen_t * pDenShiftInGrid(PDen_t * result, PDen_t * a,real dx,real dy, real dz);

/**
 * @brief apply a Gaussian filter
 */
PDen_t * pDenGaussFilter(PDen_t * result, PDen_t * a,real sigma);

/**
 * @brief apply Laplace filter 
 */
PDen_t * pDenLaplaceFilter(PDen_t * result, PDen_t * a,real sigma);

/**
 * @brief apply ramp filter 
 */
PDen_t * pDenRampFilter(PDen_t * result, PDen_t * a,real sigma);








/**
 * @brief calculate the power spectrum 
 */
int      pDenCalcPS(PDen_t * input,real *data, size_t n);

//Render
PDen_t * pDenRenderPeng2(PDen_t * result, real width, real sigma, size_t n, int * type, real * coords, real * cfactor, real *dfactor);//depreated
PDen_t * pDenRenderGauss2(PDen_t * result, real width, real sigma, size_t n, int * type, real * coords, real * cfactor, real *dfactor);//depreated

/**
 * @brief render an atomic structure into a density using peng kernel
 */
PDen_t * pDenRenderPeng(PDen_t * result, real width, real sigma, size_t n, int * type, real * coords, real * factor);

/**
 * @brief render an atomic structure into a density using simple gaussian
 */
PDen_t * pDenRenderGauss(PDen_t * result, real width, real sigma, size_t n, int * type, real * coords, real * factor);

//Gradient

/**
 * @brief calculate gradient of the density using roberts-like filter (2 x 2 x 2) at atomic positions (fast / simple)
 */
int pDenGetGradient      (PDen_t * this, PDen_t *ref, real *grad, real *x, size_t natoms);

/**
 * @brief calculate ranged gradient of the density using scharr-like filter (width x width x width) at atomic positions (slow / no advantage have been observed) 
 */
int pDenGetRangedGradient(PDen_t * this, PDen_t *ref, real *grad, real *x, size_t natoms,  real width);

/**
 * @brief calculate gradient of the density using roberts-like filter (2 x 2 x 2) at atomic positions (fast / simple) and add the value to the atomic coordiante
 */
int pDenAddGradient      (PDen_t * this, PDen_t *ref, real *x, size_t natoms,  real width);

/**
 * @brief calculate ranged gradient of the density using scharr-like filter (width x width x width) at atomic positions (slow / no advantage have been observed) and add the value to the atomic coordiante 
 */
int pDenAddRangedGradient(PDen_t * this, PDen_t *ref, real *x, size_t natoms,  real width, real weight);

// Stucture factor correction
/**
 * @Refine the parameter of Babinet's procipal
 *
 * Babinet’s principle
 * F_map =  F_cur - k_solv * F_m exp (-B_solv * |s|^2 / 4) 
 * 
 * T. D. Fenn, M. J. Schnieders and A. T. Brunger
 * Acta Cryst. (2010). D66, 1024–1031
 */
int      pDenRefineSF_Babinet (PDen_t * model, PDen_t *calc, PDen_t *mask, real *k_glob, real * B_glob, real * k_solv, real * B_solv, const real prec, const size_t steps);

/**
 * apply babinet's principle to a map (mask can be the same as input)
 */
PDen_t * pDenApplySF_Babinet (PDen_t * result, PDen_t *input, PDen_t *mask, const real k_glob, const real B_glob, const real k_solv, const real B_solv);


#ifdef __cplusplus
}
#endif
#endif /* __PEN_H */
