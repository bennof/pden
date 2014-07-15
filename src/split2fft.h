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

#ifndef __SPLIT2FFT_H
#define __SPLIT2FFT_H     0.2
#ifdef __cplusplus 
extern "C" 
#endif

/**
 * @file split2fft.h
 * @author Benjamin Falkner
 * @brief File containing a system to handle fouriertransforms on densities
 * 2nd version of splitfft 
 * uses an additional buffer to reduce alloctions
 * is not threadsafe
 * 
 * @warning internal file - not installed
 */

#include <fftw3.h> 
#include "pden.in.h"

#ifdef DOUBLE
#define FFTW_plan  fftw_plan             
#define FFTW_destroy_plan fftw_destroy_plan
#define dft_r2c    fftw_plan_dft_r2c_3d  
#define dft_c2r    fftw_plan_dft_c2r_3d  
#define ex_dft_r2c fftw_execute_dft_r2c  
#define ex_dft_c2r fftw_execute_dft_c2r  
#define complx     fftw_complex
#else              
#define FFTW_plan  fftwf_plan            
#define FFTW_destroy_plan fftwf_destroy_plan
#define dft_r2c    fftwf_plan_dft_r2c_3d 
#define dft_c2r    fftwf_plan_dft_c2r_3d 
#define ex_dft_r2c fftwf_execute_dft_r2c 
#define ex_dft_c2r fftwf_execute_dft_c2r 
#define complx     fftwf_complex
#endif



/**
 * @brief 3 split2fft structure for fft  
 *
 * split2fft data structure used for all maps of the same size (common case)
 * to reduce memory overhead and allocation
 * object used reference counting to reduce further memory 
 * please use always constuctors and deconstructors
 */
struct Split2FFT_struct {
	unsigned int cinst; /**< reference counter */
	FFTW_plan planf; /**< FFTW PLAN forward */
	FFTW_plan planb; /**< FFTW PLAN backward */
	real *switching; /**< hidden memory for swapping buffers */
	VEC3(real) mod;  /**< storing Fourie dimensions */
};

// Constructor & Destructor
/**
 * @brief split2fft constructor
 *
 * constructs a fft from an exsisting density or null
 * @param map a map of similar size or null
 * @param optimize optimization flag for FFTW (0 = no, 1 = yes) 
 * @return a reference of a split2fft object
 */
Split2FFT_t * split2FFTNew ( PDen_t * map, const int optimize);

/**
 * @brief split2fft destructor
 *
 * reduces the referenccounter of a split2fft object and deletes the object if the reference is not used anymore
 * @param this a split2fft reference
 * @return empty (null) reference
 */
Split2FFT_t * split2FFTDelete ( Split2FFT_t * this );

// Functions
/**
  * @brief execute the fft ()
  *
  * performs a fourier transform and applies normalisation only after the backtransform
  * @param this a split2fft reference
  * @param map a density map
  * @param direction (>0 forward; <0 backward)
  * @return 0
  * @deprecated this method is deprecated - use split2FFTExecuteN
  */
int split2FFTExecute ( Split2FFT_t * this, PDen_t *map, const  int direction);

/**
  * @brief execute the fft ()
  *
  * performs a fourier transform and applies normalisation 1/sqrt(N)
  * @param this a split2fft reference
  * @param map a density map
  * @param direction (>0 forward; <0 backward)
  * @return 0
  */
int split2FFTExecuteN ( Split2FFT_t * this, PDen_t *map, int direction);




#endif /* __SPLIT2FFT_H */
