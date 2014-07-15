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

#ifndef __PDEN_SFREFINE_H
#define __PDEN_SFREFINE_H
/**
 * @file pden_sfrefine.h
 * @author Benjamin Falkner
 * @brief File containing conjugated gradient methods
 *
 * This Header is only to implement new Structur Factor Optimizations
 * - typically least squares is used in this library
 * - libraray internal functions an stuff
 * - as an example see Babinet implementation
 */

#include "pden.h"
#include "tools.h"
#include "types.h"
#include "mathtools.h"

/**
 * @brief single structur factor adjustment function
 * @param s2 radius squared
 * @param calc calculted strucure factors
 * @param mask a solvent mask structure factor
 * @param x0 array of parameters 
 * @return new stucturefactor
 * @deprecated this method is deprecated to faster implementations of the algoritm
 */
typedef real (*SFfunc_t) (real s2, real calc, real mask, real x0[]);//deprecated



/**
 * @brief least square function to be optimized
 * @param X input arrays of structure factors
 * @param param0 array of parameters 
 * @return least square sum
 */
typedef real (*SFValueFunc_t) (PDen_t *X[], real param0[]);

/**
 * @brief least square funtion gradient
 * @param X input arrays of structure factors
 * @param param0 array of parameters 
 * @param grad the calculated gradient is written into this array
 */
typedef void (*SFGradFunc_t)  (PDen_t *X[], real param0[],real grad[]);


/**
 * @brief conjugated Gradients refinment of a structur factor function
 * @param X input arrays of structure factors
 * @param param0 array of parameters (input and output)
 * @param func least square function to be optimized
 * @param least square funtion gradient
 * @param number of parameters to be optimized
 * @param prec precision to stop calculation
 * @param steps max number of steps to be performed
 * @return 0
 */
int CGRefine (PDen_t *X[], real param0[], SFValueFunc_t func, SFGradFunc_t gradfunc, size_t npar, real prec, size_t steps);

#endif
