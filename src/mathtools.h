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

#ifndef __MATHTOOLS_H
#define __MATHTOOLS_H

/**
 * @file mathtools.h
 * @author Benjamin Falkner
 * @brief File containing useful math functions 
 * @warning internal file - not installed
 */

#include "types.h"
#include <math.h>

#define TWO_PI 6.2831853071795864769252866
#define PI     3.1415926535897932384626434

static real isqrt_( real number )
{
#ifdef FLOAT
	//long i;
	float x2;
	const float threehalfs = 1.5;
	union {
		float f;
		long  l;
	} y;
 
	x2 = number * 0.5F;
	y.f  = number;
	//i  = * ( long * ) &y;                      
	y.l  = 0x5f3759df - ( y.l >> 1 );               
	//y  = * ( float * ) &i;
	y.f  = y.f * ( threehalfs - ( x2 * y.f * y.f ) );   
	return y.f;
#else // double
	double x2;
	const double threehalfs = 1.5;
	union {
		double  f;
		long long  l;
	} y;

	x2 = number * 0.5;
	y.f  = number;
	//i  = * ( long long * ) &y; 
    	y.l  = 0x5fe6ec85e7de30daLL - ( y.l >> 1 );
    	//y  = * ( double * ) &i;
	y.f = y.f * ( threehalfs - ( x2 * y.f * y.f ) );
    	return y.f;
#endif
}

static inline int dft_c(int x,int NX){
  return x-(NX*(x/(NX/2+1)));
}


static inline real abscplx(cplx c)
{
	register real h;
	h = c.r*c.r+c.i*c.i;
	return h * isqrt_(h);
}



#endif
