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

#include "pden.h"
#include "mathtools.h"

real randf()
{
	return ( ( real ) rand() ) / ( real ) RAND_MAX;
}

real gaussRand ()
{
	real q, p;
	static real r1 = 0.0, r2 = 0.0;
	
	if( r2 != 0.0 ) {
		r2 = 0.0;
		return r1;
        }

	do {
		r1 = 2.0 * randf() - 1;
		r2 = 2.0 * randf() - 1;
		q = r1 * r1 + r2 * r2;
	} while (q == 0.0 || q >= 1);	
	p = sqrt(-2 * log(q) / q);
	r1 *= p;
	r2 *= p;
	return r2;
}


#define binop( sign ) \
{\
	size_t i;\
	for(i=0;i<result->n;i++)\
		result->data[i] = a->data[i] sign b->data[i];\
	return result;\
}


#define singularfunc( func ) \
{\
	size_t i;\
	for(i=0;i<result->n;i++)\
		result->data[i] = func ( a->data[i] );\
	return result;\
}



PDen_t * pDenAdd (PDen_t * result, PDen_t * a, PDen_t * b )
binop( + )

PDen_t * pDenSub (PDen_t * result, PDen_t * a, PDen_t * b )
binop( - )

PDen_t * pDenMult (PDen_t * result, PDen_t * a, PDen_t * b )
binop( * )

PDen_t * pDenDiv (PDen_t * result, PDen_t * a, PDen_t * b )
binop( / )

PDen_t * pDenScale (PDen_t * result, PDen_t * a, real s )
singularfunc( s * )

PDen_t * pDenAddNoise (PDen_t * result, PDen_t * a, real s )
singularfunc( s * (((real)rand()/((real)RAND_MAX)*2)-1.) + )

PDen_t * pDenAddGaussNoise (PDen_t * result, PDen_t * a, real s )
singularfunc( s *  gaussRand() + ) 

PDen_t * pDenSqrt (PDen_t * result, PDen_t * a )
singularfunc( sqrt )

PDen_t * pDenISqrt (PDen_t * result, PDen_t * a )
singularfunc( isqrt_ )

PDen_t * pDenMkMask (PDen_t * result, PDen_t * a, real level )
{
	size_t i;
	for(i=0;i<result->n;i++)
		result->data[i] = (level <  a->data[i] ) ? 1. : 0.;
	return result;
}


