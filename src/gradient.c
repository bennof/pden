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
#include "tools.h"
#include "mathtools.h"


int pDenGetGradient(PDen_t * this,PDen_t *ref, real *grad,real *x, size_t natoms){
	size_t q;
	size_t dim3, dim2;

	vec3 dx;  //gradients
	vec3 p;   //position
	vec3 h;   //h factor
	vec3 s;   //s factor

	debug("get gradient at atomic position");


	real r, * d, * dr;

	uvec3 gridid;

        dim2 = this->size.x;
	dim3 = this->size.y*dim2;

	d  = this->data;
	dr = ref->data;

	for(q=0;q<natoms;q++) {
		s.x = 0.;
		s.y = 0.;
		s.z = 0.;
		dx.x = 0.;
		dx.y = 0.;
		dx.z = 0.;

		//local position
		#ifdef PDEN_RENDER_NM_COORDS
		p.x = (((x[3*q  ] * 10.	- this->origin.x ) / this->apix.x));
		p.y = (((x[3*q+1] * 10. - this->origin.y ) / this->apix.y)); 
		p.z = (((x[3*q+2] * 10. - this->origin.z ) / this->apix.z));
		#else
		p.x = (((x[3*q  ] - this->origin.x ) / this->apix.x));
		p.y = (((x[3*q+1] - this->origin.y ) / this->apix.y)); 
		p.z = (((x[3*q+2] - this->origin.z ) / this->apix.z));
		#endif

		//get upper left grid point
		gridid.x = (size_t)( p.x );	
		gridid.y = (size_t)( p.y );	
		gridid.z = (size_t)( p.z );	

		//get alpha
		p.x -= gridid.x;
		p.y -= gridid.y;
		p.z -= gridid.z;

		gridid.y *= dim2;
		gridid.z *= dim3;

		debug("position %lf %lf %lf \n",p.x,p.y,p.z);

		//square walk unrolled loop
		// 0,0,0
		r = dr[gridid.z     +gridid.y     +gridid.x  ]
		  - d [gridid.z     +gridid.y     +gridid.x  ];
		h.x = isqrt_((      p.y ) * (      p.y ) + (      p.z ) * (      p.z ));
		h.y = isqrt_((      p.x ) * (      p.x ) + (      p.z ) * (      p.z ));
		h.z = isqrt_((      p.x ) * (      p.x ) + (      p.y ) * (      p.y ));

//{
//		double hh = (      p.x ) * (      p.x ) + (      p.z ) * (      p.z );
//		fprintf(stderr,"%lu:1       scale (%lf %lf %lf) => %lf %lf\n",q,h.x,h.y,h.z,hh,isqrt_(hh));
//}
		dx.x -= r * h.x;
		dx.y -= r * h.y;
		dx.z -= r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 1,0,0
		r = dr[gridid.z      +gridid.y     +gridid.x+1]
		  - d [gridid.z      +gridid.y     +gridid.x+1];
		h.x = isqrt_((      p.y ) * (      p.y ) + (      p.z ) * (      p.z ));
		h.y = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + (      p.z ) * (      p.z ));
		h.z = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + (      p.y ) * (      p.y ));

		dx.x += r * h.x;
		dx.y -= r * h.y;
		dx.z -= r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 0,1,0
		r = dr[gridid.z     +gridid.y+dim2+gridid.x  ]
		  - d [gridid.z     +gridid.y+dim2+gridid.x  ];
		h.x = isqrt_(( 1. - p.y ) * ( 1. - p.y ) + (      p.z ) * (      p.z ));
		h.y = isqrt_((      p.x ) * (      p.x ) + (      p.z ) * (      p.z ));
		h.z = isqrt_((      p.x ) * (      p.x ) + ( 1. - p.y ) * ( 1. - p.y ));

		dx.x -= r * h.x;
		dx.y += r * h.y;
		dx.z -= r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 1,1,0
		r = dr[gridid.z     +gridid.y+dim2+gridid.x+1]
		  - d[gridid.z     +gridid.y+dim2+gridid.x+1];
		h.x = isqrt_(( 1. - p.y ) * ( 1. - p.y ) + (      p.z ) * (      p.z ));
		h.y = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + (      p.z ) * (      p.z ));
		h.z = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + ( 1. - p.y ) * ( 1. - p.y ));

		dx.x += r * h.x;
		dx.y += r * h.y;
		dx.z -= r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 0,0,1
		r = dr[gridid.z+dim3+gridid.y     +gridid.x  ]
		  - d [gridid.z+dim3+gridid.y     +gridid.x  ];
		h.x = isqrt_((      p.y ) * (      p.y ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.y = isqrt_((      p.x ) * (      p.x ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.z = isqrt_((      p.x ) * (      p.x ) + (      p.y ) * (      p.y ));

		dx.x -= r * h.x;
		dx.y -= r * h.y;
		dx.z += r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 1,0,1
		r = dr[gridid.z+dim3+gridid.y     +gridid.x+1]
		  - d [gridid.z+dim3+gridid.y     +gridid.x+1];
		h.x = isqrt_((      p.y ) * (      p.y ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.y = isqrt_(( 1. - p.x ) * ( 1. - p.y ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.z = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + (      p.y ) * (      p.y ));

		dx.x += r * h.x;
		dx.y -= r * h.y;
		dx.z += r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 0,1,1
		r = dr[gridid.z+dim3+gridid.y+dim2+gridid.x  ]
		  - d [gridid.z+dim3+gridid.y+dim2+gridid.x  ];
		h.x = isqrt_(( 1. - p.y ) * ( 1. - p.y ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.y = isqrt_((      p.x ) * (      p.x ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.z = isqrt_((      p.x ) * (      p.x ) + ( 1. - p.y ) * ( 1. - p.y ));

		dx.x -= r * h.x;
		dx.y += r * h.y;
		dx.z += r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 1,1,1
		r = dr[gridid.z+dim3+gridid.y+dim2+gridid.x+1]
		  - d [gridid.z+dim3+gridid.y+dim2+gridid.x+1];
		h.x = isqrt_(( 1. - p.y ) * ( 1. - p.y ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.y = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.z = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + ( 1. - p.y ) * ( 1. - p.y ));

		dx.x += r * h.x;
		dx.y += r * h.y;
		dx.z += r * h.z;

		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		dx.x =  dx.x / s.x;
		dx.y =  dx.y / s.y;
		dx.z =  dx.z / s.z;

		#ifdef PDEN_RENDER_NM_COORDS
		grad[ 3 * q + 0 ] = dx.x * .1;
		grad[ 3 * q + 1 ] = dx.y * .1;
		grad[ 3 * q + 2 ] = dx.z * .1;
		#else
		grad[ 3 * q + 0 ] = dx.x ;
		grad[ 3 * q + 1 ] = dx.y ;
		grad[ 3 * q + 2 ] = dx.z ;
		#endif
	}
	return 0;
}

int pDenGetRangedGradient(PDen_t * this,PDen_t * ref,real *grad,real *x, size_t natoms,  real width)
{
	size_t i,j,k,q;
	size_t idxx, idx;

	vec3 dx; //gradients
	vec3 p; //position
	vec3 w; //weight
	vec3 h; //h factor
	vec3 s; //h factor

	real r, * d, * dr;
	real h1, h2;

	uvec3 gridid;
	uvec3 gridstart;
	uvec3 gridend;
	uvec3 walk;

	debug("get gradient at atomic position (radius=%lf)",width);

	//walk radius
	walk.x = ((size_t)(0.5*width/this->apix.x));
	walk.y = ((size_t)(0.5*width/this->apix.y)); 
	walk.z = ((size_t)(0.5*width/this->apix.z));


	d = this->data;
	dr = ref->data;
	for(q=0;q<natoms;q++) {
		dx.x = 0.;
		dx.y = 0.;
		dx.z = 0.;
		s.x  = 0.;
		s.y  = 0.;
		s.z  = 0.;

		//local position
		#ifdef PDEN_RENDER_NM_COORDS
		p.x = (((x[3*q  ] * 10.	- this->origin.x ) / this->apix.x));
		p.y = (((x[3*q+1] * 10. - this->origin.y ) / this->apix.y)); 
		p.z = (((x[3*q+2] * 10. - this->origin.z ) / this->apix.z));
		#else
		p.x = (((x[3*q  ] - this->origin.x ) / this->apix.x));
		p.y = (((x[3*q+1] - this->origin.y ) / this->apix.y)); 
		p.z = (((x[3*q+2] - this->origin.z ) / this->apix.z));
		#endif

		//get upper left grid point
		gridid.x = (size_t)( p.x );	
		gridid.y = (size_t)( p.y );	
		gridid.z = (size_t)( p.z );	

		//end of gridwalk
		gridend.x = gridid.x + walk.x + 1;
		gridend.y = gridid.y + walk.y + 1;
		gridend.z = gridid.z + walk.z + 1;

		//start of grid walk
		gridstart.x = gridid.x - walk.x + 1;
		gridstart.y = gridid.y - walk.y + 1;
		gridstart.z = gridid.z - walk.z + 1;

		
		for(k = gridstart.z;
		    k < gridend.z;
		    k++){
			idxx = k * this->size.x * this->size.y;
			for(j = gridstart.y;
			    j < gridend.y;
			    j++){
				idx = idxx + j * this->size.x;
				for(i = gridstart.x;
				    i < gridend.x;
				    i++){

					/*
					 * directional linear weighting
					 * using f'(x)= lim (f(x+h)-f(x)/h) = lim (f(x+a*h)-f(x-(1-a)*h)/h) with a in [0,1]
					 * for more than two points:
					 * f'(x) ~ sum ( f(x+ah)/h )  
					 * 
					 * using a gaussian for mixing directions 
					 * for x direction the kernel is:
					 * exp(r^2/x^2) <= 1.0
					 */
					w.x  = ((real)i-p.x);	
					w.y  = ((real)j-p.y);	
					w.z  = ((real)k-p.z);

					h.x  = ( w.x<0 )? -1. : 1.;  
					h.y  = ( w.y<0 )? -1. : 1.;  
					h.z  = ( w.z<0 )? -1. : 1.;  

					w.x  = ( w.x * w.x );
					w.y  = ( w.y * w.y );
					w.z  = ( w.z * w.z );

					r    = w.x+w.y+w.z;

					w.x  = ( w.x ) / r;
					w.y  = ( w.y ) / r;
					w.z  = ( w.z ) / r;

					s.x += w.x;
					s.y += w.y;
					s.z += w.z;

					h1 = dr[idx];
					h2 = d[idx];
					h1 -= h2;
					dx.x += h.x * w.x * (h1);
					dx.y += h.x * w.y * (h1);
					dx.z += h.x * w.z * (h1);

					idx++;
				}
			}
		}
		
		#ifdef PDEN_RENDER_NM_COORDS	
		grad[ 3 * q + 0 ] = 2. * dx.x / s.x * .1;
		grad[ 3 * q + 1 ] = 2. * dx.y / s.y * .1;
		grad[ 3 * q + 2 ] = 2. * dx.z / s.z * .1;
		#else
		grad[ 3 * q + 0 ] = 2. * dx.x / s.x;
		grad[ 3 * q + 1 ] = 2. * dx.y / s.y;
		grad[ 3 * q + 2 ] = 2. * dx.z / s.z;
		#endif
	}
	return 0;
}

int pDenAddGradient(PDen_t * this, PDen_t * ref, real *x, size_t natoms, real weight){
	size_t q;
	size_t dim3, dim2;

	vec3 dx; //gradients
	vec3 p;   //position
	vec3 h;   //h factor
	vec3 s;   //s factor


	real r, * d, * dr;

	uvec3 gridid;
	debug("add gradient at atomic position to atom");

        dim2 = this->size.x;
	dim3 = this->size.y*dim2;

	d  = this->data;
	dr = ref->data;

	for(q=0;q<natoms;q++) {
		s.x = 0.;
		s.y = 0.;
		s.z = 0.;
		dx.x = 0.;
		dx.y = 0.;
		dx.z = 0.;

		//local position
		#ifdef PDEN_RENDER_NM_COORDS
		p.x = (((x[3*q  ] * 10.	- this->origin.x ) / this->apix.x));
		p.y = (((x[3*q+1] * 10. - this->origin.y ) / this->apix.y)); 
		p.z = (((x[3*q+2] * 10. - this->origin.z ) / this->apix.z));
		#else
		p.x = (((x[3*q  ] - this->origin.x ) / this->apix.x));
		p.y = (((x[3*q+1] - this->origin.y ) / this->apix.y)); 
		p.z = (((x[3*q+2] - this->origin.z ) / this->apix.z));
		#endif

		//get upper left grid point
		gridid.x = (size_t)( p.x );	
		gridid.y = (size_t)( p.y );	
		gridid.z = (size_t)( p.z );	

		//get alpha
		p.x -= gridid.x;
		p.y -= gridid.y;
		p.z -= gridid.z;

		gridid.y *= dim2;
		gridid.z *= dim3;


		//square walk unrolled loop
		// 0,0,0
		r = dr[gridid.z     +gridid.y     +gridid.x  ]
		  - d [gridid.z     +gridid.y     +gridid.x  ];
		h.x = isqrt_((      p.y ) * (      p.y ) + (      p.z ) * (      p.z ));
		h.y = isqrt_((      p.x ) * (      p.x ) + (      p.z ) * (      p.z ));
		h.z = isqrt_((      p.x ) * (      p.x ) + (      p.y ) * (      p.y ));

		dx.x -= r * h.x;
		dx.y -= r * h.y;
		dx.z -= r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 1,0,0
		r = dr[gridid.z      +gridid.y     +gridid.x+1]
		  - d [gridid.z      +gridid.y     +gridid.x+1];
		h.x = isqrt_((      p.y ) * (      p.y ) + (      p.z ) * (      p.z ));
		h.y = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + (      p.z ) * (      p.z ));
		h.z = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + (      p.y ) * (      p.y ));

		dx.x += r * h.x;
		dx.y -= r * h.y;
		dx.z -= r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 0,1,0
		r = dr[gridid.z     +gridid.y+dim2+gridid.x  ]
		  - d [gridid.z     +gridid.y+dim2+gridid.x  ];
		h.x = isqrt_(( 1. - p.y ) * ( 1. - p.y ) + (      p.z ) * (      p.z ));
		h.y = isqrt_((      p.x ) * (      p.x ) + (      p.z ) * (      p.z ));
		h.z = isqrt_((      p.x ) * (      p.x ) + ( 1. - p.y ) * ( 1. - p.y ));

		dx.x -= r * h.x;
		dx.y += r * h.y;
		dx.z -= r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 1,1,0
		r = dr[gridid.z     +gridid.y+dim2+gridid.x+1]
		  - d[gridid.z     +gridid.y+dim2+gridid.x+1];
		h.x = isqrt_(( 1. - p.y ) * ( 1. - p.y ) + (      p.z ) * (      p.z ));
		h.y = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + (      p.z ) * (      p.z ));
		h.z = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + ( 1. - p.y ) * ( 1. - p.y ));

		dx.x += r * h.x;
		dx.y += r * h.y;
		dx.z -= r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 0,0,1
		r = dr[gridid.z+dim3+gridid.y     +gridid.x  ]
		  - d [gridid.z+dim3+gridid.y     +gridid.x  ];
		h.x = isqrt_((      p.y ) * (      p.y ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.y = isqrt_((      p.x ) * (      p.x ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.z = isqrt_((      p.x ) * (      p.x ) + (      p.y ) * (      p.y ));

		dx.x -= r * h.x;
		dx.y -= r * h.y;
		dx.z += r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 1,0,1
		r = dr[gridid.z+dim3+gridid.y     +gridid.x+1]
		  - d [gridid.z+dim3+gridid.y     +gridid.x+1];
		h.x = isqrt_((      p.y ) * (      p.y ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.y = isqrt_(( 1. - p.x ) * ( 1. - p.y ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.z = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + (      p.y ) * (      p.y ));

		dx.x += r * h.x;
		dx.y -= r * h.y;
		dx.z += r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 0,1,1
		r = dr[gridid.z+dim3+gridid.y+dim2+gridid.x  ]
		  - d [gridid.z+dim3+gridid.y+dim2+gridid.x  ];
		h.x = isqrt_(( 1. - p.y ) * ( 1. - p.y ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.y = isqrt_((      p.x ) * (      p.x ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.z = isqrt_((      p.x ) * (      p.x ) + ( 1. - p.y ) * ( 1. - p.y ));

		dx.x -= r * h.x;
		dx.y += r * h.y;
		dx.z += r * h.z;
		
		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		// 1,1,1
		r = dr[gridid.z+dim3+gridid.y+dim2+gridid.x+1]
		  - d [gridid.z+dim3+gridid.y+dim2+gridid.x+1];
		h.x = isqrt_(( 1. - p.y ) * ( 1. - p.y ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.y = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + ( 1. - p.z ) * ( 1. - p.z ));
		h.z = isqrt_(( 1. - p.x ) * ( 1. - p.x ) + ( 1. - p.y ) * ( 1. - p.y ));

		dx.x += r * h.x;
		dx.y += r * h.y;
		dx.z += r * h.z;

		s.x += h.x;
		s.y += h.y;
		s.z += h.z;

		dx.x =  weight * dx.x / s.x;
		dx.y =  weight * dx.y / s.y;
		dx.z =  weight * dx.z / s.z;

		#ifdef PDEN_RENDER_NM_COORDS	
		x[ 3 * q + 0 ] += dx.x * .1;
		x[ 3 * q + 1 ] += dx.y * .1;
		x[ 3 * q + 2 ] += dx.z * .1;
		#else
		x[ 3 * q + 0 ] += dx.x;
		x[ 3 * q + 1 ] += dx.y;
		x[ 3 * q + 2 ] += dx.z;
		#endif
	}
	return 0;
}

int pDenAddRangedGradient(PDen_t * this,PDen_t * ref,real *x, size_t natoms,  real width, real weight)
{
	size_t i,j,k,q;
	size_t idxx, idx;

	vec3 dx; //gradients
	vec3 p; //position
	vec3 w; //weight
	vec3 h; //h factor
	vec3 s; //h factor

	real r, * d, *dr;
	register real h1,h2;

	uvec3 gridid;
	uvec3 gridstart;
	uvec3 gridend;
	uvec3 walk;

	debug("add gradient at atomic position to atom (radius=%lf)",width);

	//walk radius
	walk.x = ((size_t)(0.5*width/this->apix.x))+1;
	walk.y = ((size_t)(0.5*width/this->apix.y))+1; 
	walk.z = ((size_t)(0.5*width/this->apix.z))+1;

	d = this->data;
	dr = ref->data;
	for(q=0;q<natoms;q++) {
		dx.x = 0.;
		dx.y = 0.;
		dx.z = 0.;
		s.x = 0.;
		s.y = 0.;
		s.z = 0.;

		//local position
		#ifdef PDEN_RENDER_NM_COORDS
		p.x = (((x[3*q  ] * 10.	- this->origin.x ) / this->apix.x));
		p.y = (((x[3*q+1] * 10. - this->origin.y ) / this->apix.y)); 
		p.z = (((x[3*q+2] * 10. - this->origin.z ) / this->apix.z));
		#else
		p.x = (((x[3*q  ] - this->origin.x ) / this->apix.x));
		p.y = (((x[3*q+1] - this->origin.y ) / this->apix.y)); 
		p.z = (((x[3*q+2] - this->origin.z ) / this->apix.z));
		#endif

		//get upper left grid point
		gridid.x = (size_t)( p.x );	
		gridid.y = (size_t)( p.y );	
		gridid.z = (size_t)( p.z );	

		//end of gridwalk
		gridend.x = gridid.x + walk.x + 1;
		gridend.y = gridid.y + walk.y + 1;
		gridend.z = gridid.z + walk.z + 1;

		//start of grid walk
		gridstart.x = gridid.x - walk.x + 1;
		gridstart.y = gridid.y - walk.y + 1;
		gridstart.z = gridid.z - walk.z + 1;

		debug("atom %lu: pos (%lf %lf %lf) \n",q,p.x,p.y,p.z);
	
		width = 1./width;
		
		for(k = gridstart.z;
		    k < gridend.z;
		    k++){
			idxx = k * this->size.x * this->size.y;
			for(j = gridstart.y;
			    j < gridend.y;
			    j++){
				idx = idxx + j * this->size.x;
				for(i = gridstart.x;
				    i < gridend.x;
				    i++){

					/*
					 * directional linear weighting
					 * using f'(x)= lim (f(x+h)-f(x)/h) = lim (f(x+a*h)-f(x-(1-a)*h)/h) with a in [0,1]
					 * for more than two points:
					 * f'(x) ~ sum ( f(x+ah)/h )  
					 * 
					 * using a gaussian for mixing directions (not working)
					 * for x direction the kernel uses a radial weight:
					 * x^2/r^2
					 */
					w.x  = ((real)i-p.x);	
					w.y  = ((real)j-p.y);	
					w.z  = ((real)k-p.z);

					h.x  = ( w.x < 0.) ? -1. : 1.;
					h.y  = ( w.y < 0.) ? -1. : 1.;
					h.z  = ( w.z < 0.) ? -1. : 1.;

					w.x  = ( w.x * w.x );
					w.y  = ( w.y * w.y );
					w.z  = ( w.z * w.z );

					r    = w.x+w.y+w.z;

					w.x  = ( w.x ) / r;
					w.y  = ( w.y ) / r;
					w.z  = ( w.z ) / r;

					s.x += (w.x);
					s.y += (w.y);
					s.z += (w.z);

					h1 = dr[idx];
					h2 = d[idx];
					h1 -= h2;
					dx.x += h.x * w.x * (h1);
					dx.y += h.x * w.y * (h1);
					dx.z += h.x * w.z * (h1);

					idx++;
				}
			}
		}

		dx.x =  weight * dx.x / s.x * width;
		dx.y =  weight * dx.y / s.y * width;
		dx.z =  weight * dx.z / s.z * width;

		debug("atom %lu: grad = (%lf %lf %lf) \n",q,dx.x,dx.y,dx.z);

		#ifdef PDEN_RENDER_NM_COORDS	
		x[ 3 * q + 0 ] += dx.x * .1;
		x[ 3 * q + 1 ] += dx.y * .1;
		x[ 3 * q + 2 ] += dx.z * .1;
		#else
		x[ 3 * q + 0 ] += dx.x;
		x[ 3 * q + 1 ] += dx.y;
		x[ 3 * q + 2 ] += dx.z;
		#endif
	}
	return 0;
}
