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

#include "pden.in.h"
#include "tools.h"

#define BUFFERSIZE 120

int pDenReadXPLOR(PDen_t * this, const char * filename, int mode)
{
	FILE * f;
	char buffer[BUFFERSIZE];
	int NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX;
	real a, b, c, alpha, beta, gamma; 
	size_t x,y,z;
	real *r, *p, h;

	f = fopen(filename,"r");
	if ( !f )
		return 1;

	if ( !fgets(buffer,BUFFERSIZE,f) ) // 1 - empty line written by the `/ ` FORTRAN format descriptor in the formatted map file
		return 2;
	
	if ( !fgets(buffer,BUFFERSIZE,f) ) // 2 - 
		return 2;

	if ( !fgets(buffer,BUFFERSIZE,f) ) // 3 - 
		return 2;

	if ( !fgets(buffer,BUFFERSIZE,f) ) // 4 - 
		return 2;

	if ( !fgets(buffer,BUFFERSIZE,f) ) // 5 - 
		return 2;

	if ( !fscanf(f,"%i %i %i %i %i %i %i %i %i\n",&NA, &AMIN, &AMAX, &NB, &BMIN, &BMAX, &NC, &CMIN, &CMAX) ) // 6 - NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX 
		return 2;
	#ifdef DOUBLE
	if ( !fscanf(f,"%lf %lf %lf %lf %lf %lf\n",&a, &b, &c, &alpha, &beta, &gamma ) ) // 7 - a, b, c, alpha, beta, gamma
	#else
	if ( !fscanf(f,"%f %f %f %f %f %f\n",&a, &b, &c, &alpha, &beta, &gamma ) ) // 7 - a, b, c, alpha, beta, gamma
	#endif
		return 2;


	this->size.x = ( size_t ) ( AMAX - AMIN + 1 );
	this->size.y = ( size_t ) ( BMAX - BMIN + 1 );
	this->size.z = ( size_t ) ( CMAX - CMIN + 1 );

	this->apix.x = ( real ) ( a / NA );
	this->apix.y = ( real ) ( b / NB );
	this->apix.z = ( real ) ( c / NC );

	this->origin.x = ( real ) ( AMIN * this->apix.x );
	this->origin.y = ( real ) ( BMIN * this->apix.y );
	this->origin.z = ( real ) ( CMIN * this->apix.z );

	this->cell.x = ( size_t ) ( AMAX );
	this->cell.y = ( size_t ) ( BMAX );
	this->cell.z = ( size_t ) ( CMAX );

	this->start.x = ( size_t ) ( AMIN );
	this->start.y = ( size_t ) ( BMIN );
	this->start.z = ( size_t ) ( CMIN );

	if ( ( alpha != 90. ) &&
	     (  beta != 90. ) &&
	     ( gamma != 90. ))
		return 2;

	if ( !fgets(buffer,BUFFERSIZE,f) ) // 8 - ZYX
		return 2;
	

	this->n = this->size.x * this->size.y * this->size.z;

	//alloc
	if ( this->data ) free( this->data );
	r = ( real * ) malloc ( ( this->size.x + 2 ) * this->size.y * this->size.z * sizeof( real ) );
	p=r;
	
	//loop reading data
	for(z=0;z<this->size.z;z++) { 
		if ( !fgets(buffer,BUFFERSIZE,f) ){ // Z - Page head
			free(r);
			return 3;
		}
		for(y=0;y<this->size.y;y++) {
			if ( !fgets(buffer,BUFFERSIZE,f) ){ // Y - Page head
				free(r);
				return 3;
			}
			for(x=0;x<this->size.x;x++) { // X - Page record
				#ifdef DOUBLE
				fscanf(f, "%lf", &h);
				#else
				fscanf(f, "%f", &h);
				#endif
				(*p++) = h;
			}
		}
	}
	fclose(f);
	this->data=r;
	return 0;
}
