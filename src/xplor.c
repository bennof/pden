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

int pDenReadXPLOR(PDen_t * this, const char * filename, const int mode)
{
	FILE * f;
	char buffer[BUFFERSIZE];
	int NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX;
	real a, b, c, alpha, beta, gamma; 
	size_t x,y,z,n;
	real *r, *p, h;

	f = fopen(filename,"r");
	if ( !f )
		return 1;

	if ( !fgets(buffer,BUFFERSIZE,f) ) // 1 - empty line written by the `/ ` FORTRAN format descriptor in the formatted map file
	{
		fclose(f);
		return 2;
	}

	if ( !fgets(buffer,BUFFERSIZE,f) ) // 2 - NTITLE
	{
		fclose(f);
		return 2;
	}

	if(!sscanf(buffer,"%lu!NTITLE",&n)){ //skip header
		fclose(f);
		return 2;
	}
	for(n;n>0;n--)
		if ( !fgets(buffer,BUFFERSIZE,f) ) // 3 - 
		{
			fclose(f);
			return 2;
		}
	if ( !fscanf(f,"%8i %8i %8i %8i %8i %8i %8i %8i %8i\n",&NA, &AMIN, &AMAX, &NB, &BMIN, &BMAX, &NC, &CMIN, &CMAX) ) // 6 - NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX 
	{
		fclose(f);
		return 2;
	}

	#ifdef DOUBLE
	if ( !fscanf(f,"%12.5lE %12.5lE %12.5lE %12.5lE %12.5lE %12.5lE\n",&a, &b, &c, &alpha, &beta, &gamma ) ) // 7 - a, b, c, alpha, beta, gamma
	#else
	if ( !fscanf(f,"%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",&a, &b, &c, &alpha, &beta, &gamma ) ) // 7 - a, b, c, alpha, beta, gamma
	#endif
	{
		fclose(f);
		return 2;
	}



	this->size.x = ( size_t ) ( AMAX - AMIN + 1 );
	this->size.y = ( size_t ) ( BMAX - BMIN + 1 );
	this->size.z = ( size_t ) ( CMAX - CMIN + 1 );

	this->apix.x = ( real ) ( a / NA );
	this->apix.y = ( real ) ( b / NB );
	this->apix.z = ( real ) ( c / NC );

	this->origin.x = ( real ) ( AMIN * this->apix.x );
	this->origin.y = ( real ) ( BMIN * this->apix.y );
	this->origin.z = ( real ) ( CMIN * this->apix.z );

	this->cell.x = ( size_t ) ( NA );
	this->cell.y = ( size_t ) ( NB );
	this->cell.z = ( size_t ) ( NC );

	this->start.x = ( size_t ) ( AMIN );
	this->start.y = ( size_t ) ( BMIN );
	this->start.z = ( size_t ) ( CMIN );

	if ( ( alpha != 90. ) &&
	     (  beta != 90. ) &&
	     ( gamma != 90. )){
		fclose(f);
		return 2;
	}

	if ( !fgets(buffer,BUFFERSIZE,f) ) // 8 - ZYX
	{
		fclose(f);
		return 2;
	}

	this->n = this->size.x * this->size.y * this->size.z;

	//alloc
	if ( this->data ) free( this->data );
	r = ( real * ) malloc ( ( this->size.x + 2 ) * this->size.y * this->size.z * sizeof( real ) );
	p=r;
	
	//loop reading data
	for(z=0;z<this->size.z;z++) { 
		if ( !fgets(buffer,BUFFERSIZE,f) ){ // Z - Page head
			fclose(f);
			free(r);
			return 3;
		}
		for(y=0;y<this->size.y;y++) {
			for(x=0;x<this->size.x;x++) { // X - Page record
				#ifdef DOUBLE
				if(!fscanf(f, "%12.5lE",&h))
				#else
				if(!fscanf(f, "%12.5E", &h))
				#endif
				{
					fclose(f);
					free(r);
					return 3;
				
				}
				(*p++) = h;
			}
		}
	}
	fclose(f);
	this->data=r;
	return 0;
}


int pDenWriteXPLOR(PDen_t * this, const char * filename, const int mode)
{
	FILE * f;
	char buffer[BUFFERSIZE];
	int NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX;
	real a, b, c, alpha, beta, gamma; 
	real *p;
	size_t x,y,z;


	f = fopen(filename,"r");
	if ( !f )
		return 1;
	
	fputs("\n",f);//1
	fputs("       2 !NTITLE                                                                \n",f); // 2
	fputs(" REMARKS FILENAME=''                                                            \n",f); // 3
	fputs(" REMARKS DATE:                          created by user:                        \n",f); // 4

        AMAX = this->start.x + this->size.x - 1;
	AMIN = this->start.x;
	NA   = this->cell.x;
	a    = this->cell.x * this->apix.x;

        BMAX = this->start.y + this->size.y - 1;
	BMIN = this->start.y;
	NB   = this->cell.y;
	b    = this->cell.y * this->apix.y;

        CMAX = this->start.z + this->size.z - 1;
	CMIN = this->start.z;
	NC   = this->cell.z;
	c    = this->cell.z * this->apix.z;


	fprintf(f,"%8i %8i %8i %8i %8i %8i %8i %8i %8i\n",NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX); // 6
	fprintf(f,"%12.5lE %12.5lE %12.5lE %12.5lE %12.5lE %12.5lE\n",a, b, c, 90.,90.,90. ); // 7 - a, b, c, alpha, beta, gamma
	fputs("ZYX                                                                             \n",f); // 2
	
	p = this->data;
	for(z=0;z<this->size.z;z++) { 
		fprintf(f,"%8i                                                                        ",z);
		for(y=0;y<this->size.y*this->size.x;y++) {
			if(!(y%6)) 
				fputs("\n",f);
			#ifdef DOUBLE
			fprintf(f, "%12.5lE", (*p++));
			#else
			fprintf(f, "%12.5E", (*p++));
			#endif
		}
		fputs("\n",f);
	}
	fputs("-9999\n",f);
	fprintf(f,"%12.5lE %12.5lEi\n",0.,0.);
	fclose(f);
}
