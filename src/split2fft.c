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

#include "tools.h"
#include "split2fft.h"
#include "mathtools.h"


// dublicate a exsiting split2fft used in new function
static Split2FFT_t * dub(Split2FFT_t * d){
	d->cinst++;
	return d;
}


Split2FFT_t * split2FFTNew ( PDen_t * map, const int optimize){
	if(map->fft) return dub(map->fft);
	else{
		Split2FFT_t * r;
		r = ( Split2FFT_t * ) malloc (sizeof(Split2FFT_t));

		r->cinst=1;
		
		r->mod.x = map->apix.x * map->apix.x;
		r->mod.y = map->apix.y * map->apix.y;
		r->mod.z = map->apix.z * map->apix.z;

		r->switching = ( real * ) malloc ( ( map->size.x + 2 ) * map->size.y * map->size.z * sizeof( real ) );

		if(optimize){
			real * hh = (real*)malloc(( map->size.x + 2 ) * map->size.y * map->size.z * sizeof( real ));
			info("Optimizing FFT");
			
			r->planf = dft_r2c(map->size.z,map->size.y,map->size.x,
					hh,(complx*)r->switching,FFTW_PATIENT);
			r->planb = dft_c2r(map->size.z,map->size.y,map->size.x,
					(complx*)hh,r->switching,FFTW_PATIENT);
			free(hh);
		}
		else {
			r->planf = dft_r2c(map->size.z,map->size.y,map->size.x,
					map->data,(complx*)r->switching,FFTW_ESTIMATE);
			r->planb = dft_c2r(map->size.z,map->size.y,map->size.x,
					(complx*)map->data,r->switching,FFTW_ESTIMATE);
		}
		return r;
	}
}


Split2FFT_t * split2FFTDelete ( Split2FFT_t * this )
{
	this->cinst--;
	if(!this->cinst){
		FFTW_destroy_plan(this->planf);
		FFTW_destroy_plan(this->planb);
		free(this->switching);	
		free(this);
	}
	return 0;
}

int split2FFTExecute(Split2FFT_t * this, PDen_t *map, const int direction)
{
	size_t i;
	real *h,N;
	if ( direction > 0 )
	{
		if( !(map->mode & PDEN_MODE_PHASE_SPACE) ) {
			ex_dft_r2c ( ( const FFTW_plan ) this->planf, map->data, (complx*) this->switching);
			h = map->data;
			map->data = this->switching;
			this->switching = h;
			map->mode |= PDEN_MODE_PHASE_SPACE;
		}

	}
	else if ( direction < 0 )
	{
		if( ( map->mode & PDEN_MODE_PHASE_SPACE ) ) {
			ex_dft_c2r ( ( const FFTW_plan ) this->planb, (complx*) map->data, this->switching);
			h = map->data;
			map->data = this->switching;
			this->switching = h;
			map->mode &= ~PDEN_MODE_PHASE_SPACE;
			N = 1. / ( ( real ) map->n);
			for(i=0;i<map->n;i++)
				map->data[i]*=N;
		}

	
	}
	return 0;
}

int split2FFTExecuteN(Split2FFT_t * this, PDen_t *map, const int direction)
{
	size_t i,n;
	real *h,N;
	if ( direction > 0 )
	{
		if( !(map->mode & PDEN_MODE_PHASE_SPACE) ) {
			ex_dft_r2c ( ( const FFTW_plan ) this->planf, map->data, (complx*) this->switching);
			h = map->data;
			map->data = this->switching;
			this->switching = h;
			map->mode |= PDEN_MODE_PHASE_SPACE;
			N = 1. / sqrt( ( real ) map->n);
			n = (map->size.x+2)*map->size.y*map->size.z;
			debug("N=%lf n=%lu",N,n);
			for(i=0;i<n;i++)
				map->data[i]*=N;
		}

	}
	else if ( direction < 0 )
	{
		if( ( map->mode & PDEN_MODE_PHASE_SPACE ) ) {
			ex_dft_c2r ( ( const FFTW_plan ) this->planb, (complx*) map->data, this->switching);
			h = map->data;
			map->data = this->switching;
			this->switching = h;
			map->mode &= ~PDEN_MODE_PHASE_SPACE;
			N = 1. / sqrt( ( real ) map->n);
			n = (map->size.x+2)*map->size.y*map->size.z;
			debug("N=%lf n=%lu",N,n);
			for(i=0;i<map->n;i++)
				map->data[i]*=N;
		}
	}
	return 0;
}



