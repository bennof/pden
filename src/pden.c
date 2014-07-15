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
#include "mathtools.h"

#include "split2fft.h"


PDen_t *pDenNew()
{
	PDen_t * r;
	r = (PDen_t*)calloc(1,sizeof(PDen_t));
	return r;
}

PDen_t *pDenNewFrom( PDen_t * reference )
{
	PDen_t * r;
	r = (PDen_t*)calloc(1,sizeof(PDen_t));
	if ( !r )
		return 0;
	memcpy( r, reference, sizeof(PDen_t));
	r->fft=split2FFTNew(reference,0);
	r->data = 0;
	return r;
}

PDen_t *pDenNewE(size_t n[],real apix[], real origin[])
{
	PDen_t * r;
	r = (PDen_t*)calloc(1,sizeof(PDen_t));
	if ( !r )
		return 0;

	r->size.x = n[0];
	r->size.y = n[1];
	r->size.z = n[2];

	r->apix.x = apix[0];
	r->apix.y = apix[1];
	r->apix.z = apix[2];

	r->origin.x = origin[0];
	r->origin.y = origin[1];
	r->origin.z = origin[2];

	r->n = r->size.x * r->size.y * r->size.z;

	r->data=0;

	return r;
}

PDen_t *pDenNewESimple(real min[], real max[],real apix, real offset)
{
	PDen_t * r;
	real h,dx,dy,dz;
	r = (PDen_t*)calloc(1,sizeof(PDen_t));
	if ( !r )
		return 0;

	//sym 
	r->apix.x = apix;
	r->apix.y = apix;
	r->apix.z = apix;

	//cubic
	dx = max[0] - min[0];
	dy = max[1] - min[1];
	h = (dx>dy)? dx : dy;
	dz = max[2] - min[2];
	if(dz>h) h=dz;

	//pos grid
	min[0] -= 0.5*( h - dx );
	min[1] -= 0.5*( h - dy );
	min[2] -= 0.5*( h - dz );

	max[0] += 0.5*( h - dx );
	max[1] += 0.5*( h - dy );
	max[2] += 0.5*( h - dz );

	r->start.x = ( size_t ) floor ( (min[0] - offset ) / r->apix.x);
	r->start.y = ( size_t ) floor ( (min[1] - offset ) / r->apix.y); 
	r->start.z = ( size_t ) floor ( (min[2] - offset ) / r->apix.z); 

	r->cell.x = ( ( size_t ) ceil ( (max[0] + offset ) / r->apix.x) );
	r->cell.y = ( ( size_t ) ceil ( (max[1] + offset ) / r->apix.y) );
	r->cell.z = ( ( size_t ) ceil ( (max[2] + offset ) / r->apix.z) );

	r->size.x = r->cell.x - r->start.x + 1;
	r->size.y = r->cell.y - r->start.y + 1;
	r->size.z = r->cell.z - r->start.z + 1;

	r->origin.x = r->start.x * r->apix.x;
	r->origin.y = r->start.y * r->apix.y;
	r->origin.z = r->start.z * r->apix.z;

	r->n = r->size.x * r->size.y * r->size.z;

	return r;
}

PDen_t * pDenNewFromMol( real coords[], size_t natoms, real apix, real offset)
{
	size_t i;
	real min[3] = {0., 0., 0.};
	real max[3] = {0., 0., 0.};

	real *c = coords;

	for (i=0;i<natoms;i++){
		if(min[0] > *c) min[0]=*c;
		if(max[0] < *c) max[0]=*c;
		c++;
		if(min[1] > *c) min[1]=*c;
		if(max[1] < *c) max[1]=*c;
		c++;
		if(min[2] > *c) min[2]=*c;
		if(max[2] < *c) max[2]=*c;
		c++;
	}
	return pDenNewESimple( min, max, apix, offset);
}

int pDenAlloc( PDen_t * this )
{	
	this->data = ( real * ) malloc( ( this->size.x + 2 ) * this->size.y * this->size.z * sizeof( real ) );
	if(!this->data)
		return 1;
	return 0;
}

PDen_t *pDenDelete ( PDen_t * this )
{
	if(this->data) free(this->data);
	free(this);
	return 0;
}

int pDenClear(PDen_t * this)
{
	this->mode=0;
	if(this->data)
		memset(this->data, 0, this->n * sizeof( real ) );
	return 0;
}

real *pDenGetData( PDen_t * this )
{
	return this->data;
}


void pDenPrint (PDen_t * this, FILE * stream)
{
	fprintf(stream,
		"\nDensity Map Info (%s):\n"
		"   Grid ........... : %7lu   %7lu   %7lu\n"
		"   Origin ......... : %7ld   %7ld   %7ld\n"
		"   Extent ......... : %7lu   %7lu   %7lu\n"
		"   Unit Cell ...... : %7.3lf   %7.3lf   %7.3lf\n"
		"   Grid Spacing ... : %7.3lf   %7.3lf   %7.3lf\n"
		"   Cartesian Origin : %7.3lf   %7.3lf   %7.3lf\n",
		(this->mode & PDEN_MODE_PHASE_SPACE) ? "Fourier" : "Real",
		this->cell.x,this->cell.y,this->cell.z,
		this->start.x,this->start.y,this->start.z,
		this->size.x,this->size.y,this->size.z,
		this->apix.x*this->cell.x,this->apix.y*this->cell.y,this->apix.z*this->cell.z,
		this->apix.x,this->apix.y,this->apix.z,
		this->origin.x,this->origin.y,this->origin.z
	);

}


int pDenMeanSd(PDen_t * this,real *mean, real *sd)
{
	register size_t i;
	register real m,v,h;
	m=0.;
	v=0.;

	for(i=0;i<this->n;i++){
		h=this->data[i];
		m+=h;
		v+=h*h;
	}

	m/=this->n;
	v/=this->n;
	v -= m*m;
	*mean = m;
	*sd = v*isqrt_(v);
	return 0;
}

real pDenCorr(PDen_t * map1, PDen_t * map2)
{
	size_t i,n;
	real *X,*Y;
	register real  Xi,Yi;
	real mx,my,sxx,syy,sxy;
	mx = my = sxx = syy = sxy = 0.;

	X = map1->data;
	Y = map2->data;
	n = map1->n;

	for( i = 0; i < n; i++){
		Xi   = X[i];
		Yi   = Y[i];
		mx  += Xi; 
		sxx += Xi*Xi;
		my  += Yi; 
		syy += Yi*Yi;
		sxy += Xi*Yi;
	}

	mx  /= n;
	my  /= n;
	sxx /= n;
	syy /= n;
	sxy /= n;
	
	sxx -= mx*mx;
	syy -= my*my;

	sxy -= mx*my;

	return sxy * isqrt_(sxx) * isqrt_(syy);
}

int pDenMapCompare(PDen_t * map1, PDen_t * map2, real *corr, real *rsr)
{
	size_t i,n;
	real *X,*Y;
	register real  Xi,Yi;
	real mx,my,sxx,syy,sxy,sumxy,difxy;

	mx = my = sxx = syy = sxy = 0.;
	sumxy=difxy=0.;

	X = map1->data;
	Y = map2->data;
	n = map1->n;

	for( i = 0; i < n; i++){
		Xi   = X[i];
		Yi   = Y[i];
		mx  += Xi; 
		sxx += Xi*Xi;
		my  += Yi; 
		syy += Yi*Yi;
		sxy += Xi*Yi;
	}

	mx  /= n;
	my  /= n;
	sxx /= n;
	syy /= n;
	sxy /= n;
	
	difxy = sxx+syy-2*sxy;
	sumxy = sxx+syy+2*sxy;

	sxx -= mx*mx;
	syy -= my*my;

	sxy -= mx*my;

	*corr = sxy * isqrt_(sxx) * isqrt_(syy);
	*rsr  = difxy/sumxy;

	return 0;
}


int pDenMapCorrError(PDen_t * map1, PDen_t * map2, real *corr, real *e_corr, size_t nboot)
{
	size_t i,j,k,n;
	real cm, cs;
	real *X,*Y;
	register real  Xi,Yi;
	real mx,my,sxx,syy,sxy;
	cm = cs = 0.0;

	X = map1->data;
	Y = map2->data;
	n = map1->n;

	for ( k = 0; k < nboot; k++) {
		mx = my = sxx = syy = sxy = 0.;

		for( i = 0; i < n; i++){
			j = rand() % n;

			Xi   = X[j];
			Yi   = Y[j];
			mx  += Xi; 
			sxx += Xi*Xi;
			my  += Yi; 
			syy += Yi*Yi;
			sxy += Xi*Yi;
		}
		mx  /= n;
		my  /= n;
		sxx /= n;
		syy /= n;
		sxy /= n;

		sxx -= mx*mx;
		syy -= my*my;

		sxy -= mx*my;

		sxy = sxy * isqrt_(sxx) * isqrt_(syy);
		
		cm += sxy;
		cs += sxy*sxy;
	}
	cm /= nboot;
	cs /= nboot;

	cs -= cm*cm;

	*corr   = cm;
	*e_corr = cs * isqrt_(cs);

	return 0;
}

int pDenNormalizeTo( PDen_t * this, real m, real s )
{
	size_t i;
	for(i=0;i<this->n;i++)
		this->data[i] = ( this->data[i] - m ) / s;
	return 0;
}

int pDenNormalize( PDen_t * this )
{
	real m,s;
	pDenMeanSd( this, &m, &s);	
	debug("Mean=%lf sd=%lf",m,s);
	pDenNormalizeTo( this, m, s);
	debug("Value: %lf\n",this->data[4]);
	return 0;
}

int pDenScaleBy(PDen_t * this, real value)
{
	size_t i;
	for(i=0;i<this->n;i++)
		this->data[i]*=value;
	return 0;
}

int PDenScaleMap(PDen_t * this, PDen_t * map)
{
	size_t i;
	for(i=0;i<this->n;i++)
		this->data[i]*=map->data[i];
	return 0;
}

int PDenMask(PDen_t * this, PDen_t * mask)
{

	return 0;
}

