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

#include "split2fft.h"  

int pDenFFT(PDen_t * this,int direction)
{
	debug("FFT (%s)",(direction>0)?"forward":"backward");
	if(!this->fft) 
		this->fft=split2FFTNew(this,1);
	return split2FFTExecute(this->fft,this,direction);
}


int pDenFFTNormal(PDen_t * this,int direction)
{
	debug("Normalized FFT (%s)",(direction>0)?"forward":"backward");
	if(!this->fft) 
		this->fft=split2FFTNew(this,1);
	return split2FFTExecuteN(this->fft,this,direction);
}

int pDenSetupFFT(PDen_t * this,int optimize)
{
	this->fft=split2FFTNew(this,optimize);
	return 0;
}

int pDenGetSizeSplitCV     (PDen_t * input, real lb1,real ub1,real lb2,real ub2, size_t *nused, size_t *nfree)
{
	size_t  k,j,i; 
	size_t  ni,nj,nk;
	ssize_t u,v;
	size_t  nu,nf;
	real    h;
	real 	modx,mody,modz;

	nu=nf=0;

	debug("split density");

	modx = input->apix.x*input->size.x;
	mody = input->apix.y*input->size.y;
	modz = input->apix.z*input->size.z;

	modx = 1./(modx*modx);
	mody = 1./(mody*mody);
	modz = 1./(modz*modz);

	//calculating inverse square values
 	lb1=1./(lb1*lb1);
	lb2=1./(lb2*lb2);
  	ub1=1./(ub1*ub1);
  	ub2=1./(ub2*ub2);

	//calculate data dimensions
  	nk = input->size.z;
  	nj = input->size.y;
  	ni = ( input->size.x / 2 + 1 );

  	for(k=0;k<nk;k++){
  	  	for(j=0;j<nj;j++){
      			for(i=0;i<ni;i++){
        			//radius
        			u = dft_c(k,nk);
        			v = dft_c(j,nj);

        			h = (double)(i*i)*(modx)
        			 	+(double)(v*v)*(mody)
         				+(double)(u*u)*(modz);

				if ((h <= ub2) && (h >= lb2)){ //if free set
					nf++;
				}
				else if((h <ub1)  && (h >= lb1)){ //if in use set 
					nu++;
				}
			}
		}
	}
	*nused = nu;
	*nfree = nf;
	return 0;
}

PDen_t * pDenSplitCV       (PDen_t * result, PDen_t * free, PDen_t * input,real lb1,real ub1,real lb2,real ub2)
{
	size_t  k,j,i; 
	size_t  ni,nj,nk;
	ssize_t u,v;
	real    *d,*f,*r,h;
	real 	modx,mody,modz;

	pDenFFTNormal(input,1);

	modx = input->apix.x*input->size.x;
	mody = input->apix.y*input->size.y;
	modz = input->apix.z*input->size.z;

	modx = 1./(modx*modx);
	mody = 1./(mody*mody);
	modz = 1./(modz*modz);

	//calculating inverse square values
 	lb1=1./(lb1*lb1);
	lb2=1./(lb2*lb2);
  	ub1=1./(ub1*ub1);
  	ub2=1./(ub2*ub2);

	debug("scaling: %lf %lf %lf",modx,mody,modz);
	debug("Bounds: %lf %lf %lf %lf",lb1,ub1,lb2,ub2);


	//calculate data dimensions
  	nk = input->size.z;
  	nj = input->size.y;
  	ni = ( input->size.x / 2 + 1 );

	f = free->data;
	d = input->data;
	r = result->data;

  	for(k=0;k<nk;k++){
  	  	for(j=0;j<nj;j++){
      			for(i=0;i<ni;i++){
        			//radius
        			u = dft_c(k,nk);
        			v = dft_c(j,nj);

        			h = (double)(i*i)*(modx)
        			 	+(double)(v*v)*(mody)
         				+(double)(u*u)*(modz);

				if ((h <= ub2) && (h >= lb2)){ //if free set
					(*f++) = (*d++);
					(*f++) = (*d++);
					(*r++) = 0.;
					(*r++) = 0.;
				}
				else if((h <ub1)  && (h >= lb1)){ //if in use set 
					(*f++) = 0.;
					(*f++) = 0.;
					(*r++) = (*d++);
					(*r++) = (*d++);
				}
				else{  // not used (0.0)
					(*f++) = 0.;
					(*f++) = 0.;
					(*r++) = 0.;
					(*r++) = 0.;
					d+=2; 
				}
			}
		}
	}
	free->mode   |= PDEN_MODE_PHASE_SPACE;
	result->mode |= PDEN_MODE_PHASE_SPACE;
	return result;
}

real pDenRFactor(PDen_t * a, PDen_t * b)
{
	size_t  k;
	size_t  ni,nj,nk,n;
	real pow1,pow2,r,c,scale, R, d;
	real *m1,*m2;

	debug("compute R factor");

	// test for phase space
	pDenFFTNormal(a,1);
	pDenFFTNormal(b,1);

	//calculate data dimensions
  	nk = a->size.z;
  	nj = a->size.y;
  	ni = ( a->size.x / 2 + 1 );
	n=nk*nj*ni;

	//scaling factor
	m1 = a->data;
	m2 = b->data;
	pow1=0.;
	pow2=0.;
  	for(k=0;k<n;k++){
		r = (*m1++);
		c = (*m1++);
		r = ( r*r / c*c );
		if ( r > 0. ){
			r = sqrt ( r );
			pow1 += r;
		}
		r = (*m2++);
		c = (*m2++);
		r = ( r*r / c*c );
		if ( r > 0. ){
			r = sqrt ( r );
			pow2 += r;
		}
	}
	scale = pow1 / pow2 ;

	// compute R Value
	R = d = 0.0;
	m1 = a->data;
	m2 = b->data;
  	for(k=0;k<nk;k++){
		r = (*m1++);
		c = (*m1++);
		r = ( r*r / c*c );
		pow1 = (r>0.)? sqrt ( r ) : 0.;
		r = (*m2++);
		c = (*m2++);
		r = ( r*r / c*c );
		pow2 = (r>0.)? sqrt ( r ) : 0.;
		R += fabs ( pow1 - pow2 * scale );
		d += pow1;		
	}
	

	return R/d;
}


#define binop( sign ) \
{\
	size_t i,n;\
	pDenFFT(a,1);\
	pDenFFT(b,1);\
	n = result->size.z*result->size.y*(result->size.x+2);\
	for(i=0;i<n;i++)\
		result->data[i] = a->data[i] sign b->data[i];\
	result->mode |=PDEN_MODE_PHASE_SPACE;\
	return result;\
}


PDen_t * pDenAddFS (PDen_t * result, PDen_t * a, PDen_t * b )
binop( + )

PDen_t * pDenSubFS (PDen_t * result, PDen_t * a, PDen_t * b )
binop( - )

