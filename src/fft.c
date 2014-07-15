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
		this->fft=split2FFTNew(this,0);
	return split2FFTExecuteN(this->fft,this,direction);
}

int pDenSetupFFT(PDen_t * this,int optimize)
{	
	if (this->fft)
		split2FFTDelete(this->fft);	
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
	//calculate data dimensions
  	const size_t nk = input->size.z;
  	const size_t nj = input->size.y;
  	const size_t ni = ( input->size.x / 2 + 1 );
	real 	modx,mody,modz;

	modx = input->apix.x*input->size.x;
	mody = input->apix.y*input->size.y;
	modz = input->apix.z*input->size.z;

	modx = 1./(modx*modx);
	mody = 1./(mody*mody);
	modz = 1./(modz*modz);


	pDenFFTNormal(input,1);

	//calculating inverse square values
 	lb1=1./(lb1*lb1);
	lb2=1./(lb2*lb2);
  	ub1=1./(ub1*ub1);
  	ub2=1./(ub2*ub2);

	debug("Bounds: %lf %lf %lf %lf",lb1,ub1,lb2,ub2);

	#ifdef __FFT_C_OPENMP 
	#pragma omp parallel 
	#endif
	{
	size_t  k,j,i,p; 
	ssize_t u,v;
	const cplx    *d;

	cplx    *f,*r;
	real    h;

	f = (cplx*)free->data;
	d = (const cplx*)input->data;
	r = (cplx*)result->data;

	#ifdef __FFT_C_OPENMP
	#pragma omp for schedule(static)
	#endif
  	for(k=0;k<nk;k++){
        	u = dft_c(k,nk);
  	  	for(j=0;j<nj;j++){
        		v = dft_c(j,nj);
      			for(i=0;i<ni;i++){
        			//radius
        			h = (double)(i*i)*(modx)
        			 	+(double)(v*v)*(mody)
         				+(double)(u*u)*(modz);


				// array position
				p = (k * nj + j) * ni + i;

				if ((h <= ub2) && (h >= lb2)){ //if free set
					f[p].r = d[p].r;
					f[p].i = d[p].i;
					r[p].r = 0.;
					r[p].i = 0.;
				}
				else if((h <ub1)  && (h >= lb1)){ //if in use set 
					r[p].r = d[p].r;
					r[p].i = d[p].i;
					f[p].r = 0.;
					f[p].i = 0.;
				}
				else{  // not used (0.0)
					r[p].r = 0.;
					r[p].i = 0.;
					f[p].r = 0.;
					f[p].i = 0.;
				}
			}
		}
	}
	} //OMP PARALLEL END	
	free->mode   |= PDEN_MODE_PHASE_SPACE;
	result->mode |= PDEN_MODE_PHASE_SPACE;
	return result;
}


real pDenRFactor(PDen_t * a, PDen_t * b)
{
	size_t  k;
	size_t  ni,nj,nk,n;
	real pow1,pow2,r,c,scale, R, d;
	const cplx *m1,*m2;

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
	m1 = (const cplx*)a->data;
	m2 = (const cplx*)b->data;
	pow1=0.;
	pow2=0.;

	R=0.;
	d=0.;


	#ifdef __FFT_C_OPENMP 
	#pragma omp parallel private(k,r,c) 
	#endif
	{
	real pow1_=0., pow2_=0.;
	real R_=0., d_=0.;
	#ifdef __FFT_C_OPENMP
	#pragma omp for schedule(static)
	#endif
  	for(k=0;k<n;k++){
		r = m1[k].r;
		c = m1[k].i;
		r = ( r*r / c*c );
		if ( r > 0. ){
			r = sqrt ( r );
			pow1_ += r;
		}
		r = m2[k].r;
		c = m2[k].i;
		r = ( r*r / c*c );
		if ( r > 0. ){
			r = sqrt ( r );
			pow2_ += r;
		}
	}
	#pragma omp atomic
	pow1 += pow1_;
	#pragma omp atomic
	pow2 += pow2_;

	#pragma omp barrier
	scale = pow1 / pow2 ;

	// compute R Value
	#ifdef __FFT_C_OPENMP
	#pragma omp for schedule(static)
	#endif
  	for(k=0;k<nk;k++){
		r = m1[k].r;
		c = m1[k].i;
		r = ( r*r / c*c );
		pow1_ = (r>0.)? sqrt ( r ) : 0.;
		r = m2[k].r;
		c = m2[k].i;
		r = ( r*r / c*c );
		pow2_ = (r>0.)? sqrt ( r ) : 0.;
		R_ += fabs ( pow1_ - pow2_ * scale );
		d_ += pow1_;		
	}

	#pragma omp atomic 
	R += R_;
	#pragma omp atomic 
	d += d_;

	} //OMP PARALLEL END
	

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



PDen_t * pDenShiftInGrid(PDen_t * result, PDen_t * a,real dx,real dy, real dz)
{
        real *d , *r; 
        size_t i,j,k,xx,yy,zz,p;
        double re[6],im[6],h;

	dx /= a->apix.x;
	dy /= a->apix.y;
	dz /= a->apix.z;

        re[0] = cos(2.0f*PI/a->size.x*dx);
        im[0] = -1.0f*sin(2.0f*PI/a->size.x*dx);

        re[1] = cos(2.0f*PI/a->size.y*dy);
        im[1] = -1.0f*sin(2.0f*PI/a->size.y*dy);

        re[2] = cos(-1.0f*PI*dy);
        im[2] = -1.0f*sin(-1.0f*PI*dy);

        re[5] = cos(2.0f*PI/a->size.z*dz);
        im[5] = -1.0f*sin(2.0f*PI/a->size.z*dz);

        re[4] = cos(-1.0f*PI*dz);
        im[4] = -1.0f*sin(-1.0f*PI*dz);

        re[4] = cos(-1.0f*PI*dz);
        im[4] = -1.0f*sin(-1.0f*PI*dz);
        xx=a->size.x/2+1;

	// test for phase space
	pDenFFTNormal(a,1);

        d = a->data;
	r = result->data;
        for(k=0;k<a->size.z;k++){
        	zz=(k+a->size.z/2)%a->size.z;
        	re[2]=re[4];
        	im[2]=im[4];

        	for(j=0;j<a->size.y;j++){
        		yy=(j+a->size.y/2)%a->size.y;
        		re[3]=re[2];
        		im[3]=im[2];
        		for(i=0;i<xx;i++){
        		        p=2*(i+yy*xx+zz*a->size.x*xx);

        		        h=d[p];
        		        r[p  ] = h*re[3]-d[p+1]*im[3];
        		        r[p+1] = d[p+1]*re[3]+h*im[3];

        		        h=re[3];
        		        re[3] = h*re[0]-im[3]*im[0];
        		        im[3] = h*im[0]+im[3]*re[0];
        		}   
        		h=re[2];
        		re[2] = h*re[1]-im[2]*im[1];
        		im[2] = h*im[1]+im[2]*re[1];
        	}   
        	h=re[4];
        	re[4] = h*re[5]-im[4]*im[5];
        	im[4] = h*im[5]+im[4]*re[5];
        }   
	result->mode |=PDEN_MODE_PHASE_SPACE;
	return result;
}

PDen_t * pDenGaussFilter(PDen_t * result, PDen_t * a,real sigma)
{

        size_t i,j,k,p;
	size_t  ni,nj,nk;
        ssize_t u,v;
        const cplx *m;
	cplx *r;
	real h;
        const real k0 = 1.0/(a->apix.x*a->apix.x*a->size.x*a->size.x);
        const real k1 = 1.0/(a->apix.y*a->apix.y*a->size.y*a->size.y);
        const real k2 = 1.0/(a->apix.z*a->apix.z*a->size.z*a->size.z);
        const real C = sqrt(2*PI)*sigma;
        const real c = -2.0*PI*PI*sigma*sigma;

	// calculate data dimensions
	nk = a->size.z;
	nj = a->size.y;
	ni = ( a->size.x / 2 + 1 );

	// test for phase space
	pDenFFTNormal(a,1);

        m = (const cplx*)a->data;
        r = (cplx*)result->data;
        for(k=0;k<nk;k++){
                v=dft_c(k,nk);
                for(j=0;j<nj;j++){
                        u=dft_c(j,nj);
                        for(i=0;i<ni;i++){
                                h=(real)(i*i)*k0+(real)(u*u)*k1+(real)(v*v)*k2;
				// array position
				p = (k * nj + j) * ni + i;

				h = C*exp(c*h);

                                r[p].r = m[p].r*h;
                                r[p].i = m[p].i*h;
                        }
                }
        }
	result->mode |=PDEN_MODE_PHASE_SPACE;
	return result;
}

PDen_t * pDenLaplaceFilter(PDen_t * result, PDen_t * a,real sigma)
{
        size_t i,j,k,p;
	size_t  ni,nj,nk;
        ssize_t u,v;
        real h; 
	const cplx *m;
	cplx *r;
        const real k0 = 1.0/(a->apix.x*a->apix.x*a->size.x*a->size.x);
        const real k1 = 1.0/(a->apix.y*a->apix.y*a->size.y*a->size.y);
        const real k2 = 1.0/(a->apix.z*a->apix.z*a->size.z*a->size.z);
	// calculate data dimensions
	nk = a->size.z;
	nj = a->size.y;
	ni = ( a->size.x / 2 + 1 );

	// test for phase space
	pDenFFTNormal(a,1);

        m = (const cplx*)a->data;
        r = (cplx*)result->data;

	for(k=0;k<nk;k++){
        	v=dft_c(k,nk);
        	for(j=0;j<nj;j++){
        		u=dft_c(j,nj);
        		for(i=0;i<ni;i++){
        			h=((real)(i*i)*k0+(real)(u*u)*k1+(real)(v*v)*k2);
				// array position
				p = (k * nj + j) * ni + i;

                                r[p].r = m[p].r*h;
                                r[p].i = m[p].i*h;
        		}
        	}
        }

	result->mode |=PDEN_MODE_PHASE_SPACE;
	return result;
}

PDen_t * pDenRampFilter(PDen_t * result, PDen_t * a,real sigma)
{
        size_t i,j,k;
        ssize_t u,v;
        real h, *m, *r;
        const real k0 = 1.0/(a->apix.x*a->apix.x*a->size.x*a->size.x);
        const real k1 = 1.0/(a->apix.y*a->apix.y*a->size.y*a->size.y);
        const real k2 = 1.0/(a->apix.z*a->apix.z*a->size.z*a->size.z);

	// test for phase space
	pDenFFTNormal(a,1);

        m = a->data;
        r = result->data;

	for(k=0;k<a->size.z;k++){
        	v=dft_c(k,a->size.z);
        	for(j=0;j<a->size.y;j++){
        	        u=dft_c(j,a->size.y);
        	        for(i=0;i<(a->size.x/2+1);i++){
        	                h=(real)(i*i)*k0+(real)(u*u)*k1+(real)(v*v)*k2;
        	                h=(h==0.0)?0.0:h*isqrt_(h);
        	                (*r++) *= h*(*m++);
        	                (*r++) *= h*(*m++);
        	        }
        	}
        }

	result->mode |=PDEN_MODE_PHASE_SPACE;
	return result;
}
