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

#include "pden_sfrefine.h"


//INTERNAL FUNCTIONS

//Babinet Functions
static real SFFunc_Babinet(const real s2, const real calc, const real mask, const real x0[])
{
	return x0[0]*exp(-x0[1] * s2/4) * (calc - mask * x0[2] * exp(-x0[3] * s2/4));
}

static real SFFunc_Babinet_dB_glob(const real s2, const real calc, const real mask, const real x0[])
{
	return (-s2/4) * x0[0]*exp(-x0[1] * s2/4) * (calc - mask * x0[2] * exp(-x0[3] * s2/4));
}

static real SFFunc_Babinet_dk_glob(const real s2, const real calc, const real mask, const real x0[])
{
	return exp(-x0[1] * s2/4) * (calc - mask * x0[2] * exp(-x0[3] * s2/4));
}

static real SFFunc_Babinet_dB_solv(const real s2, const real calc, const real mask, const real x0[])
{
	return (-s2/4) * x0[0]*exp(-x0[1] * s2/4) * -1. * mask * x0[0] * exp(-x0[3] * s2/4);
}

static real SFFunc_Babinet_dk_solv(const real s2, const real calc, const real mask, const real x0[])
{
	return x0[0]*exp(-x0[1] * s2/4) * -1. * mask * exp(-x0[3] * s2/4);
}



static real SFValueFunc_Babinet(const PDen_t *X[], const real param0[])
{

	real val=0;
	debug("get value Babinet:");

	#ifdef __BABINET_C_OPENMP 
	#pragma omp parallel reduction(+:val) 
	#endif
	{
	size_t  k,j,i,p;
	ssize_t u,v;
	size_t  ni,nj,nk;
	real Fm , Fc, Fr, s;
	real modx,mody,modz;
	real h;
	const cplx * res, * cal, *mas ;
	// calculate data dimensions
  	nk = X[0]->size.z;
  	nj = X[0]->size.y;
  	ni = ( X[0]->size.x / 2 + 1 );

	// fourier scaling factor
	modx = X[0]->apix.x*X[0]->size.x;
	mody = X[0]->apix.y*X[0]->size.y;
	modz = X[0]->apix.z*X[0]->size.z;

	modx = 1./(modx*modx);
	mody = 1./(mody*mody);
	modz = 1./(modz*modz);

	//scaling factor
	res = (const cplx *)X[0]->data;
	cal = (const cplx *)X[1]->data;
	mas = (const cplx *)X[2]->data;

	//debug("Pointer %lf %p %p",X[0]->data[56],X[1],X[2]);
	
	#ifdef __BABINET_C_OPENMP
	#pragma omp for schedule(static)
	#endif
	//element loop
  	for(k=0;k<nk;k++){
  	  	for(j=0;j<nj;j++){
      			for(i=0;i<ni;i++){
				// array position
				p = (k * nj + j) * ni + i;

        			//radius
        			u = dft_c(k,nk);
        			v = dft_c(j,nj);

        			s = (real)(i*i)*(modx)
        			   +(real)(v*v)*(mody)
         			   +(real)(u*u)*(modz);

				// Structure Factors
				Fc = abscplx(cal[p]);
				Fm = abscplx(mas[p]);
				Fr = abscplx(res[p]);			    
					    
				h = (SFFunc_Babinet(s, Fc, Fm,param0)-Fr);
				val += h*h;

				//debug("%lu: s=%lf Fc=%lf Fm=%lf Fr=%lf h=%lf \n",p,s,Fc,Fm,Fr,h);
			}
		}
	}//end element loop	
	} //OMP PARALLEL END
	debug("value = %lf(=%lf/%lu)",val/X[0]->n,val,X[0]->n);
	return val/X[0]->n;
}

static void SFGradFunc_Babinet(const PDen_t *X[], const real param0[], real grad[])
{
	real S;

	debug("get Gradient:");
	grad[0]=0.;
	grad[1]=0.;
	grad[2]=0.;
	grad[3]=0.;
	

	#ifdef __BABINET_C_OPENMP 
	#pragma omp parallel 
	#endif
	{
	size_t  k,j,i,p;
	ssize_t u,v;
	size_t  ni,nj,nk;
	real Fm , Fc, Fr, s, h;
	real modx,mody,modz;
	const cplx * res, * cal, *mas ;

	real grad_[4]={0.,0.,0.,0.};

	// calculate data dimensions
	nk = X[0]->size.z;
	nj = X[0]->size.y;
	ni = ( X[0]->size.x / 2 + 1 );
	
	// fourier scaling factor
	modx = X[0]->apix.x*X[0]->size.x;
	mody = X[0]->apix.y*X[0]->size.y;
	modz = X[0]->apix.z*X[0]->size.z;
	
	modx = 1./(modx*modx);
	mody = 1./(mody*mody);
	modz = 1./(modz*modz);
	
	//scaling factor
	res = (const cplx *) X[0]->data;
	cal = (const cplx *) X[1]->data;
	mas = (const cplx *) X[2]->data;
	#ifdef __BABINET_C_OPENMP 
	#pragma omp for schedule(static)  
	#endif
	//element loop
	for(k=0;k<nk;k++){
	  	for(j=0;j<nj;j++){
			for(i=0;i<ni;i++){
				// array position
				p = (k * nj + j) * ni + i;
	
				//radius
				u = dft_c(k,nk);
				v = dft_c(j,nj);
	
				s = (real)(i*i)*(modx)
				 	+(real)(v*v)*(mody)
	 				+(real)(u*u)*(modz);
	
				// Structure Factors
				Fc = abscplx(cal[p]);
				Fm = abscplx(mas[p]);
				Fr = abscplx(res[p]);			    
				debug("Fr=%lf  Fc=%lf  Fm=%lf",Fr,Fc,Fm);
				h = 2*(SFFunc_Babinet(s, Fc, Fm, param0)-Fr);
				grad_[0] += h * SFFunc_Babinet_dk_glob(s, Fc, Fm, param0);
				grad_[1] += h * SFFunc_Babinet_dB_glob(s, Fc, Fm, param0);
				grad_[2] += h * SFFunc_Babinet_dk_solv(s, Fc, Fm, param0);
				grad_[3] += h * SFFunc_Babinet_dB_solv(s, Fc, Fm, param0);
			}
		}
	}//end element loop	

	#ifdef __BABINET_C_OPENMP 
	#pragma omp atomic
	#endif
	grad[0]+=grad_[0];
	#ifdef __BABINET_C_OPENMP 
	#pragma omp atomic
	#endif
	grad[1]+=grad_[1];
	#ifdef __BABINET_C_OPENMP 
	#pragma omp atomic
	#endif
	grad[2]+=grad_[2];
	#ifdef __BABINET_C_OPENMP 
	#pragma omp atomic
	#endif
	grad[3]+=grad_[3];
	
	} //OMP PARALLEL END


	S = -1./(real) X[0]->n;
	grad[0] *=S;
	grad[1] *=S;
	grad[2] *=S;
	grad[3] *=S;
	debug("grad(%9.4lf %9.4lf %9.4lf %9.4lf)",grad[0],grad[1],grad[2],grad[3]);
}

//library functions
PDen_t *  pDenApplySF_Babinet (PDen_t * result, PDen_t *input, PDen_t *mask, const real k_glob, const real B_glob, const real k_solv, const real B_solv)
{

	pDenFFTNormal(input,1);
	pDenFFTNormal(mask,1);

	#ifdef __BABINET_C_OPENMP 
	#pragma omp parallel 
	#endif
	{
	size_t  k,j,i,p;
	ssize_t u,v;
	size_t  ni,nj,nk;
	real Fm , Fc, s;
	real modx,mody,modz;
	const cplx  * inp , * mas;
	 cplx * tar;
	real x0[4];

	x0[0] = k_glob;
	x0[1] = B_glob;
	x0[2] = k_solv;
	x0[3] = B_solv;

	// calculate data dimensions
  	nk = input->size.z;
  	nj = input->size.y;
  	ni = ( input->size.x / 2 + 1 );

	// fourier scaling factor
	modx = input->apix.x*input->size.x;
	mody = input->apix.y*input->size.y;
	modz = input->apix.z*input->size.z;

	modx = 1./(modx*modx);
	mody = 1./(mody*mody);
	modz = 1./(modz*modz);

	//scaling factor
	tar = (cplx *)result->data;
	inp = (const cplx *)input->data;
	mas = (const cplx *)mask->data;

	#ifdef __BABINET_C_OPENMP 
	#pragma omp for schedule(static)  
	#endif
  	for(k=0;k<nk;k++){
  	  	for(j=0;j<nj;j++){
      			for(i=0;i<ni;i++){
				// array position
				p = (k * nj + j) * ni + i;

        			//radius
        			u = dft_c(k,nk);
        			v = dft_c(j,nj);

        			s = (real)(i*i)*(modx)
        			   +(real)(v*v)*(mody)
         			   +(real)(u*u)*(modz);

				// Structure Factors
				Fc = abscplx(inp[p]);
				Fm = abscplx(mas[p]);
				Fm = SFFunc_Babinet(s,Fc,Fm,x0);

				s = Fm/Fc;

				tar[p].r = s*inp[p].r;
				tar[p].i = s*inp[p].i;
			}
		}
	}
	} //OMP PARALLEL END
	result->mode |= PDEN_MODE_PHASE_SPACE;
	return result;
}

int      pDenRefineSF_Babinet (PDen_t * model, PDen_t *calc, PDen_t *mask, real *k_glob, real * B_glob, real * k_solv, real * B_solv, const real prec, const size_t steps)
{
	int r;
	real x0[4] = {
		 1.0,
		 3.0,
		 0.35,
		46.0
	};
	PDen_t *maps[] = {
		model,
		calc,
		mask
	};

	pDenFFTNormal(model,1);
	pDenFFTNormal(calc,1);
	pDenFFTNormal(mask,1);

	if(*k_glob != 0. ) x0[0] = *k_glob;
	if(*B_glob != 0. ) x0[1] = *B_glob;
	if(*k_solv != 0. ) x0[2] = *k_solv;
	if(*B_solv != 0. ) x0[3] = *B_solv;

	info("Structure factor refinment:");
	

	say("Step      k_global  B_global  k_solvent B_solvent Value\n");
	r = CGRefine((const PDen_t **)maps,x0,SFValueFunc_Babinet, SFGradFunc_Babinet,4, prec, steps);

	if(r)
		error("Babinet stucture factor refinment failed with parameters: k_global=%lf  B_global=%lf  k_solvent=%lf B_solvent=%lf", x0[0],x0[1],x0[2],x0[3])
	else 
		info("Babinet: F_m=%.2lf*exp(-%.1lf*s^2/4)*(F_c-F_mask*%.2lf*exp(-%.1lf*s^2/4))", x0[0],x0[1],x0[2],x0[3]);

	*k_glob = x0[0];
	*B_glob = x0[1];
	*k_solv = x0[2];
	*B_solv = x0[3];


	return r;
}
