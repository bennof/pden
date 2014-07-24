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


int pDenCalcSF(PDen_t * input,real *data, size_t n)
{
	size_t  k,j,i, idx, *nc; 
	ssize_t u,v;
	size_t  ni,nj,nk;
	real    h,scale;
	real 	modx,mody,modz;
	cplx * in;
	//real scale2;

	debug("Calculate power spectrum");
	//calculate data dimensions
  	nk = input->size.z;
  	nj = input->size.y;
  	ni = ( input->size.x / 2 + 1 );
	
	pDenFFTNormal(input,1);
	//pDenFFT(input,1);

	//setup arrays
	memset(data,0,n*sizeof(real));
	nc = (size_t*) calloc(n,sizeof(size_t));

	// fourier scaling factor
	modx = input->apix.x*input->size.x;
	mody = input->apix.y*input->size.y;
	modz = input->apix.z*input->size.z;

	modx = 1./(modx*modx);
	mody = 1./(mody*mody);
	modz = 1./(modz*modz);


	//histo width scaling
	scale = 0.25 / ( input->apix.x*input->apix.x);
	scale  = 1 / sqrt(scale); 
	scale *= (real)n;

	//scale2 = 1./sqrt((real)input->n);
	in = (cplx*) input->data;
  	for(k=0;k<nk;k++){
  	  	for(j=0;j<nj;j++){
      			for(i=0;i<ni;i++){
        			//radius
        			u = dft_c(k,nk);
        			v = dft_c(j,nj);

        			h = (real)(i*i)*(modx)
        			  + (real)(v*v)*(mody)
         			  + (real)(u*u)*(modz);
				h  = sqrt(h);

				h *= scale;

				idx = (size_t) (h); // calc index

				if(idx<n){
					nc[idx]++; //count entries 
					h = abscplx(*in);
					data[idx] += h;
					debug("|F|=%lf|",h);
				}
				in++;
			}
		}
	}
	//average
	for(i=0;i<n;i++){
		data[i]	/= nc[i];
	}
	free(nc);
	return 0;
}


