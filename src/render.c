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


static real pDenGaussKernel(real r,int type, real sigma)
{
	return exp( -.5 * r / sigma );
}


real PengGFactor = 10.0;
real PengScatterFactor_A[][5] = {
{ 0.0489, 0.2091, 0.7537, 1.1420, 0.3555 }, // ANY
{ 0.0088, 0.0449, 0.1481, 0.2356, 0.0914 }, // H
{ 0.0489, 0.2091, 0.7537, 1.1420, 0.3555 }, // C
{ 0.0267, 0.1328, 0.5301, 1.1020, 0.4215 }, // N
{ 0.0365, 0.1729, 0.5805, 0.8814, 0.3121 }, // O
{ 0.0915, 0.4312, 1.0847, 2.4671, 1.0852 }, // S
{ 0.1005, 0.4615, 1.0663, 2.5854, 1.2725 }, // P
{ 0.5229, 2.2874, 4.7243, 5.0807, 5.6389 }, // BA
{ 0.1260, 0.6442, 0.8893, 1.8197, 1.2988 }, // NA
{ 0.0799, 0.3891, 1.0037, 2.3332, 1.0507 }, // CL
{ 0.2149, 0.8703, 2.4999, 2.3591, 3.0318 }, // K
{ 0.1130, 0.5575, 0.9046, 2.1580, 1.4735 }  // MG
};

real PengScatterFactor_B[][5] = {
{ 0.1140, 1.0825, 5.4281, 17.8811,  51.1341 }, // ANY
{ 0.1152, 1.0867, 4.9755, 16.5591,  43.2743 }, // H
{ 0.1140, 1.0825, 5.4281, 17.8811,  51.1341 }, // C
{ 0.0541, 0.5165, 2.8207, 10.6297,  34.3764 }, // N
{ 0.0652, 0.6184, 2.9449, 9.62980,  28.2194 }, // O
{ 0.0838, 0.7788, 4.3462, 15.5846,  44.6365 }, // S
{ 0.0977, 0.9084, 4.9654, 18.5471,  54.3648 }, // P
{ 0.1434, 1.6019, 9.4511, 42.7685, 148.4969 }, // BA
{ 0.1684, 1.7150, 8.8386, 50.8265, 147.2073 }, // NA
{ 0.0694, 0.6443, 3.5351, 12.5058,  35.8633 }, // CL
{ 0.1660, 1.6906, 8.7447, 46.7825, 165.6923 }, // K
{ 0.1356, 1.3579, 6.9255, 32.3165,  92.1138 }  // MG
};

static real pDenPengKernel(real r,int type, real sigma)
{	
	/* 
	 * Fourier transform of the scattering factor
	 * Kirkland, "Advanced Computing in Electron Microscopy"
	 * chapter C, Eq. C19.
	 */

	register real v, h1, h2, h3, h4, h0 ;
	real sf[5], sf_a[5];

	h0 = -9.8696 * r;
	
	//prelaoding
	sf[0]   = PengScatterFactor_B[type][0];
	sf[1]   = PengScatterFactor_B[type][1];
	sf[2]   = PengScatterFactor_B[type][2];
	sf[3]   = PengScatterFactor_B[type][3];
	sf[4]   = PengScatterFactor_B[type][3];
	sf_a[0] = PengScatterFactor_A[type][0];
	sf_a[1] = PengScatterFactor_A[type][1];
	sf_a[2] = PengScatterFactor_A[type][2];
	sf_a[3] = PengScatterFactor_A[type][3];
	sf_a[4] = PengScatterFactor_A[type][4];

	//unrolled loop
	// pass1
	h1 = (sf[0]);
	h2 = 1.0 / (h1 + PengGFactor);
	h3 = 1.0 / h1;
	h4 = sqrt( h3 * h3 * h3 );
	v  = sf_a[0] * h4 * exp ( h0 * h2 );

	// pass2
	h1 = (sf[1]);
	h2 = 1.0 / (h1 + PengGFactor);
	h3 = 1.0 / h1;
	h4 = sqrt( h2 * h2 * h2 );
	v += sf_a[1] * h4 * exp ( h0 * h2 );

	// pass3
	h1 = (sf[2]);
	h2 = 1.0 / (h1 + PengGFactor);
	h3 = 1.0 / h1;
	h4 = sqrt( h2 * h2 * h2 );
	v += sf_a[2] * h4 * exp ( h0 * h2 );

	// pass4
	h1 = (sf[3]);
	h2 = 1.0 / (h1 + PengGFactor);
	h3 = 1.0 / h1;
	h4 = sqrt( h2 * h2 * h2 );
	v += sf_a[3] * h4 * exp ( h0 * h2 );
	
	// pass5
	h1 = (sf[4]);
	h2 = 1.0 / (h1 + PengGFactor);
	h3 = 1.0 / h1;
	h4 = sqrt( h2 * h2 * h2 );
	v += sf_a[4] * h4 * exp ( h0 * h2 );

	return v;
}




#ifdef PDEN_RENDER_NM_COORDS 
#define getAtomPosition \
pos.x = coords[3 * atom    ] * 10. - result->origin.x ;\
pos.y = coords[3 * atom + 1] * 10. - result->origin.y ;\
pos.z = coords[3 * atom + 2] * 10. - result->origin.z ;
#else 
#define getAtomPosition \
pos.x = coords[3 * atom    ] - result->origin.x ;\
pos.y = coords[3 * atom + 1] - result->origin.y ;\
pos.z = coords[3 * atom + 2] - result->origin.z ;
#endif 



#define render(kernel) \
	size_t atom;\
	size_t idx, idxx;\
	size_t i,j,k;\
	struct {size_t x,y,z;} gridid;\
	struct {size_t x,y,z;} gridend;\
	struct {size_t x,y,z;} walk;\
	struct {real   x,y,z;} pos;\
	struct {real   x,y,z;} r;\
	real lenr;\
	real density;\
	real mass;\
\
	size_t dim2,dim3;\
	int type_;\
\
	sigma *= sigma;\
\
	/*walk radius*/\
	walk.x = ((size_t)(0.5 * width/result->apix.x))+1;\
	walk.y = ((size_t)(0.5 * width/result->apix.y))+1; \
	walk.z = ((size_t)(0.5 * width/result->apix.z))+1;\
\
\
	dim2 = result->size.x;\
	dim3 = result->size.x*result->size.y;\
	/*render loop*/\
	for(atom=0;atom<n;atom++){\
		mass = 1.0;\
		\
		getAtomPosition\
\
		/*get upper left grid point*/\
		gridid.x = (size_t)( ( pos.x ) / result->apix.x );	\
		gridid.y = (size_t)( ( pos.y ) / result->apix.y );	\
		gridid.z = (size_t)( ( pos.z ) / result->apix.z );	\
\
\
		/*calc end of walk*/ \
		gridend.x = gridid.x + walk.x;\
		gridend.y = gridid.y + walk.y;\
		gridend.z = gridid.z + walk.z;\
\
		gridid.x = gridid.x - walk.x + 1;\
		gridid.y = gridid.y - walk.y + 1;\
		gridid.z = gridid.z - walk.z + 1;\
		\
		/* calc factor */\
		mass = factor[atom]; \
		type_ = type[atom];\
\
\
		/*TODO: check size*/\
\
		for ( k = gridid.z;\
		      k < gridend.z;\
		      k++ ) {\
		      	/*slice index*/\
		      	idxx = k * dim3;\
			for ( j = gridid.y;\
			      j < gridend.y; \
			      j++ ) {\
			      	/*line index*/\
			      	idx = idxx + j * dim2+gridid.x;\
				for ( i = gridid.x;\
				      i < gridend.x; \
				      i++ ) {\
				      	/*distance vector*/\
					r.x = result->apix.x * i  - pos.x;\
					r.y = result->apix.y * j  - pos.y;\
					r.z = result->apix.z * k  - pos.z;\
					/*radius*/\
					lenr = r.x * r.x + r.y * r.y + r.z * r.z;\
\
					/*value*/\
					density = kernel(lenr,type_,sigma);\
					density *= mass;\
\
					/*add density to map and get next gridpoint*/ \
					result->data[idx++] += density;\
				}\
			}\
		}\
	}\
	return result;



PDen_t * pDenRenderPeng(PDen_t * result, real width, real sigma, size_t n, int * type, real * coords, real * factor)
{
	render(pDenPengKernel)
}


PDen_t * pDenRenderGauss(PDen_t * result, real width, real sigma, size_t n, int * type, real * coords, real * factor)
{
	render(pDenGaussKernel)
}













#define render2(kernel) \
	size_t atom; \
	size_t idx, idxx; \
	size_t i,j,k; \
	struct {size_t x,y,z;} gridid; \
	struct {size_t x,y,z;} gridend; \
	struct {size_t x,y,z;} walk; \
	struct {real   x,y,z;} pos; \
	struct {real   x,y,z;} r; \
	real lenr; \
	real density; \
	real mass; \
 \
	/*walk radius*/ \
	walk.x = ((size_t)(0.5 * width/result->apix.x))+1; \
	walk.y = ((size_t)(0.5 * width/result->apix.y))+1;  \
	walk.z = ((size_t)(0.5 * width/result->apix.z))+1; \
 \
	/*render loop*/ \
	for(atom=0;atom<n;atom++){ \
		mass = 1.0; \
		 \
		getAtomPosition\
 \
		/*get upper left grid point*/ \
		gridid.x = (size_t)( pos.x - result->origin.x / result->apix.x );	 \
		gridid.y = (size_t)( pos.y - result->origin.y / result->apix.y );	 \
		gridid.z = (size_t)( pos.z - result->origin.z / result->apix.z );	 \
 \
		/*calc end of walk*/ \
		gridend.x = gridid.x + walk.x; \
		gridend.y = gridid.y + walk.y; \
		gridend.z = gridid.z + walk.z; \
 \
		gridid.x = gridid.x - walk.x + 1; \
		gridid.y = gridid.y - walk.y + 1; \
		gridid.z = gridid.z - walk.z + 1; \
		 \
 \
		/* calc factor*/ \
		mass = cfactor[atom] * dfactor[atom];  \
 \
		/*TODO: check size*/ \
 \
		for ( k = gridid.z; \
		      k < gridend.z; \
		      k++ ) { \
		      	/*slice index*/ \
		      	idxx = k * result->size.x * result->size.y; \
			for ( j = gridid.y; \
			      j < gridend.y;  \
			      j++ ) { \
			      	/*line index*/ \
			      	idx = idxx + j * result->size.x; \
				for ( i = gridid.x; \
				      i < gridend.x;  \
				      i++ ) { \
				      	/*distance vector*/ \
					r.x = result->apix.x * i + result->origin.x - pos.x; \
					r.y = result->apix.y * i + result->origin.y - pos.y; \
					r.z = result->apix.z * i + result->origin.z - pos.z; \
					/*radius*/ \
					lenr = r.x * r.x + r.y * r.y + r.z * r.z; \
					lenr = lenr * isqrt_(lenr); \
 \
					/*value*/ \
					density = kernel(lenr,type[atom],sigma); \
					density *= mass; \
 \
					/*add density to map and get next gridpoint */ \
					result->data[idx++] += density; \
				} \
			} \
		} \
	} \
	return result; 

PDen_t * pDenRenderPeng2(PDen_t * result, real width, real sigma, size_t n, int * type, real * coords, real * cfactor, real *dfactor)
{
	render2(pDenPengKernel)
}

PDen_t * pDenRenderGauss2(PDen_t * result, real width, real sigma, size_t n, int * type, real * coords, real * cfactor, real *dfactor)
{
	render2(pDenGaussKernel)
}
