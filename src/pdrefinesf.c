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


const char help[] = 
"HELP: \n"
"PDRefineSF: refines constants for structure factor using Babinet's approximation\n\n"
"Usage: pdrefinesf <args> input\n\n"
"Arguments:\n"
" -h               print this help and exit\n"
" -c  <file>       calculated density file (mrc/xplor)\n"
" -m  <file>       mask density file (mrc/xplor) (default: calculated density file)\n"
" -o  <file>       output density file (mrc/xplor)\n"
" -ps              powerspectrum output file (plain)\n"
" -kg <float>      set global constant scaling factor (default: 1.)\n"
" -Bg <float>      set global B-factor (default: 3.0)\n"
" -ks <float>      set solvent constant scaling factor (default: 0.35 )\n"
" -Bs <float>      set solvent B-factor (default: 46.0 )\n"
" -s  <int>        set max steps (default: 100)\n"
" -p  <float>      set precision (default: 0.0001)\n"
" -norm            normalize density ( m=0.0 s=1.0)\n\n"
"Input:\n"
"mrc/xplor file as target\n";


int main(int argc, char * argv[] ) 
{
	int i = 1;
	FILE *fp;
	const char *calcname = 0;
	const char *maskname = 0;
	const char *outname = 0;
	const char *psname = 0;
	size_t n;
	real *ps;
	PDen_t * map, * calc, *mask;
	real kg=1.,ks=0.35,Bg=30.,Bs=46.;
	real prec=0.00001;
	real s;
	size_t steps=200;
	char NORMALIZE=0;

	INIT("PDRefineSF");

	map  = pDenNew();
	calc = pDenNew();

	//arguments
	ARGSTART
	  argstring(calcname,"-c")
	  argstring(maskname,"-m")
	  argstring(outname,"-o")
	  argstring(psname,"-ps")
	  argfloat(kg,"-kg")
	  argfloat(Bg,"-Bg")
	  argfloat(ks,"-ks")
	  argfloat(Bs,"-Bs")
	  argfloat(prec,"-p")
	  argint(steps,"-s")
	  argtrue(NORMALIZE,"-norm")
	ARGFINALIZE

	if(!calcname){
		die("No calculated density given (-c)");
	}

	readMap(calc,calcname);

	if(maskname){
		mask = pDenNew();
		readMapO(mask,maskname);
		if(mask->data==0){
			pDenDelete(mask);
			mask = calc;
		}
	}
	else {
		mask = calc;;
	}

	if(i<argc){
		readMap(map,argv[i])
	}
	else {
		die("No input");
	}

	if(NORMALIZE){
		pDenNormalize(map);
		pDenNormalize(calc);
		pDenNormalize(mask);
	}
	say("Babinet: F = k_global * exp(-B_global * s^2 /4) * (F_calc-k_solvent * exp(-B_solvent * s^2 /4) * F_mask)\n");
	if(pDenRefineSF_Babinet(map,calc,mask,&kg,&Bg,&ks,&Bs,prec,steps))
		error("CG not converged");
	info("Summary stucture factor refinemnt:");
	info("k_global  = %lf\tB_global  = %lf",kg,Bg);
	info("k_solvent = %lf\tB_solvent = %lf",ks,Bs);


	//apply sf
	if(outname || psname){
		pDenApplySF_Babinet(calc,calc,mask,kg,Bg,ks,Bs);
	}
	else {
		info("Arguments: -kg %lf -Bg %lf -ks %lf -Bs %lf",kg,Bg,ks,Bs);
	}

	if(psname){
		n=calc->size.x/2;
		if(!(ps = malloc(n*sizeof(real))))
			die("Malloc failed");
		pDenCalcPS(calc,ps,n);
		OUTPUTFILE(fp,psname);
		fprintf(fp,"# PDen " __PDEN_VERSION__  " (" __DATE__ " -  " __TIME__")\n");	
		fprintf(fp,"# PDen Power Spectrum\n");
		fprintf(fp,"# File:    %s\n",argv[i]);
		fprintf(fp,"# Records: %lu\n",n);

		s  =  1 / (4. * map->apix.x*map->apix.x);
		s +=  1 / (4. * map->apix.y*map->apix.y);
		s +=  1 / (4. * map->apix.z*map->apix.z);
		s  =  1 / (4. * map->apix.x*map->apix.x);

		s  = sqrt(s); 
		s /= (real)n;

		for(i=0;i<n;i++){
			fprintf(fp,"%lf \t %lf\n",((real)i)*s,ps[i]);
		}
		OUTPUTFILECLOSE(fp);	
	}

	if(outname) {
		pDenFFTNormal(map,-1);
		 writeMap(map,outname);
	}


	return 0;
}
