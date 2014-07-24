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
"PDCalcSF: calculate radial structure factor of a density maps\n\n"
"Usage: pdcalcpsf <args> input\n\n"
"Arguments:\n"
" -o               output file (plain)\n"
" -n <int>         number of intervals (default: grid size )\n"
" -norm            normalize density ( m=0.0 s=1.0)\n\n"
"Input:\n"
"mrc/xplor files\n";


int main(int argc, char * argv[] ) 
{
	int i = 1;
	const char *output = "out.dat";
	PDen_t * map;
	size_t n = 0;
	real *ps;
	FILE *fp;
	char NORMALIZE=0;

	real s;

	INIT("PDCalcPS");

	map = pDenNew();

	//arguments
	ARGSTART
	  argint(n,"-n")
	  argstring(output,"-o")
	  argtrue(NORMALIZE,"-norm")
	ARGFINALIZE


	//input & process 
	if(i<argc){
		readMap(map,argv[i])
	}
	else {
		die("No input");
	}


	if(n==0){
		n=map->size.x/2;
	}

	if(NORMALIZE)
		pDenNormalize(map);

	if(!(ps = malloc(n*sizeof(real))))
		die("Malloc failed");
	
	pDenCalcSF(map,ps,n);

	OUTPUTFILE(fp,output);
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
	return 0;
}


