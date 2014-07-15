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



const char help[] = 
"HELP: \n"
"PDApplySF: applys Babinet's approximation to modify the structure factor\n\n"
"Usage: psapplysf <args> input\n\n"
"Arguments:\n"
" -h               print this help and exit\n"
" -o  <file>       output density file (mrc/xplor)\n"
" -c  <file>       calculated density file (mrc/xplor)\n"
" -m  <file>       mask density file (mrc/xplor) (default: calculated density file)\n"
" -kg <float>      set global constant scaling factor (default: 1.)\n"
" -Bg <float>      set global B-factor (default: 3.0)\n"
" -ks <float>      set solvent constant scaling factor (default: 0.35 )\n"
" -Bs <float>      set solvent B-factor (default: 46.0 )\n"
" -norm            normalize density ( m=0.0 s=1.0)\n\n"
"Input:\n"
"ignored\n";

int main(int argc, char * argv[] ){
	int i = 1;
	const char *outname = 0;
	const char *calcname = 0;
	const char *maskname = 0;
	PDen_t * map, * calc, *mask;
	real kg=1.,ks=0.35,Bg=3.,Bs=46.;
	char NORMALIZE=0;

	calc  = pDenNew();

	//arguments
	ARGSTART
	  argstring(outname,"-o")
	  argstring(calcname,"-c")
	  argstring(maskname,"-m")
	  argfloat(kg,"-kg")
	  argfloat(Bg,"-Bg")
	  argfloat(ks,"-ks")
	  argfloat(Bs,"-Bs")
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

	if(NORMALIZE){
		pDenNormalize(calc);
		pDenNormalize(mask);
	}
	
	map = pDenNewFrom( calc);
	pDenAlloc(map);

		
	say("Babinet: F_out = k_global * exp(-B_global * s^2 /4) * (F_calc-k_solvent * exp(-B_solvent * s^2 /4) * F_mask)\n");
	
	pDenApplySF_Babinet(map,calc,mask,kg,Bg,ks,Bs);
	pDenFFTNormal(map,-1);

	if(outname){
		writeMap(map,outname);
	}

	return 0;
}
