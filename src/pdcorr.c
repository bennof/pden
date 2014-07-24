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
"PDCorr: calculates map correlation\n\n"
"Usage: pdrefinesf <args> input\n\n"
"Arguments:\n"
" -h               print this help and exit\n"
" -e               estimate correlation error\n"
" -n               bootstrapping steps (default: 100)\n"
"Input:\n"
"mrc/xplor two files\n";

int main(int argc,char *argv[])
{
	int i,mi1, mi2;
	PDen_t * map1, *map2;
	real m1,m2,s1,s2;
	real corr, e_corr;
	size_t nboot = 100;
	char BSCORR=0;

	//arguments
	ARGSTART
	argtrue(BSCORR,"-e")
	argint(nboot,"-n")
	ARGFINALIZE

	if( i < argc ) {
		map1 = pDenNew();
		readMap(map1,argv[i])
		mi1=i;
		i++;
	}

	if( i < argc ) {
		map2 = pDenNew();
		readMap(map2,argv[i])
		mi2=i;
		i++;
	}

	pDenMeanSd(map1,&m1,&s1);	
	info("Map: %s",argv[mi1]);
	info("Mean: %5.2lf \tSD:  %5.2lf ",m1,s1);
	pDenMeanSd(map2,&m2,&s2);	
	info("Map: %s",argv[mi2]);
	info("Mean: %5.2lf \tSD:  %5.2lf ",m2,s2);
	info("Correlation: %5.3lf",pDenCorr(map1, map2));
	if( BSCORR ) {
		info("Estimating correlation (bootstrapping: %li)",nboot);
        	pDenMapCorrError(map1, map2, &corr, &e_corr, nboot);
		info("Correlation: %5.3lf (+/-) %5.3lf",corr,e_corr);
	}
	return 0;
}
