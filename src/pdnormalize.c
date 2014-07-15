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
"PDnormalize: normalize density maps\n\n"
"Usage: pdnormalize <args> input\n\n"
"Arguments:\n"
" -o               output file (mrc)\n\n"
"Input:\n"
"mrc/xplor files\n";


int main(int argc, char * argv[] ) 
{
	int i=1;
	const char *output = "out.mrc";
	PDen_t * map;

	INIT("PDNormalize");

	map = pDenNew();

	//arguments
	ARGSTART
	  argstring(output,"-o")
	ARGFINALIZE

	//input & process 
	if(i<argc){
		readMap(map,argv[i])
	}
	else {
		die("No input");
	}

	pDenNormalize(map);

	writeMap(map,output);
	return 0;
}


