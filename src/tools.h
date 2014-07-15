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

#ifndef __tools_h
#define __tools_h


#ifdef DEBUG 
#define VERBOSE 1
#endif

#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#include <stdio.h>
#include <string.h>

#define error( ... ) \
{\
fprintf(stderr,"#PDEN ERROR> ");\
fprintf(stderr, __VA_ARGS__ );\
fprintf(stderr," [%s:%i]\n",__FILE__,__LINE__);\
}

#define die( ... )\
{\
fprintf(stderr,"#PDEN ERROR> ");\
fprintf(stderr, __VA_ARGS__ );\
fprintf(stderr, "[%s:%i]\n",__FILE__,__LINE__);\
fprintf(stderr,"#PDEN EXIT\n");\
exit(1);\
}

#define info( ... ) \
{\
fprintf(stdout,"#PDEN> ");\
fprintf(stdout, __VA_ARGS__ );\
fprintf(stdout,"\n");\
}


#ifdef VERBOSE
#define say( ... ) \
{\
fprintf(stdout,"#PDEN> " __VA_ARGS__ );\
}
#else
#define say( ... )
#endif

#ifdef DEBUG 
#define debug( ... ) \
{\
fprintf(stdout,"#PDEN DEBUG> ");\
fprintf(stdout, __VA_ARGS__ );\
fprintf(stdout," [%s:%i]\n",__FILE__,__LINE__);\
}
#else
#define debug( ... )
#endif



//READ MAP
#define readMap(to,name) \
{\
	info("read: %s",name);\
	if(strstr(name,".mrc")){\
		if(pDenReadMRC(to,name,0))\
			die("failed to read: %s",name);\
	}\
	else if(strstr(name,".xplor")){\
		if(pDenReadXPLOR(to,name,0))\
			 die("failed to read: %s",name);\
	}\
	else {\
		die("Unknown file format: %s (use mrc/xplor)",name);\
	}\
}

#define readMapO(to,name) \
{\
	info("read: %s",name);\
	if(strstr(name,".mrc")){\
		if(pDenReadMRC(to,name,0))\
			die("failed to read: %s",name);\
	}\
	else if(strstr(name,".xplor")){\
		if(pDenReadXPLOR(to,name,0))\
			 die("failed to read: %s",name);\
	}\
	else {\
		error("Unknown file format: %s (use mrc/xplor)",name);\
	}\
}

#define writeMap(to,name) \
{\
	info("write: %s",name);\
	if(strstr(name,".mrc")){\
		if(pDenWriteMRC(to,name,0))\
			die("failed to write: %s",name);\
	}\
	else {\
		error("Unknown file format: %s (use mrc only)",name);\
	}\
}

//OPENFILE
#define OUTPUTFILE(fp_,name) \
{\
	if(!strcmp(name,"-"))\
		fp_ = stdout;\
	else {\
		info("write: %s",name);\
		if(!(fp=fopen(name,"w")))\
			die("Faile open file: %s",name);\
	}\
}

#define OUTPUTFILECLOSE(fp_) \
{\
	if(fp_!=stdout)\
		fclose(fp_);\
}



//INIT
#define INIT(name) \
info( name " " PACKAGE_STRING  " (" __DATE__ " -  " __TIME__")")

//ARGUMENTS 

#define argfloat(val,name) \
else if(!strcmp(argv[i],name)){\
	val = atof(argv[++i]);\
}

#define argint(val,name) \
else if(!strcmp(argv[i],name)){\
	val = atoi(argv[++i]);\
}

#define argtrue(val,name) \
else if(!strcmp(argv[i],name)){\
	val=1;\
}

#define argfalse(val,name) \
else if(!strcmp(argv[i],name)){\
	val=0;\
}


#define argstring(val,name) \
else if(!strcmp(argv[i],name)){\
	val = argv[++i];\
}

#define ARGFINALIZE \
else if (argv[i][0]=='-'){\
	error("Unkown argument: %s", argv[i]);\
}\
else{ \
	break;\
}\
i++;\
}

#define ARGSTART \
while(i<argc){\
if(!strcmp(argv[i],"-h")){\
puts(help);\
exit(0);\
}\





#endif
