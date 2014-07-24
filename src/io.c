#include "pden.in.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


int pDenRead (PDen_t * this, const char * filename, const int mode)
{
	if ( strcasestr ( filename, ".xplor") )
		return pDenReadXPLOR(this,filename,mode);
	else if ( strcasestr ( filename, ".mrc") )
		return pDenReadMRC(this,filename,mode);
	else 
		return 1;
}

int pDenWrite (PDen_t * this, const char * filename, const int mode)
{
	if ( strcasestr ( filename, ".xplor") )
		return pDenWriteXPLOR(this,filename,mode);
	else if ( strcasestr ( filename, ".mrc") )
		return pDenWriteMRC(this,filename,mode);
	else 
		return 1;
}
