#include "pden.in.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


/**
 * Modes of a mrc file data section
 * @{
 */
#define MRC_MODE_CHAR                               0 /**< image : signed 8-bit bytes range -128 to 127 */
#define MRC_MODE_SHORT                              1 /**< image : 16-bit halfwords */
#define MRC_MODE_FLOAT                              2 /**< image : 32-bit reals */
#define MRC_MODE_COMPLEX_INT16                      3 /**< transform : complex 16-bit integers */
#define MRC_MODE_COMPLEX_FLOAT                      4 /**< transform : complex 32-bit reals */
#define MRC_MODE_UNKOWN                             5 /**< not used */
#define MRC_MODE_USHORT                             6 /**< image : unsigned 16-bit range 0 to 65535 */
/** @} */


/**
 * MRC header structure
 */
typedef struct MRC_Header_struct {
        int nx;                         /**< number of columns (fastest changing in map) */
        int ny;                         /**< number of rows */
        int nz;                         /**< number of sections (slowest changing in map) */
        int mode;                       /**< data type @see MRC Modes */
        int nxstart;                    /**< number of first cols in map (Default = 0) */
        int nystart;                    /**< number of first rows in map */
        int nzstart;                    /**< number of first section in map */
        int mx;                         /**< number of intervals along X */
        int my;                         /**< number of intervals along Y */
        int mz;                         /**< number of intervals along Z */
        float cella[3];                 /**< cell dimensions in angstroms */
        float cellb[3];                 /**< cell angles in degrees */
        int mapc;                       /**< axis corresp to cols (1,2,3 for X,Y,Z) */
        int mapr;                       /**< axis corresp to rows (1,2,3 for X,Y,Z) */
        int maps;                       /**< axis corresp to sections (1,2,3 for X,Y,Z) */
        float min;                      /**< minimum density value */
        float max;                      /**< maximum density value */
        float mean;                     /**< mean density value */
        int ispg;                       /**< pace group number 0 or 1 (default=0) */
        int nsymbt;                     /**< number of bytes used for symmetry data (0 or 80) */
        int extra[25];                  /**< extra space used for anything   - 0 by default */
        float origin[3];                /**< origin in X,Y,Z used for transforms */
        char map[4];                    /**< character string 'MAP ' to identify file type */
        int machst;                     /**< machine stamp */
        int rms;                        /**< rms deviation of map from mean density */
        int nlabl;                      /**< number of labels being used */
        char label[10][80];             /**< 10 80-character text labels */
} MRC_Header_t;

static void setMRCHeader(MRC_Header_t *header,size_t dim[], real apix[], real origin[], size_t cell[], size_t start[])
{
        header->nx        = ( unsigned int ) dim[0];
        header->ny        = ( unsigned int ) dim[1];
        header->nz        = ( unsigned int ) dim[2];
        header->mode      = MRC_MODE_FLOAT;
        header->nxstart   = start[0];
        header->nystart   = start[1];
        header->nzstart   = start[2];
        header->mx        = cell[0];
        header->my        = cell[1];
        header->mz        = cell[2];
        header->cella[0]  = ( float ) (apix[0] * cell[0]);
        header->cella[1]  = ( float ) (apix[1] * cell[1]);
        header->cella[2]  = ( float ) (apix[2] * cell[2]);
        header->cellb[0]  = 90.0f;
        header->cellb[1]  = 90.0f;
        header->cellb[2]  = 90.0f;
        header->mapc      = 1;
        header->mapr      = 2;
        header->maps      = 3;
        header->min       = 0.0f;
        header->max       = 0.0f;
        header->mean      = 0.0f;
        header->ispg      = 0;
        header->nsymbt    = 0;
        memset(header->extra,0,4*25);
        header->origin[0] = ( float ) origin[0];
        header->origin[1] = ( float ) origin[1];
        header->origin[2] = ( float ) origin[2];
        header->map[0]    = 'M';
        header->map[1]    = 'A';
        header->map[2]    = 'P';
        header->map[3]    = 0;
        header->machst    = 0x12345678;
        header->rms       = 0.0f;
        header->nlabl     = 0;
}

 

int pDenReadMRC(PDen_t * this, const char * filename, const int mode)
{
	int fd;
	MRC_Header_t head;
	size_t i,R,l;
	real *r , *d;
	char *c;
	short *s;
	unsigned short *u;
	#ifdef DOUBLE
	float *f;
	#endif

	

	fd = open(filename,O_RDONLY); //open file	
	if ( fd < 0 ) 
		return 1;

	if ( 1024 != read ( fd, ( void * ) &head, 1024 ) ){ // read header
		close(fd);
		return 2;
	}

	// get info from header
	this->size.x = ( size_t ) head.nx;
	this->size.y = ( size_t ) head.ny;
	this->size.z = ( size_t ) head.nz;
	
	this->apix.x = ( real ) head.cella[0] / head.mx;
	this->apix.y = ( real ) head.cella[1] / head.my;
	this->apix.z = ( real ) head.cella[2] / head.mz;

	this->origin.x = ( real ) head.origin[0];
	this->origin.y = ( real ) head.origin[1];
	this->origin.z = ( real ) head.origin[2];

	this->cell.x = ( size_t ) head.mx;
	this->cell.y = ( size_t ) head.my;
	this->cell.z = ( size_t ) head.mz;

	this->start.x = ( ssize_t ) head.nxstart;
	this->start.y = ( ssize_t ) head.nystart;
	this->start.z = ( ssize_t ) head.nzstart;

	//Chimera specific origin
	// TODO: add mode flag and correct origin
	if ( ( this->origin.x == 0. ) &&
	     ( this->origin.y == 0. ) &&
	     ( this->origin.z == 0. ) ) {
		this->origin.x = ( real ) ( head.nxstart * this->apix.x );
		this->origin.y = ( real ) ( head.nystart * this->apix.y );
		this->origin.z = ( real ) ( head.nzstart * this->apix.z );
	}

	this->n = this->size.x * this->size.y * this->size.z;

	//alloc
	if ( this->data ) free( this->data );
	r = ( real * ) malloc ( ( this->size.x + 2 ) * this->size.y * this->size.z * sizeof( real ) );
	
	switch ( head.mode ) {
	case MRC_MODE_CHAR :
                l = this->n;
                R = read ( fd, r, l );
                if ( R != l ){
                        free( r );
			close(fd);
			return 3; 
                }
                else {
		    l = this->n;
                    c = ( char * ) r;
                    c = &c[l];
                    d = &r[l];
                    for ( i = 0; i < l; i++) {
                        *d = ( ( real ) *c ) / 128.;
                        c--;
                        d--;
                    }
                }
                break;
	case MRC_MODE_SHORT:
                l = this->n * sizeof(short);
                R = read ( fd, r, l);
                if(R!=l){
                        free ( r );
			close(fd);
                        return 3;
                }
                else {
		    l = this->n;
                    s = ( short * ) r;
                    s = &s[l];
                    d = &r[l];
                    for ( i = 0; i < l ; i++) {
                        *d = ( ( real ) *s ) / 32767.;
                        d--;
                        s--;
                    }
                }
                break;
        case MRC_MODE_FLOAT:
                l = this->n * sizeof(float);
                R = read ( fd, r, l);
                if( R != l ){
                        free ( r );
			close(fd);
                        return 3;
                }
		#ifdef DOUBLE
                else {
		    l = this->n;
                    f = ( float * ) r;
                    f = &f[l-1];
                    d = &r[l-1];
                    for ( i = 0; i < l ; i++) {
                        *d = ( ( float ) *f );
                        d--;
                        f--;
                    }
		}
		#endif
                break;
        case MRC_MODE_COMPLEX_INT16:
        case MRC_MODE_COMPLEX_FLOAT:
		close(fd);
		return 3;
        case MRC_MODE_USHORT:
                l = this->n * sizeof(short);
                R = read( fd, r, l);
                if ( R != l ) {
                        free ( r );
			close(fd);
			return 3;
                }
                else {
		    l = this->n;
                    u = (unsigned short*) r;
                    u = &u[l];
                    d = &r[l];
                    for(i=0;i<l;i++){
                        *d = ( ( float ) *u ) / 65536.;
                        d--;
                        u--;
                    }
                }
                break;
        default:
		close(fd);
                return 3;
        }
	this->data = r;
	return 0;	
}



int pDenWriteMRC(PDen_t * this, const char * filename, const int mode)
{
	int fd;
	size_t size;
	#ifdef DOUBLE
	size_t i;
	#endif
	MRC_Header_t head;
	float *data;

	setMRCHeader( &head , ( size_t * ) &(this->size), ( real * ) &(this->apix), ( real * ) &(this->origin),( size_t * ) &(this->cell),( size_t * ) &(this->start));
	size=this->n*sizeof(float);


	fd = open(filename, O_WRONLY|O_CREAT|O_TRUNC,0644);
	if ( fd < 0 ) 
		return 1;

	if ( 1024 != write ( fd, ( const void * ) &head, 1024 ) ){ // write header
		close(fd);
		return 2;
	}
	
	#ifdef DOUBLE	// type conversion if double
	data = ( float * ) malloc( size );
	for(i=0;i<this->n;i++)
		data[i] = ( float ) this->data[i];
	#else
	data=this->data;
	#endif
	
	if(size!=write(fd,data,size)){ //write data
                close(fd);
                return 3;
        }

	close(fd);
	#ifdef DOUBLE
	free(data);
	#endif
	return 0;
}
