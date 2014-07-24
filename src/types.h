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

#ifndef __TYPES_H
#define __TYPES_H

#include <math.h>

/**
 * @file types.h
 * @author Benjamin Falkner
 * @brief File containing mathematical types to be used internally
 * @warning internal file - not installed
 */



/**
 * @brief 3 dimensional vector 
 *
 * Three dimensional real vector type 
 */
typedef struct {
	real    x; /**< x-coordinate */
	real    y; /**< y-coordinate */
	real    z; /**< z-coordinate */
}    vec3;

/**
 * @brief 3 dimensional vector 
 *
 * Three dimensional real vector type 
 */
typedef struct {
	vec3    x; /**< x-column */
	vec3    y; /**< y-column */
	vec3    z; /**< z-column */
}    mat3;

/**
 * @brief 3 dimensional vector 
 *
 * Three dimensional unsigned integer  vector type 
 */
typedef struct {
	size_t  x; /**< x-coordinate */ 
	size_t  y; /**< y-coordinate */
	size_t  z; /**< z-coordinate */
}   uvec3;

/**
 * @brief 3 dimensional vector 
 *
 * Three dimensional signed integer vector type 
 */
typedef struct {
	ssize_t  x; /**< x-coordinate */ 
	ssize_t  y; /**< y-coordinate */
	ssize_t  z; /**< z-coordinate */
}   ivec3;

/**
 * @brief complex number 
 *
 * Complex number type type 
 */
typedef struct {
	real    r; /**< real component */
	real    i; /**< imaginary component */
}     cplx;

/**
 * @brief quaternion number 
 *
 * Quaternion type 
 */
typedef struct {
	real    r; /**< real component */
	real    i; /**< i component */ 
	real    j; /**< j component */
	real    k; /**< k component */
} quat;


static inline real quatAbs(quat q)
{
	register float h;
	h=q.r*q.r+q.i*q.i+q.j*q.j+q.k*q.k;
	return vsqrt(h);
}

static inline quat quatAdd(quat a,quat b)
{
	register quat q;
	q.r = a.r + b.r;
	q.i = a.i + b.i;
	q.j = a.j + b.j;
	q.k = a.k + b.k;
	return q;
}

static inline quat quatSub(quat a,quat b)
{
	register quat q;
	q.r = a.r - b.r;
	q.i = a.i - b.i;
	q.j = a.j - b.j;
	q.k = a.k - b.k;
	return q;
}

static inline quat quatMult(quat a,quat b)
{
	register quat q;
	q.r = a.r * b.r - a.i * b.i - a.j * b.j - a.k * b.k;
	q.i = a.r * b.i + a.i * b.r + a.j * b.k - a.k * b.j;
	q.j = a.r * b.j - a.i * b.k + a.j * b.r + a.k * b.i;
	q.k = a.r * b.k + a.i * b.j - a.j * b.i + a.k * b.r;
	return q;
}

static inline quat quatRot(const char axis,float angle)
{
	register quat q;
	q.r = cos( angle / 2 );
	q.i = 0.0;
	q.j = 0.0;
	q.k = 0.0;

	switch (axis) {
	case 'X':
		q.i = sin( angle * 0.5 );	
		break;
	case 'Y':
		q.j = sin( angle * 0.5 );	
		break;
	case 'Z':
		q.k = sin( angle * 0.5 );	
		break;
	}
	return q;
}

static inline quat quatFromEuler(const char order[],float a, float b, float c)
{
	quat q;
	q = quatRot(order[0],a);
	q = quatMult(q,quatRot(order[1],a));
	q = quatMult(q,quatRot(order[2],a));
	return q;
}

static inline mat3 mat3FromQuat(quat q)
{
	mat3 r;
        r.x.x = 1.0f-2.0f*(q.j*q.j+q.k*q.k);
        r.x.y = 2.0f*(q.r*q.k+q.i*q.j);
        r.x.z = 2.0f*(q.i*q.k-q.r*q.j);
    
        r.y.x = 2.0f*(q.i*q.j-q.r*q.k);
        r.y.y = 1.0f-2.0f*(q.i*q.i+q.k*q.k);
        r.y.z = 2.0f*(q.r*q.i+q.j*q.k);

        r.z.x = 2.0f*(q.r*q.j+q.i*q.k);
        r.z.y = 2.0f*(q.j*q.k-q.r*q.i);
        r.z.z = 1.0f-2.0f*(q.i*q.i+q.j*q.j);
	return r;
}

#endif
