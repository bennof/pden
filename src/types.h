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


#endif
