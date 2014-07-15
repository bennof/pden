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

#include "pden_sfrefine.h"



//OPTIMIZER
static real lineSearchArmijo(const PDen_t *X[], real param0[], SFValueFunc_t func, real grad[], real val, const size_t npar){
	real lambda = 2.;
	size_t i, cc=0;
	real eval,nval; //linear estimation value, new value
	real x1[npar];
	real alpha=0.01;

	debug("Armijo rule line search loop");
	while(cc<100){
		eval = 0.0;
		for(i=0;i<npar;i++){
			x1[i] = param0[i]+lambda*grad[i]; // candidate
			//debug("%lf = %lf + %lf * %lf", x1[i],x0[i],lambda,grad[i]);
			eval += lambda*grad[i]*lambda*grad[i];
		}
		//debug("Armijo rule[%3lu]: (%lf,%lf,%lf)->(%lf,%lf,%lf)",cc,x0[0],x0[1],x0[2],x1[0],x1[1],x1[2]);

		
		eval = val - alpha * sqrt(eval);
		nval = func(X,(const real*)x1);// value at position

		debug("Armijo rule[%3lu]: %lf<=%lf (lambda=%lf)",cc,nval,eval,lambda);
		if(nval <= eval){ // success
			for(i=0;i<npar;i++){
				param0[i]=x1[i]; //update position
			}
			break;
		} 
		else {
			lambda *= .5;	
		}
		cc++;
	}
	if(cc>=100) {
		error("Armijo rule line search not converged (%lf<=%lf)(old %lf)",nval,eval,val);
		return val;

	}
	return nval;
}


//Conjugated Gradients  
int CGRefine (const PDen_t *X[], real param0[], SFValueFunc_t func, SFGradFunc_t gradfunc, const size_t npar, const real prec, const size_t steps)
{
	size_t i;
	size_t c = 0;         //loop counter
	real d[npar];         // direction
	real g[npar];         // old gradient
	real grad[npar]; 
	real val,nval,h;
	real beta;

	
	// Conjugated Gradients
	val = func(X,param0); 
	gradfunc(X,param0,grad);	
	while(c<steps){ // recursion loop
		#if defined (DEBUG) || defined (VERBOSE)
		say("%9lu ",c);
		for(i=0;i<npar;i++)
			fprintf(stdout,"%9.4lf ",param0[i]);
		fprintf(stdout,"%9.6lf\n",val);
		#endif

		//Armijo rule line search loop
                nval = lineSearchArmijo(X,param0,func,grad,val,npar);
		if (nval >= val || isnan(nval))
			return 1;


		debug("Test convergence: %lf ~= %lf (old/new)",val,nval);
		if ( fabs(val-nval) < prec)
			break;
		val=nval;

		//needed fot beta
		h = 0.;
		for(i=0;i<npar;i++){
			//debug("h=%lf grad=%lf",h,grad[i]);
			h += grad[i]*grad[i];
			g[i] = grad[i];
		}

		// gradient
		debug("calculate gradient");
		gradfunc(X,param0,grad);	
		
		debug("calculate beta");
		beta = 0.;
		// beta (Polak–Ribière)
		for(i=0;i<npar;i++){
			beta += grad[i]*(grad[i]-g[i]);
		}
		debug("beta=%lf/%lf",beta,h);
		beta = beta/h; 

                // direction
		debug("apply step (beta=%lf)",beta);
		for(i=0;i<npar;i++){
			d[i] += grad[i] + beta * d[i];	
		}

		c++;
	}
	if(c==steps)
		return 1;
	debug("CG converged\n");	
	return 0;
}



