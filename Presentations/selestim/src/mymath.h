/*
 Copyright INRA
 author: Renaud VITALIS (2013)
 
 renaud.vitalis@inra.fr
 
 This file is part of SelEstim.
 
 SelEstim is a computer program whose purpose is to is to detect
 and measure selection from gene frequency data.
 
 SelEstim is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 SelEstim is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MYMATH_H
#define _MYMATH_H

#include <float.h>                                                               /* This is to get the definition of DBL_MIN and DBL_MAX */
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dcmt0.6.1b/include/dc.h"
#include "defs.h"
#include "rand.h"

#define ISNAN(x)           (isnan(x)!=0)                                            /* Added RV 13-10-2005 from nmath.h */
#define R_FINITE(x)        isfinite(x)
#define ML_ERROR(x)
#define ML_POSINF	         (1.0 / 0.0)                                              /* Added RV 13-10-2005 from nmath.h */
#define ML_NEGINF	         ((-1.0) / 0.0)                                           /* Added RV 13-10-2005 from nmath.h */
#define ML_NAN		         (0.0 / 0.0)                                              /* Added RV 13-10-2005 from nmath.h */
#define M_PI_R	           3.141592653589793238462643383280                         /* pi */ /* Added RV 13-10-2005 from Rmath.h */

#ifndef M_LN2
#define M_LN2              0.693147180559945309417232121458                         /* ln(2) */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	     0.918938533204672741780329736406                         /* log(sqrt(2*pi)) == log(2*pi)/2 */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2     0.225791352644727432363097614947                         /* log(sqrt(pi/2)) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2	         0.301029995663981195213738894724                         /* log10(2) */
#endif

#define trunc(x)           ((x) < 0 ? (ceil((x))) : (floor((x))))                   /* Added RV 28-10-2005 */
#define ML_UNDERFLOW	     (DBL_MIN * DBL_MIN)                                      /* Added RV 13-10-2005 from nmath.h */
#define ML_ERR_return_NAN  {return ML_NAN;}                                         /* Added RV 13-10-2005 from nmath.h */
#define n_max              (100)
#define expmax	           (DBL_MAX_EXP * M_LN2)
#define ML_TREAT_psigam(_IERR_)	 \
if(_IERR_ != 0) 			           \
return ML_NAN
#define v_w_from__u1_bet(AA)		 \
v = beta * log(u1 / (1.0 - u1)); \
if (v <= expmax)			           \
w = AA * exp(v);	        	     \
else                             \
w = DBL_MAX

double gammafn(double x);
double lgammafn(double x);
double lgammacor(double x);
double chebyshev_eval(double x, const double *a, const int n);
double stirlerr(double n);
void dpsifn(double x, int n, int kode, int m, double *ans, int *nz, int *ierr);
double digamma(double x);
int	Rf_i1mach(int);
double Rf_d1mach(int);
int int_min(int x, int y);
int int_max(int x, int y);
double dbl_min(double x, double y);
double dbl_max(double x, double y);
double hypergeometric_1_F_1(double a,double b,double z);
double F_ST(int ***x,int **n,int nbr_loci,int nbr_demes);
double F_ST_poolseq(int ***x,int **n,int nbr_loci,int nbr_demes,int *pool_size);
double log_uniform_prior(double x);
double log_beta_prior(double *x,double shape1,double shape2);
double log_exponential_prior(double x,double a);
double log_inverse_gamma_prior(double x,double alpha,double beta);
int signgam;

#endif
