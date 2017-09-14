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

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  Mathlib - A Mathematical Function Library
 *  Copyright (C) 1998  Ross Ihaka
 *  Copyright (C) 2000-7 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

#include "mymath.h"

/*
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double gammafn(double x);
 *
 *  DESCRIPTION
 *
 *    This function computes the value of the gamma function.
 *
 *  NOTES
 *
 *    This function is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    The accuracy of this routine compares (very) favourably
 *    with those of the Sun Microsystems portable mathematical
 *    library.
 *
 *    MM specialized the case of  n!  for n < 50 - for even better precision
 */

double gammafn(double x)
{
  static const double gamcs[42] = {
    +.8571195590989331421920062399942e-2,
    +.4415381324841006757191315771652e-2,
    +.5685043681599363378632664588789e-1,
    -.4219835396418560501012500186624e-2,
    +.1326808181212460220584006796352e-2,
    -.1893024529798880432523947023886e-3,
    +.3606925327441245256578082217225e-4,
    -.6056761904460864218485548290365e-5,
    +.1055829546302283344731823509093e-5,
    -.1811967365542384048291855891166e-6,
    +.3117724964715322277790254593169e-7,
    -.5354219639019687140874081024347e-8,
    +.9193275519859588946887786825940e-9,
    -.1577941280288339761767423273953e-9,
    +.2707980622934954543266540433089e-10,
    -.4646818653825730144081661058933e-11,
    +.7973350192007419656460767175359e-12,
    -.1368078209830916025799499172309e-12,
    +.2347319486563800657233471771688e-13,
    -.4027432614949066932766570534699e-14,
    +.6910051747372100912138336975257e-15,
    -.1185584500221992907052387126192e-15,
    +.2034148542496373955201026051932e-16,
    -.3490054341717405849274012949108e-17,
    +.5987993856485305567135051066026e-18,
    -.1027378057872228074490069778431e-18,
    +.1762702816060529824942759660748e-19,
    -.3024320653735306260958772112042e-20,
    +.5188914660218397839717833550506e-21,
    -.8902770842456576692449251601066e-22,
    +.1527474068493342602274596891306e-22,
    -.2620731256187362900257328332799e-23,
    +.4496464047830538670331046570666e-24,
    -.7714712731336877911703901525333e-25,
    +.1323635453126044036486572714666e-25,
    -.2270999412942928816702313813333e-26,
    +.3896418998003991449320816639999e-27,
    -.6685198115125953327792127999999e-28,
    +.1146998663140024384347613866666e-28,
    -.1967938586345134677295103999999e-29,
    +.3376448816585338090334890666666e-30,
    -.5793070335782135784625493333333e-31
  };

  int i, n;
  double y;
  double sinpiy, value;

  #define ngam 22

  static const double xmin = -170.5674972726612;
  static const double xmax = 171.61447887182298;
  static const double xsml = 2.2474362225598545e-308;
  static const double dxrel = 67108864.;

  if(ISNAN(x)) return x;
    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
  if (x == 0 || (x < 0 && x == (long)x)) {
    ML_ERROR(ME_RANGE);
    return ML_NAN;
  }
  y = fabs(x);
  if (y <= 10) {
    /* Compute gamma(x) for -10 <= x <= 10
     * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
     * first of all. */
    n = x;
    if(x < 0) --n;
    y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
    --n;
    value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
    if (n == 0) return value;/* x = 1.dddd = 1+y */
    if (n < 0) {
      /* compute gamma(x) for -10 <= x < 1 */
      /* exact 0 or "-n" checked already above */
      /* The answer is less than half precision */
      /* because x too near a negative integer. */
      if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) {
        ML_ERROR(ME_PRECISION);
      }
      /* The argument is so close to 0 that the result would overflow. */
      if (y < xsml) {
        ML_ERROR(ME_RANGE);
        if(x > 0) return ML_POSINF;
        else return ML_NEGINF;
      }
      n = -n;
      for (i = 0; i < n; i++) {
        value /= (x + i);
      }
      return value;
    }
    else {
      /* gamma(x) for 2 <= x <= 10 */
      for (i = 1; i <= n; i++) {
        value *= (y + i);
      }
      return value;
    }
  }
  else {
    /* gamma(x) for	 y = |x| > 10. */
    if (x > xmax) {			/* Overflow */
      ML_ERROR(ME_RANGE);
      return ML_POSINF;
    }
    if (x < xmin) {			/* Underflow */
      ML_ERROR(ME_UNDERFLOW);
      return ML_UNDERFLOW;
    }
    if(y <= 50 && y == (int)y) { /* compute (n - 1)! */
      value = 1.;
      for (i = 2; i < y; i++) value *= i;
    }
    else { /* normal case */
      value = exp((y - 0.5) * log(y) - y + M_LN_SQRT_2PI +
      ((2*y == (int)2*y)? stirlerr(y) : lgammacor(y)));
    }
    if (x > 0) return value;
    if (fabs((x - (int)(x - 0.5))/x) < dxrel){
      /* The answer is less than half precision because */
      /* the argument is too near a negative integer. */
      ML_ERROR(ME_PRECISION);
    }
    sinpiy = sin(M_PI_R * y);
    if (sinpiy == 0) {		/* Negative integer arg - overflow */
      ML_ERROR(ME_RANGE);
      return ML_POSINF;
    }
    return -M_PI_R / (y * sinpiy * value);
  }
}


/*
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    extern int signgam;
 *    double lgammafn(double x);
 *
 *  DESCRIPTION
 *
 *    This function computes log|gamma(x)|.  At the same time
 *    the variable "signgam" is set to the sign of the gamma
 *    function.
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    The accuracy of this routine compares (very) favourably
 *    with those of the Sun Microsystems portable mathematical
 *    library.
 */

double lgammafn(double x)

{
  double ans, y, sinpiy;
  static const double xmax = 2.5327372760800758e+305;
  static const double dxrel = 1.490116119384765696e-8;

  signgam = 1;
  #ifdef IEEE_754
    if(ISNAN(x)) return x;
  #endif
  if (x < 0 && fmod(floor(-x), 2.) == 0)
    signgam = -1;
  if (x <= 0 && x == trunc(x)) { /* Negative integer argument */
    ML_ERROR(ME_RANGE);
    return ML_POSINF;/* +Inf, since lgamma(x) = log|gamma(x)| */
  }
  y = fabs(x);
  if (y <= 10) return log(fabs(gammafn(x)));
  if (y > xmax) {
    ML_ERROR(ME_RANGE);
    return ML_POSINF;
  }
  if (x > 0) { /* i.e. y = x > 10 */
  #ifdef IEEE_754
    if(x > 1e17) return(x*(log(x) - 1.));
    else
      if(x > 4934720.) return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
      else
  #endif
        return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
  }
  /* else: x < -10; y = -x */
  sinpiy = fabs(sin(M_PI_R * y));
  ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);
  if(fabs((x - trunc(x - 0.5)) * ans / x) < dxrel) {
	/* The answer is less than half precision because
	 * the argument is too near a negative integer. */
    ML_ERROR(ME_PRECISION);
  }
  return ans;
}


/*
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double lgammacor(double x);
 *
 *  DESCRIPTION
 *
 *    Compute the log gamma correction factor for x >= 10 so that
 *
 *    log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x)
 *
 *    [ lgammacor(x) is called	Del(x)	in other contexts (e.g. dcdflib)]
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    written by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *  SEE ALSO
 *
 *    Loader(1999)'s stirlerr() {in ./stirlerr.c} is *very* similar in spirit,
 *    is faster and cleaner, but is only defined "fast" for half integers.
 */

double lgammacor(double x)
{
  static const double algmcs[15] = {
    +.1666389480451863247205729650822e+0,
    -.1384948176067563840732986059135e-4,
    +.9810825646924729426157171547487e-8,
    -.1809129475572494194263306266719e-10,
    +.6221098041892605227126015543416e-13,
    -.3399615005417721944303330599666e-15,
    +.2683181998482698748957538846666e-17,
    -.2868042435334643284144622399999e-19,
    +.3962837061046434803679306666666e-21,
    -.6831888753985766870111999999999e-23,
    +.1429227355942498147573333333333e-24,
    -.3547598158101070547199999999999e-26,
    +.1025680058010470912000000000000e-27,
    -.3401102254316748799999999999999e-29,
    +.1276642195630062933333333333333e-30
    };
  double tmp;

  #define nalgm 5

  static const double xbig = 94906265.62425156;
  static const double xmax = 3.745194030963158e306;

  if (x < 10)
    ML_ERR_return_NAN
  else if (x >= xmax) {
    ML_ERROR(ME_UNDERFLOW);
    return ML_UNDERFLOW;
  }
  else if (x < xbig) {
    tmp = 10 / x;
    return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
  }
  else return 1 / (x * 12);
}


double chebyshev_eval(double x, const double *a, const int n)
{
  double b0, b1, b2, twox;
  int i;

  if (n < 1 || n > 1000) ML_ERR_return_NAN;
  if (x < -1.1 || x > 1.1) ML_ERR_return_NAN;
  twox = x * 2;
  b2 = b1 = 0;
  b0 = 0;
  for (i = 1; i <= n; i++) {
    b2 = b1;
    b1 = b0;
    b0 = twox * b1 - b2 + a[n - i];
  }
  return (b0 - b2) * 0.5;
}


/*
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 *  Merge in to R:
 *	Copyright (C) 2000, The R Core Development Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *
 *  DESCRIPTION
 *
 *    Computes the log of the error term in Stirling's formula.
 *      For n > 15, uses the series 1/12n - 1/360n^3 + ...
 *      For n <=15, integers or half-integers, uses stored values.
 *      For other n < 15, uses lgamma directly (don't use this to
 *        write lgamma!)
 *
 * Merge in to R:
 * Copyright (C) 2000, The R Core Development Team
 * R has lgammafn, and lgamma is not part of ISO C
 */

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
 *             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
 *             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
 *
 * see also lgammacor() in ./lgammacor.c  which computes almost the same!
 */

double stirlerr(double n)
{
  #define S0 0.083333333333333333333       /* 1/12 */
  #define S1 0.00277777777777777777778     /* 1/360 */
  #define S2 0.00079365079365079365079365  /* 1/1260 */
  #define S3 0.000595238095238095238095238 /* 1/1680 */
  #define S4 0.0008417508417508417508417508/* 1/1188 */

/*
  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/

  static const double sferr_halves[31] = {
    0.0, /* n=0 - wrong, place holder only */
    0.1534264097200273452913848,  /* 0.5 */
    0.0810614667953272582196702,  /* 1.0 */
    0.0548141210519176538961390,  /* 1.5 */
    0.0413406959554092940938221,  /* 2.0 */
    0.03316287351993628748511048, /* 2.5 */
    0.02767792568499833914878929, /* 3.0 */
    0.02374616365629749597132920, /* 3.5 */
    0.02079067210376509311152277, /* 4.0 */
    0.01848845053267318523077934, /* 4.5 */
    0.01664469118982119216319487, /* 5.0 */
    0.01513497322191737887351255, /* 5.5 */
    0.01387612882307074799874573, /* 6.0 */
    0.01281046524292022692424986, /* 6.5 */
    0.01189670994589177009505572, /* 7.0 */
    0.01110455975820691732662991, /* 7.5 */
    0.010411265261972096497478567, /* 8.0 */
    0.009799416126158803298389475, /* 8.5 */
    0.009255462182712732917728637, /* 9.0 */
    0.008768700134139385462952823, /* 9.5 */
    0.008330563433362871256469318, /* 10.0 */
    0.007934114564314020547248100, /* 10.5 */
    0.007573675487951840794972024, /* 11.0 */
    0.007244554301320383179543912, /* 11.5 */
    0.006942840107209529865664152, /* 12.0 */
    0.006665247032707682442354394, /* 12.5 */
    0.006408994188004207068439631, /* 13.0 */
    0.006171712263039457647532867, /* 13.5 */
    0.005951370112758847735624416, /* 14.0 */
    0.005746216513010115682023589, /* 14.5 */
    0.005554733551962801371038690  /* 15.0 */
  };
  double nn;

  if (n <= 15.0) {
    nn = n + n;
    if (nn == (int)nn) return(sferr_halves[(int)nn]);
    return(lgammafn(n + 1.) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI);
  }
  nn = n*n;
  if (n>500) return((S0-S1/nn)/n);
  if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
  if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
  return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}


/*
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    void dpsifn(double x, int n, int kode, int m,
 *		  double *ans, int *nz, int *ierr)
 *    double digamma(double x);
 *    double trigamma(double x)
 *    double tetragamma(double x)
 *    double pentagamma(double x)
 *
 *  DESCRIPTION
 *
 *    Compute the derivatives of the psi function
 *    and polygamma functions.
 *
 *    The following definitions are used in dpsifn:
 *
 *    Definition 1
 *
 *	 psi(x) = d/dx (ln(gamma(x)),  the first derivative of
 *				       the log gamma function.
 *
 *    Definition 2
 *		     k	 k
 *	 psi(k,x) = d /dx (psi(x)),    the k-th derivative
 *				       of psi(x).
 *
 *
 *    "dpsifn" computes a sequence of scaled derivatives of
 *    the psi function; i.e. for fixed x and m it computes
 *    the m-member sequence
 *
 *		  (-1)^(k+1) / gamma(k+1) * psi(k,x)
 *		     for k = n,...,n+m-1
 *
 *    where psi(k,x) is as defined above.   For kode=1, dpsifn
 *    returns the scaled derivatives as described.  kode=2 is
 *    operative only when k=0 and in that case dpsifn returns
 *    -psi(x) + ln(x).	That is, the logarithmic behavior for
 *    large x is removed when kode=2 and k=0.  When sums or
 *    differences of psi functions are computed the logarithmic
 *    terms can be combined analytically and computed separately
 *    to help retain significant digits.
 *
 *    Note that dpsifn(x, 0, 1, 1, ans) results in ans = -psi(x).
 *
 *  INPUT
 *
 *	x     - argument, x > 0.
 *
 *	n     - first member of the sequence, 0 <= n <= 100
 *		n == 0 gives ans(1) = -psi(x)	    for kode=1
 *				      -psi(x)+ln(x) for kode=2
 *
 *	kode  - selection parameter
 *		kode == 1 returns scaled derivatives of the
 *		psi function.
 *		kode == 2 returns scaled derivatives of the
 *		psi function except when n=0. In this case,
 *		ans(1) = -psi(x) + ln(x) is returned.
 *
 *	m     - number of members of the sequence, m >= 1
 *
 *  OUTPUT
 *
 *	ans   - a vector of length at least m whose first m
 *		components contain the sequence of derivatives
 *		scaled according to kode.
 *
 *	nz    - underflow flag
 *		nz == 0, a normal return
 *		nz != 0, underflow, last nz components of ans are
 *			 set to zero, ans(m-k+1)=0.0, k=1,...,nz
 *
 *	ierr  - error flag
 *		ierr=0, a normal return, computation completed
 *		ierr=1, input error,	 no computation
 *		ierr=2, overflow,	 x too small or n+m-1 too
 *			large or both
 *		ierr=3, error,		 n too large. dimensioned
 *			array trmr(nmax) is not large enough for n
 *
 *    The nominal computational accuracy is the maximum of unit
 *    roundoff (d1mach(4)) and 1e-18 since critical constants
 *    are given to only 18 digits.
 *
 *    The basic method of evaluation is the asymptotic expansion
 *    for large x >= xmin followed by backward recursion on a two
 *    term recursion relation
 *
 *	     w(x+1) + x^(-n-1) = w(x).
 *
 *    this is supplemented by a series
 *
 *	     sum( (x+k)^(-n-1) , k=0,1,2,... )
 *
 *    which converges rapidly for large n. both xmin and the
 *    number of terms of the series are calculated from the unit
 *    roundoff of the machine environment.
 *
 *  AUTHOR
 *
 *    Amos, D. E.  	(Fortran)
 *    Ross Ihaka   	(C Translation)
 *    Martin Maechler   (x < 0, and psigamma())
 *
 *  REFERENCES
 *
 *    Handbook of Mathematical Functions,
 *    National Bureau of Standards Applied Mathematics Series 55,
 *    Edited by M. Abramowitz and I. A. Stegun, equations 6.3.5,
 *    6.3.18, 6.4.6, 6.4.9 and 6.4.10, pp.258-260, 1964.
 *
 *    D. E. Amos, (1983). "A Portable Fortran Subroutine for
 *    Derivatives of the Psi Function", Algorithm 610,
 *    TOMS 9(4), pp. 494-502.
 *
 *    Routines called: Rf_d1mach, Rf_i1mach.
 */

/* From R, currently only used for kode = 1, m = 1, n in {0,1,2,3} : */
void dpsifn(double x, int n, int kode, int m, double *ans, int *nz, int *ierr)
{
    static const double bvalues[] = {	/* Bernoulli Numbers */
        1.00000000000000000e+00,
        -5.00000000000000000e-01,
        1.66666666666666667e-01,
        -3.33333333333333333e-02,
        2.38095238095238095e-02,
        -3.33333333333333333e-02,
        7.57575757575757576e-02,
        -2.53113553113553114e-01,
        1.16666666666666667e+00,
        -7.09215686274509804e+00,
        5.49711779448621554e+01,
        -5.29124242424242424e+02,
        6.19212318840579710e+03,
        -8.65802531135531136e+04,
        1.42551716666666667e+06,
        -2.72982310678160920e+07,
        6.01580873900642368e+08,
        -1.51163157670921569e+10,
        4.29614643061166667e+11,
        -1.37116552050883328e+13,
        4.88332318973593167e+14,
        -1.92965793419400681e+16
    };
    
    int i, j, k, mm, mx, nn, np, nx, fn;
    double arg, den, elim, eps, fln, fx, rln, rxsq,
	r1m4, r1m5, s, slope, t, ta, tk, tol, tols, tss, tst,
	tt, t1, t2, wdtol, xdmln, xdmy, xinc, xln = 0.0 /* -Wall */,
	xm, xmin, xq, yint;
    double trm[23], trmr[n_max + 1];
    
    *ierr = 0;
    if (n < 0 || kode < 1 || kode > 2 || m < 1) {
        *ierr = 1;
        return;
    }
    if (x <= 0.) {
        /* use	Abramowitz & Stegun 6.4.7 "Reflection Formula"
         *	psi(k, x) = (-1)^k psi(k, 1-x)	-  pi^{n+1} (d/dx)^n cot(x)
         */
        if (x == (long)x) {
            /* non-positive integer : +Inf or NaN depends on n */
            for (j=0; j < m; j++) /* k = j + n : */
                ans[j] = ((j+n) % 2) ? ML_POSINF : ML_NAN;
            return;
        }
        /* This could cancel badly */
        dpsifn(1. - x, n, /*kode = */ 1, m, ans, nz, ierr);
        /* ans[j] == (-1)^(k+1) / gamma(k+1) * psi(k, 1 - x)
         *	     for j = 0:(m-1) ,	k = n + j
         */
        
        /* Cheat for now: only work for	 m = 1, n in {0,1,2,3} : */
        if(m > 1 || n > 3) {/* doesn't happen for digamma() .. pentagamma() */
            /* not yet implemented */
            *ierr = 4; return;
        }
        x *= M_PI_R; /* pi * x */
        if (n == 0)
            tt = cos(x)/sin(x);
        else if (n == 1)
            tt = -1/pow(sin(x),2);
        else if (n == 2)
            tt = 2*cos(x)/pow(sin(x),3);
        else if (n == 3)
            tt = -2*(2*pow(cos(x),2) + 1)/pow(sin(x),4);
        else /* can not happen! */
            tt = ML_NAN;
        /* end cheat */
        
        s = (n % 2) ? -1. : 1.;/* s = (-1)^n */
        /* t := pi^(n+1) * d_n(x) / gamma(n+1)	, where
         *		   d_n(x) := (d/dx)^n cot(x)*/
        t1 = t2 = s = 1.;
        for (k=0, j=k-n; j < m; k++, j++, s = -s) {
            /* k == n+j , s = (-1)^k */
            t1 *= M_PI_R;/* t1 == pi^(k+1) */
            if(k >= 2)
                t2 *= k;/* t2 == k! == gamma(k+1) */
            if(j >= 0) /* by cheat above,  tt === d_k(x) */
                ans[j] = s*(ans[j] + t1/t2 * tt);
        }
        if (n == 0 && kode == 2) /* unused from R, but "wrong": xln === 0 :*/
            ans[0] += xln;
        return;
    } /* x <= 0 */
    
    /* else :  x > 0 */
    *nz = 0;
    xln = log(x);
    if(kode == 1 && m == 1) {/* the R case  ---  for very large x: */
        double lrg = 1/(2. * DBL_EPSILON);
        if(n == 0 && x * xln > lrg) {
            ans[0] = -xln;
            return;
        }
        else if(n >= 1 && x > n * lrg) {
            ans[0] = exp(-n * xln)/n; /* == x^-n / n  ==  1/(n * x^n) */
            return;
        }
    }
    mm = m;
    nx = int_min(-Rf_i1mach(15), Rf_i1mach(16));/* = 1021 */
    r1m5 = Rf_d1mach(5);
    r1m4 = Rf_d1mach(4) * 0.5;
    wdtol = dbl_max(r1m4, 0.5e-18); /* 1.11e-16 */
    
    /* elim = approximate exponential over and underflow limit */
    elim = 2.302 * (nx * r1m5 - 3.0);/* = 700.6174... */
    for (;;) {
        nn = n + mm - 1;
        fn = nn;
        t = (fn + 1) * xln;
        
        /* overflow and underflow test for small and large x */
        
        if (fabs(t) > elim) {
            if (t <= 0.0) {
                *nz = 0;
                *ierr = 2;
                return;
            }
        }
        else {
            if (x < wdtol) {
                ans[0] = pow(x, -n-1.0);
                if (mm != 1) {
                    for (k = 1; k < mm ; k++)
                        ans[k] = ans[k-1] / x;
                }
                if (n == 0 && kode == 2)
                    ans[0] += xln;
                return;
            }
            
            /* compute xmin and the number of terms of the series,  fln+1 */
            
            rln = r1m5 * Rf_i1mach(14);
            rln = dbl_min(rln, 18.06);
            fln = dbl_max(rln, 3.0) - 3.0;
            yint = 3.50 + 0.40 * fln;
            slope = 0.21 + fln * (0.0006038 * fln + 0.008677);
            xm = yint + slope * fn;
            mx = (int)xm + 1;
            xmin = mx;
            if (n != 0) {
                xm = -2.302 * rln - dbl_min(0.0, xln);
                arg = xm / n;
                arg = dbl_min(0.0, arg);
                eps = exp(arg);
                xm = 1.0 - eps;
                if (fabs(arg) < 1.0e-3)
                    xm = -arg;
                fln = x * xm / eps;
                xm = xmin - x;
                if (xm > 7.0 && fln < 15.0)
                    break;
            }
            xdmy = x;
            xdmln = xln;
            xinc = 0.0;
            if (x < xmin) {
                nx = (int)x;
                xinc = xmin - nx;
                xdmy = x + xinc;
                xdmln = log(xdmy);
            }
            
            /* generate w(n+mm-1, x) by the asymptotic expansion */
            
            t = fn * xdmln;
            t1 = xdmln + xdmln;
            t2 = t + xdmln;
            tk = dbl_max(fabs(t), dbl_max(fabs(t1), fabs(t2)));
            if (tk <= elim) /* for all but large x */
                goto L10;
        }
        nz++; /* underflow */
        mm--;
        ans[mm] = 0.;
        if (mm == 0)
            return;
    } /* end{for ()} */
    nn = (int)fln + 1;
    np = n + 1;
    t1 = (n + 1) * xln;
    t = exp(-t1);
    s = t;
    den = x;
    for (i=1; i <= nn; i++) {
        den += 1.;
        trm[i] = pow(den, (double)-np);
        s += trm[i];
    }
    ans[0] = s;
    if (n == 0 && kode == 2)
        ans[0] = s + xln;
    
    if (mm != 1) { /* generate higher derivatives, j > n */
        
        tol = wdtol / 5.0;
        for (j = 1; j < mm; j++) {
            t /= x;
            s = t;
            tols = t * tol;
            den = x;
            for (i=1; i <= nn; i++) {
                den += 1.;
                trm[i] /= den;
                s += trm[i];
                if (trm[i] < tols)
                    break;
            }
            ans[j] = s;
        }
    }
    return;
    
L10:
    tss = exp(-t);
    tt = 0.5 / xdmy;
    t1 = tt;
    tst = wdtol * tt;
    if (nn != 0)
        t1 = tt + 1.0 / fn;
    rxsq = 1.0 / (xdmy * xdmy);
    ta = 0.5 * rxsq;
    t = (fn + 1) * ta;
    s = t * bvalues[2];
    if (fabs(s) >= tst) {
        tk = 2.0;
        for (k = 4; k <= 22; k++) {
            t = t * ((tk + fn + 1)/(tk + 1.0))*((tk + fn)/(tk + 2.0)) * rxsq;
            trm[k] = t * bvalues[k-1];
            if (fabs(trm[k]) < tst)
                break;
            s += trm[k];
            tk += 2.;
        }
    }
    s = (s + t1) * tss;
    if (xinc != 0.0) {
        
        /* backward recur from xdmy to x */
        
        nx = (int)xinc;
        np = nn + 1;
        if (nx > n_max) {
            *nz = 0;
            *ierr = 3;
            return;
        }
        else {
            if (nn==0)
                goto L20;
            xm = xinc - 1.0;
            fx = x + xm;
            
            /* this loop should not be changed. fx is accurate when x is small */
            for (i = 1; i <= nx; i++) {
                trmr[i] = pow(fx, (double)-np);
                s += trmr[i];
                xm -= 1.;
                fx = x + xm;
            }
        }
    }
    ans[mm-1] = s;
    if (fn == 0)
        goto L30;
    
    /* generate lower derivatives,  j < n+mm-1 */
    
    for (j = 2; j <= mm; j++) {
        fn--;
        tss *= xdmy;
        t1 = tt;
        if (fn!=0)
            t1 = tt + 1.0 / fn;
        t = (fn + 1) * ta;
        s = t * bvalues[2];
        if (fabs(s) >= tst) {
            tk = 4 + fn;
            for (k=4; k <= 22; k++) {
                trm[k] = trm[k] * (fn + 1) / tk;
                if (fabs(trm[k]) < tst)
                    break;
                s += trm[k];
                tk += 2.;
            }
        }
        s = (s + t1) * tss;
        if (xinc != 0.0) {
            if (fn == 0)
                goto L20;
            xm = xinc - 1.0;
            fx = x + xm;
            for (i=1 ; i<=nx ; i++) {
                trmr[i] = trmr[i] * fx;
                s += trmr[i];
                xm -= 1.;
                fx = x + xm;
            }
        }
        ans[mm - j] = s;
        if (fn == 0)
            goto L30;
    }
    return;
    
L20:
    for (i = 1; i <= nx; i++)
        s += 1. / (x + (nx - i)); /* avoid disastrous cancellation, PR#13714 */
    
L30:
    if (kode != 2) /* always */
        ans[0] = s - xdmln;
    else if (xdmy != x) {
        xq = xdmy / x;
        ans[0] = s - log(xq);
    }
    return;
} /* dpsifn() */

double digamma(double x)
{
    double ans;
    int nz, ierr;
    if(ISNAN(x)) return x;
    dpsifn(x, 0, 1, 1, &ans, &nz, &ierr);
    ML_TREAT_psigam(ierr);
    return -ans;
}


/*
 *  Mathlib - A Mathematical Function Library
 *  Copyright (C) 1998  Ross Ihaka
 *  Copyright (C) 2000-7 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

int Rf_i1mach(int i)
{
    switch(i) {
            
        case  1: return 5;
        case  2: return 6;
        case  3: return 0;
        case  4: return 0;
            
        case  5: return CHAR_BIT * sizeof(int);
        case  6: return sizeof(int)/sizeof(char);
            
        case  7: return 2;
        case  8: return CHAR_BIT * sizeof(int) - 1;
        case  9: return INT_MAX;
            
        case 10: return FLT_RADIX;
            
        case 11: return FLT_MANT_DIG;
        case 12: return FLT_MIN_EXP;
        case 13: return FLT_MAX_EXP;
            
        case 14: return DBL_MANT_DIG;
        case 15: return DBL_MIN_EXP;
        case 16: return DBL_MAX_EXP;
            
        default: return 0;
    }
}

/* NaNs propagated correctly */

/*-- FIXME:  Eliminate calls to these
 *   =====   o   from C code when
 *	     o   it is only used to initialize "static" variables (threading)
 *  and use the DBL_... constants instead
 */

double Rf_d1mach(int i)
{
    switch(i) {
        case 1: return DBL_MIN;
        case 2: return DBL_MAX;
            
        case 3: /* = FLT_RADIX  ^ - DBL_MANT_DIG
                 for IEEE:  = 2^-53 = 1.110223e-16 = .5*DBL_EPSILON */
            return 0.5*DBL_EPSILON;
            
        case 4: /* = FLT_RADIX  ^ (1- DBL_MANT_DIG) =
                 for IEEE:  = 2^-52 = DBL_EPSILON */
            return DBL_EPSILON;
            
        case 5: return M_LOG10_2;
            
        default: return 0.0;
    }
}

double dbl_min(double x, double y)

{
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(y)) return x + y;
#endif
  return (x < y) ? x : y;
}

double dbl_max(double x, double y)

{
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(y)) return x + y;
#endif
  return (x < y) ? y : x;
}

int int_min(int x, int y)
{
    return (x < y) ? x : y;
}

int int_max(int x, int y)
{
  return (x > y) ? x : y;
}

double hypergeometric_1_F_1(double a,
                            double b,
                            double z)

{
	int j,crit;
	double a_j,s_j;
	
	a_j = 1;
	s_j = a_j;
	j = 0;
	crit = 0;
	do {
		a_j *= (a + j) / (b + j) * z / (double) (j + 1);
		++j;
		if((a_j / s_j) < PREC) ++crit; else crit = 0;
		s_j += a_j;
	}
	while (crit != 2 && j < 2000);
	if (j == 2000) {
#pragma omp critical (error)
    printf("\nProblem in computation of the hypergeometric function 1_F_1(a,b,z) with:\n");
    printf("\ta = %.3f\n",a);
    printf("\tb = %.3f\n",b);
    printf("\tz = %.3f\n",z);
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
  return(s_j);
}

double F_ST(int ***x,
            int **n,
            int nbr_loci,
            int nbr_demes)

{
  int i,j,k;
  int S_1,S_2;
  int empty_deme;
  double SSI,SSP,MSI,MSP,n_c,n_d;
  double num,den;
  double pi_hat[2];

  num = den = 0.0;
  for (j = 0; j < nbr_loci; ++j) {
    S_1 = S_2 = 0;
    empty_deme = 0;
    SSI = SSP = 0.0;
    for (k = 0; k < 2; ++k) {
      pi_hat[k] = 0.0;
    }
    for (i = 0; i < nbr_demes; ++i) {
      if (n[i][j] == 0) {
        empty_deme += 1;
        continue;
      }
      S_1 += n[i][j];
      S_2 += pow(n[i][j],2);
      for (k = 0; k < 2; ++k) {
        pi_hat[k] += (double) x[i][j][k];
        SSI += (double) x[i][j][k] - pow(x[i][j][k],2) / n[i][j];
      }
    }
    for (k = 0; k < 2; ++k) {
      pi_hat[k] /= (double) S_1;
    }
    for (i = 0; i < nbr_demes; ++i) {
      if (n[i][j] == 0) {
        continue;
      }
      for (k = 0; k < 2; ++k) {
        SSP += n[i][j] * pow(((double) x[i][j][k] / n[i][j] - pi_hat[k]),2);
      }
    }
    n_d = (nbr_demes - empty_deme);
    n_c = ((double) S_1 - (double) S_2 / S_1) / ((double) n_d - 1.0);
    MSI = SSI / ((double) S_1 - n_d);
    MSP = SSP / ((double) n_d - 1.0);
    num += (MSP - MSI);
    den += (MSP + (n_c - 1.0) * MSI);
  }
  if (den > 0) {
    return (num / den);
  }
  else {
    return ML_POSINF;
  }
}

double F_ST_poolseq(int ***x,
                    int **n,
                    int nbr_loci,
                    int nbr_demes,
                    int *pool_size)

{
  int i,j,k;
  int R_1,R_2;
  double C_1,C_1_star,n_c;
  double MSI,MSP,SSI,SSP;
  double num,den;
  double pi_hat[2];
  
  num = den = 0.0;
  for (j = 0; j < nbr_loci; ++j) {
    R_1 = R_2 = 0;
    C_1 = C_1_star = 0.0;
    SSI = SSP = 0.0;
    for (k = 0; k < 2; ++k) {
      pi_hat[k] = 0.0;
    }
    for (i = 0; i < nbr_demes; ++i) {
      C_1 += (double) n[i][j] / pool_size[i] + (double) (pool_size[i] - 1.0) / pool_size[i];
      C_1_star += (double) n[i][j] * ((double) n[i][j] / pool_size[i] + (double) (pool_size[i] - 1.0) / pool_size[i]);
      R_1 += n[i][j];
      R_2 += n[i][j] * n[i][j];
      if (n[i][j] == 0) {
        continue;
      }
      for (k = 0; k < 2; ++k) {
        pi_hat[k] += (double) x[i][j][k];
        SSI += (double) x[i][j][k] - pow(x[i][j][k],2) / n[i][j];
      }
    }
    for (k = 0; k < 2; ++k) {
      pi_hat[k] /= (double) R_1;
    }
    C_1_star /= (double) R_1;
    for (i = 0; i < nbr_demes; ++i) {
      if (n[i][j] == 0) {
        continue;
      }
      for (k = 0; k < 2; ++k) {
        SSP += n[i][j] * pow(((double) x[i][j][k] / n[i][j] - pi_hat[k]),2);
      }
    }
    n_c = ((double) R_1 - (double) R_2 / R_1) / ((double) C_1 - C_1_star);
    MSI = SSI / ((double) R_1 - C_1);
    MSP = SSP / ((double) C_1 - C_1_star);
    num += (MSP - MSI);
    den += (MSP + (n_c - 1.0) * MSI);
  }
  if (den > 0) {
    return (num / den);
  }
  else {
    return ML_POSINF;
  }
}

double log_exponential_prior(double x,                                                         // Accounts for an unbounded exponential prior
                             double a)

{
	double value;
	
	value = - x / a - log(a);
	return(value);
}

double log_uniform_prior(double x)

{
	double value;
	
	value = -log(log(max_M / min_M)) - log(x);
	return(value);
}

double log_beta_prior(double *x,
                      double a,
                      double b)

{
  double value;

  value = (a - 1.0) * log(x[0]) + (b - 1.0) * log(x[1]);
  return(value);
}

double log_inverse_gamma_prior(double x,
                               double alpha,
                               double beta)

{
 	double value;
	
	value = alpha * log(beta) - (alpha + 1) * log(x) - beta / x - lgammafn(alpha);
	return(value);
}
