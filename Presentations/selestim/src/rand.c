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

#include "rand.h"

double genrand_real1(mt_struct *seed)                                           // generates a random number on [0,1]-real-interval
{
  return genrand_mt(seed)*(1.0/4294967295.0);                                   // divided by 2^32-1
}

double genrand_real2(mt_struct *seed)                                           // generates a random number on [0,1)-real-interval
{
  return genrand_mt(seed)*(1.0/4294967296.0);                                   // divided by 2^32
}

double genrand_real3(mt_struct *seed)                                           // generates a random number on (0,1)-real-interval
{
  return (((double)genrand_mt(seed)) + 0.5)*(1.0/4294967296.0);                 // divided by 2^32
}

double exp_rand(mt_struct *seed)
{
  /* q[k-1] = sum(log(2)^k / k!)  k=1,..,n, */
  /* The highest n (here 8) is determined by q[n-1] = 1.0 */
  /* within standard precision */
  static const double q[] =
  {
    0.6931471805599453,
    0.9333736875190459,
    0.9888777961838675,
    0.9984959252914960,
    0.9998292811061389,
    0.9999833164100727,
    0.9999985691438767,
    0.9999998906925558,
    0.9999999924734159,
    0.9999999995283275,
    0.9999999999728814,
    0.9999999999985598,
    0.9999999999999289,
    0.9999999999999968,
    0.9999999999999999,
    1.0000000000000000
  };
  double a, u, ustar, umin;
  int i;
  
  a = 0.;
  /* precaution if u = 0 is ever returned */
  u = genrand_real1(seed);
  while(u <= 0.0 || u >= 1.0) u = genrand_real1(seed);
  for (;;) {
    u += u;
    if (u > 1.0)
	    break;
    a += q[0];
  }
  u -= 1.;
  
  if (u <= q[0])
    return a + u;
  
  i = 0;
  ustar = genrand_real1(seed);
  umin = ustar;
  do {
    ustar = genrand_real1(seed);
    if (ustar < umin)
	    umin = ustar;
    i++;
  } while (u > q[i]);
  return a + umin * q[0];
}

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 2000 The R Development Core Team
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
 */

/* Reference:
 * R. C. H. Cheng (1978).
 * Generating beta variates with nonintegral shape parameters.
 * Communications of the ACM 21, 317-322.
 * (Algorithms BB and BC)
 */

double rbeta(double aa,double bb,mt_struct *seed)

{
  double a, b, alpha;
  double r, s, t, u1, u2, v, w, y, z;
  int qsame;
  /* FIXME:  Keep Globals (properly) for threading */
  /* Uses these GLOBALS to save time when many rv's are generated : */
  static double beta, gamma, delta, k1, k2;
  static double olda = -1.0;
  static double oldb = -1.0;
  
  if (aa <= 0. || bb <= 0. || (!R_FINITE(aa) && !R_FINITE(bb))) ML_ERR_return_NAN;
  if (!R_FINITE(aa)) return 1.0;
  if (!R_FINITE(bb)) return 0.0;
  /* Test if we need new "initializing" */
  qsame = (olda == aa) && (oldb == bb);
  if (!qsame) { olda = aa; oldb = bb; }
  a = dbl_min(aa, bb);
  b = dbl_max(aa, bb); /* a <= b */
  alpha = a + b;
  if (a <= 1.0) {	/* --- Algorithm BC --- */
    /* changed notation, now also a <= b (was reversed) */
    if (!qsame) { /* initialize */
      beta = 1.0 / a;
      delta = 1.0 + b - a;
      k1 = delta * (0.0138889 + 0.0416667 * a) / (b * beta - 0.777778);
      k2 = 0.25 + (0.5 + 0.25 / delta) * a;
    }
    /* FIXME: "do { } while()", but not trivially because of "continue"s:*/
    for (;;) {
      u1 = genrand_real2(seed);
      u2 = genrand_real2(seed);
      if (u1 < 0.5) {
        y = u1 * u2;
        z = u1 * y;
        if (0.25 * u2 + z - y >= k1) continue;
      }
      else {
        z = u1 * u1 * u2;
        if (z <= 0.25) {
          v_w_from__u1_bet(b);
          break;
        }
        if (z >= k2) continue;
      }
      v_w_from__u1_bet(b);
      if (alpha * (log(alpha / (a + w)) + v) - 1.3862944 >= log(z)) break;
    }
    return (aa == a) ? a / (a + w) : w / (a + w);
  }
  else {		/* Algorithm BB */
    if (!qsame) { /* initialize */
      beta = sqrt((alpha - 2.0) / (2.0 * a * b - alpha));
      gamma = a + 1.0 / beta;
    }
    do {
      u1 = genrand_real2(seed);
      u2 = genrand_real2(seed);
      v_w_from__u1_bet(a);
      z = u1 * u1 * u2;
      r = gamma * v - 1.3862944;
      s = a + r - w;
      if (s + 2.609438 >= 5.0 * z) break;
      t = log(z);
      if (s > t) break;
    }
    while (r + alpha * log(alpha / (b + w)) < t);
    return (aa != a) ? b / (b + w) : w / (b + w);
  }
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2002 The R Development Core Team
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
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *  SYNOPSIS
 *
 *	#include <Rmath.h>
 *	double rbinom(double nin, double pp)
 *
 *  DESCRIPTION
 *
 *	Random variates from the binomial distribution.
 *
 *  REFERENCE
 *
 *	Kachitvichyanukul, V. and Schmeiser, B. W. (1988).
 *	Binomial random variate generation.
 *	Communications of the ACM 31, 216-222.
 *	(Algorithm BTPEC).
 */

double rbinom(double nin, double pp,mt_struct *seed)
{
  /* FIXME: These should become THREAD_specific globals : */
  
  static double c, fm, npq, p1, p2, p3, p4, qn;
  static double xl, xll, xlr, xm, xr;
  
  static double psave = -1.0;
  static int nsave = -1;
  static int m;
  
  double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
  double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
  int i,ix,k, n;
  
  if (!R_FINITE(nin)) ML_ERR_return_NAN;
  n = floor(nin + 0.5);
  if (n != nin) ML_ERR_return_NAN;
  
  if (!R_FINITE(pp) ||
      /* n=0, p=0, p=1 are not errors <TSL>*/
      n < 0 || pp < 0. || pp > 1.)	ML_ERR_return_NAN;
  
  if (n == 0 || pp == 0.) return 0;
  if (pp == 1.) return n;
  
  p = dbl_min(pp, 1. - pp);
  q = 1. - p;
  np = n * p;
  r = p / q;
  g = r * (n + 1);
  
  /* Setup, perform only when parameters change [using static (globals): */
  
  /* FIXING: Want this thread safe
   -- use as little (thread globals) as possible
   */
  if (pp != psave || n != nsave) {
    psave = pp;
    nsave = n;
    if (np < 30.0) {
	    /* inverse cdf logic for mean less than 30 */
	    qn = pow(q, (double) n);
	    goto L_np_small;
    } else {
	    ffm = np + p;
	    m = ffm;
	    fm = m;
	    npq = np * q;
	    p1 = (int)(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
	    xm = fm + 0.5;
	    xl = xm - p1;
	    xr = xm + p1;
	    c = 0.134 + 20.5 / (15.3 + fm);
	    al = (ffm - xl) / (ffm - xl * p);
	    xll = al * (1.0 + 0.5 * al);
	    al = (xr - ffm) / (xr * q);
	    xlr = al * (1.0 + 0.5 * al);
	    p2 = p1 * (1.0 + c + c);
	    p3 = p2 + c / xll;
	    p4 = p3 + c / xlr;
    }
  } else if (n == nsave) {
    if (np < 30.0)
	    goto L_np_small;
  }
  
  /*-------------------------- np = n*p >= 30 : ------------------- */
  repeat {
    u = genrand_real2(seed) * p4;
    v = genrand_real2(seed);
    /* triangular region */
    if (u <= p1) {
      ix = xm - p1 * v + u;
      goto finis;
    }
    /* parallelogram region */
    if (u <= p2) {
      x = xl + (u - p1) / c;
      v = v * c + 1.0 - fabs(xm - x) / p1;
      if (v > 1.0 || v <= 0.)
	      continue;
      ix = x;
    } else {
      if (u > p3) {	/* right tail */
	      ix = xr - log(v) / xlr;
	      if (ix > n)
          continue;
	      v = v * (u - p3) * xlr;
      } else {/* left tail */
	      ix = xl + log(v) / xll;
	      if (ix < 0)
          continue;
	      v = v * (u - p2) * xll;
      }
    }
    /* determine appropriate way to perform accept/reject test */
    k = abs(ix - m);
    if (k <= 20 || k >= npq / 2 - 1) {
      /* explicit evaluation */
      f = 1.0;
      if (m < ix) {
	      for (i = m + 1; i <= ix; i++)
          f *= (g / i - r);
      } else if (m != ix) {
	      for (i = ix + 1; i <= m; i++)
          f /= (g / i - r);
      }
      if (v <= f)
	      goto finis;
    } else {
      /* squeezing using upper and lower bounds on log(f(x)) */
      amaxp = (k / npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
      ynorm = -k * k / (2.0 * npq);
      alv = log(v);
      if (alv < ynorm - amaxp)
	      goto finis;
      if (alv <= ynorm + amaxp) {
	      /* stirling's formula to machine accuracy */
	      /* for the final acceptance/rejection test */
	      x1 = ix + 1;
	      f1 = fm + 1.0;
	      z = n + 1 - fm;
	      w = n - ix + 1.0;
	      z2 = z * z;
	      x2 = x1 * x1;
	      f2 = f1 * f1;
	      w2 = w * w;
	      if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.)
          goto finis;
      }
    }
  }
  
L_np_small:
  /*---------------------- np = n*p < 30 : ------------------------- */
  
  repeat {
    ix = 0;
    f = qn;
    u = genrand_real2(seed);
    repeat {
      if (u < f)
        goto finis;
      if (ix > 110)
        break;
      u -= f;
      ix++;
      f *= (g / ix - r);
    }
  }
finis:
  if (psave > 0.5)
    ix = n - ix;
  return (double)ix;
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000--2008 The R Core Team
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
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double rgamma(double a, double scale);
 *
 *  DESCRIPTION
 *
 *    Random variates from the gamma distribution.
 *
 *  REFERENCES
 *
 *    [1] Shape parameter a >= 1.  Algorithm GD in:
 *
 *	  Ahrens, J.H. and Dieter, U. (1982).
 *	  Generating gamma variates by a modified
 *	  rejection technique.
 *	  Comm. ACM, 25, 47-54.
 *
 *
 *    [2] Shape parameter 0 < a < 1. Algorithm GS in:
 *
 *	  Ahrens, J.H. and Dieter, U. (1974).
 *	  Computer methods for sampling from gamma, beta,
 *	  poisson and binomial distributions.
 *	  Computing, 12, 223-246.
 *
 *    Input: a = parameter (mean) of the standard gamma distribution.
 *    Output: a variate from the gamma(a)-distribution
 */

double rgamma(double a, double scale,mt_struct *seed)
{
  /* Constants : */
  static const double sqrt32 = 5.656854;
  static const double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */
  
  /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
   * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
   * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
   */
  static const double q1 = 0.04166669;
  static const double q2 = 0.02083148;
  static const double q3 = 0.00801191;
  static const double q4 = 0.00144121;
  static const double q5 = -7.388e-5;
  static const double q6 = 2.4511e-4;
  static const double q7 = 2.424e-4;
  
  static const double a1 = 0.3333333;
  static const double a2 = -0.250003;
  static const double a3 = 0.2000062;
  static const double a4 = -0.1662921;
  static const double a5 = 0.1423657;
  static const double a6 = -0.1367177;
  static const double a7 = 0.1233795;
  
  /* State variables [FIXME for threading!] :*/
  static double aa = 0.;
  static double aaa = 0.;
  static double s, s2, d;    /* no. 1 (step 1) */
  static double q0, b, si, c;/* no. 2 (step 4) */
  
  double e, p, q, r, t, u, v, w, x, ret_val;
  
  if (!R_FINITE(a) || !R_FINITE(scale) || a < 0.0 || scale <= 0.0) {
    if(scale == 0.) return 0.;
    ML_ERR_return_NAN;
  }
  
  if (a < 1.) { /* GS algorithm for parameters a < 1 */
    if(a == 0)
	    return 0.;
    e = 1.0 + exp_m1 * a;
    repeat {
	    p = e * genrand_real2(seed);
	    if (p >= 1.0) {
        x = -log((e - p) / a);
        if (exp_rand(seed) >= (1.0 - a) * log(x))
          break;
	    } else {
        x = exp(log(p) / a);
        if (exp_rand(seed) >= x)
          break;
	    }
    }
    return scale * x;
  }
  
  /* --- a >= 1 : GD algorithm --- */
  
  /* Step 1: Recalculations of s2, s, d if a has changed */
  if (a != aa) {
    aa = a;
    s2 = a - 0.5;
    s = sqrt(s2);
    d = sqrt32 - s * 12.0;
  }
  /* Step 2: t = standard normal deviate,
   x = (s,1/2) -normal deviate. */
  
  /* immediate acceptance (i) */
  t = rnorm(seed);
  x = s + 0.5 * t;
  ret_val = x * x;
  if (t >= 0.0)
    return scale * ret_val;
  
  /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
  u = genrand_real2(seed);
  if (d * u <= t * t * t)
    return scale * ret_val;
  
  /* Step 4: recalculations of q0, b, si, c if necessary */
  
  if (a != aaa) {
    aaa = a;
    r = 1.0 / a;
    q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
           + q2) * r + q1) * r;
    
    /* Approximation depending on size of parameter a */
    /* The constants in the expressions for b, si and c */
    /* were established by numerical experiments */
    
    if (a <= 3.686) {
	    b = 0.463 + s + 0.178 * s2;
	    si = 1.235;
	    c = 0.195 / s - 0.079 + 0.16 * s;
    } else if (a <= 13.022) {
	    b = 1.654 + 0.0076 * s2;
	    si = 1.68 / s + 0.275;
	    c = 0.062 / s + 0.024;
    } else {
	    b = 1.77;
	    si = 0.75;
	    c = 0.1515 / s;
    }
  }
  /* Step 5: no quotient test if x not positive */
  
  if (x > 0.0) {
    /* Step 6: calculation of v and quotient q */
    v = t / (s + s);
    if (fabs(v) <= 0.25)
	    q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
                                + a3) * v + a2) * v + a1) * v;
    else
	    q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
    
    
    /* Step 7: quotient acceptance (q) */
    if (log(1.0 - u) <= q)
	    return scale * ret_val;
  }
  
  repeat {
    /* Step 8: e = standard exponential deviate
     *	u =  0,1 -uniform deviate
     *	t = (b,si)-double exponential (laplace) sample */
    e = exp_rand(seed);
    u = genrand_real2(seed);
    u = u + u - 1.0;
    if (u < 0.0)
	    t = b - si * e;
    else
	    t = b + si * e;
    /* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
    if (t >= -0.71874483771719) {
	    /* Step 10:	 calculation of v and quotient q */
	    v = t / (s + s);
	    if (fabs(v) <= 0.25)
        q = q0 + 0.5 * t * t *
		    ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
		      + a2) * v + a1) * v;
	    else
        q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
	    /* Step 11:	 hat acceptance (h) */
	    /* (if q not positive go to step 8) */
	    if (q > 0.0) {
        w = expm1(q);
        /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
        /* if t is rejected sample again at step 8 */
        if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
          break;
	    }
    }
  } /* repeat .. until  `t' is accepted */
  x = s + 0.5 * t;
  return scale * x * x;
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000--2008 The R Core Team
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
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double rexp(double scale)
 *
 *  DESCRIPTION
 *
 *    Random variates from the exponential distribution.
 *
 */

double rexp(double scale,mt_struct *seed)
{
  if (!R_FINITE(scale) || scale <= 0.0) {
    if(scale == 0.) return 0.;
    /* else */
    ML_ERR_return_NAN;
  }
  return scale * exp_rand(seed);
}

double rnorm(mt_struct *seed)

{
	double x,y,z,r2;
	
	do {
		x = 2.0 * genrand_real3(seed) - 1.0;                                             // Choose x,y in uniform square (-1,-1) to (+1,+1) ! CORRECTION 09-12-2004
		y = 2.0 * genrand_real3(seed) - 1.0;                                             //                                                 ! CORRECTION 09-12-2004
		r2 = x * x + y * y;                                                          // See if it is in the unit circle
	}
	while (r2 >= 1.0 || r2 == 0);                                                  // Box-Muller transform                            ! CORRECTION 09-12-2004
	z = y * sqrt(-2.0 * log(r2) / r2);
	return(z);
}
