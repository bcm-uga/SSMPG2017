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

#ifndef _RAND_H
#define _RAND_H

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "dcmt0.6.1b/include/dc.h"
#include "defs.h"
#include "mymath.h"

int mt_count;
mt_struct **mtss;

double genrand_real1(mt_struct *seed);                  // generates a random number on [0,1]-real-interval
double genrand_real2(mt_struct *seed);                  // generates a random number on [0,1)-real-interval
double genrand_real3(mt_struct *seed);                  // generates a random number on (0,1)-real-interval
double exp_rand(mt_struct *seed);
double rbeta(double aa, double bb,mt_struct *seed);
double rgamma(double a, double scale,mt_struct *seed);
double rbinom(double nin, double pp,mt_struct *seed);
double rnorm(mt_struct *seed);
double rexp(double scale,mt_struct *seed);

#endif
