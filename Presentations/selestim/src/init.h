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

#ifndef _INIT_H
#define _INIT_H

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "dcmt0.6.1b/include/dc.h"
#include "defs.h"
#include "mymath.h"
#include "rand.h"

void init_parameters(markov_state_struct *current,data_struct data);
void init_moves(data_struct data,updates_struct *move);
void init_moments(data_struct data,moments_struct *postmean,moments_struct *postmsqr);
void init_updates(data_struct data,updates_struct *update,updates_struct *accept);

#endif

