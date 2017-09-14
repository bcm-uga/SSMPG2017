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

#ifndef _MOVES_H
#define _MOVES_H

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "dcmt0.6.1b/include/dc.h"
#include "defs.h"
#include "mymath.h"
#include "rand.h"

void update_counts_Gautier_unpublished(markov_state_struct *current,markov_state_struct *proposal,updates_struct move,updates_struct *update,updates_struct *accept,data_struct data);
void update_counts_Vitalis_unpublished(markov_state_struct *current,markov_state_struct *proposal,updates_struct move,updates_struct *update,updates_struct *accept,data_struct data);
void update_p(markov_state_struct *current,markov_state_struct *proposal,updates_struct move,updates_struct *update,updates_struct *accept,data_struct data);
void update_M_Gautier_unpublished(markov_state_struct *current,markov_state_struct *proposal,updates_struct move,updates_struct *update,updates_struct *accept,data_struct data);
void update_M_Vitalis_et_al_2014(markov_state_struct *current,markov_state_struct *proposal,updates_struct move,updates_struct *update,updates_struct *accept,data_struct data);
void update_pi(markov_state_struct *current,markov_state_struct *proposal,updates_struct move,updates_struct *update,updates_struct *accept,data_struct data);
void update_lambda_gibbs(markov_state_struct *current,data_struct data);
void update_lambda_metropolis(markov_state_struct *current,markov_state_struct *proposal,updates_struct move,updates_struct *update,updates_struct *accept,data_struct data);
void update_sigma(markov_state_struct *current,markov_state_struct *proposal,updates_struct move,updates_struct *update,updates_struct *accept,data_struct data);
void update_delta(markov_state_struct *current,markov_state_struct *proposal,updates_struct move,updates_struct *update,updates_struct *accept,data_struct data);
void update_kappa(markov_state_struct *current,markov_state_struct *proposal,data_struct data);
void update_beta_parameters_Vitalis_unpublished(markov_state_struct *current,updates_struct move,updates_struct *update,updates_struct *accept,data_struct data);
void update_beta_parameters_Gautier_unpublished(markov_state_struct *current,updates_struct move,updates_struct *update,updates_struct *accept,data_struct data);

#endif

