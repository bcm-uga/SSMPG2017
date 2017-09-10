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

#include "moves.h"

void update_counts_Gautier_unpublished(markov_state_struct *current,
                                       markov_state_struct *proposal,
                                       updates_struct move,
                                       updates_struct *update,
                                       updates_struct *accept,
                                       data_struct data)

{
	int test,tmp;
	int i,j,k,ij;
  int lower_bound,upper_bound,forward,reverse;
	double frequency,ratio;
  int ***counts;
  double ***p;
  int **n;
  
  counts = proposal -> value.counts;
  n = data.total_nbr_counts;
  p = current -> value.p;

#pragma omp parallel for schedule(guided) private(frequency,lower_bound,upper_bound,forward,reverse,tmp,ratio,test,k,i,j)

  for (ij = 0; ij < (data.nbr_demes * data.nbr_loci); ++ij) {
    i = ij / data.nbr_loci;
    j = ij % data.nbr_loci;
    ++update -> counts[i][j];
    ratio = 0.0;
    k = (int) (genrand_real2(mtss[omp_get_thread_num()]) * 2);
    frequency = (double) current -> value.counts[i][j][k] / n[i][j];
    if ((current -> value.counts[i][j][k] != 0) && (current -> value.counts[i][j][1 - k] != n[i][j])) {
      ratio -= lgammafn(n[i][j] + 1) - lgammafn(current -> value.counts[i][j][k] + 1) - lgammafn(current -> value.counts[i][j][1 - k] + 1) + data.reads[i][j][k] * log(frequency) + data.reads[i][j][1 - k] * log(1.0 - frequency);
    }
    lower_bound = int_max(min_counts[i][j][k],current -> value.counts[i][j][k] - (int) move.counts[i][j]);
    upper_bound = int_min(max_counts[i][j][k],current -> value.counts[i][j][k] + (int) move.counts[i][j]);
    do {
      tmp = lower_bound + (int) (genrand_real2(mtss[omp_get_thread_num()]) * (upper_bound - lower_bound + 1));
    }
    while (tmp == current -> value.counts[i][j][k]);
    counts[i][j][k] = tmp;
    counts[i][j][1 - k] = n[i][j] - tmp;
    forward = upper_bound - lower_bound;
    reverse = int_min(max_counts[i][j][k],counts[i][j][k] + (int) move.counts[i][j]) - int_max(min_counts[i][j][k],counts[i][j][k] - (int) move.counts[i][j]);
    ratio += log(forward) - log(reverse);
    frequency = (double) counts[i][j][k] / n[i][j];
    if ((counts[i][j][k] != 0) && (counts[i][j][1 - k] != n[i][j])) {
      ratio += lgammafn(n[i][j] + 1) - lgammafn(counts[i][j][k] + 1) - lgammafn(counts[i][j][1 - k] + 1) + data.reads[i][j][k] * log(frequency) + data.reads[i][j][1 - k] * log(1.0 - frequency);
    }
    ratio += (counts[i][j][k] - current -> value.counts[i][j][k]) * (log(p[i][j][k]) - log(p[i][j][1 - k]));
    test = (ratio >= 0) ? 1 : genrand_real1(mtss[omp_get_thread_num()]) < exp(ratio);
    if (test) {
      ++accept -> counts[i][j];
      for (k = 0; k < 2; ++k) {
        current -> value.counts[i][j][k] = counts[i][j][k];
      }
    }
  }
}

void update_counts_Vitalis_unpublished(markov_state_struct *current,
                                       markov_state_struct *proposal,
                                       updates_struct move,
                                       updates_struct *update,
                                       updates_struct *accept,
                                       data_struct data)

{
	int test;
	int i,j,k,ij;
  double current_freq,proposal_freq,tmp;
	double ratio;
  int ***counts;
  double ***p;
  int **n;
  
  counts = proposal -> value.counts;
  n = data.total_nbr_counts;
  p = current -> value.p;

#pragma omp parallel for schedule(guided) private(current_freq,proposal_freq,tmp,ratio,test,k,i,j)

  for (ij = 0; ij < (data.nbr_demes * data.nbr_loci); ++ij) {
    i = ij / data.nbr_loci;
    j = ij % data.nbr_loci;
    if (n[i][j] > 0) {
      ++update -> counts[i][j];
      k = (int) (genrand_real2(mtss[omp_get_thread_num()]) * 2);
      tmp = (double) current -> value.counts[i][j][k] + genrand_real2(mtss[omp_get_thread_num()]) * (2.0 * move.counts[i][j] + 1.0) - move.counts[i][j];
      if ((tmp < (max_counts[i][j][k] + 1)) && (tmp >= min_counts[i][j][k])) {
        counts[i][j][k] = (int) tmp;
      }
      else if (tmp < min_counts[i][j][k]) {
        counts[i][j][k] = (int) (2 * min_counts[i][j][k] - tmp);
      }
      else if (tmp > max_counts[i][j][k]) {
        counts[i][j][k] = ceil(2 * max_counts[i][j][k] - tmp + 1);
      }
      counts[i][j][1 - k] = n[i][j] - counts[i][j][k];
      ratio = 0.0;
      current_freq = (double) current -> value.counts[i][j][k] / n[i][j];
      if ((current -> value.counts[i][j][k] > 0) && (current -> value.counts[i][j][1 - k] > 0)) {
        ratio -= lgammafn(n[i][j] + 1) - lgammafn(current -> value.counts[i][j][k] + 1) - lgammafn(current -> value.counts[i][j][1 - k] + 1) + data.reads[i][j][k] * log(current_freq) + data.reads[i][j][1 - k] * log(1.0 - current_freq);
      }
      proposal_freq = (double) counts[i][j][k] / n[i][j];
      if ((counts[i][j][k] > 0) && (counts[i][j][1 - k] > 0)) {
        ratio += lgammafn(n[i][j] + 1) - lgammafn(counts[i][j][k] + 1) - lgammafn(counts[i][j][1 - k] + 1) + data.reads[i][j][k] * log(proposal_freq) + data.reads[i][j][1 - k] * log(1.0 - proposal_freq);
      }
      ratio += (counts[i][j][k] - current -> value.counts[i][j][k]) * (log(p[i][j][k]) - log(p[i][j][1 - k]));
      test = (ratio >= 0) ? 1 : genrand_real1(mtss[omp_get_thread_num()]) < exp(ratio);
      if (test) {
        ++accept -> counts[i][j];
        for (k = 0; k < 2; ++k) {
          current -> value.counts[i][j][k] = counts[i][j][k];
        }
      }
    }
  }
}

void update_p(markov_state_struct *current,
              markov_state_struct *proposal,
              updates_struct move,
              updates_struct *update,
              updates_struct *accept,
              data_struct data)

{
	int test;
	int i,j,k,ij;
	double tmp;
	double ratio;
	int ***x;
	double ***p;
	double *M;
	double **sigma;
	double **pi;
	int **kappa;
  
  p = proposal -> value.p;
  M = current -> value.M;
  pi = current -> value.pi;
  sigma = current -> value.sigma;
  kappa = current -> value.kappa;
  if (pooled_data) {
    x = current -> value.counts;
  }
  else {
    x = data.counts;
  }

#pragma omp parallel for schedule(guided) private(tmp,ratio,test,k,i,j)

  for (ij = 0; ij < (data.nbr_demes * data.nbr_loci); ++ij) {
    i = ij / data.nbr_loci;
    j = ij % data.nbr_loci;
    ++update -> p[i][j];
    k = (int) (genrand_real2(mtss[omp_get_thread_num()]) * 2);
    tmp = current -> value.p[i][j][k] + genrand_real2(mtss[omp_get_thread_num()]) * (2.0 * move.p[i][j]) - move.p[i][j];   // Reflecting boundaries (left side)
    if (tmp < MIN_P) tmp = 2.0 * MIN_P - tmp;
    if (tmp > (1.0 - MIN_P)) tmp = 2.0 * (1.0 - MIN_P) - tmp;                                                  // Reflecting boundaries (right side)
    proposal -> value.p[i][j][k] = tmp;
    proposal -> value.p[i][j][(1 - k)] = 1.0 - tmp;
    ratio = sigma[i][j] * (p[i][j][kappa[i][j]] - current -> value.p[i][j][kappa[i][j]]);
    for (k = 0; k < 2; ++k) {
      ratio += (x[i][j][k] + M[i] * pi[j][k] - 1.0) * (log(p[i][j][k]) - log(current -> value.p[i][j][k]));
    }
    test = (ratio >= 0) ? 1 : genrand_real1(mtss[omp_get_thread_num()]) < exp(ratio);							// Are the new parameter values accepted?
    if (test) {
      ++accept -> p[i][j];
      for (k = 0; k < 2; ++k) {
        current -> value.p[i][j][k] = proposal -> value.p[i][j][k];
      }
    }
  }
}

void update_pi(markov_state_struct *current,
               markov_state_struct *proposal,
               updates_struct move,
               updates_struct *update,
               updates_struct *accept,
               data_struct data)

{
	int test;
	int i,j,k;
	double tmp;
	double ratio;
	double ***p;
	double *M;
	double **pi;
	double **sigma;
	int **kappa;
  
  p = current -> value.p;
  M = current -> value.M;
  pi = proposal -> value.pi;
  sigma = current -> value.sigma;
  kappa = current -> value.kappa;
#pragma omp parallel for schedule(guided) private(tmp,ratio,test,i,k)
  for (j = 0; j < data.nbr_loci; ++j) {
    ++update -> pi[j];
    k = (int) (genrand_real2(mtss[omp_get_thread_num()]) * 2);
    tmp = current -> value.pi[j][k] + genrand_real2(mtss[omp_get_thread_num()]) * (2.0 * move.pi[j]) - move.pi[j];   // Reflecting boundaries (left side)
    if (tmp < MIN_PI) tmp = 2.0 * MIN_PI - tmp;
    if (tmp > (1.0 - MIN_PI)) tmp = 2.0 * (1.0 - MIN_PI) - tmp;                                                  // Reflecting boundaries (right side)
    pi[j][k] = tmp;
    pi[j][(1 - k)] = 1.0 - tmp;
    ratio = 0.0;
    for (i = 0; i < data.nbr_demes; ++i) {
      for (k = 0; k < 2; ++k) {
        ratio += M[i] * (pi[j][k] - current -> value.pi[j][k]) * log(p[i][j][k]);
      }
      proposal -> psi.log_1F1[i][j] = log(hypergeometric_1_F_1(M[i] * pi[j][kappa[i][j]],M[i],sigma[i][j]));
      proposal -> psi.Gamma_M_pi[i][j] = lgammafn(M[i] * pi[j][0]) + lgammafn(M[i] * pi[j][1]);
      ratio += (current -> psi.log_1F1[i][j] + current -> psi.Gamma_M_pi[i][j]);
      ratio -= (proposal -> psi.log_1F1[i][j] + proposal -> psi.Gamma_M_pi[i][j]);
    }
    proposal -> prior.pi[j][0] = log_beta_prior(pi[j],current -> value.beta_a,current -> value.beta_b);
    ratio += proposal -> prior.pi[j][0] - current -> prior.pi[j][0];
    test = (ratio >= 0) ? 1 : genrand_real1(mtss[omp_get_thread_num()]) < exp(ratio);							// Are the new parameter values accepted?
    if (test) {
      ++accept -> pi[j];
      for (k = 0; k < 2; ++k) {
        current -> value.pi[j][k] = pi[j][k];
      }
      for (i = 0; i < data.nbr_demes; ++i) {
        current -> psi.log_1F1[i][j] = proposal -> psi.log_1F1[i][j];
        current -> psi.Gamma_M_pi[i][j] = proposal -> psi.Gamma_M_pi[i][j];
      }
      current -> prior.pi[j][0] = proposal -> prior.pi[j][0];
    }
  }
}

void update_M_Gautier_unpublished(markov_state_struct *current,
                                  markov_state_struct *proposal,
                                  updates_struct move,
                                  updates_struct *update,
                                  updates_struct *accept,
                                  data_struct data)

{
	int i,j,k;
	int test;
  double current_F_ST;
  double proposal_F_ST;
  double tmp;
	double ratio;
  double LIMIT = 1e-3;
	double ***p;
	double *M;
	double **pi;
	double **sigma;
	int **kappa;
	
  p = current -> value.p;
  M = proposal -> value.M;
  pi = current -> value.pi;
  sigma = current -> value.sigma;
  kappa = current -> value.kappa;
#pragma omp parallel for schedule(guided) private(current_F_ST,proposal_F_ST,tmp,ratio,test,j,k)
  for (i = 0; i < data.nbr_demes; ++i) {
    ++update -> M[i];
    current_F_ST = (double) 1.0 / (1 + current -> value.M[i]);
    tmp = current_F_ST + genrand_real2(mtss[omp_get_thread_num()]) * (2 * move.M[i]) - move.M[i];   // Reflecting boundaries (left side)
    if (tmp < LIMIT) tmp = 2.0 * LIMIT - tmp;
    if (tmp > (1.0 - LIMIT)) tmp = 2.0 * (1.0 - LIMIT) - tmp;                                                  // Reflecting boundaries (right side)
    proposal_F_ST = tmp;
    M[i] = (1 - proposal_F_ST) / proposal_F_ST;
    proposal -> psi.Gamma_M[i] = lgammafn(M[i]);
    ratio = 0.0;
    for(j = 0; j < data.nbr_loci; ++j) {
      proposal -> psi.log_1F1[i][j] = log(hypergeometric_1_F_1(M[i] * pi[j][kappa[i][j]],M[i],sigma[i][j]));
      proposal -> psi.Gamma_M_pi[i][j] = lgammafn(M[i] * pi[j][0]) + lgammafn(M[i] * pi[j][1]);
      for(k = 0; k < 2; ++k) {
        ratio += pi[j][k] * (M[i] - current -> value.M[i]) * log(p[i][j][k]);
      }
      ratio += (current -> psi.log_1F1[i][j] + proposal -> psi.Gamma_M[i] + current -> psi.Gamma_M_pi[i][j]);
      ratio -= (proposal -> psi.log_1F1[i][j] + current -> psi.Gamma_M[i] + proposal -> psi.Gamma_M_pi[i][j]);
    }
    test = (ratio >= 0) ? 1 : genrand_real1(mtss[omp_get_thread_num()]) < exp(ratio);								// Are the new parameter values accepted?
    if (test) {
      ++accept -> M[i];
      current -> value.M[i] = proposal -> value.M[i];
      current -> psi.Gamma_M[i] = proposal -> psi.Gamma_M[i];
      for(j = 0; j < data.nbr_loci; ++j) {
        current -> psi.log_1F1[i][j] = proposal -> psi.log_1F1[i][j];
        current -> psi.Gamma_M_pi[i][j] = proposal -> psi.Gamma_M_pi[i][j];
      }
    }
  }
}

void update_M_Vitalis_et_al_2014(markov_state_struct *current,
                                 markov_state_struct *proposal,
                                 updates_struct move,
                                 updates_struct *update,
                                 updates_struct *accept,
                                 data_struct data)

{
	int i,j,k;
	int test;
	double reverse;
	double ratio;
	double ***p;
	double *M;
	double **pi;
	double **sigma;
	int **kappa;
	
  p = current -> value.p;
  M = proposal -> value.M;
  pi = current -> value.pi;
  sigma = current -> value.sigma;
  kappa = current -> value.kappa;
#pragma omp parallel for schedule(guided) private(reverse,ratio,test,j,k)
  for (i = 0; i < data.nbr_demes; ++i) {
    ++update -> M[i];
    reverse = move.M[i] * rnorm(mtss[omp_get_thread_num()]);												// update THETA from a log-normal draw, and
    proposal -> value.M[i] = current -> value.M[i] * exp(reverse);						// then record the new value for M.
    if (proposal -> value.M[i] < min_M || proposal -> value.M[i] > max_M) {		// If the new value lies outside the interval [min_M,max_M] then
      test = 0;																			// keep the Markov chain unchanged.
      goto end;
    }
    proposal -> psi.Gamma_M[i] = lgammafn(M[i]);
    ratio = 0.0;
    for (j = 0; j < data.nbr_loci; ++j) {
      proposal -> psi.log_1F1[i][j] = log(hypergeometric_1_F_1(M[i] * pi[j][kappa[i][j]],M[i],sigma[i][j]));
      proposal -> psi.Gamma_M_pi[i][j] = lgammafn(M[i] * pi[j][0]) + lgammafn(M[i] * pi[j][1]);
      for (k = 0; k < 2; ++k) {
        ratio += pi[j][k] * (M[i] - current -> value.M[i]) * log(p[i][j][k]);
      }
      ratio += (current -> psi.log_1F1[i][j] + proposal -> psi.Gamma_M[i] + current -> psi.Gamma_M_pi[i][j]);
      ratio -= (proposal -> psi.log_1F1[i][j] + current -> psi.Gamma_M[i] + proposal -> psi.Gamma_M_pi[i][j]);
    }
    proposal -> prior.M[i] = log_uniform_prior(M[i]);
    ratio += proposal -> prior.M[i];
    ratio -= current -> prior.M[i];
    ratio += reverse;
    test = (ratio >= 0) ? 1 : genrand_real1(mtss[omp_get_thread_num()]) < exp(ratio);								// Are the new parameter values accepted?
    end :
    if (test) {
      ++accept -> M[i];
      current -> value.M[i] = proposal -> value.M[i];
      current -> psi.Gamma_M[i] = proposal -> psi.Gamma_M[i];
      for (j = 0; j < data.nbr_loci; ++j) {
        current -> psi.log_1F1[i][j] = proposal -> psi.log_1F1[i][j];
        current -> psi.Gamma_M_pi[i][j] = proposal -> psi.Gamma_M_pi[i][j];
      }
      current -> prior.M[i] = proposal -> prior.M[i];
    }
  }
}

void update_sigma(markov_state_struct *current,
                  markov_state_struct *proposal,
                  updates_struct move,
                  updates_struct *update,
                  updates_struct *accept,
                  data_struct data)

{
	int i,j,ij;
	int test;
	double reverse;
	double ratio;
	double ***p;
	double *M;
	double **pi;
	double **sigma;
	double *delta;
	int **kappa;
  
  p = current -> value.p;
  M = current -> value.M;
  sigma = proposal -> value.sigma;
  pi = current -> value.pi;
  delta = current -> value.delta;
  kappa = current -> value.kappa;

#pragma omp parallel for schedule(guided) private(reverse,ratio,test,i,j)

  for (ij = 0; ij < (data.nbr_demes * data.nbr_loci); ++ij) {
    i = ij / data.nbr_loci;
    j = ij % data.nbr_loci;
    ++update -> sigma[i][j];                                                  // Update sigma
    reverse = move.sigma[i][j] * rnorm(mtss[omp_get_thread_num()]);										// proposal sigma from a log-normal draw, and
    proposal -> value.sigma[i][j] = current -> value.sigma[i][j] * exp(reverse);				// then record the new value for sigma.
    if (proposal -> value.sigma[i][j] < 0.0 || proposal -> value.sigma[i][j] > max_sigma) {	// If the new value lies outside the interval [0,MAX_SIGMA] then
      test = 0;																		// keep the Markov chain unchanged.
      goto end;
    }
    ratio = p[i][j][kappa[i][j]] * (sigma[i][j] - current -> value.sigma[i][j]);
    proposal -> psi.log_1F1[i][j] = log(hypergeometric_1_F_1(M[i] * pi[j][kappa[i][j]],M[i],sigma[i][j]));
    ratio += (current -> psi.log_1F1[i][j] - proposal -> psi.log_1F1[i][j]);
    proposal -> prior.sigma[i][j] = log_exponential_prior(sigma[i][j],delta[j]);	/// !!!!! was current.value.delta[i]
    ratio += (proposal -> prior.sigma[i][j] - current -> prior.sigma[i][j]);
    ratio += reverse;
    test = (ratio >= 0) ? 1 : genrand_real1(mtss[omp_get_thread_num()]) < exp(ratio);							// Are the new parameter values accepted?
    end :
    if (test) {
      ++accept -> sigma[i][j];
      current -> value.sigma[i][j] = proposal -> value.sigma[i][j];
      current -> psi.log_1F1[i][j] = proposal -> psi.log_1F1[i][j];
      current -> prior.sigma[i][j] = proposal -> prior.sigma[i][j];
    }
  }
}

void update_delta(markov_state_struct *current,
                  markov_state_struct *proposal,
                  updates_struct move,
                  updates_struct *update,
                  updates_struct *accept,
                  data_struct data)
{
	int i,j;
	int test;
	double reverse;
	double ratio;
	double **sigma;
	double *delta;
  double lambda;
  
  sigma = current -> value.sigma;
  delta = proposal -> value.delta;
	lambda = current -> value.lambda;
#pragma omp parallel for schedule(guided) private(reverse,ratio,test,i)
  for (j = 0; j < data.nbr_loci; ++j) {
    ++update -> delta[j];                                     // Update M
    reverse = move.delta[j] * rnorm(mtss[omp_get_thread_num()]);										// proposal sigma from a log-normal draw, and
    proposal -> value.delta[j] = current -> value.delta[j] * exp(reverse);				// then record the new value for sigma.
    ratio = 0.0;
    for (i = 0; i < data.nbr_demes; ++i) {
      proposal -> prior.sigma[i][j] = log_exponential_prior(sigma[i][j],delta[j]);
      ratio += (proposal -> prior.sigma[i][j] - current -> prior.sigma[i][j]);
    }
    proposal -> prior.delta[j] = log_exponential_prior(delta[j],lambda);
    ratio += (proposal -> prior.delta[j] - current -> prior.delta[j]);
    ratio += reverse;
    test = (ratio >= 0) ? 1 : genrand_real1(mtss[omp_get_thread_num()]) < exp(ratio);								// Are the new parameter values accepted?
    if (test) {
      ++accept -> delta[j];
      current -> value.delta[j] = proposal -> value.delta[j];
      current -> prior.delta[j] = proposal -> prior.delta[j];
      for (i = 0; i < data.nbr_demes; ++i) {
        current -> prior.sigma[i][j] = proposal -> prior.sigma[i][j];
      }
    }
  }
}

void update_lambda_metropolis(markov_state_struct *current,
                              markov_state_struct *proposal,
                              updates_struct move,
                              updates_struct *update,
                              updates_struct *accept,
                              data_struct data)

{
	int j;
	int test;
	double reverse;
	double ratio;
	double *delta;
	
  ++update -> lambda;
	reverse = move.lambda * rnorm(mtss[omp_get_thread_num()]);										// proposal sigma from a log-normal draw, and
	proposal -> value.lambda = current -> value.lambda * exp(reverse);				// then record the new value for sigma.
	delta = current -> value.delta;
	ratio = 0.0;
	for (j = 0; j < data.nbr_loci ; ++j) {                                                 // Update P(DELTA[i] | LAMBDA)
		proposal -> prior.delta[j] = log_exponential_prior(delta[j],proposal -> value.lambda);
		ratio += (proposal -> prior.delta[j] - current -> prior.delta[j]);
	}
	proposal -> prior.lambda = log_exponential_prior(proposal -> value.lambda,capital_lambda);
	ratio += (proposal -> prior.lambda - current -> prior.lambda);
	ratio += reverse;
	test = (ratio >= 0) ? 1 : genrand_real1(mtss[omp_get_thread_num()]) < exp(ratio);								// Are the new parameter values accepted?
  if (test) {
    ++accept -> lambda;
    current -> value.lambda = proposal -> value.lambda;
    current -> prior.lambda = proposal -> prior.lambda;
    for (j = 0; j < data.nbr_loci; ++j) {
      current -> prior.delta[j] = proposal -> prior.delta[j];
    }
  }
}

void update_lambda_gibbs(markov_state_struct *current,
                         data_struct data)

{
	int j;
  double alpha,beta;
	double *delta;
  
  delta = current -> value.delta;
  alpha = inverse_gamma_shape + data.nbr_loci;
  beta = inverse_gamma_rate;
  for (j = 0; j < data.nbr_loci ; ++j) {
		beta += delta[j];
	}
  current -> value.lambda = beta / rgamma(alpha,1.0,mtss[omp_get_thread_num()]);
  for (j = 0; j < data.nbr_loci ; ++j) {                                                 // Update P(DELTA[i] | LAMBDA)
		current -> prior.delta[j] = log_exponential_prior(delta[j],current -> value.lambda);
	}
  current -> prior.lambda = log_inverse_gamma_prior(current -> value.lambda,inverse_gamma_shape,inverse_gamma_rate);
}

void update_kappa(markov_state_struct *current,
                  markov_state_struct *proposal,
                  data_struct data)

{
	int i,j,ij;
  double A,B,rho;
	double ***p;
	double *M;
	double **pi;
	int **kappa;
	double **sigma;
  
  p = current -> value.p;
  pi = current -> value.pi;
  kappa = current -> value.kappa;
  sigma = current -> value.sigma;
  M = current -> value.M;

#pragma omp parallel for schedule(guided) private(A,B,rho,i,j)

  for (ij = 0; ij < (data.nbr_demes * data.nbr_loci); ++ij) {
    i = ij / data.nbr_loci;
    j = ij % data.nbr_loci;
    proposal -> psi.log_1F1[i][j] = log(hypergeometric_1_F_1(M[i] * pi[j][(1 - kappa[i][j])],M[i],sigma[i][j]));
    if (kappa[i][j] == 0) {
      A = (sigma[i][j] * p[i][j][kappa[i][j]]) - current -> psi.log_1F1[i][j];
      B = (sigma[i][j] * p[i][j][(1 - kappa[i][j])]) - proposal -> psi.log_1F1[i][j];
    } else {
      A = (sigma[i][j] * p[i][j][(1 - kappa[i][j])]) - proposal -> psi.log_1F1[i][j];
      B = (sigma[i][j] * p[i][j][kappa[i][j]]) - current -> psi.log_1F1[i][j];
    }
    rho = (double) 1.0 / (1.0 + exp(B - A));
    if (genrand_real1(mtss[omp_get_thread_num()]) < rho) {
      if (kappa[i][j] != 0) {
        kappa[i][j] = 0;
        current -> psi.log_1F1[i][j] = proposal -> psi.log_1F1[i][j];
      }
    } else {
      if (kappa[i][j] != 1) {
        kappa[i][j] = 1;
        current -> psi.log_1F1[i][j] = proposal -> psi.log_1F1[i][j];
      }
    }
  }
}

void update_beta_parameters_Vitalis_unpublished(markov_state_struct *current,
                                                updates_struct move,
                                                updates_struct *update,
                                                updates_struct *accept,
                                                data_struct data)

{
  int j;
	int test;
  double mu,nu,proposal_mu,proposal_nu,proposal_prior,reverse,sum_log_pi[2];
	double ratio;
  double LIMIT = 1e-3;
  double **pi;

  pi = current -> value.pi;
  nu = current -> value.beta_a + current -> value.beta_b;
  mu = current -> value.beta_a / nu;
  sum_log_pi[0] = sum_log_pi[1] = 0.0;
  for (j = 0; j < data.nbr_loci ; ++j) {
    sum_log_pi[0] += log(pi[j][0]);
    sum_log_pi[1] += log(pi[j][1]);
	}
  current -> prior.beta = (double) data.nbr_loci * (lgammafn(nu) - lgammafn(mu * nu) - lgammafn((1.0 - mu) * nu)) + (mu * nu - 1.0) * sum_log_pi[0] + ((1.0 - mu) * nu - 1.0) * sum_log_pi[1];
  ++update -> beta_mu;
  proposal_mu = mu + genrand_real2(mtss[omp_get_thread_num()]) * (2.0 * move.beta_mu) - move.beta_mu;
  if (proposal_mu < LIMIT) proposal_mu = 2.0 * LIMIT - proposal_mu;
  if (proposal_mu > (1.0 - LIMIT)) proposal_mu = 2.0 * (1.0 - LIMIT) - proposal_mu;
  proposal_prior = (double) data.nbr_loci * (lgammafn(nu) - lgammafn(proposal_mu * nu) - lgammafn((1.0 - proposal_mu) * nu)) + (proposal_mu * nu - 1.0) * sum_log_pi[0] + ((1.0 - proposal_mu) * nu - 1.0) * sum_log_pi[1];
  ratio = proposal_prior - current -> prior.beta;
  test = (ratio >= 0) ? 1 : genrand_real1(mtss[omp_get_thread_num()]) < exp(ratio);
  if (test) {
    ++accept -> beta_mu;
    mu = proposal_mu;
    current -> prior.beta = proposal_prior;
  }
  ++update -> beta_nu;
  reverse = move.beta_nu * rnorm(mtss[omp_get_thread_num()]);
  proposal_nu = nu * exp(reverse);
  proposal_prior = (double) data.nbr_loci * (lgammafn(proposal_nu) - lgammafn(mu * proposal_nu) - lgammafn((1.0 - mu) * proposal_nu)) + (mu * proposal_nu - 1.0) * sum_log_pi[0] + ((1.0 - mu) * proposal_nu - 1.0) * sum_log_pi[1];
  ratio = proposal_prior - current -> prior.beta + reverse - proposal_nu + nu;
  test = (ratio >= 0) ? 1 : genrand_real1(mtss[omp_get_thread_num()]) < exp(ratio);
  if (test) {
    ++accept -> beta_nu;
    nu = proposal_nu;
    current -> prior.beta = proposal_prior;
  }
  current -> value.beta_a = mu * nu;
  current -> value.beta_b = (1.0 - mu) * nu;
  for (j = 0; j < data.nbr_loci; ++j) {
    current -> prior.pi[j][0] = log_beta_prior(pi[j],current -> value.beta_a,current -> value.beta_b);
  }
}

void update_beta_parameters_Gautier_unpublished(markov_state_struct *current,
                                                updates_struct move,
                                                updates_struct *update,
                                                updates_struct *accept,
                                                data_struct data)

{
  int j;
  double sum_log_pi[2];
  double b_inf,b_sup,diff_log,inv_bwd,inv_fwd,mu_cur,mu_out,phi_cur,phi_out;
  double LIMIT = 1e-3;
  double **pi;
  
  pi = current -> value.pi;
  phi_cur = current -> value.beta_a + current -> value.beta_b;
  mu_cur = current -> value.beta_a / phi_cur;
  
  sum_log_pi[0] = sum_log_pi[1] = 0.0;
  for (j = 0; j < data.nbr_loci ; ++j) {
    sum_log_pi[0] += log(pi[j][0]);
    sum_log_pi[1] += log(pi[j][1]);
	}
  b_inf = dbl_max(LIMIT,(mu_cur - move.beta_mu / 2));
  b_sup = dbl_min((1.0 - LIMIT),(mu_cur + move.beta_mu / 2));
  mu_out = b_inf + genrand_real2(mtss[omp_get_thread_num()]) * (b_sup - b_inf);
  inv_fwd = b_sup - b_inf;
  inv_bwd = dbl_min((1.0 - LIMIT),(mu_out + move.beta_mu / 2)) - dbl_max(LIMIT,(mu_out - move.beta_mu / 2));
  diff_log = log(inv_fwd) - log(inv_bwd) + data.nbr_loci * (lgammafn(mu_cur * phi_cur) + lgammafn((1.0 - mu_cur) * phi_cur) - lgammafn(mu_out * phi_cur) - lgammafn((1.0 - mu_out) * phi_cur)) + phi_cur * (mu_out - mu_cur) * (sum_log_pi[0] - sum_log_pi[1]);
  ++update -> beta_mu;
  if (log(genrand_real2(mtss[omp_get_thread_num()])) > diff_log) {
    mu_out = mu_cur;
  }
  else {
    ++accept -> beta_mu;
  }
  phi_out = exp(rnorm(mtss[omp_get_thread_num()]) * move.beta_nu + log(phi_cur));
  diff_log = log(phi_out) - log(phi_cur) + data.nbr_loci * (lgammafn(mu_out * phi_cur) + lgammafn((1.0 - mu_out) * phi_cur) - lgammafn(mu_out * phi_out) - lgammafn((1.0 - mu_out) * phi_out) + lgammafn(phi_out) - lgammafn(phi_cur)) + (phi_out - phi_cur) * (mu_out * sum_log_pi[0] + (1.0 - mu_out) * sum_log_pi[1] - 1.0);
  ++update -> beta_nu;
  if (log(genrand_real2(mtss[omp_get_thread_num()])) > diff_log) {
    phi_out = phi_cur;
  }
  else {
    ++accept -> beta_nu;
  }
  current -> value.beta_a = mu_out * phi_out;
  current -> value.beta_b = (1.0 - mu_out) * phi_out;
  for (j = 0; j < data.nbr_loci; ++j) {
    current -> prior.pi[j][0] = log_beta_prior(pi[j],current -> value.beta_a,current -> value.beta_b);
  }
}
