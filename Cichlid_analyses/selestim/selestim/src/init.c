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

#include "init.h"

void init_parameters(markov_state_struct *current,
                     data_struct data)

{
	int i,j;
  int counts,sum_counts;
  double freq;
  double mean,mean_sqr,tmp,var;
	int **n;
	int ***x;
	double ***p;
	double *M;
	double **pi;
	double **sigma;
	int **kappa;
	double *delta;
  double *c;
  
	p = current -> value.p;
	M = current -> value.M;
	pi = current -> value.pi;
	sigma = current -> value.sigma;
	delta = current -> value.delta;
	kappa = current -> value.kappa;
  n = data.total_nbr_counts;
  if (pooled_data) {
    x = current -> value.counts;
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        if (data.total_nbr_reads[i][j] > 0) {
          freq = (double) data.reads[i][j][0] / data.total_nbr_reads[i][j];
          x[i][j][0] = round(freq * n[i][j]);
          if ((data.reads[i][j][0] > 0) && (x[i][j][0] == 0)) {
            x[i][j][0] = 1;
          }
          if ((data.reads[i][j][0] < data.total_nbr_reads[i][j]) && (x[i][j][0] == n[i][j])) {
            x[i][j][0] = n[i][j] - 1;
          }
          x[i][j][1] = n[i][j] - x[i][j][0];
        }
        else {
          x[i][j][0] = x[i][j][1] = 0;
        }

      }
    }
  }
  else {
    x = data.counts;
  }
  for (i = 0; i < data.nbr_demes; ++i) {
    for (j = 0; j < data.nbr_loci; ++j) {
      p[i][j][0] = (double) (x[i][j][0] + 1.0) / (n[i][j] + 2.0);
      p[i][j][1] = 1.0 - p[i][j][0];
    }
  }
	for (j = 0; j < data.nbr_loci; ++j) {
    counts = sum_counts = 0;
    for (i = 0; i < data.nbr_demes; ++i) {
      counts += x[i][j][0];
      sum_counts += n[i][j];
    }
    pi[j][0] = (double) (counts + 1.0) / (sum_counts + 2.0);
    pi[j][1] = 1.0 - pi[j][0];
	}
  c = (double *) malloc(data.nbr_demes * sizeof(double));
  for (i = 0; i < data.nbr_demes; ++i) {
    c[i] = 0.0;
    for (j = 0; j < data.nbr_loci; ++j) {
      c[i] += pow((p[i][j][0] - pi[j][0]),2) / (pi[j][0] * (1 - pi[j][0]));
    }
    c[i] /= data.nbr_loci;
  }
  for (i = 0; i < data.nbr_demes; ++i) {
		M[i] = (1 - c[i]) / c[i];
    if (M[i] < min_M) M[i] = min_M;
    if (M[i] > max_M) M[i] = max_M;
  }
  free(c);
  if (fixed_lambda) {
    current -> value.lambda = fixed_lambda_value;
  }
  else {
    if (inverse_gamma_prior) {
      current -> value.lambda = (double) inverse_gamma_rate / (inverse_gamma_shape - 1.0);
    }
    else if (exponential_prior) {
      current -> value.lambda = (double) capital_lambda;
    }
  }
  for (j = 0; j < data.nbr_loci; ++j) {
		delta[j] = current -> value.lambda;
		for (i = 0; i < data.nbr_demes; ++i) {
      kappa[i][j] = 1 - round(p[i][j][0]);
      sigma[i][j] = delta[j];
		}
	}
  if (fixed_beta) {
    current -> value.beta_a = beta_a;
    current -> value.beta_b = beta_b;
  }
  else {
    mean = mean_sqr = 0.0;
    for (j = 0; j < data.nbr_loci; ++j) {
      mean += pi[j][0];
      mean_sqr += pi[j][0] * pi[j][0];
    }
    mean /= data.nbr_loci;
    mean_sqr /= data.nbr_loci;
    var = (mean_sqr - mean * mean) * data.nbr_loci / (data.nbr_loci - 1.0);
    tmp = mean * (1.0 - mean) / var - 1.0;
    current -> value.beta_a = mean * tmp;
    current -> value.beta_b = (1.0 - mean) * tmp;
  }
  for (j = 0; j < data.nbr_loci; ++j) {
    current -> prior.pi[j][0] = log_beta_prior(pi[j],current -> value.beta_a,current -> value.beta_b);
  }
	for (i = 0; i < data.nbr_demes; ++i) {
		current -> prior.M[i] = log_uniform_prior(M[i]);
    current -> psi.Gamma_M[i] = lgammafn(M[i]);
		for (j = 0; j < data.nbr_loci; ++j) {
			current -> psi.log_1F1[i][j] = log(hypergeometric_1_F_1(M[i] * pi[j][kappa[i][j]],M[i],sigma[i][j]));
			current -> psi.Gamma_M_pi[i][j] = lgammafn(M[i] * pi[j][0]) + lgammafn(M[i] * pi[j][1]);
			current -> prior.sigma[i][j] = log_exponential_prior(sigma[i][j],delta[j]);
		}
	}
	for (j = 0; j < data.nbr_loci; ++j) {
		current -> prior.delta[j] = log_exponential_prior(delta[j],current -> value.lambda);
	}
  if (!fixed_lambda) {
    if (inverse_gamma_prior) {
      current -> prior.lambda = log_inverse_gamma_prior(current -> value.lambda,inverse_gamma_shape,inverse_gamma_rate);
    }
    else if (exponential_prior) {
      current -> prior.lambda = log_exponential_prior(current -> value.lambda,capital_lambda);
    }
  }
  else {
    current -> prior.lambda = 0.0;
  }
  if (pooled_data) {
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        binom_coeff_reads[i][j] = lgammafn(data.total_nbr_reads[i][j] + 1) - lgammafn(data.reads[i][j][0] + 1) - lgammafn(data.reads[i][j][1] + 1);
      }
    }
  }
  else {
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        binom_coeff_counts[i][j] = lgammafn(n[i][j] + 1) - lgammafn(x[i][j][0] + 1) - lgammafn(x[i][j][1] + 1);
      }
    }
  }
}

void init_moves(data_struct data,
                updates_struct *move)

{
  int i,j;
  
  for (i = 0; i < data.nbr_demes; ++i) {
    move -> M[i] = init_delta_M;
  }
  for (j = 0; j < data.nbr_loci; ++j) {
    move -> pi[j] = init_delta_pi;
    move -> delta[j] = init_delta_delta;
  }
  for (i = 0; i < data.nbr_demes; ++i) {
    for (j = 0; j < data.nbr_loci; ++j) {
      if (pooled_data) {
        if (data.total_nbr_counts[i][j] > 0) {
          move -> counts[i][j] = int_min((data.total_nbr_counts[i][j] - 1),init_delta_counts);
        }
        else {
          move -> counts[i][j] = 0;
        }
      }
      move -> p[i][j] = init_delta_p;
      move -> sigma[i][j] = init_delta_sigma;
    }
  }
  move -> lambda = init_delta_lambda;
  if (!fixed_beta) {
    move -> beta_mu = init_delta_beta_mu;
    move -> beta_nu = init_delta_beta_nu;
  }  
}

void init_moments(data_struct data,
                  moments_struct *postmean,
                  moments_struct *postmsqr)

{
  int i,j;
  
  for (i = 0; i < data.nbr_demes; ++i) {
    postmean -> M[i] = postmsqr -> M[i] = 0.0;
  }
	for (j = 0; j < data.nbr_loci; ++j) {
    postmean -> pi[j] = postmsqr -> pi[j] = 0.0;
    postmean -> delta[j] = postmsqr -> delta[j] = 0.0;
  }
  for (i = 0; i < data.nbr_demes; ++i) {
		for (j = 0; j < data.nbr_loci; ++j) {
      postmean -> kappa[i][j] = 0.0;
      postmean -> p[i][j] = postmsqr -> p[i][j] = 0.0;
      postmean -> sigma[i][j][0] = postmsqr -> sigma[i][j][0] = 0.0;
      postmean -> sigma[i][j][1] = postmsqr -> sigma[i][j][1] = 0.0;
      postmean -> sigma[i][j][2] = postmsqr -> sigma[i][j][2] = 0.0;
    }
  }
  if (pooled_data) {
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        postmean -> counts[i][j] = postmsqr -> counts[i][j] = 0.0;
      }
    }
  }
  postmean -> lambda = postmsqr -> lambda = 0.0;
  if (!fixed_beta) {
    postmean -> beta_a = postmsqr -> beta_a = 0.0;
    postmean -> beta_b = postmsqr -> beta_b = 0.0;
  }
}

void init_updates(data_struct data,
                  updates_struct *update,
                  updates_struct *accept)

{
  int i,j;

  for (i = 0; i < data.nbr_demes; ++i) {
    update -> M[i] = accept -> M[i] = 0.0;
  }
  for (j = 0; j < data.nbr_loci; ++j) {
    update -> pi[j] = accept -> pi[j] = 0.0;
    update -> delta[j] = accept -> delta[j] = 0.0;
  }
  for (i = 0; i < data.nbr_demes; ++i) {
    for (j = 0; j < data.nbr_loci; ++j) {
      update -> p[i][j] = accept -> p[i][j] = 0.0;
      update -> sigma[i][j] = accept -> sigma[i][j] = 0.0;
    }
  }
  if (pooled_data) {
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        update -> counts[i][j] = accept -> counts[i][j] = 0.0;
      }
    }
  }
  if (!fixed_lambda && exponential_prior) {
    update -> lambda = accept -> lambda = 0.0;
  }
  if (!fixed_beta) {
    update -> beta_mu = accept -> beta_mu = 0.0;
    update -> beta_nu = accept -> beta_nu = 0.0;
  }
}
