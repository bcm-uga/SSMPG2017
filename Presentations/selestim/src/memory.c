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

#include "memory.h"

void allocate_memory(data_struct data,
                     markov_state_struct *current,
                     markov_state_struct *proposal,
                     updates_struct *update,
                     updates_struct *accept,
                     updates_struct *move,
                     moments_struct *postmean,
                     moments_struct *postmsqr)

{
	int i,j;
  
  binom_coeff_counts = (double **) calloc (data.nbr_demes,sizeof(double *));
  for (i = 0; i < data.nbr_demes; ++i) {
    binom_coeff_counts[i] = (double *) calloc (data.nbr_loci,sizeof(double));
  }
	current -> value.p = (double ***) malloc(data.nbr_demes * sizeof(double **));
	current -> value.M = (double *) malloc(data.nbr_demes * sizeof(double));
	current -> value.delta = (double *) malloc(data.nbr_loci * sizeof(double));
	current -> value.pi = (double **) malloc(data.nbr_loci * sizeof(double *));
	current -> value.sigma = (double **) malloc(data.nbr_demes * sizeof(double *));
	current -> value.kappa = (int **) malloc(data.nbr_demes * sizeof(int *));
	for (j = 0; j < data.nbr_loci; ++j) {
		current -> value.pi[j] = (double *) malloc(2 * sizeof(double));
	}
	for (i = 0; i < data.nbr_demes; ++i) {
		current -> value.p[i] = (double **) malloc(data.nbr_loci * sizeof(double *));
		current -> value.sigma[i] = (double *) malloc(data.nbr_loci * sizeof(double));
		current -> value.kappa[i] = (int *) malloc(data.nbr_loci * sizeof(int));
		
	}
	for (i = 0; i < data.nbr_demes; ++i) {
		for (j = 0; j < data.nbr_loci; ++j) {
			current -> value.p[i][j] = (double *) malloc(2 * sizeof(double));
		}
	}
	update -> p = (double **) malloc(data.nbr_demes * sizeof(double *));
	accept -> p = (double **) malloc(data.nbr_demes * sizeof(double *));
	update -> M = (double *) malloc(data.nbr_demes * sizeof(double));
	accept -> M = (double *) malloc(data.nbr_demes * sizeof(double));
	update -> pi = (double *) malloc(data.nbr_loci * sizeof(double));
	accept -> pi = (double *) malloc(data.nbr_loci * sizeof(double));
	update -> delta = (double *) malloc(data.nbr_loci * sizeof(double));
	accept -> delta = (double *) malloc(data.nbr_loci * sizeof(double));
	update -> sigma = (double **) malloc(data.nbr_demes * sizeof(double *));
	accept -> sigma = (double **) malloc(data.nbr_demes * sizeof(double *));
	for (i = 0; i < data.nbr_demes; ++i) {
		update -> p[i] = (double *) malloc(data.nbr_loci * sizeof(double));
		accept -> p[i] = (double *) malloc(data.nbr_loci * sizeof(double));
		update -> sigma[i] = (double *) malloc(data.nbr_loci * sizeof(double));
		accept -> sigma[i] = (double *) malloc(data.nbr_loci * sizeof(double));
	}
	move -> p = (double **) malloc(data.nbr_demes * sizeof(double *));
	move -> M = (double *) malloc(data.nbr_demes * sizeof(double));
	move -> pi = (double *) malloc(data.nbr_loci * sizeof(double));
	move -> delta = (double *) malloc(data.nbr_loci * sizeof(double));
	move -> sigma = (double **) malloc(data.nbr_demes * sizeof(double *));
	for (i = 0; i < data.nbr_demes; ++i) {
		move -> sigma[i] = (double *) malloc(data.nbr_loci * sizeof(double));
		move -> p[i] = (double *) malloc(data.nbr_loci * sizeof(double));
	}
	proposal -> value.p = (double ***) malloc(data.nbr_demes * sizeof(double **));
	proposal -> value.pi = (double **) malloc(data.nbr_loci * sizeof(double *));
	proposal -> value.kappa = (int **) malloc(data.nbr_demes * sizeof(int *));
	proposal -> value.sigma = (double **) malloc(data.nbr_demes * sizeof(double *));
	proposal -> value.M = (double *) malloc(data.nbr_demes * sizeof(double));
	proposal -> value.delta = (double *) malloc(data.nbr_loci * sizeof(double));
	for (j = 0; j < data.nbr_loci; ++j) {
		proposal -> value.pi[j] = (double *) malloc(2 * sizeof(double));
	}
	for (i = 0; i < data.nbr_demes; ++i) {
		proposal -> value.p[i] = (double **) malloc(data.nbr_loci * sizeof(double *));
		proposal -> value.kappa[i] = (int *) malloc(data.nbr_loci * sizeof(int));
		proposal -> value.sigma[i] = (double *) malloc(data.nbr_loci * sizeof(double));
	}
	for (i = 0; i < data.nbr_demes; ++i) {
		for (j = 0; j < data.nbr_loci; ++j) {
			proposal -> value.p[i][j] = (double *) malloc(2 * sizeof(double));
		}
	}
  current -> prior.pi = (double **) malloc(data.nbr_loci * sizeof(double *));
	current -> prior.sigma = (double **) malloc(data.nbr_demes * sizeof(double *));
	current -> prior.M = (double *) malloc(data.nbr_demes * sizeof(double));
	current -> prior.delta = (double *) malloc(data.nbr_loci * sizeof(double));
  proposal -> prior.pi = (double **) malloc(data.nbr_loci * sizeof(double *));
	proposal -> prior.sigma = (double **) malloc(data.nbr_demes * sizeof(double *));
	proposal -> prior.M = (double *) malloc(data.nbr_demes * sizeof(double));
	proposal -> prior.delta = (double *) malloc(data.nbr_loci * sizeof(double));
	for (i = 0; i < data.nbr_demes; ++i) {
		current -> prior.sigma[i] = (double *) malloc(data.nbr_loci * sizeof(double));
		proposal -> prior.sigma[i] = (double *) malloc(data.nbr_loci * sizeof(double));
	}
	for (j = 0; j < data.nbr_loci; ++j) {
		current -> prior.pi[j] = (double *) malloc(2 * sizeof(double));
		proposal -> prior.pi[j] = (double *) malloc(2 * sizeof(double));
	}
	current -> psi.log_1F1 = (double **) malloc(data.nbr_demes * sizeof(double *));
	current -> psi.Gamma_M = (double *) malloc(data.nbr_demes * sizeof(double));
	current -> psi.Gamma_M_pi = (double **) malloc(data.nbr_demes * sizeof(double *));
	proposal -> psi.log_1F1 = (double **) malloc(data.nbr_demes * sizeof(double *));
	proposal -> psi.Gamma_M = (double *) malloc(data.nbr_demes * sizeof(double));
	proposal -> psi.Gamma_M_pi = (double **) malloc(data.nbr_demes * sizeof(double *));
	for (i = 0; i < data.nbr_demes; ++i) {
		current -> psi.log_1F1[i] = (double *) malloc(data.nbr_loci * sizeof(double));
		current -> psi.Gamma_M_pi[i] = (double *) malloc(data.nbr_loci * sizeof(double));
		proposal -> psi.log_1F1[i] = (double *) malloc(data.nbr_loci * sizeof(double));
		proposal -> psi.Gamma_M_pi[i] = (double *) malloc(data.nbr_loci * sizeof(double));
	}
  postmean -> M = (double *) malloc(data.nbr_demes * sizeof(double));
  postmsqr -> M = (double *) malloc(data.nbr_demes * sizeof(double));
  postmean -> delta = (double *) malloc(data.nbr_loci * sizeof(double));
  postmsqr -> delta = (double *) malloc(data.nbr_loci * sizeof(double));
  postmean -> pi = (double *) malloc(data.nbr_loci * sizeof(double));
  postmsqr -> pi = (double *) malloc(data.nbr_loci * sizeof(double));
  postmean -> kappa = (double **) malloc(data.nbr_demes * sizeof(double *));
  for (i = 0; i < data.nbr_demes; ++i) {
    postmean -> kappa[i] = (double *) malloc(data.nbr_loci * sizeof(double));
  }
  postmean -> p = (double **) malloc(data.nbr_demes * sizeof(double *));
  postmsqr -> p = (double **) malloc(data.nbr_demes * sizeof(double *));
  postmean -> sigma = (double ***) malloc(data.nbr_demes * sizeof(double **));
  postmsqr -> sigma = (double ***) malloc(data.nbr_demes * sizeof(double **));
  for (i = 0; i < data.nbr_demes; ++i) {
    postmean -> p[i] = (double *) malloc(data.nbr_loci * sizeof(double));
    postmsqr -> p[i] = (double *) malloc(data.nbr_loci * sizeof(double));
    postmean -> sigma[i] = (double **) malloc(data.nbr_loci * sizeof(double *));
    postmsqr -> sigma[i] = (double **) malloc(data.nbr_loci * sizeof(double *));
  }
  for (i = 0; i < data.nbr_demes; ++i) {
		for (j = 0; j < data.nbr_loci; ++j) {
      postmean -> sigma[i][j] = (double *) malloc(3 * sizeof(double));
      postmsqr -> sigma[i][j] = (double *) malloc(3 * sizeof(double));
    }
  }
  if (pooled_data) {
    binom_coeff_reads = (double **) calloc (data.nbr_demes,sizeof(double *));
    for (i = 0; i < data.nbr_demes; ++i) {
      binom_coeff_reads[i] = (double *) calloc (data.nbr_loci,sizeof(double));
    }
    current -> value.counts = (int ***) malloc(data.nbr_demes * sizeof(int **));
    proposal -> value.counts = (int ***) malloc(data.nbr_demes * sizeof(int **));
    move -> counts = (double **) malloc(data.nbr_demes * sizeof(double *));
    update -> counts = (double **) malloc(data.nbr_demes * sizeof(double *));
    accept -> counts = (double **) malloc(data.nbr_demes * sizeof(double *));
    postmean -> counts = (double **) malloc(data.nbr_demes * sizeof(double *));
    postmsqr -> counts = (double **) malloc(data.nbr_demes * sizeof(double *));
    for (i = 0; i < data.nbr_demes; ++i) {
      current -> value.counts[i] = (int **) malloc(data.nbr_loci * sizeof(int *));
      proposal -> value.counts[i] = (int **) malloc(data.nbr_loci * sizeof(int *));
      move -> counts[i] = (double *) malloc(data.nbr_loci * sizeof(double));
      update -> counts[i] = (double *) malloc(data.nbr_loci * sizeof(double));
      accept -> counts[i] = (double *) malloc(data.nbr_loci * sizeof(double));
      postmean -> counts[i] = (double *) malloc(data.nbr_loci * sizeof(double));
      postmsqr -> counts[i] = (double *) malloc(data.nbr_loci * sizeof(double));
      for (j = 0; j < data.nbr_loci; ++j) {
        current -> value.counts[i][j] = (int *) malloc(2 * sizeof(int));
        proposal -> value.counts[i][j] = (int *) malloc(2 * sizeof(int));
      }
    }
  }
}

void release_memory(data_struct *data,
                    markov_state_struct *current,
                    markov_state_struct *proposal,
                    updates_struct *update,
                    updates_struct *accept,
                    updates_struct *move,
                    moments_struct *postmean,
                    moments_struct *postmsqr)

{
	int i,j;
	
  for (i = 0; i < data -> nbr_demes; ++i) {
    free(binom_coeff_counts[i]);
  }
  free(binom_coeff_counts);
	for (i = 0; i < data -> nbr_demes; ++i) {
		for (j = 0; j < data -> nbr_loci; ++j) {
			free(current -> value.p[i][j]);
		}
		free(current -> value.p[i]);
	}
	free(current -> value.p);
	for (j = 0; j < data -> nbr_loci; ++j) {
		free(current -> value.pi[j]);
	}
	free(current -> value.pi);
	for (i = 0; i < data -> nbr_demes; ++i) {
		free(current -> value.sigma[i]);
		free(current -> value.kappa[i]);
	}
	free(current -> value.sigma);
	free(current -> value.kappa);
	free(current -> value.M);
	free(current -> value.delta);
	for (i = 0; i < data -> nbr_demes; ++i) {
		free(update -> p[i]);
		free(accept -> p[i]);
		free(update -> sigma[i]);
		free(accept -> sigma[i]);
	}
	free(update -> p);
	free(accept -> p);
	free(update -> sigma);
	free(accept -> sigma);
	free(update -> M);
	free(accept -> M);
	free(update -> pi);
	free(accept -> pi);
	free(update -> delta);
	free(accept -> delta);
	for (i = 0; i < data -> nbr_demes; ++i) {
		free(move -> p[i]);
		free(move -> sigma[i]);
	}
	free(move -> p);
	free(move -> sigma);
	free(move -> pi);
	free(move -> M);
	free(move -> delta);
	for (i = 0; i < data -> nbr_demes; ++i) {
		for (j = 0; j < data -> nbr_loci; ++j) {
			free(proposal -> value.p[i][j]);
		}
		free(proposal -> value.p[i]);
	}
	free(proposal -> value.p);
	for (j = 0; j < data -> nbr_loci; ++j) {
		free(proposal -> value.pi[j]);
	}
	free(proposal -> value.pi);
	for (i = 0; i < data -> nbr_demes; ++i) {
		free(proposal -> value.kappa[i]);
		free(proposal -> value.sigma[i]);
	}
	free(proposal -> value.kappa);
	free(proposal -> value.sigma);
	free(proposal -> value.M);
	free(proposal -> value.delta);
	for (j = 0; j < data -> nbr_loci; ++j) {
		free(current -> prior.pi[j]);
		free(proposal -> prior.pi[j]);
	}
  for (i = 0; i < data -> nbr_demes; ++i) {
		free(current -> prior.sigma[i]);
		free(proposal -> prior.sigma[i]);
	}
  free(current -> prior.pi);
	free(current -> prior.sigma);
	free(current -> prior.M);
	free(current -> prior.delta);
  free(proposal -> prior.pi);
	free(proposal -> prior.sigma);
	free(proposal -> prior.M);
	free(proposal -> prior.delta);
	for (i = 0; i < data -> nbr_demes; ++i) {
		free(current -> psi.log_1F1[i]);
		free(current -> psi.Gamma_M_pi[i]);
		free(proposal -> psi.log_1F1[i]);
		free(proposal -> psi.Gamma_M_pi[i]);
	}
	free(current -> psi.log_1F1);
	free(current -> psi.Gamma_M);
	free(current -> psi.Gamma_M_pi);
	free(proposal -> psi.log_1F1);
	free(proposal -> psi.Gamma_M);
	free(proposal -> psi.Gamma_M_pi);
	for (i = 0; i < data -> nbr_demes; ++i) {
		for (j = 0; j < data -> nbr_loci; ++j) {
			free(data -> counts[i][j]);
		}
	}
	for (i = 0; i < data -> nbr_demes; ++i) {
		free(data -> total_nbr_counts[i]);
		free(data -> counts[i]);
	}
	free(data -> total_nbr_counts);
	free(data -> counts);
  free(postmean -> M);
  free(postmsqr -> M);
  free(postmean -> delta);
  free(postmsqr -> delta);
  free(postmean -> pi);
  free(postmsqr -> pi);
  for (i = 0; i < data -> nbr_demes; ++i) {
    free(postmean -> kappa[i]);
  }
  free(postmean -> kappa);
  for (i = 0; i < data -> nbr_demes; ++i) {
		for (j = 0; j < data -> nbr_loci; ++j) {
      free(postmean -> sigma[i][j]);
      free(postmsqr -> sigma[i][j]);
    }
  }
  for (i = 0; i < data -> nbr_demes; ++i) {
    free(postmean -> p[i]);
    free(postmsqr -> p[i]);
    free(postmean -> sigma[i]);
    free(postmsqr -> sigma[i]);
  }
  free(postmean -> p);
  free(postmsqr -> p);
  free(postmean -> sigma);
  free(postmsqr -> sigma);
  if (pooled_data) {
    for (i = 0; i < data -> nbr_demes; ++i) {
      free(binom_coeff_reads[i]);
    }
    free(binom_coeff_reads);
    for (i = 0; i < data -> nbr_demes; ++i) {
      for (j = 0; j < data -> nbr_loci; ++j) {
        free(current -> value.counts[i][j]);
        free(proposal -> value.counts[i][j]);
      }
      free(current -> value.counts[i]);
      free(proposal -> value.counts[i]);
      free(move -> counts[i]);
      free(update -> counts[i]);
      free(accept -> counts[i]);
      free(postmean -> counts[i]);
      free(postmsqr -> counts[i]);
    }
    free(current -> value.counts);
    free(proposal -> value.counts);
    free(move -> counts);
    free(update -> counts);
    free(accept -> counts);
    free(postmean -> counts);
    free(postmsqr -> counts);
    for (i = 0; i < data -> nbr_demes; ++i) {
      for (j = 0; j < data -> nbr_loci; ++j) {
        free(data -> reads[i][j]);
        free(min_counts[i][j]);
        free(max_counts[i][j]);
      }
    }
    for (i = 0; i < data -> nbr_demes; ++i) {
      free(data -> total_nbr_reads[i]);
      free(data -> reads[i]);
      free(min_counts[i]);
      free(max_counts[i]);
    }
    free(data -> total_nbr_reads);
    free(data -> reads);
    free(data -> pool_size);
    free(min_counts);
    free(max_counts);
  }
}
