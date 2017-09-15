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

#ifndef _DEFS_H
#define _DEFS_H

#ifndef MAIN_UNIT
#define GLOBAL_VARIABLE extern
#else
#define GLOBAL_VARIABLE
#endif

#define update_M                   update_M_Vitalis_et_al_2014
#define update_counts              update_counts_Vitalis_unpublished
#define update_beta_parameters     update_beta_parameters_Vitalis_unpublished

#define VERSION              "1.1.7"
#define SCALING              1.25
#define MIN_ACCEPT           0.25
#define MAX_ACCEPT           0.40
#define PREC                 1e-15
#define MIN_P                1e-8
#define MIN_PI               1e-4
#define OUTPUTS_MCMC         "diag_mcmc.log"
#define TRACE_M              "trace_M.out"
#define TRACE_PI             "trace_pi.out"
#define TRACE_SIGMA          "trace_sigma.out"
#define TRACE_DELTA          "trace_delta.out"
#define TRACE_LAMBDA         "trace_lambda.out"
#define TRACE_KAPPA          "trace_kappa.out"
#define TRACE_BETA           "trace_beta.out"
#define TRACE_FREQ           "trace_freq.out"
#define TRACE_COUNTS         "trace_counts.out"
#define LOGFILE              "logfile.log"
#define PODFILE              "pod_"
#define SUMRY_M              "summary_M.out"
#define SUMRY_PI             "summary_pi.out"
#define SUMRY_SIGMA          "summary_sigma.out"
#define SUMRY_DELTA          "summary_delta.out"
#define SUMRY_LAMBDA         "summary_lambda.out"
#define SUMRY_KAPPA          "summary_kappa.out"
#define SUMRY_BETA           "summary_beta.out"
#define SUMRY_FREQ           "summary_freq.out"
#define SUMRY_COUNTS         "summary_counts.out"
#define KLD_QUTILE           "KLD_quantiles.out"
#define README               "readme.txt"
#define POD_SUBDIRECTORY     "calibration"
#define FALSE                0
#define TRUE                 1
#define UNDEFINED           -999.0

#define repeat for (;;)

// ---------------------
// structure definitions
// ---------------------

typedef struct {
	int nbr_loci;           // Number of loci
	int nbr_demes;          // Number of sampled subpopoulations
  int *pool_size;         // Number of genes in each pool
	int **total_nbr_counts; // Sum of allele counts (per locus and per subpopulation)
	int **total_nbr_reads;  // Sum of read counts (per locus and per subpopulation)
	int ***counts;          // Allele counts (per locus, per subpopulation, and per allele)
	int ***reads;           // Read counts (per locus, per subpopulation, and per allele)
} data_struct;

typedef struct {
  double beta;
  double beta_a;
  double beta_b;
  int ***counts;
  double ***p;
	double **pi;            // pi[i][j] gives the frequency of allele j at locus i, in the metapopulation (migrant pool)
	double *M;              // M[i] = 4Nm in the ith subpopulation
	int **kappa;            // kappa[i][j] returns an indicator of the allele which is selected (0 or 1) for the ith locus in the jth population
	double **sigma;         // sigma[i][j] = 4Ns gives the intensity of selection for the ith locus in the jth subpopulation
	double *delta;          // delta[i] gives the distribution of selection coefficients for the ith locus
	double lambda;          // Hyper-parameter for the selection coefficients
} parameters_struct;

typedef struct {
  double beta_a;
  double beta_b;
  double *pi;             // pi[i][j] gives the frequency of allele j at locus i, in the metapopulation (migrant pool)
	double *M;              // M[j] = 4Nm in the jth subpopulation
	double *delta;          // delta[i] gives the distribution of selection coefficients for the ith locus
	double **kappa;         // kappa[i][j] returns an indicator of the allele which is selected (0 or 1) for the ith locus in the jth population
  double **counts;
  double **p;
	double ***sigma;        // sigma[i][j] = 4Ns gives the intensity of selection for the ith locus in the jth subpopulation
	double lambda;          // Hyper-parameter for the selection coefficients
} moments_struct;

typedef struct {
  double **counts;
	double **p;
	double *M;
	double *pi;
	double *delta;
	double **sigma;
	double lambda;
  double beta_mu;
  double beta_nu;
} updates_struct;

typedef struct {
	double *Gamma_M;
	double **log_1F1;
	double **Gamma_M_pi;
} psi_struct;

typedef struct {
	psi_struct psi;
	parameters_struct value;
	parameters_struct prior;
} markov_state_struct;

// -----------------------------
// local functions to selestim.c
// -----------------------------

void do_pilot_runs(data_struct data,markov_state_struct *current,markov_state_struct *proposal,updates_struct *update,updates_struct *accept,updates_struct *move);
void run_mcmc(data_struct data,markov_state_struct *current,markov_state_struct *proposal,moments_struct *postmean,moments_struct *postmsqr,updates_struct *update,updates_struct *accept,updates_struct *move);

// ----------------
// GLOBAL VARIABLES
// ----------------

GLOBAL_VARIABLE unsigned int fixed_beta;
GLOBAL_VARIABLE unsigned int fixed_lambda;
GLOBAL_VARIABLE unsigned int verbose;
GLOBAL_VARIABLE unsigned int calibration;
GLOBAL_VARIABLE unsigned int calibration_only;
GLOBAL_VARIABLE unsigned int pooled_data;
GLOBAL_VARIABLE unsigned int pod_nbr_loci;
GLOBAL_VARIABLE unsigned int n_pilot;
GLOBAL_VARIABLE unsigned int l_pilot;
GLOBAL_VARIABLE unsigned int step;
GLOBAL_VARIABLE unsigned int exponential_prior;
GLOBAL_VARIABLE unsigned int inverse_gamma_prior;
GLOBAL_VARIABLE unsigned long burn_in;
GLOBAL_VARIABLE unsigned long chain_length;
GLOBAL_VARIABLE unsigned long total_chain_length;
GLOBAL_VARIABLE int init_delta_counts;
GLOBAL_VARIABLE int min_nbr_counts;
GLOBAL_VARIABLE int min_nbr_reads;
GLOBAL_VARIABLE double init_delta_p;
GLOBAL_VARIABLE double init_delta_M;
GLOBAL_VARIABLE double init_delta_pi;
GLOBAL_VARIABLE double init_delta_sigma;
GLOBAL_VARIABLE double init_delta_delta;
GLOBAL_VARIABLE double init_delta_lambda;
GLOBAL_VARIABLE double init_delta_beta_mu;
GLOBAL_VARIABLE double init_delta_beta_nu;
GLOBAL_VARIABLE double min_M;
GLOBAL_VARIABLE double max_M;
GLOBAL_VARIABLE double max_sigma;
GLOBAL_VARIABLE double default_init_delta_lambda;
GLOBAL_VARIABLE double default_init_delta_beta_mu;
GLOBAL_VARIABLE double default_init_delta_beta_nu;
GLOBAL_VARIABLE double default_capital_lambda;
GLOBAL_VARIABLE double default_inverse_gamma_shape;
GLOBAL_VARIABLE double default_inverse_gamma_rate;
GLOBAL_VARIABLE double default_beta_a;
GLOBAL_VARIABLE double default_beta_b;
GLOBAL_VARIABLE double beta_a;
GLOBAL_VARIABLE double beta_b;
GLOBAL_VARIABLE double inverse_gamma_shape;
GLOBAL_VARIABLE double inverse_gamma_rate;
GLOBAL_VARIABLE double fixed_lambda_value;
GLOBAL_VARIABLE double capital_lambda;
GLOBAL_VARIABLE double **binom_coeff_counts;
GLOBAL_VARIABLE double **binom_coeff_reads;
GLOBAL_VARIABLE const char *lambda_prior;
GLOBAL_VARIABLE const char *program_name;
GLOBAL_VARIABLE FILE *logfile;
GLOBAL_VARIABLE int ***min_counts;
GLOBAL_VARIABLE int ***max_counts;

#endif
