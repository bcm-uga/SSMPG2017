/*
 Copyright INRA
 author: Renaud VITALIS (2013)
 
 renaud.vitalis@inra.fr
 
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

#define MAIN_UNIT

#include <errno.h>
#include <float.h>
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "dcmt0.6.1b/include/dc.h"
#include "defs.h"
#include "init.h"
#include "memory.h"
#include "moves.h"
#include "mymath.h"
#include "rand.h"
#include "read.h"
#include "write.h"

int main(int argc, char *argv[]) {
  int opt = 0;
  int long_index = 0;
  static const struct option long_options[] = {
    {"file",required_argument,NULL,1},
    {"seed",required_argument,NULL,2},
    {"threads",required_argument,NULL,3},
    {"length",required_argument,NULL,4},
    {"thin",required_argument,NULL,5},
    {"burnin",required_argument,NULL,6},
    {"npilot",required_argument,NULL,7},
    {"lpilot",required_argument,NULL,8},
    {"outputs",required_argument,NULL,9},
    {"min_M",required_argument,NULL,101},
    {"max_M",required_argument,NULL,102},
    {"max_sig",required_argument,NULL,103},
    {"dlt_cnt",required_argument,NULL,104},
    {"dlt_p",required_argument,NULL,105},
    {"dlt_M",required_argument,NULL,106},
    {"dlt_pi",required_argument,NULL,107},
    {"dlt_sig",required_argument,NULL,108},
    {"dlt_del",required_argument,NULL,109},
    {"dlt_lam",required_argument,NULL,110},
    {"dlt_beta_mu",required_argument,NULL,111},
    {"dlt_beta_nu",required_argument,NULL,112},
    {"fixed_beta",no_argument,NULL,113},
    {"beta_a",required_argument,NULL,114},
    {"beta_b",required_argument,NULL,115},
    {"fixed_lambda",required_argument,NULL,116},
    {"lambda_prior",required_argument,NULL,117},
    {"invgam_shape",required_argument,NULL,118},
    {"invgam_rate",required_argument,NULL,119},
    {"captl_lambda",required_argument,NULL,120},
    {"verbose",no_argument,NULL,121},
    {"calibration",no_argument,NULL,122},
    {"calibration_only",no_argument,NULL,123},
    {"pod_nbr_loci",required_argument,NULL,124},
    {"pool",no_argument,NULL,125},
    {"help",no_argument,NULL,201},
    {"version",no_argument,NULL,202},
    {NULL,0,NULL,0}
  };
  data_struct data;
  data_struct pod;
  updates_struct update,accept,move;
  markov_state_struct current,proposal;
  moments_struct postmean,postmsqr;
  unsigned long seed       = 0;
  int n_threads            = 0;
  unsigned int pod_analysis;
  char *filename           = "data.dat";
  char *path               = "";
  char pod_path[256]       = "";
  char pod_file_name[256]  = "";
  struct stat st;
  double end,end_pilot_runs,start,start_pilot_runs;
  double duration_pilot_runs,elapsed,eta;
  struct tm *local;
  time_t t;
  
  program_name = argv[0];
  beta_a                       = UNDEFINED;
  beta_b                       = UNDEFINED;
  init_delta_beta_mu           = UNDEFINED;
  init_delta_beta_nu           = UNDEFINED;
  init_delta_lambda            = UNDEFINED;
  capital_lambda               = UNDEFINED;
  inverse_gamma_shape          = UNDEFINED;
  inverse_gamma_rate           = UNDEFINED;
  init_delta_pi                = 0.25;
  init_delta_sigma             = 2.5;
  init_delta_delta             = 0.8;
  min_M                        = 0.001;
  max_M                        = 10000;
  max_sigma                    = 700;
  init_delta_counts            = 5;
  init_delta_p                 = 0.25;
  init_delta_M                 = 0.1;
  default_init_delta_beta_mu   = 0.025;
  default_init_delta_beta_nu   = 1.0;
  default_beta_a               = 0.7;
  default_beta_b               = 0.7;
  default_init_delta_lambda    = 0.05;
  default_capital_lambda       = 1.0;
  default_inverse_gamma_shape  = 3.0;
  default_inverse_gamma_rate   = 2.0;
  fixed_beta                   = FALSE;
  inverse_gamma_prior          = TRUE;
  exponential_prior            = FALSE;
  fixed_lambda                 = FALSE;
  verbose                      = FALSE;
  calibration                  = FALSE;
  calibration_only             = FALSE;
  pooled_data                  = FALSE;
  pod_nbr_loci                 = 0;
  n_pilot                      = 25;
  l_pilot                      = 500;
  step                         = 40;
  burn_in                      = 50000;
  chain_length                 = 100000;
  while ((opt = getopt_long_only(argc, argv,"",long_options,&long_index)) != -1) {
    switch (opt) {
      case 1 :
        filename = optarg;
        break;
      case 2 :
        seed = atoi(optarg);
        break;
      case 3 :
        n_threads = atoi(optarg);
        break;
      case 4 :
        chain_length = atol(optarg);
        break;
      case 5 :
        step = atoi(optarg);
        break;
      case 6 :
        burn_in = atol(optarg);
        break;
      case 7 :
        n_pilot = atoi(optarg);
        break;
      case 8 :
        l_pilot = atoi(optarg);
        break;
      case 9 :
        path = optarg;
        break;
      case 101 :
        min_M = atof(optarg);
        break;
      case 102 :
        max_M = atof(optarg);
        break;
      case 103 :
        max_sigma = atof(optarg);
        break;
      case 104 :
        init_delta_counts = atoi(optarg);
        break;
      case 105 :
        init_delta_p = atof(optarg);
        break;
      case 106 :
        init_delta_M = atof(optarg);
        break;
      case 107 :
        init_delta_pi = atof(optarg);
        break;
      case 108 :
        init_delta_sigma = atof(optarg);
        break;
      case 109 :
        init_delta_delta = atof(optarg);
        break;
      case 110 :
        init_delta_lambda = atof(optarg);
        break;
      case 111 :
        init_delta_beta_mu = atof(optarg);
        break;
      case 112 :
        init_delta_beta_nu = atof(optarg);
        break;
      case 113 :
        fixed_beta = TRUE;
        break;
      case 114 :
        beta_a = atof(optarg);
        break;
      case 115 :
        beta_b = atof(optarg);
        break;
      case 116 :
        fixed_lambda = TRUE;
        fixed_lambda_value = atof(optarg);
        break;
      case 117 :
        lambda_prior = optarg;
        break;
      case 118 :
        inverse_gamma_shape = atof(optarg);
        break;
      case 119 :
        inverse_gamma_rate = atof(optarg);
        break;
      case 120 :
        capital_lambda = atof(optarg);
        break;
      case 121 :
        verbose = TRUE;
        break;
      case 122 :
        calibration = TRUE;
        break;
      case 123 :
        calibration_only = TRUE;
        break;
      case 124 :
        pod_nbr_loci = atoi(optarg);
        break;
      case 125 :
        pooled_data = TRUE;
        break;
      case 201 :
        print_usage();
        break;
      case 202 :
        print_version();
        exit(EXIT_SUCCESS);
      default :
        print_usage();
        exit(EXIT_FAILURE);
    }
  }
  if (optind < argc) {
    printf("Error! Unrecognized option(s) ");
    while (optind < argc)
      printf ("`%s' ", argv[optind++]);
    printf ("\n");
    print_usage();
    exit(EXIT_FAILURE);
  }
  if (n_threads < 0) {
    printf("Error! The value of option -threads has to be positive\n");
    exit(EXIT_FAILURE);
  }
  if (n_threads > omp_get_max_threads()) {
    printf("Error! The value of option -threads has to be less than the maximum number of threads available\n");
    exit(EXIT_FAILURE);
  }
  if (step > chain_length) {
    printf("Error! The value of option -step has to be less than total length of the Markov chain (-length)\n");
    exit(EXIT_FAILURE);
  }
  if (min_M < 0) {
    printf("Error! The value of option -min_M has to be positive\n");
    exit(EXIT_FAILURE);
  }
  if (init_delta_counts < 0) {
    printf("Error! The value of option -dlt_cnt has to be positive\n");
    exit(EXIT_FAILURE);
  }
  if (init_delta_p < 0) {
    printf("Error! The value of option -dlt_p has to be positive\n");
    exit(EXIT_FAILURE);
  }
  if (init_delta_p > 2.0) {
    printf("Error! The value of option -dlt_p has to be less than 2.0\n");
    exit(EXIT_FAILURE);
  }
  if (init_delta_M < 0) {
    printf("Error! The value of option -dlt_M has to be positive\n");
    exit(EXIT_FAILURE);
  }
  if (init_delta_pi < 0) {
    printf("Error! The value of option -dlt_pi has to be positive\n");
    exit(EXIT_FAILURE);
  }
  if (init_delta_pi > 2.0) {
    printf("Error! The value of option -dlt_pi has to be less than 2.0\n");
    exit(EXIT_FAILURE);
  }
  if (init_delta_sigma < 0) {
    printf("Error! The value of option -dlt_sig has to be positive\n");
    exit(EXIT_FAILURE);
  }
  if (init_delta_delta < 0) {
    printf("Error! The value of option -dlt_del has to be positive\n");
    exit(EXIT_FAILURE);
  }
  if (fixed_beta) {
    if (!(init_delta_beta_mu == UNDEFINED)) {
      printf("Warning! The value of option -dlt_beta_mu is unnecessary since the shape parameters of the beta prior distribution of pi are fixed\n");
      exit(EXIT_FAILURE);
    }
    if (!(init_delta_beta_nu == UNDEFINED)) {
      printf("Warning! The value of option -dlt_beta_nu is unnecessary since the shape parameters of the beta prior distribution of pi are fixed\n");
      exit(EXIT_FAILURE);
    }
    if (beta_a == UNDEFINED) {
      beta_a = default_beta_a;
    }
    else if (beta_a < 0) {
      printf("Error! The value of option -beta_a has to be positive\n");
      exit(EXIT_FAILURE);
    }
    if (beta_b == UNDEFINED) {
      beta_b = default_beta_b;
    }
    else if (beta_b < 0) {
      printf("Error! The value of option -beta_b has to be positive\n");
      exit(EXIT_FAILURE);
    }
  }
  else {
    if (init_delta_beta_mu == UNDEFINED) {
      init_delta_beta_mu = default_init_delta_beta_mu;
    }
    else if (init_delta_beta_mu < 0) {
      printf("Error! The value of option -dlt_beta_mu has to be positive\n");
      exit(EXIT_FAILURE);
    }
    if (init_delta_beta_nu == UNDEFINED) {
      init_delta_beta_nu = default_init_delta_beta_nu;
    }
    else if (init_delta_beta_nu < 0) {
      printf("Error! The value of option -dlt_beta_nu has to be positive\n");
      exit(EXIT_FAILURE);
    }
    if (!(beta_a == UNDEFINED)) {
      printf("Warning! The value of option -beta_a is unnecessary since the shape parameters of the beta prior distribution of pi are estimated\n");
      exit(EXIT_FAILURE);
    }
    if (!(beta_b == UNDEFINED)) {
      printf("Warning! The value of option -beta_b is unnecessary since the shape parameters of the beta prior distribution of pi are estimated\n");
      exit(EXIT_FAILURE);
    }
  }
  if (fixed_lambda) {
    if (fixed_lambda_value < 0) {
      printf("Error! The value of option -fixed_lambda has to be positive\n");
      exit(EXIT_FAILURE);
    }
    if (!(lambda_prior == NULL)) {
      printf("Warning! The value of option -lambda_prior is unnecessary since lambda is fixed\n");
      exit(EXIT_FAILURE);
    }
    if (!(init_delta_lambda == UNDEFINED)) {
      printf("Warning! The value of option -dlt_lam is unnecessary since lambda is fixed\n");
      exit(EXIT_FAILURE);
    }
    if (!(capital_lambda == UNDEFINED)) {
      printf("Warning! The value of option -captl_lambda is unnecessary since lambda is fixed\n");
      exit(EXIT_FAILURE);
    }
    if (!(inverse_gamma_shape == UNDEFINED)) {
      printf("Warning! The value of option -invgam_shape is unnecessary since lambda is fixed\n");
      exit(EXIT_FAILURE);
    }
    if (!(inverse_gamma_rate == UNDEFINED)) {
      printf("Warning! The value of option -invgam_rate is unnecessary since lambda is fixed\n");
      exit(EXIT_FAILURE);
    }
  }
  else {
    if (!(lambda_prior == NULL)) {
      if (strcmp(lambda_prior,"invgam") == 0) {
        inverse_gamma_prior = TRUE;
        exponential_prior = FALSE;
      }
      else if (strcmp(lambda_prior,"exp") == 0) {
        exponential_prior = TRUE;
        inverse_gamma_prior = FALSE;
      }
      else {
        printf("Error! Option -lambda_prior has to be either 'invgam' (default) or 'exp'\n");
        exit(EXIT_FAILURE);
      }
    }
    if (inverse_gamma_prior) {
      if (!(init_delta_lambda == UNDEFINED)) {
        printf("Warning! The value of option -dlt_lam is unnecessary since an inverse gamma prior is set on lambda (Gibbs sampling)\n");
        exit(EXIT_FAILURE);
      }
      if (!(capital_lambda == UNDEFINED)) {
        printf("Warning! The value of option -captl_lambda is unnecessary since an inverse gamma prior is set on lambda\n");
        exit(EXIT_FAILURE);
      }
      if (inverse_gamma_shape == UNDEFINED) {
        inverse_gamma_shape = default_inverse_gamma_shape;
      }
      else if (inverse_gamma_shape < 0) {
        printf("Error! The value of option -invgam_shape has to be positive\n");
        exit(EXIT_FAILURE);
      }
      if (inverse_gamma_rate == UNDEFINED) {
        inverse_gamma_rate = default_inverse_gamma_rate;
      }
      else if (inverse_gamma_rate < 0) {
        printf("Error! The value of option -invgam_rate has to be positive\n");
        exit(EXIT_FAILURE);
      }
    }
    else if (exponential_prior) {
      if (init_delta_lambda == UNDEFINED) {
        init_delta_lambda = default_init_delta_lambda;
      }
      else if (init_delta_lambda < 0) {
        printf("Error! The value of option -dlt_lam has to be positive\n");
        exit(EXIT_FAILURE);
      }
      if (!(inverse_gamma_shape == UNDEFINED)) {
        printf("Warning! The value of option -invgam_shape is unnecessary since an exponential prior is set on lambda\n");
        exit(EXIT_FAILURE);
      }
      if (!(inverse_gamma_rate == UNDEFINED)) {
        printf("Warning! The value of option -invgam_rate is unnecessary since an exponential prior is set on lambda\n");
        exit(EXIT_FAILURE);
      }
      if (capital_lambda == UNDEFINED) {
        capital_lambda = default_capital_lambda;
      }
      else if (capital_lambda < 0) {
        printf("Error! The value of option -captl_lambda has to be positive\n");
        exit(EXIT_FAILURE);
      }
    }
  }
  if (calibration && calibration_only) {
    printf("Error! Options -calibration and -calibration_only are mutually exclusive\n");
    exit(EXIT_FAILURE);
  }
  if ((pod_nbr_loci > 0) && !(calibration || calibration_only)) {
    printf("Error! Option -pod_nbr_loci should only be used the options -calibration or -calibration_only\n");
    exit(EXIT_FAILURE);
  }
  if (strcmp("",path) != 0) {
    if ((!calibration_only) && (stat(path,&st) == 0)) {
      printf("Warning! Directory %s already exists. You may not want to overwrite previous analyses\n",path);
      exit(EXIT_FAILURE);
    }
    else {
#ifdef __linux__
      mkdir(path,0755);
#elif __APPLE__
      mkdir(path,0755);
#else
      _mkdir(path);
#endif
    }
    if (stat(path,&st) != 0) {
      printf("Cannot create directory %s\n",path);
      exit(EXIT_FAILURE);
    }
  }
  open_log_file(path);
  pod_analysis = FALSE;
  t = time(NULL);
  local = localtime(&t);
  start = omp_get_wtime();
  printf("------------------------------------------------------------------------------------\n");
  fprintf(logfile,"------------------------------------------------------------------------------------\n");
  printf("%s",asctime(local));
  fprintf(logfile,"%s",asctime(local));
  printf("------------------------------------------------------------------------------------\n\n");
  fprintf(logfile,"------------------------------------------------------------------------------------\n\n");
  for (opt = 0; opt < argc; opt++) {
    printf("%s ",argv[opt]);
    fprintf(logfile,"%s ",argv[opt]);
  }
  printf("\n\n");
  fprintf(logfile,"\n\n");
  printf("This analysis was performed using selestim (version %s)\n\n",VERSION);
  fprintf(logfile,"This analysis was performed using selestim (version %s)\n\n",VERSION);
  read_data(&data,filename,pod_analysis);
  if (strcmp("",path) != 0) {
#ifdef __linux__
    if (path[strlen(path) - 1] != '/') {
      strcat(path,"/");
    }
#elif __APPLE__
    if (path[strlen(path) - 1] != '/') {
      strcat(path,"/");
    }
#else
    if (path[strlen(path) - 1] != '\\') {
      strcat(path,"\\");
    }
#endif
  }
  copy_data(data,filename,path);
  if (n_threads == 0) {
    n_threads = omp_get_max_threads();                                          // omp_get_max_threads() returns the same value whether executing from a serial or parallel region
  }
  omp_set_num_threads(n_threads);
  if (!seed) {
    seed = time(NULL);
  }
  mtss = get_mt_parameters_st(32,521,0,(n_threads - 1),4172,&mt_count);         // Initialize the MT random number generator
  if (mtss == NULL) {
    printf("Error! The random number generator could not be initialized\n");
    exit(EXIT_FAILURE);
  }
  int i;
  for (i = 0; i < mt_count; i++) {
    sgenrand_mt((seed + i),mtss[i]);
  }
  total_chain_length = burn_in + chain_length;
  allocate_memory(data,&current,&proposal,&update,&accept,&move,&postmean,&postmsqr);
  if (calibration_only) {
    read_output_files(data,path,&postmean);
  }
  else {
    print_F_ST(data);
    print_model_parameters(seed,n_threads);
    init_moves(data,&move);
    if (n_pilot > 0) {
      start_pilot_runs = omp_get_wtime();
      do_pilot_runs(data,&current,&proposal,&update,&accept,&move);
      end_pilot_runs = omp_get_wtime();
      elapsed = ((double) (end_pilot_runs - start));
      t = time(NULL);
      local = localtime(&t);
      printf("------------------------------------------------------------------------------------\n");
      printf("%s------------------------\n",asctime(local));
      printf("\tComputing time elapsed since beginning = %s\n",print_time(elapsed));
      duration_pilot_runs = end_pilot_runs - start_pilot_runs;
      eta = duration_pilot_runs * (double) total_chain_length / (n_pilot * l_pilot);
      printf("\tEstimated time until the MCMC stops    = %s\n",print_time(eta));
      printf("------------------------------------------------------------------------------------\n\n");
      fflush(stdout);
    }
    init_parameters(&current,data);
    init_moments(data,&postmean,&postmsqr);
    init_updates(data,&update,&accept);
    open_output_files(path);
    write_headers(data);
    run_mcmc(data,&current,&proposal,&postmean,&postmsqr,&update,&accept,&move);
    write_summaries(data,postmean,postmsqr,path,pod_analysis);
    close_output_files();
    print_effective_sample_size(data,postmean,postmsqr,path);
  }
  if (calibration || calibration_only) {
    printf("------------------------------------------------------------------------------------\n");
    fprintf(logfile,"------------------------------------------------------------------------------------\n");
    printf("Calibration of the Kullback-Leibler divergence using pseudo-observed data\n");
    fprintf(logfile,"Calibration of the Kullback-Leibler divergence using pseudo-observed data\n");
    printf("------------------------------------------------------------------------------------\n\n");
    fprintf(logfile,"------------------------------------------------------------------------------------\n\n");
    strcpy(pod_path,path);
    strcat(pod_path,POD_SUBDIRECTORY);
    if (strcmp("",pod_path) != 0) {
      if (stat(pod_path,&st) == 0) {
        printf("\nWarning! Directory %s already exists. You may not want to overwrite previous analyses\n\n",pod_path);
        exit(EXIT_FAILURE);
      }
      else {
#ifdef __linux__
        mkdir(pod_path,0755);
#elif __APPLE__
        mkdir(pod_path,0755);
#else
        _mkdir(pod_path);
#endif
      }
      if (stat(pod_path,&st) != 0) {
        printf("\nCannot create directory %s\n\n",pod_path);
        exit(EXIT_FAILURE);
      }
#ifdef __linux__
      strcat(pod_path,"/");
#elif __APPLE__
      strcat(pod_path,"/");
#else
      strcat(pod_path,"\\");
#endif
    }
    pod_analysis = TRUE;
    open_output_files(pod_path);
    write_pseudo_observed_data(data,postmean,pod_path,filename,pod_file_name);
    write_readme(data,filename,pod_path,pod_file_name);
    release_memory(&data,&current,&proposal,&update,&accept,&move,&postmean,&postmsqr);
    read_data(&pod,pod_file_name,pod_analysis);
    write_headers(pod);
    allocate_memory(pod,&current,&proposal,&update,&accept,&move,&postmean,&postmsqr);
    print_F_ST(pod);
    if (calibration_only) {
      print_model_parameters(seed,n_threads);
    }
    init_moves(pod,&move);
    if (n_pilot > 0) {
      start_pilot_runs = omp_get_wtime();
      do_pilot_runs(pod,&current,&proposal,&update,&accept,&move);
      end_pilot_runs = omp_get_wtime();
      elapsed = end_pilot_runs - start;
      t = time(NULL);
      local = localtime(&t);
      printf("------------------------------------------------------------------------------------\n");
      printf("%s------------------------\n",asctime(local));
      printf("\tComputing time elapsed since beginning = %s\n",print_time(elapsed));
      duration_pilot_runs = end_pilot_runs - start_pilot_runs;
      eta = duration_pilot_runs * (double) total_chain_length / (n_pilot * l_pilot);
      printf("\tEstimated time until calibration       = %s\n",print_time(eta));
      printf("------------------------------------------------------------------------------------\n\n");
      fflush(stdout);
    }
    init_parameters(&current,pod);
    init_moments(pod,&postmean,&postmsqr);
    init_updates(pod,&update,&accept);
    run_mcmc(pod,&current,&proposal,&postmean,&postmsqr,&update,&accept,&move);
    write_summaries(pod,postmean,postmsqr,pod_path,pod_analysis);
    close_output_files();
    release_memory(&pod,&current,&proposal,&update,&accept,&move,&postmean,&postmsqr);
  }
  else {
    release_memory(&data,&current,&proposal,&update,&accept,&move,&postmean,&postmsqr);
  }
  end = omp_get_wtime();
  elapsed = end - start;
  t = time(NULL);
  local = localtime(&t);
  printf("------------------------------------------------------------------------------------\n");
  fprintf(logfile,"------------------------------------------------------------------------------------\n");
  printf("%s------------------------\n",asctime(local));
  fprintf(logfile,"%s------------------------\n",asctime(local));
  printf("\tTotal computing time elapsed           = %s\n",print_time(elapsed));
  fprintf(logfile,"\tTotal computing time elapsed           = %s\n",print_time(elapsed));
  printf("------------------------------------------------------------------------------------\n\n");
  fprintf(logfile,"------------------------------------------------------------------------------------\n\n");
  printf("The program has successfully terminated.\n\n");
  fprintf(logfile,"The program has successfully terminated.\n\n");
  fclose(logfile);
  free_mt_struct_array(mtss,mt_count);
  return(EXIT_SUCCESS);
}

void do_pilot_runs(data_struct data,
                   markov_state_struct *current,
                   markov_state_struct *proposal,
                   updates_struct *update,
                   updates_struct *accept,
                   updates_struct *move)

{
  int i,j;
  unsigned int l,n;
  int missing_data,n_counts,n_p,n_M,n_beta_mu,n_beta_nu,n_pi,n_delta,n_sigma,n_lambda;
  double mean_par,mean_upd,mean_acc,min_par,min_upd,min_acc,max_par,max_upd,max_acc;
  
  n = 0;
  while (n < n_pilot) {
    init_parameters(current,data);
    init_updates(data,update,accept);
    for (l = 0; l < l_pilot; ++l) {
      if (pooled_data) {
        update_counts(current,proposal,(*move),update,accept,data);
      }
      update_p(current,proposal,(*move),update,accept,data);
      update_M(current,proposal,(*move),update,accept,data);
      if (!fixed_beta) {
        update_beta_parameters(current,(*move),update,accept,data);
      }
      update_pi(current,proposal,(*move),update,accept,data);
      if (!fixed_lambda) {
        if (inverse_gamma_prior) {
          update_lambda_gibbs(current,data);
        }
        else if (exponential_prior) {
          update_lambda_metropolis(current,proposal,(*move),update,accept,data);
        }
      }
      update_delta(current,proposal,(*move),update,accept,data);
      update_sigma(current,proposal,(*move),update,accept,data);
      update_kappa(current,proposal,data);
    }
    n_counts = n_p = n_M = n_beta_mu = n_beta_nu = n_pi = n_delta = n_sigma = n_lambda = 0;
    printf("------------------------------------------------------------------------------------\n");
    printf("Pilot run #%2d:\n",(n + 1));
    printf("------------------------------------------------------------------------------------\n");
    fprintf(logfile,"------------------------------------------------------------------------------------\n");
    fprintf(logfile,"Pilot run #%d: \n",(n + 1));
    fprintf(logfile,"------------------------------------------------------------------------------------\n");
    if (pooled_data) {
      printf("\n\tAllele counts x_ij's\n");
      fprintf(logfile,"\n\tAllele counts x_ij's\n");
      for (i = 0; i < data.nbr_demes; ++i) {
        for (j = 0; j < data.nbr_loci; ++j) {
          if (data.total_nbr_counts[i][j] > 0) {
            if (accept -> counts[i][j] / update -> counts[i][j] > MAX_ACCEPT) {
              move -> counts[i][j] = int_min((data.total_nbr_counts[i][j] - 1),move -> counts[i][j] + 1);
              ++n_counts;
            }
            if (accept -> counts[i][j] / update -> counts[i][j] < MIN_ACCEPT) {
              move -> counts[i][j] = int_max(1,move -> counts[i][j] - 1);
              ++n_counts;
            }
          }
        }
      }
      mean_upd = mean_acc = mean_par = 0.0;
      min_upd = min_acc = min_par = 4294967295.0;
      max_upd = max_acc = max_par = -4294967295.0;
      missing_data = 0;
      for (i = 0; i < data.nbr_demes; ++i) {
        for (j = 0; j < data.nbr_loci; ++j) {
          if (update -> counts[i][j] > 0) {
            mean_acc += accept -> counts[i][j] / update -> counts[i][j];
            mean_par += current -> value.counts[i][j][0];
            mean_upd += move -> counts[i][j];
            if (min_acc > accept -> counts[i][j] / update -> counts[i][j]) {min_acc = accept -> counts[i][j] / update -> counts[i][j];}
            if (max_acc < accept -> counts[i][j] / update -> counts[i][j]) {max_acc = accept -> counts[i][j] / update -> counts[i][j];}
            if (min_par > current -> value.counts[i][j][0]) {min_par = current -> value.counts[i][j][0];}
            if (max_par < current -> value.counts[i][j][0]) {max_par = current -> value.counts[i][j][0];}
            if (min_upd > move -> counts[i][j]) {min_upd = move -> counts[i][j];}
            if (max_upd < move -> counts[i][j]) {max_upd = move -> counts[i][j];}
          }
          else {
            missing_data += 1;
          }
        }
      }
      mean_upd /= (data.nbr_loci * data.nbr_demes - missing_data);
      mean_par /= (data.nbr_loci * data.nbr_demes - missing_data);
      mean_acc /= (data.nbr_loci * data.nbr_demes - missing_data);
      printf("\t\taverage value = %6.4f [%6.4f,%6.4f]\n",mean_par,min_par,max_par);
      fprintf(logfile,"\t\taverage value = %6.4f [%6.4f,%6.4f]\n",mean_par,min_par,max_par);
      printf("\t\taverage updating parameter = %6.4f [%6.4f,%6.4f]\n",mean_upd,min_upd,max_upd);
      fprintf(logfile,"\t\taverage updating parameter = %6.4f [%6.4f,%6.4f]\n",mean_upd,min_upd,max_upd);
      printf("\t\taverage acceptance rate = %6.4f [%6.4f,%6.4f]\n",mean_acc,min_acc,max_acc);
      fprintf(logfile,"\t\taverage acceptance rate = %6.4f [%6.4f,%6.4f]\n",mean_acc,min_acc,max_acc);
      printf("\t\t%d parameters have been scaled, out of %d\n",n_counts,(data.nbr_loci * data.nbr_demes));
      fprintf(logfile,"\t\t%d parameters have been scaled, out of %d\n",n_counts,(data.nbr_loci * data.nbr_demes));
    }
    printf("\n\tAllele frequencies p_ij's\n");
    fprintf(logfile,"\n\tAllele frequencies p_ij's\n");
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        if (accept -> p[i][j] / update -> p[i][j] > MAX_ACCEPT) {
          if (move -> p[i][j] < (1.0 - 2.0 * MIN_P)) {
            if ((move -> p[i][j] *= SCALING) > (1.0 - 2.0 * MIN_P)) move -> p[i][j] = (1.0 - 2.0 * MIN_P);
            ++n_p;
          }
        }
        if (accept -> p[i][j] / update -> p[i][j] < MIN_ACCEPT) {
          move -> p[i][j] /= SCALING;
          ++n_p;
        }
      }
    }
    mean_upd = mean_acc = mean_par = 0.0;
    min_upd = min_acc = min_par = 4294967295.0;
    max_upd = max_acc = max_par = -4294967295.0;
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        mean_par += current -> value.p[i][j][0];
        mean_upd += move -> p[i][j];
        mean_acc += accept -> p[i][j] / update -> p[i][j];
        if (min_upd > move -> p[i][j]) {min_upd = move -> p[i][j];}
        if (max_upd < move -> p[i][j]) {max_upd = move -> p[i][j];}
        if (min_par > current -> value.p[i][j][0]) {min_par = current -> value.p[i][j][0];}
        if (max_par < current -> value.p[i][j][0]) {max_par = current -> value.p[i][j][0];}
        if (min_acc > accept -> p[i][j] / update -> p[i][j]) {min_acc = accept -> p[i][j] / update -> p[i][j];}
        if (max_acc < accept -> p[i][j] / update -> p[i][j]) {max_acc = accept -> p[i][j] / update -> p[i][j];}
      }
    }
    mean_upd /= (data.nbr_loci * data.nbr_demes);
    mean_par /= (data.nbr_loci * data.nbr_demes);
    mean_acc /= (data.nbr_loci * data.nbr_demes);
    printf("\t\taverage value = %6.4f [%6.4f,%6.4f]\n",mean_par,min_par,max_par);
    fprintf(logfile,"\t\taverage value = %6.4f [%6.4f,%6.4f]\n",mean_par,min_par,max_par);
    printf("\t\taverage updating parameter = %6.4f [%6.4f,%6.4f]\n",mean_upd,min_upd,max_upd);
    fprintf(logfile,"\t\taverage updating parameter = %6.4f [%6.4f,%6.4f]\n",mean_upd,min_upd,max_upd);
    printf("\t\taverage acceptance rate = %6.4f [%6.4f,%6.4f]\n",mean_acc,min_acc,max_acc);
    fprintf(logfile,"\t\taverage acceptance rate = %6.4f [%6.4f,%6.4f]\n",mean_acc,min_acc,max_acc);
    printf("\t\t%d parameters have been scaled, out of %d\n",n_p,(data.nbr_loci * data.nbr_demes));
    fprintf(logfile,"\t\t%d parameters have been scaled, out of %d\n",n_p,(data.nbr_loci * data.nbr_demes));
    printf("\n\tPopulation parameters M_i's\n");
    fprintf(logfile,"\n\tPopulation parameters M_i's\n");
    for (i = 0; i < data.nbr_demes; ++i) {
      if (accept -> M[i] / update -> M[i] > MAX_ACCEPT) {
        move -> M[i] *= SCALING;
        ++n_M;
      }
      if (accept -> M[i] / update -> M[i] < MIN_ACCEPT) {
        move -> M[i] /= SCALING;
        ++n_M;
      }
    }
    mean_upd = mean_acc = mean_par = 0.0;
    min_upd = min_acc = min_par = 4294967295.0;
    max_upd = max_acc = max_par = -4294967295.0;
    for (i = 0; i < data.nbr_demes; ++i) {
      mean_par += current -> value.M[i];
      mean_upd += move -> M[i];
      mean_acc += accept -> M[i] / update -> M[i];
      if (min_upd > move -> M[i]) {min_upd = move -> M[i];}
      if (max_upd < move -> M[i]) {max_upd = move -> M[i];}
      if (min_par > current -> value.M[i]) {min_par = current -> value.M[i];}
      if (max_par < current -> value.M[i]) {max_par = current -> value.M[i];}
      if (min_acc > accept -> M[i] / update -> M[i]) {min_acc = accept -> M[i] / update -> M[i];}
      if (max_acc < accept -> M[i] / update -> M[i]) {max_acc = accept -> M[i] / update -> M[i];}
    }
    mean_upd /= data.nbr_demes;
    mean_par /= data.nbr_demes;
    mean_acc /= data.nbr_demes;
    printf("\t\taverage value = %6.4f [%6.4f,%6.4f]\n",mean_par,min_par,max_par);
    fprintf(logfile,"\t\taverage value = %6.4f [%6.4f,%6.4f]\n",mean_par,min_par,max_par);
    printf("\t\taverage updating parameter = %6.4f [%6.4f,%6.4f]\n",mean_upd,min_upd,max_upd);
    fprintf(logfile,"\t\taverage updating parameter = %6.4f [%6.4f,%6.4f]\n",mean_upd,min_upd,max_upd);
    printf("\t\taverage acceptance rate = %6.4f [%6.4f,%6.4f]\n",mean_acc,min_acc,max_acc);
    fprintf(logfile,"\t\taverage acceptance rate = %6.4f [%6.4f,%6.4f]\n",mean_acc,min_acc,max_acc);
    printf("\t\t%d parameters have been scaled, out of %d\n",n_M,data.nbr_demes);
    fprintf(logfile,"\t\t%d parameters have been scaled, out of %d\n",n_M,data.nbr_demes);
    if (!fixed_beta) {
      printf("\n\tShape parameter (a) of the prior distribution of migrant allele frequencies pi_j's\n");
      fprintf(logfile,"\n\tShape parameter (a) of the prior distribution of migrant allele frequencies pi_j's\n");
      printf("\t\tcurrent value = %6.4f\n",current -> value.beta_a);
      fprintf(logfile,"\t\tcurrent value = %6.4f\n",current -> value.beta_a);
      if (accept -> beta_mu / update -> beta_mu > MAX_ACCEPT) {
        move -> beta_mu *= SCALING;
        ++n_beta_mu;
      }
      if (accept -> beta_mu / update -> beta_mu < MIN_ACCEPT) {
        move -> beta_mu /= SCALING;
        ++n_beta_mu;
      }
      printf("\t\tupdating parameter = %6.4f\n",move -> beta_mu);
      fprintf(logfile,"\t\tupdating parameter = %6.4f\n",move -> beta_mu);
      printf("\t\taverage acceptance rate = %6.4f\n",((double) accept -> beta_mu / update -> beta_mu));
      fprintf(logfile,"\t\taverage acceptance rate = %6.4f\n",((double) accept -> beta_mu / update -> beta_mu));
      printf("\t\t%d parameters have been scaled, out of 1\n",n_beta_mu);
      fprintf(logfile,"\t\t%d parameters have been scaled, out of 1\n",n_beta_mu);
      printf("\n\tShape parameter (b) of the prior distribution of migrant allele frequencies pi_j's\n");
      fprintf(logfile,"\n\tShape parameter (b) of the prior distribution of migrant allele frequencies pi_j's\n");
      printf("\t\tcurrent value = %6.4f\n",current -> value.beta_b);
      fprintf(logfile,"\t\tcurrent value = %6.4f\n",current -> value.beta_b);
      if (accept -> beta_nu / update -> beta_nu > MAX_ACCEPT) {
        move -> beta_nu *= SCALING;
        ++n_beta_nu;
      }
      if (accept -> beta_nu / update -> beta_nu < MIN_ACCEPT) {
        move -> beta_nu /= SCALING;
        ++n_beta_nu;
      }
      printf("\t\tupdating parameter = %6.4f\n",move -> beta_nu);
      fprintf(logfile,"\t\tupdating parameter = %6.4f\n",move -> beta_nu);
      printf("\t\taverage acceptance rate = %6.4f\n",((double) accept -> beta_nu / update -> beta_nu));
      fprintf(logfile,"\t\taverage acceptance rate = %6.4f\n",((double) accept -> beta_nu / update -> beta_nu));
      printf("\t\t%d parameters have been scaled, out of 1\n",n_beta_nu);
      fprintf(logfile,"\t\t%d parameters have been scaled, out of 1\n",n_beta_nu);
    }
    printf("\n\tMigrant allele frequencies pi_j's\n");
    fprintf(logfile,"\n\tMigrant allele frequencies pi_j's\n");
    for (j = 0; j < data.nbr_loci; ++j) {
      if (accept -> pi[j] / update -> pi[j] > MAX_ACCEPT) {
        if (move -> pi[j] < (1 - 2 * MIN_PI)) {
          if ((move -> pi[j] *= SCALING) > (1 - 2 * MIN_PI)) move -> pi[j] = (1 - 2 * MIN_PI);
          ++n_pi;
        }
      }
      if (accept -> pi[j] / update -> pi[j] < MIN_ACCEPT) {
        move -> pi[j] /= SCALING;
        ++n_pi;
      }
    }
    mean_upd = mean_acc = mean_par = 0.0;
    min_upd = min_acc = min_par = 4294967295.0;
    max_upd = max_acc = max_par = -4294967295.0;
    for (j = 0; j < data.nbr_loci; ++j) {
      mean_par += current -> value.pi[j][0];
      mean_upd += move -> pi[j];
      mean_acc += accept -> pi[j] / update -> pi[j];
      if (min_upd > move -> pi[j]) {min_upd = move -> pi[j];}
      if (max_upd < move -> pi[j]) {max_upd = move -> pi[j];}
      if (min_par > current -> value.pi[j][0]) {min_par = current -> value.pi[j][0];}
      if (max_par < current -> value.pi[j][0]) {max_par = current -> value.pi[j][0];}
      if (min_acc > accept -> pi[j] / update -> pi[j]) {min_acc = accept -> pi[j] / update -> pi[j];}
      if (max_acc < accept -> pi[j] / update -> pi[j]) {max_acc = accept -> pi[j] / update -> pi[j];}
    }
    mean_upd /= data.nbr_loci;
    mean_par /= data.nbr_loci;
    mean_acc /= data.nbr_loci;
    printf("\t\taverage value = %6.4f [%6.4f,%6.4f]\n",mean_par,min_par,max_par);
    fprintf(logfile,"\t\taverage value = %6.4f [%6.4f,%6.4f]\n",mean_par,min_par,max_par);
    printf("\t\taverage updating parameter = %6.4f [%6.4f,%6.4f]\n",mean_upd,min_upd,max_upd);
    fprintf(logfile,"\t\taverage updating parameter = %6.4f [%6.4f,%6.4f]\n",mean_upd,min_upd,max_upd);
    printf("\t\taverage acceptance rate = %6.4f [%6.4f,%6.4f]\n",mean_acc,min_acc,max_acc);
    fprintf(logfile,"\t\taverage acceptance rate = %6.4f [%6.4f,%6.4f]\n",mean_acc,min_acc,max_acc);
    printf("\t\t%d parameters have been scaled, out of %d\n",n_pi,data.nbr_loci);
    fprintf(logfile,"\t\t%d parameters have been scaled, out of %d\n",n_pi,data.nbr_loci);
    if (!fixed_lambda) {
      printf("\n\tGenome-wide coefficient of selection lambda\n");
      fprintf(logfile,"\n\tGenome-wide coefficient of selection lambda\n");
      printf("\t\tcurrent value = %6.4f\n",current -> value.lambda);
      fprintf(logfile,"\t\tcurrent value = %6.4f\n",current -> value.lambda);
      if (exponential_prior) {
        if (accept -> lambda / update -> lambda > MAX_ACCEPT) {
          move -> lambda *= SCALING;
          ++n_lambda;
        }
        if (accept -> lambda / update -> lambda < MIN_ACCEPT) {
          move -> lambda /= SCALING;
          ++n_lambda;
        }
        printf("\t\tupdating parameter = %6.4f\n",move -> lambda);
        fprintf(logfile,"\t\tupdating parameter = %6.4f\n",move -> lambda);
        printf("\t\taverage acceptance rate = %6.4f\n",((double) accept -> lambda / update -> lambda));
        fprintf(logfile,"\t\taverage acceptance rate = %6.4f\n",((double) accept -> lambda / update -> lambda));
        printf("\t\t%d parameters have been scaled, out of 1\n",n_lambda);
        fprintf(logfile,"\t\t%d parameters have been scaled, out of 1\n",n_lambda);
      }
    }
    printf("\n\tLocus-specific selection coefficient delta_j's\n");
    fprintf(logfile,"\n\tLocus-specific selection coefficient delta_j's\n");
    for (j = 0; j < data.nbr_loci; ++j) {
      if (accept -> delta[j] / update -> delta[j] > MAX_ACCEPT) {
        move -> delta[j] *= SCALING;
        ++n_delta;
      }
      if (accept -> delta[j] / update -> delta[j] < MIN_ACCEPT) {
        move -> delta[j] /= SCALING;
        ++n_delta;
      }
    }
    mean_upd = mean_acc = mean_par = 0.0;
    min_upd = min_acc = min_par = 4294967295.0;
    max_upd = max_acc = max_par = -4294967295.0;
    for (j = 0; j < data.nbr_loci; ++j) {
      mean_par += current -> value.delta[j];
      mean_upd += move -> delta[j];
      mean_acc += accept -> delta[j] / update -> delta[j];
      if (min_upd > move -> delta[j]) {min_upd = move -> delta[j];}
      if (max_upd < move -> delta[j]) {max_upd = move -> delta[j];}
      if (min_par > current -> value.delta[j]) {min_par = current -> value.delta[j];}
      if (max_par < current -> value.delta[j]) {max_par = current -> value.delta[j];}
      if (min_acc > accept -> delta[j] / update -> delta[j]) {min_acc = accept -> delta[j] / update -> delta[j];}
      if (max_acc < accept -> delta[j] / update -> delta[j]) {max_acc = accept -> delta[j] / update -> delta[j];}
    }
    mean_upd /= data.nbr_loci;
    mean_par /= data.nbr_loci;
    mean_acc /= data.nbr_loci;
    printf("\t\taverage value = %6.4f [%6.4f,%6.4f]\n",mean_par,min_par,max_par);
    fprintf(logfile,"\t\taverage value = %6.4f [%6.4f,%6.4f]\n",mean_par,min_par,max_par);
    printf("\t\taverage updating parameter = %6.4f [%6.4f,%6.4f]\n",mean_upd,min_upd,max_upd);
    fprintf(logfile,"\t\taverage updating parameter = %6.4f [%6.4f,%6.4f]\n",mean_upd,min_upd,max_upd);
    printf("\t\taverage acceptance rate = %6.4f [%6.4f,%6.4f]\n",mean_acc,min_acc,max_acc);
    fprintf(logfile,"\t\taverage acceptance rate = %6.4f [%6.4f,%6.4f]\n",mean_acc,min_acc,max_acc);
    printf("\t\t%d parameters have been scaled, out of %d\n",n_delta,data.nbr_loci);
    fprintf(logfile,"\t\t%d parameters have been scaled, out of %d\n",n_delta,data.nbr_loci);
    printf("\n\tLocus- population-specific selection coefficient sigma_ij's\n");
    fprintf(logfile,"\n\tLocus- population-specific selection coefficient sigma_ij's\n");
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        if (accept -> sigma[i][j] / update -> sigma[i][j] > MAX_ACCEPT) {
          move -> sigma[i][j] *= SCALING;
          ++n_sigma;
        }
        if (accept -> sigma[i][j] / update -> sigma[i][j] < MIN_ACCEPT) {
          move -> sigma[i][j] /= SCALING;
          ++n_sigma;
        }
      }
    }
    mean_upd = mean_acc = mean_par = 0.0;
    min_upd = min_acc = min_par = 4294967295.0;
    max_upd = max_acc = max_par = -4294967295.0;
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        mean_par += current -> value.sigma[i][j];
        mean_upd += move -> sigma[i][j];
        mean_acc += accept -> sigma[i][j] / update -> sigma[i][j];
        if (min_upd > move -> sigma[i][j]) {min_upd = move -> sigma[i][j];}
        if (max_upd < move -> sigma[i][j]) {max_upd = move -> sigma[i][j];}
        if (min_par > current -> value.sigma[i][j]) {min_par = current -> value.sigma[i][j];}
        if (max_par < current -> value.sigma[i][j]) {max_par = current -> value.sigma[i][j];}
        if (min_acc > accept -> sigma[i][j] / update -> sigma[i][j]) {min_acc = accept -> sigma[i][j] / update -> sigma[i][j];}
        if (max_acc < accept -> sigma[i][j] / update -> sigma[i][j]) {max_acc = accept -> sigma[i][j] / update -> sigma[i][j];}
      }
    }
    mean_upd /= (data.nbr_loci * data.nbr_demes);
    mean_par /= (data.nbr_loci * data.nbr_demes);
    mean_acc /= (data.nbr_loci * data.nbr_demes);
    printf("\t\taverage value = %6.4f [%6.4f,%6.4f]\n",mean_par,min_par,max_par);
    fprintf(logfile,"\t\taverage value = %6.4f [%6.4f,%6.4f]\n",mean_par,min_par,max_par);
    printf("\t\taverage updating parameter = %6.4f [%6.4f,%6.4f]\n",mean_upd,min_upd,max_upd);
    fprintf(logfile,"\t\taverage updating parameter = %6.4f [%6.4f,%6.4f]\n",mean_upd,min_upd,max_upd);
    printf("\t\taverage acceptance rate = %6.4f [%6.4f,%6.4f]\n",mean_acc,min_acc,max_acc);
    fprintf(logfile,"\t\taverage acceptance rate = %6.4f [%6.4f,%6.4f]\n",mean_acc,min_acc,max_acc);
    printf("\t\t%d parameters have been scaled, out of %d\n\n",n_sigma,(data.nbr_loci * data.nbr_demes));
    fprintf(logfile,"\t\t%d parameters have been scaled, out of %d\n\n",n_sigma,(data.nbr_loci * data.nbr_demes));
    fflush(stdout);
    fflush(logfile);
    ++n;
  }
}

void run_mcmc(data_struct data,
              markov_state_struct *current,
              markov_state_struct *proposal,
              moments_struct *postmean,
              moments_struct *postmsqr,
              updates_struct *update,
              updates_struct *accept,
              updates_struct *move)

{
  int i,j;
  int cnt,tick10,tick100;
  unsigned long iter,trick;
  
  iter = 0;
  trick = (int) (total_chain_length / 100) * 100;
  tick10 = (int) (trick / 10);
  tick100 = (int) (trick / 100);
  cnt = 0;
  printf("------------------------------------------------------------------------------------\n");
  printf("Running the MCMC\n");
  printf("------------------------------------------------------------------------------------\n");
  while (iter < total_chain_length) {
    if (iter < trick) {
      if ((iter / tick10) * tick10 == iter) {
        if (cnt == 0) {
          printf("\n starting [");
        }
        else {
          printf("]\n%3d%% done [",cnt);
        }
        cnt += 10;
      }
      if ((iter / tick100) * tick100 == iter) {
        printf(".");
        fflush(stdout);
      }
    }
    iter += 1;
    if(pooled_data) {
      update_counts(current,proposal,(*move),update,accept,data);
    }
    update_p(current,proposal,(*move),update,accept,data);
    update_M(current,proposal,(*move),update,accept,data);
    if (!fixed_beta) {
      update_beta_parameters(current,(*move),update,accept,data);
    }
    update_pi(current,proposal,(*move),update,accept,data);
    if (!fixed_lambda) {
      if (inverse_gamma_prior) {
        update_lambda_gibbs(current,data);
      }
      else if (exponential_prior) {
        update_lambda_metropolis(current,proposal,(*move),update,accept,data);
      }
      else {
        exit(EXIT_FAILURE);
      }
    }
    update_delta(current,proposal,(*move),update,accept,data);
    update_sigma(current,proposal,(*move),update,accept,data);
    update_kappa(current,proposal,data);
    if (iter > burn_in) {
      if ((iter / step) * step == iter) {
        for (i = 0; i < data.nbr_demes; ++i) {
          postmean -> M[i] += current -> value.M[i];
          postmsqr -> M[i] += current -> value.M[i] * current -> value.M[i];
        }
        for (j = 0; j < data.nbr_loci; ++j) {
          postmean -> pi[j] += current -> value.pi[j][0];
          postmsqr -> pi[j] += current -> value.pi[j][0] * current -> value.pi[j][0];
          postmean -> delta[j] += current -> value.delta[j];
          postmsqr -> delta[j] += current -> value.delta[j] * current -> value.delta[j];
        }
        for (i = 0; i < data.nbr_demes; ++i) {
          for (j = 0; j < data.nbr_loci; ++j) {
            postmean -> p[i][j] += current -> value.p[i][j][0];
            postmsqr -> p[i][j] += current -> value.p[i][j][0] * current -> value.p[i][j][0];
            postmean -> kappa[i][j] += current -> value.kappa[i][j];
            postmean -> sigma[i][j][current -> value.kappa[i][j]] += current -> value.sigma[i][j];
            postmsqr -> sigma[i][j][current -> value.kappa[i][j]] += current -> value.sigma[i][j] * current -> value.sigma[i][j];
            postmean -> sigma[i][j][2] += current -> value.sigma[i][j];
            postmsqr -> sigma[i][j][2] += current -> value.sigma[i][j] * current -> value.sigma[i][j];
          }
        }
        postmean -> lambda += current -> value.lambda;
        postmsqr -> lambda += current -> value.lambda * current -> value.lambda;
        if (pooled_data) {
          for (i = 0; i < data.nbr_demes; ++i) {
            for (j = 0; j < data.nbr_loci; ++j) {
              postmean -> counts[i][j] += (double) current -> value.counts[i][j][0];
              postmsqr -> counts[i][j] += (double) current -> value.counts[i][j][0] * current -> value.counts[i][j][0];
            }
          }
        }
        if (!fixed_beta) {
          postmean -> beta_a += current -> value.beta_a;
          postmsqr -> beta_a += current -> value.beta_a * current -> value.beta_a;
          postmean -> beta_b += current -> value.beta_b;
          postmsqr -> beta_b += current -> value.beta_b * current -> value.beta_b;
        }
        write_outputs(data,(*current),(*update),(*accept),iter);
        for (i = 0; i < data.nbr_demes; ++i) {
          update -> M[i] = accept -> M[i] = 0;
        }
        for (j = 0; j < data.nbr_loci; ++j) {
          update -> pi[j] = accept -> pi[j] = 0;
          update -> delta[j] = accept -> delta[j] = 0;
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
    }
  }
  for (i = 0; i < data.nbr_demes; ++i) {
    postmean -> M[i] /= floor(chain_length / step);
    postmsqr -> M[i] /= floor(chain_length / step);
  }
  for (j = 0; j < data.nbr_loci; ++j) {
    postmean -> pi[j] /= floor(chain_length / step);
    postmsqr -> pi[j] /= floor(chain_length / step);
    postmean -> delta[j] /= floor(chain_length / step);
    postmsqr -> delta[j] /= floor(chain_length / step);
  }
  for (i = 0; i < data.nbr_demes; ++i) {
    for (j = 0; j < data.nbr_loci; ++j) {
      postmean -> p[i][j] /= floor(chain_length / step);
      postmsqr -> p[i][j] /= floor(chain_length / step);
      postmean -> sigma[i][j][0] /= floor(chain_length / step) - postmean -> kappa[i][j];
      postmsqr -> sigma[i][j][0] /= floor(chain_length / step) - postmean -> kappa[i][j];
      postmean -> sigma[i][j][1] /= postmean -> kappa[i][j];
      postmsqr -> sigma[i][j][1] /= postmean -> kappa[i][j];
      postmean -> sigma[i][j][2] /= floor(chain_length / step);
      postmsqr -> sigma[i][j][2] /= floor(chain_length / step);
    }
  }
  for (i = 0; i < data.nbr_demes; ++i) {
    for (j = 0; j < data.nbr_loci; ++j) {
      postmean -> kappa[i][j] /= floor(chain_length / step);
    }
  }
  postmean -> lambda /= floor(chain_length / step);
  postmsqr -> lambda /= floor(chain_length / step);
  if (!fixed_beta) {
    postmean -> beta_a /= floor(chain_length / step);
    postmsqr -> beta_a /= floor(chain_length / step);
    postmean -> beta_b /= floor(chain_length / step);
    postmsqr -> beta_b /= floor(chain_length / step);
  }
  if (pooled_data) {
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        postmean -> counts[i][j] /= floor(chain_length / step);
        postmsqr -> counts[i][j] /= floor(chain_length / step);
      }
    }
  }
  printf("]\n%3d%% done !\n\n",cnt);
}
