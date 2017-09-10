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

#include "write.h"

void open_log_file(char *path)

{
  char filename[256] = "";
  
  strcpy(filename,path);
  if (strcmp("",path) != 0) {
#ifdef __linux__
    strcat(filename,"/");
#elif __APPLE__
    strcat(filename,"/");
#else
    strcat(filename,"\\");
#endif
  }
  strcat(filename,LOGFILE);
  if (calibration_only) {
    if ((logfile = fopen(filename,"a")) == NULL) {
      printf("Can't open %s for output!\n",filename);
      exit(EXIT_FAILURE);
    }
  }
  else {
    if ((logfile = fopen(filename,"w")) == NULL) {
      printf("Can't open %s for output!\n",filename);
      exit(EXIT_FAILURE);
    }
  }
}

void open_output_files(char *path)

{
  char filename[256] = "";

  strcpy(filename,path);
  strcat(filename,OUTPUTS_MCMC);
  if ((outputs_mcmc = fopen(filename,"w")) == NULL) {
    printf("Can't open %s for output!\n",filename);
    exit(EXIT_FAILURE);
  }
  strcpy(filename,path);
  strcat(filename,TRACE_M);
  if ((trace_m = fopen(filename,"w")) == NULL) {
    printf("Can't open %s for output!\n",filename);
    exit(EXIT_FAILURE);
  }
  if (!fixed_lambda) {
    strcpy(filename,path);
    strcat(filename,TRACE_LAMBDA);
    if ((trace_lambda = fopen(filename,"w")) == NULL) {
      printf("Can't open %s for output!\n",filename);
      exit(EXIT_FAILURE);
    }
  }
  if (!fixed_beta) {
    strcpy(filename,path);
    strcat(filename,TRACE_BETA);
    if ((trace_beta = fopen(filename,"w")) == NULL) {
      printf("Can't open %s for output!\n",filename);
      exit(EXIT_FAILURE);
    }
  }
  if (verbose) {
    strcpy(filename,path);
    strcat(filename,TRACE_FREQ);
    if ((trace_freq = fopen(filename,"w")) == NULL) {
      printf("Can't open %s for output!\n",filename);
      exit(EXIT_FAILURE);
    }
    if (pooled_data) {
      strcpy(filename,path);
      strcat(filename,TRACE_COUNTS);
      if ((trace_counts = fopen(filename,"w")) == NULL) {
        printf("Can't open %s for output!\n",filename);
        exit(EXIT_FAILURE);
      }
    }
    strcpy(filename,path);
    strcat(filename,TRACE_PI);
    if ((trace_pi = fopen(filename,"w")) == NULL) {
      printf("Can't open %s for output!\n",filename);
      exit(EXIT_FAILURE);
    }
    strcpy(filename,path);
    strcat(filename,TRACE_KAPPA);
    if ((trace_kappa = fopen(filename,"w")) == NULL) {
      printf("Can't open %s for output!\n",filename);
      exit(EXIT_FAILURE);
    }
    strcpy(filename,path);
    strcat(filename,TRACE_SIGMA);
    if ((trace_sigma = fopen(filename,"w")) == NULL) {
      printf("Can't open %s for output!\n",filename);
      exit(EXIT_FAILURE);
    }
    strcpy(filename,path);
    strcat(filename,TRACE_DELTA);
    if ((trace_delta = fopen(filename,"w")) == NULL) {
      printf("Can't open %s for output!\n",filename);
      exit(EXIT_FAILURE);
    }
  }
}

void write_headers(data_struct data)

{
  int i,j,k;
  int nbr_spaces;
  
  fprintf(outputs_mcmc," iteration        log_Lik");
  if (pooled_data) {
    fprintf(outputs_mcmc,"      acpt_n");
  }
  fprintf(outputs_mcmc,"      acpt_p");
  fprintf(outputs_mcmc,"      acpt_M     acpt_pi");
  if (!fixed_beta) {
    fprintf(outputs_mcmc,"     acpt_mu     acpt_nu");
  }
  if (!fixed_lambda) {
    if (exponential_prior) {
      fprintf(outputs_mcmc," acpt_lambda");
    }
  }
  fprintf(outputs_mcmc,"  acpt_delta  acpt_sigma\n");
  fprintf(trace_m," iteration");
  for (i = 0; i < data.nbr_demes; ++i) {
    nbr_spaces = 15 - 5 - (floor(log10(i + 1)) + 1);
    for (k = 0; k < nbr_spaces; ++k) {
      fprintf(trace_m," ");
    }
    fprintf(trace_m,"deme_%d",(i + 1));
  }
  fprintf(trace_m,"\n");
  if (!fixed_lambda) {
    fprintf(trace_lambda," iteration      lambda\n");
  }
  if (!fixed_beta) {
    fprintf(trace_beta," iteration      beta_a      beta_b\n");
  }
  if (verbose) {
    fprintf(trace_freq,"         -");
    for (i = 0; i < data.nbr_demes; ++i) {
      nbr_spaces = 14 - 5 - (floor(log10(i + 1)) + 1);
      for (j = 0; j < data.nbr_loci; ++j) {
        for (k = 0; k < nbr_spaces; ++k) {
          fprintf(trace_freq," ");
        }
        fprintf(trace_freq,"deme_%d",(i + 1));
      }
    }
    fprintf(trace_freq,"\n");
    fprintf(trace_freq," iteration");
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        nbr_spaces = 14 - 6 - (floor(log10(j + 1)) + 1);
        for (k = 0; k < nbr_spaces; ++k) {
          fprintf(trace_freq," ");
        }
        fprintf(trace_freq,"locus_%d",(j + 1));
      }
    }
    fprintf(trace_freq,"\n");
    if (pooled_data) {
      fprintf(trace_counts,"         -");
      for (i = 0; i < data.nbr_demes; ++i) {
        nbr_spaces = 14 - 5 - (floor(log10(i + 1)) + 1);
        for (j = 0; j < data.nbr_loci; ++j) {
          for (k = 0; k < nbr_spaces; ++k) {
            fprintf(trace_counts," ");
          }
          fprintf(trace_counts,"deme_%d",(i + 1));
        }
      }
      fprintf(trace_counts,"\n");
      fprintf(trace_counts," iteration");
      for (i = 0; i < data.nbr_demes; ++i) {
        for (j = 0; j < data.nbr_loci; ++j) {
          nbr_spaces = 14 - 6 - (floor(log10(j + 1)) + 1);
          for (k = 0; k < nbr_spaces; ++k) {
            fprintf(trace_counts," ");
          }
          fprintf(trace_counts,"locus_%d",(j + 1));
        }
      }
      fprintf(trace_counts,"\n");
    }
    fprintf(trace_pi," iteration");
    for (j = 0; j < data.nbr_loci; ++j) {
      nbr_spaces = 14 - 6 - (floor(log10(j + 1)) + 1);
      for (k = 0; k < nbr_spaces; ++k) {
        fprintf(trace_pi," ");
      }
      fprintf(trace_pi,"locus_%d",(j + 1));
    }
    fprintf(trace_pi,"\n");
    fprintf(trace_kappa,"         -");
    for (i = 0; i < data.nbr_demes; ++i) {
      nbr_spaces = 14 - 5 - (floor(log10(i + 1)) + 1);
      for (j = 0; j < data.nbr_loci; ++j) {
        for (k = 0; k < nbr_spaces; ++k) {
          fprintf(trace_kappa," ");
        }
        fprintf(trace_kappa,"deme_%d",(i + 1));
      }
    }
    fprintf(trace_kappa,"\n");
    fprintf(trace_kappa," iteration");
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        nbr_spaces = 14 - 6 - (floor(log10(j + 1)) + 1);
        for (k = 0; k < nbr_spaces; ++k) {
          fprintf(trace_kappa," ");
        }
        fprintf(trace_kappa,"locus_%d",(j + 1));
      }
    }
    fprintf(trace_kappa,"\n");
    fprintf(trace_sigma,"         -");
    for (i = 0; i < data.nbr_demes; ++i) {
      nbr_spaces = 14 - 5 - (floor(log10(i + 1)) + 1);
      for (j = 0; j < data.nbr_loci; ++j) {
        for (k = 0; k < nbr_spaces; ++k) {
          fprintf(trace_sigma," ");
        }
        fprintf(trace_sigma,"deme_%d",(i + 1));
      }
    }
    fprintf(trace_sigma,"\n");
    fprintf(trace_sigma," iteration");
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        nbr_spaces = 14 - 6 - (floor(log10(j + 1)) + 1);
        for (k = 0; k < nbr_spaces; ++k) {
          fprintf(trace_sigma," ");
        }
        fprintf(trace_sigma,"locus_%d",(j + 1));
      }
    }
    fprintf(trace_sigma,"\n");
    fprintf(trace_delta," iteration");
    for (j = 0; j < data.nbr_loci; ++j) {
      nbr_spaces = 14 - 6 - (floor(log10(j + 1)) + 1);
      for (k = 0; k < nbr_spaces; ++k) {
        fprintf(trace_delta," ");
      }
      fprintf(trace_delta,"locus_%d",(j + 1));
    }
    fprintf(trace_delta,"\n");
  }
}

void close_output_files(void)

{
	fclose(outputs_mcmc);
	fclose(trace_m);
  if (verbose) {
    fclose(trace_freq);
    if (pooled_data) {
      fclose(trace_counts);
    }
    fclose(trace_pi);
    fclose(trace_kappa);
    fclose(trace_sigma);
    fclose(trace_delta);
  }
  if (!fixed_lambda) {
    fclose(trace_lambda);
  }
  if (!fixed_beta) {
    fclose(trace_beta);
  }
}

void write_outputs(data_struct data,
                   markov_state_struct current,
                   updates_struct update,
                   updates_struct accept,
                   long iter)

{
	int i,j;
  int missing_data;
	double mean;
  double lik = 0.0;
  int **n;
  int ***x;
	
  n = data.total_nbr_counts;
	fprintf(outputs_mcmc,"%10ld",iter);
  if (pooled_data) {
    x = current.value.counts;
    for (i = 0; i < data.nbr_demes; ++i) {
      lik += current.prior.M[i];
      for (j = 0; j < data.nbr_loci; ++j) {
        lik += binom_coeff_reads[i][j];
        if (x[i][j][0] > 0) {
          lik += data.reads[i][j][0] * log((double) x[i][j][0] / n[i][j]);
        }
        if (x[i][j][1] > 0) {
          lik += data.reads[i][j][1] * log((double) x[i][j][1] / n[i][j]);
        }
        lik += lgammafn(n[i][j] + 1) - lgammafn(x[i][j][0] + 1) - lgammafn(x[i][j][1] + 1);
        lik += x[i][j][0] * log(current.value.p[i][j][0]) + x[i][j][1] * log(current.value.p[i][j][1]);
        lik += current.psi.Gamma_M[i];
        lik -= current.psi.log_1F1[i][j];
        lik -= current.psi.Gamma_M_pi[i][j];
        lik += current.value.sigma[i][j] * current.value.p[i][j][current.value.kappa[i][j]];
        lik += (current.value.M[i] * current.value.pi[j][0] - 1) * log(current.value.p[i][j][0]) + (current.value.M[i] * current.value.pi[j][1] - 1) * log(current.value.p[i][j][1]);
        lik += current.prior.sigma[i][j];
      }
    }
    for (j = 0; j < data.nbr_loci; ++j) {
      lik += current.prior.delta[j];
      lik += current.prior.pi[j][0];
    }
    lik += current.prior.lambda;
  }
  else {
    x = data.counts;
    for (i = 0; i < data.nbr_demes; ++i) {
      lik += current.prior.M[i];
      for (j = 0; j < data.nbr_loci; ++j) {
        lik += binom_coeff_counts[i][j];
        lik += x[i][j][0] * log(current.value.p[i][j][0]) + x[i][j][1] * log(current.value.p[i][j][1]);
        lik += current.psi.Gamma_M[i];
        lik -= current.psi.log_1F1[i][j];
        lik -= current.psi.Gamma_M_pi[i][j];
        lik += current.value.sigma[i][j] * current.value.p[i][j][current.value.kappa[i][j]];
        lik += (current.value.M[i] * current.value.pi[j][0] - 1) * log(current.value.p[i][j][0]) + (current.value.M[i] * current.value.pi[j][1] - 1) * log(current.value.p[i][j][1]);
        lik += current.prior.sigma[i][j];
      }
    }
    for (j = 0; j < data.nbr_loci; ++j) {
      lik += current.prior.delta[j];
      lik += current.prior.pi[j][0];
    }
    lik += current.prior.lambda;
  }
  fprintf(outputs_mcmc,"%15.5e",lik);
  if (pooled_data) {
    mean = 0.0;
    missing_data = 0;
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        if (update.counts[i][j] > 0.0) {
          mean += accept.counts[i][j] / update.counts[i][j];
        }
        else {
          missing_data += 1;
        }
      }
    }
    mean /= (data.nbr_loci * data.nbr_demes - missing_data);
    fprintf(outputs_mcmc,"%12.6f",mean);
  }
  mean = 0.0;
  for (i = 0; i < data.nbr_demes; ++i) {
    for (j = 0; j < data.nbr_loci; ++j) {
      if (update.p[i][j]) {
        mean += accept.p[i][j] / update.p[i][j];
      }
    }
  }
  mean /= (data.nbr_loci * data.nbr_demes);
  fprintf(outputs_mcmc,"%12.6f",mean);
  mean = 0.0;
	for (i = 0; i < data.nbr_demes; ++i) {
		if (update.M[i]) {
			mean += accept.M[i] / update.M[i];
		}
	}
	mean /= data.nbr_demes;
	fprintf(outputs_mcmc,"%12.6f",mean);
	mean = 0.0;
	for (j = 0; j < data.nbr_loci; ++j) {
		if (update.pi[j]) {
			mean += accept.pi[j] / update.pi[j];
		}
	}
	mean /= data.nbr_loci;
	fprintf(outputs_mcmc,"%12.6f",mean);
  if (!fixed_beta) {
    if (update.beta_mu) {
      fprintf(outputs_mcmc,"%12.6f",accept.beta_mu / update.beta_mu);
    }
    else {
      fprintf(outputs_mcmc,"%12.6f",accept.beta_mu);
    }
    if (update.beta_nu) {
      fprintf(outputs_mcmc,"%12.6f",accept.beta_nu / update.beta_nu);
    }
    else {
      fprintf(outputs_mcmc,"%12.6f",accept.beta_nu);
    }
  }
  if (!fixed_lambda) {
    if (exponential_prior) {
      if (update.lambda) {
        fprintf(outputs_mcmc,"%12.6f",accept.lambda / update.lambda);
      }
      else {
        fprintf(outputs_mcmc,"%12.6f",accept.lambda);
      }
    }
  }
	mean = 0.0;
	for (j = 0; j < data.nbr_loci; ++j) {
		if (update.delta[j]) {
			mean += accept.delta[j] / update.delta[j];
		}
	}
	mean /= data.nbr_loci;
	fprintf(outputs_mcmc,"%12.6f",mean);
	mean = 0.0;
	for (i = 0; i < data.nbr_demes; ++i) {
		for (j = 0; j < data.nbr_loci; ++j) {
			if (update.sigma[i][j]) {
				mean += accept.sigma[i][j] / update.sigma[i][j];
			}
		}
	}
	mean /= (data.nbr_loci * data.nbr_demes);
	fprintf(outputs_mcmc,"%12.6f\n",mean);
  fprintf(trace_m,"%10ld",iter);
	for (i = 0; i < data.nbr_demes; ++i) {
		fprintf(trace_m,"%15.6f",current.value.M[i]);
	}
	fprintf(trace_m,"\n");
  if (!fixed_lambda) {
    fprintf(trace_lambda,"%10ld%12.6f\n",iter,current.value.lambda);
  }
  if (!fixed_beta) {
    fprintf(trace_beta,"%10ld%12.6f%12.6f\n",iter,current.value.beta_a,current.value.beta_b);
  }
  if (verbose) {
    fprintf(trace_freq,"%10ld",iter);
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        fprintf(trace_freq,"%14.6f",current.value.p[i][j][0]);
      }
    }
    fprintf(trace_freq,"\n");
    if (pooled_data) {
      fprintf(trace_counts,"%10ld",iter);
      for (i = 0; i < data.nbr_demes; ++i) {
        for (j = 0; j < data.nbr_loci; ++j) {
          fprintf(trace_counts,"%14d",current.value.counts[i][j][0]);
        }
      }
      fprintf(trace_counts,"\n");
    }
    fprintf(trace_pi,"%10ld",iter);
    for (j = 0; j < data.nbr_loci; ++j) {
      fprintf(trace_pi,"%14.6f",current.value.pi[j][0]);
    }
    fprintf(trace_pi,"\n");
    fprintf(trace_kappa,"%10ld",iter);
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        fprintf(trace_kappa,"%14d",current.value.kappa[i][j]);
      }
    }
    fprintf(trace_kappa,"\n");
    fprintf(trace_sigma,"%10ld",iter);
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < data.nbr_loci; ++j) {
        fprintf(trace_sigma,"%14.6f",current.value.sigma[i][j]);
      }
    }
    fprintf(trace_sigma,"\n");
    fprintf(trace_delta,"%10ld",iter);
    for (j = 0; j < data.nbr_loci; ++j) {
      fprintf(trace_delta,"%14.6f",current.value.delta[j]);
    }
    fprintf(trace_delta,"\n");
  }
  fflush(outputs_mcmc);
	fflush(trace_m);
  if (!fixed_lambda) {
    fflush(trace_lambda);
  }
  if (!fixed_beta) {
    fflush(trace_beta);
  }
  if (verbose) {
    fflush(trace_freq);
    if (pooled_data) {
      fflush(trace_counts);
    }
    fflush(trace_pi);
    fflush(trace_kappa);
    fflush(trace_sigma);
    fflush(trace_delta);
  }
}

void write_summaries(data_struct data,
                     moments_struct postmean,
                     moments_struct postmsqr,
                     char *path,
                     unsigned int pod_analysis)

{
  char filename[256] = "";
  int i,j,k;
  int nbr_spaces;
  int n;
  double gamma,m,p,Q;
  double k0,v,std,t,theta0;
  double *kld;
  double quantile_list[] = {0.5,0.9,0.95,0.98,0.99,0.995,0.999,0.9995,0.9999,UNDEFINED};
  
  kld = (double *) malloc(data.nbr_loci * sizeof(double));
  t = 1 / postmean.lambda;
  for (j = 0; j < data.nbr_loci; ++j) {
    m = postmean.delta[j];
    v = (postmsqr.delta[j] - pow(postmean.delta[j],2)) * data.nbr_loci / (data.nbr_loci - 1);
    k0 = pow(m,2) / v;
    theta0 = m / v;
    kld[j] = log(pow(theta0,k0) / (gammafn(k0) * t)) + (k0 - 1) * (digamma(k0) - log(theta0)) + k0 * (t - theta0) / theta0;
  }
  strcpy(filename,path);
  strcat(filename,SUMRY_M);
  if ((sumry_m = fopen(filename,"w")) == NULL) {
    printf("Can't open %s for output!\n",filename);
    exit(EXIT_FAILURE);
  }
  fprintf(sumry_m,"    deme           mean            std\n");
  for (i = 0; i < data.nbr_demes; ++i) {
    std = sqrt(postmsqr.M[i] - pow(postmean.M[i],2));
    fprintf(sumry_m,"%8d%15.6f%15.6f\n",(i + 1),postmean.M[i],std);
  }
  fclose(sumry_m);
  strcpy(filename,path);
  strcat(filename,SUMRY_PI);
  if ((sumry_pi = fopen(filename,"w")) == NULL) {
    printf("Can't open %s for output!\n",filename);
    exit(EXIT_FAILURE);
  }
  fprintf(sumry_pi,"   locus        mean         std\n");
  for (j = 0; j < data.nbr_loci; ++j) {
    std = sqrt(postmsqr.pi[j] - pow(postmean.pi[j],2));
    fprintf(sumry_pi,"%8d%12.6f%12.6f\n",(j + 1),postmean.pi[j],std);
  }
  fclose(sumry_pi);
  strcpy(filename,path);
  strcat(filename,SUMRY_LAMBDA);
  if ((sumry_lambda = fopen(filename,"w")) == NULL) {
    printf("Can't open %s for output!\n",filename);
    exit(EXIT_FAILURE);
  }
  fprintf(sumry_lambda,"        mean         std\n");
  std = sqrt(postmsqr.lambda - pow(postmean.lambda,2));
  fprintf(sumry_lambda,"%12.6f%12.6f\n",postmean.lambda,std);
  fclose(sumry_lambda);
  strcpy(filename,path);
  strcat(filename,SUMRY_FREQ);
  if ((sumry_freq = fopen(filename,"w")) == NULL) {
    printf("Can't open %s for output!\n",filename);
    exit(EXIT_FAILURE);
  }
  fprintf(sumry_freq,"   locus    deme        mean         std\n");
  for (j = 0; j < data.nbr_loci; ++j) {
    for (i = 0; i < data.nbr_demes; ++i) {
      fprintf(sumry_freq,"%8d%8d",(j + 1),(i + 1));
      std = sqrt(postmsqr.p[i][j] - pow(postmean.p[i][j],2));
      fprintf(sumry_freq,"%12.6f%12.6f\n",postmean.p[i][j],std);
    }
  }
  fclose(sumry_freq);
  if (pooled_data) {
    strcpy(filename,path);
    strcat(filename,SUMRY_COUNTS);
    if ((sumry_counts = fopen(filename,"w")) == NULL) {
      printf("Can't open %s for output!\n",filename);
      exit(EXIT_FAILURE);
    }
    fprintf(sumry_counts,"   locus    deme        mean         std\n");
    for (j = 0; j < data.nbr_loci; ++j) {
      for (i = 0; i < data.nbr_demes; ++i) {
        fprintf(sumry_counts,"%8d%8d",(j + 1),(i + 1));
        std = sqrt(postmsqr.counts[i][j] - pow(postmean.counts[i][j],2));
        fprintf(sumry_counts,"%12.6f%12.6f\n",postmean.counts[i][j],std);
      }
    }
    fclose(sumry_counts);
  }
  if (!fixed_beta) {
    strcpy(filename,path);
    strcat(filename,SUMRY_BETA);
    if ((sumry_beta = fopen(filename,"w")) == NULL) {
      printf("Can't open %s for output!\n",filename);
      exit(EXIT_FAILURE);
    }
    fprintf(sumry_beta,"        mean         std\n");
    std = sqrt(postmsqr.beta_a - pow(postmean.beta_a,2));
    fprintf(sumry_beta,"%12.6f%12.6f\n",postmean.beta_a,std);
    std = sqrt(postmsqr.beta_b - pow(postmean.beta_b,2));
    fprintf(sumry_beta,"%12.6f%12.6f\n",postmean.beta_b,std);
    fclose(sumry_beta);
  }
  strcpy(filename,path);
  strcat(filename,SUMRY_KAPPA);
  if ((sumry_kappa = fopen(filename,"w")) == NULL) {
    printf("Can't open %s for output!\n",filename);
    exit(EXIT_FAILURE);
  }
  fprintf(sumry_kappa,"   locus");
  for (i = 0; i < data.nbr_demes; ++i) {
    nbr_spaces = 12 - 5 - (floor(log10(i + 1)) + 1);
    for (k = 0; k < nbr_spaces; ++k) {
      fprintf(sumry_kappa," ");
    }
    fprintf(sumry_kappa,"deme_%d",(i + 1));
  }
  fprintf(sumry_kappa,"\n");
  for (j = 0; j < data.nbr_loci; ++j) {
    fprintf(sumry_kappa,"%8d",(j + 1));
    for (i = 0; i < data.nbr_demes; ++i) {
      fprintf(sumry_kappa,"%12.6f",postmean.kappa[i][j]);
    }
    fprintf(sumry_kappa,"\n");
  }
  fclose(sumry_kappa);
  strcpy(filename,path);
  strcat(filename,SUMRY_SIGMA);
  if ((sumry_sigma = fopen(filename,"w")) == NULL) {
    printf("Can't open %s for output!\n",filename);
    exit(EXIT_FAILURE);
  }
  fprintf(sumry_sigma,"   locus    deme        mean         std  mean_all_1   std_all_1  mean_all_2   std_all_2\n");
  for (j = 0; j < data.nbr_loci; ++j) {
    for (i = 0; i < data.nbr_demes; ++i) {
      fprintf(sumry_sigma,"%8d%8d",(j + 1),(i + 1));
      std = sqrt(postmsqr.sigma[i][j][2] - pow(postmean.sigma[i][j][2],2));
      fprintf(sumry_sigma,"%12.6f%12.6f",postmean.sigma[i][j][2],std);
      std = sqrt(postmsqr.sigma[i][j][0] - pow(postmean.sigma[i][j][0],2));
      fprintf(sumry_sigma,"%12.6f%12.6f",postmean.sigma[i][j][0],std);
      std = sqrt(postmsqr.sigma[i][j][1] - pow(postmean.sigma[i][j][1],2));
      fprintf(sumry_sigma,"%12.6f%12.6f\n",postmean.sigma[i][j][1],std);
    }
  }
  fclose(sumry_sigma);
  strcpy(filename,path);
  strcat(filename,SUMRY_DELTA);
  if ((sumry_delta = fopen(filename,"w")) == NULL) {
    printf("Can't open %s for output!\n",filename);
    exit(EXIT_FAILURE);
  }
  fprintf(sumry_delta,"   locus        mean         std         KLD\n");
  for (j = 0; j < data.nbr_loci; ++j) {
    std = sqrt(postmsqr.delta[j] - pow(postmean.delta[j],2));
    fprintf(sumry_delta,"%8d%12.6f%12.6f%12.6f\n",(j + 1),postmean.delta[j],std,kld[j]);
  }
  fclose(sumry_delta);
  if (pod_analysis) {
    strcpy(filename,path);
    strcat(filename,KLD_QUTILE);
    if ((kld_qutile = fopen(filename,"w")) == NULL) {
      printf("Can't open %s for output!\n",filename);
      exit(EXIT_FAILURE);
    }
    qsort(kld,data.nbr_loci,sizeof(double),compare);
    fprintf(kld_qutile," quantile         KLD\n");
    for (i = 0; quantile_list[i] != UNDEFINED; i++) {
      p = quantile_list[i];
      m = 1.0 - p;
      n = data.nbr_loci;
      j = floor(n * p + m);
      gamma = n * p + m - j;
      Q = (1 - gamma) * kld[j - 1] + gamma * kld[j];
      fprintf(kld_qutile,"%8.2f%%%12.6f\n",(p * 100),Q);
    }
    fclose(kld_qutile);
  }
  free(kld);
}

int compare(const void *a, const void *b)

{
  int reslt;
  
  if ( *(double*)a <  *(double*)b ) reslt = -1;
  if ( *(double*)a == *(double*)b ) reslt = 0;
  if ( *(double*)a >  *(double*)b ) reslt = 1;
  return(reslt);
}

double draw_from_psi(double sigma,
                     double M,
                     double pi)

{
  int MAX_ITER = 1e4;
  int accept;
  double hypergeometric_1F1_1,hypergeometric_1F1_2,hypergeometric_1F1_3;
  double mean,mean_sqr;
  double alpha,beta;
  double c,Y,U;
  double p;
  
  c = exp(sigma) / hypergeometric_1_F_1(M * pi,M,sigma);
  if (c < MAX_ITER) {
    accept = 0;
    do {
      Y = rbeta(M * pi,M * (1 - pi),mtss[omp_get_thread_num()]);
      U = genrand_real1(mtss[omp_get_thread_num()]);
      accept = U < exp(sigma * (Y - 1));
    } while (!accept);
    p = Y;
  }
  else {
    hypergeometric_1F1_1 = hypergeometric_1_F_1(pi * M,M,sigma);
    hypergeometric_1F1_2 = hypergeometric_1_F_1(pi * M + 1,M + 1,sigma);
    hypergeometric_1F1_3 = hypergeometric_1_F_1(pi * M + 2,M + 2,sigma);
    mean = pi * hypergeometric_1F1_2 / hypergeometric_1F1_1;
    mean_sqr = pi * (pi * M + 1) / (M + 1) * hypergeometric_1F1_3 / hypergeometric_1F1_1;
    alpha = mean * (mean_sqr - mean) / (pow(mean,2) - mean_sqr);
    beta = alpha * (1 / mean - 1);
    p = rbeta(alpha,beta,mtss[omp_get_thread_num()]);
  }
  return(p);
}

void write_pseudo_observed_data(data_struct data,
                                moments_struct postmean,
                                char *pod_path,
                                char *extended_filename,
                                char *pod_file_name)

{
  int i,j,k;
  int idx_dta_locus,idx_pod_locus,nbr_loci;
  int cnt,tick10,tick100,trick;
  int sum_counts[2];
  int sum_reads[2];
  unsigned int random_locus;
  double delta,kappa,lambda,sigma;
  double p_tilde,pi,pi_tilde;
  double rnd;
  double sum_p[2];
  double **cum_p,*p;
  int ***counts,***reads;
  char filename[256] = "";
  char tmp[256] = "";
  char *ptr;
  FILE *podfile;
  
  p = (double *) malloc(data.nbr_demes * sizeof(double));
  cum_p = (double **) malloc(data.nbr_demes * sizeof(double *));
  for (i = 0; i < data.nbr_demes; ++i) {
    cum_p[i] = (double *) malloc(2 * sizeof(double));
  }
#ifdef __linux__
  char path_separator[2] = "/";
#elif __APPLE__
  char path_separator[2] = "/";
#else
  char path_separator[2] = "\\";
#endif
  strcpy(tmp,extended_filename);
  if ((ptr = strtok(tmp,path_separator)) != NULL) {
    do {
      strcpy(filename,ptr);
    }
    while ((ptr = strtok(NULL,path_separator)) != NULL);                       // This grabs the next token, by passing NULL to strok(), and continues if it is !NULL
  }
  strcpy(pod_file_name,pod_path);
  strcat(pod_file_name,PODFILE);
  strcat(pod_file_name,filename);
  printf("Generating file `%s'...\n",pod_file_name);
  fprintf(logfile,"Generating file `%s'... ",pod_file_name);
  if (pod_nbr_loci > 0) {
    nbr_loci = pod_nbr_loci;
    random_locus = TRUE;
  }
  else {
    nbr_loci = data.nbr_loci;
    random_locus = FALSE;
  }
  trick = (int) (nbr_loci / 100) * 100;
  tick10 = (int) (trick / 10);
  tick100 = (int) (trick / 100);
  cnt = 0;
  counts = (int ***) malloc(data.nbr_demes * sizeof(int **));
  for (i = 0; i < data.nbr_demes; ++i) {
    counts[i] = (int **) malloc(nbr_loci * sizeof(int *));
    for (j = 0; j < nbr_loci; ++j) {
      counts[i][j] = (int *) malloc(2 * sizeof(int));
      counts[i][j][0] = counts[i][j][1] = 0;
    }
  }
  if (pooled_data) {
    reads = (int ***) malloc(data.nbr_demes * sizeof(int **));
    for (i = 0; i < data.nbr_demes; ++i) {
      reads[i] = (int **) malloc(nbr_loci * sizeof(int *));
      for (j = 0; j < nbr_loci; ++j) {
        reads[i][j] = (int *) malloc(2 * sizeof(int));
        reads[i][j][0] = reads[i][j][1] = 0;
      }
    }
    lambda = postmean.lambda;
    idx_pod_locus = 0;
    do {
      if (random_locus) {
        idx_dta_locus = (int) (genrand_real2(mtss[omp_get_thread_num()]) * data.nbr_loci);
      }
      else {
        idx_dta_locus = idx_pod_locus;
      }
      delta = rexp(lambda,mtss[omp_get_thread_num()]);                   // (!) rexp() takes 'scale' as an argument, not 'rate'
      if (!fixed_beta) {
        pi = rbeta(postmean.beta_a,postmean.beta_b,mtss[omp_get_thread_num()]);
        if (pi < MIN_PI) pi = MIN_PI;
        if (pi > (1.0 - MIN_PI)) pi = 1.0 - MIN_PI;
      }
      else {
        pi = postmean.pi[idx_dta_locus];
      }
      for (i = 0; i < data.nbr_demes; ++i) {
        sigma = rexp(delta,mtss[omp_get_thread_num()]);
        if (sigma > max_sigma) {
          sigma = max_sigma;
        }
        kappa = rbinom(1,postmean.kappa[i][idx_dta_locus],mtss[omp_get_thread_num()]);
        pi_tilde = pi * (1 - kappa) + (1 - pi) * kappa;
        p_tilde = draw_from_psi(sigma,postmean.M[i],pi_tilde);
        p[i] = p_tilde * (1 - kappa) + (1 - p_tilde) * kappa;
        if (i == 0) {
          cum_p[i][0] = p[i];
          cum_p[i][1] = 1.0 - p[i];
        }
        else {
          cum_p[i][0] = cum_p[i - 1][0] + p[i];
          cum_p[i][1] = cum_p[i - 1][1] + (1.0 - p[i]);
        }
        counts[i][idx_pod_locus][0] = (int) rbinom(data.total_nbr_counts[i][idx_dta_locus],p[i],mtss[omp_get_thread_num()]);
        counts[i][idx_pod_locus][1] = data.total_nbr_counts[i][idx_dta_locus] - counts[i][idx_pod_locus][0];
        reads[i][idx_pod_locus][0] = (int) rbinom(data.total_nbr_reads[i][idx_dta_locus],(double) counts[i][idx_pod_locus][0] / data.total_nbr_counts[i][idx_dta_locus],mtss[omp_get_thread_num()]);
        reads[i][idx_pod_locus][1] = data.total_nbr_reads[i][idx_dta_locus] - reads[i][idx_pod_locus][0];
      }
      sum_p[0] = cum_p[data.nbr_demes - 1][0];
      sum_p[1] = cum_p[data.nbr_demes - 1][1];
      sum_reads[0] = sum_reads[1] = 0;
      for (i = 0; i < data.nbr_demes; ++i) {
        cum_p[i][0] /= sum_p[0];
        cum_p[i][1] /= sum_p[1];
        sum_reads[0] += reads[i][idx_pod_locus][0];
        sum_reads[1] += reads[i][idx_pod_locus][1];
      }
      if (sum_reads[0] < min_nbr_reads) {
        for (i = 0; i < data.nbr_demes; ++i) {
          reads[i][idx_pod_locus][0] = 0;
        }
        for (j = 0; j < min_nbr_reads; ++j) {
          k = 0;
          rnd = genrand_real2(mtss[omp_get_thread_num()]);
          while (rnd > cum_p[k][0]) k += 1;
          reads[k][idx_pod_locus][0] += 1;
        }
        for (i = 0; i < data.nbr_demes; ++i) {
          reads[i][idx_pod_locus][1] = data.total_nbr_reads[i][idx_dta_locus] - reads[i][idx_pod_locus][0];
        }
      }
      if (sum_reads[1] < min_nbr_reads) {
        for (i = 0; i < data.nbr_demes; ++i) {
          reads[i][idx_pod_locus][1] = 0;
        }
        for (j = 0; j < min_nbr_reads; ++j) {
          k = 0;
          rnd = genrand_real2(mtss[omp_get_thread_num()]);
          while (rnd > cum_p[k][1]) k += 1;
          reads[k][idx_pod_locus][1] += 1;
        }
        for (i = 0; i < data.nbr_demes; ++i) {
          reads[i][idx_pod_locus][0] = data.total_nbr_reads[i][idx_dta_locus] - reads[i][idx_pod_locus][1];
        }
      }
      if (idx_pod_locus < trick) {
        if ((idx_pod_locus / tick10) * tick10 == idx_pod_locus) {
          if (cnt == 0) {
            printf("\n starting [");
          }
          else {
            printf("]\n%3d%% done [",cnt);
          }
          cnt += 10;
        }
        if ((idx_pod_locus / tick100) * tick100 == idx_pod_locus) {
          printf(".");
          fflush(stdout);
        }
      }
      ++idx_pod_locus;
    } while (idx_pod_locus < nbr_loci);
    if ((podfile = fopen(pod_file_name,"w")) == NULL) {
      printf("Can't open %s for output!\n",pod_file_name);
      exit(EXIT_FAILURE);
    }
    fprintf(podfile,"%d\n",data.nbr_demes);
    fprintf(podfile,"%d\n",nbr_loci);
    if (pooled_data) {
      for (i = 0; i < data.nbr_demes; ++i) {
        if (i < (data.nbr_demes - 1)) {
          fprintf(podfile,"%d\t",data.pool_size[i]);
        }
        else {
          fprintf(podfile,"%d",data.pool_size[i]);
        }
      }
      fprintf(podfile,"\n");
    }
    for (j = 0; j < nbr_loci; ++j) {
      for (i = 0; i < data.nbr_demes; ++i) {
        if (pooled_data) {
          if (i < (data.nbr_demes - 1)) {
            fprintf(podfile,"%d\t%d\t",reads[i][j][0],reads[i][j][1]);
          }
          else {
            fprintf(podfile,"%d\t%d",reads[i][j][0],reads[i][j][1]);
          }
        }
        else {
          if (i < (data.nbr_demes - 1)) {
            fprintf(podfile,"%d\t%d\t",counts[i][j][0],counts[i][j][1]);
          }
          else {
            fprintf(podfile,"%d\t%d",counts[i][j][0],counts[i][j][1]);
          }
        }
      }
      fprintf(podfile,"\n");
    }
    for (i = 0; i < data.nbr_demes; ++i) {
      for (j = 0; j < nbr_loci; ++j) {
        free(reads[i][j]);
      }
      free(reads[i]);
    }
    free(reads);
  }
  else {
    lambda = postmean.lambda;
    idx_pod_locus = 0;
    do {
      if (random_locus) {
        idx_dta_locus = (int) (genrand_real2(mtss[omp_get_thread_num()]) * data.nbr_loci);
      }
      else {
        idx_dta_locus = idx_pod_locus;
      }
      delta = rexp(lambda,mtss[omp_get_thread_num()]);                   // (!) rexp() takes 'scale' as an argument, not 'rate'
      if (!fixed_beta) {
        pi = rbeta(postmean.beta_a,postmean.beta_b,mtss[omp_get_thread_num()]);
        if (pi < MIN_PI) pi = MIN_PI;
        if (pi > (1.0 - MIN_PI)) pi = 1.0 - MIN_PI;
      }
      else {
        pi = postmean.pi[idx_dta_locus];
      }
      for (i = 0; i < data.nbr_demes; ++i) {
        sigma = rexp(delta,mtss[omp_get_thread_num()]);
        if (sigma > max_sigma) {
          sigma = max_sigma;
        }
        kappa = rbinom(1,postmean.kappa[i][idx_dta_locus],mtss[omp_get_thread_num()]);
        pi_tilde = pi * (1 - kappa) + (1 - pi) * kappa;
        p_tilde = draw_from_psi(sigma,postmean.M[i],pi_tilde);
        p[i] = p_tilde * (1 - kappa) + (1 - p_tilde) * kappa;
        if (i == 0) {
          cum_p[i][0] = p[i];
          cum_p[i][1] = 1.0 - p[i];
        }
        else {
          cum_p[i][0] = cum_p[i - 1][0] + p[i];
          cum_p[i][1] = cum_p[i - 1][1] + (1.0 - p[i]);
        }
        counts[i][idx_pod_locus][0] = (int) rbinom(data.total_nbr_counts[i][idx_dta_locus],p[i],mtss[omp_get_thread_num()]);
        counts[i][idx_pod_locus][1] = data.total_nbr_counts[i][idx_dta_locus] - counts[i][idx_pod_locus][0];
      }
      sum_p[0] = cum_p[data.nbr_demes - 1][0];
      sum_p[1] = cum_p[data.nbr_demes - 1][1];
      sum_counts[0] = sum_counts[1] = 0;
      for (i = 0; i < data.nbr_demes; ++i) {
        cum_p[i][0] /= sum_p[0];
        cum_p[i][1] /= sum_p[1];
        sum_counts[0] += counts[i][idx_pod_locus][0];
        sum_counts[1] += counts[i][idx_pod_locus][1];
      }
      if (sum_counts[0] < min_nbr_counts) {
        for (i = 0; i < data.nbr_demes; ++i) {
          counts[i][idx_pod_locus][0] = 0;
        }
        for (j = 0; j < min_nbr_counts; ++j) {
          k = 0;
          rnd = genrand_real2(mtss[omp_get_thread_num()]);
          while (rnd > cum_p[k][0]) k += 1;
          counts[k][idx_pod_locus][0] += 1;
        }
        for (i = 0; i < data.nbr_demes; ++i) {
          counts[i][idx_pod_locus][1] = data.total_nbr_counts[i][idx_dta_locus] - counts[i][idx_pod_locus][0];
        }
      }
      if (sum_counts[1] < min_nbr_counts) {
        for (i = 0; i < data.nbr_demes; ++i) {
          counts[i][idx_pod_locus][1] = 0;
        }
        for (j = 0; j < min_nbr_counts; ++j) {
          k = 0;
          rnd = genrand_real2(mtss[omp_get_thread_num()]);
          while (rnd > cum_p[k][1]) k += 1;
          counts[k][idx_pod_locus][1] += 1;
        }
        for (i = 0; i < data.nbr_demes; ++i) {
          counts[i][idx_pod_locus][0] = data.total_nbr_counts[i][idx_dta_locus] - counts[i][idx_pod_locus][1];
        }
      }
      if (idx_pod_locus < trick) {
        if ((idx_pod_locus / tick10) * tick10 == idx_pod_locus) {
          if (cnt == 0) {
            printf("\n starting [");
          }
          else {
            printf("]\n%3d%% done [",cnt);
          }
          cnt += 10;
        }
        if ((idx_pod_locus / tick100) * tick100 == idx_pod_locus) {
          printf(".");
          fflush(stdout);
        }
      }
      ++idx_pod_locus;
    } while (idx_pod_locus < nbr_loci);
    if ((podfile = fopen(pod_file_name,"w")) == NULL) {
      printf("Can't open %s for output!\n",pod_file_name);
      exit(EXIT_FAILURE);
    }
    fprintf(podfile,"%d\n",data.nbr_demes);
    fprintf(podfile,"%d\n",nbr_loci);
    for (j = 0; j < nbr_loci; ++j) {
      for (i = 0; i < data.nbr_demes; ++i) {
        if (i < (data.nbr_demes - 1)) {
          fprintf(podfile,"%d\t%d\t",counts[i][j][0],counts[i][j][1]);
        }
        else {
          fprintf(podfile,"%d\t%d",counts[i][j][0],counts[i][j][1]);
        }
      }
      fprintf(podfile,"\n");
    }
  }
  fclose(podfile);
  for (i = 0; i < data.nbr_demes; ++i) {
    for (j = 0; j < nbr_loci; ++j) {
      free(counts[i][j]);
    }
    free(counts[i]);
  }
  free(counts);
  for (i = 0; i < data.nbr_demes; ++i) {
    free(cum_p[i]);
  }
  free(cum_p);
  free(p);
  printf("]\n%3d%% done !\n\n",cnt);
  fprintf(logfile,"OK\n\n");
  printf("The pseudo-observed data consist in %d SNPs and %d sampled populations\n\n",nbr_loci,data.nbr_demes);
  fprintf(logfile,"The pseudo-observed data consist in %d SNPs and %d sampled populations\n\n",nbr_loci,data.nbr_demes);
}

void copy_data(data_struct data,
               char *extended_filename,
               char *path)

{
  int i,j;
  int file_exists;
  char fullpath[256] = "";
  char filename[256] = "";
  char tmp[256] = "";
  FILE *data_copy = NULL;
  char *ptr;
  
#ifdef __linux__
  char path_separator[2] = "/";
#elif __APPLE__
  char path_separator[2] = "/";
#else
  char path_separator[2] = "\\";
#endif
  strcpy(tmp,extended_filename);
  if ((ptr = strtok(tmp,path_separator)) != NULL) {
    do {
      strcpy(filename,ptr);
    }
    while ((ptr = strtok(NULL,path_separator)) != NULL);                       // This grabs the next token, by passing NULL to strok(), and continues if it is !NULL
  }
  strcpy(fullpath,path);
  strcat(fullpath,filename);
  if ((data_copy = fopen(fullpath,"r")) == NULL) {
    file_exists = FALSE;
  }
  else {
    file_exists = TRUE;
    fclose(data_copy);
  }
  if (!file_exists) {
    if ((data_copy = fopen(fullpath,"w")) == NULL) {
      printf("Cannot copy the data in %s directory!\n",path);
      exit(EXIT_FAILURE);
    }
    if (pooled_data) {
      fprintf(data_copy,"%d\n",data.nbr_demes);
      fprintf(data_copy,"%d\n",data.nbr_loci);
      for (i = 0; i < (data.nbr_demes - 1); ++i) {
        fprintf(data_copy,"%d\t",data.pool_size[i]);
      }
      fprintf(data_copy,"%d\n",data.pool_size[i]);
      for (j = 0; j < data.nbr_loci; ++j) {
        for (i = 0; i < data.nbr_demes; ++i) {
          if (i < (data.nbr_demes - 1)) {
            fprintf(data_copy,"%d\t%d\t",data.reads[i][j][0],data.reads[i][j][1]);
          }
          else {
            fprintf(data_copy,"%d\t%d",data.reads[i][j][0],data.reads[i][j][1]);
          }
        }
        fprintf(data_copy,"\n");
      }
    }
    else {
      fprintf(data_copy,"%d\n",data.nbr_demes);
      fprintf(data_copy,"%d\n",data.nbr_loci);
      for (j = 0; j < data.nbr_loci; ++j) {
        for (i = 0; i < data.nbr_demes; ++i) {
          if (i < (data.nbr_demes - 1)) {
            fprintf(data_copy,"%d\t%d\t",data.counts[i][j][0],data.counts[i][j][1]);
          }
          else {
            fprintf(data_copy,"%d\t%d",data.counts[i][j][0],data.counts[i][j][1]);
          }
        }
        fprintf(data_copy,"\n");
      }
    }
    fclose(data_copy);
  }
}

void write_readme(data_struct data,
                  char *file_name,
                  char *path,
                  char *pod_file_name)

{
  int i,j;
  int nbr_lines = 2;
  char fullpath[256] = "";
  FILE *readme = NULL;
  
  strcpy(fullpath,path);
  strcat(fullpath,README);
  if ((readme = fopen(fullpath,"w")) == NULL) {
    printf("Cannot copy the data in %s directory!\n",path);
    exit(EXIT_FAILURE);
  }
  fprintf(readme,"The pseudo-observed data in the file '%s' were simulated using the posterior means of the model parameters, from the analysis of the following dataset:\n\n",pod_file_name);
  fprintf(readme,"----- file '%s' begins here -----\n",file_name);
  fprintf(readme,"%d\n",data.nbr_demes);
  fprintf(readme,"%d\n",data.nbr_loci);
  if (pooled_data) {
    for (i = 0; i < data.nbr_demes; ++i) {
      if (i < (data.nbr_demes - 1)) {
        fprintf(readme,"%d\t",data.pool_size[i]);
      }
      else {
        fprintf(readme,"%d",data.pool_size[i]);
      }
    }
    fprintf(readme,"\n");
  }
  for (j = 0; j < nbr_lines; ++j) {
    for (i = 0; i < data.nbr_demes; ++i) {
      if (pooled_data) {
        if (i < (data.nbr_demes - 1)) {
          fprintf(readme,"%d\t%d\t",data.reads[i][j][0],data.reads[i][j][1]);
        }
        else {
          fprintf(readme,"%d\t%d",data.reads[i][j][0],data.reads[i][j][1]);
        }
      }
      else {
        if (i < (data.nbr_demes - 1)) {
          fprintf(readme,"%d\t%d\t",data.counts[i][j][0],data.counts[i][j][1]);
        }
        else {
          fprintf(readme,"%d\t%d",data.counts[i][j][0],data.counts[i][j][1]);
        }
      }
    }
    fprintf(readme,"\n");
  }
  fprintf(readme,"...\n");
  for (j = (data.nbr_loci - nbr_lines); j < data.nbr_loci; ++j) {
    for (i = 0; i < data.nbr_demes; ++i) {
      if (pooled_data) {
        if (i < (data.nbr_demes - 1)) {
          fprintf(readme,"%d\t%d\t",data.reads[i][j][0],data.reads[i][j][1]);
        }
        else {
          fprintf(readme,"%d\t%d",data.reads[i][j][0],data.reads[i][j][1]);
        }
      }
      else {
        if (i < (data.nbr_demes - 1)) {
          fprintf(readme,"%d\t%d\t",data.counts[i][j][0],data.counts[i][j][1]);
        }
        else {
          fprintf(readme,"%d\t%d",data.counts[i][j][0],data.counts[i][j][1]);
        }
      }
    }
    fprintf(readme,"\n");
  }
  fprintf(readme,"-----  file '%s' ends here  -----\n\n",file_name);
  if (pod_nbr_loci > 0) {
    fprintf(readme,"A total of %d loci were simulated. To that end, the posterior means of the pi parameters were sampled with replacement from the %d loci of the original data.\n",pod_nbr_loci,data.nbr_loci);
  }
  else {
    fprintf(readme,"The posterior means of the pi parameters were taken from the %d loci of the original data.\n",data.nbr_loci);
  }
  fclose(readme);
}

void print_model_parameters(unsigned long seed,
                            int n_threads)

{
  printf("------------------------------------------------------------------------------------\n");
  fprintf(logfile,"------------------------------------------------------------------------------------\n");
  if (fixed_beta) {
    printf("Beta prior distribution of pi is fixed (fixed_beta)\n");
    fprintf(logfile,"Beta prior distribution of pi is fixed (fixed_beta)\n");
    printf("    with shape parameters (beta_a,beta_b)                             = (%f,%f)\n\n",beta_a,beta_b);
    fprintf(logfile,"    with shape parameters (beta_a,beta_b)                             = (%f,%f)\n\n",beta_a,beta_b);
  }
  if (fixed_lambda) {
    printf("Parameter lambda is fixed at %.2f (fixed_lambda)\n",fixed_lambda_value);
    fprintf(logfile,"Parameter lambda is fixed at %.2f (fixed_lambda)\n",fixed_lambda_value);
  }
  else {
    if (inverse_gamma_prior) {
      printf("Prior distribution of lambda is inverse gamma (lambda_prior)\n");
      fprintf(logfile,"Prior distribution of lambda is inverse gamma (lambda_prior)\n");
      printf("    with shape parameter (invgam_shape)                               = %f\n",inverse_gamma_shape);
      fprintf(logfile,"    with shape parameter (invgam_shape)                               = %f\n",inverse_gamma_shape);
      printf("    and rate parameter (invgam_rate)                                  = %f\n",inverse_gamma_rate);
      fprintf(logfile,"    and rate parameter (invgam_rate)                                  = %f\n",inverse_gamma_rate);
    }
    if (exponential_prior) {
      printf("Prior distribution of lambda is exponential (lambda_prior)\n");
      fprintf(logfile,"Prior distribution of lambda is exponential (lambda_prior)\n");
      printf("    with rate parameter (captl_lambda)                                = %f\n",capital_lambda);
      fprintf(logfile,"    with rate parameter (captl_lambda)                                = %f\n",capital_lambda);
    }
  }
  printf("\n");
  fprintf(logfile,"\n");
  printf("Number of threads used (threads)                                      = %d\n",n_threads);
  fprintf(logfile,"Number of threads used (threads)                                      = %d\n",n_threads);
  printf("Random number generator's seed (seed)                                 = %ld\n\n",seed);
  fprintf(logfile,"Random number generator's seed (seed)                                 = %ld\n\n",seed);
  printf("Length of the burn-in period (burnin)                                 = %ld\n",burn_in);
  fprintf(logfile,"Length of the burn-in period (burnin)                                 = %ld\n",burn_in);
  printf("Run length of the Markov chain (length)                               = %ld\n",chain_length);
  fprintf(logfile,"Run length of the Markov chain (length)                               = %ld\n",chain_length);
  printf("Thinning interval (thin)                                              = %d\n",step);
  fprintf(logfile,"Thinning interval (thin)                                              = %d\n",step);
  printf("Number of MCMC samples (length / thin)                                = %.0f\n",floor(chain_length / step));
  fprintf(logfile,"Number of MCMC samples (length / thin)                                = %.0f\n",floor(chain_length / step));
  printf("Number of pilot studies (npilot)                                      = %d\n",n_pilot);
  fprintf(logfile,"Number of pilot studies (npilot)                                      = %d\n",n_pilot);
  printf("Length of each pilot study (lpilot)                                   = %d\n\n",l_pilot);
  fprintf(logfile,"Length of each pilot study (lpilot)                                   = %d\n\n",l_pilot);
  if (pooled_data) {
    printf("The data consist in read counts from pooled DNA samples (pool)\n\n");
    fprintf(logfile,"The data consist in read counts from pooled DNA samples (pool)\n\n");
  }
  printf("Lower bound of the interval for M (min_M)                             = %f\n",min_M);
	fprintf(logfile,"Lower bound of the interval for M (min_M)                             = %f\n",min_M);
	printf("Upper bound of the interval for M (max_M)                             = %.2f\n",max_M);
	fprintf(logfile,"Upper bound of the interval for M (max_M)                             = %.2f\n",max_M);
	printf("Upper bound of the interval for sigma (max_sig)                       = %.2f\n",max_sigma);
	fprintf(logfile,"Upper bound of the interval for sigma (max_sig)                       = %.2f\n",max_sigma);
  printf("Initial half window width for updates of allele counts (dlt_cnt)      = %d\n",init_delta_counts);
	fprintf(logfile,"Initial half window width for updates of allele counts (dlt_cnt)      = %d\n",init_delta_counts);
  printf("Initial half window width for updates of p (dlt_p)                    = %f\n",init_delta_p);
	fprintf(logfile,"Initial half window width for updates of p (dlt_p)                    = %f\n",init_delta_p);
	printf("Initial SD of the lognormal for updates of M (dlt_M)                  = %f\n",init_delta_M);
	fprintf(logfile,"Initial SD of the lognormal for updates of M (dlt_M)                  = %f\n",init_delta_M);
	printf("Initial half window width for updates of pi (dlt_pi)                  = %f\n",init_delta_pi);
	fprintf(logfile,"Initial half window width for updates of pi (dlt_pi)                  = %f\n",init_delta_pi);
	printf("Initial SD of the lognormal for updates of sigma (dlt_sig)            = %f\n",init_delta_sigma);
	fprintf(logfile,"Initial SD of the lognormal for updates of sigma (dlt_sig)            = %f\n",init_delta_sigma);
	printf("Initial SD of the lognormal for updates of delta (dlt_del)            = %f\n",init_delta_delta);
	fprintf(logfile,"Initial SD of the lognormal for updates of delta (dlt_del)            = %f\n",init_delta_delta);
  if (!fixed_lambda && exponential_prior) {
    printf("Initial SD of the lognormal for updates of lambda (dlt_lam)           = %f\n",init_delta_lambda);
    fprintf(logfile,"Initial SD of the lognormal for updates of lambda (dlt_lam)         = %f\n",init_delta_lambda);
  }
  if (!fixed_beta) {
    printf("Initial half window width for updates of mu (dlt_beta_mu)             = %f\n",init_delta_beta_mu);
    fprintf(logfile,"Initial half window width for updates of mu (dlt_beta_mu)             = %f\n",init_delta_beta_mu);
    printf("Initial SD of the lognormal for updates of nu (dlt_beta_nu)           = %f\n",init_delta_beta_nu);
    fprintf(logfile,"Initial SD of the lognormal for updates of nu (dlt_beta_nu)           = %f\n",init_delta_beta_nu);
  }
  if (calibration) {
    printf("\nCalibration of the Kullback-Leibler divergence (calibration)\n");
    fprintf(logfile,"\nCalibration of the Kullback-Leibler divergence (calibration)\n");
  }
  if (calibration_only) {
    printf("\nCalibration of the Kullback-Leibler divergence (calibration_only)\n");
    fprintf(logfile,"\nCalibration of the Kullback-Leibler divergence (calibration_only)\n");
  }
  if (pod_nbr_loci > 0) {
    printf("Number of loci to be simulated for calibration (pod_nbr_loci)         = %d\n",pod_nbr_loci);
    fprintf(logfile,"Number of loci to be simulated for calibration (pod_nbr_loci)         = %d\n",pod_nbr_loci);
  }
  if (verbose) {
    printf("\nPrint the traces of all parameters (verbose)\n");
    fprintf(logfile,"\nPrint the traces of all parameters (verbose)\n");
  }
  printf("------------------------------------------------------------------------------------\n\n");
  fprintf(logfile,"------------------------------------------------------------------------------------\n\n");
}

void print_F_ST(data_struct data)

{
  double fst;
  
  if (pooled_data) {
    fst = F_ST_poolseq(data.reads,data.total_nbr_reads,data.nbr_loci,data.nbr_demes,data.pool_size);

  }
  else {
    fst = F_ST(data.counts,data.total_nbr_counts,data.nbr_loci,data.nbr_demes);
  }
  printf("------------------------------------------------------------------------------------\n");
  fprintf(logfile,"------------------------------------------------------------------------------------\n");
  printf("Overall genetic differentiation (F_ST)                                = %6.4f\n",fst);
  fprintf(logfile,"Overall genetic differentiation (F_ST)                                = %6.4f\n",fst);
  printf("------------------------------------------------------------------------------------\n\n");
  fprintf(logfile,"------------------------------------------------------------------------------------\n\n");
}

void compute_autocorrelation(moments_struct postmean,
                             moments_struct postmsqr,
                             char *path)

{
  int dummy,i,l,nbr_values;
  double var;
  char X;
  char filename[256] = "";
  int lag[5] = {0,1,5,10,50};
  double autocov[5],autocorr[5];
  double *lambda;
  
  nbr_values = floor(chain_length / step);
  lambda = (double *) malloc(nbr_values * sizeof(double));
  strcpy(filename,path);
  strcat(filename,TRACE_LAMBDA);
  sumry_lambda = fopen(filename,"r");
  while(!((X = getc(sumry_lambda)) == '\n' || X == '\f' || X == '\r'));
  for (i = 0; i < nbr_values; ++i) {
    fscanf(sumry_lambda,"%d",&dummy);
    fscanf(sumry_lambda,"%lf",&lambda[i]);
  }
  fclose(sumry_lambda);
  var = postmsqr.lambda - pow(postmean.lambda,2);
  for (l = 0; l < 5; ++l) {
    autocov[l] = 0.0;
    for (i = 0; i < (nbr_values - lag[l]); ++i) {
      autocov[l] += (lambda[i] - postmean.lambda) * (lambda[i + lag[l]] - postmean.lambda);
    }
    autocov[l] /= (nbr_values - lag[l]);
    autocorr[l] = autocov[l] / (var * nbr_values / (nbr_values - lag[l]));
  }
  printf("Autocorrelation measure for the trace of the (hyper-)parameter lambda:\n");
  fprintf(logfile,"Autocorrelation measure for the trace of the (hyper-)parameter lambda:\n");
  printf("(the lag is expressed as the number of iterations):\n\n");
  fprintf(logfile,"(the lag is expressed as the number of iterations):\n\n");
  for (l = 0; l < 5; ++l) {
    printf("Lag %d\t %f\n",lag[l] * step,autocorr[l]);
    fprintf(logfile,"Lag %d\t %f\n",lag[l] * step,autocorr[l]);
  }
  printf("\n");
  fprintf(logfile,"\n");
  free(lambda);
}

void print_effective_sample_size(data_struct data,
                                 moments_struct postmean,
                                 moments_struct postmsqr,
                                 char *path)

{
  int dummy,i,j,sample_size;
  double effective_sample_size,mean,msqr,var;
  double *beta_a,*beta_b,*lambda,*lik,**M;
  char X;
  char filename[256] = "";
  
  sample_size = floor(chain_length / step);
  printf("------------------------------------------------------------------------------------\n");
  fprintf(logfile,"------------------------------------------------------------------------------------\n");
  printf("Computation of the effective sample size (ESS)\n\n");
  fprintf(logfile,"Computation of the effective sample size (ESS)\n\n");
  
  lik = (double *) malloc(sample_size * sizeof(double));
  mean = msqr = 0.0;
  strcpy(filename,path);
  strcat(filename,OUTPUTS_MCMC);
  outputs_mcmc = fopen(filename,"r");
  while(!((X = getc(outputs_mcmc)) == '\n' || X == '\f' || X == '\r'));
  for (i = 0; i < sample_size; ++i) {
    fscanf(outputs_mcmc,"%d",&dummy);
    fscanf(outputs_mcmc,"%lf",&lik[i]);
    while(!((X = getc(outputs_mcmc)) == '\n' || X == '\f' || X == '\r'));
    mean += lik[i];
    msqr += lik[i] * lik[i];
  }
  fclose(outputs_mcmc);
  mean /= sample_size;
  msqr /= sample_size;
  var = msqr - pow(mean,2);
  effective_sample_size = compute_effective_sample_size(lik,mean,var,sample_size);
  printf("\tlog posterior density                           = %f\n",effective_sample_size);
  fprintf(logfile,"\tlog posterior density                           = %f\n",effective_sample_size);
  free(lik);
  
  M = (double **) malloc(data.nbr_demes * sizeof(double *));
  for (i = 0; i < data.nbr_demes; ++i) {
    M[i] = (double *) malloc(sample_size * sizeof(double));
  }
  strcpy(filename,path);
  strcat(filename,TRACE_M);
  sumry_m = fopen(filename,"r");
  while(!((X = getc(sumry_m)) == '\n' || X == '\f' || X == '\r'));
  for (j = 0; j < sample_size; ++j) {
    fscanf(sumry_m,"%d",&dummy);
    for (i = 0; i < data.nbr_demes; ++i) {
      fscanf(sumry_m,"%lf",&M[i][j]);
    }
  }
  fclose(sumry_m);
  printf("\tparameters M                                    = (");
  fprintf(logfile,"\tparameters M                                    = (");
  for (i = 0; i < data.nbr_demes; ++i) {
    var = postmsqr.M[i] - pow(postmean.M[i],2);
    effective_sample_size = compute_effective_sample_size(M[i],postmean.M[i],var,sample_size);
    if (i < (data.nbr_demes - 1)) {
      printf("%f,",effective_sample_size);
      fprintf(logfile,"%f,",effective_sample_size);
    }
    else {
      printf("%f)\n",effective_sample_size);
      fprintf(logfile,"%f)\n",effective_sample_size);
    }
  }
  for (i = 0; i < data.nbr_demes; ++i) {
    free(M[i]);
  }
  free(M);

  if (!fixed_beta) {
    beta_a = (double *) malloc(sample_size * sizeof(double));
    beta_b = (double *) malloc(sample_size * sizeof(double));
    strcpy(filename,path);
    strcat(filename,TRACE_BETA);
    sumry_beta = fopen(filename,"r");
    while(!((X = getc(sumry_beta)) == '\n' || X == '\f' || X == '\r'));
    for (i = 0; i < sample_size; ++i) {
      fscanf(sumry_beta,"%d",&dummy);
      fscanf(sumry_beta,"%lf",&beta_a[i]);
      fscanf(sumry_beta,"%lf",&beta_b[i]);
    }
    fclose(sumry_beta);
    var = postmsqr.beta_a - pow(postmean.beta_a,2);
    effective_sample_size = compute_effective_sample_size(beta_a,postmean.beta_a,var,sample_size);
    printf("\tshape parameter (alpha) of the parameter pi     = %f\n",effective_sample_size);
    fprintf(logfile,"\tshape parameter (alpha) of the parameter pi     = %f\n",effective_sample_size);
    var = postmsqr.beta_b - pow(postmean.beta_b,2);
    effective_sample_size = compute_effective_sample_size(beta_b,postmean.beta_b,var,sample_size);
    printf("\tshape parameter (beta) of the parameter pi      = %f\n",effective_sample_size);
    fprintf(logfile,"\tshape parameter (beta) of the parameter pi      = %f\n",effective_sample_size);
    free(beta_a);
    free(beta_b);
  }

  if (!fixed_lambda) {
    lambda = (double *) malloc(sample_size * sizeof(double));
    strcpy(filename,path);
    strcat(filename,TRACE_LAMBDA);
    sumry_lambda = fopen(filename,"r");
    while(!((X = getc(sumry_lambda)) == '\n' || X == '\f' || X == '\r'));
    for (i = 0; i < sample_size; ++i) {
      fscanf(sumry_lambda,"%d",&dummy);
      fscanf(sumry_lambda,"%lf",&lambda[i]);
    }
    fclose(sumry_lambda);
    var = postmsqr.lambda - pow(postmean.lambda,2);
    effective_sample_size = compute_effective_sample_size(lambda,postmean.lambda,var,sample_size);
    printf("\t(hyper-)parameter lambda                        = %f\n",effective_sample_size);
    fprintf(logfile,"\t(hyper-)parameter lambda                        = %f\n",effective_sample_size);
    free(lambda);
  }

  printf("\nESS is a measure of how well a Markov chain is mixing. ESS represents the number\n");
  fprintf(logfile,"\nESS is a measure of how well a Markov chain is mixing. ESS represents the number\n");
  printf("of effectively independent draws from the posterior distribution that the Markov\n");
  fprintf(logfile,"of effectively independent draws from the posterior distribution that the Markov\n");
  printf("chain is equivalent to [ESS must be compared to the chain length = %d].\n",sample_size);
  fprintf(logfile,"chain is equivalent to [ESS must be compared to the chain length = %d].\n",sample_size);
  printf("\nWarning! Low ESS (due to strong autocorrelation) indicates poor mixing of the\n");
  fprintf(logfile,"\nWarning! Low ESS (due to strong autocorrelation) indicates poor mixing of the\n");
  printf("Markov chain. The ESS of the (hyper-)parameter lambda is typically lower than that\n");
  fprintf(logfile,"Markov chain. The ESS of the (hyper-)parameter lambda is typically lower than that\n");
  printf("of the other parameters. You are strongly recommended to inspect the trace of the\n");
  fprintf(logfile,"of the other parameters. You are strongly recommended to inspect the trace of the\n");
  printf("lambda parameter in the 'trace_lambda.out' file. The trace shall show relatively\n");
  fprintf(logfile,"lambda parameter in the 'trace_lambda.out' file. The trace shall show relatively\n");
  printf("good mixing (low autocorrelation, AND no decreasing trend). Otherwise, you may want\n");
  fprintf(logfile,"good mixing (low autocorrelation, AND no decreasing trend). Otherwise, you may want\n");
  printf("to increase the length of the burn-in period and/or the total length of the Markov\n");
  fprintf(logfile,"to increase the length of the burn-in period and/or the total length of the Markov\n");
  printf("chain.\n");
  fprintf(logfile,"chain.\n");
  printf("------------------------------------------------------------------------------------\n\n");
  fprintf(logfile,"------------------------------------------------------------------------------------\n\n");
}

double compute_effective_sample_size(double *chain,
                                     double mean,
                                     double var,
                                     int sample_size)

{
  int i,k;
  double effective_sample_size,sum;
  double *autocov,*autocorr;
  
  autocov = (double *) malloc((sample_size - 1) * sizeof(double));
  autocorr = (double *) malloc((sample_size - 1) * sizeof(double));
  for (k = 0; k < (sample_size - 1); ++k) {
    autocov[k] = 0.0;
    for (i = 0; i < (sample_size - k); ++i) {
      autocov[k] += (chain[i] - mean) * (chain[i + k] - mean);
    }
    autocov[k] /= (sample_size - k);
    autocorr[k] = autocov[k] / (var * sample_size / (sample_size - k));
  }
  sum = 0.0;
  k = 1;
  while ((autocorr[k] >= 0.01) && (k < (sample_size - 1))) {
    sum += autocorr[k];
    k += 1;
  }
  effective_sample_size = (double) sample_size / (1.0 + 2.0 * sum);
  free(autocov);
  free(autocorr);
  return(effective_sample_size);
}

char *print_time(double time)

{
  int total;
  int seconds = 0;
  int minutes = 0;
  int hours = 0;
  int days = 0;
  int years = 0;
  char number[8] = "\0";
  static char formatted_time[256] = "\0";
  
  total = (int) time;
  if (total < 60) {
    seconds = total;
    sprintf(number,"%d",seconds);
    strcpy(formatted_time,number);
    strcat(formatted_time," sec");
    (seconds > 1) ? strcat(formatted_time,"s.") : strcat(formatted_time,".");
  }
  else if (total < 3600) {
    minutes = floor(total / 60);
    seconds = total % 60;
    sprintf(number,"%d",minutes);
    strcpy(formatted_time,number);
    strcat(formatted_time," min");
    (minutes > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",seconds);
    strcat(formatted_time,number);
    strcat(formatted_time," sec");
    (seconds > 1) ? strcat(formatted_time,"s.") : strcat(formatted_time,".");
  }
  else if (total < 86400) {
    hours = floor(total / 3600);
    minutes = floor((total - (hours * 3600)) / 60);
    seconds = total % 60;
    sprintf(number,"%d",hours);
    strcpy(formatted_time,number);
    strcat(formatted_time," hr");
    (hours > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",minutes);
    strcat(formatted_time,number);
    strcat(formatted_time," min");
    (minutes > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",seconds);
    strcat(formatted_time,number);
    strcat(formatted_time," sec");
    (seconds > 1) ? strcat(formatted_time,"s.") : strcat(formatted_time,".");
  }
  else if (total < 31536000) {
    days = floor(total / 86400);
    hours = floor((total - (days * 86400)) / 3600);
    minutes = floor((total - (days * 86400) - (hours * 3600)) / 60);
    seconds = total % 60;
    sprintf(number,"%d",days);
    strcpy(formatted_time,number);
    strcat(formatted_time," day");
    (days > 1) ? strcat(formatted_time,"s ") : strcat(formatted_time," ");
    sprintf(number,"%d",hours);
    strcat(formatted_time,number);
    strcat(formatted_time," hr");
    (hours > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",minutes);
    strcat(formatted_time,number);
    strcat(formatted_time," min");
    (minutes > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",seconds);
    strcat(formatted_time,number);
    strcat(formatted_time," sec");
    (seconds > 1) ? strcat(formatted_time,"s.") : strcat(formatted_time,".");
  }
  else {
    years = floor(total / 31536000);
    days = floor((total - (years * 31536000)) / 86400);
    hours = floor((total - (years * 31536000) - (days * 86400)) / 3600);
    minutes = floor((total - (years * 31536000) - (days * 86400) - (hours * 3600)) / 60);
    seconds = total % 60;
    sprintf(number,"%d",years);
    strcpy(formatted_time,number);
    strcat(formatted_time," year");
    (years > 1) ? strcat(formatted_time,"s ") : strcat(formatted_time," ");
    sprintf(number,"%d",days);
    strcat(formatted_time,number);
    strcat(formatted_time," day");
    (days > 1) ? strcat(formatted_time,"s ") : strcat(formatted_time," ");
    sprintf(number,"%d",hours);
    strcat(formatted_time,number);
    strcat(formatted_time," hr");
    (hours > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",minutes);
    strcat(formatted_time,number);
    strcat(formatted_time," min");
    (minutes > 1) ? strcat(formatted_time,"s. ") : strcat(formatted_time,". ");
    sprintf(number,"%d",seconds);
    strcat(formatted_time,number);
    strcat(formatted_time," sec");
    (seconds > 1) ? strcat(formatted_time,"s.\0") : strcat(formatted_time,".");
  }
  strcat(formatted_time,"\0");
  return(formatted_time);
}

void print_usage() {
  printf("usage: %s [ options ]\n",program_name);
  printf("valid options are :\n");
  printf("-help\t\t\t print this message\n");
  printf("-version\t\t print version\n");
  printf("-file\t\t\t name of the input file (default: data.dat)\n");
  printf("-outputs\t\t directory where the outputs will be produced (default: current directory)\n");
  printf("-seed\t\t\t initial seed for the random number generator (default: computed from current time)\n");
  printf("-threads\t\t number of threads to be used (default: number of cpu available)\n");
  printf("-length\t\t\t run length of the Markov chain (default: %ld)\n",chain_length);
  printf("-thin\t\t\t thinning interval size (default: %d)\n",step);
  printf("-burnin\t\t\t length of the burn-in period (default: %ld)\n",burn_in);
  printf("-npilot\t\t\t number of pilot runs (default: %d)\n",n_pilot);
  printf("-lpilot\t\t\t length of each pilot run (default: %d)\n",l_pilot);
  printf("-pool\t\t\t option to analyse data from pooled DNA samples (default: unset)\n");
  printf("-fixed_beta\t\t option to fix the shape parameters of the beta prior distribution of pi (default: unset)\n");
  printf("-beta_a\t\t\t shape parameter of the beta prior distribution of pi (default: %.2f)\n",default_beta_a);
  printf("-beta_b\t\t\t shape parameter of the beta prior distribution of pi (default: %.2f)\n",default_beta_b);
  printf("-fixed_lambda\t\t option to fix the value of lambda (default: unset)\n");
  printf("-lambda_prior\t\t prior distribution of lambda, which can only be inverse gamma ('invgam', by default) or an exponential ('exp')\n");
  printf("-invgam_shape\t\t shape parameter of the inverse gamma prior distribution of lambda (default: %.2f)\n",default_inverse_gamma_shape);
  printf("-invgam_rate\t\t rate parameter of the inverse gamma prior distribution of lambda (default: %.2f)\n",default_inverse_gamma_rate);
  printf("-captl_lambda\t\t rate parameter of the exponential prior distribution of lambda (default: %.2f)\n",default_capital_lambda);
  printf("-min_M\t\t\t lower bound for the log-uniform prior on M (default: %.3f)\n",min_M);
  printf("-max_M\t\t\t upper bound for the log-uniform prior on M (default: %.0f)\n",max_M);
  printf("-max_sig\t\t upper bound for the exponential prior on sigma (default: %.0f)\n",max_sigma);
  printf("-dlt_cnt\t\t half window width from which updates of allele counts are randomly drawn (default: %d)\n",init_delta_counts);
  printf("-dlt_p\t\t\t half window width from which updates of p are randomly drawn (default: %.2f)\n",init_delta_p);
  printf("-dlt_M\t\t\t standard deviation of the lognormal distribution from which updates of M are drawn (default: %.2f) \n",init_delta_M);
  printf("-dlt_pi\t\t\t half window width from which updates of pi are randomly drawn (default: %.2f)\n",init_delta_pi);
  printf("-dlt_sig\t\t standard deviation of the lognormal distribution from which updates of sigma are drawn (default: %.2f)\n",init_delta_sigma);
  printf("-dlt_del\t\t standard deviation of the lognormal distribution from which updates of delta are drawn (default: %.2f)\n",init_delta_delta);
  printf("-dlt_lam\t\t standard deviation of the lognormal distribution from which updates of lambda are drawn (default: %.2f)\n",default_init_delta_lambda);
  printf("-dlt_beta_mu\t\t half window width from which updates of the beta mu parameters are drawn (default: %.2f)\n",default_init_delta_beta_mu);
  printf("-dlt_beta_nu\t\t standard deviation of the lognormal distribution from which updates of the beta nu parameters are drawn (default: %.2f)\n",default_init_delta_beta_nu);
  printf("-calibration\t\t option to generate pseudo-observed data and calibrate the Kullback-Leibler divergence\n");
  printf("-calibration_only\t option to generate pseudo-observed data and calibrate the Kullback-Leibler divergence from previous analyses\n");
  printf("-pod_nbr_loci\t\t option to specify the number of loci to be simulated for calibration (if different from the dataset)\n");
  printf("-verbose\t\t option to print the traces of all parameters (generates big output files!)\n");
  exit(EXIT_SUCCESS);
}

void print_version() {
  printf("You are running version %s\n",VERSION);
}
