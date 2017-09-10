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

#ifndef _WRITE_H
#define _WRITE_H

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dcmt0.6.1b/include/dc.h"
#include "defs.h"
#include "mymath.h"
#include "rand.h"

FILE *outputs_mcmc,*trace_m,*trace_pi,*trace_sigma,*trace_delta,*trace_lambda,*trace_kappa,*trace_beta,*trace_freq,*trace_counts;
FILE *sumry_m,*sumry_pi,*sumry_kappa,*sumry_sigma,*sumry_delta,*sumry_lambda,*sumry_beta,*sumry_counts,*sumry_freq;
FILE *kld_qutile;

void open_log_file(char *path);
void open_output_files(char *path);
void close_output_files(void);
void write_outputs(data_struct data,markov_state_struct current,updates_struct update,updates_struct accept,long iter);
void write_summaries(data_struct data,moments_struct postmean,moments_struct postmsqr,char *path,unsigned int pod_analysis);
void write_headers(data_struct data);
double draw_from_psi(double sigma,double M,double pi);
void write_pseudo_observed_data(data_struct data,moments_struct postmean,char *pod_path,char *filename,char *pod_file_name);
void copy_data(data_struct data,char *extended_filename,char *path);
void write_readme(data_struct data,char *file_name,char *path,char *pod_file_name);
void print_model_parameters(unsigned long seed,int n_threads);
int compare(const void *a, const void *b);
void print_F_ST(data_struct data);
void compute_autocorrelation(moments_struct postmean,moments_struct postmsqr,char *path);
void print_effective_sample_size(data_struct data,moments_struct postmean,moments_struct postmsqr,char *path);
double compute_effective_sample_size(double *chain,double mean,double var,int sample_size);
void print_usage();
char *print_time(double time);
void print_version();

#endif

