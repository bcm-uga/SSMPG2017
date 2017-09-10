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

#include "read.h"

void read_data(data_struct *data,
               char *filename,
               unsigned int pod_analysis)

{
	FILE *infile = NULL;
	char X;
	int i,j,k;
  int count,dummy,nloci,npops;
  int min,max,sum_counts[2],sum_reads[2];
  double mean;
  
	if ((infile = fopen(filename,"r")) == NULL) {
    printf("%s: file not found\n",filename);
    fprintf(logfile,"%s: file not found\n",filename);
    exit(EXIT_FAILURE);
	}
  if (!pod_analysis && !calibration_only) {
    printf("Checking file `%s'... ",filename);
    fprintf(logfile,"Checking file `%s'... ",filename);
  }
  data -> nbr_demes = data -> nbr_loci = 0;
	fscanf(infile,"%d",&data -> nbr_demes);
	while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));
	fscanf(infile,"%d",&data -> nbr_loci);
	while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));
  if (pooled_data) {
    npops = 0;
    do {
      dummy = -1;
      if ((X = fscanf(infile,"%i%n",&count,&dummy)) != EOF) {
        if (dummy == -1) {
          printf("\n...unexpected character in line 3... Please check the input file\n");
          fprintf(logfile,"\n...unexpected character in line 3... Please check the input file\n");
          exit(EXIT_FAILURE);
        }
      }
      npops++;
      do {
        X = getc(infile);
        if (X == '\n' || X == '\f' || X == '\r') {
          if (npops > data -> nbr_demes) {
            printf("\n...the number of pools in line 3 is larger than expected... Please check the input file\n");
            fprintf(logfile,"\n...the number of pools in line 3 is larger than expected... Please check the input file\n");
            exit(EXIT_FAILURE);
          }
          if (npops < data -> nbr_demes) {
            printf("\n...the number of pools in line 3 is lower than expected... Please check the input file\n");
            fprintf(logfile,"\n...the number of pools in line 3 is lower than expected... Please check the input file\n");
            exit(EXIT_FAILURE);
          }
          if (npops == data -> nbr_demes) {
            break;
          }
        }
      } while (X  == '\t' || X == ' ');
      fseek(infile, -1, SEEK_CUR);
    } while (X != '\n' && X != '\f' && X != '\r');
    data -> pool_size = (int *) malloc(data -> nbr_demes * sizeof(int));
   	data -> total_nbr_reads = (int **) malloc(data -> nbr_demes * sizeof(int *));
    data -> reads = (int ***) malloc(data -> nbr_demes * sizeof(int **));
    min_counts = (int ***) malloc(data -> nbr_demes * sizeof(int **));
    max_counts = (int ***) malloc(data -> nbr_demes * sizeof(int **));
    for (i = 0; i < data -> nbr_demes; ++i) {
      data -> total_nbr_reads[i]=(int *) malloc(data -> nbr_loci * sizeof(int));
      data -> reads[i]=(int **) malloc(data -> nbr_loci * sizeof(int *));
      min_counts[i] = (int **) malloc(data -> nbr_loci * sizeof(int *));
      max_counts[i] = (int **) malloc(data -> nbr_loci * sizeof(int *));
      for (j = 0; j < data -> nbr_loci; ++j) {
        data -> total_nbr_reads[i][j] = 0;
        data -> reads[i][j] = (int *) malloc(2 * sizeof(int));
        data -> reads[i][j][0] = data -> reads[i][j][1] = 0;
        min_counts[i][j] = (int *) malloc(2 * sizeof(int));
        max_counts[i][j] = (int *) malloc(2 * sizeof(int));
      }
    }
  }
	data -> total_nbr_counts = (int **) malloc(data -> nbr_demes * sizeof(int *));
	data -> counts = (int ***) malloc(data -> nbr_demes * sizeof(int **));
  for (i = 0; i < data -> nbr_demes; ++i) {
		data -> total_nbr_counts[i]=(int *) malloc(data -> nbr_loci * sizeof(int));
		data -> counts[i]=(int **) malloc(data -> nbr_loci * sizeof(int *));
    for (j = 0; j < data -> nbr_loci; ++j) {
      data -> total_nbr_counts[i][j] = 0;
      data -> counts[i][j] = (int *) malloc(2 * sizeof(int));
      data -> counts[i][j][0] = data -> counts[i][j][1] = 0;
    }
	}
  nloci = 0;
  do {
    npops = 0;
    do {
      dummy = -1;
      if ((X = fscanf(infile,"%i%n",&count,&dummy)) != EOF) {
        if (dummy == -1) {
          printf("\n...unexpected character at marker no. %d... Please check the input file\n",(nloci + 1));
          fprintf(logfile,"\n...unexpected character at marker no. %d... Please check the input file\n",(nloci + 1));
          exit(EXIT_FAILURE);
        }
      }
      npops++;
      do {
        X = getc(infile);
        if (X == '\n' || X == '\f' || X == '\r') {
          nloci++;
          if (npops > (2 * data -> nbr_demes)) {
            printf("\n...the number of SNP counts at marker no. %d is larger than expected... Please check the input file\n",nloci);
            fprintf(logfile,"\n...the number of SNP counts at marker no. %d is larger than expected... Please check the input file\n",nloci);
            exit(EXIT_FAILURE);
          }
          if (npops < (2 * data -> nbr_demes)) {
            printf("\n...the number of SNP counts at marker no. %d is lower than expected... Please check the input file\n",nloci);
            fprintf(logfile,"\n...the number of SNP counts at marker no. %d is lower than expected... Please check the input file\n",nloci);
            exit(EXIT_FAILURE);
          }
          if (npops == (2 * data -> nbr_demes)) {
            break;
          }
        }
      } while (X  == '\t' || X == ' ');
      fseek(infile, -1, SEEK_CUR);
    } while (X != '\n' && X != '\f' && X != '\r' && X != EOF);
  } while (X != EOF);
  if (nloci < data -> nbr_loci) {
    printf("\n...the number of loci (%d) is lower than expected (%d)... Please check the input file\n",nloci,data -> nbr_loci);
    fprintf(logfile,"\n...the number of loci (%d) is lower than expected (%d)... Please check the input file\n",nloci,data -> nbr_loci);
    exit(EXIT_FAILURE);
  }
  if (nloci > data -> nbr_loci) {
    printf("\n...the number of loci (%d) is larger than expected (%d)... Please check the input file\n",nloci,data -> nbr_loci);
    fprintf(logfile,"\n...the number of loci (%d) is larger than expected (%d)... Please check the input file\n",nloci,data -> nbr_loci);
    exit(EXIT_FAILURE);
  }
  rewind(infile);
	fscanf(infile,"%d",&data -> nbr_demes);
	while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));
	fscanf(infile,"%d",&data -> nbr_loci);
	while(!((X = getc(infile)) == '\n' || X == '\f' || X == '\r'));
  if (pooled_data) {
    for (i = 0; i < data -> nbr_demes; ++i) {
      fscanf(infile,"%d",&data -> pool_size[i]);
    }
  }
  if (pooled_data) {
    for (j = 0; j < data -> nbr_loci; ++j) {
      for (i = 0; i < data -> nbr_demes; ++i) {
        for (k = 0; k < 2; ++k) {
          fscanf(infile,"%d",&data -> reads[i][j][k]);
        }
      }
      for (i = 0; i < data -> nbr_demes; ++i) {
        data -> total_nbr_reads[i][j] = 0;
        for (k = 0; k < 2; ++k) {
          data -> total_nbr_reads[i][j] += data -> reads[i][j][k];
        }
        if (data -> total_nbr_reads[i][j] > 0) {
          data -> total_nbr_counts[i][j] = data -> pool_size[i];
        }
        else {
          data -> total_nbr_counts[i][j] = 0;
        }
      }
    }
    for (i = 0; i < data -> nbr_demes; ++i) {
      for (j = 0; j < data -> nbr_loci; ++j) {
        if (data -> total_nbr_counts[i][j] > 0) {
          for (k = 0; k < 2; ++k) {
            if ((data -> reads[i][j][k] > 0) && (data -> reads[i][j][k] < data -> total_nbr_reads[i][j])) {
              min_counts[i][j][k] = 1;
              max_counts[i][j][k] = data -> total_nbr_counts[i][j] - 1;
            }
            else if (data -> reads[i][j][k] == 0) {
              min_counts[i][j][k] = 0;
              max_counts[i][j][k] = data -> total_nbr_counts[i][j] - 1;
            }
            else if (data -> reads[i][j][k] == data -> total_nbr_reads[i][j]) {
              min_counts[i][j][k] = 1;
              max_counts[i][j][k] = data -> total_nbr_counts[i][j];
            }
          }
        }
        else {
          for (k = 0; k < 2; ++k) {
            min_counts[i][j][k] = 0;
            max_counts[i][j][k] = 0;
          }
        }
      }
    }
    min = 2147483647;
    for (j = 0; j < data -> nbr_loci; ++j) {
      sum_reads[0] = sum_reads[1] = 0;
      for (i = 0; i < data -> nbr_demes; ++i) {
        for (k = 0; k < 2; ++k) {
          sum_reads[k] += data -> reads[i][j][k];
        }
      }
      if (sum_reads[0] < min) min = sum_reads[0];
      if (sum_reads[1] < min) min = sum_reads[1];
      if ((sum_reads[0] == 0) && (sum_reads[1] == 0)) {
        printf("\nThere is absolutely no data at locus %d (line %d)... Please remove that locus from the input file\n",(j + 1),(j + 4));
        fprintf(logfile,"\nThere is absolutely no data at locus %d (line %d)... Please remove that locus from the input file\n",(j + 1),(j + 4));
        exit(EXIT_FAILURE);
      }
    }
    min_nbr_reads = min;
  }
  else {
    for (j = 0; j < data -> nbr_loci; ++j) {
      for (i = 0; i < data -> nbr_demes; ++i) {
        for (k = 0; k < 2; ++k) {
          fscanf(infile,"%d",&data -> counts[i][j][k]);
        }
      }
      for (i = 0; i < data -> nbr_demes; ++i) {
        data -> total_nbr_counts[i][j] = 0;
        for (k = 0; k < 2; ++k) {
          data -> total_nbr_counts[i][j] += data -> counts[i][j][k];
        }
      }
    }
    min = 2147483647;
    for (j = 0; j < data -> nbr_loci; ++j) {
      sum_counts[0] = sum_counts[1] = 0;
      for (i = 0; i < data -> nbr_demes; ++i) {
        for (k = 0; k < 2; ++k) {
          sum_counts[k] += data -> counts[i][j][k];
        }
      }
      if (sum_counts[0] < min) min = sum_counts[0];
      if (sum_counts[1] < min) min = sum_counts[1];
      if ((sum_counts[0] == 0) && (sum_counts[1] == 0)) {
        printf("\nThere is absolutely no data at locus %d (line %d)... Please remove that locus from the input file\n",(j + 1),(j + 3));
        fprintf(logfile,"\nThere is absolutely no data at locus %d (line %d)... Please remove that locus from the input file\n",(j + 1),(j + 3));
        exit(EXIT_FAILURE);
      }
    }
    min_nbr_counts = min;
  }
  if (!pod_analysis && !calibration_only) {
    printf("OK\n");
    fprintf(logfile,"OK\n");
  }
  fclose(infile);
  if (!pod_analysis && !calibration_only) {
    printf("The data consist in %d SNPs and %d sampled populations\n\n",data -> nbr_loci,data -> nbr_demes);
    fprintf(logfile,"The data consist in %d SNPs and %d sampled populations\n\n",data -> nbr_loci,data -> nbr_demes);
    if (pooled_data) {
      printf("------------------------------------------------------------------------------------\n");
      fprintf(logfile,"------------------------------------------------------------------------------------\n");
      printf("Mean number of reads (min, max) per sampled population:\n");
      fprintf(logfile,"Mean number of reads (min, max) per sampled population:\n");
      printf("------------------------------------------------------------------------------------\n");
      fprintf(logfile,"------------------------------------------------------------------------------------\n");
      for (i = 0; i < data -> nbr_demes; ++i) {
        min = 2147483647;
        max = -2147483648;
        mean = 0.0;
        for (j = 0; j < data -> nbr_loci; ++j) {
          if (data -> total_nbr_reads[i][j] > max) max = data -> total_nbr_reads[i][j];
          if (data -> total_nbr_reads[i][j] < min) min = data -> total_nbr_reads[i][j];
          mean += (double) data -> total_nbr_reads[i][j];
        }
        mean /= data -> nbr_loci;
        printf("Population no. %2d: %.2f (%d,%d)\n",(i + 1),mean,min,max);
        fprintf(logfile,"Population no. %2d: %.2f (%d,%d)\n",(i + 1),mean,min,max);
      }
      min = 2147483647;
      max = -2147483648;
      mean = 0.0;
      for (j = 0; j < data -> nbr_loci; ++j) {
        for (i = 0; i < data -> nbr_demes; ++i) {
          if (data -> total_nbr_reads[i][j] > max) max = data -> total_nbr_reads[i][j];
          if (data -> total_nbr_reads[i][j] < min) min = data -> total_nbr_reads[i][j];
          mean += (double) data -> total_nbr_reads[i][j];
        }
      }
      mean /= (data -> nbr_loci * data -> nbr_demes);
      printf("------------------------------------------------------------------------------------\n");
      fprintf(logfile,"------------------------------------------------------------------------------------\n");
      printf("Overall          : %.2f (%d,%d)\n",mean,min,max);
      fprintf(logfile,"Overall          : %.2f (%d,%d)\n",mean,min,max);
      printf("------------------------------------------------------------------------------------\n\n");
      fprintf(logfile,"------------------------------------------------------------------------------------\n\n");
    }
    else {
      printf("------------------------------------------------------------------------------------\n");
      fprintf(logfile,"------------------------------------------------------------------------------------\n");
      printf("Mean sample size (min, max) per sampled population:\n");
      fprintf(logfile,"Mean sample size (min, max) per sampled population:\n");
      printf("------------------------------------------------------------------------------------\n");
      fprintf(logfile,"------------------------------------------------------------------------------------\n");
      for (i = 0; i < data -> nbr_demes; ++i) {
        min = 2147483647;
        max = -2147483648;
        mean = 0.0;
        for (j = 0; j < data -> nbr_loci; ++j) {
          if (data -> total_nbr_counts[i][j] > max) max = data -> total_nbr_counts[i][j];
          if (data -> total_nbr_counts[i][j] < min) min = data -> total_nbr_counts[i][j];
          mean += (double) data -> total_nbr_counts[i][j];
        }
        mean /= data -> nbr_loci;
        printf("Population no. %2d: %.2f (%d,%d)\n",(i + 1),mean,min,max);
        fprintf(logfile,"Population no. %2d: %.2f (%d,%d)\n",(i + 1),mean,min,max);
      }
      min = 2147483647;
      max = -2147483648;
      mean = 0.0;
      for (j = 0; j < data -> nbr_loci; ++j) {
        for (i = 0; i < data -> nbr_demes; ++i) {
          if (data -> total_nbr_counts[i][j] > max) max = data -> total_nbr_counts[i][j];
          if (data -> total_nbr_counts[i][j] < min) min = data -> total_nbr_counts[i][j];
          mean += (double) data -> total_nbr_counts[i][j];
        }
      }
      mean /= (data -> nbr_loci * data -> nbr_demes);
      printf("------------------------------------------------------------------------------------\n");
      fprintf(logfile,"------------------------------------------------------------------------------------\n");
      printf("Overall          : %.2f (%d,%d)\n",mean,min,max);
      fprintf(logfile,"Overall          : %.2f (%d,%d)\n",mean,min,max);
      printf("------------------------------------------------------------------------------------\n\n");
      fprintf(logfile,"------------------------------------------------------------------------------------\n\n");
    }
  }
}

void read_output_files(data_struct data,
                       char *path,
                       moments_struct *postmean)

{
  int i,j;
  int max_nbr_columns,max_nbr_lines;
  double dummy;
  char X;
  char filename[150] = "";

  if(!fixed_beta) {
    strcpy(filename,path);
    strcat(filename,SUMRY_BETA);
    max_nbr_columns = 2;
    max_nbr_lines = 3;
    check_output_files(filename,max_nbr_columns,max_nbr_lines);
    sumry_beta = fopen(filename,"r");
    while(!((X = getc(sumry_beta)) == '\n' || X == '\f' || X == '\r'));
    fscanf(sumry_beta,"%lf",&(postmean -> beta_a));
    fscanf(sumry_beta,"%lf",&dummy);
    fscanf(sumry_beta,"%lf",&(postmean -> beta_b));
    fscanf(sumry_beta,"%lf",&dummy);
    fclose(sumry_beta);
  }
  
  strcpy(filename,path);
  strcat(filename,SUMRY_LAMBDA);
  max_nbr_columns = 2;
  max_nbr_lines = 2;
  check_output_files(filename,max_nbr_columns,max_nbr_lines);
  sumry_lambda = fopen(filename,"r");
  while(!((X = getc(sumry_lambda)) == '\n' || X == '\f' || X == '\r'));
  fscanf(sumry_lambda,"%lf",&(postmean -> lambda));
  fscanf(sumry_lambda,"%lf",&dummy);
  fclose(sumry_lambda);
  
  strcpy(filename,path);
  strcat(filename,SUMRY_KAPPA);
  max_nbr_columns = data.nbr_demes + 1;
  max_nbr_lines = data.nbr_loci + 1;
  check_output_files(filename,max_nbr_columns,max_nbr_lines);
  sumry_kappa = fopen(filename,"r");
  while(!((X = getc(sumry_kappa)) == '\n' || X == '\f' || X == '\r'));
  for (j = 0; j < data.nbr_loci; ++j) {
    fscanf(sumry_kappa,"%lf",&dummy);
    for (i = 0; i < data.nbr_demes; ++i) {
      fscanf(sumry_kappa,"%lf",&(postmean -> kappa[i][j]));
    }
  }
  fclose(sumry_kappa);
  
  strcpy(filename,path);
  strcat(filename,SUMRY_PI);
  max_nbr_columns = 3;
  max_nbr_lines = data.nbr_loci + 1;
  check_output_files(filename,max_nbr_columns,max_nbr_lines);
  sumry_pi = fopen(filename,"r");
  while(!((X = getc(sumry_pi)) == '\n' || X == '\f' || X == '\r'));
  for (j = 0; j < data.nbr_loci; ++j) {
    fscanf(sumry_pi,"%lf",&dummy);
    fscanf(sumry_pi,"%lf",&(postmean -> pi[j]));
    fscanf(sumry_pi,"%lf",&dummy);
  }
  fclose(sumry_pi);

  strcpy(filename,path);
  strcat(filename,SUMRY_M);
  max_nbr_columns = 3;
  max_nbr_lines = data.nbr_demes + 1;
  check_output_files(filename,max_nbr_columns,max_nbr_lines);
  sumry_m = fopen(filename,"r");
  while(!((X = getc(sumry_m)) == '\n' || X == '\f' || X == '\r'));
  for (i = 0; i < data.nbr_demes; ++i) {
    fscanf(sumry_m,"%lf",&dummy);
    fscanf(sumry_m,"%lf",&(postmean -> M[i]));
    fscanf(sumry_m,"%lf",&dummy);
  }
  fclose(sumry_m);
}

void check_output_files(char *filename,
                        int max_nbr_columns,
                        int max_nbr_lines)

{
  int nbr_columns,nbr_lines;
  double dummy;
	char X;
  FILE *file_to_check;
  
  if ((file_to_check = fopen(filename,"r")) == NULL) {
    printf("Cannot open %s for the calibration of the Kullback-Leibler divergence!\n",filename);
    exit(EXIT_FAILURE);
  }
  while(!((X = getc(file_to_check)) == '\n' || X == '\f' || X == '\r'));
  nbr_lines = 1;
  do {
    nbr_columns = 0;
    do {
      fscanf(file_to_check,"%lf",&dummy);
      nbr_columns++;
      do {
        X = getc(file_to_check);
        if (X == '\n' || X == '\f' || X == '\r') {
          nbr_lines++;
          if (nbr_columns > max_nbr_columns) {
            printf("%s: the number of columns at line no. %d is larger than expected... Cannot calibrate the Kullback-Leibler divergence\n",filename,nbr_lines);
            fprintf(logfile,"%s: the number of columns at line no. %d is larger than expected... Cannot calibrate the Kullback-Leibler divergence\n",filename,nbr_lines);
            exit(EXIT_FAILURE);
          }
          if (nbr_columns < max_nbr_columns) {
            printf("%s: the number of columns at line no. %d is lower than expected... Cannot calibrate the Kullback-Leibler divergence\n",filename,nbr_lines);
            fprintf(logfile,"%s: the number of columns at line no. %d is lower than expected... Cannot calibrate the Kullback-Leibler divergence\n",filename,nbr_lines);
            exit(EXIT_FAILURE);
          }
          if (nbr_columns == max_nbr_columns) {
            break;
          }
        }
      } while (X  == '\t' || X == ' ');
      fseek(file_to_check, -1, SEEK_CUR);
    } while (X != '\n' && X != '\f' && X != '\r' && X != EOF);
  } while (X != EOF);
  if (nbr_lines < max_nbr_lines) {
    printf("%s: the number of lines (%d) is lower than expected (%d)... Cannot calibrate the Kullback-Leibler divergence\n",filename,nbr_lines,max_nbr_lines);
    fprintf(logfile,"%s: the number of lines (%d) is lower than expected (%d)... Cannot calibrate the Kullback-Leibler divergence\n",filename,nbr_lines,max_nbr_lines);
    exit(EXIT_FAILURE);
  }
  if (nbr_lines > max_nbr_lines) {
    printf("%s: the number of lines (%d) is larger than expected (%d)... Cannot calibrate the Kullback-Leibler divergence\n",filename,nbr_lines,max_nbr_lines);
    fprintf(logfile,"%s: the number of lines (%d) is larger than expected (%d)... Cannot calibrate the Kullback-Leibler divergence\n",filename,nbr_lines,max_nbr_lines);
    exit(EXIT_FAILURE);
  }
  fclose(file_to_check);
}


