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

#ifndef _READ_H
#define _READ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"

FILE *sumry_beta,*sumry_m,*sumry_pi,*sumry_kappa,*sumry_lambda;

void read_data(data_struct *data,char *filename,unsigned int pod_analysis);
void read_output_files(data_struct data,char *path,moments_struct *postmean);
void check_output_files(char *filename,int nbr_columns,int nbr_lines);

#endif
