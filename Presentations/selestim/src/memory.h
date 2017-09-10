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

#ifndef _MEMORY_H
#define _MEMORY_H

#include <stdio.h>
#include <stdlib.h>

#include "defs.h"

void allocate_memory(data_struct data,markov_state_struct *current,markov_state_struct *proposal,updates_struct *update,updates_struct *accept,updates_struct *move,moments_struct *postmean,moments_struct *postmsqr);
void release_memory(data_struct *data,markov_state_struct *current,markov_state_struct *proposal,updates_struct *update,updates_struct *accept,updates_struct *move,moments_struct *postmean,moments_struct *postmsqr);

#endif
