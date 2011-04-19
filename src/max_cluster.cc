/**
 *  Copyright (C) 2010, 2011 Zhang Initiative Research Unit,
 *  Advance Science Institute, Riken
 *  2-1 Hirosawa, Wako, Saitama 351-0198, Japan
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *****************************************************************************/

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "DistRange.h"

int main (int argc, char** argv) {

  if (argc != 3) {
    cout << "Given the number of residues of your structure and the size of\n"
         << "memory (in bytes) that you are willing to use for clustering:\n"
         << "output the MAXIMUM NUMBER OF DECOYS YOU SHOULD TRY TO CLUSTER\n"
         << "---" << endl;
    //              0             1           2
    cout << "usage: ./max_cluster nb_residues max_mem_to_use_in_bytes" << endl;
    return 1;
  }

  float nb_res  = atof(argv[1]);
  float max_mem = atof(argv[2]);
  float dr_size = (float) sizeof(DistRange);
  float f_size  = (float) sizeof(float);
  float a_coeff = dr_size / 2;
  float b_coeff = -dr_size / 2 + nb_res * f_size;
  float c_coeff = -max_mem;
  printf("nb max PDB to cluster ~= %.2f\n",
         ((-b_coeff + sqrt(b_coeff * b_coeff - 4 * a_coeff * c_coeff)) /
          (2 * a_coeff)));
  return 0;
}
