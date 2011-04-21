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

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "SimpPDB.h"
#include "Stru.h"
#include "rmsd.h"

int mLen        = 0;
float* previous = NULL;
double** mPOS1  = NULL;
double** mPOS2  = NULL;

// WARNING: this function is not thread safe, not good for OpenMP
float ca_rmsd(float* coor1, float* coor2) {

  if (coor1 == coor2) {
    return 0.0;
  }

  double rmsd;
  int k3;

  if (coor1 == previous) {
    for (int k = 0; k < mLen; ++k) {
      k3 = k*3;
      mPOS2[k][0] = coor2[k3];
      mPOS2[k][1] = coor2[k3 + 1];
      mPOS2[k][2] = coor2[k3 + 2];
    }
  } else {
    for (int k = 0; k < mLen; ++k) {
      k3 = k*3;
      mPOS1[k][0] = coor1[k3];
      mPOS1[k][1] = coor1[k3 + 1];
      mPOS1[k][2] = coor1[k3 + 2];
      mPOS2[k][0] = coor2[k3];
      mPOS2[k][1] = coor2[k3 + 1];
      mPOS2[k][2] = coor2[k3 + 2];
    }
  }
  fast_rmsd(mPOS1, mPOS2, mLen, &rmsd);
  previous = coor1;

  return (float)rmsd;
}

int main (int argc, char** argv) {

  if (argc != 2) {
    cout << "Output on stdout the CARMSD of each PDB in the list to the\n"
         << "first one in the list (for example, you may want to put as\n"
         << "first PDB the one for the known structure).\n"
         << "---" << endl;
    //                   0            1
    cout << "usage: " << argv[0] << " PDB_list" << endl;
    return 1;
  }

  ifstream input_stream(argv[1]);
  string current_line;

  if (not input_stream.is_open()) {
    cout << "error: can't read file: " << string(argv[1]) << endl;
    exit(1);
  }

  bool initialized  = false;
  SimPDB* reference = NULL;

  while (getline(input_stream, current_line)) {

    if (not initialized) {
      reference = new SimPDB(current_line.c_str(), false);
      mLen  = reference->mNumResidue;
      mPOS1 = new double* [mLen];
      mPOS2 = new double* [mLen];
      for (int i = 0; i < mLen; ++i) {
        mPOS1[i] = new double[3];
        mPOS2[i] = new double[3];
      }
      initialized = true;
    } else {
      SimPDB* sim = new SimPDB(current_line.c_str(), mLen, false);
      //   float* coor2 = (*_read_pdbs)[j]->mCAlpha;
      cout << current_line << ":"
           << ca_rmsd(reference->mCAlpha, sim->mCAlpha) << '\n';
      delete sim;
    }
  }

  delete reference;
  for (int i = 0; i < mLen; ++i) {
    delete [] mPOS1[i];
    delete [] mPOS2[i];
  }
  delete [] mPOS1;
  delete [] mPOS2;

  return 0;
}
