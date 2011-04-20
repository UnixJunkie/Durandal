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

#include "SimpPDB.h"
#include "Stru.h"
#include "rmsd.h"

int _mLen        = 0;
float* _previous = NULL;
double** _mPOS1  = NULL;
double** _mPOS2  = NULL;

// WARNING: this function is not thread safe, not good for OpenMP
float ca_rmsd(float* coor1, float* coor2) {

  if (coor1 == coor2) {
    return 0.0;
  }

  double rmsd;
  int k3;

  if (coor1 == _previous) {
    for (int k = 0; k < _mLen; ++k) {
      k3 = k*3;
      _mPOS2[k][0] = coor2[k3];
      _mPOS2[k][1] = coor2[k3 + 1];
      _mPOS2[k][2] = coor2[k3 + 2];
    }
  } else {
    for (int k = 0; k < _mLen; ++k) {
      k3 = k*3;
      _mPOS1[k][0] = coor1[k3];
      _mPOS1[k][1] = coor1[k3 + 1];
      _mPOS1[k][2] = coor1[k3 + 2];
      _mPOS2[k][0] = coor2[k3];
      _mPOS2[k][1] = coor2[k3 + 1];
      _mPOS2[k][2] = coor2[k3 + 2];
    }
  }
  fast_rmsd(_mPOS1, _mPOS2, _mLen, &rmsd);
  _previous = coor1;

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
      _mLen  = reference->mNumResidue;
      _mPOS1 = new double* [_mLen];
      _mPOS2 = new double* [_mLen];
      for (int i = 0; i < _mLen; ++i) {
        _mPOS1[i] = new double[3];
        _mPOS2[i] = new double[3];
      }
      initialized = true;
    } else {
      SimPDB* sim = new SimPDB(current_line.c_str(), _mLen, false);
      //   float* coor2 = (*_read_pdbs)[j]->mCAlpha;
      cout << current_line << ":"
           << ca_rmsd(reference->mCAlpha, sim->mCAlpha) << '\n';
      delete sim;
    }
  }

  delete reference;
  for (int i = 0; i < _mLen; ++i) {
    delete [] _mPOS1[i];
    delete [] _mPOS2[i];
  }
  delete [] _mPOS1;
  delete [] _mPOS2;

  return 0;
}
