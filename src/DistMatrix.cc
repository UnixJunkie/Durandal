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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <utility> // provides pair

#include "DistMatrix.h"
#include "SimpPDB.h"
#include "Singleton.h"
#include "rmsd.h"

Single& single = Single::instance(); // global var to this file

bool is_shorter(pair<float, size_t> p1, pair<float, size_t> p2);
bool has_less_entropy(pair<int, size_t> p1, pair<int, size_t> p2);

DistMatrix::DistMatrix(const char* input_file, float clustering_radius,
                       clock_t start) {

  _energy_present = false;
  _start = start;
  _brute_delta_nb_undecided = 0;
  _smart_delta_nb_undecided = 0;
  _brute_delta_t = 0;
  _smart_delta_t = 0;
  _nb_references_used = 0;

  _nb_measured = 0;
  _nb_guessed  = 0;

  _DIST_RANGE_ZERO = DistRange(0.0);
  // read all PDB file names
  string current_line;
  SimPDB* firstPDB = NULL;
  bool initialized = false;
  ifstream input_stream(input_file);

  if (not input_stream.is_open()) {
    cout << "error: can't read file: " << string(input_file) << endl;
    exit(1);
  }

  int line_index = 0;
  bool first_time = true;
  bool use_URMSD = single._use_URMSD;
  while (getline(input_stream, current_line)) {

    size_t energy_index = 0;
    string pdb_name;
    char* cstr = new char[current_line.size() + 1];
    strcpy (cstr, current_line.c_str());
    const char* separator = ":";
    char* token = strtok(cstr, separator);

    while (token != NULL) {

      if (energy_index == 0) { // first is PDB filename
        pdb_name = string(token);
        _index_to_pdbname.push_back(pdb_name);
        _index_to_energy.push_back(vector<float>());
      } else { // following are any number of energies
        _energy_present = true;
        _index_to_energy[line_index].push_back(atof(token));
      }
      token = strtok(NULL, separator);
      ++energy_index;
    }
    if (first_time) {
      if (_energy_present) {
        _nb_energies = energy_index - 1; 
      } else {
        _nb_energies = 0;
      }
      first_time = false;
    } else if (_energy_present) {
      if (energy_index - 1 != _nb_energies) {
        cout << "error: line " << line_index
             << ": all energies must be present on each line,"
             << " or never" << endl;
        exit(1);
      }
    }
    if (not initialized) {
      firstPDB = new SimPDB(pdb_name.c_str(), use_URMSD);
      _mLen = firstPDB->mNumResidue;
      _mPOS1 = new double* [_mLen];
      _mPOS2 = new double* [_mLen];
      for (int i = 0; i < _mLen; ++i) {
        _mPOS1[i] = new double[3];
        _mPOS2[i] = new double[3];
      }
      initialized = true;
    }
    ++line_index;
    delete [] cstr;
  }
  input_stream.close();
  _nb_PDBs = _index_to_pdbname.size();
  _nb_pairs = _nb_PDBs * (_nb_PDBs - 1) / 2;
  _nb_undecided = _nb_pairs;
  assert(_nb_PDBs > 0);
  //cout << "nb PDBs read: " << _nb_PDBs << endl; // debug
  _read_pdbs = new vector<Stru*>(_nb_PDBs);
  (*_read_pdbs)[0] = new Stru(firstPDB);
  _all_PDBs_index.push_back(0);
  if (single._only_rank) {
    _matrix = NULL;
  } else {
    _matrix = new DistRange*[_nb_PDBs];
    _matrix[0] = NULL;
  }
  for (size_t i = 1 ; i < _nb_PDBs ; ++i) {
    (*_read_pdbs)[i] = new Stru(new SimPDB(_index_to_pdbname[i].c_str(),
                                           _mLen, use_URMSD));
    _all_PDBs_index.push_back(i);
    // create (N(N-1))/2 elements
    if (not single._only_rank) {
      _matrix[i] = new DistRange[i];
    }
  }
  _references = _all_PDBs_index;
  _clustering_radius = clustering_radius;
  _tabu_FIFO_max_size = 3;
  // at start, each PDB has all but himself nb undecided neighbors
  _nb_undecided_neighbors = vector<int>(_nb_PDBs, _nb_PDBs - 1);
  _nb_ca_rmsd_calls = 0;
  _previous_i_used = numeric_limits<std::size_t>::max();
  if (single._use_URMSD) {
    // N residues ==> N-1 Ca_i to Ca_j vectors on which to do the classic
    //                RMSD computation
    _mLen -= 1;
  }
}

DistMatrix::~DistMatrix() {
  if (single._use_URMSD) {
    // reverse dirty hack for URMSD computation to not miss some mem to free
    _mLen += 1;
  }
  for (size_t i = 0 ; i < _nb_PDBs ; ++i) {
    if (not single._only_rank) {
      delete [] _matrix[i];
    }
    delete (*_read_pdbs)[i];
  }
  for (int i = 0 ; i < _mLen ; ++i) {
    delete [] _mPOS1[i];
    delete [] _mPOS2[i];
  }
  if (not single._only_rank) {
    delete [] _matrix;
  }
  delete _read_pdbs;
  delete [] _mPOS1;
  delete [] _mPOS2;
}

DistRange& DistMatrix::get(size_t i, size_t j) {

  if (i < j) {
    return _matrix[j][i];
  } else {
    if (i == j) {
      return _DIST_RANGE_ZERO;
    } else {
      return _matrix[i][j];
    }
  }
}

float DistMatrix::get_or_measure(size_t i, size_t j) {

  float carmsd;
  DistRange& r = get(i, j);

  if (r._mini != r._maxi) {
    carmsd = ca_rmsd(i, j); // measure it
    --_nb_undecided;
    set(i, j, carmsd);      // store it in cache
  } else {
    carmsd = r._mini;
  }

  return carmsd;
}

void DistMatrix::set(size_t i, size_t j, float exact_distance) {

  if (i < j) {
    _matrix[j][i] = DistRange(exact_distance);
  } else {
    if (i == j) {
      return;
    } else {
      _matrix[i][j] = DistRange(exact_distance);
    }
  }
}

// return exact CARMSD if already stored or compute it if not yet
float DistMatrix::get_or_compute(size_t i, size_t j) {

  float carmsd;
  DistRange& r = get(i, j);

  if (r.is_exact()) {
    carmsd = r._mini;
  } else {
    carmsd = ca_rmsd(i, j);
    set(i, j, carmsd); // store for later
  }

  return carmsd;
}

// D comes from "SPICKER: [...], Zhang & Skolnick, J. Comp. Chem. 2003"
// I just use the cluster center unlike them who are using a virtual
// structure as cluster center (averaged from all cluster members)
//
// WARNING: the cluster center _MUST_ be the first in 'cluster'
float DistMatrix::compute_D(vector<int>& cluster) {

  float result = 0.0;
  size_t end   = cluster.size();

  if (end > 0) {

    size_t center  = cluster[0];
    float sum_rmsd = 0.0;

    for (size_t i = 1 ; i < end ; ++i) {
      sum_rmsd += get_or_compute(center, cluster[i]);
    }
    result = ((float)(end * end)) / ((float)(_nb_PDBs * sum_rmsd));
  }

  return result;
}

void DistMatrix::get_biggest_cluster(const vector<int>& remaining_pdbs,
                                     vector<int>& biggest_cluster_found,
                                     vector< vector<int> >&
                                     pole_position_clusters) {

  static bool first_time = true;
  vector< vector<int> > neighbors_at_d(remaining_pdbs.size());

  // each PDB is a neighbor of itself at distance 0.0
  for (size_t i = 0 ; i < remaining_pdbs.size() ; ++i) {
    neighbors_at_d[i].push_back(remaining_pdbs[i]);
  }
  size_t max_index = 0;
  size_t nb_max_neighbors = 1;
  // find other neighbors
  for (size_t i = 0 ; i < remaining_pdbs.size() ; ++i) {
    for (size_t j = i + 1 ; j < remaining_pdbs.size() ; ++j) {

      size_t fst = remaining_pdbs[i];
      size_t snd = remaining_pdbs[j];
      float carmsd = get(fst, snd)._maxi;

      if (carmsd <= _clustering_radius) {
        neighbors_at_d[i].push_back(snd);
        neighbors_at_d[j].push_back(fst);
        size_t nb_neighbors_at_d_fst = neighbors_at_d[i].size();
        size_t nb_neighbors_at_d_snd = neighbors_at_d[j].size();
        if (nb_neighbors_at_d_fst >= nb_neighbors_at_d_snd) {
          if (nb_neighbors_at_d_fst > nb_max_neighbors) {
            nb_max_neighbors = nb_neighbors_at_d_fst;
            max_index = i;
          }
        } else {
          if (nb_neighbors_at_d_snd > nb_max_neighbors) {
            nb_max_neighbors = nb_neighbors_at_d_snd;
            max_index = j;
          }
        }
      }
    }
  }
  biggest_cluster_found = neighbors_at_d[max_index];
  if (first_time) {
    first_time = false;
    // find cluster centers of same importance
    for (size_t i = 0 ; i < remaining_pdbs.size() ; ++i) {
      if (i != max_index and neighbors_at_d[i].size() == nb_max_neighbors) {
        pole_position_clusters.push_back(neighbors_at_d[i]);
      }
    }
  }
}

// WARNING: this function is not thread safe, not good for OpenMP
float DistMatrix::ca_rmsd(size_t i, size_t j) {

  if (i == j) {
    return 0.0;
  }

  float* coor1 = (*_read_pdbs)[i]->mCAlpha;
  float* coor2 = (*_read_pdbs)[j]->mCAlpha;
  double rmsd;
  int k3;

  ++_nb_measured;

  if (i == _previous_i_used) {
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

  ++_nb_ca_rmsd_calls;

  _previous_i_used = i;
  return (float)rmsd;
}

vector<int> DistMatrix::get_all_PDBs_index() {
  return _all_PDBs_index;
}

// compute all distance pairs
void DistMatrix::brute_init() {

  for (size_t i = 0 ; i < _nb_PDBs ; ++i) {

    for (size_t j = i + 1 ; j < _nb_PDBs ; ++j) {
      float carmsd = ca_rmsd(i, j);
      --_nb_undecided;
      //cout << i << ' ' << j << ' ' << carmsd << endl; // debug
      set(i, j, carmsd);
    }
  }
}

// compute remaining undecided distance pairs
void DistMatrix::brute_finish() {

  for (size_t i = 0 ; i < _nb_PDBs ; ++i) {

    for (size_t j = i + 1 ; j < _nb_PDBs ; ++j) {

      DistRange& range = get(i, j);

      if (not range._is_decided) {

        float lower_bound = 0.0;

        if (not single._no_lower_bound) {
          //cout << "using lower bound" << endl;
          lower_bound = lower_bound_carmsd(i, j);
        }
        if (lower_bound > _clustering_radius) {
          range = DistRange(lower_bound, range._maxi, _clustering_radius);
        } else {

          float upper_bound = numeric_limits<float>::max();

          if (not single._no_upper_bound) {
            //cout << "using upper bound" << endl;
            upper_bound = upper_bound_carmsd(i, j);
          }
          if (upper_bound <= _clustering_radius) {
            range = DistRange(range._mini, upper_bound, _clustering_radius);
          } else {
            set(i, j, ca_rmsd(i, j));
          }
        }
        --_nb_undecided;
      }
    }
  }
}

// compute until the matrix is decidable at _clustering_radius
void DistMatrix::smart_init() {

  float inst_smart_speed = 0.0;
  float avg_brute_speed = 0.0;
  bool switched = false;
  bool can_switch = false;

  while (_nb_undecided > 0) {

    // get distance to new reference for all other PDBs
    // then propagate this new info into the matrix
    fast_propagate(compare_to_reference());

    if (_smart_delta_t > 0 and _brute_delta_t > 0) {
      // could measure all speeds
      inst_smart_speed = ((float)_smart_delta_nb_undecided /
                          (float)_smart_delta_t);
      //cout << "smart_s: " << inst_smart_speed << endl;
      // smart speed is instantaneous, counters need to be reset each time
      _smart_delta_nb_undecided = 0;
      _smart_delta_t = 0;
      avg_brute_speed = ((float)_brute_delta_nb_undecided /
                         (float)_brute_delta_t);
      //cout << "brute_s: " << avg_brute_speed << endl;
      can_switch = true;
    }
    if (single._auto_switch and can_switch and
        (avg_brute_speed >= inst_smart_speed)) {

      int previous_nb_undecided = _nb_pairs;
      switched = true;

      //cout << "switched to brute\n";
      cout << "smart init decided:   "
           << 100*((float)(previous_nb_undecided - _nb_undecided) /
                   (float)_nb_pairs) << " %\n";
      previous_nb_undecided = _nb_undecided;
      brute_finish();
      cout << "brute_finish decided: "
           << 100*((float)(previous_nb_undecided - _nb_undecided) /
                   (float)_nb_pairs) << " %\n";
    }
  }
  if (not switched) {
    cout << "WARNING: did not switch to brute\n";
  }
  cout << "nb ref used: " << _nb_references_used << '\n';
}

size_t DistMatrix::get_next_reference() {

  size_t to_return;
  static bool first_time = true;
  static const size_t MAX_CANDIDATES = 3;
  int candidates [MAX_CANDIDATES];
  float diff_sum [MAX_CANDIDATES];
  float max_diff_sum = 0.0;
  int best_candidate_index = 0;
  int curr_candidate = 0;
  ++_nb_references_used;

  switch (single._reference_choosing_strategy) {
  case MAX_ENTROPY:
    to_return = get_highest_entropy_index();
    break;
  case SEQUENTIAL:
    to_return = pop_back(_references);
    break;
  case RANDOM:
    if (first_time) {
      random_shuffle(_references.begin(), _references.end());
      first_time = false;
    }
    to_return = pop_back(_references);
    break;
  case MAX_DIFF:
    if (first_time) {
      random_shuffle(_references.begin(), _references.end());
      first_time = false;
    }
    // examine a few candidates
    for (size_t l = 0 ; l < MAX_CANDIDATES and _references.size() > 0 ; ++l) {

      curr_candidate = candidates[l] = pop_back(_references);
      //cout << "candidate: " << curr_candidate << endl;
      diff_sum[l] = 0.0;

      if (_tabu_FIFO.size() < _tabu_FIFO_max_size) {
        // fill the tabu FIFO before using it
        best_candidate_index = l;
        break;
      } else {
        // select most different candidate to tabu-stored previous references
        for (size_t m = 0 ; m < _tabu_FIFO.size() ; ++m) {
          float carmsd = get_or_measure(curr_candidate, _tabu_FIFO[m]);
          diff_sum[l] += carmsd;
          //cout << curr_candidate << ' ' << _tabu_FIFO[m] << ' ' << carmsd
          //   << endl;
        }
        //cout << curr_candidate << " diffsum: " << diff_sum[l] << endl;
        if (diff_sum[l] > max_diff_sum) {
          max_diff_sum = diff_sum[l];
          best_candidate_index = l;
        }
      }
    }
    if (_references.size() == 0) {
      cout << "error: no more reference available" << endl;
      exit(1);
    }
    to_return = candidates[best_candidate_index];
    if (_tabu_FIFO.size() == _tabu_FIFO_max_size) {
      // keep _tabu_FIFO up to date and with constant size
      _tabu_FIFO.pop_front();
    }
    _tabu_FIFO.push_back(to_return);
    //cout << "FIFO:" << _tabu_FIFO;
    //cout << "### chose: " << to_return << endl;
    break;
  default:
    cout << "error: unsupported strategy: "
         << single._reference_choosing_strategy << endl;
    exit(1);
  }

  return to_return;
}

size_t DistMatrix::compare_to_reference() {

  clock_t previous = get_user_plus_system_times();
  int previous_nb_undecided = _nb_undecided;
  size_t current_reference = get_next_reference();
  static bool first_time = true;

  for (size_t i = 0 ; i < _nb_PDBs ; ++i) {

    if (i != current_reference) {

      DistRange& range = get(i, current_reference);
      // check if pair was already exactly measured
      if (range._mini != range._maxi) {

        float carmsd;

        if (first_time and not single._no_upper_bound) {
          //       superimposeAndReplace(REF,               MOV); // IMPORTANT
          carmsd = superimposeAndReplace(current_reference, i);
        } else {
          carmsd = ca_rmsd(current_reference, i);
        }

        if (not range._is_decided) {
          --_nb_undecided;
          _nb_undecided_neighbors[i]                 -= 1;
          _nb_undecided_neighbors[current_reference] -= 1;
          ++_brute_delta_nb_undecided;
        }
        set(i, current_reference, carmsd); // store value
      }
    }
  }
  _brute_delta_t += get_user_plus_system_times() - previous;
  _last_brute_delta = previous_nb_undecided - _nb_undecided;
  //cout << "curr_b_delta: " << _last_brute_delta << endl;
  first_time = false;

  return current_reference;
}

bool is_shorter(pair<float, size_t> p1, pair<float, size_t> p2) {
  return p1.first < p2.first;
}

bool has_less_entropy(pair<int, size_t> p1, pair<int, size_t> p2) {
  return p1.first < p2.first;
}

// return the PDB index which has the most entropy
// use to find next reference
size_t DistMatrix::get_highest_entropy_index() {

  int highest_entropy = 0;
  size_t highest_entropy_index = 0;

  for (size_t i = 0 ; i < _nb_PDBs ; ++i) {

    int current_entropy = _nb_undecided_neighbors[i];

    if (current_entropy > highest_entropy) {
      highest_entropy = current_entropy;
      highest_entropy_index = i;
    }
  }

  //cout << _nb_undecided_neighbors << endl;
  //cout << highest_entropy_index << endl;

  return highest_entropy_index;
}

void DistMatrix::update_if_necessary(size_t b, size_t c,
                                     DistRange& old_range,
                                     DistRange& new_range) {

  if ((not old_range._is_decided) and new_range._is_decided) {
    // just became decidable
    old_range = new_range;
    _nb_undecided_neighbors[b] -= 1;
    _nb_undecided_neighbors[c] -= 1;
    ++_nb_guessed;
    --_nb_undecided;
    ++_smart_delta_nb_undecided;
  }
  if (is_sharper(new_range, old_range)) {
    // just became sharper
    old_range = new_range;
  }
}

// propagate exactly measured distances
void DistMatrix::fast_propagate(size_t a) {

  clock_t previous = get_user_plus_system_times();
  int previous_nb_undecided = _nb_undecided;
  vector< pair<float, size_t> > dists_to_ref_a;

  for (size_t i = 0 ; i < _nb_PDBs ; ++i) {
    // make sure i was not a previous reference neither the current one
    if (i != a and _nb_undecided_neighbors[i] > 0) {
      // get dist to ref
      dists_to_ref_a.push_back(pair<float, size_t>(get(a, i)._mini, i));
    }
  }
  size_t end = dists_to_ref_a.size();
  if (end > 0) {
    //cout << "new_ref" << endl;
    //cout << "all: " << dists_to_ref_a << endl;
    sort(dists_to_ref_a.begin(), dists_to_ref_a.end(), is_shorter);
    //cout << "all: " << dists_to_ref_a << endl;
    // look for interesting maximums (AB + AC <= _clustering_radius)
    for (size_t i = 0 ; i < end - 1 ; ++i) { // start from short ABs
      // AB is increasing with i
      size_t b      = dists_to_ref_a[i].second;
      float ab      = dists_to_ref_a[i].first;
      float next_ac = dists_to_ref_a[i+1].first;
      //cout << "ab: " << ab << endl;
      if (ab + next_ac > _clustering_radius) {
        break; // stop increasing AB
      }
      for (size_t j = i + 1 ; j < end ; ++j) { // start from short ACs
        // AC is increasing with j
        size_t c       = dists_to_ref_a[j].second;
        float ac       = dists_to_ref_a[j].first;
        //cout << "ac: " << ac << endl;
        float new_maxi = ab + ac;
        if (new_maxi > _clustering_radius) {
          break; // stop increasing AC
        }
        DistRange& old_range = get(b, c);
        DistRange new_range  = DistRange(old_range._mini, new_maxi,
                                         _clustering_radius);
        update_if_necessary(b, c, old_range, new_range);
      }
    }
    // look for interesting minimums (AC - AB > _clustering_radius)
    float ab_min = dists_to_ref_a[0].first;
    for (size_t j = end - 1 ; j > 0 ; --j) { // start from long ACs
      // AC is decreasing with j
      size_t c = dists_to_ref_a[j].second;
      float ac = dists_to_ref_a[j].first;
      //cout << "ac: " << ac << endl;
      if (ac - ab_min <= _clustering_radius) {
        break; // stop reducing AC
      }
      for (size_t i = 0 ; i < j ; ++i) { // start from short ABs
        // AB is increasing with i
        size_t b = dists_to_ref_a[i].second;
        float ab = dists_to_ref_a[i].first;
        //cout << "ab: " << ab << endl;
        float new_mini = ac - ab;
        if (new_mini <= _clustering_radius) {
          break; // stop increasing AB
        }
        DistRange& old_range = get(b, c);
        DistRange new_range  = DistRange(new_mini, old_range._maxi,
                                         _clustering_radius);
        update_if_necessary(b, c, old_range, new_range);
      }
    }
  }
  _smart_delta_t += get_user_plus_system_times() - previous;
  _last_smart_delta = previous_nb_undecided - _nb_undecided;
  //cout << "curr_s_delta: " << _last_smart_delta << endl;
}

size_t DistMatrix::get_nb_PDBs() {
  return _nb_PDBs;
}

string& DistMatrix::resolve(size_t i) {
  return _index_to_pdbname[i];
}

float DistMatrix::resolve_energy(size_t pdb_index, size_t energy_index) {
  return _index_to_energy[pdb_index][energy_index];
}

float DistMatrix::lower_bound_carmsd(int i, int j) {

  float rev = 0;
  float v;
  Stru& stru1 = *((*_read_pdbs)[i]);
  Stru& stru2 = *((*_read_pdbs)[j]);
  stru1.init_lower_bound_carmsd(_mLen);
  stru2.init_lower_bound_carmsd(_mLen);
  float* sig1 = stru1.mSIG;
  float* sig2 = stru2.mSIG;

  for (int i = 0 ; i < _mLen ; ++i) {
    v = sig1[i] - sig2[i];
    rev += v*v;
  }
  return sqrt(rev / _mLen);
}

float DistMatrix::upper_bound_carmsd(int i, int j) {

  float rev = 0;
  int l3 = _mLen*3;
  float d;
  float* coor1 = (*_read_pdbs)[i]->mCAlpha;
  float* coor2 = (*_read_pdbs)[j]->mCAlpha;

  for (int k = 0 ; k < l3 ; ++k) {
    d = coor1[k] - coor2[k];
    rev += d * d;
  }
  return sqrt(rev / _mLen);
}

// superpose moving to reference, store it, then return CARMSD
float DistMatrix::superimposeAndReplace(size_t reference, size_t moving) {

  int k3;
  double rmsd;
  float* coor1 = (*_read_pdbs)[reference]->mCAlpha;
  float* coor2 = (*_read_pdbs)[moving]->mCAlpha;

  if (reference == _previous_i_used) {
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
  RMSDAndReplace(_mPOS1, _mPOS2, _mLen, &rmsd, coor2);

  ++_nb_ca_rmsd_calls;

  _previous_i_used = reference;
  return (float)rmsd;
}

float DistMatrix::find_clustering_distance(float cutoff_ratio,
                                           float bin_width,
                                           float* histo_median) {

  if (cutoff_ratio >= 1.0) {
    cout << "error: invalid cutoff_ratio: " << cutoff_ratio << endl;
    cout << "(must be < 1.0)" << endl;
    exit(1);
  }
  if (cutoff_ratio == 0.0) {
    cutoff_ratio = 0.05;
    cout << "find_clustering_distance: using default cutoff ratio: "
         << cutoff_ratio << endl;
  }

  float result = 0.0;
  vector<float> sampled_distances;

  cout << "sampling distances..." << endl;
  *histo_median = sample(sampled_distances);
  //cout << "dists" << endl << sampled_distances << endl << "/dists" << endl;
  int end = sampled_distances.size();
  result = sampled_distances[(int)(cutoff_ratio * (end - 1))];

  if (bin_width > 0.0) {

    // compute histogram
    float offset    = 0.0;
    float max_val   = 0.0;
    float delta     = 0.0;
    int nb_bins     = 0;
    int current_bin = 0;

    if (end > 0) {
      offset  = sampled_distances[0];
      max_val = sampled_distances[end - 1];
      delta   = max_val - offset;
      nb_bins = (int)ceil(delta / bin_width);
    }
    if (nb_bins > 0) {
      vector< pair<float, int> > bins = vector< pair<float, int> >(nb_bins);
      for (int k = 0 ; k < nb_bins ; ++k) {
        bins[k].first = (k * delta) / nb_bins + offset;
      }
      for (int k = 0 ; k < end ; ++k) {
        current_bin = (int)((sampled_distances[k] - offset) / bin_width);
        bins[current_bin].second += 1;
      }
      cerr << bins;
    }
  }

  return result;
}

// sample pairwise distances until median of the accumulated sample
// is stable enough, sampled distances are put inside sampled_distances
// the median of the sampled set is returned
float DistMatrix::sample(vector<float>& sampled_distances) {

  int sample_size      = 100;
  float last_median    = 0.0;
  float current_median = 100.0;
  size_t i, j;

  // as long as median is not stabilized
  do {
    // sample more
    for (int k = 0 ; k < sample_size ; ++k) {

      i = rand_between(0, _nb_PDBs);
      j = rand_between(0, _nb_PDBs);
      DistRange& r = get(i, j);

      if (r.is_exact()) { // was already measured
        sampled_distances.push_back(r._mini);
      } else {            // not yet measured
        float carmsd = ca_rmsd(i, j);
        --_nb_undecided;
        _nb_undecided_neighbors[i] -= 1;
        _nb_undecided_neighbors[j] -= 1;
        set(i, j, carmsd);
        sampled_distances.push_back(carmsd);
      }
    }
    sort(sampled_distances);
    last_median = current_median;
    current_median = median(sampled_distances);
    cout << " - using " << sample_size << " samples" << endl;
  } while (abs(current_median - last_median) > 0.025 * last_median);

  return current_median;
}

void DistMatrix::set_clustering_distance(float d) {
  _clustering_radius = d;
}

// use this to output RMSD to the correct structure, by putting
// the correct structure as the first in the input PDB list
void DistMatrix::only_rank() {

  for (size_t i = 1 ; i < _nb_PDBs ; ++i) {
    cerr << resolve(i) << ':' << ca_rmsd(0, i) << endl;
  }
}

bool has_lower_avg_RMSD(pair<double, size_t> p1, pair<double, size_t> p2) {
  return p1.first < p2.first;
}

void DistMatrix::display_energy(ostream& out, vector<int>& cluster) {

  vector< vector<float> > energies(_nb_energies);

  for (size_t j = 0 ; j < _nb_energies ; ++j) {

    energies.push_back(vector<float>());

    for (size_t i = 0 ; i < cluster.size() ; ++i) {
      energies[j].push_back(_index_to_energy[cluster[i]][j]);
    }
    sort(energies[j]);
    out << "min: "  << energies[j][0]
        << " max: " << energies[j][cluster.size() - 1]
        << " avg: " << average(energies[j])
        << " med: " << median(energies[j]) << endl;
  }
}

bool DistMatrix::has_energies() {
  return _energy_present;
}

size_t DistMatrix::get_nb_energies() {
  return _nb_energies;
}
