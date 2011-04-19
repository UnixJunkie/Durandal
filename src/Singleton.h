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

#ifndef SINGLETON_H
#define SINGLETON_H

#include <cstdio>
#include <deque>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sys/times.h>

using namespace std;

enum Strategy {
  MAX_DIFF,
  MAX_ENTROPY,
  RANDOM,
  SEQUENTIAL
};

class Single {

 public: // ------------------------------------------------------------------

  // some attributes are public for easier access 
  bool _brute_mode;
  bool _smart_mode;
  bool _verbose;
  bool _auto_switch;
  bool _no_lower_bound;
  bool _no_upper_bound;
  bool _output_to_file;
  bool _compute_D;
  bool _only_rank;
  bool _use_URMSD;
  Strategy _reference_choosing_strategy;

  // !!! THIS SINGLETON IMPLEMENTATION IS NOT THREAD-SAFE !!!
  static Single& instance() {
    static Single _instance ;
    return _instance ;
  }

 private: // -----------------------------------------------------------------

  // <never change their visibility> !!!
  Single() ;
  ~Single() ;
  Single(const Single&) ;
  Single& operator=(const Single&) ;
  // </never change their visibility> !!!
} ;

// some functions we put here so that everything can use them -----------------

Strategy cstring_to_strategy(const char* s);
size_t rand_between(size_t lower_bound, size_t upper_bound);
clock_t get_user_plus_system_times();
void remove(vector<int>& v, int i);

// return last element then pop it out
template <class T>
T pop_back(vector<T>& v) {

  T result = v.back();
  v.pop_back();

  return result;
}

// This is used for histogram output
// will break histogram output format if changed
template <class T>
ostream& operator<<(ostream& os, const vector<T>& v) {
  for (size_t i = 0 ; i < v.size() ; ++i) {
    os << v[i] << '\n';
  }
  return os << endl;
}

template <class T>
ostream& operator<<(ostream& os, const deque<T>& v) {
  for (size_t i = 0 ; i < v.size() ; ++i) {
    os << v[i] << '\n';
  }
  return os << endl;
}

// This is used for histogram output
// will break histogram output format if changed
template <class T, class U>
ostream& operator<<(ostream& os, const pair<T, U>& p) {
  return os << p.first << ' ' << p.second;
}

template <class T>
void sort(vector<T>& v) {
  sort(v.begin(), v.end());
}

// put in result intersection of fst and snd
// !!! BOTH FST AND SND NEED TO BE SORTED !!!
template <class T>
void intersect(vector<T>& fst,
               vector<T>& snd,
               vector<T>& result) {

  size_t fst_size = fst.size();
  size_t snd_size = snd.size();

  if (fst_size <= snd_size) { // iterate over fst
    for (size_t i = 0 ; i < fst_size ; ++i) {
      if (binary_search(snd.begin(), snd.end(), fst[i])) {
        result.push_back(fst[i]);
      }
    }
  } else { // iterate over snd
    for (size_t i = 0 ; i < snd_size ; ++i) {
      if (binary_search(fst.begin(), fst.end(), snd[i])) {
        result.push_back(snd[i]);
      }
    }
  }
}

float median(vector<float>& v);

float average(vector<float>& v);

#endif /* SINGLETON_H */
