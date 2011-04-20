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
 *  **************************************************************************/

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <limits>

#include "DistRange.h"

const float DistRange::ABSOLUTE_MIN = 0.0;

float DistRange::ABSOLUTE_MAX() {
  return numeric_limits<float>::max();
}

// default, KEEP IN SYNC WITH was_initialized !!!
DistRange::DistRange() {
  _mini = ABSOLUTE_MIN;
  _maxi = ABSOLUTE_MAX();
  _is_decided = false;
}

// KEEP IN SYNC WITH DistRange::DistRange() !!!
bool DistRange::was_initialized() {
  return (_mini != ABSOLUTE_MIN) or (_maxi != ABSOLUTE_MAX());
}

// after exact rmsd measure
DistRange::DistRange(float exact_distance) {
  _mini = _maxi = exact_distance;
  _is_decided = true;
}

// range estimation
DistRange::DistRange(float mini, float maxi, float clustering_radius) {

  if (mini > maxi) {
    printf("error: illegal range with [%.2f,%.2f]\n", mini, maxi);
    exit(1);
  }
  _mini = mini;
  _maxi = maxi;
  _is_decided = (_mini > clustering_radius) || (_maxi <= clustering_radius);
}

void DistRange::check_against(float real_dist) {

  if (_mini > real_dist) {
    cout << "_mini > real: " << _mini << ' ' << real_dist << endl;
  }
  if (_maxi < real_dist) {
    cout << "_maxi < real: " << _maxi << ' ' << real_dist << endl;
  }
}

bool is_sharper(DistRange& r1, DistRange& r2) {
  return ((r1._maxi - r1._mini) < (r2._maxi - r2._mini));
}

// return (min(r1 - r2) ; max(r1 + r2))
DistRange minus_min_lim_plus_max_lim(DistRange& r1, DistRange& r2,
                                     float clustering_radius) {

  float new_min1 = max(0.0f, r1._mini - r2._maxi);
  float new_min2 = max(0.0f, r2._mini - r1._maxi);
  float new_max  = r1._maxi + r2._maxi;

  return DistRange(min(new_min1, new_min2), new_max, clustering_radius);
}

float DistRange::get_entropy() {
  float result = _maxi - _mini;
  return min(ENTROPY_BOUND, result);
}

bool DistRange::is_exact() {
  return _mini == _maxi;
}

ostream& operator<<(ostream& os, const DistRange& d) {
  return os << d._mini << ' ' << d._maxi;
}
