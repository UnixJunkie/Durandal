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

#ifndef DIST_RANGE_H
#define DIST_RANGE_H

using namespace std;

static const float ENTROPY_BOUND = 100.0;

class DistRange {

 private:

  static const float ABSOLUTE_MIN = 0.0;
  static float ABSOLUTE_MAX(); // FBR: dirty but works

 public:

  // most atributes are public for easier access
  float _mini;
  float _maxi;
  bool  _is_decided;

  DistRange();
  DistRange(float exact_distance);
  DistRange(float mini, float maxi, float clustering_radius);
  virtual ~DistRange() {};
  void check_against(float real_dist);
  float get_entropy();
  bool is_exact();
  bool was_initialized();
};

ostream& operator<<(ostream& os, const DistRange& d);

bool is_sharper(DistRange& r1, DistRange& r2);

DistRange minus_min_lim_plus_max_lim(DistRange& r1, DistRange& r2,
                                     float clustering_radius);

#endif /* DIST_RANGE_H */
