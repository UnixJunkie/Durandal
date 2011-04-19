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
 *  **************************************************************************
 *
 *  This class implements a triple to store 2 PDB indexes along with their
 *  CARMSD distance. Intentionally avoid templates because the type
 *  declarations in the code become lengthy and the compiler messages cryptic.
 */

#ifndef TRIPLE_H
#define TRIPLE_H

#include <ostream>
#include <string>

using namespace std;

class Triple {

 public:

  size_t x;
  size_t y;
  float xy;

  Triple(size_t x, size_t y, float xy);
};

bool operator<(const Triple& fst, const Triple& snd);

ostream& operator<<(ostream& os, const Triple& t);

#endif /* TRIPLE_H */
