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

#include <vector>
#include <iostream>

#include "MiniCppUnit.h"
#include "DistMatrix.h"
#include "Singleton.h"

using namespace std;

class Tests : public TestFixture<Tests>
{

public:

  TEST_FIXTURE( Tests ) {
    TEST_CASE( test_DistMatrix_find_neighbors );
  }

  void test_DistMatrix_find_neighbors() {

    ASSERT_EQUALS( true, true );
  }

private:

};

REGISTER_FIXTURE( Tests );

