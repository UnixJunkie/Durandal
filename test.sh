#!/bin/bash

# Copyright (C) 2010, 2011 Zhang Initiative Research Unit,
# Advance Science Institute, Riken
# 2-1 Hirosawa, Wako, Saitama 351-0198, Japan

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#set -x

rm -f out current.out
./durandal.cluster_pdbs -i ./very_few_pdbs -d 2.0 --brute -o out -m -1 > \
/dev/null
egrep "^members|^pole" out > current.out
diff reference.out current.out
./durandal.cluster_pdbs -i ./very_few_pdbs -d 2.0 --smart -o out -m -1 > \
/dev/null
egrep "^members|^pole" out > current.out
diff reference.out current.out
# test clusters' energy characteristics
./durandal.cluster_pdbs -i ./very_few_pdbs_e -d 2.0 -s -o out -m -1 -v
diff out cluster_energies_reference.out
# test durandal.rank_pdbs
rm -f very_few_pdbs_rmsd_current
./durandal.rank_pdbs ./very_few_pdbs > very_few_pdbs_rmsd_current
diff very_few_pdbs_rmsd_current very_few_pdbs_rmsd_reference
