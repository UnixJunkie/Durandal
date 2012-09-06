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

The alignment expected by the Durandal software is 1 to 1.
For example, if 2 structures A and B are to be optimally superposed before
computing the CARMSD, both PDB files (A.pdb and B.pdb) must have same numbers
of alpha carbons.
The first alpha carbon in A.pdb will be matched to the first alpha carbon in
B.pdb.
The second alpha carbon in A.pdb will be matched to the second alpha carbon in
B.pdb.
And so on until the end of both files.
If it is not the case, please edit all your PDB files in advance.

IF YOU USE THIS SOFTWARE, PLEASE CITE THE CORRESPONDING PUBLICATION:

@article{Berenger2011,
author = {Berenger, Francois and Zhou, Yong and Shrestha, Rojan and Zhang, Kam Y. J.}, 
title = {Entropy-accelerated exact clustering of protein decoys}, 
volume = {27}, 
number = {7}, 
pages = {939-945}, 
year = {2011}, 
doi = {10.1093/bioinformatics/btr072}, 
URL = {http://bioinformatics.oxfordjournals.org/content/27/7/939.abstract}, 
eprint = {http://bioinformatics.oxfordjournals.org/content/27/7/939.full.pdf+html}, 
journal = {Bioinformatics},
}
