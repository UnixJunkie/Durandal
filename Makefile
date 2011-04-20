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

all:    tags
	scons opt=true

debug:
	scons debug=true

tags:
	find . -name "*.cc" >  src_files
	find . -name "*.h"  >> src_files
	cat src_files | xargs etags -a -o src/TAGS
	rm -f src_files

test: all
	./test.sh

clean:
	scons -c
	rm -rf current.out src/TAGS build very_few_pdbs_rmsd_current \
curr_avg_RMSD .sconsign.dblite out

COPTS=-c -g -O3 -DNDEBUG -W -Wall
compile:
	# in case you don't have scons installed
	mkdir -p build
	g++ -o build/DistMatrix.o ${COPTS} src/DistMatrix.cc
	g++ -o build/DistRange.o ${COPTS} src/DistRange.cc
	g++ -o build/SimpPDB.o ${COPTS} src/SimpPDB.cc
	g++ -o build/Singleton.o ${COPTS} src/Singleton.cc
	g++ -o build/Stru.o ${COPTS} src/Stru.cc
	g++ -o build/durandal.o ${COPTS} src/durandal.cc
	g++ -o build/ranker.o ${COPTS} src/ranker.cc
	g++ -o build/rmsd.o ${COPTS} src/rmsd.cc
	g++ -o durandal.rank_pdbs build/ranker.o build/Stru.o build/SimpPDB.o \
build/rmsd.o
	g++ -o durandal.cluster_pdbs build/durandal.o build/DistMatrix.o \
build/DistRange.o build/Stru.o build/SimpPDB.o build/rmsd.o build/Singleton.o
