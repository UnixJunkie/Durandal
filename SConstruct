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

import os, sys, commands

debug_mode   = True
profile_mode = False

for key, value in ARGLIST:
    if key == 'help':
        print "# production compile:"
        print "scons -j5"
        sys.exit(0)
    if key == 'debug':
        debug_mode = True
    if key == 'opt':
        debug_mode = False
    if key == 'prof':
        profile_mode = True

parse_flags_content = (" -D_USE_FAST_RMSD_ -D_SHOW_PERCENTAGE_COMPLETE_ " +
                       "-D_LARGE_DECOY_SET_ ")

link_flags = " "
common_flags = " -W -Wall "
if debug_mode:
    parse_flags_content = " -g -O0 " + common_flags
if not debug_mode:
    parse_flags_content = " -g -O3 -DNDEBUG " + common_flags
if profile_mode:
    parse_flags_content = " -pg -g -O3 -DNDEBUG " + common_flags
    link_flags = " -pg "

# dont' pollute current dir with object files
VariantDir('build', 'src', duplicate=0)

env = Environment (ENV = {'PATH' : os.environ['PATH'], # used by colorgcc
                          'TERM' : os.environ['TERM'],
                          'HOME' : os.environ['HOME']},
                   CXX = "g++",
                   CCFLAGS = parse_flags_content,
                   LINKFLAGS = link_flags)

cluster = env.Program(
    'durandal.cluster_pdbs',
    ['build/durandal.cc','build/DistMatrix.cc','build/DistRange.cc',
     'build/Stru.cc','build/SimpPDB.cc','build/rmsd.cc','build/Triple.cc',
     'build/Singleton.cc'])

rank_pdbs = env.Program('durandal.rank_pdbs',
                        ['build/ranker.cc','build/Stru.cc',
                         'build/SimpPDB.cc','build/rmsd.cc'])
