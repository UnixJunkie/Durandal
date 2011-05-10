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

#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "DistMatrix.h"
#include "Singleton.h"

using namespace std;

bool contains(int argc, char** argv, string param) {

  for (int i = 0 ; i < argc ; ++i) {

    if (param == string(argv[i])) {
      return true;
    }
  }
  return false;
}

bool contains(int argc, char** argv, string short_opt, string long_opt) {
  return (contains(argc, argv, short_opt) or
          contains(argc, argv, long_opt));
}

string get_option(int argc, char** argv, string option_flag) {

  for (int i = 0 ; i < argc ; ++i) {

    if (option_flag == string(argv[i])) {
      return string(argv[i + 1]);
    }
  }
  return string("");
}

string get_option(int argc, char** argv, string short_opt, string long_opt) {

  string res = get_option(argc, argv, short_opt);

  if (res == "") {
    res = get_option(argc, argv, long_opt);
  }

  return res;
}

int usage(char* program_name) {
  cout << "usage: " << program_name << " -i pdb_list -o file\n"
       << "                  {-b|-s} [-v]\n"
       << "STANDARD OPTIONS:\n"
       << "  --semi-auto dist_percent : choose cutoff to accept\n"
       << "                             approximately only\n"
       << "                             dist_percent smallest distances\n"
       << "                             ([0.03 .. 0.10] is reasonable)\n"
       << "                             incompatible with -d\n"
       << "  -d clustering_distance : user-chosen cutoff\n"
       << "                           incompatible with --semi-auto\n"
       << "  -b|--brute   : brute mode\n"
       << "  -m max_clusters : number of clusters to output\n"
       << "                    (default is 3, use -1 for all)\n"
       << "  -o|--output output_file : save clusters to file\n"
       << "  -s|--smart   : smart mode\n"
       << "  --stable     : output all pole position clusters\n"
       << "  -v|--verbose : verbose clusters listing\n"
       << "                 (default is index in input file only)\n"
       << "  -h|--help\n"
       << "RARELY USED OPTIONS:\n"
       << "  --compute-D  : output D of found clusters\n"
       << "  --no-lower   : don't use lower bound of CARMSD\n"
       << "  --no-upper   : don't use upper bound of CARMSD\n"
       << "  --rank       : output RMSD of each to first PDB of the list\n"
       << "################################################################"
       << "####\n"
       << "IF YOU USE THIS SOFTWARE, PLEASE CITE THE CORRESPONDING "
       << "PUBLICATION:\n"
       << "author  = {Berenger, Francois and Zhou, Yong and Shrestha, Rojan\n"
       << "           and Zhang, Kam Y. J.}\n"
       << "title   = {Entropy-accelerated exact clustering of protein "
       << "decoys}\n"
       << "journal = {Bioinformatics}\n"
       << "doi     = {10.1093/bioinformatics/btr072}\n"
       << "################################################################"
       << "####\n";
  exit(1);
}

// output cluster, remove its members from remaining_pdbs, unless this is
// a pole_position cluster
void print_cluster(ostream& out,
                   vector<int>& cluster,
                   vector<int>& remaining_pdbs,
                   DistMatrix& dm,
                   bool verbose,
                   bool pole_position,
                   bool compute_D) {

  int elt = cluster[0];
  bool has_energy = dm.has_energies();

  if (has_energy) {
    dm.display_energy(out, cluster);
  }
  if (compute_D) {
    out << "D:    " << dm.compute_D(cluster) << endl;
  }
  out << "members(" << cluster.size() << "):";
  if (verbose) {
    out << '\n' << dm.resolve(elt);
    if (has_energy) {
      for (size_t energy_index = 0 ; energy_index < dm.get_nb_energies();
           ++energy_index) {
        out << ':' << dm.resolve_energy(elt, energy_index);
      }
    }
  } else {
    out << ' ' << elt;
  }
  if (not pole_position) {
    remove(remaining_pdbs, elt);
  }
  for (size_t i = 1 ; i < cluster.size() ; ++i) {
    elt = cluster[i];
    if (verbose) {
      out << '\n' << dm.resolve(elt);
      if (has_energy) {
        for (size_t energy_index = 0 ; energy_index < dm.get_nb_energies();
             ++energy_index) {
          out << ':' << dm.resolve_energy(elt, energy_index);
        }
      }
    } else {
      out << ' ' << elt;
    }
    if (not pole_position) {
      remove(remaining_pdbs, elt);
    }
  }
  out << '\n';
}

void output_cluster(ostream& out,
                    vector<int>& cluster,
                    vector<int>& remaining_pdbs,
                    vector< vector<int> >& pole_position_clusters,
                    DistMatrix& dm,
                    bool verbose ,
                    bool stable,
                    bool compute_D) {

  if (pole_position_clusters.size() > 0) {
    if (stable) {
      out << "<pole position clusters("
          << pole_position_clusters.size() << ")>:\n";
      for (size_t i = 0 ; i < pole_position_clusters.size() ; ++i) {
        print_cluster(out, pole_position_clusters[i], remaining_pdbs,
                      dm, verbose, true, compute_D);
      }
      out << "</pole position clusters>\n";
    } else {
      out << "pole position centers("
          << pole_position_clusters.size() << "):";
      if (verbose) {
        out << '\n';
        for (size_t i = 0 ; i < pole_position_clusters.size() ; ++i) {
          out << dm.resolve(pole_position_clusters[i][0]) << '\n';
        }
      } else {
        for (size_t i = 0 ; i < pole_position_clusters.size() ; ++i) {
          out << ' ' << pole_position_clusters[i][0];
        }
        out << '\n';
      }
    }
  }
  print_cluster(out, cluster, remaining_pdbs, dm, verbose, false, compute_D);
}

int main(int argc, char** argv) {

  if (argc == 1 or contains(argc, argv, "-h", "--help")) {
    usage(argv[0]);
  }

  // seed Random Number Generator
  srand(time(NULL)); // comment   here for repeatability
  // srand(0);       // uncomment here for repeatability

  Single& single = Single::instance();

  bool stable             = contains(argc, argv, "--stable");
  bool auto_threshold     = contains(argc, argv, "--semi-auto");
  single._compute_D       = contains(argc, argv, "--compute-D");
  single._only_rank       = contains(argc, argv, "--rank");
  single._smart_mode      = contains(argc, argv, "-s", "--smart");
  single._verbose         = contains(argc, argv, "-v", "--verbose");
  single._brute_mode      = contains(argc, argv, "-b", "--brute");
  single._auto_switch     = true;
  single._no_lower_bound  = contains(argc, argv, "--no-lower");
  single._no_upper_bound  = contains(argc, argv, "--no-upper");
  single._output_to_file  = contains(argc, argv, "-o", "--output");
  single._use_URMSD       = false;

  string input_file = get_option(argc, argv, "-i");
  if (input_file == "") {
    cout << "error: -i pdb_list is mandatory" << endl;
    usage(argv[0]);
  }
  string distance = get_option(argc, argv, "-d");
  if (auto_threshold) {
    distance = string("10"); // will be overwritten later
  }
  string strategy        = get_option(argc, argv, "--strat");
  if (strategy != "") {
    single._reference_choosing_strategy = \
      cstring_to_strategy(strategy.c_str());
  } else {
    single._reference_choosing_strategy = \
      cstring_to_strategy("seq");
  }
  string output_file = get_option(argc, argv, "-o", "--output");
  int max_output = numeric_limits<int>::max();
  if (contains(argc, argv, "-m")) {
    int new_max_output = atoi(get_option(argc, argv, "-m").c_str());
    if (new_max_output != -1) {
      max_output = new_max_output; // limit set by user
    }
  } else {
    max_output = 3; // default output limit
  }
  if (not single._brute_mode and not single._smart_mode) {
    cout << "error: -b or -s is mandatory" << endl;
    return usage(argv[0]);
  }
  if (single._brute_mode and single._smart_mode) {
    cout << "error: -b XOR -s" << endl;
    return usage(argv[0]);
  }
  float clustering_distance = atof(distance.c_str());

  DistMatrix dm(input_file.c_str(), clustering_distance);

  if (auto_threshold) {

    string cutoff_str = get_option(argc, argv, "--semi-auto");
    float cutoff      = atof(cutoff_str.c_str());
    float threshold, median;

    cout << "using semi automatic threshold finding (" << cutoff
         << ')' << endl;
    threshold = dm.find_clustering_distance(cutoff, 0.0, &median);
    cout << "threshold found: " << threshold << endl;
    dm.set_clustering_distance(threshold);
  }

  if (single._only_rank) {
    dm.only_rank();
    return 0;
  }
  // will output the biggest clusters, in the order they are found
  // if two clusters have same size, only the first one is taken into account
  if (single._brute_mode) {
    dm.brute_init();
  } else { // assuming single._smart_mode
    dm.smart_init();
  }

  ofstream out (output_file.c_str());
  if (not out.is_open()) {
    cout << "error: cannot write to: " << output_file << endl;
    usage(argv[0]);
  }
  // find and print all clusters
  vector<int> remaining_pdbs = dm.get_all_PDBs_index();
  for (int printed = 0 ;
       remaining_pdbs.size() > 0 and printed < max_output ;
       ++printed) {

    vector<int> biggest_cluster;
    vector< vector<int> > pole_position_clusters;

    // find biggest
    dm.get_biggest_cluster(remaining_pdbs,
                           biggest_cluster,
                           pole_position_clusters);
    // output and remove from remaining
    output_cluster(out, biggest_cluster, remaining_pdbs,
                   pole_position_clusters, dm, single._verbose, stable,
                   single._compute_D);
  }
  out.close();

  return 0;
}
