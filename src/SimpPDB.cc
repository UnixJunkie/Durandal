/*  **************************************************************************
 *  Copyright 2009 Shuai Cheng Li and Yen Kaow Ng
 *  **************************************************************************
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
 *  BMC Bioinformatics. 2010; 11: 25.
 *  Published online 2010 January 13. doi: 10.1186/1471-2105-11-25.
 *  "Calibur: a tool for clustering large numbers of protein decoys"
 *  Shuai Cheng Li and Yen Kaow Ng
 *  **************************************************************************/

#include <cmath>
#include <iostream>
#include <fstream>

#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

using namespace std;

#include "SimpPDB.h"

const char* aa[]={"BCK","GLY","ALA","SER","CYS","VAL","THR","ILE",
                  "PRO","MET","ASP","ASN","LEU",
                  "LYS","GLU","GLN","ARG",
                  "HIS","PHE","TYR","TRP","CYX", "MSE"};

const char slc[]={'X','G','A','S','C','V','T','I',
                  'P','M','D','N','L','K','E','Q','R',
                  'H','F','Y','W','C', 'm'};

int toInt(const string& aString)
{
  static char st[20];
  int start=0;
  for(size_t i=0; i<aString.size(); i++)
    {
      if(aString[i]!=' ')
        {
          st[start++]=aString[i];
        }
    }
  st[start]='\0';
  int rev=atoi(st);
  return rev;
}

float toFloat(const string& aString)
{
  static char st[20];
  int start=0;
  for(size_t i=0; i<aString.size(); i++)
    {
      if(aString[i]!=' ')
        {
          st[start++]=aString[i];
        }
    }
  st[start]='\0';
  float rev=atof(st);
  return rev;
}


SimPDB::SimPDB(const char *aFileName, bool use_URMSD)
{
  mProteinFileName=aFileName;
  mNumResidue=LONGEST_CHAIN;
  mCAlpha=new float[3*LONGEST_CHAIN];
  read(-1, use_URMSD);
}

SimPDB::SimPDB(const char *aFileName, int len, bool use_URMSD)
{
  mProteinFileName=aFileName;
  mNumResidue=len;
  mCAlpha=new float[3*len];
  read(len, use_URMSD);
}

SimPDB::~SimPDB()
{
  delete [] mCAlpha;
}

void
SimPDB::read(int expected_count, bool use_URMSD)
{
  ifstream input(mProteinFileName);
  if(!input) {
    cerr<<"Cannot find protein file |" << mProteinFileName << "|" << endl;
    exit(0);
  }
  cout.flush();
  char buf[400];
  mNumResidue=0;
  float x, y, z;
  bool read=false;

  int prevID=-10000;
  int count=0;
  int count3=0;
  while (!input.eof()) {

    input.getline(buf, 400);
    string line = buf;

    if (line.substr(0, 3) == "TER" && read == true)
      break;
    if (line.substr(0, 6) == "ENDMDL")
      break;
    if (line.substr(0, 4) != "ATOM" && line.substr(0, 6) != "HETATM")
      continue;
    if (line.substr(13, 4) == "CA  " ||
        line.substr(13, 4) == " CA " ||
        line.substr(13, 4) == "  CA" ||
        line.substr(13, 2) == "CA") {

      if (toupper(line[21]) == 'A' ||
          toupper(line[21]) == 'C' ||
          line[21]          == ' ') {

        read = true;
        int residueID = toInt(line.substr(22, 6));
        if (residueID == prevID)
          continue;
        prevID = residueID;
        x = toFloat(line.substr(30, 8));
        y = toFloat(line.substr(38, 8));
        z = toFloat(line.substr(46, 8));
        string AAType = line.substr();
        count3 = 3 * count;
        if (expected_count == -1) {
          if (count < LONGEST_CHAIN) { // don't write after mCAlpha
            mCAlpha[count3]     = x;
            mCAlpha[count3 + 1] = y;
            mCAlpha[count3 + 2] = z;
//             cerr << "CA: "
//                  << mCAlpha[count3]     << ' '
//                  << mCAlpha[count3 + 1] << ' '
//                  << mCAlpha[count3 + 2] << endl;
          }
        } else {
          if (count < expected_count) {
            // don't store more than we expect once we know the
            // number of CAs of the first PDB
            mCAlpha[count3]     = x;
            mCAlpha[count3 + 1] = y;
            mCAlpha[count3 + 2] = z;
//             cerr << "CA: "
//                  << mCAlpha[count3]     << ' '
//                  << mCAlpha[count3 + 1] << ' '
//                  << mCAlpha[count3 + 2] << endl;
          }
        }
        ++count;
      }
    }
  }
  //cerr << "END" << endl;
  mNumResidue=count;
  input.close();
  if (expected_count != -1 and count != expected_count) {
    cout << "fatal: " << mProteinFileName << ": read " << count
         << " atoms while expecting " << expected_count << endl;
    exit(1);
  }
  if (use_URMSD) {

    for (int i = 0 ; i < mNumResidue - 1 ; ++i) {

      int i3 = 3 * i;
      // Ca_i
      float Xi = mCAlpha[i3 + 0];
      float Yi = mCAlpha[i3 + 1];
      float Zi = mCAlpha[i3 + 2];
      //cerr << "i: " << Xi << ' ' << Yi << ' ' << Zi << endl;
      // Ca_j (j = i + 1)
      float Xj = mCAlpha[i3 + 3];
      float Yj = mCAlpha[i3 + 4];
      float Zj = mCAlpha[i3 + 5];
      //cerr << "j: " << Xj << ' ' << Yj << ' ' << Zj << endl;
      // vector(Ca_i, Ca_j)
      Xi = Xj - Xi;
      Yi = Yj - Yi;
      Zi = Zj - Zi;
      //cerr << "d: " << Xi << ' ' << Yi << ' ' << Zi << endl;
      // normalize it
      float len = sqrt(Xi*Xi + Yi*Yi + Zi*Zi);
      mCAlpha[i3 + 0] = Xi / len;
      mCAlpha[i3 + 1] = Yi / len;
      mCAlpha[i3 + 2] = Zi / len;
      //cerr << mCAlpha[i3 + 0] << ' '
      //     << mCAlpha[i3 + 1] << ' '
      //     << mCAlpha[i3 + 2] << endl;
    }
    //cerr << "END" << endl;
  } else {
    // molecule centering (center of mass is at (0, 0, 0) after)
    float cx = 0;
    float cy = 0;
    float cz = 0;

    int i3=0;
    for(int i=0; i<mNumResidue; i++)
      {
        cx+=mCAlpha[i3];
        cy+=mCAlpha[i3+1];
        cz+=mCAlpha[i3+2];
        i3+=3;
      }

    cx/=mNumResidue;
    cy/=mNumResidue;
    cz/=mNumResidue;
    i3=0;
    for(int i=0; i<mNumResidue; i++)
      {
        mCAlpha[i3]-=cx;
        mCAlpha[i3+1]-=cy;
        mCAlpha[i3+2]-=cz;
        i3+=3;
      }
  }
}
