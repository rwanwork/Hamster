///////////////////////////////////////////////////////////////////////////
//  HAMSTER
//  Software for depicting microarray data sets as a set of minimum spanning
//    trees.
//  
//  Version 1.3 -- August 26, 2011
//  
//  Copyright (C) 2009-2011 by Raymond Wan, All rights reserved.
//  Contact:  r-wan@cb.k.u-tokyo.ac.jp
//  Organization:  Department of Computational Biology, Graduate School of
//                 Frontier Science, University of Tokyo and
//                 Computational Biology Research Center, AIST,
//                 Japan
//  
//  This file is part of HAMSTER.
//  
//  HAMSTER is free software; you can redistribute it and/or 
//  modify it under the terms of the GNU Lesser General Public License 
//  as published by the Free Software Foundation; either version 
//  3 of the License, or (at your option) any later version.
//  
//  HAMSTER is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//  
//  You should have received a copy of the GNU Lesser General Public 
//  License along with HAMSTER; if not, see 
//  <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////////


/*******************************************************************/
/*!
    \file io.cpp
    Additional member functions for BUILDMST class definition
      Functions for reading input in
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: io.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>  //  cerr
#include <iomanip>  //  setw
#include <fstream>  //  ifstream
#include <queue>  // priority_queue

#include <boost/tokenizer.hpp>

using namespace std;
using namespace boost;

#include "global_defn.hpp"
#include "check.hpp"
#include "heapnode.hpp"
#include "vect_spear.hpp"
#include "vect.hpp"
#include "cluster.hpp"
#include "score.hpp"
#include "build_mst.hpp"


//!  Read the microarray data file in
/*!
     The data file must be tab-separated with an experiment on each line
     (assuming the user is building MSTs on the experiments).  The first
     row and column are headers and are basically ignored.  All other
     fields must be either floating point values or the string NULL.

     Each row in the data file translates into a VECT object.
*/
bool BUILDMST::readMicroarray () {
  unsigned int m = 0;
  unsigned int n = 0;
  string str;
  vector<string> tokens;
  VECT v;

  ifstream ma_fp (getMicroarrayFn ().c_str (), ios::in);
  if (!ma_fp) {
    cerr << "==\tError:  Input file " << getMicroarrayFn () << " could not be opened!" << endl;
    return false;
  }
  //  Get the header row; if eof already, then return false...something
  //    went wrong if we can't even get the first row
  getline (ma_fp, str);
  if (ma_fp.eof ()) {
    cerr << "Error:  Failed opening microarray data for reading!" << endl;
    return false;
  }

  while (true) {
    getline (ma_fp, str);
    if (ma_fp.eof ()) {
      break;
    }

    //  The integer m is the unique ID (starting from 0) of the vector
    v = VECT (m, str);
    if (n == 0) {
      n = v.getN ();
    }
    else if (n != v.getN ()) {
      cerr << "Error:  Mismatch in vector dimensions -- " << n << " vs " << v.getN () << endl;
      return false;
    }
    data.push_back (v);
    m++;
  }
  ma_fp.close ();

  //  Set the number of rows in the data set; same as number of nodes in the first MST
  setM (m);
  setN (n);

  if (getVerbose ()) {
    cerr << left << setw (VERBOSE_WIDTH) << "==\tMicroarray dimensions:" << getM () << " by " << getN () << endl;
  }

  return true;
}


//!  Read the optional attribute file in
/*!
     The attribute file is optional.  If it is unavailable, then every
     experiment is assumed to have the default attributes (DEFAULT_COLOUR
     and DEFAULT_SHAPE, found in buildmst.h).  The file is tab-separated with
     3 fields on each line:  (name, colour, shape).
*/
bool BUILDMST::readAttr () {
  unsigned int i;
  string str;
  vector<string> tokens;

  //  Test if no attribute file was provided; the file is not required,
  //  so return a true and not a false
  if (getAttrFn ().empty ()) {
    return true;
  }

  //  Open file for input
  ifstream attr_fp (getAttrFn ().c_str (), ios::in);
  if (!attr_fp) {
    cerr << "Attribute file could not be opened!" << endl;
    return false;
  }

  i = 0;
  while (true) {
    getline (attr_fp, str);
    if (attr_fp.eof ()) {
      break;
    }

    typedef tokenizer<char_separator<char> > tokenizer;
    char_separator<char> sep ("\t");
    tokenizer tokens (str, sep);

    vector<string> values;
    for (tokenizer::iterator beg = tokens.begin (); beg != tokens.end (); ++beg) {
      values.push_back (*beg);
    }

    //  Check if the line is valid; since attribute file is not required,
    //  errors here are just warnings
    if ((values.size () != 3) && (getVerbose ())) {
      cerr << "Warning:  Insufficient tokens read in on line " << i + 1 << "!" << endl;
    }
    else {
      //  Only set the information if the names match
      if (data[i].getName () == sanitizeSampleName (values[0])) {
        data[i].setName (values[0]);
        data[i].setColour (values[1]);
        data[i].setShape (values[2]);
      }
      else {
        cerr << "Error:  Microarray and attribute file names do not match:  " << data[i].getName () << " vs " << sanitizeSampleName (values[0]) << endl;
        cerr << "Most likely the order of the sample names in the two files are different." << endl;
      }
    }
    i++;
  }

  attr_fp.close ();

  if (getVerbose ()) {
    cerr << left << setw (VERBOSE_WIDTH) << "==\tExperiment attributes read in:" << i << endl;
  }

  return true;
}
