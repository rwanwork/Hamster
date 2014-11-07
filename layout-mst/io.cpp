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
    Member functions for LAYOUTMST class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: io.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>  //  cerr
#include <fstream>  //  ifstream
#include <iomanip>  //  setprecision
#include <string>
#include <vector>

#include <cfloat>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

using namespace std;
using namespace boost;

#include "global_defn.hpp"
#include "packet.hpp"
#include "score.hpp"
#include "vertex.hpp"
#include "edge.hpp"
#include "layout_mst.hpp"


//!  Read in the nodes from file
/*!
     \param id ID of the MST (0-based)
*/
void LAYOUTMST::readNodes (unsigned int id) {
  unsigned int i = 0;
  string str;
  string fn;

  //////////////////////////////////////////////////
  //  Open the node attribute file for input
  fn = getPath () + lexical_cast<std::string>(id) + NODES_FILE_EXTENSION;

  ifstream fp (fn.c_str (), ios::in);
  if (!fp) {
    cerr << "Unexpected error:  File " + fn + " could not be opened!" << endl;
    exit (EXIT_FAILURE);
    return;
  }

  i = 0;
  while (true) {
    getline (fp, str);
    if (fp.eof ()) {
      break;
    }

    typedef tokenizer<char_separator<char> > tokenizer;
    char_separator<char> sep (" \t");
    tokenizer tokens (str, sep);

    vector<string> values;
    for (tokenizer::iterator beg = tokens.begin (); beg != tokens.end (); ++beg) {
      values.push_back (*beg);
    }

    //  Check if the line is valid; since attribute file is not required,
    //  errors here are just warnings
    if ((values.size () != 4) && (getVerbose ())) {
      cerr << "Warning [" << fn << "]:  Insufficient tokens read in on line " << i + 1 << "!" << endl;
    }
    else {
      vertices.push_back (VERTEX (values[0], values[1], values[2], lexical_cast<unsigned int>(values[3])));
    }
    i++;
  }
  fp.close ();

  return;
}


//!  Read in the edges from file
/*!
     \param id ID of the MST (0-based)
*/
double LAYOUTMST::readEdges (unsigned int id) {
  unsigned int i = 0;
  string str;
  string fn;
  double max_weight = 0.0;

  //////////////////////////////////////////////////
  //  Open the edge attribute file for input
  fn = getPath () + lexical_cast<std::string>(id) + EDGES_FILE_EXTENSION;

  ifstream fp2 (fn.c_str (), ios::in);
  if (!fp2) {
    cerr << "Unexpected error:  File " + fn + " could not be opened!" << endl;
    exit (EXIT_FAILURE);
  }

  i = 0;
  while (true) {
    getline (fp2, str);
    if (fp2.eof ()) {
      break;
    }

    typedef tokenizer<char_separator<char> > tokenizer;
    char_separator<char> sep (" \t");
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
      double weight = lexical_cast<double>(values[2]);
      edges.push_back (EDGE (values[0], values[1], weight));
      //  If weight is DBL_MAX, then it should be greater than (DBL_MAX - DBL_EPSILON).
      //  Test avoids testing of floating point.  The max_weight is the largest
      //  weight smaller than DBL_MAX.
      if ((weight > max_weight) || (weight >= (DBL_MAX - DBL_EPSILON))) {
        max_weight = weight;
      }
    }
    i++;
  }

  fp2.close ();

  return (max_weight);
}


//!  Update the node positions
/*!
     \param id ID of the MST (0-based)

     Read in the previous MST and add positions to each node.  If the current MST is 0, then immediately return.
*/
void LAYOUTMST::updateNodePositions (unsigned int id) {
  const regex expr ("\\s*(\\S+)\\s*\\[label.+pos=\"(\\d+),(\\d+)\", width=\"([\\d\\.]+)\", height=\"([\\d\\.]+)\".+");
  smatch what;

  unsigned int i = 0;
  string str;
  string fn;

  if (id == 0) {
    return;
  }

  //////////////////////////////////////////////////
  //  Open the edge attribute file for input
  fn = getPath () + lexical_cast<std::string>(id - 1) + GV_FILE_EXTENSION;

  ifstream fp2 (fn.c_str (), ios::in);
  if (!fp2) {
    cerr << "Unexpected error:  File " + fn + " could not be opened!" << endl;
    exit (EXIT_FAILURE);
  }

  i = 0;
  while (true) {
    getline (fp2, str);
    if (fp2.eof ()) {
      break;
    }

    if (regex_match (str, what, expr)) {
      string key = what[1];
      unsigned int x = lexical_cast<unsigned int>(what[2]);
      unsigned int y = lexical_cast<unsigned int>(what[3]);
      float width = lexical_cast<float>(what[4]);
      float height = lexical_cast<float>(what[5]);
      for (vector<VERTEX>::iterator iter = vertices.begin (); iter != vertices.end (); iter++) {
        if (iter -> getName () == key) {
          iter -> setX (x);
          iter -> setY (y);
          iter -> setWidth (width);
          iter -> setHeight (height);
          iter -> setUpdated (true);
          break;
        }
      }
    }
    i++;
  }

  fp2.close ();

  return;
}

