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
    \file process_scores.cpp
    Additional member functions for LAYOUTMST class definition
      Functions for reading and processing the scores
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: process_scores.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>  //  cerr
#include <fstream>  //  ifstream
#include <string>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

#include "global_defn.hpp"
#include "packet.hpp"
#include "score.hpp"
#include "vertex.hpp"
#include "edge.hpp"
#include "layout_mst.hpp"

//!  Comparison function for sorting scores into decreasing order by combined score
bool scoresGreaterCmp (SCORE x, SCORE y) {
  return (x.getCombinedScore () > y.getCombinedScore ());
}

//!  Comparison function for sorting scores into increasing order by ID
bool scoresGreaterIDCmp (SCORE x, SCORE y) {
  return (x.getID () < y.getID ());
}


//!  Read the data file of scores in
/*!
     The data file must be tab-separated with an experiment on each line
     (assuming the user is building MSTs on the experiments).  The first
     row and column are headers and are basically ignored.  All other
     fields must be either floating point values or the string NULL.

     Each row in the data file translates into a VECT object.
*/
bool LAYOUTMST::readScores () {
  string str;
  vector<string> out_tokens;

  typedef tokenizer<char_separator<char> > tokenizer;
  char_separator<char> sep (" \t");

  ifstream scores_fp (getScoresFn ().c_str (), ios::in);
  if (!scores_fp) {
    cerr << "==\tError:  Input file " << getScoresFn () << " could not be opened!" << endl;
    return false;
  }

  //  Get the number of iterations; if eof already, then return false...something
  //    went wrong if we can't even get the first row
  getline (scores_fp, str);
  if (scores_fp.eof ()) {
    cerr << "Error:  Failed opening scores file for reading!" << endl;
    return false;
  }
  tokenizer tokens (str, sep);

  setTotalIter (lexical_cast<unsigned int> (*(tokens.begin ())));

  while (true) {
    getline (scores_fp, str);
    if (scores_fp.eof ()) {
      break;
    }

    tokenizer tokens2 (str, sep);

    out_tokens.clear ();
    for (tokenizer::iterator beg = tokens2.begin (); beg != tokens2.end (); ++beg) {
      out_tokens.push_back (*beg);
    }

    if (out_tokens.size () == SCORES_FIELDS) {
       unsigned int iter = lexical_cast<unsigned int>(out_tokens[0]);
       unsigned int left = lexical_cast<unsigned int>(out_tokens[1]);
       unsigned int right = lexical_cast<unsigned int>(out_tokens[2]);
       double score1 = lexical_cast<double>(out_tokens[3]);
       double score2 = lexical_cast<double>(out_tokens[4]);
       double combined = lexical_cast<double>(out_tokens[5]);
       scores.push_back (SCORE (iter, left, right, score1, score2, combined));
    }
  }

  //  Close file pointer
  scores_fp.close ();

  return true;
}


//!  Process the scores by sorting and inserting into separate lists
/*!
     The jobs (MST IDs) are sorted by increasing score and then
     assigned to each process in a round-robin fashion.  This
     ensures that if there are k processors available, the k most
     important MSTs are processed first simultaneously by all
     processors.  Then, the next k most important MSTs are done, etc.
*/
void LAYOUTMST::processScores () {
  //  Proportion of images to generate
  double per = static_cast<double>(getPercent ()) / 100;

  //  Sort the scores by decreasing difference in scores or increasing ID, depending if getFixedPos () is true
  if (getFixedPos ()) {
    sort (scores.begin (), scores.end (), scoresGreaterIDCmp);
  }
  else {
    sort (scores.begin (), scores.end (), scoresGreaterCmp);
  }

  //  Initialize each row of the array
  for (unsigned int i = 0; i < getWorldSize (); i++) {
    all_workunits.push_back (vector <unsigned int>());
  }

  //  Create the lists
  for (unsigned int i = 0; i < getTotalIter (); i++) {
    //  Check for premature exit.  Comparing getPercent () to 100 is
    //  fine since they are both integers.  Generate at *least* "per"
    //  images, and then exit.
    if ((getPercent () != 100) && ((static_cast<double>(i) / static_cast<double>(getTotalIter ())) > per)) {
      break;
    }

    //  The list to add the value to
    unsigned int list = i % getWorldSize ();
    if (getDebug ()) {
      cerr << __FILE__ << "  [" << i % getWorldSize () << "]\t" << scores[i].getID () << endl;
    }
    all_workunits[list].push_back (scores[i].getID ());
  }

  return;
}


