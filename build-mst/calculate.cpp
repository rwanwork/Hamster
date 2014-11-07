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
    \file calculate.cpp
    Additional member functions for BUILDMST class definition
      Functions for calculations
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: calculate.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>  //  cerr, endl
#include <iomanip>  //  setw
#include <string>
#include <vector>
#include <queue>  //  priority_queue
#include <fstream>  //  ofstream

#include <cstdlib>  //  exit, EXIT_FAILURE

using namespace std;

#include "global_defn.hpp"
#include "heapnode.hpp"
#include "vect_spear.hpp"
#include "vect.hpp"
#include "cluster.hpp"
#include "score.hpp"
#include "build_mst.hpp"


//!  Initialize the distance matrix
/*!
     All of the functions that it calls must be distance (or dissimilarity)
     functions.  That is, a low value (0) indicates highly similar and a
     high value indicates highly dissimilar.

     Distances are stored twice in the matrix -- on both sides of the
     diagonal.  As distances are calculated, they are added into the
     priority queue.
*/
void BUILDMST::initializeDistances () {
  double score = 0.0;
  unsigned int total = 0;
  HEAPNODE heapnode;

  unsigned int i;
  unsigned int j;
  unsigned int end;

  //  Create the distance matrix
  unsigned int m = getM ();
  dist_matrix = new double*[m];
  for (i = 0; i < m; i++) {
    dist_matrix[i] = new double[m];
    for (j = 0; j < m; j++) {
      dist_matrix[i][j] = 0.0;
    }
  }

  end = data.size ();
  for (i = 0; i < end; i++) {
    //  No self-loops allowed in graph
    for (j = i + 1; j < end; j++) {
      //  All functions return a dissimilarity score
      switch (getDistance ()) {
        case DIST_EUC :
          score = data[i].simEuc (&data[j]);
          break;
        case DIST_MAN :
          score = data[i].simMan (&data[j]);
          break;
        case DIST_PEAR :
          score = data[i].simPear (&data[j]);
          break;
        case DIST_SPEAR :
          score = data[i].simSpear (&data[j]);
          break;
      }
      if (getDebug ()) {
        //  Print the distance we calculated
        cout << setprecision (6) << score << "\t" << i << "\t" << j << endl;
      }

      //  Add the distance to the heap; the priority queue is actually
      //  a queue of clusters and not experiment nodes.  But at this
      //  stage, they both mean the same thing; so it is fine to do this
      //  as long as the vector of ID i is also in cluster ID i.
      heapnode = HEAPNODE (data[i].getID (), data[j].getID (), score);
      pqueue.push (heapnode);

      //  Set the score on both sides of the diagonal
      dist_matrix[i][j] = score;
      dist_matrix[j][i] = score;

      total++;
    }
  }

  if (getVerbose ()) {
    cerr << left << setw (VERBOSE_WIDTH) << "==\tNumber of pairs calculated:" << total << endl;
  }

  return;
}

//!  Initialize the clusters with the data file
/*!
     Since bottom-up clustering starts off with each object in its own
     cluster, this function initializes the clusters vector.
*/
void BUILDMST::initializeClusters () {
  unsigned int M = getM ();
  unsigned int i = 0;
  CLUSTER c;

  for (i = 0; i < M; i++) {
    c = CLUSTER (i, data[i].getName (), data[i].getColour (), data[i].getShape (), data[i]);
    clusters.push_back (c);
  }

  return;
}


//!  Calculate the linkage
/*!
     For a given cluster, the linkage between it and every other cluster
     is calculated and added into the priority queue.
*/
void BUILDMST::calculateLinkage (CLUSTER &arg) {
  double score = 0.0;
  LINK_METHOD linkage = getLinkage ();
  HEAPNODE heapnode;

  unsigned int i = arg.getID ();
  unsigned int j = 0;
  unsigned int end = clusters.size ();
  for (j = 0; j < end; ++j) {
    //  Don't do a comparison if the same node
    if ((i == j) || (clusters[j].haveAncestors ())) {
      continue;
    }

    switch (linkage) {
      case LINK_SINGLE :
        score = clusters[i].linkSingle (&clusters[j], dist_matrix);
        break;
      case LINK_AVERAGE :
        score = clusters[i].linkAverage (&clusters[j], dist_matrix);
        break;
      case LINK_COMPLETE :
        score = clusters[i].linkComplete (&clusters[j], dist_matrix);
        break;
      case LINK_CENTROID :
        score = clusters[i].linkCentroid (&clusters[j], &data, getCentroid ());
        break;
    }
    heapnode = HEAPNODE (clusters[i].getID (), clusters[j].getID (), score);
    pqueue.push (heapnode);
  }

  return;
}


//!  Calculate the graph scores
/*!
     The score for the current graph configuration (based on
     intra and inter-cluster edge weights) is calculated.
*/
void BUILDMST::calculateScores (SCORE &arg) {

  switch (getScoreMethod ()) {
    case SCORE_GAPS     :
     arg.scoreGaps (getM (), &clusters, dist_matrix, getDebug ());
      break;
    case SCORE_ANOVA :
      arg.scoreANOVA (getM (), &clusters, dist_matrix, getDebug ());
      break;
    case SCORE_NASSOC :
      arg.scoreNormalizedAssoc (getM (), &clusters, dist_matrix, getDebug ());
      break;
    case SCORE_NASSOC_ORIG :
      arg.scoreNormalizedAssocOrig (getM (), &clusters, dist_matrix, getDebug ());
      break;
  }

  return;
}


//!  Normalize the scores in the scores vector
/*!
     Make two passes over the scores vector.  First, find the
     maximum combined score.  Then divide by this value and
     multiply by 100 so that they are in the range [0, 100].

     If no maximum was found, then everything has a value of 0.0.
     This indicates something might have gone wrong (i.e., data set
     too small, etc.), but we set everything to 100 instead.
*/
void BUILDMST::normalizeScores () {
  unsigned int i = 0;
  unsigned int end = 0;
  double max_combined = 0.0;
  bool found = false;

  end = scores.size ();
  for (i = 0; i < end; i++) {
    if (scores[i].getCombinedScore () > max_combined) {
      max_combined = scores[i].getCombinedScore ();
      found = true;
    }
  }

  //  Was a value larger than 0.0 found?
  if (found) {
    for (i = 0; i < end; i++) {
      scores[i].setCombinedScore (scores[i].getCombinedScore () / max_combined * 100);
    }
  }
  else {
    for (i = 0; i < end; i++) {
      scores[i].setCombinedScore (100);
    }
  }

  return;
}


//!  Print all of the scores out to a text file
bool BUILDMST::printScores (string outpath) {
  unsigned int i = 0;
  unsigned int end = 0;
  string fn = outpath + SCORES_FILENAME;

  ofstream fout (fn.c_str ());
  if (!fout) {
    cerr << "File " << fn << " could not be opened!" << endl;
    exit (EXIT_FAILURE);
  }

  fout << getM () << endl;
  end = scores.size ();
  for (i = 0; i < end; i++) {
    fout << scores[i].getID () << "\t"
        << scores[i].getLeft () << "\t"
        << scores[i].getRight () << "\t"
        << scores[i].getScore1 () << "\t"
        << scores[i].getScore2 () << "\t"
        << scores[i].getCombinedScore () << endl;
  }
  fout.close ();

  return true;
}
