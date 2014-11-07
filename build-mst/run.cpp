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
    \file run.cpp
    Member functions for BUILDMST class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: run.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>  //  cerr
#include <string>
#include <vector>
#include <queue>  //  priority_queue

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

using namespace std;

#include "global_defn.hpp"
#include "heapnode.hpp"
#include "vect_spear.hpp"
#include "vect.hpp"
#include "cluster.hpp"
#include "graph.hpp"
#include "score.hpp"
#include "build_mst.hpp"

//!  Execute the program after all parameters check out -- does the main work of the program
void BUILDMST::run () {
  //  Read the data in; return if either read produces an error
  //  The attribute file is *optional*, so not having one is not
  //    an error
  if (!readMicroarray () || !readAttr ()) {
    return;
  }

  //  Each experiment is a cluster, so we can initialize it now
  initializeClusters ();

  //  Calculate distance between experiments = clusters
  //    (at this initial stage)
  initializeDistances ();

  SCORE score = SCORE ();
  calculateScores (score);
  scores.push_back (score);

  unsigned int iter = 0;
  unsigned int M = getM ();
  for (iter = 0; iter < M; iter++) {
    //  Build a new graph, output the MST, and then destroy it
    GRAPH *g = new GRAPH (M - iter, pqueue, clusters);
    g -> printEdges (iter, clusters, getPath ());
    g -> printNodes (iter, clusters, getPath ());
    delete g;

    //  Find the next merge
    bool success = false;
    do {
      //  Dequeue
      if (pqueue.empty ()) {
        break;
      }
      HEAPNODE heapnode = pqueue.top ();
      pqueue.pop ();
      //  If not yet merged
      unsigned int left = heapnode.getLeft ();
      unsigned int right = heapnode.getRight ();
      if ((!clusters.at (left).haveAncestors ()) && (!clusters.at (right).haveAncestors ())) {
        success = true;

        if (getDebug ()) {
          cerr << "===\tAccepted:  [" << iter + 1 << "]  Size:  " << clusters.size () << "\tHeap nodes:  " << pqueue.size () << endl;

          cerr << "===\t\t(" << left << ", " << right << ")" << endl;
        }
        //  Create a new microarray vector in position clusters.size ().  The name of this merged node
        //  is clusters.size () - (number of rows in microarray).  i.e., from 0.
        CLUSTER c = CLUSTER (clusters.size (), &clusters[left], &clusters[right], getLinkage (), getM ());
        clusters.push_back (c);

        //  Calculate the similarity between this new node and every other node;
        //  add to the heap by pushing new edges on
        calculateLinkage (c);

        //  Add the score in
        score = SCORE (iter + 1, left, right);
        calculateScores (score);
        scores.push_back (score);
      }
      else {
        if (getDebug ()) {
          cerr << "===\tRejected:  (" << left << ", " << right << ")" << endl;
        }
      }
    } while (!success);
  }

  //  Normalize the scores to 0..100
  normalizeScores ();

  printScores (getPath ());

  return;
}
