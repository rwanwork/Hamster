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
    \file score.cpp
    Additional member functions for SCORE class definition
      Functions for scoring the clusters
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: score.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>
#include <string>
#include <vector>
#include <cmath>  //  fabs
#include <climits>  //  UINT_MAX
#include <cfloat>  //  DBL_MAX

using namespace std;

#include "global_defn.hpp"
#include "vect_spear.hpp"
#include "vect.hpp"
#include "cluster.hpp"
#include "score.hpp"

//!  Constructor that takes no arguments
SCORE::SCORE ()
  : id (0),
    left (UINT_MAX),
    right (UINT_MAX),
    score1 (0.0),
    score2 (0.0),
    combined (0.0)
{
}

//!  Constructor that takes three arguments
/*!
    \param arg1 ID of the node
    \param arg2 ID of the left cluster
    \param arg3 ID of the right cluster
*/
SCORE::SCORE (unsigned int arg1, unsigned int arg2, unsigned int arg3)
  : id (arg1),
    left (arg2),
    right (arg3),
    score1 (0.0),
    score2 (0.0),
    combined (0.0)
{
}

//!  Set the ID
void SCORE::setID (unsigned int arg) {
  id = arg;
}

//!  Get the ID
unsigned int SCORE::getID () const {
  return id;
}

//!  Set the ID of the left cluster
void SCORE::setLeft (unsigned int arg) {
  left = arg;
}

//!  Get the ID of the left cluster
unsigned int SCORE::getLeft () const {
  return left;
}

//!  Set the ID of the right cluster associated with this score
void SCORE::setRight (unsigned int arg) {
  right = arg;
}

//!  Get the ID of the right cluster
unsigned int SCORE::getRight () const {
  return right;
}

//!  Set score 1
void SCORE::setScore1 (double arg) {
  score1 = arg;
}

//!  Get score 1
double SCORE::getScore1 () const {
  return score1;
}

//!  Set score 2
void SCORE::setScore2 (double arg) {
  score2 = arg;
}

//!  Get score 2
double SCORE::getScore2 () const {
  return score2;
}

//!  Set combined score
void SCORE::setCombinedScore (double arg) {
  combined = arg;
}

//!  Get combined score
double SCORE::getCombinedScore () const {
  return combined;
}


//!  Calculate the MST score based on the size of the gap between the largest intra-cluster weight and the smallest inter-cluster weight
/*!
     The main data structure is a matrix called distance_to_cluster which maps a
     distance at d[i][j] to the cluster it is in; the value UINT_MAX is used
     to indicate it is not in a cluster (i.e., it is an inter-cluster relationship).

     After constructing this data structure, it is used to separate all distances
     into those that are within a cluster and those that are not.  As clustering
     is performed by increasing distances, we seek to find the size of the
     gap (in distance) between the largest intra-cluster distance and the smallest
     inter-cluster distance.  Once this is found, we set it and take the absolute
     value of their difference.
*/
void SCORE::scoreGaps (unsigned int M, vector<CLUSTER> *clusters, double **d, bool debug) {
  unsigned int i = 0;  //  Position i
  unsigned int j = 0;  //  Position j
  unsigned int max_items = 0;
  unsigned int expt_i = 0;  //  Experiment ID at position i
  unsigned int expt_j = 0;  //  Experiment ID at position j

  unsigned int k = 0;  //  Cluster ID

  bool **distance_to_cluster;  //  (i, j) TRUE means both experiments are in the same cluster

  double intra_max = 0.0;
  double inter_min = DBL_MAX;

  //  Initialize the distance_to_cluster
  distance_to_cluster = new bool*[M];
  for (i = 0; i < M; i++) {
    distance_to_cluster[i] = new bool[M];
    for (j = 0; j < M; j++) {
      distance_to_cluster[i][j] = false;
    }
  }

  //  For each cluster (with no ancestors)
  for (k = 0; k < clusters -> size (); k++) {
    if (!clusters -> at (k).haveAncestors ()) {
      vector<unsigned int> items = clusters -> at (k).getItems ();
      max_items = items.size ();

      // For every pair of items
      for (i = 0; i < max_items ; i++) {
        expt_i = items[i];
        for (j = i + 1; j < max_items; j++) {
          expt_j = items[j];

          //  Indicate that the pair are in the same cluster
          distance_to_cluster[expt_i][expt_j] = true;
          distance_to_cluster[expt_j][expt_i] = true;
        }
      }
    }
  }

  //  Use distance_to_cluster to scan the upper-triangular distance matrix
  unsigned int intra_count = 0;
  unsigned int inter_count = 0;
  for (i = 0; i < M; i++) {
    for (j = i + 1; j < M; j++) {
      if (distance_to_cluster[i][j]) {
        //  Intra-cluster score
        if (d[i][j] > intra_max) {
          intra_max = d[i][j];
          intra_count++;
        }
      }
      else {
        //  Inter-cluster score
        if (d[i][j] < inter_min) {
          inter_min = d[i][j];
          inter_count++;
        }
      }
    }
  }

  //  If either is 0 that means this is the first MST or the last one; set them to 0
  if ((intra_count == 0) || (inter_count == 0)) {
    setScore1 (0.0);
    setScore2 (0.0);
    setCombinedScore (0.0);
  }
  else {
    setScore1 (intra_max);
    setScore2 (inter_min);
    setCombinedScore (fabs (intra_max - inter_min));
  }

  if (debug) {
    cerr << "====>\t" << getScore1 () << "\t" << getScore2 () << "\t" << getCombinedScore () << endl;
  }

  return;
}


//!  Calculate the MST score based on the ANOVA of the two groups
/*!
     The distances between experiments are separated into two groups (intra-cluster
     and inter-cluster) and their sum of squares of both are calculated against
     the global mean.  That is, the mean is unchanged with each application of
     this function -- it is simply the mean of all the values in the distance
     matrix.

     The main data structure is a matrix called distance_to_cluster (see SCORE::scoreGaps).

     A vector of intra- and inter-cluster scores are kept since we make two passes
     over the data (the first to get the mean; the second to calculate the sums of
     squares).
*/
void SCORE::scoreANOVA (unsigned int M, vector<CLUSTER> *clusters, double **d, bool debug) {
  unsigned int i = 0;  //  Position i
  unsigned int j = 0;  //  Position j
  unsigned int max_items = 0;
  unsigned int expt_i = 0;  //  Experiment ID at position i
  unsigned int expt_j = 0;  //  Experiment ID at position j

  unsigned int k = 0;  //  Cluster ID

  bool **distance_to_cluster;  //  (i, j) TRUE means both experiments are in the same cluster
  vector<double> intra_scores;
  vector<double> inter_scores;

  double mean = 0.0;
  double intra_mean = 0.0;
  double inter_mean = 0.0;
  double intra_size = 0.0;
  double inter_size = 0.0;

  //  Initialize the distance_to_cluster
  distance_to_cluster = new bool*[M];
  for (i = 0; i < M; i++) {
    distance_to_cluster[i] = new bool[M];
    for (j = 0; j < M; j++) {
      distance_to_cluster[i][j] = false;
    }
  }

  //  For each cluster (with no ancestors)
  for (k = 0; k < clusters -> size (); k++) {
    if (!clusters -> at (k).haveAncestors ()) {
      vector<unsigned int> items = clusters -> at (k).getItems ();
      max_items = items.size ();

      // For every pair of items
      for (i = 0; i < max_items ; i++) {
        expt_i = items[i];
        for (j = i + 1; j < max_items; j++) {
          expt_j = items[j];

          //  Indicate that the pair are in the same cluster
          distance_to_cluster[expt_i][expt_j] = true;
          distance_to_cluster[expt_j][expt_i] = true;
        }
      }
    }
  }

  //  Use distance_to_cluster to scan the upper-triangular distance matrix
  for (i = 0; i < M; i++) {
    for (j = i + 1; j < M; j++) {
      mean += d[i][j];
      if (distance_to_cluster[i][j]) {
        //  Intra-cluster score
        intra_scores.push_back (d[i][j]);
        intra_mean += d[i][j];
        mean += d[i][j];
      }
      else {
        //  Inter-cluster score
        inter_scores.push_back (d[i][j]);
        inter_mean += d[i][j];
        mean += d[i][j];
      }
    }
  }

  //  If either is 0 that means this is the first MST or the last one; set them to 0 and leave
  if ((intra_scores.size () == 0) || (inter_scores.size () == 0)) {
    setScore1 (0.0);
    setScore2 (0.0);
    setCombinedScore (0.0);
    return;
  }

  intra_size = static_cast<double>(intra_scores.size ());
  inter_size = static_cast<double>(inter_scores.size ());
  intra_mean = intra_mean / intra_size;
  inter_mean = inter_mean / inter_size;
  mean = mean / (intra_size + inter_size);

  double ss_groups = intra_size * (intra_mean - mean) * (intra_mean - mean) +
                     inter_size * (inter_mean - mean) * (inter_mean - mean);
  double ms_groups = ss_groups;  //  k = 2; so k - 1 = 1

  double ss_total = 0.0;
  for (unsigned int i = 0; i < intra_scores.size (); i++) {
    ss_total += ((intra_scores[i] - mean) * (intra_scores[i] - mean));
  }
  for (unsigned int i = 0; i < inter_scores.size (); i++) {
    ss_total += ((inter_scores[i] - mean) * (inter_scores[i] - mean));
  }

  double sse = ss_total - ss_groups;
  double mse = sse / (intra_size + inter_size - 2);

  //  Free arrays
  for (i = 0; i < M; i++) {
    delete [] distance_to_cluster[i];
  }
  delete [] distance_to_cluster;

  setScore1 (ms_groups);
  setScore2 (mse);
  setCombinedScore (ms_groups / mse);

  if (debug) {
    cerr << "====>\t" << getScore1 () << "\t" << getScore2 () << "\t" << getCombinedScore () << endl;
  }

  return;
}

//!  Calculate the MST score based on the normalized association of Shi and Malik (2000)
/*!
     Defn:  Assoc (A, A) is the sum of all weights for all distances within cluster A;
     Assoc (A, V) is the sum of all weights for all distances from nodes in A with all
     other nodes.  If A is a cluster with only one node, then Assoc (A, A) = 0 since there
     are no self-loops.

     The main data structure is a matrix called distance_to_cluster (see SCORE::scoreGaps).

     This method performs the following steps:  (1)  Build distance_to_cluster;
     (2)  Calculate the self-association, Assoc (A, A), for each cluster;
     (3)  Calculate the all-association, Assoc (A, V), for each cluster;
     (4)  Accumulate the scores for all clusters that were "valid" [valid
     clusters are top-level clusters with no ancestors]; (5)  Adjust the
     scores by multiplying it by the number of clusters; (6)  Free memory
     and store values.

     Step (5) is not part of the definition of normalized association of Shi and
     Malik (2000).

*/
void SCORE::scoreNormalizedAssoc (unsigned int M, vector<CLUSTER> *clusters, double **d, bool debug) {
  unsigned int i = 0;  //  Position i
  unsigned int j = 0;  //  Position j
  unsigned int max_items = 0;
  unsigned int expt_i = 0;  //  Experiment ID at position i
  unsigned int expt_j = 0;  //  Experiment ID at position j

  unsigned int k = 0;  //  Cluster ID

  //  (i, j) is the cluster ID the edge is in
  //  Map distance (i, j) to its cluster
  unsigned int **distance_to_cluster;
  double *assoc_self = NULL;
  double *assoc_all = NULL;
  vector<bool> valid_cluster;

  unsigned int total_clusters = clusters -> size ();
  unsigned int num_clusters = 0;

  //  Initialize the distance_to_cluster
  distance_to_cluster = new unsigned int * [M];
  for (i = 0; i < M; i++) {
    distance_to_cluster[i] = new unsigned int[M];
    for (j = 0; j < M; j++) {
      //  Initialize with UINT_MAX to indicate inter-cluster edges
      distance_to_cluster[i][j] = UINT_MAX;
    }
  }

  //  For each cluster
  for (k = 0; k < total_clusters; k++) {
    if (!clusters -> at (k).haveAncestors ()) {
      //  Everything within this cluster
      vector<unsigned int> items = clusters -> at (k).getItems ();
      max_items = items.size ();

//        cerr << "Cluster has " << items.size () << " items." << endl;

      // For every pair of items
      for (i = 0; i < max_items ; i++) {
        expt_i = items[i];
        for (j = i; j < max_items; j++) {
          expt_j = items[j];

          //  Indicate that the pair are in the same cluster

//           cerr << "\t\tTRUE:  (" << expt_i << ", " << expt_j << ") --> " << k << endl;

          distance_to_cluster[expt_i][expt_j] = k;
          distance_to_cluster[expt_j][expt_i] = k;
        }
      }
      valid_cluster.push_back (true);
      num_clusters++;
    }
    else {
      valid_cluster.push_back (false);
    }
  }

  assoc_self = new double[total_clusters];
  assoc_all = new double[total_clusters];

  //  Initialize the arrays
  for (i = 0; i < total_clusters; i++) {
    assoc_self[i] = 0.0;
    assoc_all[i] = 0.0;
  }

  //  Scan only the upper-triangular distance_to_cluster matrix to find
  //  the intra-cluster edges  (i.e., Assoc (A, A))
  for (i = 0; i < M; i++) {
    for (j = i + 1; j < M; j++) {
      if (distance_to_cluster[i][j] != UINT_MAX) {
        //  Intra-cluster edge
//         cerr << "** Add " << d[i][j] << " at (" << i << ", " << j << ") to " << distance_to_cluster[i][j] << endl;
        assoc_self[distance_to_cluster[i][j]] += d[i][j];
      }
    }
  }

  //  Scan the entire distance_to_cluster matrix to calculate Assoc (A, V)
  for (i = 0; i < M; i++) {
    unsigned int cluster_id = distance_to_cluster[i][i];
    for (j = 0; j < M; j++) {
//       if (distance_to_cluster[i][j] == UINT_MAX) {
//         cerr << "Add " << d[i][j] << " at (" << i << ", " << j << ") to " << cluster_id << endl;
        assoc_all[cluster_id] += d[i][j];
//       }
    }
  }

//   for (i = 0; i < total_clusters; i++) {
//     if (valid_cluster[i]) {
//       assoc_all[i] += assoc_self[i];
//     }
//   }

  double score = 0.0;
  for (i = 0; i < total_clusters; i++) {
    if (valid_cluster[i]) {
//       cerr << "-->\t" << score << "\t" << assoc_self[i] << "\t" << assoc_all[i] << endl;
      score += (assoc_self[i] / assoc_all[i]);
    }
  }

  //  Adjust the scores so that we increase it based on the number of clusters
  score = score * static_cast<double> (num_clusters);

  //  Free arrays
  for (i = 0; i < M; i++) {
    delete [] distance_to_cluster[i];
  }
  delete [] distance_to_cluster;
  delete [] assoc_self;
  delete [] assoc_all;

  setScore1 (score);
  setScore2 (score);
  setCombinedScore (score);

  if (debug) {
    cerr << "====>\t" << getScore1 () << "\t" << getScore2 () << "\t" << getCombinedScore () << endl;
//     cerr << "Scoring done!" << endl;
  }

  return;
}


//!  Calculate the MST score based on the normalized association of Shi and Malik (2000)
/*!
     See SCORE::scoreNormalizedAssocOrig for a description.  Only difference is that
     the score is not multiplied by the number of clusters.
*/
void SCORE::scoreNormalizedAssocOrig (unsigned int M, vector<CLUSTER> *clusters, double **d, bool debug) {
  unsigned int i = 0;  //  Position i
  unsigned int j = 0;  //  Position j
  unsigned int max_items = 0;
  unsigned int expt_i = 0;  //  Experiment ID at position i
  unsigned int expt_j = 0;  //  Experiment ID at position j

  unsigned int k = 0;  //  Cluster ID

  //  (i, j) is the cluster ID the edge is in
  //  Map distance (i, j) to its cluster
  unsigned int **distance_to_cluster;
  double *assoc_self = NULL;
  double *assoc_all = NULL;
  vector<bool> valid_cluster;

  unsigned int total_clusters = clusters -> size ();
  unsigned int num_clusters = 0;

  //  Initialize the distance_to_cluster
  distance_to_cluster = new unsigned int * [M];
  for (i = 0; i < M; i++) {
    distance_to_cluster[i] = new unsigned int[M];
    for (j = 0; j < M; j++) {
      //  Initialize with UINT_MAX to indicate inter-cluster edges
      distance_to_cluster[i][j] = UINT_MAX;
    }
  }

  //  For each cluster
  for (k = 0; k < total_clusters; k++) {
    if (!clusters -> at (k).haveAncestors ()) {
      //  Everything within this cluster
      vector<unsigned int> items = clusters -> at (k).getItems ();
      max_items = items.size ();

//        cerr << "Cluster has " << items.size () << " items." << endl;

      // For every pair of items
      for (i = 0; i < max_items ; i++) {
        expt_i = items[i];
        for (j = i; j < max_items; j++) {
          expt_j = items[j];

          //  Indicate that the pair are in the same cluster

//           cerr << "\t\tTRUE:  (" << expt_i << ", " << expt_j << ") --> " << k << endl;

          distance_to_cluster[expt_i][expt_j] = k;
          distance_to_cluster[expt_j][expt_i] = k;
        }
      }
      valid_cluster.push_back (true);
      num_clusters++;
    }
    else {
      valid_cluster.push_back (false);
    }
  }

  assoc_self = new double[total_clusters];
  assoc_all = new double[total_clusters];

  //  Initialize the arrays
  for (i = 0; i < total_clusters; i++) {
    assoc_self[i] = 0.0;
    assoc_all[i] = 0.0;
  }

  //  Scan only the upper-triangular distance_to_cluster matrix to find
  //  the intra-cluster edges  (i.e., Assoc (A, A))
  for (i = 0; i < M; i++) {
    for (j = i + 1; j < M; j++) {
      if (distance_to_cluster[i][j] != UINT_MAX) {
        //  Intra-cluster edge
//         cerr << "** Add " << d[i][j] << " at (" << i << ", " << j << ") to " << distance_to_cluster[i][j] << endl;
        assoc_self[distance_to_cluster[i][j]] += d[i][j];
      }
    }
  }

  //  Scan the entire distance_to_cluster matrix to calculate Assoc (A, V)
  for (i = 0; i < M; i++) {
    unsigned int cluster_id = distance_to_cluster[i][i];
    for (j = 0; j < M; j++) {
//       if (distance_to_cluster[i][j] == UINT_MAX) {
//         cerr << "Add " << d[i][j] << " at (" << i << ", " << j << ") to " << cluster_id << endl;
        assoc_all[cluster_id] += d[i][j];
//       }
    }
  }

//   for (i = 0; i < total_clusters; i++) {
//     if (valid_cluster[i]) {
//       assoc_all[i] += assoc_self[i];
//     }
//   }

  double score = 0.0;
  for (i = 0; i < total_clusters; i++) {
    if (valid_cluster[i]) {
//       cerr << "-->\t" << score << "\t" << assoc_self[i] << "\t" << assoc_all[i] << endl;
      score += (assoc_self[i] / assoc_all[i]);
    }
  }

  //  Adjust the scores so that we increase it based on the number of clusters
//   score = score * static_cast<double> (num_clusters);

  //  Free arrays
  for (i = 0; i < M; i++) {
    delete [] distance_to_cluster[i];
  }
  delete [] distance_to_cluster;
  delete [] assoc_self;
  delete [] assoc_all;

  setScore1 (score);
  setScore2 (score);
  setCombinedScore (score);

  if (debug) {
    cerr << "====>\t" << getScore1 () << "\t" << getScore2 () << "\t" << getCombinedScore () << endl;
//     cerr << "Scoring done!" << endl;
  }

  return;
}
