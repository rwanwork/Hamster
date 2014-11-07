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
    \file cluster_link.cpp
    Additional member functions for CLUSTER class definition
      Functions for calculating linkages
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: cluster_link.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <string>
#include <vector>

#include <cfloat>  //  DBL_MAX

using namespace std;

#include "global_defn.hpp"
#include "vect_spear.hpp"
#include "vect.hpp"
#include "cluster.hpp"

//!  The single linkage between this cluster and another one
double CLUSTER::linkSingle (CLUSTER *other, double **d) {
  double score = 0;

  vector<unsigned int> my_items = getItems ();
  unsigned int i = 0;  //  Position i
  unsigned int max_i = my_items.size ();
  unsigned int expt_i = 0;  //  Experiment ID at position i

  vector<unsigned int> other_items = other -> getItems ();
  unsigned int j = 0;  //  Position j
  unsigned int max_j = other_items.size ();
  unsigned int expt_j = 0;  //  Experiment ID at position j

  //  Initialization to the first pair
  expt_i = my_items[0];
  expt_j = other_items[0];
  score = d[expt_i][expt_j];

  for (i = 0; i < max_i; i++) {
    expt_i = my_items[i];
    for (j = 0; j < max_j; j++) {
      expt_j = other_items[j];
      //  Take distance if smaller than the current score
      if (d[expt_i][expt_j] < score) {
        score = d[expt_i][expt_j];
      }
    }
  }

  return score;
}


//!  The average linkage between this cluster and another one
double CLUSTER::linkAverage (CLUSTER *other, double **d) {
  double score = 0;
  unsigned int count = 0;

  vector<unsigned int> my_items = getItems ();
  unsigned int i = 0;  //  Position i
  unsigned int max_i = my_items.size ();
  unsigned int expt_i = 0;  //  Experiment ID at position i

  vector<unsigned int> other_items = other -> getItems ();
  unsigned int j = 0;  //  Position j
  unsigned int max_j = other_items.size ();
  unsigned int expt_j = 0;  //  Experiment ID at position j

  for (i = 0; i < max_i; i++) {
    expt_i = my_items[i];
    for (j = 0; j < max_j; j++) {
      expt_j = other_items[j];
      if (d[expt_i][expt_j] < score) {
        score += d[expt_i][expt_j];
        count++;
      }
    }
  }

  if (count == 0) {
    score = DBL_MAX;
  }
  else {
    score = score / static_cast<double> (count);
  }

  return score;
}

//!  The complete linkage between this cluster and another one
double CLUSTER::linkComplete (CLUSTER *other, double **d) {
  double score = 0;

  vector<unsigned int> my_items = getItems ();
  unsigned int i = 0;  //  Position i
  unsigned int max_i = my_items.size ();
  unsigned int expt_i = 0;  //  Experiment ID at position i

  vector<unsigned int> other_items = other -> getItems ();
  unsigned int j = 0;  //  Position j
  unsigned int max_j = other_items.size ();
  unsigned int expt_j = 0;  //  Experiment ID at position j

  //  Initialization to the first pair
  expt_i = my_items[0];
  expt_j = other_items[0];
  score = d[expt_i][expt_j];

  for (i = 0; i < max_i; i++) {
    expt_i = my_items[i];
    for (j = 0; j < max_j; j++) {
      expt_j = other_items[j];
      //  Take distance if larger than the current score
      if (d[expt_i][expt_j] > score) {
        score = d[expt_i][expt_j];
      }
    }
  }

  return score;
}

//!  The centroid linkage between this cluster and another one
/*!
     The dissimilarity between two centroid vectors depends
     on the distance parameter chosen by the user.  That is,
     the Euclidean distance is not always used.  Always
     using the Euclidean distance may be preferred, but
     that is up to the user.
*/
double CLUSTER::linkCentroid (CLUSTER *other, vector<VECT> *data, enum DIST_METHOD distance) {
  double score = 0;
  VECT temp = other -> getCentroid ();

  switch (distance) {
    case DIST_EUC :
      score = getCentroid ().simEuc (&temp);
      break;
    case DIST_MAN :
      score = getCentroid ().simMan (&temp);
      break;
    case DIST_PEAR :
      score = getCentroid ().simPear (&temp);
      break;
    case DIST_SPEAR :
      score = getCentroid ().simSpear (&temp);
      break;
    default :
      break;
  }

  return score;
}

