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
    \file vect_dist.cpp
    Additional member functions for VECT class definition
      Functions for calculating distances between two vectors
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: vect_dist.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>  //  cerr, endl
#include <string>
#include <vector>
#include <algorithm>  //  sort

#include <cmath>  //  fabs, sqrt
#include <cfloat>  //  DBL_MAX
#include <cstdlib>  //  exit, EXIT_FAILURE

using namespace std;

#include "vect_spear.hpp"
#include "vect.hpp"

//!  The Euclidean distance between this vector and another one
/*!  Function exits if the two vectors are of different dimensions.  */
double VECT::simEuc (VECT *other) {
  unsigned int i = 0;
  double result = 0;
  double temp = 0;
  unsigned int size = 0;
  unsigned int n = nulls.size ();

  //  Ensure both rows are of the same dimensions
  if (getN () != other -> getN ()) {
    cerr << "Error:  Dimensions differ!" << endl;
    exit (EXIT_FAILURE);
  }

  //  First calculate differences between each pair of numbers
  //  and squared, provided both are non-null
  for (i = 0; i < n; i++) {
    if (!(this -> isNull (i)) && !(other -> isNull (i))) {
      temp = (this -> getExpr (i)) - (other -> getExpr (i));
      temp = temp * temp;
      result += temp;
      size++;
    }
  }

  if (size == 0) {
    result = DBL_MAX;
  }
  else {
    //  Sqrt number
    result = sqrt (result);
  }

  return (result);
}


//!  The Manhattan distance between this vector and another one
/*!  Function exits if the two vectors are of different dimensions.  */
double VECT::simMan (VECT *other) {
  unsigned int i = 0;
  double result = 0;
  unsigned int size = 0;
  unsigned int n = nulls.size ();

  //  Ensure both rows are of the same dimensions
  if (getN () != other -> getN ()) {
    cerr << "Error:  Dimensions differ!" << endl;
    exit (EXIT_FAILURE);
  }

  //  First calculate differences between each pair of numbers
  //  and squared, provided both are non-null
  for (i = 0; i < n; i++) {
    if (!(this -> isNull (i)) && !(other -> isNull (i))) {
      result += fabs ((this -> getExpr (i)) - (other -> getExpr (i)));
      size++;
    }
  }

  if (size == 0) {
    result = DBL_MAX;
  }

  return (result);
}


//!  The Pearson correlation coefficient (distance) between this vector and another one
/*!
     The Pearson correlation coefficient (r) is subtracted from 2 to obtain
     a distance whose range is [0, 2] such that 0 means two vectors are highly
     correlated.

     Function exits if the two vectors are of different dimensions.

     Note:  The calculation makes use of the population standard deviation.
*/
double VECT::simPear (VECT *other) {
  unsigned int n = nulls.size ();
  unsigned int i = 0;
  unsigned int n2 = 0;  //  Number of non-null pairs
  double result = 0;
  double sumxy = 0;
  double sumx = 0;
  double sumy = 0;
  double sumx_sqr = 0;
  double sumy_sqr = 0;
  double num = 0;
  double den1 = 0;
  double den2 = 0;

  //  Ensure both rows are of the same dimensions
  if (getN () != other -> getN ()) {
    cerr << "Error:  Dimensions differ!" << endl;
    exit (EXIT_FAILURE);
  }

  //  Calculate average of both genes
  for (i = 0; i < n; i++) {
    if (!(this -> isNull (i)) && !(other -> isNull (i))) {
      sumxy += (this -> getExpr (i)) * (other -> getExpr (i));
      sumx += (this -> getExpr (i));
      sumy += (other -> getExpr (i));
      sumx_sqr += (this -> getExpr (i)) * (this -> getExpr (i));
      sumy_sqr += (other -> getExpr (i)) * (other -> getExpr (i));
      n2++;
    }
  }

  num = (static_cast<double> (n2) * sumxy) - (sumx * sumy);
  den1 = (static_cast<double> (n2) * sumx_sqr) - (sumx * sumx);
  den2 = (static_cast<double> (n2) * sumy_sqr) - (sumy * sumy);

  //  Convert correlation to distance from 0 to 2
  if (den1 * den2 == 0) {
    //  Maximum possible distance
    result = 2.0;
  }
  else {
    result = 1 - num / (sqrt (den1 * den2));
  }

  return (result);
}


//!  The Spearman rank correlation coefficient (distance) between this vector and another one
/*!
     The Spearman rank correlation coefficient is calculated by sorting the two
     vectors by value, enumerating them separately (i.e., assign ranks), and then
     returning them to their original order.  These ranks are passed to the simPear
     function.

     Function exits if the two vectors are of different dimensions.

     The final result is a distance whose range is [0, 2] such that 0 means two
     vectors are highly correlated.
*/
double VECT::simSpear (VECT *other) {
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int k = 0;
  unsigned int n2 = 0;  //  Number of non-null pairs
  double result = 0.0;
  vector<SPEARMAN> myspears;
  vector<SPEARMAN> otherspears;

  //  Ensure both rows are of the same dimensions
  if (getN () != other -> getN ()) {
    cerr << "Error:  Dimensions differ!" << endl;
    exit (EXIT_FAILURE);
  }

  //  Add to spearman node
  for (i = 0; i < getN (); i++) {
    if (!(this -> isNull (i)) && !(other -> isNull (i))) {
      myspears.push_back (SPEARMAN (getExpr (i), n2, 0));
      otherspears.push_back (SPEARMAN (other -> getExpr (i), n2, 0));
      n2++;
    }
  }

  //  Sort vector of spearman nodes by value
  sort (myspears.begin (), myspears.end ());
  sort (otherspears.begin (), otherspears.end ());

  //  Assign rank
  for (i = 0; i < n2; i++) {
    myspears[i].setRank (i);
    otherspears[i].setRank (i);
  }

  //  Handle duplicate ranks for myspears
  i = 0;
  while (i < n2) {
    unsigned int dups = 1;
    double ranksum = myspears[i].getRank ();
    for (j = i + 1; j < n2; j++) {
      if (myspears[i].getValue () != myspears[j].getValue ()) {
        break;
      }
      else {
        ranksum += myspears[j].getRank ();
        dups++;
      }
    }
    for (k = i; k < j; k++) {
      myspears[k].setRank (ranksum / dups);
    }
    i = j;
  }

  //  Handle duplicate ranks for otherspears
  i = 0;
  while (i < n2) {
    unsigned int dups = 1;
    double ranksum = otherspears[i].getRank ();
    for (j = i + 1; j < n2; j++) {
      if (otherspears[i].getValue () != otherspears[j].getValue ()) {
        break;
      }
      else {
        ranksum += otherspears[j].getRank ();
        dups++;
      }
    }
    for (k = i; k < j; k++) {
      otherspears[k].setRank (ranksum / dups);
    }
    i = j;
  }

  //  Copy original positions to key
  for (i = 0; i < n2; i++) {
    myspears[i].copyOrigPosToKey ();
    otherspears[i].copyOrigPosToKey ();
  }

  //  Sort by original positions
  sort (myspears.begin (), myspears.end ());
  sort (otherspears.begin (), otherspears.end ());

  //  Copy SPEARMAN nodes to VECT objects so that we can apply simPear to it
  VECT myrow = VECT (myspears);
  VECT otherrow = VECT (otherspears);

  //  Calculate Pearson correlation
  result = myrow.simPear (&otherrow);

  return result;
}

