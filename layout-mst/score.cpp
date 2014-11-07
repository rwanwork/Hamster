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


#include <string>
#include <vector>

#include <climits>  //  UINT_MAX
#include <cmath>

using namespace std;

#include "global_defn.hpp"
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

//!  Constructor that takes five arguments
/*!
    \param arg1 ID of the node
    \param arg2 ID of the left cluster
    \param arg3 ID of the right cluster
    \param arg4 Score 1
    \param arg5 Score 2
    \param arg6 Combined score
*/
SCORE::SCORE (unsigned int arg1, unsigned int arg2, unsigned int arg3, double arg4, double arg5, double arg6)
  : id (arg1),
    left (arg2),
    right (arg3),
    score1 (arg4),
    score2 (arg5),
    combined (arg6)
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

//!  Overloaded operator for SCORE nodes (less than)
bool SCORE::operator< (const SCORE &arg) const {
  return (combined < arg.combined);
}

//!  Overloaded operator for SCORE nodes (greater than)
bool SCORE::operator> (const SCORE &arg) const {
  return (combined > arg.combined);
}
