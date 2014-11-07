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
    \file vect_spear.cpp
    Member functions for SPEARMAN class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: vect_spear.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/

#include <vector>
#include <iostream>

using namespace std;

#include "vect_spear.hpp"

//!  Constructor that takes no arguments
SPEARMAN::SPEARMAN ()
  : key (0),
    value (0),
    origpos (0),
    rank (0)
{
}

//!  Constructor that takes three arguments (key/value, original position, rank)
/*!
     \param v The key (the value we sort on) or the value of the node (i.e., the expression level)
     \param p The original position in the vector
     \param r The rank of this node respective to the other values
*/
SPEARMAN::SPEARMAN (double v, unsigned int p, double r) {
  key = v;
  value = v;
  origpos = p;
  rank = r;
}

//!  Set the value
void SPEARMAN::setValue (double v) {
  value = v;
}

//!  Set the original position
void SPEARMAN::setOrigPos (unsigned int p) {
  origpos = p;
}

//!  Set the rank
void SPEARMAN::setRank (double r) {
  rank = r;
}

//!  Get the value
double SPEARMAN::getValue () const {
  return value;
}

//!  Get the original position
unsigned int SPEARMAN::getOrigPos () const {
  return origpos;
}

//!  Get the rank
double SPEARMAN::getRank () const {
  return rank;
}

//!  Overloaded operator for SPEARMAN nodes (less than)
bool SPEARMAN::operator< (const SPEARMAN &arg) const {
        return (key < arg.key);
}

//!  Copy the original position to the key
void SPEARMAN::copyOrigPosToKey () {
  key = origpos;
}


