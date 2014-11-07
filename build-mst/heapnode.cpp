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
    \file heapnode.cpp
    Additional member functions for HEAPNODE class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: heapnode.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


using namespace std;

#include <climits>
#include <cfloat>

#include "heapnode.hpp"

//!  Default HEAPNODE constructor; should never be called (not needed).
HEAPNODE::HEAPNODE ()
  : left (UINT_MAX),
    right (UINT_MAX),
    score (DBL_MAX)
{
}

//!  Constructor for a HEAPNODE with three arguments.
/*!
     \param arg1 The left child
     \param arg2 The right child
     \param arg3 The score
*/
HEAPNODE::HEAPNODE (unsigned int arg1, unsigned int arg2, double arg3)
  : left (arg1),
    right (arg2),
    score (arg3)
{
}

//!  Get the left child
unsigned int HEAPNODE::getLeft () const {
  return left;
}

//!  Get the right child
unsigned int HEAPNODE::getRight () const {
  return right;
}

//!  Get the score
double HEAPNODE::getScore () const {
  return score;
}

//!  Overloaded operator for HEAPNODEs (less than)
bool HEAPNODE::operator< (const HEAPNODE &arg) const {
  return (score < arg.score);
}

//!  Overloaded operator for HEAPNODEs (greater than)
bool HEAPNODE::operator> (const HEAPNODE &arg) const {
  return (score > arg.score);
}



