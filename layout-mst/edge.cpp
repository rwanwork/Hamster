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
    \file edge.cpp
    Member functions for EDGE class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: edge.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <string>

using namespace std;

#include "edge.hpp"

//!  Default constructor that takes no arguments
EDGE::EDGE ()
  : start (""),
    end (""),
    weight (0.0)
{
}

//!  Default constructor that takes three arguments
/*!
    \param arg1 Edge start point
    \param arg2 Edge end point
    \param arg3 Edge weight
*/
EDGE::EDGE (string arg1, string arg2, double arg3)
  : start (arg1),
    end (arg2),
    weight (arg3)
{
}

//!  Set the start vertex name
void EDGE::setStart (string arg) {
  start = arg;
}

//!  Get the start vertex name
string EDGE::getStart () const {
  return start;
}

//!  Set the end vertex name
void EDGE::setEnd (string arg) {
  end = arg;
}

//!  Get the end vertex name
string EDGE::getEnd () const {
  return end;
}

//!  Set the edge weight
void EDGE::setWeight (double arg) {
  weight = arg;
}

//!  Get the edge weight
double EDGE::getWeight () const {
  return weight;
}

