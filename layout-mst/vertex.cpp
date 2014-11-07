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
    \file vertex.cpp
    Member functions for VERTEX class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: vertex.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <string>

using namespace std;

#include "vertex.hpp"

//!  Default constructor that takes no arguments
VERTEX::VERTEX ()
  : name (""),
    colour (""),
    shape (""),
    updated (false),
    components (1),
    x (0),
    y (0),
    width (0.0),
    height (0.0)
{
}

//!  Default constructor that takes four arguments
/*!
    \param arg1 Vertex name
    \param arg2 Vertex colour
    \param arg3 Vertex shape
    \param arg4 Number of components in node
*/
VERTEX::VERTEX (string arg1, string arg2, string arg3, unsigned int arg4)
  : name (arg1),
    colour (arg2),
    shape (arg3),
    updated (false),
    components (arg4),
    x (0),
    y (0),
    width (0.0),
    height (0.0)
{
}

//!  Set the vertex name
void VERTEX::setName (string arg) {
  name = arg;
}

//!  Get the vertex name
string VERTEX::getName () const {
  return name;
}

//!  Set the vertex colour
void VERTEX::setColour (string arg) {
  colour = arg;
}

//!  Get the vertex colour
string VERTEX::getColour () const {
  return colour;
}

//!  Set the vertex shape
void VERTEX::setShape (string arg) {
  shape = arg;
}

//!  Get the vertex shape
string VERTEX::getShape () const {
  return shape;
}

//!  Set the number of components
void VERTEX::setComponents (unsigned int arg) {
  components = arg;
}

//!  Get the number of components
unsigned int VERTEX::getComponents () const {
  return components;
}

//!  Change the update status of the object
void VERTEX::setUpdated (bool arg) {
  updated = arg;
}

//!  Check if the object has its coordinates updated
bool VERTEX::getUpdated () const {
  return updated;
}


//!  Set the x coordinate
void VERTEX::setX (unsigned int arg) {
  x = arg;
}

//!  Get the x coordinate
unsigned int VERTEX::getX () const {
  return x;
}

//!  Set the y coordinate
void VERTEX::setY (unsigned int arg) {
  y = arg;
}

//!  Get the y coordinate
unsigned int VERTEX::getY () const {
  return y;
}

//!  Set the node width
void VERTEX::setWidth (float arg) {
  width = arg;
}

//!  Get the node width
float VERTEX::getWidth () const {
  return width;
}

//!  Set the node height
void VERTEX::setHeight (float arg) {
  height = arg;
}

//!  Get the node height
float VERTEX::getHeight () const {
  return height;
}

