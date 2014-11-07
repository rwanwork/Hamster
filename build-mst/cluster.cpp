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
    \file cluster.cpp
    Member functions for CLUSTER class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: cluster.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>  //  cerr
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>

using namespace std;

#include <cstdlib>  //  exit, EXIT_FAILURE

#include "global_defn.hpp"
#include "vect_spear.hpp"
#include "vect.hpp"
#include "cluster.hpp"

//!  Default constructor; should never be called (not needed)
CLUSTER::CLUSTER ()
  : id (0),
    name (""),
    colour (""),
    shape (""),
    items (),
    ancestors (false),
    centroid ()
{
}

//!  Constructor for creating the initial clusters
/*!
     \param arg_id The integral ID for this cluster
     \param arg1 The name of this cluster
     \param arg2 The colour of this cluster
     \param arg3 The shape of this cluster
     \param arg4 The VECT object for this cluster of one item
*/
CLUSTER::CLUSTER (unsigned int arg_id, string arg1, string arg2, string arg3, VECT arg4)
  : id (arg_id),
    name (arg1),
    colour (arg2),
    shape (arg3),
    items (),
    ancestors (false),
    centroid (arg4)
{
  //  Add the item in; cluster of size 1
  items.push_back (arg_id);
}

//!  Constructor for merging two clusters
/*!
     \param arg_id The integral ID for this cluster
     \param arg1 The first cluster
     \param arg2 The second cluster
     \param arg3 The linkage method being used
     \param arg4 Number of experiments in microarray (for naming)

     This constructor sets the attributes of the cluster and
     if centroid linkage is used, calls formCentroid () to build
     a centroid vector.

     Note:  There is no order to clusters so arg1 and arg2 can be swapped.
*/
CLUSTER::CLUSTER (unsigned int arg_id, CLUSTER *arg1, CLUSTER *arg2, LINK_METHOD arg3, unsigned int arg4)
  : id (arg_id),
    name (""),
    colour (""),
    shape (""),
    items (),
    ancestors (false),
    centroid ()
{
  vector<unsigned int>::const_iterator start;

  vector<unsigned int> first_items = arg1 -> getItems ();
  vector<unsigned int> second_items = arg2 -> getItems ();

  //  Merge the items
  for (start = second_items.begin (); start != second_items.end (); start++) {
    first_items.push_back (*start);
  }
  items = first_items;

  //  Set the name based on its ID
  setName (string (boost::lexical_cast<std::string>(arg_id - arg4)));

  //  Set the colour of the new cluster
  string str1 = arg1 -> getColour ();
  string str2 = arg2 -> getColour ();
  if (str1 == str2) {
    setColour (str1);
  }
  else {
    setColour (DEFAULT_COLOUR);
  }

  //  Set the shape of the new cluster
  str1 = arg1 -> getShape ();
  str2 = arg2 -> getShape ();
  if (str1 == str2) {
    setShape (str1);
  }
  else {
    setShape (DEFAULT_SHAPE);
  }

  //  Both decendents have ancestors now
  arg1 -> setAncestors ();
  arg2 -> setAncestors ();

  //  Centroid linkage is being used, so we need to form a centroid vector
  if (arg3 == LINK_CENTROID) {
    formCentroid (arg1, arg2);
  }
}

//!  Get the cluster ID
unsigned int CLUSTER::getID () const {
  return id;
}

//!  Set the cluster name
void CLUSTER::setName (string arg) {
  name = arg;
}

//!  Get the cluster name
string CLUSTER::getName () const {
  return name;
}

//!  Set the cluster colour
void CLUSTER::setColour (string arg) {
  colour = arg;
}

//!  Get the cluster colour
string CLUSTER::getColour () const {
  return colour;
}

//!  Set the cluster shape
void CLUSTER::setShape (string arg) {
  shape = arg;
}

//!  Get the cluster shape
string CLUSTER::getShape () const {
  return shape;
}

//!  Get the components that make up a cluster
vector<unsigned int> CLUSTER::getItems () const {
  return items;
}

//!  Set the ancestors variable to TRUE
void CLUSTER::setAncestors () {
  ancestors = true;
}

//!  Create a centroid vector from two clusters
/*!
     The centroid vector is formed by averaging across each feature
     of the vector.  If the value of position i in any vector is NULL,
     then position i of the new column is automatically NULL such
     that the sum of any value with NULL is NULL.  The centroid
     vector is stored as a VECT object.
*/
void CLUSTER::formCentroid (CLUSTER *a, CLUSTER *b) {
  unsigned int a_len = (a -> getCentroid ()).getN ();
  unsigned int b_len = (b -> getCentroid ()).getN ();
  VECT a_vect = a -> getCentroid ();
  VECT b_vect = b -> getCentroid ();

  if (a_len != b_len) {
    cerr << "Error:  Number of columns differ! (" << a_len << ", " << b_len << ")" << endl;
    exit (EXIT_FAILURE);
  }

  //  Resize the vector since we know the size now
  centroid.resize (a_len);

  for (unsigned int i = 0; i < a_len; i++) {
    //  If either is NULL, make this one NULL too
    if (a_vect.isNull (i) || b_vect.isNull (i)) {
      centroid.putExpr (i, 0);
      centroid.putNull (i, true);
    }
    else {
      centroid.putExpr (i, (((a_len * a_vect.getExpr (i)) + (b_len * b_vect.getExpr (i))) / (a_len + b_len)));
      centroid.putNull (i, false);
    }
  }

  return;
}

//!  Get the centroid vector from this cluster
VECT CLUSTER::getCentroid () const {
  return centroid;
}


