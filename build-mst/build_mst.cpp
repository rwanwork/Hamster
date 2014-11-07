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
    \file build_mst.cpp
    Member functions for BUILDMST class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: build_mst.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>
#include <string>
#include <vector>
#include <queue>  //  priority_queue

#include <cstdlib>  //  NULL and srand

using namespace std;

#include "global_defn.hpp"
#include "check.hpp"
#include "heapnode.hpp"
#include "vect_spear.hpp"
#include "vect.hpp"
#include "cluster.hpp"
#include "score.hpp"
#include "build_mst.hpp"

//!  Constructor that takes no arguments
BUILDMST::BUILDMST ()
  : debug_flag (false),
    verbose_flag (false),
    distance (DIST_EUC),
    linkage (LINK_SINGLE),
    scoring (SCORE_GAPS),
    centroid (DIST_EUC),
    attr_fn (""),
    microarray_fn (""),
    path (""),
    M (0),
    N (0)
{
  //  Set the random seed using the current time
  srand (time (NULL));
}


//!  Set whether or not debugging output is required
void BUILDMST::setDebug (bool arg) {
  debug_flag = arg;
}

//!  Get the debug setting
bool BUILDMST::getDebug () const {
  return debug_flag;
}

//!  Set whether or not verbose output is required
void BUILDMST::setVerbose (bool arg) {
  verbose_flag = arg;
}

//!  Get the verbose setting
bool BUILDMST::getVerbose () const {
  return verbose_flag;
}

//!  Set the distance method
void BUILDMST::setDistance (DIST_METHOD arg) {
  distance = arg;
}

//!  Get the distance setting
DIST_METHOD BUILDMST::getDistance () const {
  return distance;
}

//!  Set the linkage method
void BUILDMST::setLinkage (LINK_METHOD arg) {
  linkage = arg;
}

//!  Get the linkage setting
LINK_METHOD BUILDMST::getLinkage () const {
  return linkage;
}

//!  Set the scoring method
void BUILDMST::setScoreMethod (SCORE_METHOD arg) {
  scoring = arg;
}

//!  Get the score method setting
SCORE_METHOD BUILDMST::getScoreMethod () const {
  return scoring;
}

//!  Set the distance method for centroid linkage
void BUILDMST::setCentroid (DIST_METHOD arg) {
  centroid = arg;
}

//!  Get the distance setting for centroid linkage
DIST_METHOD BUILDMST::getCentroid () const {
  return centroid;
}

//!  Set the microarray filename
void BUILDMST::setMicroarrayFn (string arg) {
  string tmp = sanitizeFilename (arg);

  microarray_fn = "";
  if (tmp.length () != 0) {
    microarray_fn = tmp;
  }
}

//!  Get the microarray filename
string BUILDMST::getMicroarrayFn () const {
  return microarray_fn;
}

//!  Set the attribute filename
void BUILDMST::setAttrFn (string arg) {
  string tmp = sanitizeFilename (arg);

  attr_fn = "";
  if (tmp.length () != 0) {
    attr_fn = tmp;
  }
}

//!  Get the attribute filename
string BUILDMST::getAttrFn () const {
  return attr_fn;
}

//!  Set the output path (the path where files will be written to)
void BUILDMST::setPath (string arg) {
  string tmp = sanitizePath (arg);

  path = "";
  if (tmp.length () != 0) {
    path = tmp;
  }
}
//!  Get the output path
string BUILDMST::getPath () const {
  return path;
}

//!  Set M -- the number of objects (experiments or rows) in the data set
void BUILDMST::setM (unsigned int arg) {
  M = arg;
}

//!  Get the number of objects in the data set
unsigned int BUILDMST::getM () const {
  return M;
}

//!  Set N -- the number of attributes (probes or columns) in the data set
void BUILDMST::setN (unsigned int arg) {
  N = arg;
}

//!  Get the number of attributes in the data set
unsigned int BUILDMST::getN () const {
  return N;
}

