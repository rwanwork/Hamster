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
    \file layout_mst.cpp
    Member functions for LAYOUTMST class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: layout_mst.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <string>
#include <vector>

#include <cstdlib>  //  NULL and srand

using namespace std;

#include "LayoutMSTConfig.hpp"

#if HAVE_MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
using boost::mpi::environment;
using boost::mpi::communicator;
#endif

#include "global_defn.hpp"
#include "check.hpp"
#include "packet.hpp"
#include "score.hpp"
#include "vertex.hpp"
#include "edge.hpp"
#include "layout_mst.hpp"

//!  Default constructor that takes no arguments
LAYOUTMST::LAYOUTMST ()
  : debug_flag (false),
    verbose_flag (false),
    preview_flag (false),
    fixed_pos (false),
    path (""),
    url (""),
    width (0.0),
    height (0.0),
    pwidth (0.0),
    pheight (0.0),
    outtype (FILETYPE_UNSET),
    percent (100),
    dpi (0),
    spline_flag (false),
    fontsize (0),
    total_iter (0),
    scores_fn (""),
    scores (),
    all_workunits (),
    my_workunit (),
#if HAVE_MPI
    env (NULL),
    comm (NULL),
#endif
    rank (0),
    world_size (0)
{
  //  Set the random seed using the current time
  srand (time (NULL));
}

//////////////////////////////////////////////////
//!  Set whether or not debugging output is required
void LAYOUTMST::setDebug (bool arg) {
  debug_flag = arg;
}

//!  Get the debug setting
bool LAYOUTMST::getDebug () const {
  return debug_flag;
}

//!  Set whether or not verbose output is required
void LAYOUTMST::setVerbose (bool arg) {
  verbose_flag = arg;
}

//!  Get the verbose setting
bool LAYOUTMST::getVerbose () const {
  return verbose_flag;
}

//!  Set whether or not preview images are required
void LAYOUTMST::setPreview (bool arg) {
  preview_flag = arg;
}

//!  Get the preview setting
bool LAYOUTMST::getPreview () const {
  return preview_flag;
}

//!  Set whether or not node positions are fixed using the previous MST
void LAYOUTMST::setFixedPos (bool arg) {
  fixed_pos = arg;
}

//!  Get the fixed node position setting
bool LAYOUTMST::getFixedPos () const {
  return fixed_pos;
}

//!  Set the input/output path (the path where files are)
void LAYOUTMST::setPath (string arg) {
  string tmp = sanitizePath (arg);

  path = "";
  if (tmp.length () != 0) {
    path = tmp;
  }
}

//!  Get the output path
string LAYOUTMST::getPath () const {
  return path;
}

//!  Set the URL for image maps
void LAYOUTMST::setURL (string arg) {
  string tmp = sanitizeURL (arg);

  url = "";
  if (tmp.length () != 0) {
    url = arg;
  }
}
//!  Get the URL
string LAYOUTMST::getURL () const {
  return url;
}

//!  Set the width of the images
void LAYOUTMST::setWidth (double arg) {
  if ((arg < MIN_DIM) || (arg > MAX_DIM)) {
    arg = DEFAULT_WIDTH;
  }
  width = arg;
}

//!  Get the width of the images
double LAYOUTMST::getWidth () const {
  return width;
}

//!  Set the height of the images
void LAYOUTMST::setHeight (double arg) {
  if ((arg < MIN_DIM) || (arg > MAX_DIM)) {
    arg = DEFAULT_HEIGHT;
  }
  height = arg;
}

//!  Get the height of the images
double LAYOUTMST::getHeight () const {
  return height;
}

//!  Set the width of the preview images
void LAYOUTMST::setPWidth (double arg) {
  if ((arg < MIN_DIM) || (arg > MAX_DIM)) {
    arg = DEFAULT_PWIDTH;
  }
  pwidth = arg;
}

//!  Get the width of the preview images
double LAYOUTMST::getPWidth () const {
  return pwidth;
}

//!  Set the height of the preview images
void LAYOUTMST::setPHeight (double arg) {
  if ((arg < MIN_DIM) || (arg > MAX_DIM)) {
    arg = DEFAULT_PHEIGHT;
  }
  pheight = arg;
}

//!  Get the height of the preview images
double LAYOUTMST::getPHeight () const {
  return pheight;
}

//!  Set the type of output image
void LAYOUTMST::setOuttype (enum FILETYPE arg) {
  outtype = arg;
}

//!  Get the type of output image
enum FILETYPE LAYOUTMST::getOuttype () const {
  return outtype;
}

//!  Set the percent
void LAYOUTMST::setPercent (unsigned int arg) {
  //  Percent must be between 1 and 100
  if ((arg == 0) || (arg > 100)) {
    arg = DEFAULT_PERCENT;
  }
  percent = arg;
}

//!  Get the percent
unsigned int LAYOUTMST::getPercent () const {
  return percent;
}

//!  Set the DPI
void LAYOUTMST::setDPI (unsigned int arg) {
  if ((arg < MIN_DPI) || (arg > MAX_DPI)) {
    arg = DEFAULT_DPI;
  }
  dpi = arg;
}

//!  Get the DPI
unsigned int LAYOUTMST::getDPI () const {
  return dpi;
}

//!  Set the spline flag
void LAYOUTMST::setSpline (bool arg) {
  spline_flag = arg;
}

//!  Get the spline flag
bool LAYOUTMST::getSpline () const {
  return spline_flag;
}

//!  Set the font size
void LAYOUTMST::setFontsize (unsigned int arg) {
  if ((arg < MIN_FONTSIZE) || (arg > MAX_FONTSIZE)) {
    arg = DEFAULT_FONTSIZE;
  }
  fontsize = arg;
}

//!  Get the font size
unsigned int LAYOUTMST::getFontsize () const {
  return fontsize;
}

//////////////////////////////////////////////////
//!  Set the total number of iterations
void LAYOUTMST::setTotalIter (unsigned int arg) {
  total_iter = arg;
}

//!  Get the total number of iterations
unsigned int LAYOUTMST::getTotalIter () const {
  return total_iter;
}

//!  Set the scores filename
void LAYOUTMST::setScoresFn (string arg) {
  string tmp = sanitizeFilename (arg);

  scores_fn = "";
  if (tmp.length () != 0) {
    scores_fn = tmp;
  }
}

//!  Get the scores filename
string LAYOUTMST::getScoresFn () const {
  return scores_fn;
}

//////////////////////////////////////////////////
//!  Set the workunits for the current process to do
void LAYOUTMST::setMyWorkunit (vector<unsigned int> arg) {
  my_workunit = arg;
}

//!  Get the workunits for the current process to do
vector<unsigned int> LAYOUTMST::getMyWorkunit () const {
  return my_workunit;
}

//!  Set the workunits for the child process arg1 to do (used by the primary processor only)
/*!
     \param arg1 The rank of the process that will get this workunit
     \param arg2 The workunit
*/
void LAYOUTMST::setAllWorkunits (unsigned int arg1, vector<unsigned int> arg2) {
  all_workunits[arg1] = arg2;
}

//!  Get the workunits for the child process arg1 to do from the primary processor
/*!
     \param arg The rank of the process that will get this workunit
*/
vector<unsigned int> LAYOUTMST::getAllWorkunits (unsigned int arg) const {
  return all_workunits[arg];
}

//////////////////////////////////////////////////
#if HAVE_MPI
//!  Set the MPI environment
void LAYOUTMST::setEnv (environment *arg) {
  env = arg;
}

//!  Get the MPI environment
environment *LAYOUTMST::getEnv () const {
  return env;
}

//!  Set the MPI communicator
void LAYOUTMST::setComm (communicator *arg) {
  comm = arg;
}

//!  Get the MPI communicator
communicator *LAYOUTMST::getComm () const {
  return comm;
}
#endif

//!  Set the rank for the current process
void LAYOUTMST::setRank (unsigned int arg) {
  rank = arg;
}

//!  Get the rank for the current process
unsigned int LAYOUTMST::getRank () const {
  return rank;
}

//!  Get the total number of processes
void LAYOUTMST::setWorldSize (unsigned int arg) {
  world_size = arg;
}

//!  Set the total number of processes
unsigned int LAYOUTMST::getWorldSize () const {
  return world_size;
}


