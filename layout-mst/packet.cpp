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
    \file packet.cpp
    Member functions for PACKET class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: packet.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <string>
#include <vector>

#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

using namespace std;

#include "global_defn.hpp"
#include "packet.hpp"

//!  Default constructor that takes no arguments
PACKET::PACKET ()
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
    my_workunit ()
{
}

//!  Set whether or not debugging output is required
void PACKET::setDebug (bool arg) {
  debug_flag = arg;
}

//!  Get the debug setting
bool PACKET::getDebug () const {
  return debug_flag;
}

//!  Set whether or not verbose output is required
void PACKET::setVerbose (bool arg) {
  verbose_flag = arg;
}

//!  Get the verbose setting
bool PACKET::getVerbose () const {
  return verbose_flag;
}

//!  Set whether or not verbose output is required
void PACKET::setPreview (bool arg) {
  preview_flag = arg;
}

//!  Get the verbose setting
bool PACKET::getPreview () const {
  return preview_flag;
}

//!  Set whether or not node positions are fixed using the previous MST
void PACKET::setFixedPos (bool arg) {
  fixed_pos = arg;
}

//!  Get the fixed node position setting
bool PACKET::getFixedPos () const {
  return fixed_pos;
}

//!  Set the path to the files
void PACKET::setPath (string arg) {
  path = arg;
}

//!  Get the path
string PACKET::getPath () const {
  return path;
}

//!  Set the URL for image maps
void PACKET::setURL (string arg) {
  url = arg;
}

//!  Get the URL
string PACKET::getURL () const {
  return url;
}

//!  Set the width of the images
void PACKET::setWidth (double arg) {
  width = arg;
}

//!  Get the width of the images
double PACKET::getWidth () const {
  return width;
}

//!  Set the height of the images
void PACKET::setHeight (double arg) {
  height = arg;
}

//!  Get the height of the images
double PACKET::getHeight () const {
  return height;
}

//!  Set the width of the preview images
void PACKET::setPWidth (double arg) {
  pwidth = arg;
}

//!  Get the width of the preview images
double PACKET::getPWidth () const {
  return pwidth;
}

//!  Set the height of the preview images
void PACKET::setPHeight (double arg) {
  pheight = arg;
}

//!  Get the height of the preview images
double PACKET::getPHeight () const {
  return pheight;
}

//!  Set the type of output image
void PACKET::setOuttype (enum FILETYPE arg) {
  outtype = arg;
}

//!  Get the type of output image
enum FILETYPE PACKET::getOuttype () const {
  return outtype;
}

//!  Set the percent
void PACKET::setPercent (unsigned int arg) {
  percent = arg;
}

//!  Get the percent
unsigned int PACKET::getPercent () const {
  return percent;
}

//!  Set the DPI
void PACKET::setDPI (unsigned int arg) {
  dpi = arg;
}

//!  Get the DPI
unsigned int PACKET::getDPI () const {
  return dpi;
}

//!  Set the spline flag
void PACKET::setSpline (bool arg) {
  spline_flag = arg;
}

//!  Get the spline flag
bool PACKET::getSpline () const {
  return spline_flag;
}

//!  Set the font size
void PACKET::setFontsize (unsigned int arg) {
  fontsize = arg;
}

//!  Get the font size
unsigned int PACKET::getFontsize () const {
  return fontsize;
}

//!  Set the workunits for the current process to do
void PACKET::setMyWorkunit (vector<unsigned int> arg) {
  my_workunit = arg;
}

//!  Get the workunits for the current process to do
vector<unsigned int> PACKET::getMyWorkunit () const {
  return my_workunit;
}

