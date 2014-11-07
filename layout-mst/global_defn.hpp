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
    \file global_defn.hpp
    Header file for global definitions
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: global_defn.hpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#ifndef GLOBAL_DEFN_HPP
#define GLOBAL_DEFN_HPP

//!  Spacing for aligning the verbose output (in characters).
#define VERBOSE_WIDTH 45

//!  The default filename for the configuration file.
#define CFG_FILENAME "layout-mst.cfg"

//!  The number of fields in the file of scores
#define SCORES_FIELDS 6

//!  Set to 1 to enable large nodes (set LARGE_NODE_INCREASE and LARGE_NODE_THRESH); OFF by default!
#define ENABLE_LARGE_NODE 0

//!  Size of node increase (as a multiplier)
#define LARGE_NODE_INCREASE 2

//!  Threshold for when a node increases by LARGE_NODE_INCREASE; 0.25 means that if a node has 25% of the total number of experiments, enlarge it
#define LARGE_NODE_THRESH 0.25

//!  File extension for the file of edges
#define EDGES_FILE_EXTENSION ".edges"

//!  File extension for the file of nodes
#define NODES_FILE_EXTENSION ".nodes"

//!  File extension for the intermediate GraphViz file
#define DOT_FILE_EXTENSION ".dot"

//!  File extension for the GraphViz file
#define GV_FILE_EXTENSION ".graphviz"

//!  File extension for the preview GraphViz file
#define GVPV_FILE_EXTENSION "-pv.graphviz"

//!  Default base URL
#define DEFAULT_URL "http://localhost/"

//!  Default image width (cm); A4
#define DEFAULT_WIDTH 21.0

//!  Default image height (cm); A4
#define DEFAULT_HEIGHT 29.7

//!  Default preview image width (cm); A10
#define DEFAULT_PWIDTH 2.6

//!  Default preview image height (cm); A10
#define DEFAULT_PHEIGHT 3.7

//!  Minimum dimension (cm)
#define MIN_DIM 2

//!  Maximum dimension (cm)
#define MAX_DIM 30

//!  Default percentage is everything (100); this is an unsigned int
#define DEFAULT_PERCENT 100

//!  Default resolution (DPI)
#define DEFAULT_DPI 96

//!  Minimum resolution (DPI)
#define MIN_DPI 48

//!  Maximum resolution (DPI)
#define MAX_DPI 900

//!  Default fontsize
#define DEFAULT_FONTSIZE 12

//!  Minimum fontsize
#define MIN_FONTSIZE 8

//!  Maximum fontsize
#define MAX_FONTSIZE 24

//!  Conversion from cm to inches, which is used by GraphViz
#define SCALE_FACTOR 2.54

//!  The type of output file
enum FILETYPE {
  /*! File type not yet set */ FILETYPE_UNSET,
  /*! Final Graphviz file */ FILETYPE_GV,
  /*! Intermediate Graphviz/DOT file */ FILETYPE_DOT,
  /*! PNG file */ FILETYPE_PNG,
  /*! SVG file */ FILETYPE_SVG,
  /*! Client-side map file */ FILETYPE_CMAP,
  /*! Postscript file */ FILETYPE_PS
};

#endif

