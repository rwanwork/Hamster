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
#define VERBOSE_WIDTH 35

//!  The default filename for the configuration file.
#define CFG_FILENAME "build-mst.cfg"

//!  The default filename for the list of scores
#define SCORES_FILENAME "summary.txt"

//!  File extension for the file of edges
#define EDGES_FILE_EXTENSION ".edges"

//!  File extension for the file of nodes
#define NODES_FILE_EXTENSION ".nodes"

//!  The default node colour.
#define DEFAULT_COLOUR "gray"

//!  The default node shape.
#define DEFAULT_SHAPE "ellipse"

//!  The numerical place-holder for a NULL expression; value does not matter
#define NULL_EXPR 0

//!  The distance method used
enum DIST_METHOD {
  /*! Euclidean distance */ DIST_EUC,
  /*! Manhattan distance */ DIST_MAN,
  /*! Pearson correlation coefficient */ DIST_PEAR,
  /*! Spearman rank correlation coefficient */ DIST_SPEAR
};

//!  The linkage method used
enum LINK_METHOD {
  /*! Single linkage */ LINK_SINGLE,
  /*! Average linkage */ LINK_AVERAGE,
  /*! Complete linkage */ LINK_COMPLETE,
  /*! Centroid linkage */ LINK_CENTROID
};

//!  The scoring method used
enum SCORE_METHOD {
  /*! Score based on gap size */ SCORE_GAPS,
  /*! Score based on ANOVA */ SCORE_ANOVA,
  /*! Score based on the normalized association of Shi and Malik (2000) (modified) */ SCORE_NASSOC,
  /*! Score based on the normalized association of Shi and Malik (2000) (original) */ SCORE_NASSOC_ORIG,
};

#endif

