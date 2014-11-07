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
    \file score.hpp
    Header file for SCORE class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: score.hpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#ifndef SCORE_HPP
#define SCORE_HPP

/*!
     A SCORE node keeps track of the intra and inter-cluster scores
     for a particular merge.  In addition to these scores, information
     about the merge (its unique integral ID and the IDs of the two
     clusters that were merged) are also kept track of
*/
class SCORE {
  public:
    SCORE ();
    SCORE (unsigned int arg1, unsigned int arg2, unsigned int arg3);

    //  Mutators
    void setID (unsigned int arg);
    void setLeft (unsigned int arg);
    void setRight (unsigned int arg);
    void setScore1 (double arg);
    void setScore2 (double arg);
    void setCombinedScore (double arg);

    //  Accessors
    unsigned int getID () const;
    unsigned int getLeft () const;
    unsigned int getRight () const;
    double getScore1 () const;
    double getScore2 () const;
    double getCombinedScore () const;

    void scoreGaps (unsigned int M, vector<CLUSTER> *clusters, double **d, bool debug);
    void scoreANOVA (unsigned int M, vector<CLUSTER> *clusters, double **d, bool debug);
    void scoreNormalizedAssoc (unsigned int M, vector<CLUSTER> *clusters, double **d, bool debug);
    void scoreNormalizedAssocOrig (unsigned int M, vector<CLUSTER> *clusters, double **d, bool debug);
  private:
    //!  The merge ID, numbered from 0
    unsigned int id;
    //!  The left cluster in the merge
    unsigned int left;
    //!  The right cluster in the merge
    unsigned int right;
    //!  The intra-cluster score (within-cluster) or the mean square for groups
    double score1;
    //!  The inter-cluster score (between-cluster) or the mean square error
    double score2;
    //!  The combined score calculated by either subtracting or dividing score1 and score2
    double combined;
};

#endif


