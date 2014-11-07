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
    \file cluster.hpp
    Header file for CLUSTER class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: cluster.hpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#ifndef CLUSTER_HPP
#define CLUSTER_HPP

/*!
     The CLUSTER class represents a set of one or more experiments.
     Like experiments, it has attributes like a name, colour, and shape.
     These attributes are inherited from the experiments that form the
     cluster.  Colours and shapes are their default values if the
     experiments that make it up are a mix of different attributes.  The
     name is simply a string representation of their integral IDs.

     Clusters also have a vector of the experiments that make it up.  We
     call them "items".  The vector of items has no order.

     If a cluster has an ancestor (ancestor = TRUE), then that means it is
     contained within a larger cluster.  When this program completes, every
     cluster has an ancestor except for the very last one.

     Clusters are formed using constructors from either a single experiment
     or from two clusters.

     Since linkages are calculated between clusters, their functions are
     part of this class.
*/
class CLUSTER {
  public:
    CLUSTER ();
    CLUSTER (unsigned int arg_id, string arg1, string arg2, string arg3, VECT arg4);
    CLUSTER (unsigned int arg_id, CLUSTER *arg1, CLUSTER *arg2, LINK_METHOD arg3, unsigned int arg4);

    //  Accessors
    void setName (string arg);
    void setColour (string arg);
    void setShape (string arg);
    unsigned int getID () const;
    string getName () const;
    string getColour () const;
    string getShape () const;
    vector<unsigned int> getItems () const;

    //  Other functions
    inline bool haveAncestors () {
      return ancestors;
    }
    void setAncestors ();
    void formCentroid (CLUSTER *a, CLUSTER *b);
    VECT getCentroid () const;

    //  Linkage functions  [cluster_link.cpp]
    double linkSingle (CLUSTER *other, double **d);
    double linkAverage (CLUSTER *other, double **d);
    double linkComplete (CLUSTER *other, double **d);
    double linkCentroid (CLUSTER *other, vector<VECT> *data, enum DIST_METHOD distance);
  private:
    //!  Numerical ID for this cluster.
    unsigned int id;
    //!  The name of this cluster; essentially the string representation of the id.
    string name;
    //!  The colour of this cluster when drawn.
    string colour;
    //!  The shape of this cluster when drawn.
    string shape;
    //!  The experiments (components) that make up this cluster.
    vector<unsigned int> items;
    //!  Does this cluster have an ancestor?
    /*!  i.e., if set to TRUE, then a larger cluster has been formed that includes
         this cluster; default value FALSE for all clusters.  */
    bool ancestors;
    //!  The vector that represents the centroid of this cluster
    VECT centroid;
};

#endif

