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
    \file build_mst.hpp
    Header file for BUILDMST class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: build_mst.hpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#ifndef BUILD_MST_HPP
#define BUILD_MST_HPP

/*!
     The BUILDMST class is the main class for this program.  The main driver
     (in main.cpp) creates an instance of this class as the first task.  It is
     also the last class destroyed before the program exits.
*/
class BUILDMST {
  public:
    BUILDMST ();

    //  Functions related to parameters  [parameters.cpp]
    bool processOptions (int argc, char *argv[]);
    void initSettings ();
    bool checkSettings ();
    
    //  Main part of the class  [run.cpp]
    void run ();

    //  Data file I/O  [io.cpp]
    bool readMicroarray ();
    bool readAttr ();

    //  Calculate distances or clusters  [calculate.cpp]
    void initializeDistances ();
    void initializeClusters ();
    void calculateLinkage (CLUSTER &arg);
    void calculateScores (SCORE &arg);

    //  Normalize and print the scores  [calculate.cpp]
    void normalizeScores ();
    bool printScores (string outpath);

    //  Accessors  [build_mst.cpp]
    void setDebug (bool arg);
    bool getDebug () const;
    void setVerbose (bool arg);
    bool getVerbose () const;

    void setDistance (DIST_METHOD arg);
    DIST_METHOD getDistance () const;
    void setLinkage (LINK_METHOD arg);
    LINK_METHOD getLinkage () const;
    void setScoreMethod (SCORE_METHOD arg);
    SCORE_METHOD getScoreMethod () const;
    void setCentroid (DIST_METHOD arg);
    DIST_METHOD getCentroid () const;

    void setMicroarrayFn (string arg);
    string getMicroarrayFn () const;
    void setAttrFn (string arg);
    string getAttrFn () const;
    void setPath (string arg);
    string getPath () const;
    void setM (unsigned int arg);
    unsigned int getM () const;
    void setN (unsigned int arg);
    unsigned int getN () const;
  private:
    //!  Set to true if debug output is required; false otherwise
    bool debug_flag;
    //!  Set to true if verbose output is required; false otherwise
    bool verbose_flag;
    //!  Distance method
    enum DIST_METHOD distance;
    //!  Linkage method
    enum LINK_METHOD linkage;
    //!  MST scoring method
    enum SCORE_METHOD scoring;
    //!  Distance method for centroid linkage
    enum DIST_METHOD centroid;
    //!  Attribute filename
    string attr_fn;
    //!  Microarray filename
    string microarray_fn;
    //!  Output path
    string path;

    //!  Number of rows
    unsigned int M;
    //!  Number of columns
    unsigned int N;

    //!  The vector of clusters; grows from M to at most (M + M - 1) entries
    vector<CLUSTER> clusters;

    //!  The priority queue, implemented as a heap
    priority_queue<HEAPNODE, std::vector<HEAPNODE>, greater<HEAPNODE> > pqueue;

    /*!  The vector of scores; each position in this array represents
    the score of the MST from one merge step.  */
    vector<SCORE> scores;

    //!  Original microarray data with each row as an entry in this vector
    vector<VECT> data;
    //!  Distance matrix of size (M * M)
    double **dist_matrix;
};

#endif
