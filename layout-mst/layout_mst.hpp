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
    \file layout_mst.hpp
    Header file for LAYOUTMST class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: layout_mst.hpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#ifndef LAYOUT_MST_HPP
#define LAYOUT_MST_HPP

/*!
     The LAYOUTMST class is the main class for this program.  The main driver
     (in main.cpp) creates an instance of this class as the first task.  It is
     also the last class destroyed before the program exits.
*/
class LAYOUTMST {
  public:
    LAYOUTMST ();

    //  Functions related to parameters  [parameters.cpp]
    bool processOptions (int argc, char *argv[]);
    void initSettings ();
#if HAVE_MPI
    void initSettings (environment *arg1, communicator *arg2);
#endif
    bool checkSettings ();

    //  Run the main part of the program  [run.cpp]
    void sendOKFail (bool arg);
    bool recvOKFail ();
    bool runPrimary ();
    void runAllProcessors ();

    //  Read in and process scores  [process_scores.cpp]
    bool readScores ();
    void processScores ();

    //  MPI-related functions  [transmit.cpp]
    void expandPacket (PACKET *p);
    void broadcastData ();
    void receiveData ();

    //  GraphViz-related functions  [graphviz.cpp]
    void callGraphViz (unsigned int id, enum FILETYPE in_ftype, enum FILETYPE out_ftype);
    void process (unsigned int id);
    void printGraphViz (unsigned int id, double width, double height, string extension, vector<VERTEX> vertices, vector<EDGE> edges);

    //  IO-related functions  [io.cpp]
    void readNodes (unsigned int id);
    double readEdges (unsigned int id);
    void updateNodePositions (unsigned int id);

    //  Accessors  [buildmst.cpp]
    void setDebug (bool arg);
    bool getDebug () const;
    void setVerbose (bool arg);
    bool getVerbose () const;
    void setPreview (bool arg);
    bool getPreview () const;
    void setFixedPos (bool arg);
    bool getFixedPos () const;

    void setPath (string arg);
    string getPath () const;
    void setURL (string arg);
    string getURL () const;

    void setWidth (double arg);
    double getWidth () const;
    void setHeight (double arg);
    double getHeight () const;
    void setPWidth (double arg);
    double getPWidth () const;
    void setPHeight (double arg);
    double getPHeight () const;

    void setOuttype (enum FILETYPE);
    enum FILETYPE getOuttype () const;

    void setPercent (unsigned int arg);
    unsigned int getPercent () const;

    void setDPI (unsigned int arg);
    unsigned int getDPI () const;
    void setSpline (bool arg);
    bool getSpline () const;
    void setFontsize (unsigned int arg);
    unsigned int getFontsize () const;

    void setTotalIter (unsigned int arg);
    unsigned int getTotalIter () const;

    void setScoresFn (string arg);
    string getScoresFn () const;

    void setMyWorkunit (vector<unsigned int> arg);
    vector<unsigned int> getMyWorkunit () const;
    void setAllWorkunits (unsigned int arg1, vector<unsigned int> arg2);
    vector<unsigned int> getAllWorkunits (unsigned int arg) const;

#if HAVE_MPI
    void setEnv (environment *arg);
    environment *getEnv () const;
    void setComm (communicator *arg);
    communicator *getComm () const;
#endif

    void setRank (unsigned int arg);
    unsigned int getRank () const;
    void setWorldSize (unsigned int arg);
    unsigned int getWorldSize () const;
  private:
    //////////////////////////////////////////////////
    //  Values that are sent to child processes
    //!  Set to true if debug output is required; false otherwise
    bool debug_flag;
    //!  Set to true if verbose output is required; false otherwise
    bool verbose_flag;
    //!  Set to true if smaller, preview images are required; false otherwise
    bool preview_flag;
    //!  Set to true if each MSTs initial node position is fixed using the previous MST (as a result, if MPI is in use, this value is forced to false)
    bool fixed_pos;

    //!  Input and output path
    string path;
    //!  URL for image maps
    string url;
    //!  Width of the images
    double width;
    //!  Height of the images
    double height;
    //!  Width of the preview images
    double pwidth;
    //!  Height of the preview images
    double pheight;
    //!  Type of output file
    enum FILETYPE outtype;

    //!  Percentage of images to generate
    unsigned int percent;

    //!  The resolution (dots per inch) to use in the MST
    unsigned int dpi;

    //!  Set to true if lines should not cross in the MST (takes more time)
    bool spline_flag;

    //!  Size of the fonts to use
    unsigned int fontsize;

    //!  Record the vertex attributes of the current workunit
    vector<VERTEX> vertices;

    //!  Record the edges of the current workunit
    vector<EDGE> edges;

    //////////////////////////////////////////////////
    //  Values that are kept on the primary server
    //!  Total number of iterations
    unsigned int total_iter;

    //!  Scores filename
    string scores_fn;

    //!  The vector of scores
    vector<SCORE> scores;

    //!  All of the workunits as stored by the main processor only (not used by child processors)
    vector< vector<unsigned int> > all_workunits;

    //!  The workunit for processing by the current processor
    vector< unsigned int > my_workunit;

    //////////////////////////////////////////////////
    //  Open MPI values
#if HAVE_MPI
    //!  The MPI environment; NULL if MPI not available
    environment *env;

    //!  The MPI communicator; NULL if MPI not available
    communicator *comm;
#endif

    //!  Rank of this process
    unsigned int rank;

    //!  Total number of processes
    unsigned int world_size;
};

#endif


