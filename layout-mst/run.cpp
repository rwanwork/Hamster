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
    \file run.cpp
    Member functions for LAYOUTMST class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: run.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>  //  cerr
#include <string>
#include <vector>

using namespace std;

#include "LayoutMSTConfig.hpp"

#if HAVE_MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
using boost::mpi::environment;
using boost::mpi::communicator;
#endif

#include "global_defn.hpp"
#include "packet.hpp"
#include "score.hpp"
#include "vertex.hpp"
#include "edge.hpp"
#include "layout_mst.hpp"

//!  Send the initial signal to all processes to tell them to continue or not
/*!
     \param arg Indicate whether or not sub-processes should proceed (true - yes; false - no)
*/
void LAYOUTMST::sendOKFail (bool arg) {

#if HAVE_MPI
  for (unsigned int i = 1; i < getWorldSize (); i++) {
    getComm () -> send (i, 0, arg);
  }
#endif

  return;
}


//!  Receive the initial signal from the primary process to determine whether to continue or not
bool LAYOUTMST::recvOKFail () {
  bool status = false;

#if HAVE_MPI
  getComm () -> recv (0, 0, status);
#endif

  return (status);
}


//!  Execute the part of the program that is for the primary processor only
bool LAYOUTMST::runPrimary () {
  //  Read the data in; return if an error occurs
  if (!readScores ()) {
    sendOKFail (false);
    return false;
  }

  processScores ();

  //  Tell all processes that everything is fine
  sendOKFail (true);

  broadcastData ();

  return true;
}


//!  Execute the part of the program that is for all processors (including the primary)
void LAYOUTMST::runAllProcessors () {
  unsigned int i = 0;

  receiveData ();

  if (getDebug ()) {
    cerr << "==\tWorkunit size:  " << getMyWorkunit ().size () << endl;
  }

  for (i = 0; i < getMyWorkunit ().size (); i++) {
    unsigned int id = getMyWorkunit ()[i];
    process (id);
    if (HAVE_GRAPHVIZ) {
      if (getFixedPos ()) {
        callGraphViz (id, FILETYPE_DOT, FILETYPE_GV);
      }
      callGraphViz (id, FILETYPE_GV, getOuttype ());
    }
  }

  return;
}

