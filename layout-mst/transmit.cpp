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
    \file transmit.cpp
    Member functions for LAYOUTMST class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: transmit.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <string>
#include <vector>

#include <boost/serialization/vector.hpp>

using namespace std;
using namespace boost;

#if HAVE_MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
using boost::mpi::environment;
using boost::mpi::communicator;
#endif

#include "global_defn.hpp"
#include "score.hpp"
#include "packet.hpp"
#include "vertex.hpp"
#include "edge.hpp"
#include "layout_mst.hpp"

//!  Expand the information in a PACKET object and set the private variables
void LAYOUTMST::expandPacket (PACKET *p) {
  setDebug (p -> getDebug ());
  setVerbose (p -> getVerbose ());
  setPreview (p -> getPreview ());

  setPath (p -> getPath ());
  setURL (p -> getURL ());
  setWidth (p -> getWidth ());
  setHeight (p -> getHeight ());
  setPWidth (p -> getPWidth ());
  setPHeight (p -> getPHeight ());
  setOuttype (p -> getOuttype ());

  setPercent (p -> getPercent ());
  setDPI (p -> getDPI ());
  setSpline (p -> getSpline ());
  setFontsize (p -> getFontsize ());

  setMyWorkunit (p -> getMyWorkunit ());
}


//!  Construct PACKET objects for each processor and send it to them
/*!
     PACKET objects are delivered to each process; in case MPI is
     unavailable, then all of the workunits are simply copied from
     one vector (all_workunits[0]) to another (my_workunit)
*/
void LAYOUTMST::broadcastData () {
#if HAVE_MPI
  PACKET p;

  p.setDebug (getDebug ());
  p.setVerbose (getVerbose ());
  p.setPreview (getPreview ());
  p.setPath (getPath ());
  p.setURL (getURL ());
  p.setWidth (getWidth ());
  p.setHeight (getHeight ());
  p.setPWidth (getPWidth ());
  p.setPHeight (getPHeight ());
  p.setOuttype (getOuttype ());
  p.setPercent (getPercent ());
  p.setDPI (getDPI ());
  p.setSpline (getSpline ());
  p.setFontsize (getFontsize ());

  //  Send to each process
  for (unsigned int i = 0; i < getWorldSize (); i++) {
    if (getDebug ()) {
      cerr << "==\tBroadcasting data to " << i << " of " << getWorldSize () << " whose size is " << getAllWorkunits (i).size () << endl;
    }

    //  Copy the processor's workunits to the packet
    p.setMyWorkunit (getAllWorkunits (i));
    getComm () -> send (i, 1, p);
  }
#else
    my_workunit = all_workunits[0];
#endif

  return;
}


//!  Receive PACKET objects from the primary processor
/*!
     A PACKET object is received from the primary processor; in case MPI is
     unavailable, then nothing is done.  Instead, the copying was
     accomplished by broadcastData () already.
*/
void LAYOUTMST::receiveData () {
#if HAVE_MPI
  PACKET p;

  //  Receive from the main process (0)
  getComm () -> recv (0, 1, p);
  expandPacket (&p);

  if (getDebug ()) {
    cerr << "==\tI am process " << getRank() << " of " << getWorldSize () << " and I have " << getMyWorkunit ().size ()
         << " jobs." << endl;
  }
#endif

  return;
}

