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
    \file main.cpp
    main () function
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: main.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>  //  cerr, endl
#include <string>
#include <vector>

#include <cstdlib>  //  EXIT_SUCCESS

using namespace std;

#include "global_defn.hpp"
#include "packet.hpp"
#include "score.hpp"
#include "vertex.hpp"
#include "edge.hpp"
#include "layout_mst.hpp"

//!  The main () function of the program
/*!
     Create a LAYOUTMST object and then uses it to read in the parameters from
     the file and the command line.  If all the settings check out, then run
     the main program.
*/
int main (int argc, char *argv[]) {
#if HAVE_MPI
  environment *env = new environment (argc, argv);
  communicator *comm = new communicator ();
#endif
  LAYOUTMST *hamster = new LAYOUTMST ();

#if HAVE_MPI
  hamster -> initSettings (env, comm);
#else
  hamster -> initSettings ();
#endif
  if (hamster -> getRank () == 0) {
    //  Process the options from the command line and configuration file
    if (!hamster -> processOptions (argc, argv)) {
      //  Something has gone wrong; kill all other processes
      hamster -> sendOKFail (false);
      return (EXIT_FAILURE);
    }

    //  Check if the parameters are all ok; if not, send a termination signal
    if (!hamster -> checkSettings ()) {
      //  Something has gone wrong; kill all other processes
      hamster -> sendOKFail (false);
      return (EXIT_FAILURE);
    }

    //  If there is a problem reading in the scores file, we should exit
    if (!hamster -> runPrimary ()) {
      return (EXIT_FAILURE);
    }
  }
  else {
    //  If a "kill" signal is received, then terminate
    if (hamster -> recvOKFail () == false) {
      return (EXIT_SUCCESS);
    }
  }

  //  All processes, including the main one run this
  hamster -> runAllProcessors ();

  return (EXIT_SUCCESS);
}


