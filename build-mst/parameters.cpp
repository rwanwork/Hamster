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
    \file parameters.cpp
    Member functions for BUILDMST class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: parameters.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>
#include <iomanip>  //  setw
#include <fstream>
#include <string>
#include <vector>
#include <queue>  //  priority_queue

#include <boost/program_options.hpp>

using namespace std;
using namespace boost;
namespace po = boost::program_options;

#include "BuildMSTConfig.hpp"
#include "global_defn.hpp"
#include "heapnode.hpp"
#include "vect_spear.hpp"
#include "vect.hpp"
#include "cluster.hpp"
#include "score.hpp"
#include "build_mst.hpp"


//!  Initialize settings to provide default values
/*!
     We initialize values which are not initialized by the processOptions function.
*/
void BUILDMST::initSettings () {
  setM (0);
  setN (0);

  return;
}


//!  Process options from the command line and the configuration file CFG_FILENAME
/*!
     This function makes use of Boost's program_options for handling
     arguments on the command line and in options in a configuration
     file whose format resembles .ini files.

     Initially, boolean and enumerated values are given default values.
     Then, the available options are set up, with default values for
     string and numeric types.  The description of the options are
     recorded.

     Next, the command line options are read, followed by the configuration
     file options.  The command line options take priority over the
     configuration file ones.  Then, the options are processed, one-by-one.

     All of this is encapsulated within a try...catch block.
*/
bool BUILDMST::processOptions (int argc, char *argv[]) {
  //  Initialize default values
  setDebug (false);
  setVerbose (false);
  setDistance (DIST_EUC);
  setLinkage (LINK_SINGLE);
  setScoreMethod (SCORE_GAPS);
  setCentroid (DIST_EUC);

  //  Provide a help if no arguments provided
  if (argc == 1) {
    cout << "Execute " << argv[0] << " with the --help parameter for a list of options." << endl;
    return false;
  }

  try {
    //  Options that are allowed only on the command line
    po::options_description program_only ("Program options");
    program_only.add_options()
      ("version,v", "Print version")
      ("help,h", "This help message")
      ;

    //  Options that are allowed on the command line and in the configuration file
    po::options_description config ("Configuration");
    config.add_options()
      ("debug", "Turn debugging on")
      ("verbose", "Turn verbose output on")
      ("path", po::value<string>() -> default_value ("./"), "Input/output path")
      ("distance", po::value<string>(), "Distance method [ euclidean* | manhattan | pearson | spearman ]")
      ("linkage", po::value<string>(), "Linkage method [ single* | average | complete | centroid ]")
      ("scoring", po::value<string>(), "Scoring method [ gaps* | anova | nassoc | nassoc_orig ]")
      ("centroid", po::value<string>(), "Centroid distance method [ euclidean* | manhattan | pearson | spearman ]")
      ("attr", po::value<string>(), "Attribute filename")
      ;

    //  Hidden options that are allowed on both the command line and the configuration
    //  file, but will be hidden from the user
    po::options_description hidden ("Hidden options");
    hidden.add_options()
      ("microarray", po::value<string>(), "Microarray filename")
      ;

    po::options_description cmdline_options;
    cmdline_options.add(program_only).add(config).add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config).add(hidden);

    po::options_description visible ("Allowed options");
    visible.add(program_only).add(config);

    po::positional_options_description p;
    p.add("microarray", -1);

    po::variables_map vm;
    store (po::command_line_parser (argc, argv).options (cmdline_options).positional(p).run (), vm);

    ifstream cfg_fp (CFG_FILENAME, ios::in);
    if (!cfg_fp) {
      cerr << "==\tWarning:  Configuration file " << CFG_FILENAME << " could not be opened!" << endl;
    }
    else {
      store (parse_config_file (cfg_fp, config_file_options), vm);
      notify (vm);
    }

    //  Handle each option
    if (vm.count ("version")) {
      cout << "Build MST version " << BuildMST_VERSION_MAJOR << "." << BuildMST_VERSION_MINOR << ":  " << __DATE__ <<  " (" << __TIME__ << ")" << endl;
      return false;
    }

    if (vm.count ("help")) {
      cout << visible << endl;
      cout << "Name of input file provided on the command line with no parameter." << endl;
      cout << "* indicates default values." << endl;
      return false;
    }

    if (vm.count ("debug")) {
      setDebug (true);
    }

    if (vm.count ("verbose")) {
      setVerbose (true);
    }

    if (vm.count ("path")) {
      setPath (vm["path"].as<string>());
    }

    if (vm.count ("distance")) {
      string distance_tmp = vm["distance"].as<string>();
      if (distance_tmp == "euclidean") {
        setDistance (DIST_EUC);
      }
      else if (distance_tmp == "manhattan") {
        setDistance (DIST_MAN);
      }
      else if (distance_tmp == "pearson") {
        setDistance (DIST_PEAR);
      }
      else if (distance_tmp == "spearman") {
        setDistance (DIST_SPEAR);
      }
      else {
        cerr << "The argument to --distance was not recognized:  " << distance_tmp << endl;
        return false;
      }
    }

    if (vm.count ("linkage")) {
      string linkage_tmp = vm["linkage"].as<string>();
      if (linkage_tmp == "single") {
        setLinkage (LINK_SINGLE);
      }
      else if (linkage_tmp == "average") {
        setLinkage (LINK_AVERAGE);
      }
      else if (linkage_tmp == "complete") {
        setLinkage (LINK_COMPLETE);
      }
      else if (linkage_tmp == "centroid") {
        setLinkage (LINK_CENTROID);
      }
      else {
        cerr << "The argument to --linkage was not recognized:  " << linkage_tmp << endl;
        return false;
      }
    }

    if (vm.count ("scoring")) {
      string scoring_tmp = vm["scoring"].as<string>();
      if (scoring_tmp == "gaps") {
        setScoreMethod (SCORE_GAPS);
      }
      else if (scoring_tmp == "anova") {
        setScoreMethod (SCORE_ANOVA);
      }
      else if (scoring_tmp == "nassoc") {
        setScoreMethod (SCORE_NASSOC);
      }
      else if (scoring_tmp == "nassoc_orig") {
        setScoreMethod (SCORE_NASSOC_ORIG);
      }
      else {
        cerr << "The argument to --scoring was not recognized:  " << scoring_tmp << endl;
        return false;
      }
    }

    if (vm.count ("centroid")) {
      string centroid_tmp = vm["centroid"].as<string>();
      if (centroid_tmp == "euclidean") {
        setCentroid (DIST_EUC);
      }
      else if (centroid_tmp == "manhattan") {
        setCentroid (DIST_MAN);
      }
      else if (centroid_tmp == "pearson") {
        setCentroid (DIST_PEAR);
      }
      else if (centroid_tmp == "spearman") {
        setCentroid (DIST_SPEAR);
      }
      else {
        cerr << "The argument to --centroid was not recognized:  " << centroid_tmp << endl;
        return false;
      }
    }

    if (vm.count ("attr")) {
      setAttrFn (vm["attr"].as<string>());
    }

    if (vm.count ("microarray")) {
      setMicroarrayFn (vm["microarray"].as<string>());
    }
  }
  catch(std::exception& e) {
    cout << e.what() << "\n";
    return false;
  }

  return true;
}


//!  Check the settings to ensure they are valid.
/*!
     If no distance or linkage is set, then Euclidean and single linkage are
     assigned by default.  If --verbose has been set, then the options that the
     user chose are printed out to STDERR.
*/
bool BUILDMST::checkSettings () {
  string str;
  unsigned int len;

  if (getMicroarrayFn ().empty ()) {
    cerr << "==\tError:  Microarray filename required!" << endl;
    return false;
  }

  if (!getPath ().empty ()) {
    //  Ensure the path is terminated by a "/"
    str = getPath ();
    len = str.length () - 1;
    if (str.at (len) != '/') {
      setPath (getPath () + "/");
    }
  }

  if (getVerbose ()) {
    cerr << left << setw (VERBOSE_WIDTH) << "==\tDistance method:";
    switch (getDistance ()) {
      case DIST_EUC    : cerr << "Euclidean distance";
        break;
      case DIST_MAN    : cerr << "Manhattan distance";
        break;
      case DIST_PEAR   : cerr << "Pearson correlation";
        break;
      case DIST_SPEAR  : cerr << "Spearman correlation";
        break;
    }
    cerr << endl;

    cerr << left << setw (VERBOSE_WIDTH) << "==\tLinkage method:";
    switch (getLinkage ()) {
      case LINK_SINGLE  :  cerr << "Single linkage";
        break;
      case LINK_AVERAGE :  cerr << "Average linkage";
        break;
      case LINK_COMPLETE:  cerr << "Complete linkage";
        break;
      case LINK_CENTROID:
        cerr << "Centroid linkage";
        cerr << endl;
        cerr << left << setw (VERBOSE_WIDTH) << "==\tCentroid distance method:";
        switch (getCentroid ()) {
          case DIST_EUC    : cerr << "Euclidean distance";
            break;
          case DIST_MAN    : cerr << "Manhattan distance";
            break;
          case DIST_PEAR   : cerr << "Pearson correlation";
            break;
          case DIST_SPEAR  : cerr << "Spearman correlation";
            break;
        }
        break;
    }
    cerr << endl;

    cerr << left << setw (VERBOSE_WIDTH) << "==\tScoring method:";
    switch (getScoreMethod ()) {
      case SCORE_GAPS         : cerr << "Gap-based";
        break;
      case SCORE_ANOVA        : cerr << "ANOVA";
        break;
      case SCORE_NASSOC       : cerr << "Normalized association (modified)";
        break;
      case SCORE_NASSOC_ORIG  : cerr << "Normalized association (original)";
        break;
    }
    cerr << endl;

    cerr << left << setw (VERBOSE_WIDTH) << "==\tMicroarray filename:" << getMicroarrayFn () << endl;
    cerr << left << setw (VERBOSE_WIDTH) << "==\tAttribute filename:";
    if (getAttrFn ().empty ()) {
      cerr << "N/A" << endl;
    }
    else {
      cerr << getAttrFn () << endl;
    }

    cerr << left << setw (VERBOSE_WIDTH) << "==\tOutput path:";
    if (getPath ().empty ()) {
      cerr << "Current directory" << endl;
    }
    else {
      cerr << getPath () << endl;
    }
  }

  return true;
}


