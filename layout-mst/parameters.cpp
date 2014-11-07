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
    Member functions for LAYOUTMST class definition
    
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

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;
using namespace boost;

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
bool LAYOUTMST::processOptions (int argc, char *argv[]) {
  //  Initialize values
  setDebug (false);
  setVerbose (false);
  setFixedPos (false);
  setSpline (false);
  setPreview (false);
  setOuttype (FILETYPE_PNG);

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
      ("fixedpos", "Attempt to fix the position of nodes")
      ("path", po::value<string>() -> default_value ("./"), "Input/output path")
      ("url", po::value<string>() -> default_value (DEFAULT_URL), "Base URL for image maps")
      ("filetype", po::value<string>(), "Output file type [ png* | svg | ps | map | dot ]")
      ("width", po::value<double>() -> default_value (DEFAULT_WIDTH), "Image width (cm)")
      ("height", po::value<double>() -> default_value (DEFAULT_HEIGHT), "Image height (cm)")
      ("percent", po::value<unsigned int>() -> default_value (DEFAULT_PERCENT), "Percentage of MSTs to generate (int)")
      ("dpi", po::value<unsigned int>() -> default_value (DEFAULT_DPI), "Resolution (dots per inch)")
      ("spline", "Prevent overlap of edges in MSTs")
      ("fontsize", po::value<unsigned int>() -> default_value (DEFAULT_FONTSIZE), "Fontsize")
      ("preview", "Show preview images")
      ("pwidth", po::value<double>() -> default_value (DEFAULT_PWIDTH), "Preview image width (cm)")
      ("pheight", po::value<double>() -> default_value (DEFAULT_PHEIGHT), "Preview image height (cm)")
      ;

    //  Hidden options that are allowed on both the command line and the configuration
    //  file, but will be hidden from the user
    po::options_description hidden ("Hidden options");
    hidden.add_options()
      ("scores", po::value<string>(), "Filename of scores")
      ;

    po::options_description cmdline_options;
    cmdline_options.add(program_only).add(config).add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config).add(hidden);

    po::options_description visible ("Allowed options");
    visible.add(program_only).add(config);

    po::positional_options_description p;
    p.add("scores", -1);

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
      cout << "Layout MST version " << LayoutMST_VERSION_MAJOR << "." << LayoutMST_VERSION_MINOR << ":  " << __DATE__ <<  " (" << __TIME__ << ")" << endl;
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

    if (vm.count ("fixedpos")) {
      setFixedPos (true);
    }

    if (vm.count ("path")) {
      setPath (vm["path"].as<string>());
    }

    if (vm.count ("url")) {
      setURL (vm["url"].as<string>());
    }

    if (vm.count ("scores")) {
      setScoresFn (vm["scores"].as<string>());
    }

    if (vm.count ("filetype")) {
      string filetype_tmp = vm["filetype"].as<string>();
      if (filetype_tmp == "png") {
        setOuttype (FILETYPE_PNG);
      }
      else if (filetype_tmp == "svg") {
        setOuttype (FILETYPE_SVG);
      }
      else if (filetype_tmp == "ps") {
        setOuttype (FILETYPE_PS);
      }
      else if (filetype_tmp == "map") {
        setOuttype (FILETYPE_CMAP);
      }
      else if (filetype_tmp == "dot") {
        setOuttype (FILETYPE_DOT);
      }
      else {
        cerr << "The argument to --filetype was not recognized:  " << filetype_tmp << endl;
        return false;
      }
    }

    if (vm.count ("width")) {
      setWidth (vm["width"].as<double>());
    }

    if (vm.count ("height")) {
      setHeight (vm["height"].as<double>());
    }

    if (vm.count ("percent")) {
      setPercent (vm["percent"].as<unsigned int>());
    }

    if (vm.count ("dpi")) {
      setDPI (vm["dpi"].as<unsigned int>());
    }

    if (vm.count ("spline")) {
      setSpline (true);
    }

    if (vm.count ("fontsize")) {
      setFontsize (vm["fontsize"].as<unsigned int>());
    }

    if (vm.count ("preview")) {
      setPreview (true);
    }

    if (vm.count ("pwidth")) {
      setPWidth (vm["pwidth"].as<double>());
    }

    if (vm.count ("pheight")) {
      setPHeight (vm["pheight"].as<double>());
    }
  }
  catch(std::exception& e) {
    cout << e.what() << "\n";
    return false;
  }

  return true;
}

//!  Initialize settings when there is no MPI available
/*!
     If MPI is unavailable, then the rank of this process is 0 and
     the size of the "world" is 1.
*/
void LAYOUTMST::initSettings () {
  setRank (0);
  setWorldSize (1);

  return;
}


#if HAVE_MPI
//!  Initialize settings to provide default values when MPI is available
/*!
     \param arg1 MPI environment
     \param arg2 MPI communicator

     We initialize enumerated types only here since strings have a default
     value already (the empty string).
*/
void LAYOUTMST::initSettings (environment *arg1, communicator *arg2) {
  setEnv (arg1);
  setComm (arg2);

  setRank (getComm () -> rank ());
  setWorldSize (getComm () -> size ());

  return;
}
#endif


//!  Check the settings to ensure they are valid.
/*!
     If no distance or linkage is set, then Euclidean and single linkage are
     assigned by default.  If --verbose has been set, then the options that the
     user chose are printed out to STDERR.
*/
bool LAYOUTMST::checkSettings () {
  string str;
  unsigned int len;

  if (getScoresFn ().empty ()) {
    cerr << "==\tError:  Filename for scores file required!" << endl;
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

  if (!getURL ().empty ()) {
    //  Ensure the path is terminated by a "/"
    str = getURL ();
    len = str.length () - 1;
    if (str.at (len) != '/') {
      setURL (getURL () + "/");
    }
  }

#if HAVE_MPI
  //  Not allowed if MPI is available and in use with > 1 processor, since each MST depends on the results of the previous one
  if (getWorldSize () != 1) {
    setFixedPos (false);
  }
#endif

  if (getVerbose ()) {
    cerr << left << setw (VERBOSE_WIDTH) << "==\tFix node positions:";
    if (getFixedPos ()) {
      cerr << "Yes" << endl;
    }
    else {
      cerr << "No" << endl;
    }

    cerr << left << setw (VERBOSE_WIDTH) << "==\tPath:";
    cerr << getPath () << endl;

    cerr << left << setw (VERBOSE_WIDTH) << "==\tURL:";
    cerr << getURL () << endl;

    cerr << left << setw (VERBOSE_WIDTH) << "==\tScores filename:" << getScoresFn () << endl;

    cerr << left << setw (VERBOSE_WIDTH) << "==\tOutput file type:";
    switch (getOuttype ()) {
      case FILETYPE_PNG :   cerr << "PNG" << endl;
                           break;
      case FILETYPE_SVG :   cerr << "SVG" << endl;
                           break;
      case FILETYPE_CMAP :  cerr << "Client-side map" << endl;
                           break;
      case FILETYPE_PS :    cerr << "Postscript" << endl;
                           break;
      case FILETYPE_DOT :    cerr << "DOT (Graphviz text)" << endl;
                           break;
      default :  cerr << "NOT SET!" << endl;
                 exit (EXIT_FAILURE);
    }

    cerr << left << setw (VERBOSE_WIDTH) << "==\tImage dimensions in cm (w x h):  "
         << getWidth () << " x " << getHeight () << endl;

    cerr << left << setw (VERBOSE_WIDTH) << "==\tPercentage of MSTs to generate:";
    cerr << getPercent () << endl;

    cerr << left << setw (VERBOSE_WIDTH) << "==\tResolution (dpi):";
    cerr << getDPI () << endl;

    cerr << left << setw (VERBOSE_WIDTH) << "==\tPrevent edge overlap:";
    if (getSpline ()) {
      cerr << "Yes" << endl;
    }
    else {
      cerr << "No" << endl;
    }

    cerr << left << setw (VERBOSE_WIDTH) << "==\tFontsize:";
    cerr << getFontsize () << endl;

    cerr << left << setw (VERBOSE_WIDTH) << "==\tShow preview images:";
    if (getPreview ()) {
      cerr << "Yes" << endl;
      cerr << left << setw (VERBOSE_WIDTH) << "==\tPreview image dimensions in cm (w x h):  "
           << getPWidth () << " x " << getPHeight () << endl;
    }
    else {
      cerr << "No" << endl;
    }

    cerr << left << setw (VERBOSE_WIDTH) << "==\tEnable large nodes:";
#if ENABLE_LARGE_NODE
      cerr << "Yes" << endl;
      cerr << left << setw (VERBOSE_WIDTH) << "==\t  Increase multipler:";
      cerr << LARGE_NODE_INCREASE << endl;
      cerr << left << setw (VERBOSE_WIDTH) << "==\t  Increase threshold:";
      cerr << LARGE_NODE_THRESH << endl;
#else
      cerr << "No" << endl;
#endif

    cerr << left << setw (VERBOSE_WIDTH) << "==\tMPI in use:";
#if HAVE_MPI
      cerr << "Yes" << endl;
      cerr << left << setw (VERBOSE_WIDTH) << "==\t  MPI World size:";
      cerr << getWorldSize () << endl;
#else
      cerr << "No" << endl;
#endif
    cerr << left << setw (VERBOSE_WIDTH) << "==\tGraphviz available:";
#if HAVE_GRAPHVIZ
      cerr << "Yes" << endl;
#else
      cerr << "No" << endl;
#endif
  }

  return true;
}

