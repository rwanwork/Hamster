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
    \file graphviz.cpp
    Member functions for LAYOUTMST class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: graphviz.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <iostream>  //  cerr
#include <fstream>  //  ifstream
#include <iomanip>  //  setprecision
#include <string>
#include <vector>

#include <cfloat>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

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

//!  Call GraphViz by executing it as a shell command
/*!
     This function calls GraphViz (the program defined in
     config.h as GRAPHVIZ_NEATO_PATH) using "system".  If GraphViz has
     not been installed, then this function essentially does nothing.

     \param id ID of the output file we are generating
     \param in_ftype Type of the input file (only 2 possibilities:  DOT or GV)
     \param out_ftype Type of the output file
*/
void LAYOUTMST::callGraphViz (unsigned int id, enum FILETYPE in_ftype, enum FILETYPE out_ftype) {
  int result;
  string in_fn = " " + getPath () + lexical_cast<std::string>(id);
  string inpv_fn = " " + getPath () + lexical_cast<std::string>(id);
  string out_fn = " >" + getPath () + lexical_cast<std::string>(id);
  string outpv_fn = " >" + getPath () + lexical_cast<std::string>(id);
  string param;
  string cmd;

  if (in_ftype == FILETYPE_DOT) {
    in_fn = in_fn + DOT_FILE_EXTENSION;
    inpv_fn = inpv_fn + DOT_FILE_EXTENSION;
  }
  else if (in_ftype == FILETYPE_GV) {
    in_fn = in_fn + GV_FILE_EXTENSION;
    inpv_fn = inpv_fn + GVPV_FILE_EXTENSION;
  }
  else {
    cerr << "Invalid input file type:  " << in_ftype << endl;
    exit (1);
  }

  switch (out_ftype) {
    case FILETYPE_GV :
      out_fn += ".graphviz";
      outpv_fn += "-pv.graphviz";
      param = " -Tdot";
      break;
    case FILETYPE_DOT :
      out_fn += ".dot";
      outpv_fn += "-pv.dot";
      param = " -Tdot";
      break;
    case FILETYPE_CMAP :
      out_fn += ".map";
      outpv_fn += "-pv.map";
      param = " -Tcmapx";
      break;
    case FILETYPE_PS :
      out_fn += ".eps";
      outpv_fn += "-pv.eps";
      param = " -Tps";
      break;
    case FILETYPE_SVG :
      out_fn += ".svg";
      outpv_fn += "-pv.svg";
      param = " -Tsvg";
      break;
    default :
      out_fn += ".png";
      outpv_fn += "-pv.png";
      param = " -Tpng";
      break;
  }

  cmd = "GRAPHVIZ_NEATO_PATH" + in_fn + param + out_fn + "\n";
  result = system (cmd.c_str ());

  if (getPreview ()) {
    cmd = "GRAPHVIZ_NEATO_PATH" + inpv_fn + param + outpv_fn + "\n";
    result = system (cmd.c_str ());
  }

  return;
}


//!  Process a single MST
/*!
     \param id ID of the MST (0-based)

     This function processes (creates) a single MST, which is assigned
     a unique integral id.  For each MST, the following is done:

     1)  Read in the corresponding nodes file.
     2)  Read in the corresponding edges file.
     3)  Read in the previous MST's nodes and edges.
     4)  Normalize edge weights by the largest weight smaller than DBL_MAX.
     5)  Output this information in GraphViz format by calling printGraphViz ().
*/
void LAYOUTMST::process (unsigned int id) {
  unsigned int i = 0;
  string str;
  string fn;
  double max_weight = 0.0;

//   setUpdatedPositions (false);

  //  Reset the vertices and edges
  vertices.clear ();
  edges.clear ();

  readNodes (id);
  max_weight = readEdges (id);

  if ((getFixedPos ()) && (HAVE_GRAPHVIZ) && (getWorldSize () == 1)) {
    updateNodePositions (id);
  }

  //////////////////////////////////////////////////
  //  Normalize the edge weights to values between 0 and 1.
  //  A weight that is DBL_MAX (i.e., infinity, the "magic value"
  //  returned by some distance functions) is forced to 1.0.
  for (i= 0; i < edges.size (); i++) {
    double old_weight = edges[i].getWeight ();
    if (old_weight >= (DBL_MAX - DBL_EPSILON)) {
      //  The maximum edge weight is 1.0 after normalization
      edges[i].setWeight (1.0);
    }
    else {
      edges[i].setWeight (old_weight / max_weight);
    }
  }

  //////////////////////////////////////////////////
  //  Output DOT or GV file
  if (getFixedPos ()) {
    printGraphViz (id, getWidth (), getHeight (), DOT_FILE_EXTENSION, vertices, edges);
  }
  else {
    printGraphViz (id, getWidth (), getHeight (), GV_FILE_EXTENSION, vertices, edges);

    //  Output preview GraphViz file
    if (getPreview ()) {
      printGraphViz (id, getPWidth (), getPHeight (), GVPV_FILE_EXTENSION, vertices, edges);
    }
  }

  return;
}


//!  Output the graph information in GraphViz format
/*!
     The function has been generalized so that it can print the
     actual MSTs or the preview (smaller) ones depending on the
     arguments given to it.
*/
void LAYOUTMST::printGraphViz (unsigned int id, double width, double height, string extension, vector<VERTEX> vertices, vector<EDGE> edges) {
  string fn = getPath () + lexical_cast<std::string>(id) + extension;

#if LARGE_NODE_ENABLE
  unsigned int total_expts = 0;
  //  Find out the total number of experiments
  for (unsigned int i = 0; i < vertices.size (); i++) {
    total_expts += vertices[i].getComponents ();
  }
#endif

  ofstream fout (fn.c_str ());

  fout << "digraph test {" << endl;
  unsigned int origprec = fout.precision ();
  fout << "  graph [center=true, dpi=\"" << getDPI () << "\", ratio=\"fill\", overlap=false";
  fout << setprecision (1) << ", size=\""
      << fixed << static_cast<double> (width) / SCALE_FACTOR << ","
      << fixed << static_cast<double> (height) / SCALE_FACTOR << "\"";
  if (getSpline ()) {
    fout << ", splines=true";
  }
  else {
    fout << ", splines=false";
  }
  fout << "];" << endl;

  fout << setprecision (origprec);
  fout << "  node [label=\"\\N\", style=\"filled\"];" <<  endl;
  fout << "  edge [dir=none];" << endl;

  //  Print the nodes out
  for (unsigned int i = 0; i < vertices.size (); i++) {
    fout << "    \"" << vertices[i].getName () << "\" [label=\"" << vertices[i].getName () << "\"";
    fout << ", URL=\"" << getURL () << id << ".html?type=node&id=" << vertices[i].getName () << "\"";
    fout << ", color=\"" << vertices[i].getColour () << "\"";
    fout << ", shape=\"" << vertices[i].getShape () << "\"";
    fout << ", fontsize=" << getFontsize () << ", target=_top";
    if ((vertices[i].getUpdated ()) && (getFixedPos ())) {
      fout << ", pos=\"" << vertices[i].getX () << "," << vertices[i].getY () << "\"";
#if LARGE_NODE_ENABLE
      if (static_cast<double>(vertices[i].getComponents ()) / static_cast<double>(total_expts) > LARGE_NODE_THRESH) {
        vertices[i].setWidth (1.25 * vertices[i].getWidth ());
        vertices[i].setHeight (1.25 * vertices[i].getHeight ());
      }
#endif
      fout << ", width=\"" << vertices[i].getWidth () << "\"";
      fout << ", height=\"" << vertices[i].getHeight () << "\"";
    }
    fout << "];" << endl;
  }

  //  Print the edges out
  for (unsigned int i = 0; i < edges.size (); i++) {
    fout << "    \"" << edges[i].getStart ()
         << "\" -> \"" << edges[i].getEnd ()
         << "\" [URL=\"" << getURL () << id
         << ".html?type=edge&from=" << edges[i].getStart ()
         << "&to=" << edges[i].getEnd ()
         << "\", length=\"" << edges[i].getWeight ()
         << "\", target=_top];" << endl;
  }
  fout << "}" << endl;
  fout.close ();

  return;
}

