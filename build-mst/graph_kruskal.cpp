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
    \file graph_kruskal.cpp
    Functions for GRAPH class definition related to Kruskal's algorithm
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: graph_kruskal.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <fstream>  //  ofstream
#include <vector>
#include <iostream>
#include <queue>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

#include "global_defn.hpp"
#include "vect_spear.hpp"
#include "vect.hpp"
#include "cluster.hpp"
#include "heapnode.hpp"
#include "graph.hpp"


//!  Constructor for GRAPH object with three parameters
/*!
     \param M Number of experiments (rows) for the current graph (decreases by 1 with each iteration)
     \param pqueue Priority queue of potential clusters
     \param clusters Vector of clusters

     The constructor builds the graph and then calculates
     the MST using Boost Graph Library's implementation of
     Kruskal's algorithm.
*/
GRAPH::GRAPH (unsigned int M, priority_queue<HEAPNODE, std::vector<HEAPNODE>, greater<HEAPNODE> > pqueue, vector<CLUSTER> clusters)
  : g (),
    edges (),
    weights (),
    weights_pmap (),
    spanning_tree_edges (),
    spanning_tree_vertices ()
{
  unsigned int i = 0;
  HEAPNODE heapnode;

  EdgePair *edge_array = NULL;
  double *weights_array = NULL;

  priority_queue<HEAPNODE, std::vector<HEAPNODE>, greater<HEAPNODE> > pqueue2;

  //  Number of edges and number of weights is the same
  size_t num_edges = ((M * M) - M) / 2;
  edge_array = new EdgePair[num_edges];
  weights_array = new double[num_edges];

  while (!pqueue.empty ()) {
    heapnode = pqueue.top ();
    //  If left and right don't have parents, then not merged yet
    unsigned int left = heapnode.getLeft ();
    unsigned int right = heapnode.getRight ();
    double score = heapnode.getScore ();
    if ((!clusters[left].haveAncestors ()) && (!clusters[right].haveAncestors ())) {
      edge_array[i] = EdgePair (left, right);
      weights_array[i] = score;
      i++;
      pqueue2.push (heapnode);
    }
    pqueue.pop ();
  }

  //  Update the priority queue with only pairs that have no ancestors
  pqueue = pqueue2;

  //  Create the graph using the edge iterator constructor or add_edge ();
  //  adjacency_matrix can only accept add_edge ()!!
  g = new adjGraph (edge_array, edge_array + num_edges, weights_array, M);
  weights_pmap = get (boost::edge_weight, *g);

  //  Calculate the MST using Kruskal's algorithm
  boost::kruskal_minimum_spanning_tree (*g, back_inserter (spanning_tree_edges));

  //  Deallocate memory of the arrays
  delete [] edge_array;
  delete [] weights_array;

  return;
}


//!  Print the edges of the MST out to file
void GRAPH::printEdges (unsigned int id, vector<CLUSTER> clusters, string outpath) {
  string src;
  string dest;
  string fn = outpath + lexical_cast<std::string>(id) + EDGES_FILE_EXTENSION;

  ofstream fout (fn.c_str ());
  if (!fout) {
    cerr << "File " << fn << " could not be opened!" << endl;
    exit (EXIT_FAILURE);
  }

  vector < Edge >::iterator ei;
  for (ei = spanning_tree_edges.begin(); ei != spanning_tree_edges.end(); ++ei) {
    src = clusters[boost::source(*ei, *g)].getName ();
    dest = clusters[boost::target(*ei, *g)].getName ();
    //  Ensure src < dest
    if (src > dest) {
      src.swap (dest);
    }
    fout << src << "\t"
         << dest << "\t"
         << weights_pmap[*ei] << endl;
  }
  fout.close ();

  return;
}


//!  Print the nodes of the MST out to file
void GRAPH::printNodes (unsigned int id, vector<CLUSTER> clusters, string outpath) {
  unsigned int i = 0;
  string src;
  string dest;
  string fn = outpath + lexical_cast<std::string>(id) + NODES_FILE_EXTENSION;

  ofstream fout (fn.c_str ());
  if (!fout) {
    cerr << "File " << fn << " could not be opened!" << endl;
    exit (EXIT_FAILURE);
  }

  //  Print out clusters which have no parents
  for (i = 0; i < clusters.size (); i++) {
    //  If no parents, then print it out
    if (!clusters[i].haveAncestors ()) {
      fout << clusters[i].getName () << "\t"
           << clusters[i].getColour () << "\t"
           << clusters[i].getShape () << "\t"
           << clusters[i].getItems ().size () << endl;
    }
  }

  fout.close ();

  return;
}

