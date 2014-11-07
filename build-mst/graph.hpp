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
    \file graph.hpp
    Header file for GRAPH class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: graph.hpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#ifndef GRAPH_HPP
#define GRAPH_HPP


//!  An undirected graph represented as an adjacency list
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property < boost::edge_weight_t, double > > adjGraph;

//!  A vertex in the adjacency list graph
typedef boost::graph_traits< adjGraph >::vertex_descriptor Vertex;

//!  An edge in the adjacency list graph
typedef boost::graph_traits< adjGraph >::edge_descriptor Edge;

//!  A pair of edges represented as a pair of integers
typedef pair<int, int> EdgePair;

//!  Edge weight for undirected graphs
//   typedef typename UndirectedGraph::edge_property_type Weight;
//   typedef typename property_map < UndirectedGraph, edge_weight_t >::type weight = get(edge_weight, digraph);

/*!
     A GRAPH object includes an undirected graph represented as an adjacency list
     and the corresponding MST representation.

     Edges and edge weights are kept in two separate vectors in a one-to-one
     relationship such that the weight of edge i is in the weights vector at
     position i.

     In addition to creating the original graph and its corresponding MST,
     there are functions for printing the graph and its nodes to separate
     files.
*/
class GRAPH {
  public:
    GRAPH (unsigned int M, priority_queue<HEAPNODE, std::vector<HEAPNODE>, greater<HEAPNODE> > pqueue, vector<CLUSTER> clusters);
    void printEdges (unsigned int id, vector<CLUSTER> clusters, string outpath);
    void printNodes (unsigned int id, vector<CLUSTER> clusters, string outpath);
  private:
    //!  The undirected graph represented as an adjacency list
    adjGraph *g;
    //!  The set of edges
    vector<EdgePair> edges;
    //!  The set of edge weights
    vector<double> weights;
    //!  A property map holding the edge weights for the graph
    boost::property_map < adjGraph, boost::edge_weight_t >::type weights_pmap;
    //!  The MST as a vector of edges
    vector < Edge > spanning_tree_edges;
    //!  The MST as a vector of vertices [unused]
    vector < Vertex > spanning_tree_vertices;
};

#endif


