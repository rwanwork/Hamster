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
    \file edge.hpp
    Header file for EDGE class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: edge.hpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#ifndef EDGE_HPP
#define EDGE_HPP

/*!
     An EDGE contains three attributes:  start, end, and weight.  As the
     graph is undirected, there assignment to start and end is arbitrary.
     The edge weight is a floating point value.
*/
class EDGE {
  public:
    EDGE ();
    EDGE (string arg1, string arg2, double arg3);

    void setStart (string arg);
    string getStart () const;
    void setEnd (string arg);
    string getEnd () const;
    void setWeight (double arg);
    double getWeight () const;
  private:
    //!  Edge start
    string start;

    //!  Edge end
    string end;

    //!  Edge weight
    double weight;
};

#endif


