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
    \file heapnode.hpp
    Header file for HEAPNODE class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: heapnode.hpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#ifndef HEAPNODE_HPP
#define HEAPNODE_HPP

/*!
     HEAPNODEs represent potential clusters that are formed from two
     existing clusters.  We call these existing clusters "children".
     One child is the "left" child and the other child is the "right"
     child.  Of course, there is no meaning to this labeling since
     there is no notion of "left" or "right".

     A HEAPNODE becomes an actual cluster if (a) it appears
     at the top of the priority queue, and (b) both of its child
     clusters have not already been used to form another cluster earlier
     in the agglomeration process.
*/
class HEAPNODE {
	public:
		HEAPNODE ();
		//  left, right, and score
		HEAPNODE (unsigned int arg1, unsigned int arg2, double arg3);
		bool operator< (const HEAPNODE &arg) const;
		bool operator> (const HEAPNODE &arg) const;

		//  Accessors
		unsigned int getLeft () const;
		unsigned int getRight () const;
		double getScore () const;
	private:
    //!  The left child.
		unsigned int left;
    //!  The right child.
		unsigned int right;
    //!  The score indicating the dissimilarity between the two children.
		double score;
};

#endif
