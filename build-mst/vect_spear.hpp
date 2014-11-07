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
    \file vect_spear.hpp
    Header file for SPEARMAN class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: vect_spear.hpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#ifndef VECT_SPEAR_HPP
#define VECT_SPEAR_HPP

/*!
     A SPEARMAN node is used to keep track of the expression value,
     original rank (position), and final, sorted rank.

     The nodes are sorted based on the value in the key.  For
     example, if we are sorting based on the expression values,
     then they are copied to the key and a sort function is applied
     to the array.
*/
class SPEARMAN {
  public:
    SPEARMAN ();
    SPEARMAN (double v, unsigned int p, double r);

    //  Mutators
    void setValue (double v);
    void setOrigPos (unsigned int p);
    void setRank (double r);

    //  Accessors
    double  getValue () const;
    unsigned int getOrigPos () const;
    double getRank () const;
    bool operator< (const SPEARMAN &arg) const;
    void copyOrigPosToKey ();
  private:
    //!  The sort key (either the expression value or the original position in the vector)
    double key;
    //!  The expression level of the node
    double value;
    //!  The original position in the vector
    unsigned int origpos;
    //!  The rank of node within the vector
    double rank;
};

#endif

