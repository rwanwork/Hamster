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
    \file vect.hpp
    Header file for VECT class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: vect.hpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#ifndef VECT_HPP
#define VECT_HPP

/*!
     The VECT class represents a vector (or row) in a microarray
     data file.  Each instance has attributes such as a name,
     colour, and shape.  Also, they have a vector of expression
     levels and a vector of null values.  The lengths of these
     two vectors are identical and match one-to-one.  That is,
     if the value in position (column) i is NULL, then the
     expression level is a 0 (and should not be used).

     Since dissimilarity functions operate between two vectors,
     these functions are part of this class but are in their own
     file.
*/
class VECT {
  public:
    VECT ();
    VECT (unsigned int arg1, string arg2);
    VECT (vector<SPEARMAN> values);
    VECT (const VECT &src);

//     const VECT &operator= (const VECT &rhs);

    //  Mutators
    void setID (unsigned int arg);
    void setName (string arg);
    void setColour (string arg);
    void setShape (string arg);

    //  Accessors
    unsigned int getID () const;
    string getName () const;
    string getColour () const;
    string getShape () const;

    //!  Get the expression level at position i
    inline double getExpr (unsigned int i) const {
      return exprs[i];
    }

    //!  Get the NULL flag at position i
    inline bool getNull (unsigned int i) const {
      return nulls[i];
    }

    //!  Test if the expression level in position i is NULL
    inline bool isNull (unsigned int i) {
      return nulls[i];
    }

    //!  Put the expression level at position i
    inline void putExpr (unsigned int pos, double value) {
      exprs[pos] = value;
    }

    //!  Put the NULL value (true or false) at position i
    inline void putNull (unsigned int pos, bool value) {
      nulls[pos] = value;
    }

    //!  Get the size of the vector
    /*!
         Note that since exprs and nulls are the same size,
         we could have also returned the size of exprs.
    */
    inline unsigned int getN () const {
      return nulls.size ();
    }

    void resize (unsigned int arg);

    //  Dissimilarity functions  [vect_dist.cpp]
    double simEuc (VECT *other);
    double simMan (VECT *other);
    double simPear (VECT *other);
    double simSpear (VECT *other);
  private:
    //!  ID of this vector (or row)
    unsigned int id;
    //!  Name of the experiment corresponding to this vector
    string name;
    //!  Colour of the experiment corresponding to this vector
    string colour;
    //!  Shape of the experiment corresponding to this vector
    string shape;
    //!  Vector of expression levels
    vector<double> exprs;
    //!  Vector of NULL levels as a boolean flag (TRUE = NULL expression level)
    vector<bool> nulls;
};

#endif
