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
    \file vect.cpp
    Member functions for VECT class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: vect.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <string>
#include <vector>
#include <iostream>  //  cerr, endl

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

#include "global_defn.hpp"
#include "check.hpp"
#include "vect_spear.hpp"
#include "vect.hpp"

//!  Default constructor that takes no arguments
VECT::VECT ()
  : id (0),
    name (""),
    colour (DEFAULT_COLOUR),
    shape (DEFAULT_SHAPE),
    exprs (),
    nulls ()
{
}

//!  Constructor that takes two arguments
/*!
     \param arg1 ID of the vector
     \param arg2 The row from the microarray data file, represented as a string
*/
VECT::VECT (unsigned int arg1, string arg2)
  : id (arg1),
    name (""),
    colour (DEFAULT_COLOUR),
    shape (DEFAULT_SHAPE),
    exprs (),
    nulls ()
{
  typedef tokenizer<char_separator<char> > tokenizer;
  char_separator<char> sep ("\t");  //  Must be tab-separated
  tokenizer tokens (arg2, sep);

  //  Grab the name of the vector (first column)
  tokenizer::iterator beg = tokens.begin ();
  setName (*beg);
  beg++;

  for (; beg != tokens.end (); ++beg) {
    if (*beg == "NULL") {
      exprs.push_back (NULL_EXPR);
      nulls.push_back (true);
    }
    else {
      try {
        exprs.push_back (lexical_cast<double>(*beg));
      }
      catch (bad_lexical_cast &) {
        cerr << "\nUnexpected error in reading in microarray data.\nExpecting the keyword \"NULL\" or a double value, but found this instead:  " << *beg << "\nExiting...\n\n";
        exit (EXIT_FAILURE);
      }
      nulls.push_back (false);
    }
//       cout << *beg << endl;
  }
}

//!  Constructor that takes a single argument
/*!
     \param values A vector of SPEARMAN nodes

     This constructor copies the values in the SPEARMAN nodes into a
     VECT object.  It assumes that there are NULL values since every
     value is a rank.  NULL values have an expression level of NULL_EXPR
     and would have all been pushed to one side.
*/
VECT::VECT (vector<SPEARMAN> values)
  : id (0),
    name (""),
    colour (DEFAULT_COLOUR),
    shape (DEFAULT_SHAPE),
    exprs (),
    nulls ()
{
  unsigned int i = 0;
  unsigned int size = values.size ();

  for (i = 0; i < size; i++) {
    exprs.push_back (values[i].getRank ());
    nulls.push_back (false);
  }

  return;
}

//!  Copy constructor for a VECT object
/*!
     \param src The source that we are copying from

     This function is used to copy the two vectors an element at a time.
     The function receives a reference to prevent infinite recursion.
*/
VECT::VECT (const VECT &src)
  : id (src.id),
    name (src.name),
    colour (src.colour),
    shape (src.shape),
    exprs (),
    nulls ()
{
   const unsigned int size = src.getN ();

   //  Copy the two arrays one element at a time
   for (unsigned int i = 0; i < size; i++) {
     exprs.push_back (src.getExpr (i));
     nulls.push_back (src.getNull (i));
   }
}

//!  Set the ID
void VECT::setID (unsigned int arg) {
  id = arg;
}

//!  Get the ID
unsigned int VECT::getID () const {
  return id;
}

//!  Set the name
void VECT::setName (string arg) {
  name = sanitizeSampleName (arg);
}

//!  Get the name
string VECT::getName () const {
  return name;
}

//!  Set the colour
void VECT::setColour (string arg) {
  colour = arg;
}

//!  Get the colour
string VECT::getColour () const {
  return colour;
}

//!  Set the shape
void VECT::setShape (string arg) {
  shape = arg;
}

//!  Get the shape
string VECT::getShape () const {
  return shape;
}

//!  Resize the vector
/*!
     In order to resize the VECT object, we have to keep the lengths
     of both exprs and nulls the same at all times; so we change
     both simultaneously here.
*/
void VECT::resize (unsigned int arg) {
  exprs.resize (arg);
  nulls.resize (arg);
}

