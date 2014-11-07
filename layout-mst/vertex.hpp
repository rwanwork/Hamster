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
    \file vertex.hpp
    Header file for VERTEX class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: vertex.hpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#ifndef VERTEX_HPP
#define VERTEX_HPP

/*!
     A VERTEX contains three attributes:  name, colour, and shape.
*/
class VERTEX {
  public:
    VERTEX ();
    VERTEX (string arg1, string arg2, string arg3, unsigned int arg4);

    void setName (string arg);
    string getName () const;
    void setColour (string arg);
    string getColour () const;
    void setShape (string arg);
    string getShape () const;
    void setComponents (unsigned int arg);
    unsigned int getComponents () const;
    void setUpdated (bool arg);
    bool getUpdated () const;

    void setX (unsigned int arg);
    unsigned int getX () const;
    void setY (unsigned int arg);
    unsigned int getY () const;
    void setWidth (float arg);
    float getWidth () const;
    void setHeight (float arg);
    float getHeight () const;
  private:
    //!  Vertex name
    string name;

    //!  Vertex colour
    string colour;

    //!  Vertex shape
    string shape;

    //!  Updated with node positions?
    bool updated;

    //!  Number of components
    unsigned int components;

    //!  x position
    unsigned int x;

    //!  y position
    unsigned int y;

    //!  Node width
    float width;

    //!  Node height
    float height;
};

#endif


