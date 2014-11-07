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
    \file check.cpp
    General non-class fucntions for checking strings
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: check.cpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#include <string>

#include <cctype>  //  isalnum

using namespace std;

#include "check.hpp"

//!  Sanitize a filename by disallowing [/\]
/*!
     We allow all characters since the filename is not as dangerous.
     However, we do not allow the slash or backslash characters to
     prevent a change in path and any form of escaping.
*/
string sanitizeFilename (string arg) {
  string str;

  for (unsigned int i = 0; i < arg.length (); i++) {
    //  Note that we do not allow the backslash
    if ((arg[i] != '/') || (arg[i] != '\\')) {
      str += arg[i];
    }
  }

  return str;
}

//!  Sanitize a path by allowing alphanumeric characters and [./]
/*!
     Note that the backslash character has been purposely excluded to prevent
     any escaping.  This will cause problems to the Windows' family of
     operating systems and should be added in, if required.
*/
string sanitizePath (string arg) {
  string str;

  for (unsigned int i = 0; i < arg.length (); i++) {
    //  Note that we do not allow the backslash, so this could
    //  cause problems with Windows' family of operating systems
    if (isalnum (arg[i]) || (arg[i] == '/') || (arg[i] == '.')) {
      str += arg[i];
    }
  }

  return str;
}

//!  Sanitize a URL by allowing alphanumeric characters and [/:.]
/*!
     Note that the backslash character has been purposely excluded to prevent
     any escaping.
*/
string sanitizeURL (string arg) {
  string str;

  for (unsigned int i = 0; i < arg.length (); i++) {
    //  Note that we do not allow the backslash
    if (isalnum (arg[i]) || (arg[i] == '/') || (arg[i] == ':') || (arg[i] == '.')) {
      str += arg[i];
    }
  }

  return str;
}

