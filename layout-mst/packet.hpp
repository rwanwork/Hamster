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
    \file packet.hpp
    Header file for PACKET class definition
    
    \author Raymond Wan (r-wan@cb.k.u-tokyo.ac.jp)
    \par Organizations
          - Department of Computational Biology, Graduate School of
            Frontier Science, University of Tokyo
          - Computational Biology Research Center, AIST, Japan
          
    $Id: packet.hpp 3 2011-08-25 10:19:50Z rwan $

*/
/*******************************************************************/


#ifndef PACKET_HPP
#define PACKET_HPP

/*!
     A PACKET represents the object with all of the variables that each
     child process needs from the main, primary process to do its work.
     They have been collected together into a class to facilitate the
     Boost Serialization class.

     Many of the variables are duplicated from the LAYOUTMST class.
*/
class PACKET {
  public:
    PACKET ();

    void setDebug (bool arg);
    bool getDebug () const;
    void setVerbose (bool arg);
    bool getVerbose () const;
    void setPreview (bool arg);
    bool getPreview () const;
    void setFixedPos (bool arg);
    bool getFixedPos () const;

    void setPath (string arg);
    string getPath () const;
    void setURL (string arg);
    string getURL () const;

    void setWidth (double arg);
    double getWidth () const;
    void setHeight (double arg);
    double getHeight () const;
    void setPWidth (double arg);
    double getPWidth () const;
    void setPHeight (double arg);
    double getPHeight () const;

    void setOuttype (enum FILETYPE);
    enum FILETYPE getOuttype () const;

    void setPercent (unsigned int arg);
    unsigned int getPercent () const;

    void setDPI (unsigned int arg);
    unsigned int getDPI () const;
    void setSpline (bool arg);
    bool getSpline () const;
    void setFontsize (unsigned int arg);
    unsigned int getFontsize () const;

    void setMyWorkunit (vector<unsigned int> arg);
    vector<unsigned int> getMyWorkunit () const;
  private:
    //////////////////////////////////////////////////
    //  Values that are sent to child processes
    //!  Set to true if debug output is required; false otherwise
    bool debug_flag;
    //!  Set to true if verbose output is required; false otherwise
    bool verbose_flag;
    //!  Set to true if smaller, preview images are required; false otherwise
    bool preview_flag;
    //!  Set to true if each MSTs initial node position is fixed using the previous MST (as a result, if MPI is in use, this value is forced to false)
    bool fixed_pos;

    //!  Path to files
    string path;
    //!  URL for image maps
    string url;
    //!  Width of the images
    double width;
    //!  Height of the images
    double height;
    //!  Width of the preview images
    double pwidth;
    //!  Height of the preview images
    double pheight;
    //!  Type of output file
    enum FILETYPE outtype;

    //!  Percent of images to generate
    unsigned int percent;

    //!  The resolution (dots per inch) to use in the MST
    unsigned int dpi;

    //!  Set to true if lines should not cross in the MST (takes more time)
    bool spline_flag;

    //!  Size of the fonts to use
    unsigned int fontsize;

    //!  The workunit that the current process has to worry about
    vector< unsigned int > my_workunit;
#if HAVE_MPI
    friend class boost::serialization::access;
    template<class Archive>
    void serialize (Archive &ar, const unsigned int version) {
      ar & debug_flag;
      ar & verbose_flag;
      ar & preview_flag;
      ar & fixed_pos;
      ar & path;
      ar & url;
      ar & width;
      ar & height;
      ar & pwidth;
      ar & pheight;
      ar & outtype;
      ar & percent;
      ar & dpi;
      ar & spline_flag;
      ar & fontsize;
      ar & my_workunit;
    }
#endif
};

#endif


