###########################################################################
##  HAMSTER
##  Software for depicting microarray data sets as a set of minimum spanning
##    trees.
##  
##  Version 1.3 -- August 26, 2011
##  
##  Copyright (C) 2009-2011 by Raymond Wan, All rights reserved.
##  Contact:  r-wan@cb.k.u-tokyo.ac.jp
##  Organization:  Department of Computational Biology, Graduate School of
##                 Frontier Science, University of Tokyo and
##                 Computational Biology Research Center, AIST,
##                 Japan
##  
##  This file is part of HAMSTER.
##  
##  HAMSTER is free software; you can redistribute it and/or 
##  modify it under the terms of the GNU Lesser General Public License 
##  as published by the Free Software Foundation; either version 
##  3 of the License, or (at your option) any later version.
##  
##  HAMSTER is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU Lesser General Public License for more details.
##  
##  You should have received a copy of the GNU Lesser General Public 
##  License along with HAMSTER; if not, see 
##  <http://www.gnu.org/licenses/>.
###########################################################################


############################################################
##  CMake file
##
##  Raymond Wan
##  Organizations
##    - Department of Computational Biology, Graduate School of
##      Frontier Science, University of Tokyo
##    - Computational Biology Research Center, AIST, Japan
##
##  $Id: CMakeLists.txt 3 2011-08-25 10:19:50Z rwan $
############################################################
CMAKE_MINIMUM_REQUIRED (VERSION 2.8)


##  Project name
PROJECT (hamster CXX)

########################################
##  Include externally defined functions

IF (NOT DEFINED ${MAIN_BINARY_PATH})
  SET (MAIN_BINARY_PATH ${CMAKE_CURRENT_BINARY_DIR})
ENDIF (NOT DEFINED ${MAIN_BINARY_PATH})

IF (NOT DEFINED ${MAIN_SRC_PATH})
  SET (MAIN_SRC_PATH ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF (NOT DEFINED ${MAIN_SRC_PATH})

SET (CMAKE_MODULE_PATH ${MAIN_SRC_PATH}/cmake)
INCLUDE (global)
INCLUDE (add_subdirectory_once)

##  Set the top-level path and make it as a possible include directory
SET (TOPLEVEL_PATH ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "Top-level directory")
INCLUDE_DIRECTORIES (${TOPLEVEL_PATH})

##  Create a directory for the Doxygen-generated documentation
FILE (MAKE_DIRECTORY ${MAIN_BINARY_PATH}/doc/)

########################################
##  Add dependencies and directories

INCLUDE_DIRECTORIES (build-mst)
ADD_SUBDIRECTORY_ONCE (build-mst ${CMAKE_CURRENT_BINARY_DIR}/build-mst)

INCLUDE_DIRECTORIES (layout-mst)
ADD_SUBDIRECTORY_ONCE (layout-mst ${CMAKE_CURRENT_BINARY_DIR}/layout-mst)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/CMakeLists.txt)
  ADD_SUBDIRECTORY (doc)
endif ()


########################################
##  Set initial compilation flags

##  Set compiler flags based on global variable
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")


############################################################
##  CTest
############################################################

ENABLE_TESTING ()


############################################################
##  Installation
############################################################

INSTALL (FILES doc/documentation.html DESTINATION doc)
INSTALL (FILES doc/AUTHORS DESTINATION doc)
INSTALL (FILES doc/COPYING DESTINATION doc)
INSTALL (FILES doc/COPYING.DOC DESTINATION doc)
INSTALL (FILES doc/ChangeLog DESTINATION doc)
INSTALL (FILES doc/README DESTINATION doc)
INSTALL (FILES doc/README.html DESTINATION doc)
INSTALL (FILES doc/VERSION DESTINATION doc)
INSTALL (FILES data/build-mst.cfg DESTINATION data)
INSTALL (FILES data/layout-mst.cfg DESTINATION data)
INSTALL (FILES data/sample.attr DESTINATION data)
INSTALL (FILES data/sample.data DESTINATION data)


############################################################
##  CPack
############################################################

if (EXISTS ${MAIN_SRC_PATH}/cmake/cpack.cmake)
  Include (${MAIN_SRC_PATH}/cmake/cpack.cmake)
  Include (CPack)
endif ()
