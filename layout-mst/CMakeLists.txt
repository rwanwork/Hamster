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
##  $Id$
############################################################
CMAKE_MINIMUM_REQUIRED (VERSION 2.8)


########################################
##  Include externally defined functions

##  Set the main path and where additional CMake modules are
IF (NOT DEFINED MAIN_SRC_PATH)
  SET (MAIN_SRC_PATH ${CMAKE_CURRENT_SOURCE_DIR}/..)
  SET (CMAKE_MODULE_PATH ${MAIN_SRC_PATH}/cmake)
ENDIF (NOT DEFINED MAIN_SRC_PATH)

INCLUDE (global)
INCLUDE (add_subdirectory_once)
INCLUDE (FindGraphviz)

##  Set the top-level path if it isn't set already
IF (NOT DEFINED TOPLEVEL_PATH)
  SET (TOPLEVEL_PATH ${CMAKE_CURRENT_BINARY_DIR})
ENDIF (NOT DEFINED TOPLEVEL_PATH)

########################################
##  Set up the software

##  Project name, enabled C also because of bzip2 which seems to need it
PROJECT (layout-mst CXX C)

##  Software version number
SET (LayoutMST_VERSION_MAJOR ${LAYOUTMST_VERSION_MAJOR})
SET (LayoutMST_VERSION_MINOR ${LAYOUTMST_VERSION_MINOR})

##  List of files, excluding the main driver
SET (LAYOUTMST_SRCFILES 
  check.cpp
  edge.cpp
  packet.cpp
  io.cpp
  layout_mst.cpp
  main.cpp
  process_scores.cpp
  score.cpp
  transmit.cpp
  vertex.cpp
  parameters.cpp
  run.cpp
  graphviz.cpp
)


########################################
##  Detect OpenMP, MPI, and Graphviz -- must be before the creation of the configuration file

FIND_PACKAGE (MPI)
IF (MPI_FOUND)
  SET (HAVE_MPI 1)
ENDIF (MPI_FOUND)

FIND_PACKAGE (OpenMP)
IF (OPENMP_FOUND)
  SET (HAVE_OPENMP 1)
ENDIF (OPENMP_FOUND)

FIND_PACKAGE (Graphviz)
IF (GRAPHVIZ_FOUND)
  SET (HAVE_GRAPHVIZ 1)
  SET (GRAPHVIZ_NEATO_PATH ${GRAPHVIZ_NEATO_EXECUTABLE})
ENDIF (GRAPHVIZ_FOUND)


########################################
##  Create configuration file

##  Configure a header file to pass some of the CMake settings
##  to the source code
CONFIGURE_FILE (
  "${PROJECT_SOURCE_DIR}/LayoutMSTConfig.hpp.in"
  "${TOPLEVEL_PATH}/LayoutMSTConfig.hpp"
  )

##  Add the binary tree to the search path for include files
INCLUDE_DIRECTORIES ("${TOPLEVEL_PATH}")


########################################
##  Create the targets

##  Create the test executable
ADD_EXECUTABLE (layout-mst ${LAYOUTMST_SRCFILES})
INSTALL (TARGETS layout-mst DESTINATION bin)


########################################
##  Set initial compilation flags

##  Set compiler flags based on global variable
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")


########################################
##  Set up for Boost

##  Set the module path to get the latest version of it
SET (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake)

##  Set the root path for Boost from the environment variable
# SET (Boost_DEBUG OFF)
SET (BOOST_ROOT $ENV{BOOST_ROOT})
SET (Boost_NO_SYSTEM_PATHS true)
SET (Boost_USE_STATIC_LIBS   ON)
SET (Boost_USE_MULTITHREADED ON)
SET (Boost_USE_STATIC_RUNTIME OFF)
FIND_PACKAGE (Boost 1.42.0 REQUIRED COMPONENTS regex graph mpi program_options serialization)

IF (Boost_FOUND)
  IF (TARGET layout-mst)
    LINK_DIRECTORIES (${Boost_LIBRARY_DIRS})
    INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES (layout-mst ${Boost_LIBRARIES})
  ENDIF (TARGET layout-mst)
ENDIF ()


########################################
##  Set various values based on libraries found

IF (MPI_FOUND)
  INCLUDE_DIRECTORIES (${MPI_INCLUDE_PATH})
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
  SET (LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
  SET (CMAKE_CXX_COMPILER "${MPI_COMPILER}")
ENDIF (MPI_FOUND)

IF (OPENMP_FOUND)
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF (OPENMP_FOUND)


########################################
##  Documentation

##  Add a target to generate Doxygen documentation
FIND_PACKAGE (Doxygen)
IF (DOXYGEN_FOUND)
  CONFIGURE_FILE (${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  IF (NOT TARGET layout-mst-doc)
    ADD_CUSTOM_TARGET (layout-mst-doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating API documentation with Doxygen" VERBATIM)
    INSTALL (DIRECTORY ${MAIN_BINARY_PATH}/doc/layout-mst DESTINATION doc OPTIONAL)
  ENDIF (NOT TARGET layout-mst-doc)
ENDIF (DOXYGEN_FOUND)


############################################################
##  CTest
############################################################

# ENABLE_TESTING ()

