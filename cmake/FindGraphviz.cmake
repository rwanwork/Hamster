# - This module looks for the various Graphviz executables.  Please see
# http://www.graphviz.org/ for more information.
#
# This modules defines the following variables:
#
#   GRAPHVIZ_FOUND              = At least one Graphviz executable was found
#   GRAPHVIZ_CIRCO_EXECUTABLE   = Path to the circo command
#   GRAPHVIZ_DOT_EXECUTABLE     = Path to the dot command
#   GRAPHVIZ_FDP_EXECUTABLE     = Path to the fdp command
#   GRAPHVIZ_NEATO_EXECUTABLE   = Path to the neato command
#   GRAPHVIZ_SFDP_EXECUTABLE    = Path to the sfdp command
#   GRAPHVIZ_TWOPI_EXECUTABLE   = Path to the twopi command
#
#
#=============================================================================

#
# Find each Graphviz executable
#

LIST (APPEND myPaths "$ENV{ProgramFiles}/Graphviz 2.21/bin" 
                     "C:/Program Files/Graphviz 2.21/bin" 
                     "$ENV{ProgramFiles}/ATT/Graphviz/bin" 
                     "C:/Program Files/ATT/Graphviz/bin" 
                     "[HKEY_LOCAL_MACHINE\\SOFTWARE\\ATT\\Graphviz;InstallPath]/bin"
                     "/Applications/Graphviz.app/Contents/MacOS"
                     "/Applications/Doxygen.app/Contents/Resources"
                     "/Applications/Doxygen.app/Contents/MacOS")

SET (GRAPHVIZ_FOUND "NO")
FIND_PROGRAM (GRAPHVIZ_CIRCO_EXECUTABLE
  NAMES circo
  PATHS ${myPaths}
  DOC "Graphviz circo tool (http://www.graphviz.org/)"
)

FIND_PROGRAM (GRAPHVIZ_DOT_EXECUTABLE
  NAMES dot
  PATHS ${myPaths}
  DOC "Graphviz dot tool (http://www.graphviz.org/)"
)

FIND_PROGRAM (GRAPHVIZ_FDP_EXECUTABLE
  NAMES fdp
  PATHS ${myPaths}
  DOC "Graphviz fdp tool (http://www.graphviz.org/)"
)

FIND_PROGRAM (GRAPHVIZ_NEATO_EXECUTABLE
  NAMES neato
  PATHS ${myPaths}
  DOC "Graphviz neato tool (http://www.graphviz.org/)"
)

FIND_PROGRAM (GRAPHVIZ_SFDP_EXECUTABLE
  NAMES sfdp
  PATHS ${myPaths}
  DOC "Graphviz sfdp tool (http://www.graphviz.org/)"
)

FIND_PROGRAM (GRAPHVIZ_TWOPI_EXECUTABLE
  NAMES twopi
  PATHS ${myPaths}
  DOC "Graphviz twopi tool (http://www.graphviz.org/)"
)


IF (GRAPHVIZ_CIRCO_EXECUTABLE OR GRAPHVIZ_DOT_EXECUTABLE OR
    GRAPHVIZ_FDP_EXECUTABLE OR GRAPHVIZ_NEATO_EXECUTABLE OR
    GRAPHVIZ_SFDP_EXECUTABLE OR GRAPHVIZ_TWOPI_EXECUTABLE)
  SET (GRAPHVIZ_FOUND "YES")
ENDIF ()


