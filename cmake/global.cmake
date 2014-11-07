############################################################
##  CMake file defining global variables
##
##  Raymond Wan
##  Organizations
##    - Department of Computational Biology, Graduate School of
##      Frontier Science, University of Tokyo
##    - Computational Biology Research Center, AIST, Japan
##
##  $Id$
############################################################

##  Software versions
SET (HAMSTER_VERSION_MAJOR 1)
SET (HAMSTER_VERSION_MINOR 3)

SET (BUILDMST_VERSION_MAJOR ${HAMSTER_VERSION_MAJOR})
SET (BUILDMST_VERSION_MINOR ${HAMSTER_VERSION_MINOR})

SET (LAYOUTMST_VERSION_MAJOR ${HAMSTER_VERSION_MAJOR})
SET (LAYOUTMST_VERSION_MINOR ${HAMSTER_VERSION_MINOR})


##  Testing compilation flags, some of which are suggested by the Valgrind 3.3 book
# SET (MY_CXX_FLAGS "-pedantic -Wno-long-long -g -fno-inline -O0 -Wall")
##  Release compilation flags, suggested by the Valgrind 3.3 book
SET (MY_CXX_FLAGS "-O3 -Wall")

