############################################################
##  CMake file defining add_subdirectory_once
##
##  Raymond Wan
##  Organizations
##    - Department of Computational Biology, Graduate School of
##      Frontier Science, University of Tokyo
##    - Computational Biology Research Center, AIST, Japan
##
##  $Id$
############################################################


##  Source code for this function from 
##    Rolf Eike Beer on the CMake mailing list
##    (2010/12/02)
function (Add_Subdirectory_Once SUBDIRECTORY)
  get_filename_component (FULLPATH ${SUBDIRECTORY} REALPATH)

  GET_PROPERTY (_INCLUDED_DIRS GLOBAL PROPERTY ADD_SUBDIRECTORY_ONCE_INCLUDED)
  LIST (FIND _INCLUDED_DIRS "${FULLPATH}" _USED_INDEX)

  if (_USED_INDEX EQUAL -1)
    SET_PROPERTY (GLOBAL PROPERTY ADD_SUBDIRECTORY_ONCE_INCLUDED "${_INCLUDED_DIRS}" "${FULLPATH}")
    if(${ARGC} EQUAL 1)
      add_subdirectory (${SUBDIRECTORY})
    else(${ARGC} EQUAL 1)
      add_subdirectory (${SUBDIRECTORY} ${ARGV1})
    endif (${ARGC} EQUAL 1)
  endif (_USED_INDEX EQUAL -1)
endfunction (Add_Subdirectory_Once)



