# - Try to find libcdd
# Once done this will define
#  CLP_FOUND - System has CLP
#  CLP_INCLUDE_DIRS - The CLP include directories
#  CLP_LIBRARIES - The libraries needed to use CLP
#  CLP_DEFINITIONS - Compiler switches required for using CLP

# /usr/include/coin, /usr/lib/libClp.so

find_path(CLP_INCLUDE_DIR coin/ClpSimplex.hpp
          HINTS ${CLP_INCLUDEDIR}
          PATH_SUFFIXES CLP )

find_library(CLP_LIBRARY NAMES libclp.so
             HINTS ${CLP_LIBDIR} ${CLP_LIBRARY_DIRS} )

set(CLP_LIBRARIES ${CLP_LIBRARY} )
set(CLP_INCLUDE_DIRS ${CLP_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CDD_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CLP  DEFAULT_MSG
                                  CLP_LIBRARY CLP_INCLUDE_DIR)

mark_as_advanced(CLP_INCLUDE_DIR CLP_LIBRARY )
