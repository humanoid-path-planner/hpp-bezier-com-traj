# - Try to find libcdd
# Once done this will define
#  CDD_FOUND - System has CDD
#  CDD_INCLUDE_DIRS - The CDD include directories
#  CDD_LIBRARIES - The libraries needed to use CDD
#  CDD_DEFINITIONS - Compiler switches required for using CDD


find_path(CDD_INCLUDE_DIR cdd/cdd.h
          HINTS ${CDD_INCLUDEDIR} /usr/include
          PATH_SUFFIXES CDD )

find_library(CDD_LIBRARY NAMES libcdd.so
             HINTS ${CDD_LIBDIR} ${CDD_LIBRARY_DIRS} /usr/lib/libcdd.so )

set(CDD_LIBRARIES ${CDD_LIBRARY} )
set(CDD_INCLUDE_DIRS ${CDD_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CDD_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CDD  DEFAULT_MSG
                                  CDD_LIBRARY CDD_INCLUDE_DIR)

mark_as_advanced(CDD_INCLUDE_DIR CDD_LIBRARY )
