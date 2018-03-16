# - Try to find libcdd
# Once done this will define
#  CDL_FOUND - System has CDL
#  CDL_INCLUDE_DIRS - The CDL include directories
#  CDL_LIBRARIES - The libraries needed to use CDL
#  CDL_DEFINITIONS - Compiler switches required for using CDL


find_path(CDL_INCLUDE_DIR centroidal-dynamics-lib/centroidal_dynamics.h
          HINTS ${CDL_INCLUDEDIR} /usr/include
          PATH_SUFFIXES CDL )

find_library(CDL_LIBRARY NAMES libcentroidal-dynamics-lib.so
             HINTS ${CDL_LIBDIR} ${CDL_LIBRARY_DIRS} /usr/lib/libcentroidal-dynamics-lib.so )

set(CDL_LIBRARIES ${CDL_LIBRARY} )
set(CDL_INCLUDE_DIRS ${CDL_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CDL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CDL  DEFAULT_MSG
                                  CDL_LIBRARY CDL_INCLUDE_DIR)

mark_as_advanced(CDL_INCLUDE_DIR CDL_LIBRARY )
