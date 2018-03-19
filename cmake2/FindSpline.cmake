# - Try to find spline
# Once done this will define
#  SPLINE_FOUND - System has SPLINE
#  SPLINE_INCLUDE_DIRS - The SPLINE include directories
#  SPLINE_DEFINITIONS - Compiler switches required for using SPLINE

# /usr/include/coin, /usr/lib/libSPLINE.so

find_path(SPLINE_INCLUDE_DIR spline/bezier_curve.h
          HINTS ${SPLINE_INCLUDEDIR}
          PATH_SUFFIXES SPLINE )

#~ find_library(SPLINE_LIBRARY NAMES libSPLINE.so
             #~ HINTS ${SPLINE_LIBDIR} ${SPLINE_LIBRARY_DIRS} )

set(SPLINE_INCLUDE_DIRS ${SPLINE_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CDD_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SPLINE  DEFAULT_MSG
                                  SPLINE_INCLUDE_DIR)

mark_as_advanced(SPLINE_INCLUDE_DIR)
