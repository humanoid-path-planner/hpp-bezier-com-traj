/*
 * Copyright 2017, LAAS-CNRS
 * Author: Steve Tonneau
 */

#ifndef _BEZIER_COM_TRAJ_LIB_CONFIG_HH
#define _BEZIER_COM_TRAJ_LIB_CONFIG_HH

// Package version (header).
# define BEZIER_COM_TRAJ_VERSION "UNKNOWN"

// Handle portable symbol export.
// Defining manually which symbol should be exported is required
// under Windows whether MinGW or MSVC is used.
//
// The headers then have to be able to work in two different modes:
// - dllexport when one is building the library,
// - dllimport for clients using the library.
//
// On Linux, set the visibility accordingly. If C++ symbol visibility
// is handled by the compiler, see: http://gcc.gnu.org/wiki/Visibility
# if defined _WIN32 || defined __CYGWIN__
// On Microsoft Windows, use dllimport and dllexport to tag symbols.
#  define BEZIER_COM_TRAJ_DLLIMPORT __declspec(dllimport)
#  define BEZIER_COM_TRAJ_DLLEXPORT __declspec(dllexport)
#  define BEZIER_COM_TRAJ_DLLLOCAL
# else
// On Linux, for GCC >= 4, tag symbols using GCC extension.
#  if __GNUC__ >= 4
#   define BEZIER_COM_TRAJ_DLLIMPORT __attribute__ ((visibility("default")))
#   define BEZIER_COM_TRAJ_DLLEXPORT __attribute__ ((visibility("default")))
#   define BEZIER_COM_TRAJ_DLLLOCAL  __attribute__ ((visibility("hidden")))
#  else
// Otherwise (GCC < 4 or another compiler is used), export everything.
#   define BEZIER_COM_TRAJ_DLLIMPORT
#   define BEZIER_COM_TRAJ_DLLEXPORT
#   define BEZIER_COM_TRAJ_DLLLOCAL
#  endif // __GNUC__ >= 4
# endif // defined _WIN32 || defined __CYGWIN__

# ifdef BEZIER_COM_TRAJ_STATIC
// If one is using the library statically, get rid of
// extra information.
#  define BEZIER_COM_TRAJ_DLLAPI
#  define BEZIER_COM_TRAJ_LOCAL
# else
// Depending on whether one is building or using the
// library define DLLAPI to import or export.
#  ifdef BEZIER_COM_TRAJ_EXPORTS
#   define BEZIER_COM_TRAJ_DLLAPI BEZIER_COM_TRAJ_DLLEXPORT
#  else
#   define BEZIER_COM_TRAJ_DLLAPI BEZIER_COM_TRAJ_DLLIMPORT
#  endif // BEZIER_COM_TRAJ_EXPORTS
#  define BEZIER_COM_TRAJ_LOCAL BEZIER_COM_TRAJ_DLLLOCAL
# endif // BEZIER_COM_TRAJ_STATIC

#endif //_BEZIER_COM_TRAJ_LIB_CONFIG_HH
