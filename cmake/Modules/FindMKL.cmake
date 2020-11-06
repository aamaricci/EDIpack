# The following are set after the configuration is done:
#  MKL_FOUND        -  system has MKL
#  MKL_ROOT_DIR     -  path to the MKL base directory
#
#
# Sample usage:
#    find_package(MKL REQUIRED)
#
# AUTHOR
# Adriano Amaricci (adriano.amaricci.AT.sissa.it)



SET(CMAKE_FIND_DEBUG_MODE 0)

# Find MKL ROOT
FIND_PATH(MKL_ROOT_DIR NAMES include/mkl.h include/mkl.fi  PATHS $ENV{MKLROOT})

# Convert symlinks to real paths

GET_FILENAME_COMPONENT(MKL_ROOT_DIR ${MKL_ROOT_DIR} REALPATH)

IF (NOT MKL_ROOT_DIR)
  IF (MKL_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find MKL: please set environment variable {MKLROOT}")
  ELSE()
    UNSET(MKL_ROOT_DIR CACHE)
  ENDIF()
  RETURN()
ENDIF()
  


INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_ROOT_DIR)

MARK_AS_ADVANCED(MKL_ROOT_DIR)


if (CMAKE_FIND_DEBUG_MODE)
  MESSAGE(STATUS "Found MKL:${MKL_ROOT_DIR}")
endif()

