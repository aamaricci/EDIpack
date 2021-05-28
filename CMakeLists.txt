##################################################
# DMFTED PROJECT
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)
PROJECT(edipack Fortran)
SET(VERSION 1.0.1)

# Define the library name
SET(LANCEDLIB edipack)


MESSAGE(STATUS "OS: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
SET(BUILD_TYPE "RELEASE" CACHE STRING "Build type")
STRING(TOUPPER "${BUILD_TYPE}" BUILD_TYPE)
MESSAGE(STATUS "BUILD_TYPE=${BUILD_TYPE}")
SET(CMAKE_BUILD_TYPE ${BUILD_TYPE} CACHE INTERNAL "Build type" FORCE)

# Add our local CMAKE modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)

# Check that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
  MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

IF( (${CMAKE_Fortran_COMPILER_ID} MATCHES Intel) OR (${CMAKE_Fortran_COMPILER_ID} MATCHES GNU))
  MESSAGE(STATUS "Fortran Compiler id   = ${CMAKE_Fortran_COMPILER_ID}")
  MESSAGE(STATUS "Fortran Compiler ver. = ${CMAKE_Fortran_COMPILER_VERSION}")
ELSEIF()
  MESSAGE(FATAL_ERROR "Unsupported Fortran compiler (use Intel or GNU). Try export FC=<your FC compiler> ")
ENDIF()  

############################################################
# Set Fortran options based on BUILD_TYPE and FC ID
############################################################
# -mcmodel=large  this is to remove the 2Gb limit of virtual memory allocation
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
  SET(CMAKE_Fortran_MODDIR_FLAG   "-J")
  SET(CMAKE_Fortran_FLAGS         "-cpp -ffree-line-length-none -fPIC -w ")
  IF(CMAKE_Fortran_COMPILER_VERSION LESS 10.0.0)
    SET(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wno-argument-mismatch")
  ELSE()
    SET(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  ENDIF()
  SET(CMAKE_Fortran_FLAGS_TESTING "-O2 -funroll-loops")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -p -g -Wsurprising -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -fbounds-check")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3   -funroll-loops")   
elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  SET(CMAKE_Fortran_MODDIR_FLAG   "-module ") #remember the ending white space here 
  SET(CMAKE_Fortran_FLAGS         "-fpp")
  SET(CMAKE_Fortran_FLAGS_TESTING "-O2 -ftz")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-p -O0 -g -fpe1 -warn -debug extended -traceback -check all,noarg_temp_created")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ftz")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES G95)
  SET(CMAKE_Fortran_MODDIR_FLAG   "-fmod=")
  SET(CMAKE_Fortran_FLAGS         "-cpp")
  SET(CMAKE_Fortran_FLAGS_TESTING "-O1  -fsloppy-char")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -ftrace=full")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  
elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  SET(CMAKE_Fortran_MODDIR_FLAG   "-module ")
  SET(CMAKE_Fortran_FLAGS         "-")
  SET(CMAKE_Fortran_FLAGS         "")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Mframe")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mcmodel=medium -fast -Munroll")
endif()

IF( "${BUILD_TYPE}" MATCHES "DEBUG")
  MESSAGE(STATUS "Fortran Compiler options = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
ELSEIF("${BUILD_TYPE}" MATCHES "TESTING")
  MESSAGE(STATUS "Fortran Compiler options = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_TESTING}")
ELSEIF("${BUILD_TYPE}" MATCHES "RELEASE")
  MESSAGE(STATUS "Fortran Compiler options = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
ENDIF()


# Set some options the user may choose
OPTION(USE_MPI "Set MPI pre-compilation flag" ON)
OPTION(VERBOSE "Build a verbose Makefile" OFF)
OPTION(WITH_BLAS_LAPACK "Use internal Blas/Lapack" OFF)
IF(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE ON)
ELSE(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE OFF)
ENDIF(VERBOSE)



############################################################
# Look for MPI if needed
############################################################
# reSet compiler to MPI is required
IF(USE_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  IF(MPI_Fortran_FOUND)
    SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
    MESSAGE(STATUS "Set Fortran compiler to ${Yellow}${CMAKE_Fortran_COMPILER}${ColourReset}, ID=${CMAKE_Fortran_COMPILER_ID}")
    SET(MPI_CPP "MPI")		#pre-processing option
  ELSE()
    MESSAGE(FATAL_ERROR "MPI Found but No MPI-Fortran compiler can be determined.")    
  ENDIF()
ELSE(USE_MPI)
  SET(MPI_CPP "")
ENDIF(USE_MPI)

#set MPI preprocessing defition:
IF(USE_MPI)
  ADD_DEFINITIONS(-D_MPI)
ELSE(USE_MPI)
  ADD_DEFINITIONS(-D_)
ENDIF(USE_MPI)

IF( "${BUILD_TYPE}" MATCHES "DEBUG")
  ADD_DEFINITIONS(-D_DEBUG)
ENDIF()

############################################################
# Check for .git and extract the current SHA1 revision.
############################################################
SET(GIT_SHA1 0)
SET(GIT_BRANCH "master")
SET(GIT_TAG ${VERSION})
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  SET(GIT_REPO TRUE)
  FIND_PACKAGE(Git)
  IF(GIT_FOUND)
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_SHA1
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "${Yellow}Git SHA1: ${ColourReset}${GIT_SHA1}" )
    #
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_BRANCH
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "${Yellow}Git branch: ${ColourReset}${GIT_BRANCH}" )
    #
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} describe --tags
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_TAG
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE(STATUS "${Yellow}Git tag: ${ColourReset}${GIT_TAG}")
  ENDIF()
ELSE()
  SET(GIT_REPO FALSE)
ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
SET(VERSION ${GIT_TAG})    	#Update version





############################################################
# Look for required libraries: scifor
# (includes Blas/Lapack/Scalapack/Quadpack/Minpack/FFTpack/Arpack/...) 
############################################################
# Find Scifor
FIND_PACKAGE(PkgConfig)
PKG_SEARCH_MODULE(SCIFOR REQUIRED scifor)
IF(${SCIFOR_FOUND})
  MESSAGE(STATUS "SCIFOR found at: ${SCIFOR_PREFIX}")
ELSE()
  MESSAGE(FATAL_ERROR "
Can not find SciFortran in this system. Check SciFor installation. Try adding the provided PKG-CONFIG file to the environment:
$ export PKG_CONFIG_PATH=<path_to_scifor_dir>/etc:$PKG_CONFIG_PATH")
ENDIF()
IF("${SCIFOR_LDFLAGS}" MATCHES "blas" OR "lapack")
  FIND_PACKAGE(BLAS REQUIRED)
  FIND_PACKAGE(LAPACK REQUIRED)
ELSEIF("${SCIFOR_LDFLAGS}" MATCHES "mkl")
  FIND_PACKAGE(MKL REQUIRED)
ENDIF()
INCLUDE_DIRECTORIES(BEFORE ${SCIFOR_INCLUDE_DIRS})










############################################################
# Define the actual files and folders that make up the build
############################################################
SET(USER_HOME $ENV{HOME})
SET(USER $ENV{USER})

STRING(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" FC_ID)

SET(LIB_PATH  "$ENV{HOME}/opt/${PROJECT_NAME}/${FC_ID}")
SET(FULL_VER "${FC_ID}")
IF("${CMAKE_BUILD_TYPE}" MATCHES "DEBUG")
  STRING(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)
  SET(LIB_PATH  "${LIB_PATH}/${BUILD_TYPE}")
  SET(FULL_VER  "${FULL_VER}/${BUILD_TYPE}")
ENDIF()

IF(GIT_REPO)
  IF(NOT GIT_BRANCH MATCHES "master")
    SET(LIB_PATH  "${LIB_PATH}/${GIT_BRANCH}")
    SET(FULL_VER  "${FULL_VER}/${GIT_BRANCH}")
  ENDIF()
ENDIF()

SET(FULL_VER  "${FULL_VER}/${VERSION}")
SET(PREFIX  "${LIB_PATH}/${VERSION}" CACHE PATH "Prefix prepended to install directories")
SET(CMAKE_INSTALL_PREFIX "${PREFIX}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)



# Define some directories
SET(LIB_SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB_ETC ${CMAKE_SOURCE_DIR}/etc)
SET(LIB_ENV ${CMAKE_SOURCE_DIR}/etc/environment_modules)
SET(LIB_TMP_INC ${CMAKE_BINARY_DIR}/include)
SET(LIB_TMP_ETC ${CMAKE_BINARY_DIR}/etc)
SET(LIB_TMP_VER ${CMAKE_BINARY_DIR}/version)
SET(LIB_TARGET_LIB ${CMAKE_INSTALL_PREFIX}/lib)
SET(LIB_TARGET_INC ${CMAKE_INSTALL_PREFIX}/include)
SET(LIB_TARGET_ETC ${CMAKE_INSTALL_PREFIX}/etc)
SET(LIB_TARGET_DIR ${CMAKE_INSTALL_PREFIX})

# Have the .mod files placed in the INC folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB_TMP_INC})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_Fortran_MODULE_DIRECTORY})




ADD_SUBDIRECTORY(${LIB_SRC}/BATH)

ADD_SUBDIRECTORY(${LIB_SRC}/IO)

ADD_SUBDIRECTORY(${LIB_SRC}/DIAG)

ADD_SUBDIRECTORY(${LIB_SRC}/GF)

ADD_SUBDIRECTORY(${LIB_SRC}/CHI)

ADD_SUBDIRECTORY(${LIB_SRC})


ADD_LIBRARY(edipack STATIC
  $<TARGET_OBJECTS:INIT_LIB>
  $<TARGET_OBJECTS:BATH_LIB>
  $<TARGET_OBJECTS:IO_LIB>
  $<TARGET_OBJECTS:DIAG_LIB>
  $<TARGET_OBJECTS:GF_LIB>
  $<TARGET_OBJECTS:CHI_LIB>
  $<TARGET_OBJECTS:GLOB_LIB>
  )
ADD_DEPENDENCIES(IO_LIB
  BATH_LIB
  )


ADD_DEPENDENCIES(GLOB_LIB
  INIT_LIB
  BATH_LIB
  IO_LIB
  DIAG_LIB
  GF_LIB
  CHI_LIB
  )



SET(PERMISSION_777
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_WRITE GROUP_READ GROUP_EXECUTE
  WORLD_WRITE WORLD_READ WORLD_EXECUTE)



#Build the Edipack module for environment modules
FILE(REMOVE_RECURSE ${LIB_TMP_ETC}/modules/)
SET(TMP_MODULE_NAME "${FULL_VER}" CACHE PATH "Prefix prepended to install directories")
SET(TMP_ENV_MODULE_FILE ${LIB_TMP_ETC}/modules/${PROJECT_NAME}/${TMP_MODULE_NAME})
SET(TMP_VER_MODULE_FILE ${LIB_TMP_ETC}/modules/${PROJECT_NAME}/${FC_ID}/.version)
CONFIGURE_FILE(${LIB_ENV}/module.in ${TMP_ENV_MODULE_FILE} @ONLY)
CONFIGURE_FILE(${LIB_ENV}/version.in ${TMP_VER_MODULE_FILE} @ONLY)

#Build the PKG-CONFIG file
SET(TMP_PKCONFIG_FILE ${LIB_TMP_ETC}/${PROJECT_NAME}.pc)
CONFIGURE_FILE( ${LIB_ETC}/${PROJECT_NAME}.pc.in ${TMP_PKCONFIG_FILE} @ONLY)

#Build the user config script (sourced in user shell config file, i.e. .bashrc)
SET(CONFIG_SCRIPT edipack_config.sh)
SET(TMP_CONFIGVARS_FILE ${LIB_TMP_ETC}/${CONFIG_SCRIPT})
CONFIGURE_FILE(${LIB_ETC}/config_edipack.sh.in ${TMP_CONFIGVARS_FILE} @ONLY)


FILE(WRITE  ${LIB_TMP_VER}  "${VERSION}\n")

INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${LIB_TARGET_INC})

INSTALL(TARGETS edipack DESTINATION ${LIB_TARGET_LIB})

INSTALL(DIRECTORY ${LIB_TMP_ETC}/ DESTINATION ${LIB_TARGET_ETC})

INSTALL(FILES ${TMP_CONFIGVARS_FILE} DESTINATION ${LIB_TARGET_ETC}/
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(FILES ${TMP_VER_FILE} DESTINATION ${LIB_TARGET_DIR} 
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(FILES ${LIB_TARGET_ETC}/${PROJECT_NAME}.pc DESTINATION $ENV{HOME}/.pkgconfig.d/
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(DIRECTORY ${LIB_TARGET_ETC}/modules/ DESTINATION $ENV{HOME}/.modules.d)


MESSAGE( STATUS "${Red}Library version:${ColourReset} ${VERSION}")
MESSAGE( STATUS "${Red}Library will be installed in:${ColourReset} ${CMAKE_INSTALL_PREFIX}")
MESSAGE( STATUS "
>> ${Red}TO CONCLUDE INSTALLATION${ColourReset} <<
Compile with:
$ make
Install with:
$ make install

Uninstall with:
$ make uninstall

Add library to your OS: 
pick ONE method below  [add it in your bash profile to make it permanent (e.g. ~/.bashrc)]:
${Yellow}Method 1: use the provided ${PROJECT_NAME} environment module${ColourReset}:
   $ module use $HOME/.modules.d
   $ module load ${PROJECT_NAME}/${FC_PLAT}

${Yellow}Method 2: source the config script${ColourReset}:
   $ source ${LIB_TARGET_ETC}/${CONFIG_SCRIPT}

${Yellow}Method 3: use pkg-config with the provided ${PROJECT_NAME}.pc${ColourReset}:
   $ export PKG_CONFIG_PATH=${LIB_TARGET_ETC}/:$PKG_CONFIG_PATH
   $ pkg-config --cflags --libs ${PROJECT_NAME} (to get lib info)
")


# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_MODULE_PATH}/distclean.cmake
)

# Uninstall target
if(NOT TARGET uninstall)
  CONFIGURE_FILE(
    "${CMAKE_MODULE_PATH}/ConfigUninstall.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
    "${PROJECT_NAME}/${TMP_MODULE_NAME}" "${PROJECT_NAME}.pc" )
  
ENDIF()
