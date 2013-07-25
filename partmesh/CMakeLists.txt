##-------------------------------------------------------
##          
##                  Partition Tool Project
##
##                By Wenqing Wang
##
##                     02.2012        
##-------------------------------------------------------


cmake_minimum_required(VERSION 2.4.5)

SET (VERSION_MAJOR 2)
SET (VERSION_MINOR 0)
SET (VERSION_PATCH 0)

  if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
    endif(COMMAND cmake_policy)

IF (WIN32)
	SET ( CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE TYPE INTERNAL FORCE )
    ADD_DEFINITIONS(-DWIN)	
ENDIF (WIN32)

IF (UNIX)
	string( LENGTH "$ENV{DESTDIR}" DESTDIRLEN )
	IF ( ${DESTDIRLEN} GREATER 0 )
		message( STATUS "DESTDIR specified ($ENV{DESTDIR}). Overriding 
		CMAKE_INSTALL_PREFIX (${CMAKE_INSTALL_PREFIX})" )
		set( CMAKE_INSTALL_PREFIX "" CACHE PATH "Use DESTDIR instead of CMAKE_INSTALL_PREFIX" FORCE )
	ENDIF()
ENDIF (UNIX)

SET(PROJECT_NAME partmesh)
# Project name
PROJECT( ${PROJECT_NAME} )
INCLUDE( Cmake/base.cmake )

# Option
OPTION(USE_METIS_SOURCE "Use METIS source code" OFF)

# Set build directories
# Binaries are created in /bin and libraries in /lib
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )
SET( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )

# Relative paths are used in the created projects
# May not work correcty due to the CMake docs
#SET( CMAKE_USE_RELATIVE_PATHS TRUE )

IF (MSVC)
	ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
	ADD_DEFINITIONS(-D_CRT_NONSTDC_NO_WARNINGS)
ENDIF (MSVC)


### For GNU C/CXX.
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -Wall -Wextra")
        ADD_DEFINITIONS(
#		-pedantic
		-fno-nonansi-builtins
	)
ENDIF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)

IF(USE_METIS_SOURCE) 
	MESSAGE( STATUS "The METIS package is copyrighted by the Regents of the University of
Minnesota." )
	MESSAGE( STATUS "Please read the license of the METIS package carefully before you use the METIS." )

	add_definitions(-DUSE_METIS_SOURCE -DUSE_GKREGEX)

	SET(METIS_PATH "${CMAKE_SOURCE_DIR}/metis-5.0.2")
	
	set(GKLIB_PATH "${METIS_PATH}/GKlib" CACHE PATH "path to GKlib")
	set(SHARED FALSE CACHE BOOL "build a shared library")

	if(MSVC)
		set(METIS_INSTALL FALSE)
	else()
		set(METIS_INSTALL TRUE)
	endif()

	# Configure libmetis library.
	if(SHARED)
		set(METIS_LIBRARY_TYPE SHARED)
	else()
		set(METIS_LIBRARY_TYPE STATIC)
	endif(SHARED)

	include(${GKLIB_PATH}/GKlibSystem.cmake)
	# Add include directories.
	include_directories(${GKLIB_PATH})
	include_directories(${METIS_PATH}/include)
	# Recursively look for CMakeLists.txt in subdirs.
	add_subdirectory("${METIS_PATH}/include")
	add_subdirectory("${METIS_PATH}/libmetis")
ENDIF(USE_METIS_SOURCE)


##Sow all compilation messages
SET( CMAKE_VERBOSE_MAKEFILE on )

#Options
# 
#OPTION(USE_PETSC "Use PETSc libs for parallel computing" OFF)

###SET(ROOT_PATH ..)

INCLUDE_DIRECTORIES(
.
Mesh
MatrixVec
)

#ADD_DEFINITIONS(
##
#		)

ADD_SUBDIRECTORY( MatrixVec )
ADD_SUBDIRECTORY( Mesh )
ADD_SUBDIRECTORY( METIS2ogs )



		


