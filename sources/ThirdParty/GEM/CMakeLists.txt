# Source files

# Recursively collect all header files from the current directory
file(GLOB HEADERS *.h)

# Recursively collect all source files from the current directory
file(GLOB SOURCES *.cpp)

# Create the library
add_library( GEM STATIC 
	     ${SOURCES}
             ${HEADERS} 
)
include_directories(
        ${CMAKE_SOURCE_DIR}/GEM
)


##set_target_properties(GEM PROPERTIES LINKER_LANGUAGE CXX)

if (WIN32)
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPMGEMPLUGIN")
endif ()

if (UNIX)
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPMGEMPLUGIN -D__unix")
endif ()
