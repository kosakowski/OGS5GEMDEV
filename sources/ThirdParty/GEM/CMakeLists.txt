# Source files


if(USE_SPDLOG_PRECOMPILED)
   if(NOT TARGET spdlog::spdlog)
       find_package(spdlog CONFIG REQUIRED)
   endif()
endif()



# Recursively collect all header files from the current directory
file(GLOB HEADERS *.h)

# Recursively collect all source files from the current directory
file(GLOB SOURCES *.cpp)

# Create the library
add_library( GEM STATIC 
	     ${SOURCES}
             ${HEADERS} 
)
# add fmt, thermofun und spdlog  library
target_link_libraries( GEM fmt )
#target_link_libraries( GEM ThermoFun)
target_link_libraries( GEM spdlog)

include_directories(
        ${CMAKE_SOURCE_DIR}/GEM
)




##set_target_properties(GEM PROPERTIES LINKER_LANGUAGE CXX)

if (WIN32)
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPMGEMPLUGIN")
endif ()

if (UNIX)
	set( CMAKE_CXX_FLAGS "-O3 -DNDEBUG -DIPMGEMPLUGIN -D__unix -DUSE_NLOHMANNJSON -std=c++17")
#	set( CMAKE_CXX_FLAGS "-O3 -DNDEBUG -DIPMGEMPLUGIN -D__unix -DUSE_NLOHMANNJSON -DUSE_THERMOFUN -std=c++17")
endif ()
