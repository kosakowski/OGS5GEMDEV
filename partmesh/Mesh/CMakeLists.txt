# Source files

set( MESH_SOURCES
      Grain.h
      Grain.cpp    
      Node.h
	  Node.cpp
      Edge.h
      Edge.cpp
      Elem.h
      Elem.cpp
      Mesh.h
      Mesh.cpp
)

#include_directories(
#	.
#	../MatrixVec
#)


# Create the library
add_library(Mesh STATIC ${MESH_SOURCES} )

target_link_libraries(
	Mesh
	Math
)

