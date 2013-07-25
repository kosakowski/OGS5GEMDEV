# Source files

set( Math_SOURCES
      Matrix.h
      Matrix.cpp    
      SymMatrix.h
      SymMatrix.cpp
      vec.h
      vec.cpp
)

#include_directories(
#	.
#	../Mesh
#)


# Create the library
add_library(Math STATIC ${Math_SOURCES} )

target_link_libraries(
	Math
)

