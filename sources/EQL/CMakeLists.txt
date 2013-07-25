# Source files
SET ( SOURCES 
	cacint.h
	cacint.cpp
	eqlink.h
	eqlink.cpp
)

# Create the library
ADD_LIBRARY ( EQL STATIC ${SOURCES} )

INCLUDE_DIRECTORIES (
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/Mathlib
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/FEM
)
