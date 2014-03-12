##-------------------------------------------------------
##          
##    Tool to convert GIS shape data into FEM data
##                By Wenqing Wang
##
##                     03.2010        
##-------------------------------------------------------

INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/UTL/MSHGEOTOOLS
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO	
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/MSHGEOTOOLS
)
  
SET(GIS2FEM_SOURCES main.cpp misc.cpp misc.h )
    
add_executable(fem_toolkits  ${GIS2FEM_SOURCES})

TARGET_LINK_LIBRARIES( fem_toolkits
	MSH
	MSHGEOTOOLS
	FileIO
)	