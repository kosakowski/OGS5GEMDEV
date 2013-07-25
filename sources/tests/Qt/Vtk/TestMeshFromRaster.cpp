/*
 * testMeshFromRaster.cpp
 * 2012/01/13 KR Initial implementation
 */

#include "gtest.h"
#include "TestHelperFunctions.h"

#include "GridAdapter.h"
#include "VtkMeshConverter.h"
#include "VtkGeoImageSource.h"

#include <sstream>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>

TEST(VTK, TestMeshFromRaster)
{
	QString fileName = getTestdataInputDir();
	fileName += "testMeshFromRaster.asc";
	vtkSmartPointer<VtkGeoImageSource> geo_image = vtkSmartPointer<VtkGeoImageSource>::New();
	geo_image->readImage(fileName);
	vtkSmartPointer<vtkImageData> image = geo_image->GetOutput();
	image->Update();
	
	double origin[3];
	geo_image->getOrigin(origin);

	GridAdapter* grid = VtkMeshConverter::convertImgToMesh(image, origin, geo_image->getSpacing(), MshElemType::TRIANGLE, UseIntensityAs::ELEVATION);
		
	// Correct number of nodes?
	ASSERT_EQ((size_t)626, grid->getNodes()->size());
	
	// Correct number of elements?
	ASSERT_EQ((size_t)1082, grid->getElements()->size());

	FileIO::OGSMeshIO meshIO;
	meshIO.setMesh(grid->getCFEMesh());
	meshIO.setPrecision(6);
	meshIO.setFormat(ios::fixed);
	
	compareToReference(meshIO.writeToString().c_str(), QString("testMeshFromRaster_result.msh"));
}
