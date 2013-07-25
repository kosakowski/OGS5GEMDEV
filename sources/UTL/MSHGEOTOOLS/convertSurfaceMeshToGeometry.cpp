/*
 * convertSurfaceMeshToGeometry.cpp
 *
 *  Created on: May 30, 2012
 *      Author: fischeth
 */

// STL
#include <vector>
#include <string>

// FileIO
#include "OGSIOVer4.h"
#include "XmlIO/XmlGmlInterface.h"

// FileIO/GeoIO
#include "GeoIO/Gmsh2GeoIO.h"
// FileIO/MeshIO
#include "GMSHInterface.h"

// GEO
#include "GEOObjects.h"
#include "Point.h"
#include "PolylineVec.h"
#include "ProjectData.h"

int main (int argc, char* argv[])
{
	if (argc < 5) {
		std::cout << "Usage: " << argv[0] << " --mesh gmsh_mesh --geometry new_geo_file" << std::endl;
		return -1;
	}

	GEOLIB::GEOObjects* geo (new GEOLIB::GEOObjects);

	// *** read GMSH mesh
	std::string gmsh_msh_fname(argv[2]);
	FileIO::Gmsh2GeoIO::loadMeshAsGeometry(gmsh_msh_fname, geo);

	ProjectData* project_data (new ProjectData);
	project_data->setGEOObjects (geo);
	FileIO::XmlGmlInterface xml_out (project_data, "OpenGeoSysGLI.xsd");

	std::string geo_fname(argv[4]);
	xml_out.setNameForExport(gmsh_msh_fname);
	xml_out.writeToFile(geo_fname);

	// get Surface for writing a OpenGeoSys TIN
	std::vector<GEOLIB::Surface*> const& sfcs(*(geo->getSurfaceVec(gmsh_msh_fname)));
	GEOLIB::Surface const*const sfc(sfcs[0]);
	std::ofstream out((geo_fname+".tin").c_str());
	for (size_t k(0); k<sfc->getNTriangles(); k++) {
		GEOLIB::Triangle const& tri_k(*((*sfc)[k]));
		for (size_t j(0); j<3; j++) {
			out << *(tri_k.getPoint(j)) << " " << std::flush;
		}
		out << std::endl;
	}
	out.close();

	delete project_data;
}
