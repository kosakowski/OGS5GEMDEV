/*
 * mainExtractSurface.cpp
 *
 *  Created on: Jan 26, 2011
 *      Author: TF
 */

#include "ExtractSurface.h"

// FileIO
#include "OGSIOVer4.h"
#include "XmlIO/XmlGmlInterface.h"

// GEO
#include "GEOObjects.h"
#include "Polygon.h"
#include "PolylineVec.h"
#include "ProjectData.h"

// MSH
#include "msh_lib.h" // for FEMRead
#include "msh_mesh.h"

// STL
#include <string>
Problem* aproblem = NULL;

int main (int argc, char* argv[])
{
	if (argc < 5)
	{
		std::cout << "Usage: " << argv[0] <<
		" --mesh ogs_meshfile --geometry ogs_geometry" << std::endl;
		return -1;
	}

	// *** read mesh
	std::string tmp (argv[1]);
	if (tmp.find ("--mesh") == std::string::npos)
	{
		std::cout << "could not extract mesh file name" << std::endl;
		return -1;
	}

	tmp = argv[2];
	std::string file_base_name (tmp);
	if (tmp.find (".msh") != std::string::npos)
		file_base_name = tmp.substr (0, tmp.size() - 4);

	std::vector<MeshLib::CFEMesh*> mesh_vec;
	FEMRead(file_base_name, mesh_vec);
	if (mesh_vec.empty())
	{
		std::cerr << "could not read mesh from file " << std::endl;
		return -1;
	}
	MeshLib::CFEMesh* mesh (mesh_vec[mesh_vec.size() - 1]);

	// *** read geometry
	tmp = argv[3];
	if (tmp.find ("--geometry") == std::string::npos)
	{
		std::cout << "could not extract geometry file name" << std::endl;
		return -1;
	}

	GEOLIB::GEOObjects* geo (new GEOLIB::GEOObjects);
	tmp = argv[4];
	std::string unique_name;
	std::vector<std::string> error_strings;
	FileIO::readGLIFileV4(tmp, geo, unique_name, error_strings);

	// *** get Polygon
	const std::vector<GEOLIB::Polyline*>* plys (geo->getPolylineVec (unique_name));
	if (!plys)
	{
		std::cout << "could not get vector of polylines" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	MeshLib::ExtractSurface extract_surface (mesh);
	const size_t n_plys (plys->size());

	ProjectData* project_data (new ProjectData);
	project_data->setGEOObjects (geo);
	FileIO::XmlGmlInterface xml_out (project_data, "OpenGeoSysGLI.xsd");
	for (size_t k(0); k < n_plys; k++) {
		std::vector<GEOLIB::Point*>* sfc_pnts (new std::vector<GEOLIB::Point*>);
		std::vector<GEOLIB::Surface*>* surfaces (new std::vector<GEOLIB::Surface*>);

		if ( ((*plys)[k])->isClosed() ) {
			GEOLIB::Polygon polygon (*((*plys)[k]));
			GEOLIB::Surface* surface(extract_surface.extractSurface(polygon, *sfc_pnts));
			surfaces->push_back (surface);
			std::cout << "number of triangles in surface " << k << ": " << surface->getNTriangles() << std::endl;
		}

		if (! surfaces->empty()) {
			std::string fname ("PointsForSurface");
			fname += number2str (k);
			geo->addPointVec (sfc_pnts, fname);
			geo->addSurfaceVec (surfaces, fname);

			std::string out_fname ("Surface");
			out_fname += number2str (k);
			out_fname += ".gml";
			xml_out.setNameForExport(fname);
			xml_out.writeToFile(out_fname);
		} else {
			std::cout << "could not create any surface" << std::endl;
		}
	}

	delete mesh;
	delete project_data;
}
