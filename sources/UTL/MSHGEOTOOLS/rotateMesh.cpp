/*
 * rotateMesh.cpp
 *
 *  Created on: Jan 30, 2012
 *      Author: TF
 */

#include <ctime>

// FileIO
#include "MeshIO/OGSMeshIO.h"
#include "XmlIO/XmlGmlInterface.h"

// GeoLib
#include "AxisAlignedBoundingBox.h"
#include "GEOObjects.h"
#include "ProjectData.h"

// MathLib
#include "AnalyticalGeometry.h"
#include "Vector3.h"

// MSH
#include "msh_lib.h" // for FEMRead

// we need this for using the xml functions of Qt
#include <QApplication>

int main (int argc, char* argv[])
{
	// Creating a non-gui (console) Qt application
	QApplication app(argc, argv, false);

	if (argc < 9) {
		std::cout << "Usage: " << argv[0] << " --mesh-in meshfile --mesh-out mesh-out-file "
						<< std::flush;
		std::cout << "--geo-in geo-file --geo-out geo-out-file [--upside-down true/false] "
						<< std::endl;
		return -1;
	}

	// *** parsing mesh options from command line
	// *** parsing mesh in file name
	std::string tmp(argv[1]);
	if (tmp.find("--mesh-in") == std::string::npos) {
		std::cout << "could not find switch for reading mesh file name" << std::endl;
		return -1;
	}
	tmp = argv[2];
	std::string file_base_name(tmp);
	if (tmp.find(".msh") != std::string::npos) file_base_name = tmp.substr(0, tmp.size() - 4);

	// *** parsing mesh output name
	tmp = argv[3];
	if (tmp.find("--mesh-out") == std::string::npos) {
		std::cout << "could not find switch for file name for writing the mesh" << std::endl;
		return -1;
	}
	std::string mesh_out_fname(argv[4]);

	// *** parsing geometry options from command line
	tmp = argv[5];
	if (tmp.find("--geo-in") == std::string::npos) {
		std::cout << "could not find switch for reading geometry file name" << std::endl;
		return -1;
	}
	std::string geo_fname_in(argv[6]);

	// *** read output name
	tmp = argv[7];
	if (tmp.find("--geo-out") == std::string::npos) {
		std::cout << "could not find switch for file name for writing the geometry" << std::endl;
		return -1;
	}
	std::string geo_fname_out(argv[8]);

	// *** read mesh
	std::vector<MeshLib::CFEMesh*> mesh_vec;
	FEMRead(file_base_name, mesh_vec);
	if (mesh_vec.empty()) {
		std::cerr << "could not read mesh from file " << std::endl;
		return -1;
	}
	MeshLib::CFEMesh* mesh (mesh_vec[mesh_vec.size() - 1]);
	mesh->setNumberOfNodesFromNodesVectorSize();

	// *** read geometry
	GEOLIB::GEOObjects* geo_objs (new GEOLIB::GEOObjects);
	ProjectData* project_data (new ProjectData);
	project_data->setGEOObjects (geo_objs);

	std::string schema_name("./OpenGeoSysGLI.xsd");
	FileIO::XmlGmlInterface xml(project_data, schema_name);
	xml.readFile(QString::fromStdString (geo_fname_in));
	std::vector<std::string> original_geo_names;
	geo_objs->getGeometryNames(original_geo_names);
	std::vector<GEOLIB::Point*> const& original_geo_pnts(*(geo_objs->getPointVec(original_geo_names[0])));
	const size_t n_geo_pnts(original_geo_pnts.size());
	// create a copy of the original geometric points
	std::vector<GEOLIB::Point*> *geo_pnts(new std::vector<GEOLIB::Point*>);
	for (size_t k(0); k<n_geo_pnts; k++) {
		geo_pnts->push_back(new GEOLIB::Point(original_geo_pnts[k]->getData()));
	}

	// *** transfer mesh nodes to points
	// to make a copy is for this small program useless, but it is more clear and
	// if this functions will become part of official OGS we are ready to use it
	MeshLib::CFEMesh mesh_copy(*mesh);

	std::vector<GEOLIB::Point*> pnts;
	const size_t n_nodes(mesh_copy.GetNodesNumber(false));
	std::vector<MeshLib::CNode*>& nodes(const_cast<std::vector<MeshLib::CNode*>& >(mesh_copy.getNodeVector()));
	for (size_t k(0); k<n_nodes; k++) {
		pnts.push_back(new GEOLIB::Point(nodes[k]->getData()));
	}

	MathLib::Vector plane_normal(0.0, 0.0, 0.0);
	double d(0.0);
	MathLib::getNewellPlane(pnts, plane_normal, d);
	MathLib::Vector geo_plane_normal(plane_normal);
	MathLib::rotatePointsToXZ(plane_normal, pnts);
	MathLib::rotatePointsToXZ(geo_plane_normal, *geo_pnts);

	bool switch_upside_down(true);
	if (argc > 10) {
		tmp = argv[9];
		if (tmp.find("--upside-down") != std::string::npos) {
			tmp = argv[10];
			if (tmp.find("true") != std::string::npos || tmp.find("TRUE") != std::string::npos) {
				switch_upside_down = true;
			} else {
				switch_upside_down = false;
			}
		}
	}

	if (switch_upside_down) {
		GEOLIB::AABB aabb(&pnts);
		double val(aabb.getMinPoint()[2] + aabb.getMaxPoint()[2]);
		for (size_t k(0); k<n_nodes; k++) {
			(*(pnts[k]))[2] = val - (*(pnts[k]))[2];
		}
		for (size_t k(0); k<n_geo_pnts; k++) {
			(*((*geo_pnts)[k]))[2] = val - (*((*geo_pnts)[k]))[2];
		}
	}

	double mean_val_y((*pnts[0])[1]);
	double min_val_y((*pnts[0])[1]), max_val_y((*pnts[0])[1]);
	for (size_t k(1); k<n_nodes; k++) {
		mean_val_y += (*pnts[k])[1];
		if ((*pnts[k])[1] < min_val_y) min_val_y = (*pnts[k])[1];
		if (max_val_y < (*pnts[k])[1]) max_val_y = (*pnts[k])[1];
	}
	mean_val_y /= n_nodes;

	double varianz (MathLib::fastpow(mean_val_y-(*pnts[0])[1],2));
	for (size_t k(1); k<n_nodes; k++) {
		varianz += MathLib::fastpow(mean_val_y-(*pnts[k])[1],2);
	}
	std::cout << "statistical data of y coordinate:" << std::endl;
	std::cout << "\tmean value: " << mean_val_y << std::endl;
	std::cout << "\tminimal value: " << min_val_y << std::endl;
	std::cout << "\tmaximal value: " << max_val_y << std::endl;
	std::cout << "\tvarianz: " << varianz << std::endl;
	std::cout << "\tstandard deviation: " << sqrt(varianz) << std::endl << std::endl;

	size_t answer(0);
	while (answer < 1 || answer > 3) {
		std::cout << "Would should I do? Please give the number of the alternative!" << std::endl;
		std::cout << "\t1 Leave the y coordinate for every mesh node untouched." << std::endl;
		std::cout << "\t2 Choose the mean value as y coordinate for every mesh node." << std::endl;
		std::cout << "\t3 Set y = 0 for every mesh node." << std::endl;
		std::cin >> answer;
	}

	if (answer == 1) {
		for (size_t k(0); k<n_nodes; k++) {
			nodes[k]->SetCoordinates(pnts[k]->getData());
		}
	} else {
		if (answer == 2) {
			for (size_t k(0); k<n_nodes; k++) {
				(*pnts[k])[1] = mean_val_y;
				nodes[k]->SetCoordinates(pnts[k]->getData());
			}
			for (size_t k(0); k<n_geo_pnts; k++) {
				(*(*geo_pnts)[k])[1] = mean_val_y;
			}
		} else {
			for (size_t k(0); k<n_nodes; k++) {
				(*pnts[k])[1] = 0.0;
				nodes[k]->SetCoordinates(pnts[k]->getData());
			}
			for (size_t k(0); k<n_geo_pnts; k++) {
				(*(*geo_pnts)[k])[1] = 0.0;
			}
		}
	}

	std::cout << "writing rotated mesh to " << mesh_out_fname << " ... " << std::flush;
	FileIO::OGSMeshIO mesh_io;
	mesh_io.setMesh(&mesh_copy);
	mesh_io.writeToFile (mesh_out_fname);
	std::cout << "done" << std::endl;

	geo_objs->addPointVec(geo_pnts, geo_fname_out);
	geo_objs->getGeometryNames(original_geo_names);

	// copy polylines
	std::vector<GEOLIB::Polyline*> const& original_plys(*(geo_objs->getPolylineVec(original_geo_names[0])));
	std::vector<GEOLIB::Polyline*>* ply_copies(new std::vector<GEOLIB::Polyline*>);
	const size_t n_plys(original_plys.size());
	for (size_t k(0); k<n_plys; k++) {
		GEOLIB::Polyline* ply(new GEOLIB::Polyline(*geo_pnts));
		const size_t n_pnts_in_ply(original_plys[k]->getNumberOfPoints());
		for (size_t j(0); j<n_pnts_in_ply; j++) {
			ply->addPoint(original_plys[k]->getPointID(j));
		}
		ply_copies->push_back(ply);
	}
	geo_objs->addPolylineVec(ply_copies, geo_fname_out);

	// copy surfaces
	std::vector<GEOLIB::Surface*> const*const original_sfcs((geo_objs->getSurfaceVec(original_geo_names[0])));
	if (original_sfcs) {
		std::vector<GEOLIB::Surface*>* sfc_copies(new std::vector<GEOLIB::Surface*>);
		const size_t n_sfcs(original_sfcs->size());
		for (size_t k(0); k<n_sfcs; k++) {
			GEOLIB::Surface* sfc(new GEOLIB::Surface(*geo_pnts));
			const size_t n_tris_in_sfc((*original_sfcs)[k]->getNTriangles());
			for (size_t j(0); j<n_tris_in_sfc; j++) {
				GEOLIB::Triangle const& tri(*(*((*original_sfcs)[k]))[j]);
				sfc->addTriangle(tri[0], tri[1], tri[2]);
			}
			sfc_copies->push_back(sfc);
		}
		geo_objs->addSurfaceVec(sfc_copies, geo_fname_out);
	}

	xml.setPrecision(2);
	xml.setNameForExport(original_geo_names[1]);
	xml.writeToFile(geo_fname_out);

	delete project_data;

	return 0;
}
