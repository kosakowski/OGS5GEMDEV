/*
 * createLayeredMesh.cpp
 *
 *  Created on: Apr 27, 2012
 *      Author: TF
 */

#include <ctime>

// FileIO
#include "MeshIO/OGSMeshIO.h"

// MSH
#include "msh_lib.h" // for FEMRead

// we need this for using the xml functions of Qt
#include <QApplication>

#include "Qt/DataView/MshLayerMapper.h"

int main (int argc, char* argv[])
{
	// Creating a non-gui (console) Qt application
	QApplication app(argc, argv, false);

	if (argc < 3) {
		std::cout << "Usage: " << argv[0] << " --mesh-in meshfile "
						<< std::flush;
		std::cout << "--number-of-layers number " << std::flush;
		std::cout << std::endl;
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

	// *** parsing number of search points
	tmp = argv[3];
	size_t n(0);
	if (tmp.find("--number-of-layers") == std::string::npos) {
		std::cout << "could not find switch for reading number of layers " << std::endl;
		return -1;
	}
	n = static_cast<size_t>(atoi(argv[4]));

	// *** read mesh
	std::vector<MeshLib::CFEMesh*> mesh_vec;
	FEMRead(file_base_name, mesh_vec);
	if (mesh_vec.empty()) {
		std::cerr << "could not read mesh from file " << std::endl;
		return -1;
	}
	MeshLib::CFEMesh* mesh (mesh_vec[mesh_vec.size() - 1]);
	mesh->setNumberOfNodesFromNodesVectorSize();
	mesh->ConstructGrid();

	double thickness (10);
	MeshLib::CFEMesh* layered_mesh (MshLayerMapper::CreateLayers(mesh, n, thickness));

	std::string mesh_out_fname("LayeredMesh.msh");
	FileIO::OGSMeshIO mesh_io;
	mesh_io.setMesh(layered_mesh);
	mesh_io.writeToFile (mesh_out_fname);

	return 0;
}
