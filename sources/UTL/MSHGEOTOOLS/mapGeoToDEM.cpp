/**
 * \file mapGeoToDEM.cpp
 * 2011/12/19 KR Initial implementation
 */

#include "ProjectData.h"
#include "GEOObjects.h"
#include "FileFinder.h"
#include "VtkRaster.h"
#include "XmlIO/XmlGmlInterface.h"
#include "XmlIO/XmlStnInterface.h"

#include <QApplication>

#include <fstream>
#include <iostream>

float* img_data(NULL);	// pixel information
double origin_x(0), origin_y(0), cellsize(0); // image origin + pixel size
size_t width(0), height(0); // image dimensions

float getElevation(size_t x, size_t y)
{
	if ((x<origin_x) || (x>origin_x+(width*cellsize)) || (y<origin_y) || (y>origin_y+(height*cellsize)))
		return 0;

	size_t x_index = static_cast<size_t>((x-origin_x)/cellsize);
	size_t y_index = static_cast<size_t>((y-origin_y)/cellsize);

	return img_data[2*(y_index*width+x_index)];
}

int mapGeometry(const std::string &geo_name)
{
	std::cout << "Mapping " << geo_name << std::endl;
	ProjectData project;
	GEOLIB::GEOObjects* geo_objects = new GEOLIB::GEOObjects();
	project.setGEOObjects(geo_objects);

	FileFinder fileFinder;
	fileFinder.addDirectory(".");
	fileFinder.addDirectory(std::string(SOURCEPATH).append("/FileIO"));

	FileIO::XmlGmlInterface xml(&project, fileFinder.getPath("OpenGeoSysGLI.xsd"));
	if (xml.readFile(QString::fromStdString(geo_name)))
	{
		std::vector<std::string> names;
		geo_objects->getGeometryNames(names);
		std::vector<GEOLIB::Point*> *points = const_cast<std::vector<GEOLIB::Point*>*>(geo_objects->getPointVec(names[0]));
		size_t nPoints (points->size());
		for (size_t j=0; j<nPoints; j++)
		{
			GEOLIB::Point* pnt = (*points)[j];
			(*pnt)[2] = getElevation((*pnt)[0],(*pnt)[1]);
		}

		std::string new_geo_name = geo_name.substr(0, geo_name.length()-4) + "_elevation.gml";
		xml.setNameForExport(names[0]);
		xml.writeToFile(new_geo_name);

		std::cout << "New file \"" << new_geo_name << " successfully written." << std::endl;
		std::cout << std::endl;

		return 1;
	}

	std::cout << "Error: Could not open geometry file..." << std::endl;
	return -1;
}

int mapStations(const std::string &geo_name)
{
	std::cout << "Mapping " << geo_name << std::endl;
	ProjectData project;
	GEOLIB::GEOObjects* geo_objects = new GEOLIB::GEOObjects();
	project.setGEOObjects(geo_objects);

	FileFinder fileFinder;
	fileFinder.addDirectory(".");
	fileFinder.addDirectory(std::string(SOURCEPATH).append("/FileIO"));

	FileIO::XmlStnInterface xml(&project, fileFinder.getPath("OpenGeoSysSTN.xsd"));
	if (xml.readFile(QString::fromStdString(geo_name)))
	{
		bool is_borehole(false);
		std::vector<std::string> names;
		geo_objects->getStationVectorNames(names);
		std::vector<GEOLIB::Point*> *points = const_cast<std::vector<GEOLIB::Point*>*>(geo_objects->getStationVec(names[0]));
		if (static_cast<GEOLIB::StationBorehole*>((*points)[0])->type() == GEOLIB::Station::BOREHOLE)
			is_borehole = true;
		size_t nPoints (points->size());
		for (size_t j=0; j<nPoints; j++)
		{
			GEOLIB::Point* pnt = (*points)[j];

			if (!is_borehole)
				(*pnt)[2] = getElevation((*pnt)[0],(*pnt)[1]);
			else
			{
				double offset(getElevation((*pnt)[0],(*pnt)[1]) - (*pnt)[2]);
				GEOLIB::StationBorehole* borehole = static_cast<GEOLIB::StationBorehole*>(pnt);
				const std::vector<GEOLIB::Point*> layers = borehole->getProfile();
				size_t nLayers = layers.size();
				for (size_t k=0; k<nLayers; k++)
				{
					GEOLIB::Point* layer_pnt = layers[k];
					(*layer_pnt)[2] = (*layer_pnt)[2] + offset;
				}
			}
		}

		std::string new_geo_name = geo_name.substr(0, geo_name.length()-4) + "_elevation.stn";
		xml.setNameForExport(names[0]);
		xml.writeToFile(new_geo_name);

		std::cout << "New file \"" << new_geo_name << " successfully written." << std::endl;
		std::cout << std::endl;

		return 1;
	}

	std::cout << "Error: Could not open geometry file..." << std::endl;
	return -1;
}

int main (int argc, char* argv[])
{
	QApplication app(argc, argv, false);

	if (argc != 3)
	{
		std::cout << "Changes the z-Coordinates of the geometric objects in the gml- or stn-files according to a given DEM." << std::endl;
		std::cout << std::endl;
		std::cout << "Usage: " << argv[0] << " <geo-file.gml> <DEM-file.asc>" << std::endl;
		return -1;
	}

	bool isList(false);
	std::string geo_name = argv[1];
	std::string dem_name = argv[2];
	std::string gml_name("");

	if ((geo_name.substr(geo_name.length()-4, 4).compare(".gml") != 0) &&
		(geo_name.substr(geo_name.length()-4, 4).compare(".stn") != 0) &&
		(geo_name.substr(geo_name.length()-4, 4).compare(".lst") != 0))
	{
		std::cout << "Error: Parameter 1 should be a gml- or stn-file" << std::endl;
		std::cout << "Usage: " << argv[0] << " <geo-file.gml> <DEM-file.asc>" << std::endl;
		std::cout << std::endl;
		return -1;
	}

	if (dem_name.substr(dem_name.length()-4, 4).compare(".asc") != 0)
	{
		std::cout << "Error: Parameter 2 should be an asc-file" << std::endl;
		std::cout << "Usage: " << argv[0] << " <geo-file.gml> <DEM-file.asc>" << std::endl;
		return -1;
	}

	if (geo_name.substr(geo_name.length()-4, 4).compare(".lst") == 0)
		isList = true;

	img_data = VtkRaster::loadDataFromASC(dem_name, origin_x, origin_y, width, height, cellsize);

	if (img_data != NULL)
	{
		// map list of geometries to the same DEM
		if (isList)
		{
			std::ifstream in(geo_name.c_str());
			while (!in.eof())
			{
				in >> gml_name;

				if (gml_name.substr(gml_name.length()-4, 4).compare(".gml") == 0)
					mapGeometry(gml_name);
				else if (gml_name.substr(gml_name.length()-4, 4).compare(".stn") == 0)
					mapStations(gml_name);
				else
					std::cout << "File extension for " << gml_name << " unknown." << std::endl;
			}
			return 1;
		}
		// map only one geometry
		else
		{
			if (geo_name.substr(geo_name.length()-4, 4).compare(".gml") == 0)
				mapGeometry(geo_name);
			else if (geo_name.substr(geo_name.length()-4, 4).compare(".stn") == 0)
				mapStations(geo_name);
		}
	}
	else
	{
		std::cout << "Error: Could not read DEM file..." << std::endl;
		return -1;
	}
	return -1;
}



