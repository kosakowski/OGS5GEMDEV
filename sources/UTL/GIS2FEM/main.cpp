// expre_new_Operator.cpp
// compile with: /EHsc
#include <cmath>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>

#include <float.h>
#include <iomanip>
#include <string>
#include <vector>

//#include <time.h>

#include "misc.h"
#include "msh_mesh.h"

std::string FileName;
std::string FilePath; //WW

using namespace std;
int main(int argc, char* argv[])
{
	const int max_size = 1028;
	char str1[max_size];
	char str2[max_size];

	int option = 1;
	CFEMesh* a_mesh = NULL;
	fstream ofile;
	string ofname;

	cout << "\t|================================================| " << endl;
	cout << "\t|                                                | " << endl;
	cout << "\t|        Toolkit for hydraulic modeling          | " << endl;
	cout << "\t|                                                | " << endl;
	cout << "\t|                 By    WW                       | " << endl;
	cout << "\t|                                                | " << endl;
	cout << "\t|    Argument:                                   | " << endl;
	cout << "\t|      option and file name (no extension)       | " << endl;
	cout << "\t|    Option:                                     | " << endl;
	cout << "\t|     1. Generation OGS Neumman BC from          | " << endl;
	cout << "\t|        raster files of GIS                     | " << endl;
	cout << "\t|     2. Top surface integration for 3D mesh     | " << endl;
	cout << "\t|     3. Convert GIS raster cells into FE mesh   | " << endl;
	cout << "\t|                                                | " << endl;
	cout << "\t|================================================| " << endl;
	cout << "\tInput file name: ";

	if(argc > 1)
	{
		strcpy(str2,argv[1]);
		strcpy(str1,argv[2]);
	}
	else
		scanf(" %s %s%*[^\n]%*c",str2, str1);

	sscanf(str2, "%d", &option);

	FileName = str1;
	basic_string <char>::size_type indexChWin, indexChLinux;
	indexChWin = indexChLinux = 0;
	indexChWin = FileName.find_last_of('\\');
	indexChLinux = FileName.find_last_of('/');
	//
	if(indexChWin != string::npos)
		FilePath = FileName.substr(0,indexChWin) + "\\";
	else if(indexChLinux != string::npos)
		FilePath = FileName.substr(0,indexChLinux) + "/";

	if(option != 3)
	{
		FEMRead(FileName);
		a_mesh = fem_msh_vector[0];
	}

	switch(option)
	{
	case 1:
		a_mesh->mHM2NeumannBC();
		break;
	case 2:
		a_mesh->TopSurfaceIntegration();
		break;
	case 3:
		a_mesh = new CFEMesh();
		a_mesh->ConvertShapeCells(FileName + ".asc");
		ofname = FileName + ".msh";
		ofile.open(ofname.c_str(), ios::out | ios::trunc);
		a_mesh->Write(&ofile);
		break;

	default:

		break;
	}
	return 0;
}
