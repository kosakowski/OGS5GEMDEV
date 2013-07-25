
#include "misc.h"

#include "msh_mesh.h"
vector<MeshLib::CFEMesh*> fem_msh_vector;

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2010 WW  Based on FEMRead in MSHLib
**************************************************************************/
bool FEMRead(string file_base_name)
{
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	CFEMesh* m_fem_msh = NULL;
	char line[1024];
	string sub_line;
	string line_string;
	ios::pos_type position;
	//========================================================================
	// File handling
	bool msh_file_binary = false;
	string msh_file_name_bin = file_base_name + "_binary" + ".msh";
	string msh_file_name_ascii = file_base_name + ".msh";
	ifstream msh_file_bin;
	ifstream msh_file_ascii;

	msh_file_bin.open(msh_file_name_bin.c_str(),ios::binary | ios::in);
	if(msh_file_bin.good())
		msh_file_binary = true;

	//----------------------------------------------------------------------
	cout << "MSHRead: ";

	cout << "ASCII file" << endl;
	msh_file_ascii.open(msh_file_name_ascii.data(),ios::in);
	if (!msh_file_ascii.good())
	{
		cout << "Opening MSH file failed" << endl;
		return false;
	}
	//----------------------------------------------------------------------
	//========================================================================
	// Keyword loop
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	while(!msh_file_ascii.eof())
	{
		msh_file_ascii.getline(line,1024);
		line_string = line;
		if(line_string.find("#STOP") != string::npos)
			return true;
		//..................................................................
		if(line_string.find("#FEM_MSH") != string::npos) // keyword found
		{
			m_fem_msh = new CFEMesh();
			position = m_fem_msh->Read(&msh_file_ascii);
			fem_msh_vector.push_back(m_fem_msh);
			msh_file_ascii.seekg(position,ios::beg);
		} // keyword found
	} // eof

	msh_file_ascii.close();

	//========================================================================
	return true;
}
