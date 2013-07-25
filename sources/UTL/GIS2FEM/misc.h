#ifndef misc_INC
#define misc_INC
#include <vector>
using namespace std;

namespace  MeshLib
{class CFEMesh;
}
using  MeshLib::CFEMesh;
extern vector<MeshLib::CFEMesh*> fem_msh_vector;

extern bool FEMRead(string);
#endif