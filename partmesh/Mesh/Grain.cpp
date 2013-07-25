
#include "Grain.h"

//------------------------------------------------------
//   Topology definition of geometrical element.
//    WW. 10.01.2005
//    WW 03.2012
//------------------------------------------------------
namespace Mesh_Group
{

//-----------------------------------------------------
//    WW. 06.2005
//0. Base class
Grain::Grain(const int id)
{
   index = id;
   mark = false;
}

}

