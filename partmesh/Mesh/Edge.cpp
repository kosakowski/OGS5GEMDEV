#include "Edge.h"

#include <iomanip>

#include "Node.h"


//------------------------------------------------------
//   Topology definition of geometrical element.
//    WW. 10.01.2005
//------------------------------------------------------
namespace Mesh_Group
{

using namespace std;

//    WW. 06.2005
//-----------------------------------------------------
//2. Edge
Edge::Edge(const int Index, bool quadr)
   :Grain(Index)
{
   quadratic = quadr;
   index = Index;
   // Assume that each edge has three nodes
   nodes_of_edges.resize(3);
   for(int i=0; i<3; i++)
      nodes_of_edges[i] = NULL;
}
Edge::~Edge()
{
   nodes_of_edges.resize(0);
}

//    WW. 06.2005
void Edge::operator = (Edge& ed)
{
   index = ed.index;
   mark = ed.mark;
   for(int i=0; i<nodes_of_edges.Size(); i++)
      nodes_of_edges[i] = ed.nodes_of_edges[i];
}
//    WW. 06.2005
bool Edge::operator == (Edge& ed)
{
   int identical;

   // Compare two ends
   identical=0;
   for(int i=0; i<2; i++)
   {
      if(nodes_of_edges[i] == ed.nodes_of_edges[i])
         identical++;
   }
   if(identical==2)
      return true;

   identical=0;
   for(int i=0; i<2; i++)
   {
      if(nodes_of_edges[1-i] == ed.nodes_of_edges[i])
         identical++;
   }
   if(identical==2)
      return true;

   return false;
}

//    WW. 06.2005
// Output
void Edge::Write(ostream& osm) const
{
   osm<<"Edge: "<< index<<endl;
   for(int i=0; i<nodes_of_edges.Size(); i++)
   {
      osm<<"Node: "<< i<<endl;
      nodes_of_edges[i]->Write(osm);
   }
   osm<<endl;
}

}// Namespace

