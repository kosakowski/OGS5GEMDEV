#ifndef Edge_INC
#define Edge_INC

#include "vec.h"

#include "Grain.h"

//------------------------------------------------------
//   Topology declartion of geometrical element.
//   WW. 06.2005
//   WW. 02.2012
//------------------------------------------------------


namespace Mesh_Group
{
class Node;
class Elem;
//3.  Edge declaration
class Edge:public Grain
{
   public:
      Edge(const int Index, bool quadr=false);
      ~Edge();

      // Access members
      void setNodes( Math_Group::vec<Node*>& Nodes)
      {
         for(int i=0; i<3; i++)  nodes_of_edges[i] = Nodes[i];
      }
      void setNode(const int index,  Node *aNode)
      {
         nodes_of_edges[index] = aNode;
      }
      void getNodes( Math_Group::vec<Node*>& Nodes)
      {
         for(int i=0; i<3; i++)  Nodes[i] = nodes_of_edges[i];
      }
      Node* getNode(const int l_index)
      {
         return nodes_of_edges[l_index];
      }

      void setOrder(const bool order)
      {
         quadratic = order;
      }
 
      // Operator
      void operator = (Edge& edg);
      bool operator == (Edge& edg);

      // Output
      void Write(std::ostream& osm = std::cout) const;
   private:
      // High order
      bool quadratic;

      Math_Group::vec<Node*>  nodes_of_edges;
      friend class Elem;
};

} //end namespace


#endif
