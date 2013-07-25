#ifndef node_INC
#define node_INC

#include<vector>

#include "Grain.h"

//------------------------------------------------------
//   Topology declartion of geometrical element.
//   WW. 06.2005
//   WW. 02.2012
//------------------------------------------------------
namespace Mesh_Group
{
class Edge;
class Elem;
class Mesh;

//2.  Node declaration
class Node:public Grain
{
   public:
      Node(const int Index):Grain(Index)
      {
         Coordinate = new double[3];
      }
      Node(const int Index, const double x,
           const double y, const double z=0.0);
      ~Node()
      {
         delete [] Coordinate;
         Coordinate = NULL;
         ElementsRelated.resize(0);
      }

      // Operator
      void operator = (const Node& n);
      bool operator == (const Node & n);

      // Change members;
      // By component
      void setX(const double argX)
      {
         Coordinate[0] = argX;
      }
      void setY(const double argY)
      {
         Coordinate[1] = argY;
      }
      void setZ(const double argZ)
      {
         Coordinate[2] = argZ;
      }
      void SetCoordinates(const double* argCoord);

      // Access to members
      double X() const
      {
         return Coordinate[0];
      }
      double Y() const
      {
         return Coordinate[1];
      }
      double Z() const
      {
         return Coordinate[2];
      }
      double *getCoordinates() const
      {
         return  Coordinate;
      }

      void setLocalIndex(const long l_index)
      {
         local_index = l_index;
      }
      long getLocalIndex() const
      {
         return local_index;
      }

      // Output
      void Write(std::ostream& os = std::cout) const;
      void WriteCoordinates(std::ostream& os = std::cout) const;

   private:
      double *Coordinate;
      long local_index; // For domain decomposition
      std::vector<long>  ElementsRelated;
      std::vector<long>  NodesRelated;
      friend class Mesh_Group::Edge;
      friend class Mesh_Group::Elem;
      friend class Mesh_Group::Mesh;

};

}

#endif
