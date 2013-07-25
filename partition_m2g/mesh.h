#ifndef rf_msh_fem_INC
#define rf_msh_fem_INC

#include<string>
#include<vector>
#include<iostream>
#include "matrix_class.h"

//------------------------------------------------------ 
//   Topology declartion of geometrical element.
//   WW. 06.2005
//------------------------------------------------------ 
using namespace std;

namespace Mesh_Group
{
 enum BCType  { INTERIOR, DIRCHLET, NEUMANN, CAUCHY, BOUNDARY };
 //enum Element_Type  {Line, Quad, Hex, Tri, Tet, Pris };
 using Math_Group::SymMatrix;
 using Math_Group::vec;

//1.  Mesh declaration
class Grain
{
   public:
      Grain(const int id);
      virtual  ~Grain() {}
      // Operator
      virtual void operator = (const Grain & g) {}
      virtual bool operator == (const Grain & g) {return false;}

      // Set members
      void secBC(const char BC_type) {boundayC = BC_type;}
      void SetOrder(const bool order) {quadratic = order;}
      void Marking(const bool state) {mark = state;}
      // Get members
      long GetIndex() const {return index;} 
      bool GetStatus() const {return mark;}
      bool Dirichlet() const  { return (BCType(boundayC) == DIRCHLET); }
      bool Neumann()   const  { return (BCType(boundayC) == NEUMANN);   }
      bool Cauchy ()   const  { return (BCType(boundayC) == CAUCHY);    }
      bool onBoundary() const { return (BCType(boundayC) == BOUNDARY);  }
      bool Interior() const { return (BCType(boundayC) == INTERIOR);  }

      // Output
      virtual void output(ostream& os=cout) const {};
   protected:
      long index;
      char boundayC;
      // Towards special purpose, 
      // e.g. marked to be refined or active 
      bool mark;
      // High order
      bool quadratic;
    
      // delimitor
      string deli;

};
    
//2.  Node declaration
class Node:public Grain
{
   public:
						Node(const int Index):Grain(Index) {}
      Node(const int Index, const double x, 
           const double y, const double z=0.0);
      ~Node() {ElementsBelonged.resize(0);}

      // Operator
      void operator = (const Node& n);
      bool operator == (const Node & n);

      // Change members;
      // By component
      void setX(const double argX)  { Coordinate[0] = argX;}
      void setY(const double argY)  { Coordinate[1] = argY;}
      void setZ(const double argZ)  { Coordinate[2] = argZ;}
      void SetCoordinates(const double* argCoord); 

      // Access to members
      double X() const {return Coordinate[0];}
      double Y() const {return Coordinate[1];}
      double Z() const {return Coordinate[2];}
      void Coordinates(double *xyz) const
      {  for(int i=0; i<3; i++)  xyz[i] = Coordinate[i];} 
	  void SetLocalIndex(const long l_index) {local_index = l_index; } 
	  long GetLocalIndex() const {return local_index;}
     
      // Output
      void Write(ostream& os=cout) const;

      vector<long>  ElementsBelonged;

   private:
      double Coordinate[3];
      long local_index; // For domain decomposition
      friend class Edge;
      friend class Elem;
};

//3.  Edge declaration
class Edge:public Grain
{
   public:
      Edge(const int Index, bool quadr=false);
      ~Edge(); 

      // Access members
      void setNodes( vec<Node*>& Nodes)
        { for(int i=0; i<3; i++)  nodes_of_edges[i] = Nodes[i]; }
      void getNodes( vec<Node*>& Nodes) 
        { for(int i=0; i<3; i++)  Nodes[i] = nodes_of_edges[i]; }										
      // Operator
      void operator = (Edge& edg);
      bool operator == (Edge& edg);
      // Output
      void Write(ostream& osm=cout) const;
   private:
      vec<Node*>  nodes_of_edges;
      friend class Elem;
};

//3.  Element declaration
class Elem:public Grain
{
   public:
      Elem(const int Index);
      ~Elem();

      // Operator
      // virtual void operator = (const Elem& elem);

      // Access to members
      int getNodesNumber() const {return nnodes;}
      int getNodesNumberHQ() const {return nnodesHQ;}
      int getEdgesNumber() const{return nedges;}
      int getFacesNumber() const {return nfaces;}
      int getElementType() const {return ele_Type;}
      int getPatchIndex() const {return PatchIndex;}
      int Dim() const {return ele_dim;}
      string getName() const;
      void setJacobian(const SymMatrix& Jac) {Jacobian = Jac;}
      void setVolume(const double Vol) {Volume = Vol;}
      double getVolume() const {return Volume;}

   	  // Nodes
      void getNodeIndeces(vec<long>&  node_index) const 
	            {for (int i=0; i< (int) nodes_index.Size();i++)
                              node_index[i]= nodes_index[i]; } 
      long getNodeIndex(const int loc_lndex) const {  return nodes_index[loc_lndex];}
      void setNodes(vec<Node*>&  ele_nodes, const bool ReSize=false);
      void getNodes(vec<Node*>&  ele_nodes) 
         { for (int i=0; i< (int) nodes.Size();i++) ele_nodes[i]= nodes[i]; }
      void MarkingNodes(bool merker) 
	    { for (int i=0; i< (int) nodes.Size();i++) nodes[i]->Marking(merker); }
      void SetLocalNodeIndex(const int li, const long n_lindex) 
	     {  nodes[li]->local_index = n_lindex; }
      long GetLocalNodeIndex(const int li) const 
	     {  return nodes[li]->local_index; }
         		   
     // Edges
      void getEdges(vec<Edge*>&  ele_edges) 
        {for (int i=0; i<nedges; i++) ele_edges[i]= edges[i];} 
      void setEdges(vec<Edge*>&  ele_edges) 
        {for (int i=0; i<nedges; i++) edges[i]= ele_edges[i];} 
						
	  void setEdges_Orientation(vec<int>&  ori_edg) 
         {for (int i=0; i<nedges; i++) edges_orientation[i]= ori_edg[i];} 
      // Neighbors
      void setNeighbors(vec<Elem*>&  ele_neighbors) 
         { for (int i=0; i< nfaces;i++) neighbors[i] = ele_neighbors[i];}
      void setNeighbor(const int LocalIndex, Elem* ele_neighbor) 
         { neighbors[LocalIndex] = ele_neighbor;}
      void getNeighbors(vec<Elem*>&  ele_neighbors)  
         {for (int i=0; i< nfaces;i++) ele_neighbors[i]= neighbors[i];} 
      //Domain partition
      long GetDomNodeIndex(const int loc_index) { return locnodes_index[loc_index];}
	  void setDomNodeIndex(const int loc_index, const long dom_nindex) { locnodes_index[loc_index]=dom_nindex;}
      void AllocateLocalIndexVector() {locnodes_index.resize(nodes_index.Size());}
      void setDomainIndex(const int dom) {sub_dom = dom;} 
      int GetDomainIndex() const { return sub_dom;} 
      // Local indicis
      void GetLocalIndices_EdgeNodes(const int Edge, int *EdgeNodes);
      int GetElementFaceNodes(const int Face, int *FacesNode);

      // Faces: Face, local face index
      Elem( const int Index, Elem* onwer, const int Face);
      int GetFaceType();
      
      // Output
      void Read(istream& is=cin, int fileType=1);
      void WriteIndex(ostream& os=cout) const;
      void Write_index(ostream& os=cout) const;
      void WriteAll(ostream& os=cout) const;
      void WriteNeighbors(ostream& os=cout) const;
   private:
      int nnodes;
      int nnodesHQ;
      int ele_dim;         // Dimension of element
      int nfaces;
      int nedges;
      int sub_dom;
	  //
      double Volume;
      SymMatrix Jacobian;     
      int PatchIndex;
	  // Element type
      // 1 Line, 2 Quad, 3 Hex, 4 Tri, 5 Tet, 6 Pris 
      int ele_Type;
      Elem* Owner;
      vec<long>   nodes_index;
      vec<long>   locnodes_index;
      vec<Node*>  nodes;
      vec<Edge*>  edges;
      vec<int>  edges_orientation;
      vec<Elem*>  neighbors;
 
      // Private methods
      int GetElementFaces1D(int *FaceNode);
      int GetElementFacesTri(const int Face, int *FaceNode);
      int GetElementFacesQuad(const int Face, int *FaceNode);
      int GetElementFacesHex(const int Face, int *FaceNode);
      int GetElementFacesTet(const int Face, int *FaceNode);
      int GetElementFacesPri(const int Face, int *FaceNode);

};

}; //end namespace


// The following can be members of grid class
extern long NodesNumber_Linear;
extern long NodesNumber_Quadratic;
// All nodes
extern vector<Mesh_Group::Node*> NodesVector;
// All edges
extern vector<Mesh_Group::Edge*> EdgeVector;
// All surface feces
extern vector<Mesh_Group::Elem*> SurfaceFaces;
// All elements 
extern vector<Mesh_Group::Elem*> ElementsVector;

extern void ConstructGrid( const bool quadratic=false);
extern void ConstructDomain(char *fname,  const int num_parts);
// Test read
extern void ReadGrid(istream& is=cin);
extern void ReadGridGeoSys(istream& is=cin);
// Test dextructor
extern void Realese();


#endif
