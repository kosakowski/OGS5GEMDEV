#ifndef Elem_INC
#define Elem_INC

#include<vector>

#include"vec.h"
#include"SymMatrix.h"
#include "Grain.h"
#include "Node.h"


//#define BUILD_MESH_EDGE

//------------------------------------------------------
//   Topology declartion of geometrical element.
//   WW. 06.2005
//   WW. 02.2012
//------------------------------------------------------

namespace Mesh_Group
{
class Node;
class Edge;
class Mesh;
enum ElemType {line, quadri, hex, tri, tet, prism, pyramid}; 


//3.  Element declaration
class Elem:public Grain
{
   public:
      Elem(const int Index);
      // Faces: Face, local face index
      Elem( const int Index, Elem* onwer, const int Face);

      ~Elem();

	  void Init();
      // Operator
      // virtual void operator = (const Elem& elem);

      // Access to members
      int getNodesNumber() const
      {
         return nnodes;
      }
      int getNodesNumber(bool quad) const
      {
         if(quad) return nnodesHQ;
         else return nnodes;
      }
      int getNodesNumberHQ() const
      {
         return nnodesHQ;
      }
      int getEdgesNumber() const
      {
         return nedges;
      }
      int getFacesNumber() const
      {
         return nfaces;
      }
      int getElementType() const
      {
         return ele_Type;
      }
      int getPatchIndex() const
      {
         return PatchIndex;
      }
      int Dim() const
      {
         return ele_dim;
      }
      std::string getName() const;

      // Nodes
      void getNodeIndeces(long  *node_index) const
      {
         for (int i=0; i< (int) nodes.Size(); i++)
			 node_index[i]= nodes[i]->index;
      }
      long getNodeIndex(const int loc_lndex) const
      {
		  return nodes[loc_lndex]->index;
      }
      void setNodes(Math_Group::vec<Node*>&  ele_nodes, const bool ReSize=false);
      void getNodes(Math_Group::vec<Node*>&  ele_nodes)
      {
         for (int i=0; i< (int) nodes.Size(); i++)
            ele_nodes[i]= nodes[i];
      }
      Node* getNode(const int i)
      {
         return nodes[i];
      }
      void MarkingNodes(bool maker);
      //
      void setLocalNodeIndex(const int li, const long n_lindex);
      //
      long getLocalNodeIndex(const int li) const;

#ifdef BUILD_MESH_EDGE
      // Edges
      void getEdges(Math_Group::vec<Edge*>&  ele_edges)
      {
         for (int i=0; i<nedges; i++) ele_edges[i]= edges[i];
      }
      Edge* getEdge(const int index)
      {
         return edges[index];
      }
      void setEdges(Math_Group::vec<Edge*>&  ele_edges)
      {
         for (int i=0; i<nedges; i++) edges[i]= ele_edges[i];
      }

      void setEdges_Orientation(Math_Group::vec<int>&  ori_edg)
      {
         for (int i=0; i<nedges; i++) edges_orientation[i]= ori_edg[i];
      }
#endif
      // Neighbors
      void setNeighbors(Math_Group::vec<Elem*>&  ele_neighbors)
      {
         for (int i=0; i< nfaces; i++) neighbors[i] = ele_neighbors[i];
      }
      void setNeighbor(const int LocalIndex, Elem* ele_neighbor)
      {
         neighbors[LocalIndex] = ele_neighbor;
      }
      void getNeighbors(Math_Group::vec<Elem*>&  ele_neighbors)
      {
         for (int i=0; i< nfaces; i++) ele_neighbors[i]= neighbors[i];
      }
      //Domain partition
      long getDomNodeIndex(const int loc_index)
      {
         return locnodes_index[loc_index];
      }
      void setDomNodeIndex(const int loc_index, const long dom_nindex)
      {
         locnodes_index[loc_index] = dom_nindex;
      }
      void AllocateLocalIndexVector()
      {
         locnodes_index.resize(nodes.Size());
      }
      void setDomainIndex(const int dom)
      {
         sub_dom = dom;
      }
      int getDomainIndex() const
      {
         return sub_dom;
      }
      // Local indicis
      void getLocalIndices_EdgeNodes(const int Edge, int *EdgeNodes);
      int getElementFaceNodes(const int Face, int *FacesNode);


      int getFaceType();

      void setOrder(const bool order)
      {
         quadratic = order;
      }

      // Output
	  void Read(std::istream& is, Mesh_Group::Mesh *mesh, int fileType);
      void WriteIndex(std::ostream& os = std::cout) const;
      void WriteGmsh(std::ostream& os, const int sdom_idx = 0) const;
      void WriteGSmsh(std::ostream& os, bool quad = false) const;
      void WriteSubDOM(std::ostream& os, const long node_id_shift, bool quad = false) const;
      void WriteVTK_Type(std::ostream& os, bool isquad) const;
      void Write_index(std::ostream& os = std::cout) const;
      void WriteAll(std::ostream& os = std::cout) const;
      void WriteNeighbors(std::ostream& os = std::cout) const;
   private:
      // High order
      bool quadratic;

      int nnodes;
      int nnodesHQ;
      int ele_dim;         // Dimension of element
      int nfaces;
      int nedges;
      int sub_dom;
      int no_faces_on_surface;
      //

	  int PatchIndex;
      // Element type
      // 1 Line, 2 Quad, 3 Hex, 4 Tri, 5 Tet, 6 Pris
      ElemType ele_Type;
      Elem* Owner;
      //Math_Group::vec<long>   nodes_index;
      Math_Group::vec<long>   locnodes_index;
      Math_Group::vec<Node*>  nodes;

#ifdef BUILD_MESH_EDGE
      Math_Group::vec<Edge*>  edges;
      Math_Group::vec<int>  edges_orientation;
#endif
      Math_Group::vec<Elem*>  neighbors;
      //vec<Elem*>  sons;
  
      int nnodes_gl; //> number of ghost nodes for linear element 
      std::vector<int>  ghost_nodes;

      // Private methods
      int getElementFaces1D(int *FaceNode);
      int getElementFacesTri(const int Face, int *FaceNode);
      int getElementFacesQuad(const int Face, int *FaceNode);
      int getElementFacesHex(const int Face, int *FaceNode);
      int getElementFacesTet(const int Face, int *FaceNode);
      int getElementFacesPri(const int Face, int *FaceNode);
      int getElementFacesPyramid(const int Face, int *FaceNode);
      friend class Mesh;

};

} //end namespace

#endif
