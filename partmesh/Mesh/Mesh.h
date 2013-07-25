#ifndef Mesh_INC
#define Mesh_INC

#include<string>
#include<vector>
#include<iostream>

//------------------------------------------------------
//   Topology declartion of geometrical element.
//   WW. 06.2005
//------------------------------------------------------

namespace Mesh_Group
{
class Node;
class Edge;
class Elem;


/*!
   \class Mesh
*/
class Mesh
{
   public:
      Mesh(bool quad = false);
      ~Mesh();


      void ReadGrid(std::istream& is = std::cin);
      void ReadGridGeoSys(std::istream& is = std::cin);

      void Write2METIS(std::ostream& os);
      void WriteVTK_Nodes(std::ostream& os);
      void WriteVTK_Nodes(std::ostream& os, std::vector<Node*>& nod_vec, const size_t start = 0);
      void WriteVTK_Elements_of_Subdomain(std::ostream& os, std::vector<Elem*>& ele_vec,
                                          const int sbd_index, const long node_shift = 0);

      void ConstructSubDomain_by_Elements(const std::string fname,  const int num_parts, const bool osdom);
      void ConstructSubDomain_by_Nodes(const std::string fname, const std::string fpath, const std::string mat_fname,
		                               const int num_parts, const bool is_quad, const bool osdom);

      void ConnectedNodes(bool quadratic);
      void ConnectedElements2Node(bool quadratic=false);

      void ConstructGrid();
      void GenerateHighOrderNodes();

      void setOrder(const bool is_quad)
      {
         useQuadratic = is_quad;
      }
   private:
      // The following can be members of grid class
      long NodesNumber_Linear;
      long NodesNumber_Quadratic;
      bool useQuadratic;
      bool axisymmetry;

      // Coordinate indicator
      // 1:  X component only
      // 12: Y component only
      // 13: Z component only
      // 2:  X, Y component
      // 23:  X, Z component
      // 3:  X, Y, Z component
      int coordinate_system;
      int max_ele_dim;

      // All nodes
      std::vector<Node*> node_vector;
      // All edges
#ifdef BUILD_MESH_EDGE
      std::vector<Edge*> edge_vector;
#endif
      // All surface feces
#ifdef BUILD_MESH_FACE
      std::vector<Elem*> face_vector;
#endif
      // All elements
      std::vector<Elem*> elem_vector;

      long msh_no_line;
      long msh_no_quad;
      long msh_no_hexs;
      long msh_no_tris;
      long msh_no_tets;
      long msh_no_pris;
      long msh_no_pyra;
      int msh_max_dim;

	  friend class Elem;
};

}

#endif