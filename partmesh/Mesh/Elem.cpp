#include "Elem.h"

#include <iomanip>

#include "Node.h"
#include "Edge.h"
#include "Mesh.h"

//------------------------------------------------------
//   Topology definition of geometrical element.
//    WW. 10.01.2005
//------------------------------------------------------
namespace Mesh_Group
{

using namespace std;
using namespace Math_Group;


int EdgeLocalNodeIndex [] =
{
    // Line, 3 entries
	0, 1, 2,
	// quadri, 4 edges, 12 entries. sh: 3
	0, 1, 4, 
    1, 2, 5,
	2, 3, 6,
	3, 0, 7, 
	// hex, 12 edges, 36 entries. sh: 15
	0, 1, 8, 
	1, 2, 9,
	2, 3, 10,
	3, 0, 11, 
	4, 5, 12, 
	5, 6, 13,
	6, 7, 14,
	7, 4, 15, 
	4, 0, 16, 
	5, 1, 17,
	6, 2, 18,
	7, 3, 19,
	//tri, 3 edges, 9 entries. sh: 51     
	0, 1, 3, 
    1, 2, 4,
	2, 0, 5,
	//tet, 6 edges, 18 entries. sh: 60
	0, 1, 4, 
    1, 2, 5,
	2, 0, 6,
	0, 3, 7, 
    1, 3, 8,
	2, 3, 9,
	//prism, 9 edges, 27 entries. sh: 78 
	0, 1, 6, 
    1, 2, 7,
	2, 0, 8,
	3, 4, 9, 
    4, 5, 10,
	5, 3, 11,
	0, 3, 12,
	1, 4, 13,
	2, 5, 15,
	//pyramid, 8 edges, 24 entries. sh: 105
	0, 1, 5,
	1, 2, 6,
	2, 3, 7,
	3, 0, 8,
	4, 0, 9,
	4, 1, 10,
	4, 2, 11,
	4, 3, 12
};
int EdgeLocalIndexArrayElemShift [] =
{  
	0, //Line
	3, //quad
	15, //hex
	51, //tri
	60, //tet
	78, //prism
	106 //pyramid
};


//-----------------------------------------------------
//2. Mesh
//    WW. 06.2005
Elem::Elem(const int Index):Grain(Index)
{

   nnodes = 0;
   nnodesHQ = 0;
   ele_dim = 1;         // Dimension of element
   PatchIndex = 0;
   //
   quadratic = false;

   Owner = NULL;
}
// Face element
//    WW. 06.2005
Elem::  Elem( const int Index,  Elem* onwer, const int Face):
   Grain(Index), Owner(onwer)
{
   int i, n;
   static int faceIndex_loc[10];
   static int edgeIndex_loc[10];
//  Owner = onwer;
   n = Owner->getElementFaceNodes(Face, faceIndex_loc);
   switch(Owner->ele_Type)
   {
      case line:  // 1-D bar element
         break;
      case quadri: // 2-D quadrilateral element
         ele_Type = line;
         break;
      case hex: // 3-D hexahedral element
         ele_Type = quadri;
         break;
      case tri:  // 2-D triagular element
         ele_Type = line;
         break;
      case tet:  // 3-D tetrahedral element
         ele_Type = tri;
         break;
      case prism:
         if(Face<2)
         {
            ele_Type = tri;
         }
         else
         {
            ele_Type = quadri;
         }
         break; // 3-D prismatic element
      case pyramid:
         if(Face == 0) // Bottom
         {
            ele_Type = quadri;
         }
         else 
         {
            ele_Type = tri;
         }
         break; // 3-D prismatic element
	  default:
		 break;
   }

   Init();
   PatchIndex =  Owner->PatchIndex;
   quadratic = Owner->quadratic;
   nodes.resize(n);

   for(i=0; i<n; i++)
   {
      nodes[i] = Owner->nodes[faceIndex_loc[i]];
   }

#ifdef BUILD_MESH_EDGE
   int j, k, ne;

   // Face edges
   ne = Owner->getEdgesNumber();
   edges.resize(nnodes);
   edges_orientation.resize(nnodes);
   edges_orientation = 1;
   for(i=0; i<nnodes; i++)
   {
      k = (i+1)%nnodes;
      for(j=0; j<ne; j++)
      {
         Owner->getLocalIndices_EdgeNodes(j, edgeIndex_loc);
         if( (faceIndex_loc[i]==edgeIndex_loc[0]&&
               faceIndex_loc[k]==edgeIndex_loc[1])||
               (faceIndex_loc[i]==edgeIndex_loc[1]&&
                faceIndex_loc[k]==edgeIndex_loc[0]) )
         {
            edges[i] = Owner->edges[j];
            if(faceIndex_loc[i]==edgeIndex_loc[1]&&
                  faceIndex_loc[k]==edgeIndex_loc[0] )
               edges_orientation[i] = -1;
            break;
         }
      }
   }
#endif
}

//    WW. 06.2005
Elem::~Elem()
{
   locnodes_index.resize(0);
   nodes.resize(0);
#ifdef BUILD_MESH_EDGE
   edges.resize(0);
#endif
   neighbors.resize(0);
   ghost_nodes.resize(0);
}

void Elem::Init()
{
   // 2 Element configuration
   switch(ele_Type)
   {
      case line:
         nnodes = 2;
         nnodesHQ = 3;
         ele_dim = 1;
         nfaces = 2;
         nedges = 0;
         break;
      case quadri:
         nnodes = 4;
         nnodesHQ = 9;
         ele_dim = 2;
         nfaces = 4;
         nedges = 4;
         break;
      case hex:
         nnodes = 8;
         nnodesHQ = 20;
         ele_dim = 3;
         nfaces = 6;
         nedges = 12;
         break;
      case tri:
         nnodes = 3;
         nnodesHQ = 6;
         ele_dim = 2;
         nfaces = 3;
         nedges = 3;
         break;
      case tet:
         nnodes = 4;
         nnodesHQ = 10;
         ele_dim = 3;
         nfaces = 4;
         nedges = 6;
         break;
      case prism:
         nnodes = 6;
         nnodesHQ = 15;
         ele_dim = 3;
         nfaces = 5;
         nedges = 9;
         break;
      case pyramid:
         nnodes = 5;
         nnodesHQ = 14;
         ele_dim = 3;
         nfaces = 5;
         nedges = 8;
         break;
	  default:
		break;	
   }
}

//    WW. 06.2005
string Elem::getName() const
{

   switch(ele_Type)
   {
      case line:
         return "line";
      case quadri:
         return "quad";
      case hex:
         return "hex";
      case tri:
         return "tri";
      case tet:
         return "tet";
      case prism:
         return "pris";
      case pyramid:
         return "pyra";
      default:
         return "none";
   }
}


void Elem::setLocalNodeIndex(const int li, const long n_lindex)
{
   nodes[li]->local_index = n_lindex;
}

long Elem::getLocalNodeIndex(const int li) const
{
   return nodes[li]->local_index;
}

//    WW. 06.2005
void Elem::Read(istream& is,  Mesh_Group::Mesh *mesh, int fileType)
{
   //fileType=0: msh
   //fileType=1: rfi
   //fileType=2: gmsh
   //fileType=3: GMS
   //fileType=4: SOL
   int idummy, et;
   string buffer, name;
   idummy=et=-1;
//   is.ignore(numeric_limits<int>::max(), '\n');
   //----------------------------------------------------------------------
   // 1 Reading element type data
   switch(fileType)
   {
         //....................................................................
      case 0: // msh
         is>>index>>PatchIndex;
         is>>buffer;
         if(buffer.find("-1")!=string::npos)
            is>>name;
         else
            name = buffer;
         if(name.find("line")!=string::npos)
            ele_Type = line;
         else if(name.find("quad")!=string::npos)
            ele_Type = quadri;
         else if(name.find("hex")!=string::npos)
            ele_Type = hex;
         else if(name.find("tri")!=string::npos)
            ele_Type = tri;
         else if(name.find("tet")!=string::npos)
            ele_Type = tet;
         else if(name.find("pri")!=string::npos)
            ele_Type = prism;
         else if(name.find("pyra")!=string::npos)
            ele_Type = pyramid;
         break;
         //....................................................................
      case 1: // rfi
         is>>index>>PatchIndex>>name;
         if(name.find("line")!=string::npos)
            ele_Type = line;
         else if(name.find("quad")!=string::npos)
            ele_Type = quadri;
         else if(name.find("hex")!=string::npos)
            ele_Type = hex;
         else if(name.find("tri")!=string::npos)
            ele_Type = tri;
         else if(name.find("tet")!=string::npos)
            ele_Type = tet;
         else if(name.find("pri")!=string::npos)
            ele_Type = prism;
         else if(name.find("pyra")!=string::npos)
            ele_Type = pyramid;
         break;
         //....................................................................
      case 2: // gmsh
         int gmsh_patch_index;
         is>>index>>et>>gmsh_patch_index>>idummy>>nnodes;
         PatchIndex = gmsh_patch_index-1; //OK
         switch(et)
         {
            case 1:
               ele_Type = line;
               break;
            case 2:
               ele_Type = tri;
               break;
            case 3:
               ele_Type = quadri;
               break;
            case 4:
               ele_Type = tet;
               break;
            case 5:
               ele_Type = hex;
               break;
            case 6:
               ele_Type = prism;
               break;
            case 7:
               ele_Type = pyramid;
               break;
			default:
			   break;		
         }
         index--;
         break;
         //....................................................................
      case 3: // GMS
         ele_Type = tri;
         break;
         //....................................................................
      case 4: // SOL
         ele_Type = tri;
         break;
   }
   //----------------------------------------------------------------------
   // 2 Element configuration
   Init();
   nodes.resize(nnodes);
   //----------------------------------------------------------------------
   // 3 Reading element node data
   long nidx = 0;
   switch(fileType)
   {
         //....................................................................
      case 0: // msh
         for(int i=0; i<nnodes; i++)
         {   is >> nidx;		  
             nodes[i] = mesh->node_vector[nidx];
         }
         break;
         //....................................................................
      case 1: // rfi
         for(int i=0; i<nnodes; i++)
         {   is >> nidx;		  
             nodes[i] = mesh->node_vector[nidx];
         }
         break;
         //....................................................................
      case 2: // gmsh
         for(int i=0; i<nnodes; i++)
         {   is >> nidx;		  
             nodes[i] = mesh->node_vector[nidx-1];
         }
         break;
         //....................................................................
      case 3: // GMS
         for(int i=0; i<nnodes; i++)
         {   is >> nidx;		  
             nodes[i] = mesh->node_vector[nidx-1];
         }
         break;
         //....................................................................
      case 4: // SOL
         for(int i=0; i<nnodes; i++)
         {   is >> nidx;		  
             nodes[i] = mesh->node_vector[nidx-1];
         }
         is >> PatchIndex;
         break;
   }
   is>>ws;
   //----------------------------------------------------------------------
   // Initialize topological properties
   neighbors.resize(nfaces);
   for(int i=0; i<nfaces; i++)
      neighbors[i] = NULL;

#ifdef BUILD_MESH_EDGE 
   edges.resize(nedges);
   edges_orientation.resize(nedges);
   for(int i=0; i<nedges; i++)
   {
      edges[i] = NULL;
      edges_orientation[i] = 1;
   }
#endif
}
//  WW. 03.2009
void Elem::WriteGmsh(ostream& os,  const int sdom_idx) const
{
   //int igeo=14;

   int et=1;
   int ntags=3;
   string deli = " ";

   int nn = nnodes;
   if(quadratic)
   {
      nn = nnodesHQ;
      switch(ele_Type)
      {
         case line:
            et = 8;
            break;    //Line
         case quadri:
            et = 10;
            break;    //Quad
         case hex:
            et = 12;
            break;    //Hex
         case tri:
            et = 9;
            break;    //Tri
         case tet:
            et = 11;
            break;    //Tet
         case prism:
            et = 18;
            break;    //Pris
         case pyramid:
            et = 14;
            break;    //Pris
      }
   }
   else
   {
      switch(ele_Type)
      {
         case line:
            et = 1;
            break;    //Line
         case quadri:
            et = 3;
            break;    //Quad
         case hex:
            et = 5;
            break;    //Hex
         case tri:
            et = 2;
            break;    //Tri
         case tet:
            et = 4;
            break;    //Tet
         case prism:
            et = 6;
            break;    //Pris
         case pyramid:
            et = 7;
            break;    //Pris
      }
   }
   os<<index+1<<deli<<et<<deli<<ntags<<deli<<PatchIndex+1<<deli<<PatchIndex+1<<deli<<sdom_idx<<deli;
   for(int i=0; i<nn; i++)
      os<<nodes[i]->index + 1<<deli;
   os<<endl;
}

//  WW. 03.2009
void Elem::WriteGSmsh(ostream& os, bool quad) const
{
   string ename;
   string deli = " ";

   int nn = getNodesNumber(quad);

   switch(ele_Type)
   {
      case line:
         ename = "line";
         break;
      case quadri:
         ename = "quad";
         break;
      case hex:
         ename = "hex";
         break;
      case tri:
         ename = "tri";
         break;
      case tet:
         ename = "tet";
         break;
      case prism:
         ename = "pris";
         break;
      case pyramid:
         ename = "pyra";
         break;
   }
   os<<index<<deli<<PatchIndex<<deli<<ename<<deli;
   for(int i=0; i<nn; i++)
   {
//      nodes_index[i] = nodes[i]->getIndex();
      os<<nodes[i]->getIndex()<<deli;
   }
   os<<endl;
}

//  WW. 03.2012
void Elem::WriteSubDOM(ostream& os, const long node_id_shift, bool quad) const
{
   int nn = getNodesNumber(quad);

   os<<PatchIndex<<" "<<ele_Type+1<<" "<<nn<<" ";
   for(int i=0; i<nn; i++)
   {
//      nodes_index[i] = nodes[i]->getIndex();
      os<<nodes[i]->getIndex()-node_id_shift<<" ";
   }
   os<<endl;
}

//  WW. 02.2012
void Elem::WriteVTK_Type(ostream& os,  bool isquad) const
{
   if(!isquad)
   {
      switch(ele_Type)
      {
         case line:
            os<< "3  "<<endl;
            break;
         case quadri:
            os<< "9  "<<endl;
            break;
         case hex:
            os<< "12 "<<endl;
            break;
         case tri:
            os<< "5  "<<endl;
            break;
         case tet:
            os<< "10 "<<endl;
            break;
         case prism:
            os<< "13 "<<endl;
            break;
         case pyramid:
            os<< "14 "<<endl;
            break;
      }
   }
   else
   {
      switch(ele_Type)
      {
         case line:
            os<< "21  "<<endl;
            break;
         case quadri:
            os<< "23  "<<endl;
            break;
         case hex:
            os<< "25 "<<endl;
            break;
         case tri:
            os<< "22  "<<endl;
            break;
         case tet:
            os<< "24 "<<endl;
            break;
         default:
            cout<<"Warning: quadratic element is not available for prism and pyramid elements"<<endl;				
            break;
      }
   }
}

//    WW. 06.2005
void Elem::WriteIndex(ostream& os) const
{
    string deli = " ";
   os<<index<<deli<<PatchIndex<<deli<<getName()<<deli;
   for(int i=0; i<nnodes; i++)
      os<<nodes[i]->index<<deli;
   os<<endl;
}
void Elem::Write_index(ostream& os) const
{
    string deli = " ";
   if(nodes.Size()>0)
   {
      for(int i=0; i<nnodes; i++)
         os<<nodes[i]->index+1<<deli;
   }
   else
   {
      for(int i=0; i<nnodes; i++)
         os<<nodes[i]->index + 1<<deli;
   }
   os<<endl;
}
//    WW. 06.2005
void Elem::WriteAll(ostream& os) const
{
    string deli = " ";
   os<<index<<deli<<PatchIndex<<deli<<getName()<<deli;
   //if(index==0)
   os<<"Index X Y Z: "<<endl;
   for(int i=0; i<nodes.Size(); i++)
   {
      const Node *anode = nodes[i];  
      os<<anode->index
        <<deli<<anode->X()
        <<deli<<anode->Y()
        <<deli<<anode->Z()<<endl;
   }
}

void Elem::WriteNeighbors(ostream& os) const
{
   os<<"Neighbors of "<<index<<endl;
   for(int i=0; i<nfaces; i++)
      neighbors[i]->WriteAll(os);
   os<<"End neighbors of "<<index<<endl<<endl;;
}

void Elem::MarkingNodes(bool maker)
{
   int SizeV = nnodes;
   if(quadratic) SizeV = nnodesHQ;

   for (int i=0; i< SizeV; i++)
   {
       nodes[i]->Marking(maker);
   }
}


//    WW. 06.2005
void Elem::setNodes(vec<Node*>&  ele_nodes, const bool ReSize)
{
   int SizeV = nnodes;
   if(quadratic) SizeV = nnodesHQ;
   if(ReSize)
   {
      nodes.resize(SizeV);
   }
   for (int i=0; i< SizeV; i++)
   {
      nodes[i] = ele_nodes[i];
   }
}

//   WW  06.2005
//   WW  08.2012
void  Elem::getLocalIndices_EdgeNodes(const int Edge, int *EdgeNodes)
{
   const int start_index  = EdgeLocalIndexArrayElemShift[ele_Type] + Edge*3; 	
   EdgeNodes[0] = EdgeLocalNodeIndex[start_index];
   EdgeNodes[1] = EdgeLocalNodeIndex[start_index + 1];
   EdgeNodes[2] = EdgeLocalNodeIndex[start_index + 2];
}

/**************************************************************************
GetElementFaceNodes
Task: Get local indeces of an element face nodes
Return: number of nodes of a face
Programing:
06/2004 WW
**************************************************************************/
int Elem::getElementFaces1D(int *FaceNode)
{
   FaceNode[0] = 0;
   FaceNode[1] = 1;
   return 2;
}
/**************************************************************************
GetElementFaceNodesTri
Task: Get local indeces of a traingle element face nodes
Augs.:
        const int Face :  Local index of element face
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
Return: number of nodes of a face
Programing:
09/2004 WW
**************************************************************************/
int Elem::getElementFacesTri(const int Face, int *FaceNode)
{
   if(!quadratic)
   {
      FaceNode[0] = Face;
      FaceNode[1] = (Face+1)%3;
      return 2;
   }
   else
   {
      FaceNode[0] = Face;
      FaceNode[1] = (Face+1)%3;
      FaceNode[2] = Face+3;
      return 3;
   }
}

/**************************************************************************
GetElementFaceNodesQuad
Task: Get local indeces of a quadralateral element face nodes
Augs.:
        const int Face :  Local index of element face
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
Return: number of nodes of a face
Programing:
09/2004 WW
**************************************************************************/
int Elem::getElementFacesQuad(const int Face, int *FaceNode)
{
   if(!quadratic)
   {
      FaceNode[0] = Face;
      FaceNode[1] = (Face+1)%4;
      return 2;
   }
   else
   {
      FaceNode[0] = Face;
      FaceNode[1] = (Face+1)%4;
      FaceNode[2] = Face+4;
      return 3;
   }
}

/**************************************************************************
GetElementFaceNodesHex
Task: Get local indeces of a hexahedra element face nodes
Augs.:
        const int Face :  Local index of element face
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
Return: number of nodes of a face
Programing:
09/2004 WW
**************************************************************************/
int Elem::getElementFacesHex(const int Face, int *FaceNode)
{
   int nn=4, k = 0;
   if(quadratic) nn = 8;
   switch(Face)
   {
      case 0:
         for(k=0; k<4; k++)
            FaceNode[k] = k;
         if(quadratic)
         {
            for(k=0; k<4; k++)
               FaceNode[k+4] = k+8;
         }
         break;
      case 1:
         for(k=0; k<4; k++)
            FaceNode[k] = k+4;
         if(quadratic)
         {
            for(k=0; k<4; k++)
               FaceNode[k+4] = k+12;
         }
         break;
      case 2:
         FaceNode[0] = 0;
         FaceNode[1] = 4;
         FaceNode[2] = 5;
         FaceNode[3] = 1;
         if(quadratic)
         {
            FaceNode[4] = 16;
            FaceNode[5] = 12;
            FaceNode[6] = 17;
            FaceNode[7] = 8;
         }
         break;
      case 3:
         FaceNode[0] = 1;
         FaceNode[1] = 5;
         FaceNode[2] = 6;
         FaceNode[3] = 2;
         if(quadratic)
         {
            FaceNode[4] = 17;
            FaceNode[5] = 13;
            FaceNode[6] = 18;
            FaceNode[7] = 9;
         }

         break;
      case 4:
         FaceNode[0] = 2;
         FaceNode[1] = 6;
         FaceNode[2] = 7;
         FaceNode[3] = 3;
         if(quadratic)
         {
            FaceNode[4] = 18;
            FaceNode[5] = 14;
            FaceNode[6] = 19;
            FaceNode[7] = 10;
         }
         break;
      case 5:
         FaceNode[0] = 0;
         FaceNode[1] = 3;
         FaceNode[2] = 7;
         FaceNode[3] = 4;
         if(quadratic)
         {
            FaceNode[4] = 11;
            FaceNode[5] = 19;
            FaceNode[6] = 15;
            FaceNode[7] = 16;
         }
         break;
   }
   return nn;
}


/**************************************************************************
GetElementFaceNodesTet
Task: Get local indeces of a Tedrahedra element face nodes
Augs.:
        const int Face :  Local index of element face
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
Return: number of nodes of a face
Programing:
09/2004 WW
**************************************************************************/
int Elem::getElementFacesTet(const int Face, int *FaceNode)
{
   int nn=3;
   if(quadratic) nn =6;
   switch(Face)
   {
      case 0:
         FaceNode[0] = 1;
         FaceNode[1] = 2;
         FaceNode[2] = 3;
         if(quadratic)
         {
            FaceNode[3] = 5 ;
            FaceNode[4] = 8;
            FaceNode[5] = 7;
         }
         break;
      case 1:
         FaceNode[0] = 3;
         FaceNode[1] = 2;
         FaceNode[2] = 0;
         if(quadratic)
         {
            FaceNode[3] = 8 ;
            FaceNode[4] = 6;
            FaceNode[5] = 9;
         }
         break;
      case 2:
         FaceNode[0] = 1;
         FaceNode[1] = 3;
         FaceNode[2] = 0;
         if(quadratic)
         {
            FaceNode[3] = 7 ;
            FaceNode[4] = 9;
            FaceNode[5] = 4;
         }
         break;
      case 3:
         FaceNode[0] = 0;
         FaceNode[1] = 2;
         FaceNode[2] = 1;
         if(quadratic)
         {
            FaceNode[3] = 6 ;
            FaceNode[4] = 5;
            FaceNode[5] = 4;
         }
         break;

   }
   return nn;
}


/**************************************************************************
GetElementFaceNodesPri
Task: Get local indeces of a prismal element face nodes
Augs.:
        const int Face :  Local index of element face
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
Return: number of nodes of a face
Programing:
09/2004 WW
**************************************************************************/
int Elem::getElementFacesPri(const int Face, int *FaceNode)
{
   int nn=3, k = 0;
   switch(Face)
   {
      case 0:
         nn = 3;
         for(k=0; k<3; k++)
            FaceNode[k] = k;
         if(quadratic)
         {
            for(k=0; k<3; k++)
               FaceNode[k+3] = k+6;
            nn = 6;
         }
         break;
      case 1:
         for(k=0; k<3; k++)
            FaceNode[k] = k+3;
         nn = 3;
         if(quadratic)
         {
            for(k=0; k<3; k++)
               FaceNode[k+3] = k+9;
            nn = 6;
         }
         break;
      case 2:
         FaceNode[0] = 1;
         FaceNode[1] = 2;
         FaceNode[2] = 5;
         FaceNode[3] = 4;
         nn = 4;
         if(quadratic)
         {
            FaceNode[4] = 7 ;
            FaceNode[5] = 14;
            FaceNode[6] = 10;
            FaceNode[7] = 13;
            nn = 8;
         }
         break;
      case 3:
         FaceNode[0] = 5;
         FaceNode[1] = 2;
         FaceNode[2] = 0;
         FaceNode[3] = 3;
         nn = 4;
         if(quadratic)
         {
            FaceNode[4] = 14 ;
            FaceNode[5] =  8;
            FaceNode[6] = 12;
            FaceNode[7] = 10;
            nn = 8;
         }
         break;
      case 4:
         FaceNode[0] = 0;
         FaceNode[1] = 1;
         FaceNode[2] = 4;
         FaceNode[3] = 3;
         nn = 4;
         if(quadratic)
         {
            FaceNode[4] = 6 ;
            FaceNode[5] = 13;
            FaceNode[6] = 9;
            FaceNode[7] = 12;
            nn = 8;
         }
         break;

   }
   return nn;
}

// 08.2012. WW
int Elem::getElementFacesPyramid(const int Face, int *FaceNode)
{
   int nn = 3;

   switch(Face)
   {
      case 0:
         nn= 4;
         FaceNode[0] = 0;
         FaceNode[1] = 1;
         FaceNode[2] = 2;
         FaceNode[3] = 3;
         if(quadratic)
         {
            FaceNode[4] = 5 ;
            FaceNode[5] = 6;
            FaceNode[6] = 7;
            FaceNode[7] = 8;
            nn = 8;
         }
         break;
      case 1:
         nn = 3;
         FaceNode[0] = 0;
         FaceNode[1] = 1;
         FaceNode[2] = 4;
         if(quadratic)
         {
            nn = 6;
            FaceNode[3] = 5 ;
            FaceNode[4] = 10;
            FaceNode[5] = 9;
         }
         break;
      case 2:
         nn = 3;
         FaceNode[0] = 1;
         FaceNode[1] = 2;
         FaceNode[2] = 4;
         if(quadratic)
         {
            nn = 6;
            FaceNode[3] = 6 ;
            FaceNode[4] = 11;
            FaceNode[5] = 10;
         }
         break;
      case 3:
         nn = 3;
         FaceNode[0] = 2;
         FaceNode[1] = 3;
         FaceNode[2] = 4;
         if(quadratic)
         {
            nn = 6;
            FaceNode[3] = 7 ;
            FaceNode[4] = 12;
            FaceNode[5] = 11;
         }
         break;
      case 4:
         nn = 3;
         FaceNode[0] = 3;
         FaceNode[1] = 0;
         FaceNode[2] = 4;
         if(quadratic)
         {
            nn = 6;
            FaceNode[3] = 8 ;
            FaceNode[4] = 9;
            FaceNode[5] = 12;
         }
         break;
   }
   return nn;
}



/**************************************************************************
GetElementFaces
Task: set element faces (Geometry)
Augs.:
        const int Face :  Local index of element face
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes

Programing:
09/2004 WW
**************************************************************************/
int Elem::getElementFaceNodes(const int Face, int *FacesNode)
{
   switch(ele_Type)
   {
      case line:  // 1-D bar element
         return getElementFaces1D(FacesNode);
      case quadri: // 2-D quadrilateral element
         return getElementFacesQuad(Face, FacesNode);
      case hex: // 3-D hexahedral element
         return getElementFacesHex(Face, FacesNode);
      case tri:  // 2-D triagular element
         return getElementFacesTri(Face, FacesNode);
      case tet:  // 3-D tetrahedral element
         return getElementFacesTet(Face, FacesNode);
      case prism:
         return getElementFacesPri(Face, FacesNode);
	  case pyramid:	
		 return getElementFacesPyramid(Face, FacesNode);
	  default:
        break;
   }
   return 0;
}



}//end namespace

