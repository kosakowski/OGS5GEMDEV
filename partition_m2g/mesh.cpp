#include <iomanip>
#include <limits>
#include <cstdlib> // For GCC 4.3.2
#include <string.h>  // For GCC 4.3.2
//#include <sstream>
#include "mesh.h"
//#include "itoa.h"  // Only for GCC

using namespace std;

long NodesNumber_Linear =0;
long NodesNumber_Quadratic =0;



//------------------------------------------------------ 
//   Topology definition of geometrical element.
//    WW. 10.01.2005
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
   quadratic = false;
   deli = "  ";
 }

//1.  Node declaration
 //    WW. 06.2005
Node:: Node(const int Index, const double x, 
            const double y, const double z):Grain(Index)
{
  
   Coordinate[0] = x;
   Coordinate[1] =y;
   Coordinate[2] =z;
   boundayC = 'I';
  
}
//    WW. 06.2005
void Node::operator = (const Node& n)
{
   boundayC = n.boundayC; 
   index = n.index;
   mark = n.mark;
   Coordinate[0] = n.Coordinate[0]; 
   Coordinate[1] = n.Coordinate[1]; 
   Coordinate[2] = n.Coordinate[2]; 
}
//    WW. 06.2005
bool Node::operator == (const Node& n)
{
   if(index == n.index)
      return true;
   else 
	  return false;
}
//    WW. 06.2005
// Output
void Node::Write(ostream& osm) const
{
    osm.setf(ios::scientific, ios::floatfield);
    setw(14);
    osm.precision(14);
    osm<<index<<deli
       <<Coordinate[0]<<deli
       <<Coordinate[1]<<deli
       <<Coordinate[2]<<endl;
}
//    WW. 06.2005
// Set
void Node::SetCoordinates(const double* argCoord)
{
	  Coordinate[0] = argCoord[0];
  	Coordinate[1] = argCoord[1];
  	Coordinate[2] = argCoord[2];
}
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
   boundayC = ed.boundayC; 
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
   Volume = 0.0;
   
   Owner = NULL;
}
//    WW. 06.2005
Elem::~Elem()
{
   nodes_index.resize(0);
   locnodes_index.resize(0);
   nodes.resize(0);
   edges.resize(0);
   neighbors.resize(0);
}

//    WW. 06.2005
Elem::  Elem( const int Index,  Elem* onwer, const int Face):
  Grain(Index), Owner(onwer)
{
   int i, j, k, n, ne; 
   static int faceIndex_loc[10];
   static int edgeIndex_loc[10];
 //  Owner = onwer;
   n = Owner->GetElementFaceNodes(Face, faceIndex_loc);
   switch(Owner->ele_Type)
   {
       case 1:  // 1-D bar element
           break;          
       case 2: // 2-D quadrilateral element
           nnodes = 2;
           nnodesHQ = 3;    
           ele_dim = 1;
           ele_Type = 1;
           nfaces = 2;
           nedges = 0;
           break;           
       case 3: // 3-D hexahedral element 
           nnodes = 4;
           nnodesHQ = 8;      
           ele_dim = 2;
           ele_Type = 2;
           nfaces = 4;
           nedges = 4;
           break;           
       case 4:  // 2-D triagular element 
           nnodes = 2;
           nnodesHQ = 3;    
           ele_dim = 1;
           ele_Type = 1;
           nfaces = 2;
           nedges = 0;
           break;           
       case 5:  // 3-D tetrahedral element 
           nnodes = 3;
           nnodesHQ = 6;
           ele_dim = 2;
           ele_Type = 4;
           nfaces = 3;
           nedges = 3;
           break;           
       case 6: 
          if(Face<2) 
          {
              nnodes = 3;
              nnodesHQ = 6;
              ele_dim = 2;
              ele_Type = 4;
              nfaces = 3;
              nedges = 3;
           }
           else  
           {
              nnodes = 4;
              nnodesHQ = 8;      
              ele_dim = 2;
              ele_Type = 2;
              nfaces = 4;
              nedges = 4;              
           }
           break; // 3-D prismatic element 
    }

    PatchIndex =  Owner->PatchIndex;
    quadratic = Owner->quadratic;
    nodes_index.resize(n);
    nodes.resize(n);

    boundayC='B';
    for(i=0; i<n; i++)
    {
       nodes_index[i] =
                  Owner->nodes_index[faceIndex_loc[i]];
       nodes[i] = Owner->nodes[faceIndex_loc[i]];
							nodes[i]->boundayC = 'B';
    }
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
           Owner->GetLocalIndices_EdgeNodes(j, edgeIndex_loc);    
           if( (faceIndex_loc[i]==edgeIndex_loc[0]&&
                faceIndex_loc[k]==edgeIndex_loc[1])||															
               (faceIndex_loc[i]==edgeIndex_loc[1]&&
                faceIndex_loc[k]==edgeIndex_loc[0]) )
           {
               edges[i] = Owner->edges[j];
               if(faceIndex_loc[i]==edgeIndex_loc[1]&&
                      faceIndex_loc[k]==edgeIndex_loc[0] )
               edges_orientation[i] = -1; 
               edges[i]->boundayC = 'B';
               break;
            } 
         }
    } 
}

//    WW. 06.2005
string Elem::getName() const
{
  
   switch(ele_Type)
   {
     case 1:
       return "line";
       break;
     case 2:
       return "quad";
       break;
     case 3:
       return "hex";
       break;
     case 4:
       return "tri";
       break;
     case 5:
       return "tet";
       break;
     case 6:
       return "pris";
       break;
     default:
       return "none";
       break;
   }
}
//    WW. 06.2005
void Elem::Read(istream& is, int fileType)
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
  switch(fileType){
    //....................................................................
    case 0: // msh
      is>>index>>PatchIndex;
      is>>buffer;
	  if(buffer.find("-1")!=string::npos)
         is>>name;
	  else
	    name = buffer;
      if(name.find("line")!=string::npos)
         ele_Type = 1;
      else if(name.find("quad")!=string::npos)
         ele_Type = 2;
      else if(name.find("hex")!=string::npos)
         ele_Type = 3;
      else if(name.find("tri")!=string::npos)
         ele_Type = 4;
      else if(name.find("tet")!=string::npos)
         ele_Type = 5;
      else if(name.find("pri")!=string::npos)
         ele_Type = 6;
      break;
    //....................................................................
    case 1: // rfi
	  is>>index>>PatchIndex>>name;
      if(name.find("line")!=string::npos)
         ele_Type = 1;
      else if(name.find("quad")!=string::npos)
         ele_Type = 2;
      else if(name.find("hex")!=string::npos)
         ele_Type = 3;
      else if(name.find("tri")!=string::npos)
         ele_Type = 4;
      else if(name.find("tet")!=string::npos)
         ele_Type = 5;
      else if(name.find("pri")!=string::npos)
         ele_Type = 6;
      break;
    //....................................................................
    case 2: // gmsh
      int gmsh_patch_index; 
      is>>index>>et>>gmsh_patch_index>>idummy>>nnodes;
      PatchIndex = gmsh_patch_index-1; //OK
      switch(et)
      {
         case 1: ele_Type = 1; break;
         case 2: ele_Type = 4; break;
         case 3: ele_Type = 2; break;
         case 4: ele_Type = 5; break;
         case 5: ele_Type = 3; break;
         case 6: ele_Type = 6; break;
      }
	  index--;
      break;
    //....................................................................
    case 3: // GMS
      ele_Type = 4;
      break;
    //....................................................................
    case 4: // gmsh
      ele_Type = 4;
      break;
  }
  //----------------------------------------------------------------------
  // 2 Element configuration
  switch(ele_Type)
   {
      case 1:
         nnodes = 2;
         nnodesHQ = 3;    
         ele_dim = 1;
         ele_Type = 1;
         nfaces = 2;
         nedges = 0;
         break;
      case 2:
         nnodes = 4;
         nnodesHQ = 9;      
         ele_dim = 2;
         ele_Type = 2;
         nfaces = 4;
         nedges = 4;
         break;
      case 3:
         nnodes = 8;
         nnodesHQ = 20;
         ele_dim = 3;
         nfaces = 6;
         nedges = 12;
         ele_Type = 3;
         break;
      case 4:
         nnodes = 3;
         nnodesHQ = 6;
         ele_dim = 2;
         ele_Type = 4;
         nfaces = 3;
         nedges = 3;
         break;
      case 5:
         nnodes = 4;
         nnodesHQ = 10;
         ele_dim = 3;
         ele_Type = 5;
         nfaces = 4;
         nedges = 6;
         break;
      case 6:
         nnodes = 6;
         nnodesHQ = 15;
         ele_dim = 3;
         ele_Type = 6;
         nfaces = 5;
         nedges = 9;
		 break;
   }
   nodes_index.resize(nnodes);
  //----------------------------------------------------------------------
  // 3 Reading element node data
  switch(fileType){
    //....................................................................
    case 0: // msh
      for(int i=0; i<nnodes; i++)
        is>>nodes_index[i];
      break;
    //....................................................................
    case 1: // rfi
      for(int i=0; i<nnodes; i++)
        is>>nodes_index[i];
      break;
    //....................................................................
    case 2: // gmsh
      for(int i=0; i<nnodes; i++){
        is>>nodes_index[i];
        nodes_index[i] -= 1;
      }
      break;
    //....................................................................
    case 3: // GMS
      for(int i=0; i<nnodes; i++){
        is>>nodes_index[i];
        nodes_index[i] -= 1;
      }
      break;
    //....................................................................
    case 4: // SOL
      for(int i=0; i<nnodes; i++){
        is>>nodes_index[i];
        nodes_index[i] -= 1;
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
  edges.resize(nedges);    
  edges_orientation.resize(nedges);
  for(int i=0; i<nedges; i++)
  {
    edges[i] = NULL;
	edges_orientation[i] = 1;
  }
}
//    WW. 06.2005
void Elem::WriteIndex(ostream& os) const
{
    os<<index<<deli<<PatchIndex<<deli<<getName()<<deli;
    for(int i=0; i<nnodes; i++)
      os<<nodes_index[i]<<deli;
	   os<<endl;
}
void Elem::Write_index(ostream& os) const
{
    for(int i=0; i<nnodes; i++)
      os<<nodes_index[i]+1<<deli;
    os<<endl;
}
//    WW. 06.2005
void Elem::WriteAll(ostream& os) const
{
    os<<index<<deli<<PatchIndex<<deli<<getName()<<deli;
    //if(index==0) 
    os<<"Index X Y Z: "<<endl;
				for(int i=0; i<nodes.Size(); i++)
    {
       os<<nodes_index[i]
       <<deli<<nodes[i]->X()
       <<deli<<nodes[i]->Y()
       <<deli<<nodes[i]->Z()<<endl;
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
       	nodes_index.resize(SizeV);
    }
    for (int i=0; i< SizeV;i++)
    {
       nodes[i] = ele_nodes[i];
      nodes_index[i] = nodes[i]->GetIndex();       
    }
}



void Elem::WriteNeighbors(ostream& os) const
{
    os<<"Neighbors of "<<index<<endl;
    for(int i=0; i<nfaces; i++)
      neighbors[i]->WriteAll(os);
    os<<"End neighbors of "<<index<<endl<<endl;;
}
//    WW. 06.2005
void  Elem::GetLocalIndices_EdgeNodes(const int Edge, int *EdgeNodes)
{
	switch(ele_Type)
	{
       case 1: 
           break; // 1-D bar element 
       case 2: // 2-D quadrilateral element 
          EdgeNodes[0] = Edge;
          EdgeNodes[1] = (Edge+1)%4;
          break;             
       case 3: // 3-D hexahedral element
          if(Edge<8)
          {
             EdgeNodes[0] = Edge;
             EdgeNodes[1] = (Edge+1)%4+4*(int)(Edge/4);
          }
          else 
          {
             EdgeNodes[0] = Edge%4;
             EdgeNodes[1] = Edge%4+4;
          }
          break;  
       case 4:  // 2-D triagular element 
          EdgeNodes[0] = Edge;
          EdgeNodes[1] = (Edge+1)%3;
          break;
       case 5:  // 3-D tetrahedra
          if(Edge<3)
          {
             EdgeNodes[0] = Edge;
             EdgeNodes[1] = (Edge+1)%3;
          }
          else
          {
             EdgeNodes[0] = 3;
             EdgeNodes[1] = (Edge+1)%3;
          }

          break;  
       case 6: // 3-D prismatic element
          if(Edge<6)
          {
             EdgeNodes[0] = Edge;
             EdgeNodes[1] = (Edge+1)%3+3*(int)(Edge/3);
          }
          else 
          {
             EdgeNodes[0] = Edge%3;
             EdgeNodes[1] = Edge%3+3;
          }
          break;  
	}
}


/**************************************************************************
GetElementFaceNodes
Task: Get local indeces of an element face nodes
Return: number of nodes of a face
Programing:
06/2004 WW  
**************************************************************************/
int Elem::GetElementFaces1D(int *FaceNode)
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
int Elem::GetElementFacesTri(const int Face, int *FaceNode)
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
       FaceNode[1] = Face+3;
       FaceNode[2] = (Face+1)%3;
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
int Elem::GetElementFacesQuad(const int Face, int *FaceNode)
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
       FaceNode[1] = Face+4;
       FaceNode[2] = (Face+1)%4;
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
int Elem::GetElementFacesHex(const int Face, int *FaceNode)
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
int Elem::GetElementFacesTet(const int Face, int *FaceNode)
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
int Elem::GetElementFacesPri(const int Face, int *FaceNode)
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
int Elem::GetElementFaceNodes(const int Face, int *FacesNode)
{
   switch(ele_Type)
   {
       case 1:  // 1-D bar element
           return GetElementFaces1D(FacesNode);
           break;          
       case 2: // 2-D quadrilateral element
           return GetElementFacesQuad(Face, FacesNode);
           break;           
       case 3: // 3-D hexahedral element 
           return GetElementFacesHex(Face, FacesNode);
           break;           
       case 4:  // 2-D triagular element 
           return GetElementFacesTri(Face, FacesNode);
           break;           
       case 5:  // 3-D tetrahedral element 
           return GetElementFacesTet(Face, FacesNode);
           break;           
       case 6: 
           return GetElementFacesPri(Face, FacesNode);
           break; // 3-D prismatic element 
    }
    return 0;
}



}//end namespace


// Functions
using Mesh_Group::Grain;
using Mesh_Group::Node;
using Mesh_Group::Edge;
using Mesh_Group::Elem;
using Math_Group::vec;

vector<Mesh_Group::Node*> NodesVector;
// All surface nodes
vector<Mesh_Group::Elem*> SurfaceFaces;
// All edges
vector<Mesh_Group::Edge*> EdgeVector;
// All elements 
vector<Mesh_Group::Elem*> ElementsVector;



// Read grid for test purpose
void ReadGrid(istream& is)
{
   long i, ne, nn, counter;
   int ibuff;
   double x,y,z;
   string buffer;
//   is.seekg(position);
   // Read description
   is>>buffer>>ws;
   // Read numbers of nodes and elements
   is>>ibuff>>nn>>ne>>ws;
   if(nn==0||ne==0)
   {
       cout<<"Error: number of elements or nodes is zero"<<endl;
       abort();
   }
   
   // Read Nodes
   counter = 0;
   for(i=0; i<nn; i++)
   {
      is>>ibuff>>x>>y>>z>>ws;
      Node* newNode = new Node(ibuff,x,y,z);
      NodesVector.push_back(newNode);            
      counter++;
   }
   if(counter!=nn)
   {
       cout<<"Error: number nodes do not match"<<endl;
       abort();
   }
   NodesNumber_Linear = nn;
   NodesNumber_Quadratic = nn;
      
   // Read Elements
   counter = 0;
   for(i=0; i<ne; i++)
   { 
      Elem* newElem = new Elem(i);
      newElem->Read(is);
      ElementsVector.push_back(newElem);           
      counter++;
   }     
   if(counter!=ne)
   {
       cout<<"Error: number elements do not match"<<endl;
       abort();
   }

//   position = is.tellg();
}

// Construct grid
// 
void ConstructGrid( const bool quadratic)
{
   int counter;	 
   int i, j, k, ii, jj, m0, m, n0, n;
			int nnodes0, nedges0, nedges;
   long e, ei, ee,  e_size,  e_size_l;
   bool done;
   double x0,y0,z0;

   int edgeIndex_loc0[2];
   int edgeIndex_loc[2];
   int faceIndex_loc0[10];
   int faceIndex_loc[10];
   vec<Node*> e_nodes0(20);
   vec<long> node_index_glb(20);
   vec<long> node_index_glb0(20);
   vec<int> Edge_Orientation(15);
   vec<Edge*> Edges(15);
   vec<Edge*> Edges0(15);
   vec<Elem*> Neighbors(15);
   vec<Elem*> Neighbors0(15);

   vec<Node*> e_edgeNodes0(3);
   vec<Node*> e_edgeNodes(3);
   Elem* thisElem0=NULL;
   Elem* thisElem=NULL;

   //Elem->nodes not initialized

   e_size = (long)ElementsVector.size();  
   Edge_Orientation = 1;

   // Set neighbors
   for(e=0; e<e_size; e++)
   {
       thisElem0 = ElementsVector[e];   
       thisElem0->getNodeIndeces(node_index_glb0);
       for(i=0; i<thisElem0->getNodesNumber(); i++)
       {
          done = false;
          for(j=0; j<(int)NodesVector[node_index_glb0[i]]
					         ->ElementsBelonged.size(); j++)
          {
            if(e==NodesVector[node_index_glb0[i]]
                            ->ElementsBelonged[j])
              done = true;
              break;
          }
          if(!done)  
          NodesVector[node_index_glb0[i]]->ElementsBelonged.push_back(e);
      }
   }

   // Compute neighbors and edges
   for(e=0; e<e_size; e++)
   {
       thisElem0 = ElementsVector[e];   
       nnodes0 = thisElem0->getNodesNumber(); // Number of nodes for linear element
       thisElem0->getNodeIndeces(node_index_glb0);
       thisElem0->getNeighbors(Neighbors0);
       for(i=0; i<nnodes0; i++) // Nodes
         e_nodes0[i] = NodesVector[node_index_glb0[i]];  
       m0 = thisElem0->getFacesNumber();
	      // neighbors
       for(i=0; i<m0; i++) // Faces
       {
          if(Neighbors0[i])
               continue;
          n0 = thisElem0->GetElementFaceNodes(i, faceIndex_loc0);
          done = false;  
          for(k=0;k<n0;k++)
          {    
             e_size_l = (long)e_nodes0[faceIndex_loc0[k]]->ElementsBelonged.size();         
             for(ei=0; ei<e_size_l; ei++)
             {
                ee = e_nodes0[faceIndex_loc0[k]]->ElementsBelonged[ei];   
                if(ee==e) continue;
                thisElem = ElementsVector[ee];   
                thisElem->getNodeIndeces(node_index_glb);
                thisElem->getNeighbors(Neighbors);
                m = thisElem->getFacesNumber();

                for(ii=0; ii<m; ii++) // Faces
                {
                   n = thisElem->GetElementFaceNodes(ii, faceIndex_loc);
                   if(n0!=n) continue;
                   counter = 0;
                   for(j=0; j<n0; j++) 
                   {             
                      for(jj=0; jj<n; jj++) 
                      {
                         if(node_index_glb0[faceIndex_loc0[j]]
                             ==node_index_glb[faceIndex_loc[jj]])    
                         {
                             counter++;
                             break;
                         }
                      }
                   }
                   if(counter==n)
                   {
                       Neighbors0[i] = thisElem;
                       Neighbors[ii] = thisElem0;  
                       thisElem->setNeighbor(ii, thisElem0);
                       done = true;
                       break;                 
                   }
                }
                if(done) break;
             }
             if(done) break;
          }
       }
       thisElem0->setNeighbors(Neighbors0);						
       // 
       // Edges
       nedges0 = thisElem0->getEdgesNumber();
       thisElem0->getEdges(Edges0);
       for(i=0; i<nedges0; i++)
       { 
          thisElem0->GetLocalIndices_EdgeNodes(i, edgeIndex_loc0);    
          // Check neighbors 
          done = false; 
          for(k=0;k<2;k++)
          {    
             e_size_l = (long)e_nodes0[edgeIndex_loc0[k]]->ElementsBelonged.size();         
             for(ei=0; ei<e_size_l; ei++)
             {
                ee = e_nodes0[edgeIndex_loc0[k]]->ElementsBelonged[ei];   
                if(ee==e) continue;
                thisElem = ElementsVector[ee];                   
                thisElem->getNodeIndeces(node_index_glb);
                nedges = thisElem->getEdgesNumber();
                thisElem->getEdges(Edges);
                // Edges of neighbors
                for(ii=0; ii<nedges; ii++)
                { 
                    thisElem->GetLocalIndices_EdgeNodes(ii, edgeIndex_loc);
                    if((  node_index_glb0[edgeIndex_loc0[0]]==node_index_glb[edgeIndex_loc[0]]
                        &&node_index_glb0[edgeIndex_loc0[1]]==node_index_glb[edgeIndex_loc[1]])
				                 ||(  node_index_glb0[edgeIndex_loc0[0]]==node_index_glb[edgeIndex_loc[1]]
                        &&node_index_glb0[edgeIndex_loc0[1]]==node_index_glb[edgeIndex_loc[0]]) )
                     {
                         if(Edges[ii])
                         {
                            Edges0[i] = Edges[ii]; 
                            Edges[ii]->getNodes(e_edgeNodes); 
                            if(  node_index_glb0[edgeIndex_loc0[0]]==e_edgeNodes[1]->GetIndex()
                             && node_index_glb0[edgeIndex_loc0[1]]==e_edgeNodes[0]->GetIndex())
			                             Edge_Orientation[i] = -1; 
                            if(quadratic)  // Get middle node
                            {
                               node_index_glb0[nnodes0] = e_edgeNodes[2]->GetIndex();
                               e_nodes0[nnodes0] = e_edgeNodes[2];
                               nnodes0++;
                            }
                            done = true;
                            break;
                        }
		            }                   
                } //  for(ii=0; ii<nedges; ii++)
                if(done) break;
             } // for(ei=0; ei<e_size_l; ei++)
             if(done) break;
          }//for(k=0;k<2;k++)
          if(!done) // new edges and new node
          {
              Edges0[i] = new Edge((long)EdgeVector.size()); 
              Edges0[i]->SetOrder(quadratic); 
              e_edgeNodes0[0] = e_nodes0[edgeIndex_loc0[0]];
              e_edgeNodes0[1] = e_nodes0[edgeIndex_loc0[1]];
              if(quadratic)  // new node: middle point of edges
              {
                  e_edgeNodes0[2] = new Node((long)NodesVector.size());
                  e_edgeNodes0[2]->setX(0.5*(e_edgeNodes0[0]->X()+e_edgeNodes0[1]->X()));                
                  e_edgeNodes0[2]->setY(0.5*(e_edgeNodes0[0]->Y()+e_edgeNodes0[1]->Y()));                
                  e_edgeNodes0[2]->setZ(0.5*(e_edgeNodes0[0]->Z()+e_edgeNodes0[1]->Z()));    
                  NodesVector.push_back(e_edgeNodes0[2]);
                  node_index_glb0[nnodes0] = e_edgeNodes0[2]->GetIndex();
                  e_nodes0[nnodes0] = e_edgeNodes0[2];
                  nnodes0++;
              }
             Edges0[i]->setNodes(e_edgeNodes0); 
             EdgeVector.push_back(Edges0[i]);		               
          } // new edges
   	  } //  for(i=0; i<nedges0; i++)
      //
      if(quadratic&&thisElem0->getElementType()==2) // Quadrilateral
      {
         x0=y0=z0=0.0;
         Node* newNode = new Node((long)NodesVector.size());
         e_nodes0[nnodes0] = newNode;
         nnodes0 = thisElem0->getNodesNumber();
         for(i=0; i<nnodes0; i++) // Nodes
         {
            x0 += e_nodes0[i]->X();	
            y0 += e_nodes0[i]->Y();	
            z0 += e_nodes0[i]->Z();	
         }         
         x0 /= (double)nnodes0;
         y0 /= (double)nnodes0;
         z0 /= (double)nnodes0;
         newNode->setX(x0);
         newNode->setY(y0);
         newNode->setZ(z0);
         NodesVector.push_back(newNode);         
      }     
      // Set edges and nodes
      thisElem0->SetOrder(quadratic);
      thisElem0->setEdges_Orientation(Edge_Orientation); 
      thisElem0->setEdges(Edges0); 
      // Resize is true
      thisElem0->setNodes(e_nodes0, true);						
   }// Over elements

   // Set faces on surfaces
   for(e=0; e<e_size; e++)
   {
       thisElem0 = ElementsVector[e];   
       thisElem0->getNodeIndeces(node_index_glb0);
       thisElem0->getNeighbors(Neighbors0);
       m0 = thisElem0->getFacesNumber();

       // Check face on surface
       for(i=0; i<m0; i++) // Faces
       {		  
          if(Neighbors0[i])
             continue;
          Elem* newFace = new Elem((long)SurfaceFaces.size(), thisElem0, i);
          SurfaceFaces.push_back(newFace);
          Neighbors0[i] = newFace;        
       }
       thisElem0->setNeighbors(Neighbors0);						
   }
   NodesNumber_Quadratic= (long)NodesVector.size(); 
}

// Read grid for test purpose
void ReadGridGeoSys(istream& is)
{
  string sub_line;
  string line_string;
  bool new_keyword = false;
  string hash("#");
  string sub_string,sub_string1;
  long i, ibuff;
  long no_elements;
  long no_nodes;
  double x,y,z;
  Node* newNode = NULL;
  Elem* newElem = NULL;
  //========================================================================
  // Keyword loop
  while (!new_keyword) {
    //if(!GetLineFromFile(line,fem_file)) 
    //  break;
    //line_string = line;
    getline(is, line_string);
    if(is.fail()) 
      break;
    /*	 
    if(line_string.find(hash)!=string::npos)
	{
      new_keyword = true;
      break;
    }
    */
    //....................................................................
    //....................................................................
    if(line_string.find("$NODES")!=string::npos) { // subkeyword found
      is  >> no_nodes>>ws;
      for(i=0;i<no_nodes;i++){
         is>>ibuff>>x>>y>>z>>ws;
         newNode = new Node(ibuff,x,y,z);
         NodesVector.push_back(newNode);            
      }
      continue;
    }
    //....................................................................
    if(line_string.find("$ELEMENTS")!=string::npos) { // subkeyword found
      is >> no_elements>>ws;
      for(i=0;i<no_elements;i++){
         newElem = new Elem(i);
         newElem->Read(is, 0);
		 ElementsVector.push_back(newElem);
      }
      continue;
    }
  }
  //========================================================================
}


void ConstructDomain(char *fname, const int num_parts)
{
   char str[1028];
   char stro[1028];
   char str_buf[3];
   int dom;
   int max_dom;
   int k,kk;
   long i,j;
   //
   sprintf(str_buf, "%d",num_parts);

   /////////////////vector<int> node_dom;
   strcpy(str,fname);
   strcpy(stro,fname);
   strcat(str,".mesh.epart.");
   strcat(str,str_buf);
   strcat(stro,".");
   strcat(stro,str_buf);
   strcat(stro,"ddc");
   //namef = ".mesh.epart."; //+str_buf;
   ifstream part_in;
   fstream part_out;
   part_out.open(stro, ios::out );
   // Output for gmsh
   strcpy(stro,fname);
   strcat(stro,"_gmsh.msh");
   fstream gmsh_out;
   gmsh_out.open(stro, ios::out );
   gmsh_out<<"$NOD"<<endl;
   gmsh_out<<NodesVector.size()<<endl;
   for(i=0; i<(long)NodesVector.size(); i++)
   {
     gmsh_out<<i+1<<" ";
     gmsh_out<<NodesVector[i]->X()<<" ";
     gmsh_out<<NodesVector[i]->Y()<<" ";
     gmsh_out<<NodesVector[i]->Z()<<endl;
   }
   gmsh_out<<"$ENDNOD"<<endl; 
   gmsh_out<<"$ELM"<<endl; 
   gmsh_out<<(long)ElementsVector.size()<<endl; 
   //
   part_in.open(str);
   if(!part_in.is_open())
   {
       cerr<<("Error: cannot open .epart file . It may not exist !");
       abort();
   } 

   max_dom=0;
   int et;
   for(i=0; i<(long)ElementsVector.size(); i++)
   {
      part_in>>dom>>ws;
	  ElementsVector[i]->setDomainIndex(dom);
//      ElementsVector[i]->AllocateLocalIndexVector(); 
      if(dom>max_dom) max_dom = dom;
      // GMSH output
      switch(ElementsVector[i]->getElementType())
      {
         case 1: et = 1; break;
         case 4: et = 2; break;
         case 2: et = 3; break;
         case 5: et = 4; break;
         case 3: et = 5; break;
         case 6: et = 6; break;
      }
      gmsh_out<<i+1<<" "<<et<<" "<<dom+1<<" "<<dom+1<<" ";
      gmsh_out<<ElementsVector[i]->getNodesNumber()<<"  ";
      for(k=0; k<ElementsVector[i]->getNodesNumber(); k++)
        gmsh_out<<ElementsVector[i]->getNodeIndex(k)+1<<" ";
      gmsh_out<<endl;
   }
   max_dom++;
   gmsh_out<<"$ENDELM"<<endl; 
   gmsh_out.close();
   part_in.close();
   //
   /*
   strcpy(str,fname);
   strcat(str,".npart");
   ifstream npart_in;

   npart_in.open(str);
   if(!npart_in.is_open())
   {
       cerr<<("Error: cannot open .npart file . It may not exist !");
       abort();
   } 
   for(i=0; i<(long)NodesVector.size(); i++)
   {
      npart_in>>dom>>ws;
	  node_dom.push_back(dom);
   }
   npart_in.close();
   */
   //Output ddc file
   // long *nod_dom = new long[max_dom];
   long *ele_dom = new long[max_dom];
   for(k=0; k<max_dom; k++)
   {
      ele_dom[k]=0;
      //nod_dom[k]=0;
      for(j=0; j<(long)ElementsVector.size(); j++)
      {
	     if(ElementsVector[j]->GetDomainIndex()==k)
             ele_dom[k] += 1;
      } 
      /*
      for(j=0; j<(long)NodesVector.size(); j++)
      {
	     if(node_dom[j]==k)
             nod_dom[k] += 1;
      } 
      */
   }

   Elem *ele=0;
   bool done = false;
   long n_index=0;
   vector<int> nodes_dom;
   vector<int> eles_dom;

   //TEST
   Node *nod = 0;
   //


   //
   for(k=0; k<max_dom; k++)
   {      	   
	  part_out<<"#DOMAIN "<<k<<endl;   
	  part_out<<"$ELEMENTS "<<ele_dom[k]<<endl; 
	  nodes_dom.clear(); 
	  eles_dom.clear();
      for(j=0; j<(long)ElementsVector.size(); j++)
	  {
         ele = ElementsVector[j];
         for(kk=0; kk<ele->getNodesNumber(); kk++)
		 {
            ele->SetLocalNodeIndex(kk, -1);
            ele->AllocateLocalIndexVector();
            ele->setDomNodeIndex(kk, -1);
		 }
      }
      for(j=0; j<(long)ElementsVector.size(); j++)
      {
         ele = ElementsVector[j]; 
		 //ele->AllocateLocalIndexVector();
	     if(ele->GetDomainIndex()==k)   
		 {
            for(kk=0; kk<ele->getNodesNumber(); kk++)
			{  
               done = false;
               n_index = ele->GetLocalNodeIndex(kk);
               if(n_index>-1)
               {
                  ele->setDomNodeIndex(kk, n_index);                    
                  done = true;
               }
 			   if(!done)
			   { 
				   ele->setDomNodeIndex(kk, (long)nodes_dom.size()); //For test output
                   ele->SetLocalNodeIndex(kk, (long)nodes_dom.size());
				   nodes_dom.push_back(ele->getNodeIndex(kk));
			   } 
			}
            part_out<<ele->GetIndex()<<endl;
			eles_dom.push_back(ele->GetIndex()); //TEST OUT
		 }
      }        
	  part_out<<"$NODES_INNER "<<(long)nodes_dom.size()<<endl;
      for(j=0; j<(long)nodes_dom.size(); j++)
          part_out<<nodes_dom[j]<<endl;
      /*
	  part_out<<"$NODES_INNER "<<nod_dom[k]<<endl;
      for(j=0; j<(long)NodesVector.size(); j++)
      {
	     if(node_dom[j]==k)
            part_out<<NodesVector[j]->GetIndex()<<endl;
      } 
	  */
      //TEST OUT
      //_itoa(k,stro, 10);
      sprintf(stro, "%d",k);

      strcpy(str,fname);
      string aa = str;
      string bb= stro;
      string name_f = aa+bb+".msh";
      fstream test_out;
      test_out.open(name_f.c_str(), ios::out );

      //GMSH test_out<<"$NOD"<<endl;
 	  //GMSH test_out<<(long)nodes_dom.size()<<endl;
      test_out<<"#0#0#0#1#0.0#0#################################################################"<<endl;
      test_out<<"0 "<<(long)nodes_dom.size()<<" "<<(long)eles_dom.size()<<endl;
      for(j=0; j<(long)nodes_dom.size(); j++)
	  {
          nod = NodesVector[nodes_dom[j]];
        //GMSH  test_out<<j+1<<"  "
          test_out<<j<<"  "
                  << nod->X()<<"  "<< nod->Y()<<"  "<< nod->Z() <<endl;
	  }
      //GMSH test_out<<"$ENDNOD"<<endl;
      //GMSH test_out<<"$ELE"<<endl;
	  //GMSG test_out<<(long)eles_dom.size()<<endl;
      for(j=0; j<(long)eles_dom.size(); j++)
	  {
          ele = ElementsVector[eles_dom[j]]; 
          /* //GMSH
          test_out<<j+1<<"  ";
          int e_t=1;       
          switch(ele->getElementType())
          {
            case 1: e_t = 1; break;
            case 4: e_t = 2; break;
            case 2: e_t = 3; break;
            case 5: e_t = 4; break;
            case 3: e_t = 5; break;
            case 6: e_t = 6; break;
          }
          test_out<<e_t<<" 1 23 "<<ele->getNodesNumber()<<" ";

          for(kk=0; kk<ele->getNodesNumber(); kk++)
             test_out<< ele->GetDomNodeIndex(kk)+1<<" ";
          */
		  test_out<<j<<"  "<<ele->getPatchIndex()<<" "<<ele->getName()<<" ";
          for(kk=0; kk<ele->getNodesNumber(); kk++)
             test_out<< ele->GetDomNodeIndex(kk)<<" ";
          test_out<<endl;
	  }
     //GMSH test_out<<"$ENDELE"<<endl;
      test_out.clear();
	  test_out.close();    
   }
   part_out<<"#STOP "<<endl;  
   part_out.clear();
   part_out.close();

   //
   delete ele_dom;
   //delete nod_dom;
}


// Test realse
extern void Realese()
{
			long i;
			for(i=0; i<(long)NodesVector.size(); i++)
     delete NodesVector[i];
			NodesVector.clear();
			for(i=0; i<(long)EdgeVector.size(); i++)
     delete EdgeVector[i];
			EdgeVector.clear();
			for(i=0; i<(long)SurfaceFaces.size(); i++)
     delete SurfaceFaces[i];
			SurfaceFaces.clear();

			for(i=0; i<(long)ElementsVector.size(); i++)
     delete ElementsVector[i];
			ElementsVector.clear();
}

