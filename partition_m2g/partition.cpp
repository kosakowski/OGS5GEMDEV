// expre_new_Operator.cpp
// compile with: /EHsc
#include <cmath>

#include <stdio.h>
#include <cstdlib>  // For GCC 4.3.2
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <vector>

#include "matrix_class.h"
#include "mesh.h"

#include <time.h>

using namespace std;

using Math_Group::SymMatrix;
using Math_Group::Matrix;
using Math_Group::Vec;
using Math_Group::vec;

int main(int argc, char* argv[])
{
  ifstream rfi_in;
  fstream rfi_out;


  const int MAX_ZEILE = 1028; 
  char str0[MAX_ZEILE];
  char str1[MAX_ZEILE];
  char str2[MAX_ZEILE];
  bool quad = false;
  int aug=0;
  //ios::pos_type position;
 
  cout<<"\t|=====================================| \n"<<endl;
  cout<<"\t|     Grid partition converter        | \n"<<endl;
  cout<<"\t|          (From METIS)               | \n"<<endl;
  cout<<"\t|                                     | \n"<<endl;
  cout<<"\t|    University of Tuebingen (WW@ZAG) | \n"<<endl;
  cout<<"\t|                                     | \n"<<endl;
  cout<<"\t|    File name (no extension)         | \n"<<endl;
  cout<<"\t|    Integer:                         | \n"<<endl;
  cout<<"\t|      0: Linear, do not convert      | \n"<<endl;
  cout<<"\t|      1: Quadratic, do not convert   | \n"<<endl;
  cout<<"\t|      2: Write METIS element data    | \n"<<endl;
  cout<<"\t|     -n: Negative non-zero number    | \n"<<endl;
  cout<<"\t|         Convert! n: number of parts.| \n"<<endl;
  cout<<"\t|=====================================| \n"<<endl;
  cout<<"\tInput file name: ";

  if(argc>1) 
  {
      strcpy(str1,argv[1]);
      strcpy(str2,argv[2]);
  }
  else 
    scanf(" %s  %s%*[^\n]%*c",str1, str2);

  aug = atoi(str2);
  if(aug==1) quad = true;

  strcpy(str0,str1);
  strcpy(str2,str1);
  strcat(str1,".msh");
  strcat(str2,".mesh");

  rfi_in.open(str1);
  if(aug>0) rfi_out.open(str2, ios::out );
  if(!rfi_in.is_open())
  {
      cerr<<("Error: cannot open msh file . It may not exist !");
      abort();
  }

  bool rfiMesh = true;
  string line_string;
  getline(rfi_in,line_string); // The first line
  if(line_string.find("#FEM_MSH")!=string::npos)
    rfiMesh = false;
  if(line_string.find("GeoSys-MSH")!=string::npos) 
    rfiMesh = false;
  rfi_in.seekg(0L,ios::beg);

  if(rfiMesh)
     ReadGrid(rfi_in);
  else
     ReadGridGeoSys(rfi_in);
  clock_t start, finish;
  start = clock();

  ConstructGrid(quad);
  if(aug<0) ConstructDomain(str0, abs(aug));
 
  finish = clock();
  cout<<"CPU time elapsed in deformation process: "
      <<(double)(finish - start) / CLOCKS_PER_SEC<<"s"<<endl;
   
  if(aug==2)
  {  
     rfi_out<<(long)ElementsVector.size()<<" ";
     int e_type =0;
	 switch(ElementsVector[0]->getElementType())    
	 {
        case 1: cout<<"Not for 1D element"<<endl; abort(); 
        case 2: e_type =4; break; 
        case 3: e_type =3; break;  
        case 4: e_type =1; break;  
        case 5: e_type =2; break; 
        case 6: cout<<"Not for prismal element"<<endl; abort(); 
	 } 
     rfi_out<<e_type<<endl;
     for(long i=0; i<(long)ElementsVector.size(); i++)
        ElementsVector[i]->Write_index(rfi_out);
  }

  Realese();
  return 0;

}
