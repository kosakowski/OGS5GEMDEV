/*========================================================================
 GeoSys - class Matrix (Definition)
          class vec  
 Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation. 
 Function:   See the definition below
 programming:
  22/08/2004  WW  
==========================================================================*/
/// Matrix
#include <iomanip>

#include "matrix_class.h"
#ifdef NewSparseMatrix
#include "msh_mesh.h"
#endif

namespace Math_Group{
// Constructors
Matrix::Matrix(const int rows, const int cols)
{
    
    if(rows*cols>0)
    {
      Sym = false;
      nrows = rows;
      ncols = cols;
      nrows0 = rows;
      ncols0 = ncols;
      size = nrows*ncols;
      data = new double[size];
      for(int i=0; i<size; i++) data[i] = 0.0;
    }
}
Matrix::Matrix()
{
     Sym = false;
     nrows = 0;
     ncols = 0;
     nrows0 = 0;
     ncols0 = 0;
     size = 0;
     data = 0;
}
Matrix::Matrix(const Matrix& m)
{
	Sym = m.Sym;
	nrows = m.nrows;
	ncols = m.ncols;
	nrows0 = m.nrows0;
	ncols0 = m.ncols0;
	size = m.size;
    data = new double[size];
    for(int i=0; i<size; i++) data[i] = 0.0;
}

void Matrix::resize(const int rows, const int cols)
{
 
   if(size>0) { 
      delete [] data;
      data = NULL;
   }
     
   if(rows*cols>0)
   {
      Sym = false;
      nrows = rows;
      ncols = cols;
      nrows0 = rows;
      ncols0 = ncols;
      size = nrows*ncols;
      data = new double[size];
      for(int i=0; i<size; i++) data[i] = 0.0;
   }
}

Matrix::~Matrix()
{
    delete [] data;
    data = NULL;
}

//----------------------------------------------
#ifdef OverLoadNEW_DELETE
void* Matrix::operator new(size_t sz) {
  //printf("operator new: %d Bytes\n", sz);
  void* m = malloc(sz);
  if(!m) puts("out of memory");
  return m;
}

void Matrix::operator delete(void* m) {

  Matrix* mm = static_cast<Matrix*>(m); 
  free(mm); 
}
#endif
//----------------------------------------------
//
void Matrix::operator = (const double a)
{
    for(int i=0; i<size; i++) data[i] = a;
}
void Matrix::operator *= (const double a)
{
    for(int i=0; i<size; i++) data[i] *= a;
}
void Matrix::operator += (const double a)
{
    for(int i=0; i<size; i++) data[i] += a;
}
//
void Matrix::operator = (const Matrix& m)
{
#ifdef gDEBUG    
    if(nrows!=m.Rows()||ncols!=m.Cols())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
#endif
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
         data[i*ncols+j] = m(i,j);
}

//
void Matrix::operator += (const Matrix& m)
{
#ifdef gDEBUG    
   if(nrows!=m.Rows()||ncols!=m.Cols())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
#endif
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
         data[i*ncols+j] += m(i,j);
}

//
void Matrix::operator -= (const Matrix& m)
{
#ifdef gDEBUG    
    if(nrows!=m.Rows()||ncols!=m.Cols()) //Assertion, will be removed
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
#endif
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
         data[i*ncols+j] -= m(i,j);
}
//
void Matrix::GetTranspose(Matrix& m)
{
 #ifdef gDEBUG    
    if(ncols!=m.Rows()&&nrows!=m.Cols())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
 #endif

   for(int i=0; i<m.Rows(); i++)
       for(int j=0; j<m.Cols(); j++)
       {
//          m(i,j) = data[j*ncols+i];
           m(i,j) = (*this)(j,i);
       }

}
//
// m_results = this*m. m_results must be initialized
void Matrix::multi(const Matrix& m, Matrix& m_result, const double fac)
{
 #ifdef gDEBUG    
    if(ncols!=m.Rows()&&nrows!=m_result.Rows()&&m.Cols()!=m_result.Cols())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
 #endif
    for(int i=0; i<m_result.Rows(); i++) 
    {
       for(int j=0; j<m_result.Cols(); j++) 
       { 
           if(Sym&&(j>i)) continue;
           // m_result(i,j) = 0.0;
           for(int k=0; k<ncols; k++)
//            m_result(i,j) += fac*data[i*ncols+k]*m(k,j);
              m_result(i,j) += fac*(*this)(i,k)*m(k,j);
           
       }
    }
}

//
// m_results = this*m1*m2. m_results must be  initialized
void Matrix::multi(const Matrix& m1, const Matrix& m2, Matrix& m_result)
{
  #ifdef gDEBUG    
    if(ncols!=m1.Rows()&&m1.Cols()!=m2.Rows()
       &&m2.Cols()!=m_result.Cols()&&nrows!=m_result.Rows())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
 #endif
    for(int i=0; i<m_result.Rows(); i++) 
    {
       for(int j=0; j<m_result.Cols(); j++) 
       { 
           if(Sym&&(j>i)) continue;
           //m_result(i,j) = 0.0;
           for(int k=0; k<ncols; k++)
           {
              for(int l=0; l<m2.Rows(); l++) 
//                m_result(i,j) += data[i*ncols+k]*m1(k,l)*m2(l,j);
                m_result(i,j) += (*this)(i,k)*m1(k,l)*m2(l,j);
           }
       }
    }
}
// vec_result = This*vec. vec_result must be  initialized
void Matrix::multi(const double *vec, double *vec_result, const double fac)
{
    for(int i=0; i<nrows; i++)
	{
       for(int j=0; j<ncols; j++)
//         vec_result[i] += fac*data[i*ncols+j]*vec[j];
         vec_result[i] += fac*(*this)(i,j)*vec[j];
    }
}

double& Matrix::operator() (const int i, const int j) const
{
 #ifdef gDEBUG    
    if(i>=nrows||j>=ncols)
    {
        cout<<"\n Index exceeds the size of the matrix"<<endl;
        abort();
    }
 #endif
    return data[i*ncols+j]; 
} 
void  Matrix::LimitSize(const int nRows, const int nCols)
{
 #ifdef gDEBUG    
    if(nRows>nrows0||nCols>ncols0)
    {
        cout<<"\n Given size exceeds the original size of the matrix"<<endl;
        abort();
    }
 #endif
    nrows = nRows;
	ncols = nCols;
    size = nrows*ncols;
}


/**************************************************************************
MathLib-Method: 
Task: 
Programing:
08/2004 WW Implementation
02/2005 WW Change name
**************************************************************************/
void Matrix::Write(ostream& os)
{
    //os<<"============================================="<<endl;
    //os<<"Rows: "<<Rows()<<"  Columns: "<<Cols()<<endl;
 
    os.setf(ios::scientific,ios::floatfield);
    os.precision(12);

    for(int i=0; i<nrows; i++)
	{
       os<< "| ";
       for(int j=0; j<ncols; j++)
         os<<(*this)(i,j)<<" ";
       os<< "| "<<endl;
    }
    os<<endl;
    //os<<"============================================="<<endl;
    //os<<endl;     
}

/**************************************************************************
MathLib-Method: 
Task: 
Programing:
01/2006 WW Implementation
**************************************************************************/
void Matrix::Write_BIN(fstream& os)
{
    for(int i=0; i<size; i++) 
      os.write((char*)(&data[i]), sizeof(data[i]));
}
/**************************************************************************
MathLib-Method: 
Task: 
Programing:
01/2006 WW Implementation
**************************************************************************/
void Matrix::Read_BIN(fstream& is)
{
    for(int i=0; i<size; i++) 
      is.read((char*)(&data[i]), sizeof(data[i]));
}


//-----------------------------------------------------
// Symmetrical matrix 
SymMatrix::SymMatrix(const int dim):Matrix(0)
{
    Sym = true; 
    nrows = ncols = dim;    
    size = (int)nrows*(nrows+1)/2;
    data = new double[size];
    nrows0 = ncols0 = dim;
    for(int i=0; i<size; i++) data[i] = 0.0;
}

SymMatrix::SymMatrix():Matrix(0)
{
     Sym = true;
     nrows = 0;
     ncols = 0;
     nrows0 = 0;
     ncols0 = 0;
     size = 0;
     data = 0;
}
SymMatrix::SymMatrix(const SymMatrix& m):Matrix(0)
{
	Sym = m.Sym;
	nrows = m.nrows;
	ncols = m.ncols;
	nrows0 = m.nrows0;
	ncols0 = m.ncols0;
	size = m.size;
    data = new double[size];
    for(int i=0; i<size; i++) data[i] = 0.0;
}

void SymMatrix::resize(const int dim)
{
 
   if(size>0)
   {
	  delete [] data;
      data = NULL;
   }
     
    Sym = true; 
    nrows = ncols = dim;    
    size = (int)nrows*(nrows+1)/2;
    data = new double[size];
    nrows0 = ncols0 = dim;
    for(int i=0; i<size; i++) data[i] = 0.0;
}


//----------------------------------------------
#ifdef OverLoadNEW_DELETE
void* SymMatrix::operator new(size_t sz) {
  //printf("operator new: %d Bytes\n", sz);
  void* m = malloc(sz);
  if(!m) puts("out of memory");
  return m;
}
#endif


void SymMatrix::operator = (const double a)
{
    for(int i=0; i<size; i++) data[i] = a;
}
void SymMatrix::operator *= (const double a)
{
    for(int i=0; i<size; i++) data[i] *= a;
}
void SymMatrix::operator += (const double a)
{
    for(int i=0; i<size; i++) data[i] += a;
}


//
void SymMatrix::operator = (const SymMatrix& m)
{
#ifdef gDEBUG    
    if(nrows!=m.Rows()||ncols!=m.Cols())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
#endif
    int id=0;
    for(int i=0; i<nrows; i++)
    {
       for(int j=0; j<ncols; j++)
       {
          if(j>i) continue;
          id = (int)i*(i+1)/2+j; // temporary  
          data[id] = m(i,j);
       }
    }
}

//
void SymMatrix::operator += (const SymMatrix& m)
{
#ifdef gDEBUG    
   if(nrows!=m.Rows())
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
#endif
    int id=0;
    for(int i=0; i<nrows; i++)
    {
       for(int j=0; j<ncols; j++)
       {
          if(j>i) continue;
          id = (int)i*(i+1)/2+j; // temporary  
          data[id] += m(i,j);
       }
    }
}

//
void SymMatrix::operator -= (const SymMatrix& m)
{
#ifdef gDEBUG    
    if(nrows!=m.Rows()) //Assertion, will be removed
    {
        cout<<"\n The sizes of the two matrices are not matched"<<endl;
        abort();
    }
#endif
    int id=0;
    for(int i=0; i<nrows; i++)
    {
       for(int j=0; j<ncols; j++)
       {
          if(j>i) continue;
          id = (int)i*(i+1)/2+j; // temporary  
          data[id] -= m(i,j);
       }
    }

}
//
double& SymMatrix::operator() (const int i, const int j) const
{
 #ifdef gDEBUG    
    if(i>=nrows||j>=nrows)
    {
        cout<<"\n Index exceeds the size of the matrix"<<endl;
        abort();
    }
 #endif
    
    int id=0;
    if(i>=j)
       id = (int)i*(i+1)/2+j; // temporary
    else
       id = (int)j*(j+1)/2+i; // temporary
    return data[id]; 
} 

void  SymMatrix::LimitSize(const int dim)
{
 #ifdef gDEBUG    
    if(dim>nrows0)
    {
        cout<<"\n Given size exceeds the original size of the matrix"<<endl;
        abort();
    }
 #endif
    nrows = ncols = dim;
    size = (int)nrows*(nrows+1)/2;
}



/*========================================================================
MathLib-Method: 
Task:       Carry out vector operation 
Function:   See the declaration below
programming:
 05/2005  WW  
==========================================================================*/
//1.
template<class T>  vec<T>:: vec(const int argSize):size(argSize)
{		
   entry = new T[argSize];
 #ifdef gDEBUG    
    if (!entry)  {
       cout << "\n*** failed in memory allocatiing for vec "; 
       abort();
    }
#endif
}	
template<class T>  vec<T>:: vec(const vec<T>& v)
{		
   size = v.Size();
   resize(size);

 #ifdef gDEBUG    
    if (!entry)  {
       cout << "\n*** failed in memory allocatiing for vec "; 
       abort();
    }
#endif
   for (int i=0; i<size; i++)
	  entry[i] = v.entry[i];   

}


template<class T> vec<T>:: ~vec()
{		
   delete[] entry;
   entry = 0;
}	

template<class T>  void vec<T>:: resize(const int argSize)
{    
   if(size>0)
   {
	   delete[] entry;
       entry = NULL; 
   }
   size = argSize;
   entry = new T[argSize];
#ifdef gDEBUG    
    if (!entry)  {
       cout << "\n*** failed in memory allocatiing for vec "; 
       abort();
    }
#endif
}	

template<class T> void vec<T>:: operator = (const vec<T>& v)
{
#ifdef gDEBUG    
    if (size!=v.Size())  {
       cout << "\n*** Sizes do not match in vec "; 
       abort();
    }
#endif
    for (int i=0; i<size; i++) entry[i] = v[i];
}

template<class T> void vec<T>:: Write(ostream& os) const
{
   for (int i=0; i<size; i++)
       os<< entry[i]<<"  ";
   os<<endl;
}


//2.
vec<void*>:: vec (const int argSize):size(argSize)
{		
   entry = new void*[argSize];
 #ifdef gDEBUG    
    if (!entry)  {
       cout << "\n*** failed in memory allocatiing for vec "; 
       abort();
    }
#endif
}

vec<void*>:: vec (const vec<void*>& v)
{		
   size = v.Size();
   resize(size);

 #ifdef gDEBUG    
    if (!entry)  {
       cout << "\n*** failed in memory allocatiing for vec "; 
       abort();
    }
#endif
   for (int i=0; i<size; i++)
	  entry[i] = v.entry[i];   

}

vec<void*>::~vec()
{
   delete[] entry;
   entry = 0;
}
void vec<void*>:: resize(const int argSize)
{
    if(size>0)
	{
        delete[] entry;
        entry=NULL;
	}
    size = argSize;
    entry = new void*[argSize];
#ifdef gDEBUG    
    if (!entry)  {
       cout << "\n*** failed in memory allocatiing for vec "; 
       abort();
    }
#endif
}

void vec<void*>:: Write(ostream& os) const
{
   for (int i=0; i<size; i++)
       os<< entry[i]<<"  ";
   os<<endl;
}

void vec<void*>:: operator = (const vec<void*>& v)
{
#ifdef gDEBUG    
    if (size!=v.Size())  {
       cout << "\n*** Sizes do not match in vec "; 
       abort();
    }
#endif
	for (int i=0; i<size; i++) entry[i] = v.entry[i];
}

//3.
template<class T>  void vec<T*>:: operator = (const vec<T*>& v)
{
#ifdef gDEBUG    
    if (size!=v.Size())  {
       cout << "\n*** Sizes do not match in vec "; 
       abort();
    }
#endif
	for (int i=0; i<size; i++) entry[i] = v.entry[i];
}

////////////////////////////////////////////////////////////
//#define NewSparseMatrix
#ifdef NewSparseMatrix
/*\!
   Create sparse matrix table
   01/2006 WW
*/
/*
SparseTable::SparseTable(CFEMesh *a_mesh, bool symmetry)
{
         
}
*/
CSparseMatrix::CSparseMatrix(CFEMesh *a_mesh)
{
     
   

}





// 06/2006 WW
double& CSparseMatrix::operator() (const long i, const long j) const
{
 #ifdef gDEBUG    
    if(i>=Size_row_index||j>=Size_row_index)
    {
        cout<<"\n Index exceeds the dimension of the matrix"<<endl;
        abort();
    }
 #endif    
    long ii, jj, k, new_id, count;
    ii = i;
    jj = j;
    if(symmetry)
    {
       if(ii>jj)
       {
          k = ii;
          ii = jj;
          jj = k; 
       }       
    }
    new_id = row_index_o2new[ii];
    for (k = 0; k < max_column; k++)
    {
       count += new_id;
       if(entry_column[count]==jj)
          return entry[count];  // Found the entry  
	   count += column_size[k]; 
	}
    return entry[count]; 
} 

// 06/2006 WW
void CSparseMatrix::multiVec(double *vec_aug, const double fac)
{
    long i, k, ii, jj, counter;
    for(i=0; i<dimension; i++)
    {
       vec_buffer[i] = fac*vec_aug[i];
       vec_aug[i] = 0.0;
    }
    counter=0;
    for (k = 0; k < max_column; k++)
    {
       for (i = 0; i < column_size[k]; i++)
       {          
          ii = row_index[i];  
          jj=entry_column[counter];
		  vec_aug[ii] += entry[counter]*vec_buffer[jj];
          if(symmetry&&(ii!=jj))
             vec_aug[jj] += entry[counter]*vec_buffer[ii];
		  counter++;
       }         
    }
}
#endif
///////////////////////////////////////////////////////////


}// Namespace

using Math_Group::vec;
using Math_Group::SymMatrix;

template class vec<int>;       
template class vec<long>;       
template class vec<double>;	

// End of class Matrix
//==========================================================================
