/*========================================================================
 GeoSys - class Matrix (Declaration)
 Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation.
 Function:   See the declaration below
 programming:
  22/08/2004  WW
==========================================================================*/
#ifndef Matrix_INC
#define Matrix__INC

#include<iostream>
#include<fstream>

//#define OverLoadNEW_DELETE

namespace Math_Group
{
class Matrix
{
   public:
      Matrix(const int rows, const int cols=1);
      Matrix();
      explicit Matrix(const Matrix& m);
      //
      void resize(const int rows, const int cols=1);
      //
      virtual ~Matrix();
//----------------------------------------------
#ifdef OverLoadNEW_DELETE
      // Allocate memory
      void* operator new(size_t sz);
      // Release memory
      void operator delete(void* m);
#endif
//----------------------------------------------
      // Operators
      void operator = (const double a);
      void operator *= (const double a);
      void operator += (const double a);
      void operator = (const Matrix& m);
      void operator += (const Matrix& m);
      void operator -= (const Matrix& m);

      void GetTranspose(Matrix& m);

      // vec_result = This*vec. vec_result must be initialized
      void multi(const double *vec, double *vec_result, const double fac=1.0);
      // m_result = this*m. m_result must be initialized
      void multi(const Matrix& m, Matrix& m_result, const double fac=1.0);
      // m_result = this*m1*m2. m_result must be initialized
      void multi(const Matrix& m1, const Matrix& m2, Matrix& m_result);

      // Access to members
      virtual double& operator() (const int i, const int j=0) const;
      void LimitSize(const int nRows, const int nCols=1);

      int Rows() const
      {
         return nrows;
      }
      int Cols() const
      {
         return ncols;
      }
      int Size() const
      {
         return size;
      }

      // Print
      void Write(std::ostream& os = std::cout);
      void Write_BIN(std::fstream& os);
      void Read_BIN(std::fstream& is);
   protected:
      double *data;
      int nrows, nrows0;
      int ncols, ncols0;
      int size;
      bool Sym;
};
typedef Matrix Vec;
}
//==========================================================================
#endif
