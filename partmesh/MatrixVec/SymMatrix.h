/*========================================================================
 GeoSys - class Matrix (Declaration)
 Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation.
 Function:   See the declaration below
 programming:
  22/08/2004  WW
==========================================================================*/
#ifndef SymMatrix_INC
#define SymMatrix_INC

#include "Matrix.h"

//#define OverLoadNEW_DELETE

namespace Math_Group
{

// Symmetrical matrix. 12-01-2005. WW
class SymMatrix:public Matrix
{
   public:
      SymMatrix(const int dim);
      SymMatrix();
      explicit SymMatrix(const SymMatrix& m);

      void resize(const int dim);

      ~SymMatrix() {}

//----------------------------------------------
#ifdef OverLoadNEW_DELETE
      // Allocate memory
      void* operator new(size_t sz);
#endif
//----------------------------------------------

      // Operators

      void LimitSize(const int dim);

      // Access to members
      double& operator() (const int i, const int j) const;
};

}
//==========================================================================
#endif
