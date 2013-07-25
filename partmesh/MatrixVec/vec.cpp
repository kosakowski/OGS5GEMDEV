/*========================================================================
 GeoSys - class my_vector (Declaration)
 Task:       Carry out vector operation
 Function:   See the declaration below
 programming:
 05/2005  WW
==========================================================================*/
#include "vec.h"

namespace Math_Group
{
using namespace std;

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
   if (!entry)
   {
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
   if (!entry)
   {
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
   if (!entry)
   {
      cout << "\n*** failed in memory allocatiing for vec ";
      abort();
   }
#endif
}

template<class T> void vec<T>:: operator = (const vec<T>& v)
{
#ifdef gDEBUG
   if (size!=v.Size())
   {
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
   if (!entry)
   {
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
   if (!entry)
   {
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
   if (!entry)
   {
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
   if (size!=v.Size())
   {
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
   if (size!=v.Size())
   {
      cout << "\n*** Sizes do not match in vec ";
      abort();
   }
#endif
   for (int i=0; i<size; i++) entry[i] = v.entry[i];
}



}// Namespace

using Math_Group::vec;

template class vec<int>;
template class vec<long>;
template class vec<double>;

// End of class Matrix
//==========================================================================
