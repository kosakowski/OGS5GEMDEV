/*========================================================================
 GeoSys - class my_vector (Declaration)
 Task:       Carry out vector operation
 Function:   See the declaration below
 programming:
 05/2005  WW
==========================================================================*/
#ifndef vec_INC
#define vec_INC

#include<iostream>


namespace Math_Group
{

template<class T> class vec
{
   public:
      vec(const int argSize);
      vec()
      {
         size=0;
      }
      explicit vec(const vec<T>& v);

      virtual ~vec();
      // Operator
      virtual void operator = (T v)
      {
         for (int i=0; i<size; i++) entry[i] = v;
      }
      virtual void operator = (const vec<T>&);
      virtual void resize(const int newh);
      virtual T& operator[] (int i)
      {
         return (T&) entry[i];
      }
      virtual const T& operator[] (int i) const
      {
         return (const T&) entry[i];
      }
      virtual int Size() const
      {
         return size;
      }

      T* Entry()
      {
         return entry;
      }
      T* Entry()  const
      {
         return entry;
      }

      virtual void Write(std::ostream& os = std::cout) const;
   protected:
      T* entry;
      int size;

};
template<> class vec<void*>
{
   public:
      vec(const int argSize);
      vec()
      {
         size=0;
      }
      explicit vec(const vec<void*>& v);

      virtual ~vec();
      // Operator
      void operator = (void* v)
      {
         for (int i=0; i<size; i++) entry[i] = v;
      }
      void operator = (const vec<void*>& v);
      void*& operator[] (int i)
      {
         return entry[i];
      }
      const void*& operator[] (int i) const
      {
         return (const void*&) entry[i];
      }

      // Access to memebers
      void** Entry()
      {
         return entry;
      }
      const void** Entry()  const
      {
         return (const void**)entry;
      }

      virtual void resize(const int newh);
      virtual int Size() const
      {
         return size;
      }
      virtual void Write(std::ostream& os = std::cout) const;

   protected:
      void** entry;
      int size;
};

template<class T> class vec<T*> : public vec<void*>
{
   public:
      vec(const int Size) : vec<void*>(Size) { }
      vec()               : vec<void*>()     { }
      explicit vec(const vec<T*>& v) : vec<void*>(v) { }

      ~vec() { }

      // Operator
      void operator = (T* v)
      {
         for (int i=0; i<size; i++) entry[i] = v;
      }
      void operator = (const vec<T*>& v);
      T*& operator[] (int i)
      {
         return (T*&) entry[i];
      }
      const T*& operator[] (int i) const
      {
         return (const T*&) entry[i];
      }

      T** Entry()
      {
         return entry;
      }
      T** Entry()  const
      {
         return (const T**)entry;
      }


};

}
//==========================================================================

#endif
