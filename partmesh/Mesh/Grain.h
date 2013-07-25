#ifndef Grain_INC
#define Grain_INC

#include<string>
#include<iostream>

//------------------------------------------------------
//   Topology declartion of geometrical element.
//   WW. 06.2005
//   WW  02.2012
//------------------------------------------------------

namespace Mesh_Group
{
//enum Element_Type  {Line, Quad, Hex, Tri, Tet, Pris };

//1.  Mesh declaration
class Grain
{
   public:
      Grain(const int id);
      virtual  ~Grain() {}
      // Operator
      virtual void operator = (const Grain & g) {}
      virtual bool operator == (const Grain & g)
      {
         return false;
      }

      void Marking(const bool state)
      {
         mark = state;
      }
      // Get members
      long getIndex() const
      {
         return index;
      }
      bool getStatus() const
      {
         return mark;
      }

      // Output
      virtual void output(std::ostream& os = std::cout) const {};
   protected:
      long index;
      // Towards special purpose,
      // e.g. marked to be refined or active
      bool mark;
};

} //end namespace

#endif
