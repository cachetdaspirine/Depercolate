#ifndef CG_H
#define CG_H

#include "nr3.h"
#include "Ham.h"
/* #include "Parameter.h" */
/* #include "Fiber.h" */
/* #include "Lattice.h" */


class Fiber;
class Lattice;
class Parameter;


class CG
{
 public:
  CG(Parameter* parameter);
  void Evolv(Lattice* lattice, Fiber* fiber);
  void set_Gamma(double gamx,double gamy);
  
 private:
  int Z,Lx,Ly,Lz,D;
  double eps,gammax,gammay,k;
  string Geometry;
  
  
};


#endif
