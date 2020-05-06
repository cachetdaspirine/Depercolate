//
//  Lattice.h
//  Simple_String_Lattice
//
// Created by AAAAAAAAAAAAAAAIIIIIIIIIIIIIIIGGGGGGGGGGGHHHHHHHHHHTTTTTTT


#ifndef Lattice_h
#define Lattice_h
/*-----------------------------include classique-----------------------*/
#include <iostream>
#include <cmath>
#include <fstream>
#include <blitz/array.h>
#include <vector>
/*---------------------------------------------------------------------*/

/*------------------------------include d objet------------------------*/
#include "Parameter.h"
#include "Fiber.h"
/*---------------------------------------------------------------------*/

class Fiber;


using namespace blitz;

class Lattice
{    
public:
  Lattice(Parameter* parameter);
  
  void Print( int n);
  void Evolv( Fiber* fiber);
  /* void Evolv_Relax(Fiber* fiber); */

  void Geometrie();
  void Square();
  int get_coord(int i, int j, int k);
  Array<int,3> get_coord();
  void p_coord(int i, int j, int k);
  void m_coord(int i, int j, int k);
  void set_coord(int i, int j, int k, int cord);
  double get_Grid(int i, int j, int k,int l);
  double get_Grid_rot(int i, int j, int k,int l);
  void set_Grid(int i, int j, int k, int l, double pos);
  void Pull(int i, int j, int k,int l, double dl);
  void set_busy(int i, int j, int k, int aight);
  int get_busy(int i, int j, int k);
  void set_fixe(int i, int j, int k, int fix);
  double dist(int i, int j, int k, int l, int m, int n);
  int get_fixe(int i,int j,int k);
  Array<int,3> get_fixe();
  Array<double,4> get_Grid();
  void set_Grid(Array<double,4> CopGrid);
  void Triangular(Parameter* parameter);
private:
  Array<int,3> busy;
  Array<double,4> grid;
  Array<int,3> fixe;
  Array<int,3> coord;
  int D,Lx,Ly,Lz,rot;
  double dt,eps;
  string path;
  double Gammay,Gammax;
  
};
#endif
