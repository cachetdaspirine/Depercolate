//
//  Fiber.h
//  Simple_String_Lattice
//
// Created by AAAAAAAAAAAAAAAIIIIIIIIIIIIIIIGGGGGGGGGGGHHHHHHHHHHTTTTTTT

#ifndef Fiber_h
#define Fiber_h
/*-----------------------------include classique-----------------------*/
#include <iostream>
#include <cmath>
#include <fstream>
#include <blitz/array.h>
#include <map>
#include <random>
/*---------------------------------------------------------------------*/

/*------------------------------include d objet------------------------*/
#include "Parameter.h"
#include "Lattice.h"
/*---------------------------------------------------------------------*/


class Lattice;

using namespace blitz;

class Fiber
{    
public:
  Fiber(Parameter* parameter, Lattice* lattice);
  void Print(int n);/*(Array<double,1> E1,Array<double,1> E2,Array<double,1> E3, int i, int j, int k );*/
  void PrintStress(int n);
  void PrintLength(int n);
  void CloseStress();
  /* void CloseFiber(); */
  void OpenStress();
  /* void OpenFiber(int n); */
  void Evolv(Lattice* lattice, int n);
  void Printer(Lattice* lattice,int n);
  Array<double,1> Compute(Lattice* lattice, int i, int j, int k);
  void CheckCoord(Lattice* lattice);
  double get_Force(int i, int j, int k,int l);
  Array<double,4> get_Force();
  int get_Crack();
  bool get_Dilute(int i,int j, int k);
  void set_Dilute(int i, int j, int k, bool dil);
  double get_Length(int i, int j, int k) const;
  void IterateTime(double dt);
  std::vector<int> get_BreakingFiber(double& time);
  void NewFiber(int index);
  void set_Length(int i, int j, int k,double L);
  double Hdist(int i, int j, int k, int l, int m, int n, Lattice* lattice);
  double dist(int i, int j, int k, int l, int m, int n, Lattice* lattice);
  double planedist(int i, int j, int k, int l, int m, int n, Lattice* lattice);
  Array<double,1> vect(int i, int j, int k, int l,int m, int n, Lattice* lattice);
  Array<double,1> Uvect(int i, int j, int k, int l,int m, int n, Lattice* lattice);
  Array<double,1> vect_rot(int i, int j, int k, int l,int m, int n, Lattice* lattice);
  Array<double,1> Uvect_rot(int i, int j, int k, int l,int m, int n, Lattice* lattice);
  void Clean(int i, int j, int k, int l, int m, int n,Lattice* lattice);
  void Clean(int i, int j, int k,Lattice* lattice);
  void rClean(int i, int j, int k,Lattice* lattice);
  void cut(int i, int j, int k, int o,Lattice* lattice);
  void CreateHole(Lattice* lattice);
  Array<int,1> Num(int i, int j, int k, int l, int m, int n);
  void set_Gamma(Parameter* param);
  double PickTime();
  void CloseCheck();

  
 private:

  bool generate();
  std::map<double,std::vector<int>> Time;
  ofstream outStress,output,outLength,CheckDistrib;
  int Hx,Hy,Hz,dDump;
  double RH;
  Parameter* param;
  Array<bool,3> dilute;
  Array<double,4> force;
  Array<double,3> Length;
  int D,Lx,Ly,Lz,Z;
  Array<double,1> size;
  Array<double,4> printer;
  double kel,eps,epsLim,Gammax,Gammay,Pdil,Pcrack;
  int crack,NCrackLim;
  string path,pathStress,Geometry;
  ofstream test;

  std::normal_distribution<double>* Distrib;
};



#endif
