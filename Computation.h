//
// Computation.h
// Simple_String_Network
//
// Created by AAAAAAAAAAAAAAAIIIIIIIIIIIIIIIGGGGGGGGGGGHHHHHHHHHHTTTTTTT

#ifndef Computation_h
#define Computation_h

#include <vector>
#include <iostream>
#include <fstream>

class Parameter;
class Lattice;
class Fiber;



class Computation
{
 public:
  Computation(Parameter* parameter);
  void Time_Evolv(Parameter* parameter, Lattice* lattice, Fiber* fibre);
  //double PickTime();
  std::vector<int> PickFiber(Parameter* parameter,Fiber* fiber,bool cutable);
 private:
  int NFiber;
  std::ofstream test;

};

#endif
