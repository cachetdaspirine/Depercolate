//
// Parameter.hpp
// Simple_String_Network
// Created by AAAAAAAAAAAAAAAAIIIIIIIIIIIIGGGGGGGGGHHHHHHHHHHTTTTTTTTTTT

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <blitz/array.h>

#include "Parameter.h"


using namespace blitz;

Parameter::Parameter()
{
  ifstream input("data.in");
  string Trash;
  
  getline(input,Trash);
  getline(input,Trash);
  getline(input,Trash);

  /*-----------------------Parametre de serie-------------------------*/
  input>>Trash>>D;
  HolePos.resize(D); 
  input>>Trash>>SNum;
  input>>Trash>>dnEvolvStress;
  input>>Trash>>dNTot;
  input>>Trash>>ddGamma;
  input>>Trash>>dHoleSize;
  input>>Trash>>dP;
  input>>Trash>>dseed;
  /*------------------------------------------------------------------*/
  getline(input,Trash);
  getline(input,Trash);
  /*-------------------Parameter d une simulation---------------------*/
  input>>Trash>>LatticeType;
  input>>Trash>>NTot;
  input>>Trash>>TimeStep;
  input>>Trash>>dDump;
  input>>Trash>>nEvolvStress;
  input>>Trash>>seed;
  /*------------------------------------------------------------------*/

  /*---------------------Parametres du reseaux------------------------*/
  input>>Trash>>Lx;
  input>>Trash>>Ly;
  input>>Trash>>Lz; // doit valoir au moins 1
  input>>Trash>>k;
  input>>Trash>>Epsilon;
  input>>Trash>>Gammax;
  input>>Trash>>Gammay;
  input>>Trash>>dGammax;
  input>>Trash>>dGammay;
  input>>Trash>>P;
  input>>Trash>>EpsilonLim;
  input>>Trash>>PCrack;
  input>>Trash>>NCrackLim;
  input>>Trash>>HoleSize;
  input>>Trash>>HolePos(0)>>HolePos(1)>>HolePos(2);
  input>>Trash>>CG;
  input>>Trash>>W0;
  input>>Trash>>Geometry;
  input>>Trash>>Ncut;
  input>>Trash>>rot;
  input>>Trash>>Dump;
  input>>Trash>>drapeaux;
  input>>Trash>>alpha;
  input>>Trash>>Tau;
  input.close();
  generator=new std::default_random_engine(seed);
}
std::default_random_engine* Parameter::g_generator(){
  return generator;
  
}
bool Parameter::get_Restore()
{
  return Restore;
}
int Parameter::get_mobile()
{
  return mobile;
}
int Parameter::get_rot()
{
  return rot;
}

int Parameter::flag()
{
  return drapeaux;
}

// int Parameter::get_Boundary()
// {
//   return Boundary;
// }
int Parameter::get_Dump()
{
  return Dump;
}
double Parameter::get_Njump()
{
  return Njump;
}
double Parameter::get_Ncut()
{
  return Ncut;
}
string Parameter::get_Geometry()
{
  return Geometry;
}
double Parameter::get_W0()
{
  return W0;
}
int Parameter::get_CG()
{
  return CG;
}
double Parameter::get_Gammay()
{
  return Gammay;
}
double Parameter::get_dGammax()
{
  return dGammax;
}
double Parameter::get_dGammay()
{
  return dGammay;
}
Array<int,1> Parameter::get_HolePos()
{
  return HolePos;
}
int Parameter::get_NCrackLim()
{
  return NCrackLim;
}
int Parameter::get_SNum()
{
  return SNum;
}
int Parameter::get_LatticeType()
{
  return LatticeType;
}
int Parameter::get_Lx()
{
  return Lx;
}
std::string Parameter::get_PathSaveForce()
{
  return PathSaveForce;
}
std::string Parameter::get_PathSaveFiber()
{
  return PathSaveFiber;
}
void Parameter::set_PathSaveFiber(std::string path)
{
  PathSaveFiber=path;
}
std::string Parameter::get_PathSaveLattice()
{
  return PathSaveLattice;
}
void Parameter::set_PathSaveForce(std::string path)
{
  PathSaveForce=path;
}
void Parameter::set_PathSaveLattice(std::string path)
{
  PathSaveLattice=path;
}

std::string Parameter::get_PathSaveStress()
{
  return PathSaveStress;
}
void Parameter::set_PathSaveStress(std::string path)
{
  PathSaveStress=path;
}
int Parameter::get_Ly()
{
  return Ly;}
int Parameter::get_Lz()
{
  return Lz;
}
int Parameter::get_seed()
{
  return seed;
}
int Parameter::get_NTot()
{
  return NTot;
}
double Parameter::get_TimeStep()
{
  return TimeStep;
}
int Parameter::get_dDump()
{
  return dDump;
}
int Parameter::get_nEvolvStress()
{
  return nEvolvStress;
}
double Parameter::get_k()
{
  return k;
}
double Parameter::get_Epsilon()
{
  return Epsilon;
}
double Parameter::get_Gammax()
{
  return Gammax;
}
double Parameter::get_P()
{
  return P;
}
double Parameter::get_EpsilonLim()
{
  return EpsilonLim;
}
double Parameter::get_PCrack()
{
  return PCrack;
}
double Parameter::get_HoleSize()
{
  return HoleSize;
}
double Parameter::g_alpha(){return alpha;}
double Parameter::g_Tau(){return Tau;}
int Parameter::get_D()
{
  return D;
}


void Parameter::AdjustParam()
{
  nEvolvStress+=dnEvolvStress;
  NTot +=dNTot;
  dGammax+=ddGamma;
  dGammay+=ddGamma;
  HoleSize+=dHoleSize;
  P+=dP;
  seed+=dseed;
  
}
void Parameter::set_Gamma()
{
  Gammax+=dGammax;
  Gammay+=dGammay;
}
