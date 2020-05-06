//
// Parameter.hpp
// Simple_String_Network
// Created by AAAAAAAAAAAAAAAAIIIIIIIIIIIIGGGGGGGGGHHHHHHHHHHTTTTTTTTTTT

#ifndef Parameter_h
#define Parameter_h

#include <iostream>
#include <blitz/array.h>
#include <vector>
#include <random>


class Parameter
{
 public :
  
  Parameter();
  bool get_Restore();
  int flag();
  double get_W0();
  int get_mobile();
  int get_rot();
  int get_CG();
  int get_SNum();
  int get_LatticeType(); //1=triangulaire
  int get_Lx();
  int get_Ly();
  int get_Lz();
  int get_seed();
  int get_NCrackLim();
  int get_NTot();
  int get_D();
  int get_Dump();
  int get_Boundary();
  double get_TimeStep();
  int get_dDump();
  int get_nEvolvStress();
  double get_Njump();
  double get_Ncut();
  double get_ddGamma();
  double get_k();
  double get_Epsilon();
  double get_Gammax();
  double get_Gammay();
  double get_dGammax();
  double get_dGammay();
  double get_P();
  double get_EpsilonLim();
  double get_PCrack();
  double get_HoleSize();
  double g_Tau();
  double g_alpha();
  blitz::Array<int,1> get_HolePos();
  std::string get_Geometry();
  void AdjustParam();
  void set_Gamma();
  std::string get_PathSaveFiber();
  std::string get_PathSaveLattice();
  std::string get_PathSaveStress();
  std::string get_PathSaveForce();
  void set_PathSaveForce(std::string path);
  void set_PathSaveFiber(std::string path);
  void set_PathSaveLattice(std::string path);
  void set_PathSaveStress(std::string path);
  std::default_random_engine* g_generator();
  
  
 private :
  /*-----------------------Parametre de serie-----------------------------*/
  std::string PathSaveFiber;
  std::string PathSaveLattice;
  std::string PathSaveStress;
  std::string PathSaveForce;
  std::default_random_engine* generator;
  bool Restore;
  double Ncut,Njump;
  int rot,mobile;
  int Dump,Boundary;
  int CG;
  int D,drapeaux;
  int SNum;
  int dnEvolvStress;
  double W0;
  double dNTot;
  double ddGamma;
  double dHoleSize;
  blitz::Array<int,1> HolePos;
  double dP;
  int  dseed;
    
  /*----------------------------------------------------------------------*/

  /*----------------Parametres d initialisation d une simulation----------*/
  int LatticeType;
  int NTot;
  double TimeStep;
  int dDump;
  int nEvolvStress;
  int seed;
  /*----------------------------------------------------------------------*/

  /*------------------------Parametre du reseaux--------------------------*/
  int Lx,Ly,Lz;
  double k;
  double Epsilon;
  double Gammax,Gammay;
  double dGammax,dGammay;
  double P;
  double EpsilonLim;
  double PCrack;
  int NCrackLim;
  double HoleSize;
  double Tau,alpha;
  std::string Geometry;
  /*----------------------------------------------------------------------*/
};

#endif
