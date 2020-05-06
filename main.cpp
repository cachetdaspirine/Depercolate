//
//  main.cpp
//  simple_string_network_Sim
//
//Created by AAAAAAAAAAAIIIIIIIIIIIIGGGGGGHHHHHHTTTTT!!!!


/*----------------------include technique-----------------------*/
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <random>
/*--------------------------------------------------------------*/

/*---------------------include output/input---------------------*/
#include "Parameter.h"
#include "Computation.h"
#include <time.h>
/*--------------------------------------------------------------*/

/*---------------------include selectif-------------------------*/
#include "Lattice.h"
#include "Fiber.h"

/*--------------------------------------------------------------*/

using namespace::std;


void OnlyNetwork(Parameter* parameter, Computation* compute, string arg1)
{
  /*---------------------Boucle de la serie----------------------*/
  /*------------------loop over different series-----------------*/
  for (int S=0; S<parameter->get_SNum();S++)
    {
      /*----------gestion des fichiers d entre/sortie------------*/
      /*---------------manage the input/output file--------------*/
      string MkFoldSerie;
      MkFoldSerie="mkdir Res/simulation"+arg1+"/S"+to_string(S);
      system(MkFoldSerie.c_str());
      
      parameter->set_PathSaveLattice("Res/simulation"+arg1+"/S"+to_string(S)+"/Lattice/");
      system(("mkdir Res/simulation"+arg1+"/S"+to_string(S)+"/Lattice").c_str());
      
      parameter->set_PathSaveStress("Res/simulation"+arg1+"/S"+to_string(S)+"/Stress/");
      system(("mkdir Res/simulation"+arg1+"/S"+to_string(S)+"/Stress").c_str());
      
      
      parameter->set_PathSaveFiber("Res/simulation"+arg1+"/S"+to_string(S) +"/Fiber");
      system(("mkdir Res/simulation"+arg1+"/S"+to_string(S)+"/Fiber").c_str());
      
      parameter->set_PathSaveForce("Res/simulation"+arg1+"/S"+to_string(S));
	
      
      /*---------------------------------------------------------*/
      /*--------------------allocate pointer---------------------*/
      Lattice* lattice=new Lattice(parameter);
      lattice->Square();
      Fiber* fiber=new Fiber(parameter, lattice);
      /*---------------------------------------------------------*/
      /*------------------Compute the main loop------------------*/
      compute->Time_Evolv(parameter,lattice,fiber);	  
      /*--------------------delete pointeur---------------------*/
      delete lattice;
      delete fiber;
      /*--------------------------------------------------------*/
      /*--ajuste la valeur des parametre pour la serie suivant--*/
      /*--adjust the value of the parameter for the next series-*/     
      parameter->AdjustParam();
    }
}
 
int main (int argc, char* argv[])
{
  float temps;
  clock_t t1,t2;
  t1=clock();
  Parameter* parameter = new Parameter();
  Computation* compute = new Computation(parameter);

  
  
  /*-----------gestion des parametres d un serie---------------*/
  string arg1(argv[1]);
  system(("rm -r Res/simulation"+arg1).c_str());
  system(("mkdir Res/simulation"+arg1).c_str());
  system(("cp data.in Res/simulation"+arg1).c_str());

  /*-----------------------------------------------------------*/

  
  /*--------lecture des parametre d une simulation-----------*/

  srand(parameter->get_seed());
  

  OnlyNetwork(parameter,compute,arg1);

  delete parameter;
  delete compute;
  t2=clock();
  temps=(float) (t2-t1)/CLOCKS_PER_SEC;
  cout<<"time of execution: "<<temps<<endl;
  return 0;
}
