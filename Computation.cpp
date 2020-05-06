//
// Computation.cpp
// Simple_String_Network
//
// Created by AAAAAAAAAAAAAAAIIIIIIIIIIIIIIIGGGGGGGGGGGHHHHHHHHHHTTTTTTT

#include "Computation.h"
/*-----------------------------include classique-----------------------*/
#include <iostream>
#include <cmath>
#include <fstream>
/*---------------------------------------------------------------------*/

/*------------------------------include d objet------------------------*/
#include "Parameter.h"
#include "Lattice.h"
#include "Fiber.h"
#include "CG.h"
/*---------------------------------------------------------------------*/


//using namespace blitz;
using namespace std;
Computation::Computation(Parameter* parameter){NFiber=0;test.open("test.res");}

// {{{ Time Evolv no Cell
void Computation::Time_Evolv(Parameter* parameter, Lattice* lattice, Fiber* fiber)
{
  for(int i=0;i<parameter->get_Lx();i++){for(int j=0;j<parameter->get_Ly();j++){for(int o=0;o<3;o++)
      {
	if(fiber->get_Dilute(3*i+o,j,0)){NFiber++;}
      }}}
  // cout<<NFiber<<endl;
  // string trash;
  // cin>>trash;
  double tot_sim_time=parameter->get_NTot();
  //double dt=parameter->get_TimeStep();
  ofstream OutForce;
  OutForce.open((parameter->get_PathSaveForce()+"/MForce.res").c_str(),ios::out|ios::trunc);
  OutForce<<"time n force"<<endl;
  fiber->CreateHole(lattice);
  fiber->OpenStress();
  fiber->PrintStress(0);
  lattice->Print(0);
  double time(0);

  // un ptit coup de gradient conjugue pour le debut

  CG* cg=new CG(parameter);
  
  if(parameter->get_CG()==1)
    {
      cg->set_Gamma(parameter->get_Gammax(),parameter->get_Gammay());
      cg->Evolv(lattice,fiber);
    }
  double F0(0);
  /*-------------------------boucle de temps----------------------------*/
  for(int n=1;n<=parameter->get_NTot();n++)
    {
      parameter->set_Gamma();
      fiber->set_Gamma(parameter);
      cg->set_Gamma(parameter->get_Gammax(), parameter->get_Gammay());
      int ii(0),jj(0),oo,p(0);
      double dt(0);
      vector<int> Index(fiber->get_BreakingFiber(dt));
      //cout<<"dt="<<dt<<endl;
       //-------------------------------------------------------------
      //cut the bond for which time has come to be cut
      int Ncuted(0);
      for(auto& it : Index)
	{
	  jj=it/(3*parameter->get_Lx());
	  //cout<<"Index="<<Index;
	  //Index-=jj*parameter->get_Lx()*3;
	  ii=(it%(3*parameter->get_Lx()))/3;
	  oo=(it%(3*parameter->get_Lx()))%3;
	  //cout<<"ii="<<ii<<"jj="<<jj<<"oo="<<oo<<endl;
	  fiber->cut(ii,jj,0,oo,lattice);
	  test<<ii<<" "<<jj<<" "<<oo<<endl;
	  NFiber--;
	  Ncuted++;
	}

      if(n==1)
	{
	  Array<double,1> Mforce(2);
	  Mforce=0.;
	  for(int i=1;i<parameter->get_Lx()-1;i++)
	    {
	      Mforce+=abs(fiber->Compute(lattice,i,1,0));
	      Mforce+=abs(fiber->Compute(lattice,i,parameter->get_Ly()-1,0));
	    }
	  Mforce=Mforce/(6*(parameter->get_Lx()-2));
	  F0=sqrt(pow(Mforce(0),2)+pow(Mforce(1),2));
	}

	  
      if(n%parameter->get_dDump()==0)
	{
	  cg->Evolv(lattice,fiber);
	  if(parameter->get_Dump()==1)
	    {
	      lattice->Print(n);
	      fiber->PrintStress(n);
	      fiber->Printer(lattice,n);
	      fiber->PrintLength(n);
	    }
	  //---------------------------output--------------------------------

	  //Compute the average force on the boundaries :
	  Array<double,1> Mforce(2);
	  Mforce=0.;
	  for(int i=1;i<parameter->get_Lx()-1;i++)
	    {
	      Mforce+=abs(fiber->Compute(lattice,i,1,0));
	      Mforce+=abs(fiber->Compute(lattice,i,parameter->get_Ly()-1,0));
	    }
	  Mforce=Mforce/(6*(parameter->get_Lx()-2));
	  if(n==1){F0=sqrt(pow(Mforce(0),2)+pow(Mforce(1),2));}
	  OutForce<<time<<" "<<n<<" "<<sqrt(pow(Mforce(0),2)+pow(Mforce(1),2))/F0<<endl;
	}
      //-----------Adjust the time
      time+=dt;
      //cout<<"adjust time"<<endl;
      
      fiber->IterateTime(dt);

      //-------------------------
      for(int Nnew=0;Nnew<Ncuted;Nnew++)
	{
      //cout<<"add a new fiber now"<<endl;
      vector<int> VectIndex=PickFiber(parameter,fiber,false);
      //cout<<"fiber picked"<<endl;
      ii=VectIndex[0];jj=VectIndex[1];oo=VectIndex[2];
      fiber->set_Dilute(3*ii+oo,jj,0,true);
      fiber->NewFiber(ii*3+parameter->get_Lx()*3*jj+oo);
      NFiber++;
      int i1,i2,i3(0),i4,i5,i6(0),l,pair(0);
      //l=ii%3;
      l=oo;
      i1=ii;
      i2=jj;
      if(jj%2==0){pair=1;}
      if(l==0){i4=i1-1+pair;i5=i2+1;}
      if(l==1){i4=i1+pair;i5=i2+1;}
      if(l==2){i4=i1+1;i5=i2;}
      if(i4==parameter->get_Lx()-1 || i4==0 || i5==parameter->get_Ly()-1)
	{
	  cout<<"ii/jj/oo : "<<ii<<"/"<<jj<<"/"<<oo<<endl;
	  cout<<"i1/i2/i4/i5/l"<<i1<<"/"<<i2<<"/"<<i4<<"/"<<i5<<"/"<<l<<"\n";
	  string trash;
	  cin>>trash;
	}
      //cout<<"restore"<<endl;
      //restored with a new restlength
      fiber->set_Length(3*ii+oo,jj,0,sqrt(pow(lattice->get_Grid_rot(i1,i2,i3,0)-lattice->get_Grid_rot(i4,i5,i6,0),2)+pow(lattice->get_Grid_rot(i1,i2,i3,1)-lattice->get_Grid_rot(i4,i5,i6,1),2)));
	}
      //------------------------Time-------------------------------------
      cout<<"time = "<<time<<" n= "<<n<<endl;      

      
    }
  /*--------------------------------------------------------------------*/
  /*-------------------------------Close--------------------------------*/
  fiber->CloseStress();
  OutForce.close();
  
}
// }}}
// double Computation::PickTime()
// {
//   double Rd(static_cast<double>(rand() % 10000)/10000);
//   return log(1/(1-Rd))/static_cast<double>(NFiber);
//   //return (pow((1-Rd),1/(1-2))-1)/static_cast<double>(NFiber);
// }
vector<int> Computation::PickFiber(Parameter* parameter, Fiber* fiber,bool cutable)
{
  // if cutable is true then we look for a bond that we can cut, otherwise we look for a slot for a bond
  vector<int> Index(3);
  int ii,jj,oo,p;
  bool Yes(true);
  if(cutable)
    {
      int count(0);
      while(Yes)
       	{
	  if(count>=10000){cout<<"Too much iteration cannot find a cutable bond, count="<<count<<endl;}
	  ii=rand() % (parameter->get_Lx()-2);
	  jj=rand() % (parameter->get_Ly()-2);
	  ii++;
	  jj++;
	  p=1-(jj%2);
	  oo=rand() % 3;
	  if(fiber->get_Dilute(3*ii+oo,jj,0) & fiber->get_Length(3*ii+oo,jj,0)==1)
	    {
	      if(ii!=parameter->get_Lx()-2 & jj!=parameter->get_Ly()-2 & ii!=1){Yes=false;}
	      if(ii==parameter->get_Lx()-2 & p==1 & oo==0 & jj!=parameter->get_Ly()-2){Yes=false;}
	      if(ii==parameter->get_Lx()-2 & p==0 & oo!=2 & jj!=parameter->get_Ly()-2){Yes=false;}
	      if(jj==parameter->get_Ly()-2 & oo==2 & ii!=parameter->get_Lx()-2){Yes=false;}
	      if(ii==1 & p==1 & jj!=parameter->get_Ly()-2){Yes=false;}
	    }
	  count++;
	}
      Index[0]=ii;
      Index[1]=jj;
      Index[2]=oo;
    }
  else
    {
      int count(0);
      while(Yes)
       	{
	  if(count>=10000){cout<<"Too much iteration cannot find a free bond, count="<<count<<endl;exit(0);}
	  ii=rand() % (parameter->get_Lx()-2);
	  jj=rand() % (parameter->get_Ly()-2);
	  ii++;
	  jj++;
	  p=1-(jj%2);
	  oo=rand() % 3;
	  if(fiber->get_Dilute(3*ii+oo,jj,0)){}
	  else
	    {
	      if(ii!=parameter->get_Lx()-2 & jj!=parameter->get_Ly()-2 & ii!=1){Yes=false;}
	      if(ii==parameter->get_Lx()-2 & p==1 & oo==0 & jj!=parameter->get_Ly()-2){Yes=false;}
	      if(ii==parameter->get_Lx()-2 & p==0 & oo!=2 & jj!=parameter->get_Ly()-2){Yes=false;}
	      if(jj==parameter->get_Ly()-2 & oo==2 & ii!=parameter->get_Lx()-2){Yes=false;}
	      if(ii==1 & p==1 & jj!=parameter->get_Ly()-2){Yes=false;}
	    }
	  count++;	  
	}
      Index[0]=ii;
      Index[1]=jj;
      Index[2]=oo;
      //test<<ii<<" "<<jj<<" "<<oo<<endl;
    }
  return Index;
}
