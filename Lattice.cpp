//
//  Lattice.cpp
//  Simple_String_Lattice
//
// Created by AAAAAAAAAAAAAAAIIIIIIIIIIIIIIIGGGGGGGGGGGHHHHHHHHHHTTTTTTT
#include "Lattice.h"
/*-----------------------------include classique-----------------------*/
#include <iostream>
#include <cmath>
#include <fstream>
#include <blitz/array.h>
#include <random>
/*---------------------------------------------------------------------*/

/*------------------------------include d objet------------------------*/
#include "Parameter.h"
#include "Fiber.h"
/*---------------------------------------------------------------------*/
using namespace blitz;
// {{{ constructor
Lattice::Lattice(Parameter* parameter)
{
  Lx=parameter->get_Lx();
  Ly=parameter->get_Ly();
  Lz=parameter->get_Lz();
  D=parameter->get_D();
  dt=parameter->get_TimeStep();
  eps=parameter->get_Epsilon();
  path=parameter->get_PathSaveLattice();
  Gammay=parameter->get_Gammay();
  Gammax=parameter->get_Gammax();
  rot=parameter->get_rot();
  /*----------------Cree un vecteur de la bonne dimension--------------*/
  grid.resize(Lx,Ly,Lz,3);
  coord.resize(Lx,Ly,Lz);
  coord=0;
  fixe.resize(Lx,Ly,Lz);
  busy.resize(Lx,Ly,Lz);
  busy=0;
  /*-------------------------------------------------------------------*/
  // {{{ Geometry
  if(parameter->get_LatticeType()==1){Triangular(parameter);}
  fixe=1;
  // }}}  
}

// }}}

// {{{ Geometry

void Lattice::Geometrie()
{
  /*-------------------creer un reseau de la bonne geometrie-----------*/

  double gammax(0),gammay(0),gammaxx(0),gammaxy(0);
  gammax=Gammay*0.25882;
  gammay=Gammay*0.96593;
  gammaxx=Gammax*0.96593;
  gammaxy=Gammax*0.25882;
  cout<<"circle"<<endl;
  for(int i=0; i<Lx;i++)
    {
      for(int j=0;j<Ly;j++)
	{
	  for(int k=0;k<Lz;k++)
	    {
	      if(dist(i,j,k,Lx/2,Ly/2,Lz/2)<0.866*eps*(Ly/2.-1) && dist(i,j,k,Lx/2,Ly/2,Lz/2)>=0.866*eps*(Ly/2.-2))
		{
		      
		  fixe(i,j,k)=0;
		  //		  grid(i,j,k,0)=grid(i,j,k,0)+Gammax*Lx;
		  //grid(i,j,k,1)=grid(i,j,k,1)+Gammay*Ly;
		  
		  grid(i,j,k,0)=grid(i,j,k,0)+gammax*(get_Grid(i,j,k,1))+gammaxx*get_Grid(i,j,k,0);
		  grid(i,j,k,1)=grid(i,j,k,1)+gammay*(get_Grid(i,j,k,1))+gammaxy*get_Grid(i,j,k,0);
		}
	    }
	}
    }   
}
void Lattice::Square()
{
  /*-------------------creer un reseau de la bonne geometrie-----------*/

  //double gammax(0),gammay(0),gammaxx(0),gammaxy(0);
  // gammax=Gammay*0.25882;
  // gammay=Gammay*0.96593;
  // gammaxx=Gammax*0.96593;
  // gammaxy=Gammax*0.25882;
  cout<<"square"<<endl;
  int F(0);
  for(int i=0; i<Lx;i++)
    {
      for(int j=0;j<Ly;j++)
	{
	  for(int k=0;k<Lz;k++)
	    {
	      if(i==Lx-2 )
	      	{
	      	  //
		  //fixe(i,j,k)=0;
		  //grid(i,j,k,1)=grid(i,j,k,1)+Gammax*Lx;
		  
	      	  // grid(i,j,k,1)=grid(i,j,k,1)+Gammay*(get_Grid(i,j,k,1));
	      	  //grid(i,j,k,0)=grid(i,j,k,0)+Gammax*Lx;

	      	  //grid(i,j,k,0)=grid(i,j,k,0)+gammax*(get_Grid(i,j,k,1))+gammaxx*get_Grid(i,j,k,0);
	      	  //grid(i,j,k,1)=grid(i,j,k,1)+gammay*(get_Grid(i,j,k,1))+gammaxy*get_Grid(i,j,k,0);
	      	}
	      else if(i==1)
	      	{
	      	  //fixe(i,j,k)=0;
	      	  //grid(i,j,k,0)=grid(i,j,k,0)-Gammax*Lx;
	      	}
	      if(j==Ly-1)
		{
		  grid(i,j,k,0)=grid(i,j,k,0)-Lx*Gammay;
		  grid(i,j,k,1)=grid(i,j,k,1)+Lx*Gammax;
		  fixe(i,j,k)=0;
		}
	      else if(j==0)
		{
		  grid(i,j,k,0)=grid(i,j,k,0)+Gammay*Lx;
		  grid(i,j,k,1)=grid(i,j,k,1)-Gammax*Lx;
		  fixe(i,j,k)=0;
		}
	    }
	}
    }
}
/*-------------------------------------------------------------------*/

// }}}

// {{{ Triangle
/*---------------------------Lattice triangulaire--------------------*/
void Lattice::Triangular(Parameter* parameter)
{     
  double Xpair(0), Ypair(0);
  for(int i=0;i<Lx;i++)
    {
      for(int j=0;j<Ly;j++)
  	{
  	  if(j%2==0){Xpair=1;}
  	  else{Xpair=0;}
  	  for(int k=0;k<Lz;k++)
  	    {
  	      if(k%2==0){Ypair=1;}
  	      else{Ypair=0;}
		  
  	      grid(i,j,k,0)=i*parameter->get_Epsilon()+Xpair/2.;
  	      grid(i,j,k,1)=j*parameter->get_Epsilon()*0.866+Ypair/2.;
  	      if(D>2)
  		{
  		  grid(i,j,k,2)=k*parameter->get_Epsilon()*0.866;
  		}
  	    }
  	}
    } 
  
  
}
/*-------------------------------------------------------------------*/
// }}}

// {{{ Evolution

/*---------------------Evolv Lattice (no Cell)-----------------------*/
void Lattice::Evolv(Fiber* fiber)
{
 
  for(int i=0;i<Lx;i++)
    {
      for(int j=0; j<Ly;j++)
	{
	  for(int k=0; k<Lz;k++)
	    {
	      for(int d=0;d<D;d++)
		{
		  grid(i,j,k,d)=grid(i,j,k,d)+fiber->get_Force(i,j,k,d)*dt*fixe(i,j,k);
		}
	    }
	}
    }
}
/*-------------------------------------------------------------------*/
// void Lattice::Evolv_Relax(Fiber* fiber)
// {
//   vector<double> project(D);
//   project[0]=0.965926;
//   project[1]=0.258819;
//   if(D==3){project[2]=0;}   
//   for(int i=0;i<Lx;i++)
//     {
//       for(int j=0; j<Ly;j++)
// 	{
// 	  for(int k=0; k<Lz;k++)
// 	    {
// 	      for(int d=0;d<D;d++)
// 		{
// 		  grid(i,j,k,d)=grid(i,j,k,d)+10*fiber->get_Force(i,j,k,d)*dt*project[d];
// 		}
// 	    }
// 	}
//     }

// }
/*---------------------Evolv Lattice (with Cell)-----------------------*/
/*-------------------------------------------------------------------*/

// }}}

// {{{ get truc

/*--------------------------get truc---------------------------------*/
Array<int,3> Lattice::get_coord()
{
  return coord;
}
int Lattice::get_coord(int i, int j, int k)
{
  return coord(i,j,k);
}
double Lattice::get_Grid(int i, int j, int k, int l)
{
  if(rot==0)
    {
      return grid(i,j,k,l);
    }
  else
    {
      if(l==0)
	{
	  return (grid(i,j,k,0)-Lx/2.*eps)*0.965926-(grid(i,j,k,1)-Ly/2*eps*0.866)*0.258819;
	}
      if(l==1)
	{
	  //return grid(i,j,k,l);
	  return (grid(i,j,k,0)-Lx/2.*eps)*0.258819+(grid(i,j,k,1)-Ly/2.*eps*0.866)*0.965926;
	}
      if(l==2)
	{
	  return grid(i,j,k,l);
	}
    }
}
void Lattice::set_Grid(Array<double,4> CopGrid)
{
  grid=CopGrid;
}
double Lattice::get_Grid_rot(int i, int j, int k, int l)
{
  return grid(i,j,k,l);
}
int Lattice::get_fixe(int i,int j,int k)
{
  return fixe(i,j,k);
}
Array<int,3> Lattice::get_fixe()
{
  return fixe;
}
void Lattice::set_fixe(int i, int j, int k,int fix)
{
  fixe(i,j,k)=fix;
}
Array<double,4> Lattice::get_Grid()
{
  return grid;
}
/*-------------------------------------------------------------------*/

// }}}

// {{{ Print

/*----------------------Print et Close-------------------------------*/
void Lattice::Print( int n)
{
  ofstream output((path+"network"+to_string(n)+".res").c_str());
  for(int i=0;i<Lx;i++)
    {
      for(int j=0;j<Ly;j++)
	{
	  for(int k=0;k<Lz;k++)
	    {
	      for(int l=0;l<D;l++)
		{
		  output<<get_Grid(i,j,k,l)<<" ";
		}
	      output<<i<<" "<<j<<" "<<k<<endl;
	    }
	}
    }
  output.close();
}
/*-------------------------------------------------------------------*/

// }}}

// {{{ coordination
void Lattice::p_coord(int i, int j, int k)
{
  coord(i,j,k)++;
}
void Lattice::m_coord(int i, int j, int k)
{
  coord(i,j,k)--;
}
void Lattice::set_coord(int i, int j, int k, int cord)
{
  coord(i,j,k)=cord;
}
void Lattice::set_Grid(int i, int j, int k, int l, double pos)
{
  grid(i,j,k,l)=pos;
}
// }}}

// {{{ La vie qu on mène
void Lattice::Pull(int i, int j, int k,int l, double dl)
{
  grid(i,j,k,l)=grid(i,j,k,l)+dl;
}
double Lattice::dist(int i, int j, int k, int l,int m, int n)
{
  double distance(0);

  for(int o=0;o<D;o++)
    {
      distance+=pow(grid(i,j,k,o)-grid(l,m,n,o),2);
    }
    
  distance=sqrt(distance);
  return distance;
}
void Lattice::set_busy(int i, int j, int k, int aight)
{
  busy(i,j,k)=aight;
}
int Lattice::get_busy(int i, int j, int k)
{
  return busy(i,j,k);
}
// }}}
