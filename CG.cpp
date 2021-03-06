#include "Parameter.h"
#include "Fiber.h"
#include "Lattice.h"
#include "CG.h"
// #include <blitz/array.h>
#include <vector>

CG::CG(Parameter* parameter)
{
  if(parameter->get_LatticeType()==1){Z=6;}
  Lx=parameter->get_Lx();
  Ly=parameter->get_Ly();
  Lz=parameter->get_Lz();
  D=parameter->get_D();
  eps=parameter->get_Epsilon();
  gammax=parameter->get_Gammax();
  gammay=parameter->get_Gammay();
  k=parameter->get_k();
  Geometry=parameter->get_Geometry();
}
void CG::set_Gamma(double gamx,double gamy)
{
  gammax=gamx;
  gammay=gamy;
}

void CG::Evolv(Lattice* lattice, Fiber* fiber)
{
  Ham ham;
  /*------------------Put all data from originel vector to VecDoub--------------*/
  //VecDoub Position(D*Lx*Ly*Lz);

  
    
  /*---------------------------------------------------------------------------*/
  /*------------------------------Perform the minimization---------------------*/
  
  /*--------------------------------Set the parameter of ham-------------------*/
  // ham.fixe.resize(Lx,Ly,Lz);
  // ham.fixe=lattice->get_fixe();
  //ham.Length.resize(Z/2*Lx*Ly);
  ham.lien.resize(2*Lx*Ly);
  ham.lattice=lattice;
  ham.Lx=Lx;
  ham.Ly=Ly;
  ham.Lz=Lz;
  //ham.eps=eps;
  //for(int i=0;i<Lx;i++){for(int j=0;j<Ly;j++){for(int k=0;k<3;k++){ham.Length[3*Lx*j+3*i+k]=Fiber->get_Length(3*i+k,j,0);}}}
  ham.D=D;
  ham.gammax=gammax;
  ham.gammay=gammay;
  ham.k=k;
  ham.Z=Z;
  ham.size[0]=Lx*(1+gammax)*eps;
  ham.size[1]=Ly*(1+gammay)*0.866*eps;
  int NDOF(0);
  for(int j=0;j<Ly;j++)
    {	
      int pair(0);
      if(j%2==0){pair=1;}
      for(int i =0;i<Lx;i++)
	{
	  // if(i==0 || i==Lx-1 || j==0 || j==Ly-1)
	  //   {
	  //     if(fiber->get_Dilute(Z/2*i,j,0)==1)
	  // 	{
	  // 	  ham.Bound.push_back(i);
	  // 	  ham.Bound.push_back(j);
	  // 	  ham.Bound.push_back(ham.r(i+pair-1,Lx));
	  // 	  ham.Bound.push_back(ham.r(j+1,Ly));
	  // 	}
	  //     if(fiber->get_Dilute(Z/2*i+1,j,0)==1)
	  // 	{
	  // 	  ham.Bound.push_back(i);
	  // 	  ham.Bound.push_back(j);
	  // 	  ham.Bound.push_back(ham.r(i+pair,Lx));
	  // 	  ham.Bound.push_back(ham.r(j+1,Ly));
	  // 	}
	  //     if(fiber->get_Dilute(Z/2*i+2,j,0)==1)
	  // 	{
	  // 	  ham.Bound.push_back(i);
	  // 	  ham.Bound.push_back(j);
	  // 	  ham.Bound.push_back(ham.r(i+1,Lx));
	  // 	  ham.Bound.push_back(j);
	  // 	}
	  //   }
	  // else
	  //   {
	      if(lattice->get_fixe(i,j,0)==0)
		{
		  if(fiber->get_Dilute(Z/2*i,j,0)==1 && lattice->get_fixe(i+pair-1,j+1,0)==1 )
		    {
		      ham.fixe.push_back(i);
		      ham.fixe.push_back(j);
		      ham.fixe.push_back(i+pair-1);
		      ham.fixe.push_back(j+1);
		      ham.Lengthfixe.push_back(fiber->get_Length(Z/2*i,j,0));
		    }
		  if(fiber->get_Dilute(Z/2*i+1,j,0)==1 && lattice->get_fixe(i+pair,j+1,0)==1)
		    {			
		      ham.fixe.push_back(i);
		      ham.fixe.push_back(j);
		      ham.fixe.push_back(i+pair);
		      ham.fixe.push_back(j+1);
		      ham.Lengthfixe.push_back(fiber->get_Length(Z/2*i+1,j,0));
		    }
		  if(fiber->get_Dilute(Z/2*i+2,j,0)==1 && lattice->get_fixe(i+1,j,0)==1)
		    {
		      ham.fixe.push_back(i);
		      ham.fixe.push_back(j);
		      ham.fixe.push_back(i+1);
		      ham.fixe.push_back(j);
		      ham.Lengthfixe.push_back(fiber->get_Length(Z/2*i+2,j,0));
		    }
		}
	      else
		{
		  NDOF++;
		  if(fiber->get_Dilute(Z/2*i,j,0)==1 )
		    {
		      if(lattice->get_fixe(i+pair-1,j+1,0)!=0)
			{
			  ham.DOF.push_back(i);
			  ham.DOF.push_back(j);
			  ham.DOF.push_back(i+pair-1);
			  ham.DOF.push_back(j+1);
			  ham.Length.push_back(fiber->get_Length(Z/2*i,j,0));
			}
		      else
			{
			  ham.rfixe.push_back(i);
			  ham.rfixe.push_back(j);
			  ham.rfixe.push_back(i+pair-1);
			  ham.rfixe.push_back(j+1);
			  ham.Lengthrfixe.push_back(fiber->get_Length(Z/2*i,j,0));
			}
		    }
		  if(fiber->get_Dilute(Z/2*i+1,j,0)==1)
		    {
		      if(lattice->get_fixe(i+pair,j+1,0)!=0)
			{
			  ham.DOF.push_back(i);
			  ham.DOF.push_back(j);
			  ham.DOF.push_back(i+pair);
			  ham.DOF.push_back(j+1);
			  ham.Length.push_back(fiber->get_Length(Z/2*i+1,j,0));
			}
		      else
			{
			  ham.rfixe.push_back(i);
			  ham.rfixe.push_back(j);
			  ham.rfixe.push_back(i+pair);
			  ham.rfixe.push_back(j+1);
			  ham.Lengthrfixe.push_back(fiber->get_Length(Z/2*i+1,j,0));
			}
		    }
		  if(fiber->get_Dilute(Z/2*i+2,j,0)==1)
		    {
		      if(lattice->get_fixe(i+1,j,0)!=0)
			{
			  ham.DOF.push_back(i);
			  ham.DOF.push_back(j);
			  ham.DOF.push_back(i+1);
			  ham.DOF.push_back(j);
			  ham.Length.push_back(fiber->get_Length(Z/2*i+2,j,0));
			}
		      else
			{
			  ham.rfixe.push_back(i);
			  ham.rfixe.push_back(j);
			  ham.rfixe.push_back(i+1);
			  ham.rfixe.push_back(j);
			  ham.Lengthrfixe.push_back(fiber->get_Length(Z/2*i+2,j,0));
			}
		    }
		}
	      //}
	}
    }
  cout<<"number of degree of freedom : "<<(ham.DOF.size())/4<<endl;
  cout<<"number of fixed nodes : "<<(ham.rfixe.size()+ham.fixe.size())/4<<endl;
  cout<<"number of Length : "<<(ham.Length.size())<<endl;
  cout<<"number of Lengthfixe : "<<(ham.Lengthfixe.size()+ham.Lengthrfixe.size())<<endl;

  VecDoub Position(2*NDOF);
  cout<<"NDOF="<<NDOF<<endl;
  int index(0);
  for(int i=0;i<Lx;i++)
    {
      for(int j=0; j<Ly;j++)
	{
	  for(int k=0;k<Lz;k++)
	    {
	      // if(lattice->get_fixe(i,j,0)==1)
	      // 	{
	      //Position[2*i+2*j*Lx+2*k*Ly*Lx]=lattice->get_Grid_rot(i,j,k,0);
		  //Position[2*i+j*2*Lx+k*2*Ly*Lx+1]=lattice->get_Grid_rot(i,j,k,1);
		  //if(D>2){Position[2*i+2*j*Lx+2*k*Ly*Lx+2]=lattice->get_Grid_rot(i,j,k,2);}
	       // 	}
	       // else
	       // 	{
	       // 	  Position[2*i+2*j*Lx]=0.;
	       // 	  Position[2*i+2*j*Lx+1]=0.;
	       // 	  cout<<i<<" "<<j<<endl;
	       // 	}
	      if(lattice->get_fixe(i,j,0)==1)
		{
		  //if(i!=0 || i!=Lx-1 || j!=0 || j!=Ly-1)
		  //  {
		  //cout<<2*NDOF<<" "<<index<<endl;
		      Position[index]=lattice->get_Grid_rot(i,j,k,0);
		      Position[index+1]=lattice->get_Grid_rot(i,j,k,1);
		      ham.lien[2*i+2*Lx*j]=index;
		      ham.lien[2*i+2*Lx*j+1]=index+1;
		      index+=2;
		      //  }
		}
	    }
	}
    }
  // ofstream aight("ouai.res");
  // for(int h=0;h<ham.DOF.size()/4;h++)
  //   {
  //     aight<<DOF(4*h
  //   }
  // aight.close();
  //cout<<Position[0]<<" "<<Position[1]<<" "<<ham.lien[2*9+2*Lx*44]<<" "<<ham.lien[2*9+2*Lx*44+1]<< endl;
  /*--------------------------------------------------------------------------*/
  cout<<"nrg before "<<ham(Position)<<endl;
  Frprmn<Ham> frprmn(ham);
  Position=frprmn.minimize(Position);
  //delete ham.dilute;
  cout<<"nrg after "<<ham(Position)<<endl;
  /*---------------------------------------------------------------------------*/
  /*-------------------------re-transfert the data back------------------------*/
  for(int i=0;i<Lx;i++)
    {
      for(int j=0; j<Ly;j++)
	{
	  if(lattice->get_fixe(i,j,0)==1)
	     {
	       lattice->set_Grid(i,j,0,0,Position[ham.lien[2*i+2*j*Lx]]);
	       lattice->set_Grid(i,j,0,1,Position[ham.lien[2*i+2*j*Lx+1]]);
	      }
	   //else
	     //{
	      //cout<<lattice->get_Grid(i,j,0,0)<<" "<<lattice->get_Grid(i,j,0,1)<<endl;
	      //}
	}
    }
  /*---------------------------------------------------------------------------*/
}
