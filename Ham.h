#ifndef Ham_H
#define Ham_H

#include "Frprmn.h"
#include <vector>
#include "Lattice.h"




struct Ham
{
  

  int Lx,Ly,Lz,D,Z;
  double gammax,gammay,k;
  std::vector<double> Lengthfixe,Lengthrfixe,Length;
  double size[2];
  std::vector<int> DOF,fixe,rfixe;
  //std::vector<int> Bound;
  std::vector<int> lien;
  bool* dilute;
  Lattice* lattice;
  
  int r(int a, int b)
  {
    if(b<0){cout<<"Jsuis deçu..."<<endl; exit(0);}
    while(a<0){a+=b;}
    //if(a!=a || b!=b){cout<<"r"<<endl; exit(0);}
    return a%b;
  }
  
  double r(double a, double b)
  {
    /* if(abs(a)<abs(a-b) && abs(a)<abs(a+b)){return a;} */
    /* if(abs(a-b)<abs(a) && abs(a-b)<abs(a+b)){return a-b;} */
    /* else{return a+b;} */
    return a;
  }
  
  double dist(int i, int j, int l,int m, VecDoub_I &x)
  {
    double distance(0);

    distance+=pow(x[lien[2*i+j*2*Lx]]-x[lien[2*l+m*2*Lx]],2);
    distance+=pow(x[lien[2*i+j*2*Lx+1]]-x[lien[2*l+m*2*Lx+1]],2);
    
    /* if(lien[2*i+j*2*Lx]==0){cout<<"dist 1"<<i<<" "<<j<<endl;} */
    /* if(lien[2*i+j*2*Lx+1]==0){cout<<"dist 2"<<i<<" "<<j<<endl;} */
    /* if(lien[2*l+m*2*Lx]==0){cout<<"dist 3"<<l<<" "<<m<<endl;} */
    /* if(lien[2*l+m*2*Lx+1]==0){cout<<"dist 4"<<l<<" "<<m<<endl;} */
    
    //if(distance!=distance){cout<<"dist is NaN "<<x[2*l+m*2*Lx]<<" "<<x[2*l+m*2*Lx+1] <<" "<<l<< " "<<m <<endl; exit(0);}
    //if(i==0 || i==Lx){distance+=Lx*eps*gamma;}
    return sqrt(distance);
    
  }
  /* double fdist(int l, int m, int i, int j, VecDoub_I &x) */
  /* { */
  /*   double distance(0); */
  /*   distance+=pow(lattice->get_Grid_rot(i,j,0,0)-x[2*l+m*2*Lx],2); */
  /*   distance+=pow(lattice->get_Grid_rot(i,j,0,1)-x[2*l+m*2*Lx+1],2); */
  /*   return sqrt(distance);  */
  /* } */
  double fdist(int i, int j, int l,int m, VecDoub_I &x)
  {
    double distance(0);
    distance+=pow(lattice->get_Grid_rot(i,j,0,0)-x[lien[2*l+m*2*Lx]],2);
    distance+=pow(lattice->get_Grid_rot(i,j,0,1)-x[lien[2*l+m*2*Lx+1]],2);
    /* if(lien[2*l+m*2*Lx]==0){cout<<"fdist 1"<<l<<" "<<m<<endl;} */
    /* if(lien[2*l+m*2*Lx+1]==0){cout<<"fdist 2"<<l<<" "<<m<<endl;} */

    return sqrt(distance);
  }
  double rfdist(int i, int j, int l,int m, VecDoub_I &x)
  {
    double distance(0);
    distance+=pow(x[lien[2*i+j*2*Lx]]-lattice->get_Grid_rot(l,m,0,0),2);
    distance+=pow(x[lien[2*i+j*2*Lx+1]]-lattice->get_Grid_rot(l,m,0,1),2);

    
    /* if(lien[2*i+j*2*Lx]==0){cout<<"rfdist 1"<<i<<" "<<j<<endl;} */
    /* if(lien[2*i+j*2*Lx+1]==0){cout<<"rfdist 2"<<i<<" "<<j<<endl;} */
    return sqrt(distance);
  }
  /* double rfdist(int l, int m, int i,int j, VecDoub_I &x) */
  /* { */
  /*   double distance(0); */
  /*   distance+=pow(x[2*i+j*2*Lx]-lattice->get_Grid_rot(l,m,0,0),2); */
  /*   distance+=pow(x[2*i+j*2*Lx+1]-lattice->get_Grid_rot(l,m,0,1),2); */
  /*   return sqrt(distance); */
  /* } */
  double rdist(int i, int j, int l, int m, VecDoub_I &x)
  {
        double distance(0);
	distance+=pow(r(x[lien[2*i+j*2*Lx]]-x[lien[2*l+m*2*Lx]],size[0]),2);
	distance+=pow(r(x[lien[2*i+j*2*Lx+1]]-x[lien[2*l+m*2*Lx+1]],size[1]),2);
	return sqrt(distance);	
  }

  
  Doub operator() (VecDoub_I &x)
  {
    /*
    Doub Hamilt(0);
    for(int j=0;j<Ly;j++)
      {
	int pair(0);
	if(j%2==0){pair=1;}
	for (int i=0;i<Lx;i++)
	  {
	    int p1,p2,p3,p4,p5;;
	    p1=r(i+pair,Lx);
	    p2=r(j+1,Ly);
	    p3=r(i+1,Lx);
	    p5=r(i+pair-1,Lx);
	    if(dilute[Z/2*i+Z/2*j*Ly]==1){Hamilt+=0.5*k*pow(dist(i,j,p1,p2,x)-eps,2);}
	    if(dilute[Z/2*i+1+Z/2*j*Ly]==1){Hamilt+=0.5*k*pow(dist(i,j,p5,p2,x)-eps,2);}
	    if(dilute[Z/2*i+2+Z/2*j*Ly+1]==1){Hamilt+=0.5*k*pow(dist(i,j,p3,j,x)-eps,2);}
	  }
      }

    return Hamilt;
    */
    
    Doub Hamilt(0);
    //#pragma omp parallel for
    for(int i=0;i<DOF.size()/4;i++)
      {
	Hamilt+=0.5*k*pow(dist(DOF[4*i],DOF[4*i+1],DOF[4*i+2],DOF[4*i+3],x)-Length[i],2);
      }
    /* for(int i=0;i<Bound.size()/4;i++) */
    /*   { */
    /* 	//cout<<"ntm"<<endl; */
    /* 	Hamilt+=0.5*k*pow(rdist(Bound[4*i],Bound[4*i+1],Bound[4*i+2],Bound[4*i+3],x)-eps,2); */
    /*   } */
    for(int i=0;i<fixe.size()/4;i++)
      {
    	Hamilt+=0.5*k*pow(fdist(fixe[4*i],fixe[4*i+1],fixe[4*i+2],fixe[4*i+3],x)-Lengthfixe[i],2);
      }
    for(int i=0;i<rfixe.size()/4;i++)
      {
    	Hamilt+=0.5*k*pow(rfdist(rfixe[4*i],rfixe[4*i+1],rfixe[4*i+2],rfixe[4*i+3],x)-Lengthrfixe[i],2);
      }
    return Hamilt;
    
  }

/* double Hamilt(VecDoub_I &x) */
/*   { */
/*     double Hamilt(0); */


/* #pragma omp parallel for */
/*     for(int j=0;j<Ly;j++) */
/*       { */
/* 	int pair(0); */
/* 	if(j%2==0){pair=1;} */
/* 	for (int i=0;i<Lx;i++) */
/* 	  { */
/* 	    int p1,p2,p3,p4,p5;; */
/* 	    p1=r(i+pair,Lx); */
/* 	    p2=r(j+1,Ly); */
/* 	    p3=r(i+1,Lx); */
/* 	    p5=r(i+pair-1,Lx); */


/* 	    Hamilt+=0.5*k*pow(dist(i,j,p1,p2,x)-eps,2); */
/* 	    Hamilt+=0.5*k*pow(dist(i,j,p5,p2,x)-eps,2);	  */
/* 	    Hamilt+=0.5*k*pow(dist(i,j,p3,j,x)-eps,2); */
/* 	  } */
/*       } */
/*     return Hamilt; */

/*   } */

  /*--------------------ANCIENNE DF----------------------*/
  /*  
  void df(VecDoub_I &x, VecDoub_O &deriv)
  {
#pragma omp parallel for
    for( int j=0;j<Ly;j++)
      {
  	int pair(0);
  	if(j%2==0){pair=1;}
  	for( int i=0; i<Lx;i++)
  	  {
  	    int p1,p2,p3,p5,p6,p7;
  	    p1=r(i+pair,Lx);
  	    p2=r(j+1,Ly);
  	    p3=r(i+1,Lx);
  	    p5=r(i+pair-1,Lx);
  	    p6=r(j-1,Ly);
  	    p7=r(i-1,Lx);

  	    deriv[2*i+2*j*Lx]=0;
	    deriv[2*i+2*j*Lx+1]=0;
  	    double a(0);
	    if(dilute[Z/2*i+1+Z/2*j*Ly]==1)
	      {
		a=dist(i,j,p1,p2,x);
		if(a!=0)
		  {
		    deriv[2*i+j*2*Lx]+=k*(a-eps)*r(x[2*i+j*2*Lx]-x[2*p1+Lx*2*p2],size[0])/a;
		    deriv[2*i+j*2*Lx+1]+=k*(a-eps)*r(x[2*i+j*2*Lx+1]-x[2*p1+Lx*2*p2+1],size[1])/a;
		  }
	      }
	    if(dilute[Z/2*i+Z/2*j*Ly]==1)
	      {
		a=dist(i,j,p5,p2,x);
		if(a!=0)
		  {
		    deriv[2*i+j*2*Lx]+=k*(a-eps)*r(x[2*i+j*2*Lx]-x[2*p5+p2*2*Lx],size[0])/a;
		    deriv[2*i+j*2*Lx+1]+=k*(a-eps)*r(x[2*i+j*2*Lx+1]-x[2*p5+p2*2*Lx+1],size[1])/a;
		  }
	      }
	    if(dilute[Z/2*i+2+Z/2*j*Ly]==1)
	      {
		a=dist(i,j,p3,j,x);
		if(a!=0)
		  {
		    deriv[2*i+j*2*Lx]+=k*(a-eps)*r(x[2*i+j*2*Lx]-x[2*p3+j*2*Lx],size[0])/a;
		    deriv[2*i+j*2*Lx+1]+=k*(a-eps)*r(x[2*i+j*2*Lx+1]-x[2*p3+j*2*Lx+1],size[1])/a;
		  }
	      }
	    if(dilute[Z/2*p1+Z/2*p6*Ly] ==1)
	      {
		a=dist(i,j,p1,p6,x);
		if(a!=0)
		  {
		    deriv[2*i+j*2*Lx]+=k*(a-eps)*r(x[2*i+j*2*Lx]-x[2*p1+p6*2*Lx],size[0])/a;
		    deriv[2*i+j*2*Lx+1]+=k*(a-eps)*r(x[2*i+j*2*Lx+1]-x[2*p1+p6*2*Lx+1],size[1])/a;
		  }
	      }
	    if(dilute[Z/2*p5+1+Z/2*p6*Ly] ==1)
	      {
		a=dist(i,j,p5,p6,x);
		if(a!=0)
		  {
		    deriv[2*i+j*2*Lx]+=k*(a-eps)*r(x[2*i+j*2*Lx]-x[2*p5+p6*2*Lx],size[0])/a;
		    deriv[2*i+j*2*Lx+1]+=k*(a-eps)*r(x[2*i+j*2*Lx+1]-x[2*p5+p6*2*Lx+1],size[1])/a;
		  }
	      }
	    if(dilute[Z/2*p7+2+Z/2*j*Ly] ==1)
	      {
		a=dist(i,j,p7,j,x);
		if(a!=0)
		  {
		    deriv[2*i+j*2*Lx]+=k*(a-eps)*r(x[2*i+j*2*Lx]-x[2*p7+j*2*Lx],size[0])/a;
		    deriv[2*i+j*2*Lx+1]+=k*(a-eps)*r(x[2*i+j*2*Lx+1]-x[2*p7+j*2*Lx+1],size[1])/a;
		  }
	      }
  	  }
      }

  }
  */

  /*--------------------------------------------------------------------------*/
  
/* void df(VecDoub_I &x, VecDoub_O &deriv) */
/*   { */
/*     VecDoub xdx(deriv.size()); */
/*     double dx=0.001; */
/*     xdx=deriv; */
/*     //#pragma omp parallel for */
/*     for(int j=0; j<Ly;j++) */
/*       { */
/* 	for( int i=0;i<Lx;i++) */
/* 	  { */
/* 	    xdx[2*i+2*j*Lx]+=dx; */
/* 	    deriv[2*i+2*j*Lx]=(Hamilt(xdx)-Hamilt(x))/dx; */
/* 	    xdx[2*i+2*j*Lx]-=dx; */
/* 	    xdx[2*i+2*j*Ly+1]+=dx; */
/* 	    deriv[2*i+2*j*Lx+1]=(Hamilt(xdx)-Hamilt(x))/dx; */
/* 	    xdx[2*i+2*j*Ly+1]-=dx; */
/* 	  } */
/*       } */
/*   } */


  /*------------------------------------nouvelle def------------------------------*/
  
void df(VecDoub_I &x, VecDoub_O &deriv)
  {
    double a(0);
    int i,j,l,m;
    for(int h=0;h<deriv.size();h++){deriv[h]=0;}
    for(int h=0;h<DOF.size()/4;h++)
      {
	i=DOF[4*h];
	j=DOF[4*h+1];
	l=DOF[4*h+2];
	m=DOF[4*h+3];
	a=dist(i,j,l,m,x);
	if(a!=0)
	  {	      
	    deriv[lien[2*i+j*2*Lx]]+=k*(a-Length[h])*(x[lien[2*i+j*2*Lx]]-x[lien[2*l+m*2*Lx]])/a;
	    deriv[lien[2*i+j*2*Lx+1]]+=k*(a-Length[h])*(x[lien[2*i+j*2*Lx+1]]-x[lien[2*l+m*2*Lx+1]])/a; 
	    deriv[lien[2*l+m*2*Lx]]+=k*(a-Length[h])*(x[lien[2*l+m*2*Lx]]-x[lien[2*i+j*2*Lx]])/a;
	    deriv[lien[2*l+m*2*Lx+1]]+=k*(a-Length[h])*(x[lien[2*l+m*2*Lx+1]]-x[lien[2*i+j*2*Lx+1]])/a;
	    /* if(lien[2*i+j*2*Lx]==0){cout<<"DOF 1"<<i<<" "<<j<<endl;} */
	    /* if(lien[2*i+j*2*Lx+1]==0){cout<<"DOF 2"<<i<<" "<<j<<endl;} */
	    /* if(lien[2*l+m*2*Lx]==0){cout<<"DOF 3"<<l<<" "<<m<<endl;} */
	    /* if(lien[2*l+m*2*Lx+1]==0){cout<<"DOF 4"<<l<<" "<<m<<endl;} */
	  }
      }
    for(int h=0;h<fixe.size()/4;h++)
      {

	/* l=fixe[4*h]; */
	/* m=fixe[4*h+1]; */
	/* i=fixe[4*h+2]; */
	/* j=fixe[4*h+3]; */
	i=fixe[4*h];
	j=fixe[4*h+1];
	l=fixe[4*h+2];
	m=fixe[4*h+3];
	a=fdist(i,j,l,m,x);
	if(a!=0)
	  {
	    deriv[lien[2*l+m*2*Lx]]+=k*(a-Lengthfixe[h])*(x[lien[2*l+m*2*Lx]]-lattice->get_Grid_rot(i,j,0,0))/a;
	    deriv[lien[2*l+m*2*Lx+1]]+=k*(a-Lengthfixe[h])*(x[lien[2*l+m*2*Lx+1]]-lattice->get_Grid_rot(i,j,0,1))/a;

	    /* if(lien[2*l+m*2*Lx]==0){cout<<"fixe 1 "<<l<<" "<<m<<endl;} */
	    /* if(lien[2*l+m*2*Lx+1]==0){cout<<"fixe 2 "<<l<<" "<<m<<endl;} */
	    /* deriv[2*i+j*2*Lx]=1/0.; */
	    /* deriv[2*i+2*Lx*j+1]=1/0.; */
	  }
      }
    for(int h=0;h<rfixe.size()/4;h++)
      {
	/* l=rfixe[4*h]; */
	/* m=rfixe[4*h+1]; */
	/* i=rfixe[4*h+2]; */
	/* j=rfixe[4*h+3]; */
	i=rfixe[4*h];
	j=rfixe[4*h+1];
	l=rfixe[4*h+2];
	m=rfixe[4*h+3];
	a=rfdist(i,j,l,m,x);
	if(a!=0)
	  {
	    deriv[lien[2*i+j*2*Lx]]+=k*(a-Lengthrfixe[h])*(x[lien[2*i+j*2*Lx]]-lattice->get_Grid_rot(l,m,0,0))/a;
	    deriv[lien[2*i+j*2*Lx+1]]+=k*(a-Lengthrfixe[h])*(x[lien[2*i+j*2*Lx+1]]-lattice->get_Grid_rot(l,m,0,1))/a;	    
	    /* if(lien[2*i+j*2*Lx]==0){cout<<"rfixe 1"<<i<<" "<<j<<endl;} */
	    /* if(lien[2*i+j*2*Lx+1]==0){cout<<"rfixe 2"<<i<<" "<<j<<endl;} */
	    /* deriv[2*l+m*2*Lx]=1/0.; */
	    /* deriv[2*l+m*2*Lx+1]=1/0.; */
	  }
      }
    /* for(int h=0;h<Bound.size()/4;h++) */
    /*   { */
    /* 	i=Bound[4*h]; */
    /* 	j=Bound[4*h+1]; */
    /* 	l=Bound[4*h+2]; */
    /* 	m=Bound[4*h+3]; */
    /* 	a=rdist(i,j,l,m,x); */
    /* 	if(a!=0) */
    /* 	  {	       */
    /* 	    deriv[lien[2*i+j*2*Lx]]+=k*(a-eps)*r(x[lien[2*i+j*2*Lx]]-x[lien[2*l+m*2*Lx]],size[0])/a; */
    /* 	    deriv[lien[2*i+j*2*Lx+1]]+=k*(a-eps)*r(x[lien[2*i+j*2*Lx+1]]-x[lien[2*l+m*2*Lx+1]],size[1])/a;  */
    /* 	    deriv[lien[2*l+m*2*Lx]]+=k*(a-eps)*r(x[lien[2*l+m*2*Lx]]-x[lien[2*i+j*2*Lx]],size[0])/a; */
    /* 	    deriv[lien[2*l+m*2*Lx+1]]+=k*(a-eps)*r(x[lien[2*l+m*2*Lx+1]]-x[lien[2*i+j*2*Lx+1]],size[1])/a; */
    /* 	  } */
    /*   } */
  }

  /*--------------------------------------------------------------------------------*/
};

#endif
