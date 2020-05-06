#/usr/bin/ipython

import pylab
from Folder import *
from Network import *
import os

time=100
for
SimNum=6
repertoir="/data/H.Le/simulation"+str(SimNum)+"/S/network/"
k=1.1


Fold=Folder(Folder=repertoir,ParameterFile="Parameters.log",ResFile="Positions",Ext=".res")
numero=time/(Fold.dt*Fold.pas)
N=Network(filename=Fold.OpenDataFile(numero),time=numero*Fold.dt,sizex=Fold.Lx,sizey=Fold.Ly,sizez=Fold.Lz,distance=Fold.eps,k=k)
N.tabl()

N.StressProfile(direction=1,hole="false",mean="false",large=0,Nbande=10)
