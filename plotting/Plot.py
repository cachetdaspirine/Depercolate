#/usr/bin/ipython

import pylab
from Network import *
from Folder import *
import os

#------------------Liste des donnees d'entre---------------
SimNum=1
repertoir="/data/H.Le/simulation"+str(SimNum)+"/S/network/"
NumFile=len(os.listdir(repertoir))
NumPlot=2

Zoom=2
Border=0

#----------------------------------------------------------
Fold=Folder(Folder=repertoir,ParameterFile="Parameters.log",ResFile="Positions",Ext=".res")
N=Network(filename=Fold.OpenDataFile(NumPlot),time=NumPlot*Fold.dt,sizex=Fold.Lx,sizey=Fold.Ly,sizez=Fold.Lz,distance=Fold.eps)
N.plot2D_fiber(cracking=Fold.PCrack)
N.AdjustWindows(zoom=Zoom,border=Border)
N.show()
