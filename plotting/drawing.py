#!/usr/bin/python

import pylab
from Network import *
from Folder import *
from Cell import *
import os


show=0
nbeg=0
save=1
tVid=20
SimNum=1
Sbeg=0
Nserie=1
zoom=1.7
follow=0
Quality=5 #qualite de la video en M


for S in range(Sbeg,Sbeg+Nserie):
    os.system("mkdir ../Res/simulation"+str(SimNum)+"/S"+str(S)+"/Image")
    Fold=Folder(Folder="../Res/simulation"+str(SimNum)+"/",Serie="S"+str(S)+"/",ParameterFile="data.in",NodeFile="Lattice/network",FiberFile="Fiber/fiber",Ext=".res")
    n=1
    if(Fold.evolvcell==1):
        CellFold=[Folder(Folder="../Res/simulation"+str(SimNum)+"/",Serie="S"+str(S)+"/", ParameterFile="data.in", CellFile="Cell/cell"+str(i)+".res",VectFile="Cell/force",Ext=".res") for i in range(Fold.Ncell)]
        C=[Cell(CellFold[i]) for i in range (Fold.Ncell)]
        limite=Fold.Ntot
    else:
        limite=Fold.Ntot#*Fold.dDump
    while n<limite :
        print(n)
        if(Fold.evolvcell==1):
            for i in range(Fold.Ncell):
                C[i].UnPlotCell()
        #try :
        N=Network(folder=Fold,N=n)
        #except :
        #    print("Break!")
        #    n+=1
        #    continue
            #break
        N.plot2D_fiber(cracking=Fold.PCrack/100)
        #-----------------------------------a modifier-----------
        if(Fold.evolvcell==1):
        #    try:
            for i in range(Fold.Ncell):
            #C[i].OpenVector(n,i)
                C[i].PlotCell(n)
            #C[i].PlotVect(n)
             #except:
              #  break
        #-------------------------------------------------
        numero=str(n)
        if len(numero)==1 :
            numero="000"+numero
        if len(numero)==2 :
            numero="00"+numero
        if len(numero)==3 :
            numero="0"+numero
        if(Fold.evolvcell==1 and follow==1):
            N.AdjustWindows(zoom=zoom,Cell=1,Xcell=C[0].Xcell,Ycell=C[0].Ycell)
        #else:
            #N.AdjustWindows(zoom=zoom,Xcell=37.4, Ycell=55.42, Cell=1)
            #N.AdjustWindows(zoom=zoom,Cell=1,Xcell=0,Ycell=0)
        if save==1 and n>=nbeg:
            N.save("../Res/simulation"+str(SimNum)+"/S"+str(S)+"/Image/network_"+numero+".png")
        if show==1 and n>=nbeg:
            N.show()
        #-----------------------a modifier--------------------------
        #---------------------------------------------------------------
        N.close()
        del(N)
        n+=1    

    #os.system("cd saving && convert network_*.png anime.gif && cp anime.gif /home/h/H.Le/image")
    # print("cp ffmpeg /home/hugo/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/Image/ && cd /home/hugo/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/Image && ./ffmpeg -r "+str(n//tVid)+" -i network_%04d.png -q:a 0 -q:v 0 vid.avi")
    # os.system("cp ffmpeg /home/hugo/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/Image/ && cd /home/hugo/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/Image && ./ffmpeg -r "+str(n//tVid)+" -i network_%04d.png -q:a 0 -q:v 0 vid.avi")

    #print("cp ffmpeg /home/hugo/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/Image/ && cd /home/hugo/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/Image && ./ffmpeg -r "+str(n//tVid)+" -i network_%04d.png -vcodec mpeg4 -b "+str(Quality)+"M vid.avi")
    #os.system("cp ffmpeg /home/hugo/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/Image/ && chmod +x /home/hugo/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/Image/ffmpeg && cd /home/hugo/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/Image && ./ffmpeg -r "+str(n//tVid)+" -i network_%04d.png -vcodec mpeg4 -b "+str(Quality)+"M vid.mp4")
    os.system("ffmpeg -r "+str(n//tVid)+" -i ../Res/simulation"+str(SimNum)+"/S"+str(S)+"/Image/network_%04d.png -vcodec mpeg4 -b "+str(Quality)+"M ../Res/simulation"+str(SimNum)+"/S"+str(S)+"/Image/vid.mp4")
    #os.system("cp /home/hugo/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/Image/vid.mp4 /home/hugo/H.Le/simulation"+str(SimNum)+"/videos/vid"+str(S)+".mp4")
    #if(S!=1):
    #    os.system("rm -r /home/hugo/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/Image")

