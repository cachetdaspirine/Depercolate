
#/usr/bin/ipython

import pylab
from Network import *
from Folder import *
import os
#nom du fichier d entree
# folder="./network/"
# param="parameters.log"
# file="Positions"
# extension=".res"

entier=1
Nserie=10
Nimage=50
numshow="none"
SimNum=1
tVid=10
os.system("rm -r cracking")
os.system("mkdir cracking")
for S in range(0,Nserie):
    if S==0:
        Fichier="/data/H.Le/simulation"+str(SimNum)+"/S/network/"
        Fichier2="/data/H.Le/simulation"+str(SimNum)+"/S/"
    else:
        Fichier="/data/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/network"+str(S)+"/"
        Fichier2="/data/H.Le/simulation"+str(SimNum)+"/S"+str(S)+"/"
    Fold=Folder(Folder=Fichier,ParameterFile="Parameters.log",ResFile="Positions",Ext=".res")
    Fold2=Folder(Folder=Fichier2, ResFile="Stress",Ext=".res")
    print(Fold2.OpenStressFile())
    
    rm= "cd ./cracking && rm hole" +str(S)
    mkdir = "cd ./cracking && mkdir hole"+str(S)
    os.system(rm)
    os.system(mkdir)
    NumFile=len(os.listdir(Fichier))
    if entier==1:
        Nimage=NumFile
    
    for num in range(0,Nimage):
        n=NumFile-num-2
        print(n)
        time=n*Fold.dt
        try:
            N=Network(filename=Fold.OpenDataFile(n),time=time,sizex=Fold.Lx,sizey=Fold.Ly,sizez=Fold.Lz,distance=Fold.eps)
        except:
            print("fin de la boucle pour "+str(n))
            break
        N.plot2D_fiber(cracking=Fold.PCrack)
        N.AdjustWindows(zoom=1.8,border=0)
        if num == numshow:
            N.show()
        #for i in range(0,Fold.pas):
        numero=str(n)#+i)
        if len(numero)==1 :
            numero="000"+numero
        if len(numero)==2 :
            numero="00"+numero
        if len(numero)==3 :
            numero="0"+numero
        
        N.save("./cracking/hole"+str(S)+"/"+numero+".png")
        N.close()
        del(N)
        time=time-Fold.dt*Fold.pas
        n+=-1
    print("cp ffmpeg /home/h/H.Le/image/cracking/hole"+str(S)+" && cd cracking/hole"+str(S)+" && ./ffmpeg -r "+str(int(Nimage//tVid))+" -i %04d.png -q:a 0 -q:v 0 vid.avi")
    os.system("cp ffmpeg /home/h/H.Le/image/cracking/hole"+str(S)+" && cd cracking/hole"+str(S)+" && ./ffmpeg -r "+str(int(Nimage//tVid))+" -i %04d.png -q:a 0 -q:v 0 vid.avi")

