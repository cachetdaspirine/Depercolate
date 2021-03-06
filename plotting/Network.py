#This module contain all the propreties of the object Network.

from Folder import *
#import mayavi.mlab
import numpy as np
import pylab
import os

from matplotlib.colors import LinearSegmentedColormap

cdict = {'blue':   ((0.0,  0.9,0.9),
                    (0.5,  0.4, 0.4),
                    (1.0,  0.1, 0.1)),

         'green': ((0.0,  0.5, 0.5),
                   (0.5 , 1, 1),
                   (1.0,  0.3, 0.3)),

         'alpha': ((0.0,  1, 1),
                   (0.5 , 0.8, 0.8),
                   (1.0,  1, 1)),

         'red':  ((0.0,  0.4, 0.4),
                   (0.5,  0.5, 0.5),
                   (1.0,  0.9,0.9)),
}
cm = LinearSegmentedColormap('my_colormap', cdict, 1024)
fig, ax = pylab.subplots()


# Here we define the class Network
class Network:
    
    #first the constructor, it need the name/path of the file, the type of simulation
    #we want to draw, and the time if we want to creat an animation.
    def __init__(self, folder,N=0):
        #we start by opening the file, the format is :
        #x y z neigh1x, neigh1y, neigh1z, neigh2x, ... , neighNy
        #neigh[i]j are three integers which give the position of the ith neighboor.
        self.f=folder
        #we start by opening the fiber and the node
        if(self.f.evolvcell==1):
            Node=open(folder.OpenNetworkFile(N),'r')
            Fiber=open(folder.OpenFiberFile(N),'r')
        else:
            Node=open(folder.OpenNetworkFile(N),'r')
            Fiber=open(folder.OpenFiberFile(N),'r')
            self.Length=np.loadtxt(folder.OpenLengthFile(N),dtype=float)
        
        
        self.DataNode=[i.split(" ") for i in Node.read().split("\n")]
        self.DataFiber=[i.split(" ") for i in Fiber.read().split("\n")]
        Node.close()
        Fiber.close()

        #six list needed to plot the file
        self.X=list()
        self.Y=list()
        self.Z=list()
        self.U=list()
        self.V=list()
        self.W=list()
        self.Rlength=list()


        
    def plot2D_fiber(self,cracking=1):
        
        #we will build line by line the "vector field" assiciated to the results file.
        
        if self.f.Lz > 1 :
            print("Lz != 0 impossible to plot 2d, you certainly want to plot a 3d network")

       
        for j in range(0,self.f.Ly):
            if(j%2==0):
                pair=1
            else:
                pair=0
            for k in range(0,self.f.Lz):
                for i in range(0,self.f.Lx):
                    self.X.append(float(self.DataNode[i+j*self.f.Ly+k*self.f.Lz][0]))
                    self.Y.append(float(self.DataNode[i+j*self.f.Ly+k*self.f.Lz][1]))
                    self.U.append(float(self.DataFiber[i+j*self.f.Ly+k*self.f.Lz][0]))
                    self.V.append(float(self.DataFiber[i+j*self.f.Ly+k*self.f.Lz][1]))
                    self.Rlength.append(self.Length[i+j*self.f.Ly,0])
                    
                    self.X.append(float(self.DataNode[i+j*self.f.Ly+k*self.f.Lz][0]))
                    self.Y.append(float(self.DataNode[i+j*self.f.Ly+k*self.f.Lz][1]))
                    self.U.append(float(self.DataFiber[i+j*self.f.Ly+k*self.f.Lz][2]))
                    self.V.append(float(self.DataFiber[i+j*self.f.Ly+k*self.f.Lz][3]))
                    self.Rlength.append(self.Length[i+j*self.f.Ly,1])

                    self.X.append(float(self.DataNode[i+j*self.f.Ly+k*self.f.Lz][0]))
                    self.Y.append(float(self.DataNode[i+j*self.f.Ly+k*self.f.Lz][1]))
                    self.U.append(float(self.DataFiber[i+j*self.f.Ly+k*self.f.Lz][4]))
                    self.V.append(float(self.DataFiber[i+j*self.f.Ly+k*self.f.Lz][5]))
                    self.Rlength.append(self.Length[i+j*self.f.Ly,2])

        Norm=(pylab.hypot(self.U,self.V)-self.Rlength)
        pylab.quiver(self.X,self.Y,self.U,self.V,Norm ,scale = 1.0,angles='xy',scale_units = 'xy',width = 0.005,minlength=0.,headlength=0.,headaxislength=0.,headwidth=0.,alpha=1,edgecolor='k',cmap=cm)
        #pylab.quiver(self.X,self.Y,self.U,self.V,Norm ,scale = 1.0,angles='xy',scale_units = 'xy',width = 0.005,minlength=0.,headlength=10,headaxislength=10,headwidth=10,alpha=1,edgecolor='k',cmap=cm)
        #pylab.quiver(self.X,self.Y,self.U,self.V ,scale = 1.0,angles='xy',scale_units = 'xy',width = 0.005,minlength=0.,headlength=0.,headaxislength=0.,headwidth=0.,alpha=1,edgecolor='k',cmap=cm)
        #pylab.clim(0.,0.5)

        #pylab.clim(-0.15,0.20)

        pylab.clim(-0.05,0.05)

        #set the color scale visible or not
        #pylab.gca().axis('off')
        #set the axis visible or not
        #pylab.gca().patch.set_visible(False)

        #pylab.colorbar()

        #pylab.quiver(self.Xlim,self.Ylim,self.Ulim,self.Vlim ,scale = 1.0,angles='xy',scale_units = 'xy',width = 0.0003,minlength=0.,lw=1,headlength=0.,headaxislength=0.,headwidth=0.,alpha=1,edgecolor='k')
       

    def show(self):
        pylab.show()
    def save(self,name):
        pylab.savefig(name)
    def close(self):
        pylab.close()
    def AdjustWindows(self,Cell=0,Xcell=25,Ycell=25,NX=10,NY=10,zoom=2,border=0):
        if Cell==1:
            Xcentre=Xcell
            Ycentre=Ycell
        else :
            Xcentre=0#self.f.Lx/2*self.f.eps
            Ycentre=0#self.f.Ly/2*self.f.eps*0.866
        Xlarge=self.f.Lx*self.f.eps/zoom
        Ylarge=self.f.Ly*0.866*self.f.eps/zoom

        pylab.xlim(Xcentre-Xlarge,Xcentre+Xlarge)
        pylab.ylim(Ycentre-Ylarge,Ycentre+Ylarge)

    def tabl(self):
    #Ici on cree un tableau a deux entres Lx*Ly,14 dans lequel il y a la position X,Y du noeud i*j, puis X1
    #Y1 du premier voisin ect...
        self.table = [[ 0 for _ in range(16)] for _ in range(len(self.line))]
        #self.table=[[0]*16]*len(self.line)
        for k,vertex in enumerate(self.line) :
            try:
                inf=[float(i) for i in vertex.split(" ")]
            except:
                pass
            for j,data in enumerate(inf):
                if j<len(inf)-2:
                    self.table[k][j]=data

    def StressProfile(self,direction=1,hole="false",mean="false",large=0,Nbande=0):
        #cree le profil de stress selon une coupe 'direction' (1=x/2=y) pour des bandes de largeur 'large'
        #par default il n y a pas de troue dans le reseaux. 
        #mean permet de faire une moyenne sur toute la largeur du reseaux
        #possibilite de creer N bandes plutot
        os.system("rm -r /home/h/H.Le/image/StressProfile")
        os.system("mkdir /home/h/H.Le/image/StressProfile")

        #il faut avoir cree un tableau avec tout ce qu'il faut

        #On cree ensuite N categorie selon la bande d'une certaine largeur

        if direction==1:
            LargeurTot=self.Lx
            ecart=0.5*self.distance
        if direction==2:
            LargeurTot=self.Ly
            ecart=0.577*self.distance
            
            
        if Nbande==0:
            Nbande=LargeurTot//large
        if large==0:
            large=LargeurTot*self.distance/(2*Nbande)
            #On parcours tous les noeuds, et on les places dans des listes selon leur position
        classi=[[]for _ in range(Nbande)]
        for num,info in enumerate(self.table):
            n=int(info[direction-1]//large)
            if n>=Nbande:
                n=Nbande-1
            classi[n].append(num)

        #maintenant pour chaque bande on calcul le profil de contrainte
        #i.e la contrainte applique a chaque noeuds selon 
        for num,info in enumerate(classi):
            sigm=[[0 for _ in range(4)]for _ in range(self.Lx)]
            ResFile=open("/home/h/H.Le/image/StressProfile/"+str(num),'w')
            for NumNoeud in info:
                x=self.table[NumNoeud][0]
                y=self.table[NumNoeud][1]
                sigmx=0
                sigmy=0
                ind=0
                #Ici on calcul le stress applique a un noeud
                for i in range(1,6):
                    xv=self.table[NumNoeud][2*i]
                    yv=self.table[NumNoeud][2*i+1]
                    if xv!=0 or yv!=0:
                        ind+=1
                    dist=((xv-x)**2+(yv-y)**2)**(1/2)
                    sigmx+=self.k*(dist-self.distance/2.)*(x-xv)**2
                    sigmy+=self.k*(dist-self.distance/2.)*(y-yv)**2
                if ind!=0:
                    sigmx=sigmx/ind
                    sigmy=sigmy/ind
                else :
                    sigmx=0
                    sigmy=0
                if self.table[NumNoeud][2-direction]<=1:
                    indice=0
                elif self.table[NumNoeud][2-direction]>self.Lx*self.distance/2 :
                    indice=99
                else :
                    indice=int(self.table[NumNoeud][2-direction]//ecart)                   
                print(str(LargeurTot)+" "+str(self.table[NumNoeud][2-direction])+" "+str(indice))
                sigm[indice][0]=self.table[NumNoeud][2-direction]
                sigm[indice][1]+=sigmx
                sigm[indice][2]+=sigmy
                sigm[indice][3]+=1
                
            for oui in sigm:
                if oui[3]!=0:
                    ResFile.write(str(oui[0]*2/self.distance)+" "+str(oui[1]/float(oui[3]))+" "+str(oui[2]/float(oui[3]))+"\n")

# class Cell:
    
#     def __init__(self,filename,radius,distance):
        
        
#         self.CellResFile=open(filename, 'r')
#         self.R=radius
#         self.distance=distance
    
#         contenue=self.CellResFile.read()
        
#         lines=contenue.split("\n")

#         del(lines[0])
#         del(lines[0])

#         self.line2=[i.split(" ") for i in lines]
#         for c,i in enumerate(self.line2):
#             # print(i)
#             # print(i[0])
#             # print(i[1])
#             # print(i[2])
#             for cc,k in enumerate(i):
#                 try:
#                     #float(line[c][cc])
#                     self.line2[c][cc]=float(line2[c][cc])
#                 except:
#                     pass
#                     #del(self.line[c][cc])
#         # for c,i in enumerate(line):
#         #     for cc,k in enumerate(i):
                
                
#         #donne=[[float(k) for k in i.split(" ")] for i in line]
    
#     def SetCoord(self,n):
#         self.Xcell=float(self.line2[n][0])#*self.distance*0.577
#         self.Ycell=float(self.line2[n][1])#*self.distance*0.5
#     def AddCell(self):
#         #print(self.Xcell)
#         #print(self.Ycell)
#         self.circle=pylab.Circle((self.Xcell,self.Ycell),radius=self.R,facecolor='r')
#         pylab.gca().add_patch(self.circle)
#     def Reset(self):
#         self.circle.remove()

        
                
    
            
                
               
            
    
        
        
        

    
