from Folder import *
import pylab

class Cell:
    def __init__(self, folder):
        
        self.f=folder
        CellResFile=open(self.f.OpenCellFile(), 'r')    
        CellPos1=[i.split(" ") for i in CellResFile.read().split("\n")]
        del(CellPos1[-1])
        self.CellPos=[[float(j) for j in i] for i in CellPos1]
    def PlotCell(self,n):
        R=0.5*self.f.eps
        self.Xcell=self.CellPos[n-1][0]
        self.Ycell=self.CellPos[n-1][1]
        self.circle=pylab.Circle((self.Xcell,self.Ycell),radius=R,facecolor='r')
        pylab.gca().add_patch(self.circle)
    def UnPlotCell(self):
        try:
            self.circle.remove()
        except:
            pass
        #pylab.close()
    def OpenVector(self,n,i):
        VectorResFile=open(self.f.OpenVectFile(n,i), 'r')
        Vect1=[i.split(" ") for i in VectorResFile.read().split("\n")]
        del(Vect1[-1])
        self.Vect=[[float(j) for j in i] for i in Vect1]
    def PlotVect(self,n):
        x=list()
        y=list()
        u=list()
        v=list()
        for i in self.Vect:
            x.append(i[0])
            y.append(i[1])
            u.append(i[3])
            v.append(i[4])
        pylab.quiver(x,y,u,v,scale=10000.,angles='xy',scale_units='xy',width=0.005,edgecolor='k')
        
