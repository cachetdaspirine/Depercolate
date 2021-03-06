import mayavi.mlab
from pylab import *
from numpy import *
import scipy.ndimage
from scipy.spatial import Delaunay
import numpy.linalg
import numpy.ma
from DeformationData import *

 
import sys
sys.path.append("/home/ronceray/Dropbox/WORK/MYTOOLKITS/") 
import mytoolkit

######## PYLAB OUTPUT ############################

fig_size = [12.,7.5]
params = {'axes.labelsize': 8,
          'text.fontsize':   6,
          'legend.fontsize': 7,
          'xtick.labelsize': 6,
          'ytick.labelsize': 6,
          'text.usetex': False,
          'figure.figsize': fig_size,
          }
rcParams.update(params)

pointsize = 3

fig = figure(1)
clf()
 
fig.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.99,
                wspace=0.05, hspace=0.05)

ax = subplot(111)

# Convert from lattice coordinates to cartesian (works for FCC and TRI2D)
Lattice_to_cartesian =  (0.5**0.5)*array( [[0.,1.,1.],\
                       [1.,0.,1.],\
                       [1.,1.,0.] ]  )
# Project so that (i,j) -> (x,y)
Cartesian_to_xy =  array( [[0.          , 2./(6.**0.5) , 1./(3.**0.5) ],\
                           [1./(2.**0.5),-1./(6.**0.5) , 1./(3.**0.5) ],\
                           [1./(2.**0.5), 1./(6.**0.5) ,-1./(3.**0.5) ] ]  )

Passage = dot(Lattice_to_cartesian,Cartesian_to_xy)

def fcc_to_cartesian(vec):
    return dot(vec,Passage)



def to_polar_coordinates(x,y):
    R = (x**2+y**2)**0.5
    theta = arccos(x/R) if x >= 0. else -arccos(x/R)     
    return (R,theta)


def radial_decomposition(coords,mesh = 1., filter_list = None):
    cX,cY = sum(X[0] for X in coords)/len(coords),sum(X[1] for X in coords)/len(coords)
    return [ to_polar_coordinates(X[0]-cX,X[1]-cY) for X in coords ]

# Colors :

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

#edges_cm = LinearSegmentedColormap('Edges', cdict)
edges_cm = cm.RdYlGn
edges_cm = cm.Spectral
edges_cm = cm.RdYlBu
edges_cm = cm.coolwarm
edges_cm = cm.copper


def edge_cmap_default(x,alpha_fac):
    return edges_cm(x)


def edge_cmap(x,alpha_fac):
    r,g,b,a = cm.RdYlBu(x)
    f = 2*(x*(1-x))**2
    r *=  (1-2*f)
    g *=  (1-f)
    b *=  (1-2*f)
    return array((r,g,b,a))

"""
def edge_cmap(x,alpha_fac):
    r = max(0.8*(x-0.5),0.)+0.2
    g = 0.5*(1-4*abs(x-0.5)**2)+0.5*max(0.5-x,0.)
    alpha = 1. - alpha_fac*(1-4*(x-0.5)**2)
    b = max(1.*(0.5-x),0.)
    return (r,g,b,alpha)
"""

def edge_cmap_mayavi(alpha_fac,brightness = 1):
    return array([[floor(255*(1- brightness*(1-a))) for a in edge_cmap(u/255.,alpha_fac)] for u in range(256)])



class Network:
    """ Reads XXXnetwork.dat files from C++ program.

    This class reads the network file, is able to display the network
    and to deform it using results stored int XXXspatial.dat.
    """ 

    def __init__(self,ID,direc = "../RESULTS/",center=True):
        # Initializer : read network file.
 
        lines = [L.split() for L in open(direc+ID+"network.dat").readlines()]
        header = lines.pop(0)

        self.direc = direc
        self.ID = ID
        self.lattice = header[0]
        self.w,self.h,self.d = int(header[1]),int(header[2]),int(header[3])
        self.Nsites,self.Nedges = int(header[4]),int(header[5])
        self.edge_cmap = edge_cmap
        
        # Sites are stored as lattice coordinates (i,j,k,index) :
        self.sites = [ array([ int(i) for i in pos]) for pos in lines[0:self.Nsites]]
        # Edges are stored as tuples of sites :
        self.edges = [ tuple([ int(i) for i in indices[:2]]) for indices in lines[self.Nsites+1:]]
        # Neighbourhood set: sites->edges map
        
        self.edge_reference_vectors = [ array([ double(i) for i in indices[2:5]]) for indices in lines[self.Nsites+1:]]
        self.spring_constants = [  (1. if len(indices) < 6 else float(indices[5])) for indices in lines[self.Nsites+1:]] # condition for compatibility with pre-3.2 data

        self.coordinates = [ array([0.,0.,0.]) for s in self.sites ]
        self.cursor = 0
        self.maxgamma = 0.
        self.D = DeformationData(ID)
        self.read_deformation_data()

        self.plist,self.kappalist = mytoolkit.load_columns(direc+ID+"results.dat",(3,4))
        self.p = self.plist[0]

        self.center()

    def edge_cmap_mayavi(self,alpha_fac,brightness = 1):
        return array([[floor(255*(1- brightness*(1-a))) for a in self.edge_cmap(u/255.,alpha_fac)] for u in range(256)])

    def read_deformation_data(self):
        """ Read ***spatial.dat to get nonaffine deformation. Beware,
        that may be slow and memory-consuming. """

        lines = [L.split() for L in open(self.direc+self.ID+"spatial.dat").readlines()]
        data = [ [eval(tup) for tup in L[2:] ] for L in lines]
        self.gamma_list = [ array(eval(L[1])) for L in lines ]
        self.maxgamma = max( amax(abs(g)) for g in self.gamma_list) + 1.e-10

        self.all_coordinates = [ [ array(d[:3]) for d in dat   ] for dat in data]
        self.all_energies = [ array([ d[3] for d in dat   ]) for dat in data]
        self.all_nonaffinities = [ array([ d[4] for d in dat   ]) for dat in data]
        self.all_nonaffine_vectors = [ [ array(d[5:8]) for d in dat   ] for dat in data]
        self.connectivities = [ array([ d[8] for d in dat   ]) for dat in data][0]
        self.all_bending = [ array([ d[9] for d in dat   ]) for dat in data]

        self.cursor = 0

        self.stress_scale = 1. 
        self.maxE = max(max([max(E) for E in self.all_energies]),1.e-250)
        self.maxNA = max(max([max(NA) for NA in self.all_nonaffinities]),1.e-250)
        self.read_motors_data()
        self.read_forces()
        self.push(0)

    def rotate(self,rotation_matrix,center=array([0,0,0])):
        self.all_coordinates = [ [ dot(vec-center,rotation_matrix)+center for vec in coords   ] for coords in self.all_coordinates]
        self.all_forces = [ [ dot(vec,rotation_matrix) for vec in coords   ] for coords in self.all_forces]
        self.push(0)

    def center(self,cen=None):
        for i,L in enumerate(self.all_coordinates):
            if cen is None:
                c = L[0]*0.
                for v in L:
                    c+=v
                c /= len(L)
            else:
                c = cen
            self.all_coordinates[i] = [v-c for v in L]
        self.push(0)

    def rotate2d(self,angle,center=array([0,0,0]),degrees=True):
        if degrees:
            theta = angle * pi / 180.
        else:
            theta = angle
        M = array([[ cos(theta), sin(theta), 0. ],
                   [ -sin(theta),cos(theta), 0. ],
                   [ 0.,0., 1 ]])
        self.rotate(M,center)
        self.push(0) 

    def read_motors_data(self):
        """ Read active unit info."""
        lines = [L.split() for L in open(self.direc+self.ID+"dipoles.dat").readlines()]
        self.all_dipoles_indices = [ [  u  for coords in L for u in eval(coords)[:2]  ] for L in lines]
        self.all_dipoles_energies = [ [ eval(coords)[2] for coords in L   ] for L in lines]
        self.dipoles_energies = (self.all_dipoles_energies[0] if len(self.all_dipoles_indices)>0 else [])
        self.dipoles_indices = (self.all_dipoles_indices[0] if len(self.all_dipoles_indices)>0 else [])

        try:
            self.all_processive_motors = [eval( (L[:-3]+']' if len(L)>3 else L)) for L in open(self.direc+self.ID+"processive_motors.dat").readlines()]
        except:
            self.all_processive_motors = [ [] for i in self.all_dipoles_indices]
        self.processive_motors = self.all_processive_motors[0]

    def read_forces(self):
        """ Read ***forces.dat to get elastic info."""
        lines = [L.split() for L in open(self.direc+self.ID+"forces.dat").readlines()]
        self.all_forces = [ [ array(eval(coords)) for coords in L  ] for L in lines]
        self.forces = self.all_forces[0]
        # Edge stresses :
        lines = [L.split() for L in open(self.direc+self.ID+"edge_stress.dat").readlines()]
        self.all_edge_stresses = [ [ array(eval(coords)[1]) for coords in L  ] for L in lines]
        self.edge_stresses = self.all_edge_stresses[0]
        self.all_edge_forces = [ [ array(eval(coords)[0]) for coords in L  ] for L in lines]
        self.edge_forces = self.all_edge_stresses[0]


 
    def push(self, n = 1, reset = False):
        self.cursor += n
        if reset :
            self.cursor = n

        self.coordinates = self.all_coordinates[self.cursor%len(self.all_coordinates)]
        self.energies = self.all_energies[self.cursor%len(self.all_coordinates)]
        self.nonaffinities = self.all_nonaffinities[self.cursor%len(self.all_coordinates)]
        self.dipoles_energies = (self.all_dipoles_energies[self.cursor%len(self.all_coordinates)] if len(self.all_dipoles_indices)>0 else [])
        self.dipoles_indices = (self.all_dipoles_indices[self.cursor%len(self.all_coordinates)] if len(self.all_dipoles_indices)>0 else [])
        self.processive_motors = self.all_processive_motors[self.cursor%len(self.all_coordinates)]
        self.nonaffine_vectors = self.all_nonaffine_vectors[self.cursor%len(self.all_coordinates)]
        self.forces = self.all_forces[self.cursor%len(self.all_coordinates)]
        self.edge_stresses = self.all_edge_stresses[self.cursor%len(self.all_coordinates)]
        self.edge_forces = self.all_edge_forces[self.cursor%len(self.all_coordinates)]
        self.gamma = self.gamma_list[self.cursor%len(self.all_coordinates)]
        self.bending = self.all_bending[self.cursor%len(self.all_coordinates)]
        

    def quiver_NA(self, slice_index = 0, show_copies = True, clear = True, amplification = 0., normalized = False, norm_power = 0.5,remove_average=True):
        site_coords = zip(*[ R[:3] for R in self.coordinates])
        X,Y = array(site_coords[0]),array(site_coords[1])
        arrows_coords = zip(*(self.nonaffine_vectors))
        
        U = array(arrows_coords[0])
        V = array(arrows_coords[1])
        if remove_average:
            U -= sum(U)/len(U)
            V -= sum(V)/len(V)

 
    def quiver_NA(self, slice_index = 0, show_copies = True, clear = True, amplification = 0., normalized = False, norm_power = 0.5,remove_average=True,headwidth=0.1):
        site_coords = zip(*[ R[:3] for R in self.coordinates])
        X,Y = array(site_coords[0]),array(site_coords[1])
        arrows_coords = zip(*(self.nonaffine_vectors))
        
        U = array(arrows_coords[0])
        V = array(arrows_coords[1])
        if remove_average:
            U -= sum(U)/len(U)
            V -= sum(V)/len(V)
        
        X,Y = X-U,Y-V
        ampli = amplification
        if normalized:
            U,V = U/(U**2 + V**2)**norm_power,  V/(U**2 + V**2)**norm_power
            ampli = 0.9 / max(U**2+V**2)**0.5

        if amplification == 0. and not normalized:
            ampli = self.stress_scale

        quiver(X,Y,ampli*U,ampli*V,scale = 1.0,units = 'xy',width = headwidth,headwidth = 2.0,headlength = 2.0,color = 'g')
        frame1 = gca()
        frame1.axes.set_aspect('equal')
        frame1.axes.get_xaxis().set_visible(False)
        frame1.axes.get_yaxis().set_visible(False)
        show()
 
    def plot_forces(self, clear = True, norm_passive = None, norm_active = None, amplification = 0.,amplification_active = 1.,color_active='tomato',color_passive="k",width=0.3):
        if amplification == 0.:
            amplification = 2.

        if norm_passive is None :
            norm_passive = max( [sum([ f**2 for f in F]) for i,F in enumerate(self.forces) if i not in self.dipoles_indices]) ** 0.5
        print "Force norm (passive) : ", norm_passive
        site_coords = zip(*[ R[:3] for i,R in enumerate(self.coordinates) if i not in self.dipoles_indices])
        X,Y = array(site_coords[0]),array(site_coords[1])
        
        arrows_coords = zip(*([F for i,F in enumerate(self.forces) if i not in self.dipoles_indices]))        
        U = array(arrows_coords[0])/norm_passive
        V = array(arrows_coords[1])/norm_passive
        #quiver(X,Y,amplification*U,amplification*V,scale = 1.0,units = 'xy',width = width,color =color_passive,minlength=0.)
        quiver(X,Y,amplification*U,amplification*V,scale = 1.0,units = 'xy',width = width,edgecolor ='k',facecolor=color_passive,minlength=0.,linewidth=0)

        if len(self.dipoles_indices)>0 or len(self.processive_motors)>0 :
            
            point_forces_coordinates = [ self.coordinates[i] for i in self.dipoles_indices ] + \
                           [ self.coordinates[M['center']]+array(PF[0]) for M in self.processive_motors for PF in M['crossing forces'] ]
            point_forces = [ self.forces[i] for i in self.dipoles_indices ] + \
                           [ array(PF[1]) for M in self.processive_motors for PF in M['crossing forces'] ]

            if norm_active is None :
                norm_active = max( norm(f) for f in point_forces )

            print "Force norm (active) : ", norm_active
            Xa,Ya = array([ R[0] for R in point_forces_coordinates]),array([ R[1] for R in point_forces_coordinates])
            Ua,Va = array([ F[0] for F in point_forces]) * (amplification_active/norm_active),array([ F[1] for F in point_forces]) * (amplification_active/norm_active)
            quiver(Xa,Ya,Ua,Va,scale = 1.0,units = 'xy',width = 0.6,color = color_active,minlength=0.,edgecolor ='k',lw=1)
 
        xlim(min(X)-3.5,max(X)+3.5)
        ylim(min(Y)-3.5,max(Y)+3.5)

        frame1 = gca()
        frame1.axes.set_aspect('equal')
        frame1.axes.get_xaxis().set_visible(False)
        frame1.axes.get_yaxis().set_visible(False)
        show()
        return norm_passive,norm_active
        
    def plot_stress(self, clear = True, amplification = 0., norm_power = 0.5,on_edge=False):
        site_coords = zip(*[ R[:3] for R in self.coordinates])
        X,Y = array(site_coords[0]),array(site_coords[1])
        U0,U1,V0,V1,ev = [],[],[],[],[]
        
        for i,S in enumerate(self.edge_stresses):
            (s1,s2) = self.edges[i]
            x,y = 0.5*(X[s1]+X[s2]),0.5*(Y[s1]+Y[s2])
            eigvals,eigvecs = linalg.eigh(S[:2,:2].T + S[:2,:2])
            if sum(eigvals) != 0.:
                print S
                for j in range(2):
                    U0.append( x - eigvecs[0,j]/2.)
                    U1.append(  eigvecs[0,j])
                    V0.append( y - eigvecs[1,j]/2.)
                    V1.append( eigvecs[1,j])
                    ev.append(eigvals[j])
        ev = array(ev)
        norm = max(abs(ev))
        values = sign(ev) * (abs(ev)/norm)**norm_power
        colors = [ ( min(1., 1.+ V ), min(1.- V, 1.+ V) , min(1., 1.- V)) for V in values ] 
        quiver(U0,V0,U1,V1,scale = 1.0,units = 'xy',width = 0.2,headwidth = 1.,headlength = .0 , color = colors)
        frame1 = gca()
        frame1.axes.set_aspect('equal')
        frame1.axes.get_xaxis().set_visible(False)
        frame1.axes.get_yaxis().set_visible(False)
        show()




 
        

    def show(self, show_copies = True, clear = True, plot_sites = False, width = 0.2, color_elongation = False, alpha_fac = 0., color_max = 0.,gray_exponent=1.,buckled_links_only = False,buckled_links_threshold = 0.5,plot_buckled_midedges=False,xmin=0,xmax=1):
        """ Display the system, or if 3D a z-fixed slice of the system (3D todo) """
        if clear:
            cla()        
        # X,Y = edge origins ; U,V : edge vectors ; C : colors
        X,Y,U,V,C = [],[],[],[],[]
        max_spring_constant = ( color_max if color_max != 0. else max(self.spring_constants) )

        if buckled_links_only:
            # Plot buckled mid edges :
            self.find_buckling_zone_bending(threshold=0.005)
            S = set(self.buckled_midedges) # search in a set is faster

            for edge_index,(i1,i2) in enumerate(self.edges):
                if ((i1 in S) or (i2 in S)):
                    vec1,vec2 = self.coordinates[i1][:2],self.coordinates[i2][:2]
                    if numpy.linalg.norm(vec2-vec1) < 2.5 + 5.*self.maxgamma :
                        X.append(vec1[0])
                        Y.append(vec1[1])
                        U.append(vec2[0] - vec1[0])
                        V.append(vec2[1] - vec1[1])
            (r,g,b,a) = self.edge_cmap(0,0)

            quiver(X,Y,U,V, scale = 1.0,units = 'xy',width =  width,  headwidth = 1.0,headlength = 0.0, headaxislength=0. , color = (r,g,b,1-alpha_fac),edgecolor=(r,g,b,1-alpha_fac) )
 
        else:

            for edge_index,(i1,i2) in enumerate(self.edges):
                s1,s2 = self.sites[i1],self.sites[i2]
                vec1,vec2 = self.coordinates[i1][:2],self.coordinates[i2][:2]
                if numpy.linalg.norm(vec2-vec1) > 5.5 + 5.*self.maxgamma :
                    # winding edge                   
                    pass
                else:
                    # Bulk case :

                    if color_elongation:
                        affine_d = dot(self.edge_reference_vectors[edge_index],self.gamma)
                        L_minus_one = dot( (affine_d + self.nonaffine_vectors[i2] - self.nonaffine_vectors[i1]).T, 2*self.edge_reference_vectors[edge_index] +(affine_d + self.nonaffine_vectors[i2] - self.nonaffine_vectors[i1]) )/(1. + numpy.linalg.norm(vec2-vec1))
                        N = self.stress_scale * L_minus_one
                        x = (1. + N/sqrt(1+N**2) )/2.

                    else :
                        N = self.stress_scale * self.edge_stresses[edge_index].trace()
                        x = (1. + N/sqrt(1+N**2) )/2.
                        r = max(0.8*(x-0.5),0.)+0.2
                        g = 0.5*(1-4*abs(x-0.5)**2)+0.5*max(0.5-x,0.)
                        alpha = 1. - alpha_fac*(1-4*(x-0.5)**2)
                        b = max(1.*(0.5-x),0.)
                        #C.append( (r,g,b,alpha) )

                    if xmin<=x<=xmax: 
                        C.append( self.edge_cmap(x,alpha_fac) )
                        X.append(vec1[0])
                        Y.append(vec1[1])
                        U.append(vec2[0] - vec1[0])
                        V.append(vec2[1] - vec1[1])
                    #C.append( (0.05,0.4,0.2,1.) )
            quiver(X,Y,U,V, scale = 1.0,units = 'xy',width =  width,  headwidth = 1.0,headlength = 0.0, headaxislength=0. , color = C,edgecolor=(0.05,0.4,0.2,1.) )

 
        x_points,y_points = array([v[0] for i,v in enumerate(self.coordinates) ]),array([v[1] for i,v in enumerate(self.coordinates)  ])
        # Plot sites :
        if plot_sites:
            quiver(x_points,y_points,0.*x_points,0.*y_points, scale = 1.0,units = 'xy',width =  width, headwidth = 1.0,headlength = 0.0, headaxislength=0. , color = 'w')



        # Center the image :
        x0 = sum(x_points)/len(x_points)
        y0 = sum(y_points)/len(y_points)
        xmin = x0 - (0.9 + self.maxgamma) * self.w
        xmax = x0 + (0.9 + self.maxgamma) * self.w
        ymin = y0 - 0.7 * self.h
        ymax = y0 + 0.8 * self.h

        xlim(-1.2*(x0-min(x_points))+x0,1.2*(max(x_points)-x0)+x0)
        ylim(-1.2*(y0-min(y_points))+y0,1.2*(max(y_points)-y0)+y0)

        frame1 = gca()
        frame1.axes.set_aspect('equal')
        frame1.axes.get_xaxis().set_visible(False)
        frame1.axes.get_yaxis().set_visible(False)
        show()



               
 
    def plot_motors(self,color="tomato",size=None,zorder=10):
        for M in self.processive_motors:
            if size is None:
                size = M['radius']/2
            #circle = Circle(self.coordinates[M['center']][:2], M['radius'], facecolor=(1,0,0,0.3), edgecolor='0.5')
            #gca().add_patch(circle)
            circle = Circle(self.coordinates[M['center']][:2], size, facecolor=color, edgecolor='none',zorder=zorder)
            gca().add_patch(circle)





    def plot_sites_3D(self,size=0.8):
        #fast display of a 3D point list as spheres.

        #Option values allows to display different colors for different
        #subsets, by default all colors are different.
        f = mayavi.mlab.gcf()
        f.scene.background = (1.,1.,1.)
        f.scene.parallel_projection = True

        values = [ e for ind,e in enumerate(self.energies)   ]   
        points = [p[:3] for ind,p in enumerate(self.coordinates) ]    # a copy : do not modify the original point_list !!
        x,y,z = tuple(zip(*points))
        mayavi.mlab.points3d(x, y, z ,values, scale_factor = size, colormap = 'summer', scale_mode = 'none', resolution = 16)
        mayavi.mlab.show()

    def plot_edges_3D(self,max_sites = 10000.,cmap = "RdBu",mode='cylinder',resolution=8,transparent=False,lw = 10,alpha_fac=1.,brightness=1.,buckled_cmap = False, zmin=-1e100, zmax=1e100 , slice_direction=array([0,0,1.]), minvalue=0,maxvalue=1,clear=True):
        """ Plot the edges between points whose position is given in
        input."""
        if clear:
            mayavi.mlab.clf()
        f = mayavi.mlab.gcf()
        f.scene.background = (1.,1.,1.)
        f.scene.parallel_projection = True

        
        points_positions = [ R[:3] for R in self.coordinates]

        # X,Y = edge origins ; U,V : edge vectors ; C : colors
        X,Y,Z,U,V,W,C = [],[],[],[],[],[],[]

        if buckled_cmap:
            # Plot buckled mid edges :
            self.find_buckling_zone_bending(threshold=0.005)
            S = set(self.buckled_midedges) # search in a set is faster
        
        for edge_index,(i1,i2) in enumerate(self.edges):
            s1,s2 = self.sites[i1],self.sites[i2]
            if max(max(s1),max(s2))<max_sites:
                vec1,vec2 = points_positions[i1],points_positions[i2]
                if zmin < slice_direction.dot(vec1) < zmax and zmin < slice_direction.dot(vec2) < zmax :
                    if numpy.linalg.norm(vec2-vec1) > 3. + 5.*self.maxgamma :
                        # winding edge 
                        pass
                    else:
                        # Bulk case :
                        """
                        affine_d = dot(self.edge_reference_vectors[edge_index],self.gamma)
                        L_minus_one = dot( (affine_d + self.nonaffine_vectors[i2] - self.nonaffine_vectors[i1]).T, 2*self.edge_reference_vectors[edge_index] +(affine_d + self.nonaffine_vectors[i2] - self.nonaffine_vectors[i1]) )/(1. + numpy.linalg.norm(vec2-vec1))
                        N = self.stress_scale * L_minus_one
                        x = (1. + N/sqrt(1+N**2) )/2.
                        r = 0.85*x
                        g = 0.85*(1-4*(x-0.5)**2)
                        b = 0.85*(1-x)
                        """
                        
                        N = self.stress_scale * self.edge_stresses[edge_index].trace()
                        x =   0.5*(1+ N/sqrt(1+N**2) )

                        if buckled_cmap or minvalue <= x <= maxvalue:

                            X.append(vec1[0])
                            Y.append(vec1[1])
                            Z.append(vec1[2])
                            U.append(vec2[0] - vec1[0])
                            V.append(vec2[1] - vec1[1])
                            W.append(vec2[2] - vec1[2]) 
                        if buckled_cmap:
                            if ((i1 in S) or (i2 in S)):
                                C.append(0)
                            else:
                                C.append(0.5)
                        elif   minvalue <= x <= maxvalue:
                            C.append( x ) 


        C = array(C)
        vecs = mayavi.mlab.quiver3d(X,Y,Z,U,V,W, line_width = lw, scalars= C, mode=mode, scale_factor = 1.,scale_mode = 'vector', resolution = resolution,vmin=0,vmax=1,colormap=cmap,transparent=transparent)
        vecs.glyph.color_mode = 'color_by_scalar'
        vecs.module_manager.scalar_lut_manager.lut.table = self.edge_cmap_mayavi(alpha_fac,brightness)
        mayavi.mlab.show()
        mayavi.mlab.draw()
        return vecs

    def plot_forces_3D(self, clear = True, norm_passive = None, norm_active = None, amplification = 0.,amplification_active = 1.):
        if amplification == 0.:
            amplification = 2.

        if norm_passive is None :
            norm_passive = max( [sum([ f**2 for f in F]) for i,F in enumerate(self.forces) if i not in self.dipoles_indices]) ** 0.5

        # Truncate mini-forces to save memory:
        passive_indices = [i for i,F in enumerate(self.forces) if i not in self.dipoles_indices and norm(F)/norm_passive > 1e-4 ]

        site_coords = zip(*[ self.coordinates[i] for i in passive_indices])
        X,Y,Z = array(site_coords[0]),array(site_coords[1]),array(site_coords[2])

        arrows_coords = zip(*[ self.forces[i] for i in passive_indices])
        U = array(arrows_coords[0])/norm_passive
        V = array(arrows_coords[1])/norm_passive
        W = array(arrows_coords[2])/norm_passive
        print "Force norm (passive) : ", norm_passive,"   --   points represented:",len(U),"/",len(self.coordinates)
        

        if len(self.dipoles_indices)>0 or len(self.processive_motors)>0 :
            
            point_forces_coordinates = [ self.coordinates[i] for i in self.dipoles_indices ] + \
                           [ self.coordinates[M['center']]+array(PF[0]) for M in self.processive_motors for PF in M['crossing forces'] ]
            point_forces = [ self.forces[i] for i in self.dipoles_indices ] + \
                           [ array(PF[1]) for M in self.processive_motors for PF in M['crossing forces'] ]

            if norm_active is None :
                norm_active = max( norm(f) for f in point_forces )

            print "Force norm (active) : ", norm_active
            Xa,Ya,Za = array([ R[0] for R in point_forces_coordinates]),array([ R[1] for R in point_forces_coordinates]),array([ R[2] for R in point_forces_coordinates])
            Ua,Va,Wa = array([ F[0] for F in point_forces]) * (amplification_active/norm_active),array([ F[1] for F in point_forces]) * (amplification_active/norm_active),array([ F[2] for F in point_forces]) * (amplification_active/norm_active)
            vecs = mayavi.mlab.quiver3d(Xa,Ya,Za,Ua,Va,Wa, line_width = 5., mode='cone', scale_factor = 2.,scale_mode = 'vector', resolution = 16 ,color = (1.0, 0.38823529411764707, 0.2784313725490196))


        f = mayavi.mlab.gcf()
        f.scene.background = (1.,1.,1.)
        f.scene.parallel_projection = True 
        vecs = mayavi.mlab.quiver3d(X,Y,Z,U,V,W, line_width = 5., mode='cone', scale_factor = 2.,scale_mode = 'vector', resolution = 8 ,color = (0.,0.,0.))
        vecs.glyph.color_mode = 'color_by_scalar'
        mayavi.mlab.show()


    def plot_motors_3D(self,radius=None,color=(1.0, 0.38823529411764707, 0.2784313725490196),opacity=1.):
        X,Y,Z,R= [],[],[],[]
        for M in self.processive_motors:
            #circle = Circle(self.coordinates[M['center']][:2], M['radius']/2, facecolor=(1,0,0,1.), edgecolor='0.5')
            X.append(self.coordinates[M['center']][0])
            Y.append(self.coordinates[M['center']][1])
            Z.append(self.coordinates[M['center']][2])
            if radius is None:
                R.append(M['radius']/2)
            else:
                R.append(radius)
    
        mayavi.mlab.points3d(X,Y,Z,R, scale_factor = 1., color = color, scale_mode = 'scalar', resolution = 32,opacity=opacity)

        mayavi.mlab.show()



    def find_buckling_zone(self,threshold = 0.5,compressive_only=True):
        self.buckled_edges = []

        for edge_index,(i1,i2) in enumerate(self.edges):
            s = self.edge_stresses[edge_index].trace()
            if (s < -threshold) or ((not compressive_only) and abs(s)>threshold) :
                self.buckled_edges.append(edge_index)
 
    def find_buckling_zone_bending(self,threshold = 0.005):
        self.buckled_midedges = []
        for i,s in enumerate(self.bending):
            if (s > threshold):
                self.buckled_midedges.append(i)


    def buckling_zone_analysis(self,threshold):
        self.BZ_msr_buckled_list = []
        self.BZ_msr_all_list = []
        self.BZ_rmax_all_list = []
        self.BZ_rmax_buckled_list = []
        self.cursor = 0
        self.push(0)
        for n in self.all_coordinates:
            self.find_buckling_zone(threshold)
            self.BZ_msr_buckled_list.append(self.msr_buckled)
            self.BZ_msr_all_list.append(self.msr_all)
            self.BZ_rmax_buckled_list.append(self.rmax_buckled)
            self.BZ_rmax_all_list.append(self.rmax_all)
            self.push()
            


