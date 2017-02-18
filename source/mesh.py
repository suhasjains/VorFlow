import numpy as np
from scipy.spatial import *

class Mesh(N,L_x,L_y,BCs):
	
	def __init__(self):
		# Don't want globals in mesh.py otherwise prohibitive to module testing!
		self.N = N;
		self.L_x = L_x;
		self.L_y = L_y;
		self.is_periodic = np.all(BCs == 0);
		
		# Initialise mesh dataspace
		self.site = np.random.rand(N,2);
                self.site[:,0] = L_x * self.site[:,0]; self.site[:,1] = L_y * self.site[:,1];
		self.area = np.zeros(N);
                self.n_neighbor = np.zeros(N);
		# These will be lists of arrays, i.e. neighbour[i][j]
		self.neighbor = range(N);
                self.length = range(N);
                self.face = range(N);
                self.face_center = range(N);
                self.grad_area = range(N);
                self.grad_area_t = range(N);
                self.is_boundary = range(N);
		

	def update_mesh(self,data,dt):
		# Forward Euler advance data by dt
		for i in range(self.N):
			self.site[i,0] = self.site[i,0] + dt * data.u_vel[i];
			self.site[i,1] = self.site[i,1] + dt * data.v_vel[i];
		
		# Regenerate mesh
		self.generate_mesh();


	def generate_mesh(self):
		# Extend periodic domain
		if self.is_periodic:
			# Tile the sites NOTE: This is very inefficient! Later try using just a few critical edge sites.
			tiled_site = np.zeros([9*N,2]);
			tiled_site[:,0] = np.concatenate((self.site[:,0],self.site[:,0]-L_x,self.site[:,0],self.site[:,0]+L_x,self.site[:,0]-L_x,self.site[:,0]+L_x,self.site[:,0]-L_x,self.site[:,0],self.site[:,0]+L_x),axis=1)
			tiled_site[:,1] = np.concatenate((self.site{:,0],self.site[:,1]-L_x,self.site[:,1]-L_x,self.site[:,1]-L_x,self.site[:,1],self.site[:,1],self.site[:,1]+L_x,self.site[:,1]+L_x,self.site[:,1]+L_x),axis=1)		

		# Call Voronoi Mesher
		if self.is_periodic:
			voronoi = Voronoi(tiled_site);	
		else:
			voronoi = Voronoi(self.site);

		# Set the connectivity (self.neighbor)
                for i in range(N):
                        self.n_neighbor[i] = np.sum(voronoi.ridge_points[0:N,:] == i);
            
                    for i in range(N):
                        here = np.where(voronoi.ridge_points[0:N,:] == i); # Finds the site indices of points i
                        self.neighbor[i] = np.zeros(self.n_neighbor[i]);
                        for j in range(self.n_neighbor[i]):
                                # Pick out the index which is across from i
                                self.neighbor[i][j] = voronoi.ridge_points[here[0][j], not here[1][j]];
                
                # Delete last 8/9 of voronoi

		# Calculate each of the properties sequentially (in individual for loops):
		
			
			# area	
		
			# length

			# face

			# face_center

			# gred_area

			# grad_area_t

			# is_boundary
	
	
	
	

	













