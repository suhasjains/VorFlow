import numpy as np
from scipy.spatial import *

class Mesh(N,L_x,L_y):
	
	def __init__(self):
		# Don't want globals in mesh.py otherwise prohibitive to module testing!
		self.N = N;
		self.L_x = L_x;
		self.L_y = L_y;
		
		# Initialise mesh dataspace
		self.site = np.random.rand(N,2);
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
		# Call scipy meshing
		this = Voronoi(points);	
		# ...periodicity stuff:
			# Connectivity across periodic bounds
			# Remove the points outside domain
	
	
		# Calculates all properties here in individual for loops...
	
			# site

			# area

			# n_neighbor
		
			# neighbor
		
			# length

			# face

			# face_center

			# gred_area

			# grad_area_t

			# is_boundary
	
	
	
	

	













