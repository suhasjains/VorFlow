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
				self.voronoi = 0; # Stores the Voronoi Mesh Scipy object
				

		def update_mesh(self,data,dt):
				# Forward Euler advance data by dt
				for i in range(self.N):
						self.site[i,0] = self.site[i,0] + dt * data.u_vel[i];
						self.site[i,1] = self.site[i,1] + dt * data.v_vel[i];
				
				if is_periodic:
						for i in range(self.N):
								# Check if cell has been advected across the periodic boundary
							if self.site[i,0] < 0.:
									self.site[i,0] += self.L_x;
							if self.site[i,0] >= self.L_x:
									self.site[i,0] -= self.L_x;
							if self.site[i,1] < 0.:
									self.site[i,1] += self.L_y;
							if self.site[i,1] >= self.L_y:
									self.site[i,0] -= self.L_y;
	
				# Regenerate mesh
				self.generate_mesh();


		def generate_mesh(self):
				# Extend periodic domain
				if self.is_periodic:
						# Tile the sites NOTE: This is very inefficient! Later try using just a few critical edge sites. However, this will require rewriting the periodic connectivity algorithm. (AW)
						tiled_site = np.zeros([9*N,2]);
						tiled_site[:,0] = np.concatenate((self.site[:,0],self.site[:,0]-L_x,self.site[:,0],self.site[:,0]+L_x,self.site[:,0]-L_x,self.site[:,0]+L_x,self.site[:,0]-L_x,self.site[:,0],self.site[:,0]+L_x),axis=1)
						tiled_site[:,1] = np.concatenate((self.site{:,0],self.site[:,1]-L_y,self.site[:,1]-L_y,self.site[:,1]-L_y,self.site[:,1],self.site[:,1],self.site[:,1]+L_y,self.site[:,1]+L_y,self.site[:,1]+L_y),axis=1)		

				# Call Voronoi Mesher
				if self.is_periodic:
						voronoi = Voronoi(tiled_site);
				else:
						voronoi = Voronoi(self.site);

				# Set the connectivity
				for i in range(N):
						self.n_neighbor[i] = np.sum(voronoi.ridge_points[0:N,:] == i);

				# Neighbours
				for i in range(N):
						here = np.where(voronoi.ridge_points[0:N,:] == i); # Finds the site indices of point i
						self.neighbor[i] = np.zeros(self.n_neighbor[i]);
						for j in range(self.n_neighbor[i]):
								# Pick out the index which is across from i
								self.neighbor[i][j] = voronoi.ridge_points[here[0][j], not here[1][j]];
				

				# Calculate each of the properties sequentially (in individual for loops):
				# Only for first N sites!

				# Length
				for i in range(N):
						self.length[i] = np.zeros(self.n_neighbor[i],2);
						for j in range(self.n_neighbor[i]):
								self.length[i][j,:] = self.site[self.neighbor[i][j]] - self.site[i];

				# Face
				for i in range(N):
						ridge_indices = np.where(voronoi.ridge_points[0:N,:] == i)[0]; # Finds the ridge indices of point i
						self.face[i] = np.zeros(self.n_neighbor[i]);
						for j in range(self.n_neighbor[i]):
								vertex_indices = voronoi.ridge_vertices[ridge_indices[j]];
								f = voronoi.vertices[vertex_indices[0],:] - voronoi.vertices[vertex_indices[1],:];
								self.face[i][j] = np.sqrt(f[0]**2 + f[1]**2);


				# Area
				
				
				# FaceCentre
				
				
				# GradArea
				
				
				# GradAreaT


				# isBoundary
				for i in range(N):
						self.is_boundary(i) = False; # Extend this after Infinite domain works.


				# Fix neighbours and sew edge connectivity together
				if self.is_periodic:
						for j in range(self.n_neighbor[i]):
								escaped = False;
								where = np.zeros(2);
								site_neighbour = self.site[self.neighbor[i][j],:];
								if site_neighbour[0] < 0.:
										escaped = True;
										where[0] = -1;
								if site_neighbour[0] >= self.L_x:
										escaped = True;
										where[0] = +1;
								if site_neighbour[1] < 0.:
										escaped = True;
										where[1] = -1;
								if site_neighbour[1] >= self.L_y:
										escaped = True;
										where[1] = +1;
						
						# Find missing soul mate <3
						if escaped:
								moveback = where[0]*1 + where[1]*3 + 5;
								if moveback > 4: moveback -= 1;
								self.neighbor[i][j] -= moveback;




		













