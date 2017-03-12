import numpy as np
import timeit
from scipy.spatial import *
from shapely.geometry import Polygon, Point
from sympy.geometry import *


class Mesh:
	
		def __init__(self,N,L_x,L_y,BCs,mesh_type):
				# Don't want globals in mesh.py otherwise prohibitive to module testing!
				self.N = N;
				self.L_x = L_x;
				self.L_y = L_y;
				self.is_periodic = np.all(BCs == 0);    # BCs = 0 is periodic
				if ~self.is_periodic:
					print "Non-periodic boundary is still in testing phase!\n\n\n"
				
				self.mesh_type = mesh_type;
				
				# Initialise mesh dataspace
				if mesh_type=="random": # Make random sites
						self.site = np.random.rand(N,2);
						self.site[:,0] = L_x * self.site[:,0];
						self.site[:,1] = L_y * self.site[:,1];

				elif mesh_type=="cartesian": # Make Cartesian grid (Square N, L_x=L_y)
						self.site = np.zeros([N,2]);
						sqrtN = int(np.sqrt(N));
						for i in range(sqrtN):
								for j in range(sqrtN):
										self.site[i*sqrtN+j,0] =  i;

						self.site[:,1] = np.tile(np.arange(sqrtN),sqrtN);
						self.site[:,0] = L_x*self.site[:,0]/sqrtN + L_x/(2.*sqrtN);
						self.site[:,1] = L_y*self.site[:,1]/sqrtN + L_y/(2.*sqrtN);

				elif mesh_type=="nonuniform": # Make non-uniform grid (Square N, L_x=L_y)
						self.site = np.zeros([N,2]);
						sqrtN = int(np.sqrt(N));
						for i in range(sqrtN):
								for j in range(sqrtN):
										self.site[i*sqrtN+j,0] = i;

						self.site[:,1] = np.tile(np.arange(sqrtN),sqrtN);
						#self.site[:,0] = L_x*(1. - np.cos((L_x*self.site[:,0]/sqrtN + L_x/(2.*sqrtN))*np.pi/(2.*L_x)));
						#self.site[:,1] = L_y*(1. - np.cos((L_y*self.site[:,1]/sqrtN + L_y/(2.*sqrtN))*np.pi/(2.*L_x)));
						self.site[:,0] = L_x*((self.site[:,0]/sqrtN + L_x/(2.*sqrtN))**2);
						self.site[:,1] = L_y*((self.site[:,1]/sqrtN + L_y/(2.*sqrtN))**2);
				else:
						raise ValueError("You can't spell, fool...")



				self.centroid = np.zeros([N,2]);
				self.area = np.zeros(N);
				self.N_neighbor = np.zeros(N,dtype=int);
				# These will be lists of arrays, i.e. neighbour[i][j]
				self.neighbor = range(N);
				self.length = range(N);
				self.face = range(N);
				self.face_center = range(N);
				self.grad_area = range(N);
				self.grad_area_t = range(N);
				self.is_boundary = range(N);
                                #This will be a list of dictionaries
                                self.boundary = range(N);

				#self.boundary_cutting_ridge = np.zeros([N,2]);
				#self.xBoundary_points = np.zeros([N,2]);
				#self.yBoundary_points = np.zeros([N,2]);
				self.voronoi = 0; # Make dataspace for storing Scipy Voronoi object for easy plotting
				
				self.generate_mesh();
				

		def update_mesh(self,data,dt):
				# Ensure data is valid
				if np.any(np.abs(data.u_vel) > 1.e16) or np.any(np.abs(data.v_vel) > 1.e16):
						raise ValueError("Do you really think this flow is incompressible? Really? It's pretty fast FYI.")	
				# Forward Euler advance data by dt
				for i in range(self.N):
						self.site[i,0] = self.site[i,0] + dt * data.u_vel[i];
						self.site[i,1] = self.site[i,1] + dt * data.v_vel[i];
				
				if self.is_periodic:
						for i in range(self.N):
								# Check if cell has been advected across the periodic boundary
								if self.site[i,0] < 0.:
										self.site[i,0] += self.L_x;
								elif self.site[i,0] >= self.L_x:
										self.site[i,0] -= self.L_x;
								if self.site[i,1] < 0.:
										self.site[i,1] += self.L_y;
								elif self.site[i,1] >= self.L_y:
										self.site[i,1] -= self.L_y;
	
				# Regenerate mesh
				self.generate_mesh();


		def generate_mesh(self):
				# Extend periodic domain
				
				profile = False;
				
				metatic = timeit.default_timer()
				
				if self.is_periodic:
						# Tile the sites -- NOTE: This is very inefficient! Later try using just a few critical edge sites.
						# However, this will require rewriting the periodic connectivity algorithm. (AW)
						tiled_site = np.zeros([9*self.N,2]);
						# Ordered tiling for easy indexing later
						tiled_site[:,0] = np.concatenate((self.site[:,0],
														  self.site[:,0]-self.L_x,
														  self.site[:,0],
														  self.site[:,0]+self.L_x,
														  self.site[:,0]-self.L_x,
														  self.site[:,0]+self.L_x,
														  self.site[:,0]-self.L_x,
														  self.site[:,0],
														  self.site[:,0]+self.L_x))
						tiled_site[:,1] = np.concatenate((self.site[:,1],
														  self.site[:,1]-self.L_y,
														  self.site[:,1]-self.L_y,
														  self.site[:,1]-self.L_y,
														  self.site[:,1],
														  self.site[:,1],
														  self.site[:,1]+self.L_y,
														  self.site[:,1]+self.L_y,
														  self.site[:,1]+self.L_y))
				
                                else:
						# Tile the sites -- NOTE: This is very inefficient! Later try using just a few critical edge sites.
						# However, this is the only option for creating boundaries using this library. (SSJ)
						tiled_site = np.zeros([9*self.N,2]);
						# Differently ordered tiling for easy indexing later  and to create a boundary effect
						tiled_site[:,0] = np.concatenate((self.site[:,0],
														  -self.site[:,0],
														  self.site[:,0],
														  2.*self.L_x - self.site[:,0],
														  -self.site[:,0],
														  2.*self.L_x - self.site[:,0],
														  -self.site[:,0],
														  self.site[:,0],
                                                                                                                  2.*self.L_x - self.site[:,0]))
						tiled_site[:,1] = np.concatenate((self.site[:,1],
														 - self.site[:,1],
														 - self.site[:,1],
														 - self.site[:,1],
														  self.site[:,1],
														  self.site[:,1],
														 2.*self.L_y - self.site[:,1],
														 2.*self.L_y - self.site[:,1],
														 2.*self.L_y - self.site[:,1]))

				# Call Voronoi Mesher
				voronoi = Voronoi(tiled_site,qhull_options='Qbb Qc');
			#	if self.is_periodic:
			#			voronoi = Voronoi(tiled_site,qhull_options='Qbb Qc');
			#	else:
			#			voronoi = Voronoi(self.site,qhull_options='Qbb Qc');
				
				self.voronoi = voronoi; # Store for use later when plotting
				
				
				if profile: tic = timeit.default_timer()



				# determine boundary cutting ridges and boundary points
			#	if not self.is_periodic:
			#		P1 = Point(0.,0.);
			#		P2 = Point(self.L_x,0.);
			#		P3 = Point(self.L_x,self.L_y);	
			#		P4 = Point(0.,self.L_y);	
		
			#		domain = Polygon(P1,P2,P3,P4);
	
			#		South = Segment(P1,P2);
			#		East = Segment(P2,P3);
			#		North = Segment(P3,P4);
			#		West = Segment(P4,P1);	

			#		#print voronoi.ridge_points;
			#	
			#		boundary_ridge_index = np.zeros(self.N, dtype=int);

			#		x = np.nditer(voronoi.ridge_points, flags=['multi_index'])
			#		while not x.finished:
			#			vertex_indices = voronoi.ridge_vertices[x.multi_index[0]];
			#			vertex1 = Point(voronoi.vertices[vertex_indices[0]][0],voronoi.vertices[vertex_indices[0]][1]);
			#			vertex2 = Point(voronoi.vertices[vertex_indices[1]][0],voronoi.vertices[vertex_indices[1]][1]);
			#			if (vertex_indices[0] == -1) and domain.encloses_point(vertex2):
			#				#print str(voronoi.vertices[vertex_indices[1],:]) + str(vertex_indices);			
			#				node = vertex2;			

			#				site_0 = Point(voronoi.points[voronoi.ridge_points[x.multi_index[0]][x.multi_index[1]]][0],voronoi.points[voronoi.ridge_points[x.multi_index[0]][x.multi_index[1]]][1]);
			#				site_1 = Point(voronoi.points[voronoi.ridge_points[x.multi_index[0]][1-x.multi_index[1]]][0],voronoi.points[voronoi.ridge_points[x.multi_index[0]][1-x.multi_index[1]]][1]);
			#				site_site_line = Line(site_0,site_1);
			#				boundary_vertex_line = site_site_line.perpendicular_line(node);

			#				#print "Yes"

			#				dist = float(10**6);
			#				if intersection(North,boundary_vertex_line):
			#					boundary_point = intersection(North,boundary_vertex_line)[0];
			#					self.xBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.x);	
			#					self.yBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.y);	
			#					#print str(float(boundary_point.x)) + " " + str(float(boundary_point.y));
			#					dist = node.distance(boundary_point);
			#					self.boundary_cutting_ridge[x[0],boundary_ridge_index[x[0]]] = 1.;	#North
			#				if intersection(South,boundary_vertex_line):
			#					boundary_point = intersection(South,boundary_vertex_line)[0];
			#					if node.distance(boundary_point) < dist:
			#						self.xBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.x);	
			#						self.yBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.y);	
			#						#print str(float(boundary_point.x)) + " " + str(float(boundary_point.y));
			#						dist = node.distance(boundary_point);
			#						self.boundary_cutting_ridge[x[0],boundary_ridge_index[x[0]]] = 2.; #South
			#				if intersection(East,boundary_vertex_line):
			#					boundary_point = intersection(East,boundary_vertex_line)[0];
			#					#print "yes east"
			#					#print str(float(boundary_point.x)) + " " + str(float(boundary_point.y));
			#					if node.distance(boundary_point) < dist:
			#						self.xBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.x);	
			#						self.yBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.y);	
			#						#print str(float(boundary_point.x)) + " " + str(float(boundary_point.y));
			#						dist = node.distance(boundary_point);
			#						self.boundary_cutting_ridge[x[0],boundary_ridge_index[x[0]]] = 3.; #East
			#				if intersection(West,boundary_vertex_line):
			#					boundary_point = intersection(West,boundary_vertex_line)[0];
			#					#print "yes west"
			#					#print str(float(boundary_point.x)) + " " + str(float(boundary_point.y));
			#					if node.distance(boundary_point) < dist:
			#						self.xBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.x);	
			#						self.yBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.y);	
			#						#print str(float(boundary_point.x)) + " " + str(float(boundary_point.y));
			#						dist = node.distance(boundary_point);
			#						self.boundary_cutting_ridge[x[0],boundary_ridge_index[x[0]]] = 4.; #West

			#				print str(self.xBoundary_points[x[0],boundary_ridge_index[x[0]]]) + " " + str(self.yBoundary_points[x[0],boundary_ridge_index[x[0]]]);
			#			
			#				#print self.boundary_cutting_ridge[boundary_ridge_index[x[0]]];
			#				
			#				#print boundary_ridge_index[x[0]];
			#				#boundary_ridge_index[x[0]] += 1;


			#		#	elif (not domain.encloses_point(vertex1)) or (not domain.encloses_point(vertex2)):
					#		#print str(voronoi.vertices[vertex_indices[1],:]) + str(vertex_indices);			
					#		if (not domain.encloses_point(vertex1)):
					#			node = vertex2;		
					#		elif (not domain.encloses_point(vertex2)):
					#			node = vertex1;
				
					#		site_0 = Point(voronoi.points[voronoi.ridge_points[x.multi_index[0]][x.multi_index[1]]][0],voronoi.points[voronoi.ridge_points[x.multi_index[0]][x.multi_index[1]]][1]);
					#		site_1 = Point(voronoi.points[voronoi.ridge_points[x.multi_index[0]][1-x.multi_index[1]]][0],voronoi.points[voronoi.ridge_points[x.multi_index[0]][1-x.multi_index[1]]][1]);
					#		site_site_line = Line(site_0,site_1);
					#		boundary_vertex_line = site_site_line.perpendicular_line(node);

					#		#print "Yes"

					#		dist = float(10**6);
					#		if intersection(North,boundary_vertex_line):
					#			boundary_point = intersection(North,boundary_vertex_line)[0];
					#			self.xBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.x);	
					#			self.yBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.y);	
					#			#print str(float(boundary_point.x)) + " " + str(float(boundary_point.y));
					#			dist = node.distance(boundary_point);
					#			self.boundary_cutting_ridge[x[0],boundary_ridge_index[x[0]]] = 1.;	#North
					#		if intersection(South,boundary_vertex_line):
					#			boundary_point = intersection(South,boundary_vertex_line)[0];
					#			if node.distance(boundary_point) < dist:
					#				self.xBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.x);	
					#				self.yBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.y);	
					#				#print str(float(boundary_point.x)) + " " + str(float(boundary_point.y));
					#				dist = node.distance(boundary_point);
					#				self.boundary_cutting_ridge[x[0],boundary_ridge_index[x[0]]] = 2.; #South
					#		if intersection(East,boundary_vertex_line):
					#			boundary_point = intersection(East,boundary_vertex_line)[0];
					#			#print "yes east"
					#			#print str(float(boundary_point.x)) + " " + str(float(boundary_point.y));
					#			if node.distance(boundary_point) < dist:
					#				self.xBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.x);	
					#				self.yBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.y);	
					#				#print str(float(boundary_point.x)) + " " + str(float(boundary_point.y));
					#				dist = node.distance(boundary_point);
					#				self.boundary_cutting_ridge[x[0],boundary_ridge_index[x[0]]] = 3.; #East
					#		if intersection(West,boundary_vertex_line):
					#			boundary_point = intersection(West,boundary_vertex_line)[0];
					#			#print "yes west"
					#			#print str(float(boundary_point.x)) + " " + str(float(boundary_point.y));
					#			if node.distance(boundary_point) < dist:
					#				self.xBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.x);	
					#				self.yBoundary_points[x[0],boundary_ridge_index[x[0]]] = float(boundary_point.y);	
					#				#print str(float(boundary_point.x)) + " " + str(float(boundary_point.y));
					#				dist = node.distance(boundary_point);
					#				self.boundary_cutting_ridge[x[0],boundary_ridge_index[x[0]]] = 4.; #West

					#		print str(self.xBoundary_points[x[0],boundary_ridge_index[x[0]]]) + " " + str(self.yBoundary_points[x[0],boundary_ridge_index[x[0]]]);
						
							#print self.boundary_cutting_ridge[boundary_ridge_index[x[0]]];
							
							#print boundary_ridge_index[x[0]];
				#			boundary_ridge_index[x[0]] += 1;



				#		x.iternext()
				
						

				# Set the connectivity

				# Number of neighbors - O(N)
				self.N_neighbor = np.zeros(self.N, dtype=int)
				for x in np.nditer(voronoi.ridge_points):
						if x < self.N:
								self.N_neighbor[x] += 1;

				if profile:
						toc = timeit.default_timer()
						print '0: '+'{:.2e}'.format((toc-tic)/self.N)+' s'
						tic = timeit.default_timer()


				# Neighbors - O(N)
				neighbor_index = np.zeros(self.N, dtype=int);
				for i in range(self.N):
						self.neighbor[i] = np.zeros(self.N_neighbor[i],dtype=int);
                                                self.boundary[i] = {};

				x = np.nditer(voronoi.ridge_points, flags=['multi_index'])
				while not x.finished:
						if x[0] < self.N:
								self.neighbor[x[0]][neighbor_index[x[0]]] = voronoi.ridge_points[x.multi_index[0]][1-x.multi_index[1]];
        
                                                                if not self.is_periodic:
                                                                     neighbor = self.neighbor[x[0]][neighbor_index[x[0]]]
                                                                     neighbor_point = voronoi.points[neighbor]
                                                                     if neighbor_point[0] < 0:
                                                                         self.boundary[x[0]][neighbor_index[x[0]]] = "West"
                                                                     elif neighbor_point[0] > self.L_x:
                                                                         self.boundary[x[0]][neighbor_index[x[0]]] = "East"
                                                                     elif neighbor_point[1] < 0:
                                                                         self.boundary[x[0]][neighbor_index[x[0]]] = "South"
                                                                     elif neighbor_point[1] > self.L_y:
                                                                         self.boundary[x[0]][neighbor_index[x[0]]] = "North"
                                                                     else:
                                                                         self.boundary[x[0]][neighbor_index[x[0]]] = "Internal"

								neighbor_index[x[0]] += 1;
						
						x.iternext()
				
				if profile:
						toc = timeit.default_timer()
						print '1: '+'{:.2e}'.format((toc-tic)/self.N)+' s'
						tic = timeit.default_timer()
			    

				# Calculate each of the properties sequentially (in individual for loops):
				# Only for first N sites!

				# Length - O(N)
				for i in range(self.N):
						self.length[i] = np.zeros([self.N_neighbor[i],2]);
						for j in range(self.N_neighbor[i]):
								self.length[i][j,:] = voronoi.points[self.neighbor[i][j]] - voronoi.points[i];

				if profile:
						toc = timeit.default_timer()
						print '2: '+'{:.2e}'.format((toc-tic)/self.N)+' s'
						tic = timeit.default_timer()


				# Face - O(N)
				neighbor_index = np.zeros(self.N,dtype=int);
				for i in range(self.N):
						self.face[i] = np.zeros(self.N_neighbor[i]);

				x = np.nditer(voronoi.ridge_points, flags=['multi_index'])
				while not x.finished:
						if x[0] < self.N:
								vertex_indices = voronoi.ridge_vertices[x.multi_index[0]];
								f = voronoi.vertices[vertex_indices[0],:] - voronoi.vertices[vertex_indices[1],:];
								self.face[x[0]][neighbor_index[x[0]]] = np.sqrt(f[0]**2 + f[1]**2)
								neighbor_index[x[0]] += 1;

						x.iternext()

				if profile:
						toc = timeit.default_timer()
						print '3: '+'{:.2e}'.format((toc-tic)/self.N)+' s'
						tic = timeit.default_timer()


				# Area & Centroid - O(N)
				for i in range(self.N):
						region_index = voronoi.point_region[i]; # Finds the region index of point i
						vertex_indices = voronoi.regions[region_index]; # Finds indices of vertices ordered around region of site i
						vertices = voronoi.vertices[vertex_indices]; # Listed sometimes clockwise, sometimes counter-clockwise...
						self.area[i] = 0.;
						N_vertices = len(vertex_indices);
						for j in range(N_vertices):
								self.area[i] += 0.5 * (vertices[j][0] * vertices[(j+1)%N_vertices][1] -
													   vertices[(j+1)%N_vertices][0] * vertices[j][1]);
						
						self.centroid[i,:] = 0.;
						for j in range(N_vertices):
								self.centroid[i,0] += ((vertices[j][0] + vertices[(j+1)%N_vertices][0]) *
																(vertices[j][0] * vertices[(j+1)%N_vertices][1] -
																 vertices[(j+1)%N_vertices][0] * vertices[j][1]));
								self.centroid[i,1] += ((vertices[j][1] + vertices[(j+1)%N_vertices][1]) *
																(vertices[j][0] * vertices[(j+1)%N_vertices][1] -
																 vertices[(j+1)%N_vertices][0] * vertices[j][1]));
				
						self.centroid[i,:] /= 6.*self.area[i];
						self.area[i] = np.abs(self.area[i]);

				if profile:
						toc = timeit.default_timer()
						print '4: '+'{:.2e}'.format((toc-tic)/self.N)+' s'
						tic = timeit.default_timer()
				
				
				# FaceCentre - O(N)
				neighbor_index = np.zeros(self.N, dtype=int);
				for i in range(self.N):
						self.face_center[i] = np.zeros((self.N_neighbor[i],2));

				x = np.nditer(voronoi.ridge_points, flags=['multi_index'])
				while not x.finished:
						if x[0] < self.N:
								vertex_indices = voronoi.ridge_vertices[x.multi_index[0]];
								f = 0.5 * (voronoi.vertices[vertex_indices[0],:] + voronoi.vertices[vertex_indices[1],:]);
								self.face_center[x[0]][neighbor_index[x[0]],:] = f[:] - self.site[x[0],:];
								neighbor_index[x[0]] += 1;
									
						x.iternext()

				if profile:
						toc = timeit.default_timer()
						print '5: '+'{:.2e}'.format((toc-tic)/self.N)+' s'
						tic = timeit.default_timer()


				# GradArea - O(N)
				for i in range(self.N):
						self.grad_area[i] = np.zeros((self.N_neighbor[i]+1,2));
						for j in range(self.N_neighbor[i]):
								X = self.site[i,:];
								Y = tiled_site[self.neighbor[i][j],:]; # Must be the periodic extension for the distances to work out
								XY = np.sqrt(np.sum(np.square(Y-X)));
								T = 0.5 * (X + Y) - self.face_center[i][j,:] - X; # Tangent vector
								self.grad_area[i][j,:] = self.face[i][j] * (0.5*(Y - X) + T) / XY

				
				
				# GradAreaT 
				for i in range(self.N):
						self.grad_area_t[i] = np.zeros((self.N_neighbor[i]+1,2));
						for j in range(self.N_neighbor[i]):
								neighbor_j = self.neighbor[i][j]%self.N; # Absolute index of the jth neighbour to i (Note wraps around)
                                                                #print neighbor_j;
								relative_i = np.where(self.neighbor[neighbor_j]%self.N == i)[0][0]; # Relative neighbour index of i from cell j
								self.grad_area_t[i][j,:] = self.grad_area[neighbor_j][relative_i,:];
						
						for j in range(self.N_neighbor[i]): # Find dAi/dXi (Self)
								self.grad_area_t[i][-1,:] -= self.grad_area[i][j,:];
						
						self.grad_area[i][-1,:] = self.grad_area_t[i][-1,:];

				if profile:
						toc = timeit.default_timer()
						print '6: '+'{:.2e}'.format((toc-tic)/self.N)+' s'
						tic = timeit.default_timer()

				# isBoundary - O(N)
                                if self.is_periodic:
				    for i in range(self.N):
						self.is_boundary[i] = False;
                                else:
				    for i in range(self.N):
                                        b = self.boundary[i]
					self.is_boundary[i] = False; 
                                        for j in range(self.N_neighbor[i]):
                                            if not b[j] == "Internal":
					        self.is_boundary[i] = True; 
                                   

				# Fix neighbours and sew edge connectivity together <3
				if self.is_periodic:
						for i in range(self.N):
								for j in range(self.N_neighbor[i]):
										self.neighbor[i][j] = self.neighbor[i][j]%self.N; # Well that's embarrassing

				if profile:
						toc = timeit.default_timer()
						print '7: '+'{:.2e}'.format((toc-tic)/self.N)+' s'
						tic = timeit.default_timer()


				metatoc = timeit.default_timer()
				print 'Meshing Complete: '+'{:.2e}'.format(metatoc-metatic)+' s'







