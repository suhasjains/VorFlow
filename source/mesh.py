import numpy as np

class Mesh(N,L_x,L_y):
	
	def __init__(self):
		self.N = N; # Don't want globals in mesh.py otherwise prohibitive to module testing!
		self.L_x = L_x;
		self.L_y = L_y;
		
		
	# properties:
		# Domain limits Lx, Ly...



def initialize_mesh(N,Lx,Ly):
	# Constructor


def generate_mesh():
	# Base mesh generator
	
	# Call scipy meshing
	
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
	
	
	
	

	



def update_mesh(data):










