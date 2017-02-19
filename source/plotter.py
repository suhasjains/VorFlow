import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.spatial import voronoi_plot_2d



def plot_mesh(mesh):
		# Plots just the mesh
		voronoi_plot_2d(mesh.voronoi);
		plt.xlabel(r'$x$')
		plt.ylabel(r'$y$')
		plt.title('Voronoi Mesh')
		plt.xlim([0.,mesh.L_x]);
		plt.ylim([0.,mesh.L_y]);
		plt.axis('scaled');
		plt.show();


def make_frame(mesh,field,name):
		# Plots the mesh and colours
		voronoi_plot_2d(mesh.voronoi);
		for i in range(mesh.N):
				region_index = mesh.voronoi.point_region[i]; # Index of Voronoi region corresponding to site i
				region = mesh.voronoi.regions[region_index];
				if not -1 in region:
						polygon = [mesh.voronoi.vertices[k] for k in region]
						plt.fill(*zip(*polygon),c=field,cmap=plt.cm.RdBu); # Colours in the ith polygon with corresponding colour data from cmap(field)
		
		plt.xlabel(r'$x$')
		plt.ylabel(r'$y$')
		plt.title('Voronoi Mesh and '+name)
		plt.xlim([0.,mesh.L_x]);
		plt.ylim([0.,mesh.L_y]);
		plt.axis('scaled');
		plt.show();


