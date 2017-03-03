import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.spatial import voronoi_plot_2d



def plot_mesh(mesh):
		# Plots just the mesh
		voronoi_plot_2d(mesh.voronoi)
		plt.xlabel(r'$x$')
		plt.ylabel(r'$y$')
		plt.title('Voronoi Mesh')
		plt.axis('scaled')
		plt.xlim([0.,mesh.L_x])
		plt.ylim([0.,mesh.L_y])
		plt.show()


def make_frame(mesh,field,name):
		# Plots the mesh and colours
		voronoi_plot_2d(mesh.voronoi)
		cmap = cm.get_cmap('RdBu')
		
		# Scale field
		field_scaled = (field - np.min(field)) / (np.max(field) - np.min(field))
		
		for i in range(9*mesh.N):
				region_index = mesh.voronoi.point_region[i] # Index of Voronoi region corresponding to site i
				region_vertex_index = mesh.voronoi.regions[region_index]
				if not -1 in region_vertex_index:
						polygon = [mesh.voronoi.vertices[k] for k in region_vertex_index]
						plt.fill(*zip(*polygon),c=cmap(field_scaled[i%mesh.N])) # Colours in the ith polygon with corresponding colour data from cmap(field)
		
		plt.xlabel(r'$x$')
		plt.ylabel(r'$y$')
		plt.title('Voronoi Mesh and '+name)
		plt.axis('scaled')
		plt.xlim([0.,mesh.L_x])
		plt.ylim([0.,mesh.L_y])
		plt.show()


