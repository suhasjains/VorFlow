import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.spatial import voronoi_plot_2d



def plot_mesh(mesh,ax):
		## Plots just the mesh
		make_frame(mesh,np.zeros(mesh.N),'Nothing More',ax,True)


def make_frame(mesh,field,name,ax,plotMesh=False):
		## Plots the mesh and colours
		
		# Remove previous plot
		fig = plt.gcf();
		ax.clear()
		# Plots the mesh and colours
		if plotMesh: voronoi_plot_2d(mesh.voronoi,ax)
		
		cmap = cm.get_cmap('RdBu')
		for i in range(mesh.N):
				region_index = mesh.voronoi.point_region[i] # Index of Voronoi region corresponding to site i
				region_vertex_index = mesh.voronoi.regions[region_index]
				if not -1 in region_vertex_index:
						polygon = [mesh.voronoi.vertices[k] for k in region_vertex_index]
						im=plt.fill(*zip(*polygon),color=cmap(field[i%mesh.N])) # Colours in the ith polygon with corresponding colour data from cmap(field)
		
		plt.xlabel(r'$x$')
		plt.ylabel(r'$y$')
		plt.title('Voronoi Mesh and '+name)
		plt.axis('scaled')
		plt.xlim([0.,mesh.L_x])
		plt.ylim([0.,mesh.L_y])
		mappable = cm.ScalarMappable(cmap=cmap)
		mappable.set_array(field)
		if len(fig.axes) > 1:
				cax = fig.axes[1]
				plt.colorbar(mappable,ax=ax,cax=cax)
		else:
				plt.colorbar(mappable,ax=ax)

		plt.show()
		plt.pause(1e-6)






