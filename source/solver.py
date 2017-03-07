import numpy as np
import scipy.sparse.linalg as lp
import scipy.sparse as sp
import timeit

class Data:

        def __init__(self,N):
                self.u_vel = np.zeros(N)
                self.v_vel = np.zeros(N)
                self.press = np.zeros(N)

def build_matrices(mesh):
	# Initialize matrices
	N = mesh.N
	Dx = sp.coo_matrix((N,N))
	Dy = sp.coo_matrix((N,N))
	L = sp.coo_matrix((N,N))
	Gx = sp.coo_matrix((N,N))
	Gy = sp.coo_matrix((N,N))
	type = 1
	# Populate matrices
	for i in range(N):
			#print i
			#print "N",mesh.face[i]
			for j in range(mesh.N_neighbor[i]):
					k = mesh.neighbor[i][j]
					# Comment out first, don't delete
					if type == 0:
							mfac = mesh.face[i][j]
							mlen = np.sqrt(mesh.length[i][j,0]**2 + mesh.length[i][j,1]**2)
							coeff = mfac/mlen/mesh.area[i]
							# Div, x
							Dx[i,k] = coeff*(mesh.length[i][j,0]-mesh.face_center[i][j,0])
							Dx[i,i] -= coeff*(mesh.length[i][j,0]-mesh.face_center[i][j,0])
							# Div, y
							Dy[i,k] = coeff*(mesh.length[i][j,1]-mesh.face_center[i][j,1])
							Dy[i,i] -= coeff*(mesh.length[i][j,1]-mesh.face_center[i][j,1])
							# Laplacian
							L[i,k] = coeff
							L[i,i] -= coeff
							# Grad, x
							Gx[i,k] = coeff*mesh.face_center[i][j,0]
							Gx[i,i] -= coeff*mesh.face_center[i][j,0]
							# Grad, y
							Gy[i,k] = coeff*mesh.face_center[i][j,1]
							Gy[i,i] -= coeff*mesh.face_center[i][j,1]
					## Potential new matrices
					else:
							mlen = np.sqrt(np.sum(np.square(mesh.length[i][j])))
							# Div, x
							Dx[i,k] = mesh.grad_area[i][j][0] / mesh.area[i]
							# Div, y
							Dy[i,k] = mesh.grad_area[i][j][1] / mesh.area[i]
							# Laplacian
							L[i,k] = mesh.face[i][j]/mlen/mesh.area[i]
							L[i,i] -= mesh.face[i][j]/mlen/mesh.area[i]
							## Grad, x
							Gx[i,k] = -mesh.grad_area_t[i][j][0] / mesh.area[i]
							# Grad, y
							Gy[i,k] = -mesh.grad_area_t[i][j][1] / mesh.area[i]
			if type == 1:
					Dx[i,i] += mesh.grad_area[i][-1][0] / mesh.area[i]
					Dy[i,i] += mesh.grad_area[i][-1][1] / mesh.area[i]
					Gx[i,i] -= mesh.grad_area_t[i][-1][0] / mesh.area[i]
					Gy[i,i] -= mesh.grad_area_t[i][-1][1] / mesh.area[i]
	# If periodic, connectivity handled in mesh.py
	# If Dirichlet or Neumann, adjust BCs (to be done? if they are not handled by mesh.py)
	# Sparsify matrices
	# (This may need to be done simultaneous with construction above if this procedure is memory-limited)
	Dx.tocsr()
	Dy.tocsr()
	L.tocsr()
	Gx.tocsr()
	Gy.tocsr()
	return Dx, Dy, L, Gx, Gy


def time_step(mesh,data,dt,nu):
	
		metatic = timeit.default_timer()
	
		Dx, Dy, L, Gx, Gy = build_matrices(mesh);
		N = mesh.N

		solver_type = 'CrankNicolson'

		if (solver_type == 'CrankNicolson'):
				# Solution for u by pressure projection
				#I = np.identity((N))
				I = sp.eye(N)
				# Implicit projection method
				# Solve for u_star
				A1 = I - nu*dt*L/(2.)
				A2 = I + nu*dt*L/(2.)
				#rhs_u = -0.5*dt*np.dot(Gx, data.press) + np.dot(A2, data.u_vel)
				rhs_u = -0.5*dt*Gx.dot(data.press) + A2.dot(data.u_vel)
				#rhs_v = -0.5*dt*np.dot(Gy, data.press) + np.dot(A2, data.v_vel)
				rhs_v = -0.5*dt*Gy.dot(data.press) + A2.dot(data.v_vel)
				#print(sp.sparse.issparse(A1))
				#u_star = np.linalg.solve(A1, rhs_u)
				#v_star = np.linalg.solve(A1, rhs_v)
				#print('1')
				u_star,B = lp.gmres(A1, rhs_u)
				#print('2')
				v_star,B = lp.gmres(A1, rhs_v)
				#print('3')
				# Vel star star
				#u_star_star = u_star + dt*0.5*np.dot(Gx, data.press)
				u_star_star = u_star + dt*0.5*Gx.dot(data.press)
				#v_star_star = v_star + dt*0.5*np.dot(Gy, data.press)
				v_star_star = v_star + dt*0.5*Gy.dot(data.press)
				#print('3a')
				# Pressure correction
				#lhsPressure_x = dt*np.dot(Dx, Gx)
				lhsPressure_x = dt*Dx*Gx
				#lhsPressure_y = dt*np.dot(Dy, Gy)
				lhsPressure_y = dt*Dy*Gy
				#rhsPressure_u = 2.*np.dot(Dx, u_star_star)
				rhsPressure_u = 2.*Dx.dot(u_star_star)
				#rhsPressure_v = 2.*np.dot(Dy, v_star_star)
				rhsPressure_v = 2.*Dy.dot(v_star_star)
				#print('3b')
				lhsPressure = lhsPressure_x + lhsPressure_y
				rhsPressure = rhsPressure_u + rhsPressure_v
				#print(sp.sparse.issparse(lhsPressure))
				#P_tild, res, ra, s = np.linalg.lstsq(lhsPressure, rhsPressure)
				#print('4')
				P_tild, B = lp.gmres(lhsPressure, rhsPressure)
				#print('5')
				# Update velocity and pressure
				#GPx = np.dot(Gx, P_tild)
				GPx = Gx.dot(P_tild)
				#GPy = np.dot(Gy, P_tild)
				GPy = Gy.dot(P_tild)
				data.u_vel = u_star_star - 0.5*dt*GPx
				data.v_vel = v_star_star - 0.5*dt*GPy
				data.press = P_tild
		elif (solver_type == 'FEuler'):
				# Solution for u by pressure projection
				#I = np.identity((N))
				I = sp.eye(N)
				A1 = nu*L
				# Explicit
				# Pressure correction
				#lhsPressure_x = np.dot(Dx, Gx)
				lhsPressure_x = Dx*Gx
				#lhsPressure_y = np.dot(Dy, Gy)
				lhsPressure_y = Dy*Gy
				#Hx = np.dot(A1, data.u_vel)
				Hx = A1.dot(data.u_vel)
				#Hy = np.dot(A1, data.v_vel)
				Hy = A1.dot(data.v_vel)
				#rhsPressure_u = np.dot(Dx, Hx)
				rhsPressure_u = Dx.dot(Hx)
				#rhsPressure_v = np.dot(Dy, Hy)
				rhsPressure_v = Dy.dot(Hy)
				lhsPressure = lhsPressure_x + lhsPressure_y
				rhsPressure = rhsPressure_u + rhsPressure_v
				#P_tild, res, ra, s = np.linalg.lstsq(lhsPressure, rhsPressure)
				P_tild, B = lp.gmres(lhsPressure, rhsPressure)
				#P_tild = np.linalg.solve(lhsPressure, rhsPressure)
				#res_norm = np.linalg.norm(res)
				#print(res_norm)
				# Solve for u_star
				#rhs_u = -dt*np.dot(Gx, data.press) + dt*np.dot(A1, data.u_vel)
				rhs_u = -dt*Gx.dot(data.press) + dt*A1.dot(data.u_vel)
				#rhs_v = -dt*np.dot(Gy, data.press) + dt*np.dot(A1, data.v_vel)
				rhs_v = -dt*Gy.dot(data.press) + dt*A1.dot(data.v_vel)
				data.u_vel = data.u_vel + rhs_u
				data.v_vel = data.v_vel + rhs_v
				data.press = P_tild

		elif solver_type == 'BEuler': # AW
				# Solve using Backward Euler (no dt restriction)
				# Inefficient, but robust.

				#I = np.identity((N))
				I = sp.eye(N)
				
				#u_star = np.linalg.solve(I - nu * dt * L, data.u_vel);
				u_star = lp.gmres(I - nu * dt * L, data.u_vel);
				#v_star = np.linalg.solve(I - nu * dt * L, data.v_vel);
				v_star = lp.gmres(I - nu * dt * L, data.v_vel);
			
				#Div = np.dot(Dx,u_star) + np.dot(Dy,v_star);
				Div = Dx.dot(u_star) + Dy.dot(v_star);
				#DG = np.dot(Dx,Gx) + np.dot(Dy,Gy); # == L ???  Nope...
				DG = Dx*Gx + Dy*Gy; # == L ???  Nope...
			
				#q = np.linalg.solve(DG,Div);
				q, B = lp.gmres(DG,Div);
			
				#data.u_vel = u_star - np.dot(Gx,q);
				data.u_vel = u_star - Gx.dot(q);
				#data.v_vel = v_star - np.dot(Gy,q);
				data.v_vel = v_star - Gy.dot(q);
		
		metatoc = timeit.default_timer()
		print 'Solving Complete: '+'{:.2e}'.format(metatoc-metatic)+' s'
		
		return data
