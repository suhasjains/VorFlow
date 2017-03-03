import numpy as np


class Data:

        def __init__(self,N):
                self.u_vel = np.zeros(N)
                self.v_vel = np.zeros(N)
                self.press = np.zeros(N)


def time_step(mesh):
	# Initialize matrices
	N = mesh.N
	Dx = np.zeros((N,N))
	Dy = np.zeros((N,N))
	L = np.zeros((N,N))
	Gx = np.zeros((N,N))
	Gy = np.zeros((N,N))
	# Populate matrices
	for i in range(N):
		#print i
		#print "N",mesh.face[i]
		for j in range(mesh.N_neighbor[i]):
			k = mesh.neighbor[i][j]
			mfac = mesh.face[i][j]
			mlen = np.sqrt(mesh.length[i][j,0]**2 + mesh.length[i][j,1]**2)
			coeff = mfac/mlen/mesh.area[i]
			# Div, x
			Dx[i][k] += coeff*(mesh.length[i][j,0]-mesh.face_center[i][j,0])
			# Div, y
			Dy[i][k] += coeff*(mesh.length[i][j,1]-mesh.face_center[i][j,1])
			# Laplacian
			L[i][k] += coeff
			L[i][i] -= coeff
			# Grad, x
			Gx[i][k] += coeff*mesh.face_center[i][j,0]
			# Grad, y
			Gy[i][k] += coeff*mesh.face_center[i][j,1]
	# If periodic, connectivity handled in mesh.py
	# If Dirichlet or Neumann, adjust BCs (to be done? if they are not handled by mesh.py)
	return Dx, Dy, L, Gx, Gy
			

def solve(data, Dx, Dy, L, Gx, Gy, dt, nu):
	solver_type = 'explicit'
	if (solver_type == 'proj'):
		# Solution for u by pressure projection
		n = len(data.u_vel)
		I = np.identity((n))
		# Implicit projection method
		# Solve for u_star
		A1 = I - nu*dt*L/(2.)
		A2 = I + nu*dt*L/(2.)
		rhs_u = -0.5*dt*np.dot(Gx, data.press) + np.dot(A2, data.u_vel)
		rhs_v = -0.5*dt*np.dot(Gy, data.press) + np.dot(A2, data.v_vel)
		u_star = np.linalg.solve(A1, rhs_u)
		v_star = np.linalg.solve(A1, rhs_v)
		# Vel star star
		u_star_star = u_star + dt*0.5*np.dot(Gx, data.press)
		v_star_star = v_star + dt*0.5*np.dot(Gy, data.press)
		# Pressure correction
		lhsPressure_x = dt*np.dot(Dx, Gx)
		lhsPressure_y = dt*np.dot(Dy, Gy)
		rhsPressure_u = 2.*np.dot(Dx, u_star_star)
		rhsPressure_v = 2.*np.dot(Dy, v_star_star)
		lhsPressure = lhsPressure_x + lhsPressure_y
		rhsPressure = rhsPressure_u + rhsPressure_v
		P_tild, res, ra, s = np.linalg.lstsq(lhsPressure, rhsPressure)
		# Update velocity and pressure
		GPx = np.dot(Gx, P_tild)
		GPy = np.dot(Gy, P_tild)
		data.u_vel = u_star_star - 0.5*dt*GPx
		data.v_vel = v_star_star - 0.5*dt*GPy
		data.press = P_tild
	elif (solver_type == 'explicit'):
		# Solution for u by pressure projection
                n = len(data.u_vel)
                I = np.identity((n))
		A1 = nu*L
		# Explicit
	 	# Pressure correction
                lhsPressure_x = np.dot(Dx, Gx)
                lhsPressure_y = np.dot(Dy, Gy)
		Hx = np.dot(A1, data.u_vel)
		Hy = np.dot(A1, data.v_vel)
                rhsPressure_u = np.dot(Dx, Hx)
                rhsPressure_v = np.dot(Dy, Hy)
                lhsPressure = lhsPressure_x + lhsPressure_y
                rhsPressure = rhsPressure_u + rhsPressure_v
                P_tild, res, ra, s = np.linalg.lstsq(lhsPressure, rhsPressure)
                #P_tild = np.linalg.solve(lhsPressure, rhsPressure)
		res_norm = np.linalg.norm(res)
		#print(res_norm)
		# Solve for u_star
                rhs_u = -dt*np.dot(Gx, data.press) + dt*np.dot(A1, data.u_vel)
                rhs_v = -dt*np.dot(Gy, data.press) + dt*np.dot(A1, data.v_vel)
                data.u_vel = data.u_vel + rhs_u
		data.v_vel = data.v_vel + rhs_v
		data.press = P_tild
	return data
