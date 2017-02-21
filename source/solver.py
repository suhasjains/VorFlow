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
		print i
		print "N",mesh.face[i]
		for j in range(mesh.N_neighbor[i]):
			k = mesh.neighbor[i][j]
			mfac = mesh.face[i][j]
			mlen = np.sqrt(mesh.length[i][j,0]**2 + mesh.length[i][j,1]**2)
			coeff = mfac/mlen
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
			

def solve(data, Dx, Dy, L, Gx, Gy, Re, dt):
	# Solution for u by pressure projection
	n = len(data.u_vel)
	I = np.identity((n))
	# Implicit projection method
	# Solve for u_star
	A1 = I - L/(2.*Re)
	A2 = I + L/(2.*Re)
	rhs_u = -0.5*dt*np.dot(Gx, data.press) + np.dot(A2, data.u_vel)
	rhs_v = -0.5*dt*np.dot(Gy, data.press) + np.dot(A2, data.v_vel)
	u_star = np.linalg.solve(A1, rhs_u)
	v_star = np.linalg.solve(A1, rhs_v)
	# Pressure correction
	lhsPressure_x = dt*np.dot(Dx, Gx)
	lhsPressure_y = dt*np.dot(Dy, Gy)
	rhsPressure_u = np.dot(Dx, u_star)
	rhsPressure_v = np.dot(Dy, v_star)
	lhsPressure = lhsPressure_x + lhsPressure_y
	rhsPressure = rhsPressure_u + rhsPressure_v
	P_tild, res, ra, s = np.linalg.lstsq(lhsPressure, rhsPressure)
	# Update velocity and pressure
	GPx = np.dot(Gx, P_tild)
	GPy = np.dot(Gy, P_tild)
	data.u_vel = u_star - 0.5*dt*GPx
	data.v_vel = v_star - 0.5*dt*GPy
	data.press = data.press + P_tild
	return data
