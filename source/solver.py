import numpy as np


class Data:

        def __init__(self,N):
                self.u_vel = np.zeros(N);
                self.v_vel = np.zeros(N);
                self.press = np.zeros(N);


def solve(data, D, L, G, Re, dt):
	# Solution for u by pressure projection
	n = len(data.u_vel)
	I = np.identity((n))
	# Semi implicit projection method
	# Solve for u_star
	A1 = I - L/(2.*Re)
	A2 = I + L/(2.*Re)
	rhs = -np.dot(G, data.press) + np.dot(A2, data.u_vel)
	u_star = np.linalg.solve(A1, rhs)
	# Pressure correction
	lhsPressure = np.dot(D, G)
	rhsPressure = np.dot(D, u_star)
	P_tild, res, ra, s = np.linalg.lstsq(lhsPressure, rhsPressure)
	# Update velocity and pressure
	GP = np.dot(G, P_tild)
	data.u_vel = u_star - dt*GP
	data.press = data.press + P_tild
	return data




