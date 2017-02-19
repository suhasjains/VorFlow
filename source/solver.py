import numpy as np

# Parameters
# N = 10
dt = 0.5
Re = 1000

class Data:

        def __init__(self,N):
                self.u_vel = np.zeros(N);
                self.v_vel = np.zeros(N);
                self.press = np.zeros(N);



def solve(data):
	# Solution for u by pressure projection
	# Hardcoded parameters for testing
	u = data.u_vel
	P = data.press
	n = len(u)
	D = np.zeros((n, n))
	L = np.zeros((n, n))
	G = np.zeros((n, n))
	I = np.identity((n))
	# Semi implicit projection method
	# Solve for u_star
	A1 = I - L/(2.*Re)
	A2 = I + L/(2.*Re)
	rhs = -np.dot(G, P) + np.dot(A2, u)
	u_star = np.linalg.solve(A1, rhs)
	# Pressure correction
	lhsPressure = np.dot(D, G)
	rhsPressure = np.dot(D, u_star)
	P_tild, res, ra, s = np.linalg.lstsq(lhsPressure, rhsPressure)
	# Update velocity and pressure
	GP = np.dot(G, P_tild)
	u = u_star - dt*GP
	P = P + P_tild
	print(u)
	print(P)
	return Data
N = 10
data = Data(N)
solve(data)
