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
	Dx = sp.dok_matrix((N,N))
	Dy = sp.dok_matrix((N,N))
	L = sp.dok_matrix((N,N))
	Gx = sp.dok_matrix((N,N))
	Gy = sp.dok_matrix((N,N))
	type = 1
	# Populate matrices
	for i in range(N):
			#print i
			#print "N",mesh.face[i]
			for j in range(mesh.N_neighbor[i]):
					k = mesh.neighbor[i][j]
					if type == 0:
							mfac = mesh.face[i][j]
							mlen = np.sqrt(mesh.length[i][j,0]**2 + mesh.length[i][j,1]**2)
							coeff = mfac/mlen/mesh.area[i]
							# Div, x
							if k < N:
								Dx[i,k] += coeff*(mesh.length[i][j,0]-mesh.face_center[i][j,0])
							Dx[i,i] -= coeff*(mesh.length[i][j,0]-mesh.face_center[i][j,0])
							# Div, y
							if k < N:
								Dy[i,k] += coeff*(mesh.length[i][j,1]-mesh.face_center[i][j,1])
							Dy[i,i] -= coeff*(mesh.length[i][j,1]-mesh.face_center[i][j,1])
							# Laplacian
							if k < N:
								L[i,k] += coeff
							L[i,i] -= coeff
							# Grad, x
							if k < N:
								Gx[i,k] += coeff*mesh.face_center[i][j,0]
							Gx[i,i] -= coeff*mesh.face_center[i][j,0]
							# Grad, y
							if k < N:
								Gy[i,k] += coeff*mesh.face_center[i][j,1]
							Gy[i,i] -= coeff*mesh.face_center[i][j,1]
					else:
							mlen = np.sqrt(np.sum(np.square(mesh.length[i][j])))
							if k < N:
							# Div, x
								Dx[i,k] += mesh.grad_area[i][j][0] / mesh.area[i]
							# Div, y
								Dy[i,k] += mesh.grad_area[i][j][1] / mesh.area[i]
							# Laplacian
								L[i,k] += mesh.face[i][j]/mlen/mesh.area[i]
							L[i,i] -= mesh.face[i][j]/mlen/mesh.area[i]
							## Grad, x
							if k < N:
								Gx[i,k] -= mesh.grad_area_t[i][j][0] / mesh.area[i]
							# Grad, y
								Gy[i,k] -= mesh.grad_area_t[i][j][1] / mesh.area[i]
			if type == 1:
					Dx[i,i] += mesh.grad_area[i][-1][0] / mesh.area[i]
					Dy[i,i] += mesh.grad_area[i][-1][1] / mesh.area[i]
					Gx[i,i] -= mesh.grad_area_t[i][-1][0] / mesh.area[i]
					Gy[i,i] -= mesh.grad_area_t[i][-1][1] / mesh.area[i]
	# If periodic, connectivity handled in mesh.py
	# If Dirichlet or Neumann, adjust BCs (to be done? if they are not handled by mesh.py)
	# Sparsify matrices
	# (This may need to be done simultaneous with construction above if this procedure is memory-limited)
	Dx = sp.csr_matrix(Dx)
	Dy = sp.csr_matrix(Dy)
	L = sp.csr_matrix(L)
	Gx = sp.csr_matrix(Gx)
	Gy = sp.csr_matrix(Gy)
	
	return Dx, Dy, L, Gx, Gy



def build_rhs(mesh,u,v,p,BCu,BCuvals,BCv,BCvvals):
	# Initialize vectors
	N = mesh.N
	rhsDxu = np.zeros(N) 
	rhsDyv = np.zeros(N) 
	rhsLu = np.zeros(N) 
	rhsLv = np.zeros(N)
	rhsLp = np.zeros(N) 
	rhsGxp = np.zeros(N) 
	rhsGyp = np.zeros(N) 
	type = 1
	uHold = 0
	vHold = 0
	pHold = 0
	if ~mesh.is_periodic:
		for i in range(N):
				for j in range(mesh.N_neighbor[i]):
					# Correct ghost cell velocities and pressures (BC order West, East, South, North) (1=Dirichlet, 2=Neumann)
					if mesh.boundary[i][j] == "West":
						if BCu[0] == 1:
							uHold = 2*BCuvals[0]-u[i]
							pHold = p[i]
						#else if BCu[0] == 2:
						if BCv[0] == 1:
							vHold = 2*BCvvals[0]-v[i]
						#else if BCv[0] == 2:
					elif mesh.boundary[i][j] == "East":
						if BCu[1] == 1:
							uHold = 2*BCuvals[1]-u[i]
							pHold = p[i]
						#else if BCu[1] == 2:
						if BCv[1] == 1:
							vHold = 2*BCvvals[1]-v[i]
						#else if BCv[1] == 2:
					elif mesh.boundary[i][j] == "South":
						if BCu[2] == 1:
							uHold = 2*BCuvals[2]-u[i]
						#else if BCu[2] == 2:
						if BCv[2] == 1:
							vHold = 2*BCvvals[2]-v[i]
							pHold = p[i]
						#else if BCv[2] == 2:
					elif mesh.boundary[i][j] == "North":
						if BCu[3] == 1:
							uHold = 2*BCuvals[3]-u[i]
						#else if BCu[3] == 2:
						if BCv[3] == 1:
							vHold = 2*BCvvals[3]-v[i]
							pHold = p[i]
						#else if BCv[3] == 2:
					# Now the vectors
					if type == 0:
							mfac = mesh.face[i][j]
							mlen = np.sqrt(mesh.length[i][j,0]**2 + mesh.length[i][j,1]**2)
							coeff = mfac/mlen/mesh.area[i]
							# Div, x
							rhsDxu[i] += coeff*(mesh.length[i][j,0]-mesh.face_center[i][j,0])*uHold
							# Div, y
							rhsDyv[i] += coeff*(mesh.length[i][j,1]-mesh.face_center[i][j,1])*vHold
							# Laplacian
							rhsLu[i] += coeff*uHold
							rhsLv[i] += coeff*vHold
							rhsLp[i] += coeff*pHold
							# Grad, x
							rhsGxp[i] += coeff*mesh.face_center[i][j,0]*pHold
							# Grad, y
							rhsGyp[i] += coeff*mesh.face_center[i][j,1]*pHold
					else:
							mlen = np.sqrt(np.sum(np.square(mesh.length[i][j])))
							# Div, x
							rhsDxu[i] += mesh.grad_area[i][j][0] / mesh.area[i] * uHold
							# Div, y
							rhsDyv[i] += mesh.grad_area[i][j][1] / mesh.area[i] * vHold
							# Laplacian
							rhsLu[i] += mesh.face[i][j]/mlen/mesh.area[i] * uHold
							rhsLv[i] += mesh.face[i][j]/mlen/mesh.area[i] * vHold
							rhsLp[i] += mesh.face[i][j]/mlen/mesh.area[i] * pHold
							## Grad, x
							rhsGxp[i] -= mesh.grad_area_t[i][j][0] / mesh.area[i] * pHold
							# Grad, y
							rhsGyp[i] -= mesh.grad_area_t[i][j][1] / mesh.area[i] * pHold
	
	return rhsDxu, rhsDyv, rhsLu, rhsLv, rhsLp, rhsGxp, rhsGyp



def time_step(mesh,data,dt,nu,BCu=[0,0,0,0],BCuvals=[0,0,0,0],BCv=[0,0,0,0],BCvvals=[0,0,0,0]):
	
		Dx, Dy, L, Gx, Gy = build_matrices(mesh);
		N = mesh.N
	
		metatic = timeit.default_timer()	

		solver_type = 'BEuler'

		if (solver_type == 'CrankNicolson'):
				rhsDxu, rhsDyv, rhsLu, rhsLv, rhsLp, rhsGxp, rhsGyp = build_rhs(mesh,data.u_vel,data.v_vel,data.press,BCu,BCuvals,BCv,BCvvals)
				# Solution for u by pressure projection
				#I = np.identity((N))
				I = sp.eye(N)
				# Implicit projection method
				# Solve for u_star
				A1 = I - nu*dt*L/(2.)
				A2 = I + nu*dt*L/(2.)
				#rhs_u = -0.5*dt*np.dot(Gx, data.press) + np.dot(A2, data.u_vel)
				rhs_u = -0.5*dt*Gx.dot(data.press) + A2.dot(data.u_vel) - 0.5*dt*rhsGxp + nu*dt*rhsLu
				#rhs_v = -0.5*dt*np.dot(Gy, data.press) + np.dot(A2, data.v_vel)
				rhs_v = -0.5*dt*Gy.dot(data.press) + A2.dot(data.v_vel) - 0.5*dt*rhsGyp + nu*dt*rhsLv
				#print(sp.sparse.issparse(A1))
				#u_star = np.linalg.solve(A1, rhs_u)
				#v_star = np.linalg.solve(A1, rhs_v)
				#print('1')
				#u_star,B = lp.gmres(A1, rhs_u)
				u_star = lp.spsolve(A1, rhs_u)
				#print('2')
				#v_star,B = lp.gmres(A1, rhs_v)
				v_star = lp.spsolve(A1, rhs_v)
				#print('3')
				# Vel star star
				#u_star_star = u_star + dt*0.5*np.dot(Gx, data.press)
				u_star_star = u_star + dt*0.5*Gx.dot(data.press) + dt*0.5*rhsGxp
				#v_star_star = v_star + dt*0.5*np.dot(Gy, data.press)
				v_star_star = v_star + dt*0.5*Gy.dot(data.press) + dt*0.5*rhsGyp
				#print('3a')
				rhsDxu, rhsDyv, rhsLu, rhsLv, rhsLp, rhsGxp, rhsGyp = build_rhs(mesh,u_star_star,v_star_star,data.press,BCu,BCuvals,BCv,BCvvals)
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
				rhsPressure = rhsPressure_u + rhsPressure_v + 2.*rhsDxu + 2.*rhsDyv - dt*rhsLp
				#print(sp.sparse.issparse(lhsPressure))
				#P_tild, res, ra, s = np.linalg.lstsq(lhsPressure, rhsPressure)
				#print('4')
				#P_tild, B = lp.gmres(lhsPressure, rhsPressure)
				P_tild = lp.spsolve(lhsPressure, rhsPressure)
				#print('5')
				rhsDxu, rhsDyv, rhsLu, rhsLv, rhsLp, rhsGxp, rhsGyp = build_rhs(mesh,u_star_star,v_star_star,P_tild,BCu,BCuvals,BCv,BCvvals)
				# Update velocity and pressure
				#GPx = np.dot(Gx, P_tild)
				GPx = Gx.dot(P_tild)
				#GPy = np.dot(Gy, P_tild)
				GPy = Gy.dot(P_tild)
				data.u_vel = u_star_star - 0.5*dt*GPx - 0.5*dt*rhsGxp
				data.v_vel = v_star_star - 0.5*dt*GPy - 0.5*dt*rhsGyp
				data.press = P_tild
		elif (solver_type == 'FEuler'):
				# NON PERIODIC BC NOT IMPLEMENTED
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
				#P_tild, B = lp.gmres(lhsPressure, rhsPressure)
				P_tild = lp.spsolve(lhsPressure, rhsPressure)
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
				rhsDxu, rhsDyv, rhsLu, rhsLv, rhsLp, rhsGxp, rhsGyp = build_rhs(mesh,data.u_vel,data.v_vel,data.press,BCu,BCuvals,BCv,BCvvals)
				
				I = sp.eye(N,format='csr')
	
				u_star = lp.spsolve(I - nu * dt * L, data.u_vel + nu*dt*rhsLu)
				v_star = lp.spsolve(I - nu * dt * L, data.v_vel + nu*dt*rhsLv)
				
				rhsDxu, rhsDyv, rhsLu, rhsLv, rhsLp, rhsGxp, rhsGyp = build_rhs(mesh,u_star,v_star,data.press,BCu,BCuvals,BCv,BCvvals)
				
				#Div = Dx.dot(u_star) + Dy.dot(v_star) + rhsDxu + rhsDyv - rhsLp
				Div = Dx.dot(u_star) + Dy.dot(v_star) 
				DG = Dx*Gx + Dy*Gy; # == L ???  Nope...
				q = lp.spsolve(DG,Div);
				
				rhsDxu, rhsDyv, rhsLu, rhsLv, rhsLp, rhsGxp, rhsGyp = build_rhs(mesh,u_star,v_star,q,BCu,BCuvals,BCv,BCvvals)
				data.u_vel = u_star - Gx.dot(q) - rhsGxp;
				data.v_vel = v_star - Gy.dot(q) - rhsGyp;
				
				try:
					data.tracer;
				except AttributeError:
					toodle=0;
				else: # Well then diffuse the tracer already
					nuTracer = 1.e-5
					Visc = I - nuTracer * dt * L;
					#VDivX = dt * sp.csr_matrix((data.u_vel,(range(N),range(N))))*Dx; 
					#VDivY = dt * sp.csr_matrix((data.v_vel,(range(N),range(N))))*Dy;
					data.tracer = lp.spsolve(Visc, data.tracer) # Just diffusive
				
				data.press = q/dt;
		elif solver_type == 'BEulerStab': #AW
				# Solves with STABILISED Projection method (eq 5.9 Borgers)
				
				I = sp.eye(N,format='csr')
				
				u_star = lp.spsolve(I - nu * dt * L, data.u_vel);
 				v_star = lp.spsolve(I - nu * dt * L, data.v_vel);
				
				# Stable iterative projection:
				rho = np.abs(lp.eigs(L,1,which='LM',tol=0.01)[0][0]) # Spectral radius
				#print rho
				omega = 0.8/rho;
				
				tol = 1e-8; # Relative tolerance for divergence
				u_J = u_star;
				v_J = v_star;
				q_J = u_star * 0.;
				i_outer = 0;
				k = 0;
				while np.linalg.norm( Dx.dot(u_J) + Dy.dot(v_J), 2) > tol and i_outer < 1000:
					# Iterate until uJ ~ Pu
					i_outer += 1;
					k = 0;
					u_jk = u_J;
					v_jk = v_J;
					
					q_jkP1 = lp.spsolve(L,Dx.dot(u_jk)+Dy.dot(v_jk));
					u_jkP1 = u_jk - Gx.dot(q_jkP1);
					v_jkP1 = v_jk - Gy.dot(q_jkP1);
					
					while (np.sqrt(np.linalg.norm( u_jkP1 ,2)**2 + np.linalg.norm( v_jkP1 ,2)**2)) >= (np.sqrt(np.linalg.norm( u_J ,2)**2 + np.linalg.norm( v_J ,2)**2) + 1.e-14) and k < 1000:
						k += 1;
						#print k
						#print (np.linalg.norm( u_jkP1 ,2) - np.linalg.norm( u_J ,2))
						#print (np.linalg.norm( v_jkP1 ,2) - np.linalg.norm( v_J ,2))
						#print (np.sqrt(np.linalg.norm( u_jkP1 ,2)**2 + np.linalg.norm( v_jkP1 ,2)**2) - np.sqrt(np.linalg.norm( u_J ,2)**2 + np.linalg.norm( v_J ,2)**2))
						u_jk = u_jkP1;
						v_jk = v_jkP1;
						q_jk = q_jkP1;
						
						
						q_jkP1 = L.dot(q_jk) * omega;
						u_jkP1 = u_jk - Gx.dot(q_jkP1);
						v_jkP1 = v_jk - Gy.dot(q_jkP1);
					
					#print k
					u_J = u_jkP1;
					v_J = v_jkP1;
					q_J = q_jkP1;
				
				if i_outer >= 999 or k >= 999:
					print np.linalg.norm( Dx.dot(u_J) + Dy.dot(v_J), 2)
				
				data.u_vel = u_J;
				data.v_vel = v_J;
				
				try:
					data.tracer;
 				except AttributeError:
 					toodle=0;
 				else: # Well then diffuse the tracer already
 					nuTracer = 1.e-8;
 					Visc = I - nuTracer * dt * L;
 					data.tracer = lp.spsolve(Visc, data.tracer) # Just diffusive
 				
				data.press = q_J/dt;
		
		metatoc = timeit.default_timer()
		print 'Solving Complete: '+'{:.2e}'.format(metatoc-metatic)+' s'
		
		return data
