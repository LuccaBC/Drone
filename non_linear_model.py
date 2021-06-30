# Projeto de Controle Digital
# Parte 2 - movimento Z

from control.matlab import *
import control as ctrl
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig

def non_linear_sis():

	# Constantes

	w = 0
	z = 0
	z_dot = 0
	z_2dot = 0

	# Motor

	tau = (R*J)/(Km**2) 

	A = -(1/tau)

	B = -d/(eta*R**3*J)

	C = 1/(Km*tau)

	# Constantes de tempo 

	timestep = 0.02
	vec_time = np.linspace(0, tf, int(tf/timestep)+1)

	# Vetores

	vect_w = [w]
	vect_z = [z]
	vect_zdot = [z_dot]
	vect_z2dot = [z_2dot]

	# Loop

	for count_t, timestamp in enumerate(vec_time[1:]):

		Xbar = np.matmul(ssA-np.matmul(ssB,ssK)-np.matmul(ssL,ssC),Xbar)
		U = np.matmul(-K,Xbar)+Nbar*ref

		w_dot = A*w + B*w**2 + C*U

		# Integração de w

		w += w_dot*Ts 

		# Sistema

		z_2dot = -g -(d/m)*z_dot**2 + 3*(b/m)*w**2

		z_dot += z_2dot*Ts 

		z += z_2dot*Ts

		# Adicionando valores aos vetores

		vect_w.append(w)
		vect_z.append(z)
		vect_zdot.append(z_dot)
		vect_z2dot.append(z_2dot)

		xbar = [[z],[z_dot],[w]]

	return None