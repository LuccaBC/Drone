# Projeto de Controle Digital
# Parte 1 - movimento Z

from control.matlab import *
import control as ctrl
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig
from non_linear_model import * 

#Enviroment input
#São Carlos -----------
#height = 856 #meters above seal level
#latitude = -22.0104
#humidity = 0.39

###São José do Rio Preto --------------
#height = 489
#latitude =  -20.8202
#humidity = 0.3
##
###Ribeirão Preto ---------------
#height = 546.8
#latitude =  -21.1036
#humidity = 0.36
##
###Recife---------------
height = 4
latitude =  -8.0403
humidity = 0.85

#Atmospheric calculation

#constants
vapor_density = humidity*12.83/1000 #kg/m^3
M_dry = 0.0289652
M_water = 0.018016
R_gas = 8.31446
T0 = 288.15
L = 0.0065
p0 = 101325

#equations
temp = T0 - L*height
lat_r = np.deg2rad(latitude)
g = 9.780327*(1+0.0053024*(np.sin(lat_r))**2 - 0.0000058*(np.sin(2*lat_r))**2)
#Buck equation for Vapour pressure
p_vapour = 0.61121*np.exp((18.678 - (temp-273.15)/234.5)*((temp-273.15)/(257.14+(temp-273.15))))
#dry pressure
p_dry = p0*(1 - L*height/T0)**(9.81*M_dry/(R_gas*L))

rho_dry = p_dry*M_dry/(R_gas*temp)
hum_mass = vapor_density/rho_dry
partial_dry = (1 - hum_mass) * p_dry
partial_wet = hum_mass * p_vapour * 100000
air_density = (partial_dry*M_dry+partial_wet*M_water)/(R_gas*temp)

# Sistema

# Constantes do problema

diameter = 11 #inch
radius = diameter*25.4/2000
Ts = 0.02  # Amostragem
w_max = 9380 # Vel. max. rpm
c_tu = 0.1059
A = np.pi*radius**2
A_drone = 0.021
m = 1
Ohm0 = np.pi/30*(0.4*w_max) # convertido para rad/s
ct_br = (8/(np.pi ** 3))*c_tu
cd = 0.8
dz0 = 0.1

# Parametros motor

km = 12.1*(10**(-3))
Cq = 0.8
B = 5.3*(10**(-5))
J = 1*(10**(-4))
Rmotor = 0.275
c_pu = 0.0407
cq_br = (8/(np.pi ** 4))*c_pu

# Constantes do ss

b = (1/2)*(ct_br*air_density*A*radius**2)/m
d_corpo = (1/2)*(cd*air_density*A_drone)
d_motor = (1/2)*(cq_br*air_density*A*radius**3)

# Constantes motor

Am = -(km**2/(Rmotor*J) + (2*d_motor*Ohm0)/J)
Bm = km/(Rmotor*J)

# Constantes dos estados
const_x1 = (d_corpo*dz0**2)/m
const_x2 = (3*b*Ohm0**2)/m +g
const_x3 = ((d_motor*Ohm0**2)/J)

vet_const = [const_x1, const_x2, const_x1]
print("Veotr de constantes:\n", vet_const)
# Espaço de estados - Continuo

Ac = [[0,1,0],[0,(-2*d_corpo*dz0)/m,(6*b*Ohm0)/m],[0,0,Am]]
Bc = [[0],[0],[Bm]]
Cc = [1,0,0]
Dc = [0]

#Atraso z**-1
#Ac = [[0,1,0,0],[0,(-2*d_corpo*dz0)/m,(6*b*Ohm0)/m,0],[0,0,Am,0],[0,0,0,0]]
#Bc = [[0],[0],[Bm],[1]]

#constant_sis = -(6*ct_br*air_density*A*((w0*radius)**2)/m + g) 
#print("Constante:", constant_sis)

cont_ss = ss(Ac,Bc,Cc,Dc)
print("Sistema Continuo:\n", cont_ss)

# Discretização

disc_ss = c2d(cont_ss, Ts)
print("Sistema Discreto:\n", disc_ss)


print(disc_ss.pole())
# Teste

yout, T = step(cont_ss)
plt.step(T, yout, where='post')
plt.show()

# Vel
C_vel = [0,1,0]
cont_ss_vel = ss(Ac,Bc,C_vel,Dc)

yout_vel, T_vel = step(cont_ss_vel)
plt.step(T_vel, yout_vel, where='post')
plt.show()

# Controlador
# Ms = 25%
# Ts de 1% = 10s
# Acc max = 10,2 cm/s^2

zeta = 0.4037
#zeta = 0.5037
ts_1 = 10

wn = 4.6/(zeta*ts_1)
wd = wn*np.sqrt(1-zeta**2)

# polos continuos
p_c1 = complex(-zeta*wn,wd)
p_c2 = complex(-zeta*wn,-wd)
print("polos Continuos:\n", np.round(p_c1,4), "e" , np.round(p_c2,4))

# polos discretos
p_d1 = np.exp(p_c1*Ts)
p_d2 = np.exp(p_c2*Ts)
print("polos Discretos:\n", np.round(p_d1,4), "e", np.round(p_d2,4))

# Valor inicial - gravidade medida pela IMU ?
X0 = [[-1], [0], [0],[0]]

#termo de "ruido" do resto da expansão de taylor - ficou bem ruim
aux_c = [0,0,0,0]
cM = np.array(aux_c)
mT = np.eye(4)
#new_ss, mT = ctrl.reachable_form(disc_ss)
new_ss = disc_ss
new_ss_aux = np.concatenate((disc_ss.A.A,[[0],[0],[0]]),axis = 1).T
new_ss_atz_A = np.concatenate((new_ss_aux,[[0],[0],[0],[0]]),axis = 1).T
new_ss_atz_B = np.concatenate((new_ss.B.T,[[1]]),axis = 1).T
init_state = np.matmul(mT, X0)
# colocando os polos
K = place(new_ss_atz_A, new_ss_atz_B, [p_d1, p_d2, 0.6309566, 0.62]) # Controlabilidade forma canonica
PHI = new_ss_atz_A - np.matmul(new_ss_atz_B,K)
print('K=', np.round(K, 2))

#step controle sem ref
final_time = 10
offset = 0
print('Matriz nova controle \n', PHI)
cl_ss = ss(PHI, [[0], [0], [0],[0]], [1,0,0,0], new_ss.D, Ts)
print('A=', np.round(cl_ss.A, 2))
time_d = linspace(0, int(final_time/Ts)*Ts, int(final_time/Ts)+1)
yout, tout, xout = initial(cl_ss, time_d, init_state, return_x=True)
xx = np.matmul(np.linalg.inv(mT), xout.T)
#print('xx=\n', xx)
plt.step(tout, (xx[0, :].T)+offset, 'r', where='post', label='modal')
plt.step(tout, (xx[1, :].T)+offset, 'b', where='post', label='modal')
plt.show()

'''
#step normal do controlado sem ref
yout, T = step(cl_ss, time_d, init_state)
plt.step(T, yout, where='post')
plt.show()
'''

#Observador
#disc_ss
p_o1 = np.exp(10*p_c1*Ts)
p_o2 = np.exp(10*p_c2*Ts)
X0 = [[0], [0], [0]]
Ltr = place(disc_ss.A.T, disc_ss.C.T,[p_o1, p_o2, 0.6309566])
L = Ltr.T
newA = disc_ss.A - np.matmul(L,disc_ss.C)

print('matriz L do controlador: \n', L)
#pondo tudo juntao
#deve-se expandir a matriz 
K = place(disc_ss.A, disc_ss.B,[p_d1, p_d2, 0.6309566])
all_A = np.concatenate((disc_ss.A.A,-np.matmul(disc_ss.B,K)), axis=1)
aux_A = np.concatenate((np.matmul(L, disc_ss.C),newA - np.matmul(disc_ss.B,K)),axis=1)
all_A = np.concatenate((all_A, aux_A), axis=0)
print('print A com controlador :\n', all_A)
#####################################
#novo sistema somente com o observador
teste_ss = ss(all_A,np.zeros([6,1]),[1,0,0,0,0,0],[0],Ts)
#######################################
#acrescentando a ref
all_up = np.concatenate((disc_ss.A.A, disc_ss.B.A), axis=1)
all_dw = np.concatenate((disc_ss.C.A,disc_ss.D.A),  axis=1)
all_2 = np.concatenate((all_up, all_dw), axis=0)
Nrf = np.matmul(np.linalg.inv(all_2), [[1], [0], [0], [1]]) # -> [[1], [0], [0], [1]]
intall2 = np.linalg.inv(all_2)

#para provar mal condicionamento da matriz
condicional_M = np.linalg.norm(all_2,2)*np.linalg.norm(intall2,2)
#como foi mal condicionada pode-se utiilizar o matmul em ->[[1], [0], [0], [1]]


print("inversa da A:\n",np.linalg.inv(all_2))
print('Nrf\n',Nrf,'\nall_2\n',all_2)

print('inv da all_2', np.round(intall2,1))
Nx = Nrf[0:3]
Nu = Nrf[3]
N =  Nu + np.matmul(K, Nx)
#N = -np.linalg.inv(np.matmul([1,0,0], np.matmul(np.linalg.inv(disc_ss.A - np.matmul(disc_ss.B,K)),disc_ss.B)))
print("N:\n",N)
auxB = np.matmul(disc_ss.B, N)
newB = np.concatenate((auxB, auxB), axis=0)
print('matriz auxB\n', auxB)
aug_ss2 = ss(all_A, newB, [1,0,0,0,0,0], [0], Ts)
print('printando a bagaça final\n',aug_ss2)

mT_c = np.eye(6)
X0_c = [[0], [0], [0], [0], [0], [0]]
init_state_c = np.matmul(mT_c, X0_c)
final_time_c = 10
time_d_c = linspace(0, int(final_time_c/Ts)*Ts, int(final_time_c/Ts)+1)

yout_cc, Tcc = step(aug_ss2, time_d_c, init_state_c)
plt.step(Tcc, yout_cc, where='post')
plt.show()

'''
#grafico somente do controle com observador
yout, tout, xout = initial(teste_ss, time_d_c, init_state_c, return_x=True)
xx_c = np.matmul(np.linalg.inv(mT_c), xout.T)
#print('xx=\n', xx)
plt.step(tout, (xx_c[0, :].T), 'r', where='post', label='modal')
plt.show()
'''

#info do step
info_do_step_ref = stepinfo(aug_ss2)
print('stepinfo aug_ss2:\n',info_do_step_ref)
#FT do sistema discrtizado
G_c = ss2tf(aug_ss2)



'''
#Acrescentando ruido do barometro
#vetor de ruido ate 0.5
Bnoise = np.random.rand(1,len(time_d_c))
mT_n = np.eye(6)
X0_n = [[0], [0], [0], [0], [0], [0]]
init_state_n = np.matmul(mT_n, X0_n)
final_time_n = 10
yout_noise, Tnoise = step(aug_ss2, time_d_c, init_state_n)
plt.step(Tnoise, yout_noise, where='post')
plt.show()
'''

#Parte 1.2

# Constantes

w = 0
z = 0
z_dot = 0
z_2dot = 0
U = 0

# Motor

tau = (Rmotor*J)/(km**2) 

A = -(1/tau)

B = -d_motor/J

C = 1/(km*tau)

# Constantes de tempo 
tf = 10
timestep = 0.001
Ts = 0.02
vec_time = np.linspace(0, tf, int(tf/timestep)+1)
vec_time_plot = np.linspace(0, tf, int(tf/Ts))
# Vetores

vect_w = [w]
vect_z = [z]
vect_zdot = [z_dot]
vect_z2dot = [z_2dot]
xbar = [[0],[0],[0]]
Xbar = [[0],[0],[0]]
ref = 1

# Teste
cont = 0
vect_aux_u = [0]
vect_Nref = [0]

# Loop

for count_t, timestemp in enumerate(vec_time[1:]):

	if cont == 20:
		aux_xbar1 = np.matmul(disc_ss.B.A,K)
		aux_xbar2 = np.matmul(L,disc_ss.C.A)
		aux_xbar3 = -aux_xbar1-aux_xbar2
		Xbar = np.matmul(disc_ss.A.A+aux_xbar3,Xbar) + L*(z - ref) + auxB*ref

		xbar = Xbar + [[const_x1], [const_x2], [const_x3]]

		aux_u = np.matmul(-K,Xbar)
        
        #Aproximandamente a porcentagem dada por 'Desired Q' * tensão da bateria = 4.656
		U = aux_u[0,0] + 4.656 #4.8 #N[0,0]*ref + 4
		vect_aux_u.append(aux_u[0,0])
		vect_Nref.append(N[0,0]*ref)

		cont = 0

	# MIN DE u É 0,6 q é 5% da vbat
	# max é 95% de 12v 11,4
	w_dot = A*w + B*w**2 + C*U

	# Integração de w

	w += w_dot*timestep

	# Sistema
    
	z_2dot =  -g - (d_corpo/m)*z_dot**2 + 6*(b/m)*w**2

	z_dot += z_2dot*timestep 

	z += z_dot*timestep #+ 0.5*z_2dot*Ts**2

	# Adicionando valores aos vetores

	vect_w.append(w)
	vect_z.append(z)
	vect_zdot.append(z_dot)
	vect_z2dot.append(z_2dot)
	cont += 1

plt.plot(vec_time, vect_z)
plt.show()
#Parte 2
#tentando acrescentar o esp de estados dos angulos
#Conforme PDF Lucas de Paula
l = 0.072456 #metade do lado do quadrado que froma o drone - A no witeboard que usamos



