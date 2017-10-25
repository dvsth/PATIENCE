#########################################################################
# Have you heard the tragedy of Frog Leaper the Wise? It's not a story 	
# a programmer would tell you. Once upon a time, in the days of yore,	
# when PATIENCE was still in development, many brave souls made their 		
# efforts to simulate the dangerous and elusive magnetic field. Euler
# was the bold venturer, ever ready to fight in the name of science. 
# He drew first blood from the simulator, but alas, he could only 
# simulate electric fields with his inferior first-order method. 
# Runge-Kutta came next, taking the baton from 
# Euler as champions of accuracy, pioneering with their second-order
# scheme. It was not enough. Weep, O, you programmers, for under
# Runge-Kutta it was that the dark side started to gain victory.
# Blood was spilled by the megabytes in the form of energy leakage,
# for though our champions were accurate, their valiant scheme 
# failed to conserve the energy that was the lifeblood of the magnetic
# field. Particles spiraled out of control, struggling to maintain a 
# circular path but eventually fading away into the depths of magnitude
# where even the bravest dare not venture.
# 
# Enter Boris.
#########################################################################

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
import ConfigParser

#--------------------- fetch simulation settings

# open the config file to read parameters
config = ConfigParser.ConfigParser()
config.read("configtv.ini")

#analytical mode
analysis_on = config.getint("analysis" , "analysis_on")

# field data
B = [ config.getfloat("fields" , "B_x") , config.getfloat("fields" , "B_y") , config.getfloat("fields" , "B_z") ]
E = [ config.getfloat("fields" , "E_x") , config.getfloat("fields" , "E_y") , config.getfloat("fields" , "E_z") ] 

#time-variance conditions
Omega_E = config.getfloat("fields" , "Omega_E")
Omega_B = config.getfloat("fields" , "Omega_B")
phase_diff = config.getfloat("fields" , "Phase_difference")

# particle data
mass = config.getfloat("particle" , "mass")
charge = config.getfloat("particle" , "charge")

# boundary conditions
pos_0 = [ config.getfloat("boundary-conditions","x_0") , config.getfloat("boundary-conditions","y_0") , config.getfloat("boundary-conditions","z_0") ]
vel_0 = [ config.getfloat("boundary-conditions","v_x_0") , config.getfloat("boundary-conditions","v_y_0") , config.getfloat("boundary-conditions","v_z_0") ]

pos_init = pos_0
vel_init = vel_0
	
#time info
steps = config.getint("time-step","steps")
t_final = config.getfloat("time-step","t_final")
h = t_final/steps

#storage variables for plotting
plot_pos_x = list()
plot_pos_y = list()
plot_pos_z = list()
#plot_time = list()
#plot_field = list()

#--------------------- run the simulation
larmor_radius = (mass * np.linalg.norm(vel_0))**2 / (charge *np.linalg.norm(np.cross(vel_0 , B)))
#analysis mode
if analysis_on == 1 :
	plot_pos_x_aly = list()
	plot_pos_x_aly = list()
	larmor_radius = (mass * np.linalg.norm(vel_0))**2 / (charge *np.linalg.norm(np.cross(vel_0 , B)))
	
#calculate the magnetic field as a function of time
def b_field ( t ) :
	return np.multiply(B, np.sin(Omega_B*t))
	return B

def e_field ( t ) :
	return np.multiply(E, np.sin(Omega_E*t + phase_diff))

#calculate the acceleration
def acceleration ( EField , BField , velocity ) :
	a_0 = np.add(np.multiply( charge/mass , np.cross( vel_0 , BField ) ), np.multiply( charge/mass , EField ) )
	return a_0
	
#calculate the position
def positionUpdate ( h , x_prev , v_prev , v_new ) :
	return np.add(x_prev , np.multiply(h/2.0 , np.add(v_new , v_prev)))
	
#algorithm bootstrapping step
a_0 = acceleration ( E , B , vel_0 )
v_nminushalf = np.subtract(vel_0 , np.multiply(h/2. , a_0) )
vel_0 = v_nminushalf

for step in range(steps) :
	
	#calculate total time elapsed
	time = h*step 
	
	#get instantaneous fields
	B_inst = b_field(time)
	E_inst = e_field(time)
	
	#boris variables
	t_vector = np.multiply(charge * h / mass * 0.5 , B_inst)
	s_vector = np.divide( np.multiply(2. , t_vector) , 1 + (np.linalg.norm(t_vector)**2))
	
	#boris first pass
	v_minus = np.add( vel_0 , np.multiply( h/2. , np.multiply( charge/mass , E_inst ) ) )
	v_nminushalf = np.subtract( v_minus , np.multiply( charge/mass , E_inst ) )
	
	#perform the B-rotation
	v_prime = np.add( v_minus , np.cross( v_minus , t_vector ) )
	
	#boris second pass
	v_plus = np.add( v_minus , np.cross( v_prime , s_vector ) )
	v_nplushalf = np.add( v_plus , np.multiply( charge/mass , E_inst ) )
	
	#update position
	pos_1 = positionUpdate( h , pos_0 , v_nminushalf , v_nplushalf )
	
	#add data to plotting array
	plot_pos_x.append(pos_0[0])
	plot_pos_y.append(pos_0[1])
	plot_pos_z.append(pos_0[2])
	
	#carry over variables for next iteration
	pos_0 = pos_1
	vel_0 = v_nplushalf

#--------------------- generate the plot
fig = plt.figure()

fig.suptitle('Trajectory in EM field \n computed using Boris pusher \n in ' + str(steps) + " steps" , fontweight='bold' , fontsize=10)

ax = fig.add_subplot(111, projection='3d')
ax.text2D(0.05, 0.95, "Omega values:\nE: " + str(Omega_E) + "\nB: " + str(Omega_B), transform=ax.transAxes)

ax.scatter(plot_pos_x[0] , plot_pos_y[0] , plot_pos_z[0] , color='r' , s=200,marker='o' , label='initial position')
ax.plot(plot_pos_x, plot_pos_y, plot_pos_z, color='b', label='trajectory')
ax.scatter(plot_pos_x[-1] , plot_pos_y[-1] , plot_pos_z[-1] , color='r' , s=200,marker='x' , label='final position')

#!--TODO--!
	#def animate(i):
		#ax.plot(plot_pos_x[0 : i], plot_pos_y[0 : i], plot_pos_z[0 : i])

		
	#dplot = fig.add_subplot(122)
	#dplot.plot(plot_time, plot_pos_z, 'b-', label='position')
	#dplot.plot(plot_time, plot_field, 'r-', label='field')
	#plt.legend(loc='best' , fontsize=20)
	#ax.set_xlim(-100)
	#ax.set_ylim(-99)
	#ax.set_zlim(98)

	#anim = FuncAnimation(fig, animate, frames=1)

	#ax.text(pos_init[0] + 1, pos_init[1] + 1, pos_init[2] + 1, 'starting position' + str(pos_init))

	#fig.savefig('myimage.png', format='png', dpi=1200)

ax.set_xlabel('x-position (m) | E = ' + str(E[0]) +  ' | B = ' + str(B[0]) , fontweight='light' , fontsize=15)
ax.set_ylabel('y-position (m) | E = ' + str(E[1]) +  ' | B = ' + str(B[1]) , fontweight='light' , fontsize=15)
ax.set_zlabel('z-position (m) | E = ' + str(E[2]) +  ' | B = ' + str(B[2]) , fontweight='light' , fontsize=15)
ax.legend()

#ax.set_aspect(1)
plt.legend(loc='best' , fontsize=15)


plt.show()
