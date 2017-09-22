#
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
# And then there was Boris.
#
# enter Boris

import numpy as np
import matplotlib.pyplot as plt
import ConfigParser

#--------------------- fetch simulation settings

# open the config file to read parameters
config = ConfigParser.ConfigParser()
config.read("config.ini")

#analytical mode
analysis_on = config.getint("analysis" , "analysis_on")

# field data
B = [ config.getfloat("field" , "B_x") , config.getfloat("field" , "B_y") , config.getfloat("field" , "B_z") ] 

# particle data
mass = config.getfloat("particle" , "mass")
charge = config.getfloat("particle" , "charge")

# boundary conditions
pos_0 = [ config.getfloat("boundary-conditions","x_0") , config.getfloat("boundary-conditions","y_0") , config.getfloat("boundary-conditions","z_0") ]
vel_0 = [ config.getfloat("boundary-conditions","v_x_0") , config.getfloat("boundary-conditions","v_y_0") , config.getfloat("boundary-conditions","v_z_0") ]

#time info
steps = config.getint("time-step","steps")
t_final = config.getfloat("time-step","t_final")
h = t_final/steps

#storage variables for plotting
plot_pos_x = list()
plot_pos_y = list()

#boris variables
t_vector = np.multiply(charge * h / mass * 0.5 , B)
s_vector = np.divide( np.multiply(2. , t_vector) , 1 + (np.linalg.norm(t_vector)**2))


#--------------------- run the simulation

#analysis mode
if analysis_on == 1 :
	plot_pos_x_aly = list()
	plot_pos_x_aly = list()
	larmor_radius = (mass * np.linalg.norm(vel_0))**2 / (charge*np.linalg.norm(np.cross(vel_0 , B)))
	vel_init = vel_0

#calculate the magnetic field as a function of time
def field ( t ) :
	return B
	
#calculate the acceleration
def acceleration ( BField , velocity ) :
	return np.multiply(charge/mass , np.cross( velocity , BField ))
	
#calculate the velocity	
def velocityUpdate ( h , x_prev , v_prev , t ) :
	return np.add(v_prev ,  np.multiply( h , acceleration( field(t) , v_prev )))
	
#calculate the position
def positionUpdate ( h , x_prev , v_prev , v_new ) :
	return np.add(x_prev , np.multiply(h/2.0 , np.add(v_new , v_prev)))
	
#algorithm
for step in range(steps) :
	
	a_0 = np.multiply( charge/mass , np.cross( vel_0 , B ) )
	
	v_minus = np.add( vel_0 , np.multiply( h/2. , a_0 ) )
	v_prime = np.add( v_minus , np.cross( v_minus , t_vector ) )
	v_plus = np.add( v_minus , np.cross( v_prime , s_vector ) )
	
	vel_0 = np.divide( np.add( v_minus, v_plus ) , 2.0 )

	plot_pos_x.append(vel_0[0])
	plot_pos_y.append(vel_0[1])
		
#--------------------- generate the plot
fig = plt.figure()

fig.suptitle('Trajectory for a proton \n computed using Boris pusher \n in ' + str(steps) + " steps" , fontweight='bold')
pltPos = fig.add_subplot(111)

pltPos.grid(b=True, color='k', linestyle='--')
pltPos.set_xlabel('x-velocity (m/s)' , fontweight='bold')
pltPos.set_ylabel('y-velocity (m/s)' , fontweight='bold')

pltPos.plot(plot_pos_x , plot_pos_y , 'b-', label='Boris path')
#pltPos.scatter(plot_pos_x[0] ,plot_pos_y[0], color='k' , s=100,marker='o' , label='initial position')
#pltPos.scatter(plot_pos_x[-1] , plot_pos_y[-1], color='k' , s=100,marker='x' , label='final position')
	
#analysis mode
if analysis_on == 1 :
	t_arr = np.linspace(0.0, t_final, steps)
	Bmag = np.linalg.norm(B)
	Omega_L = charge*Bmag/mass
	x_ana = vel_init[0]*np.cos(Omega_L*t_arr) + vel_init[1]*np.sin(Omega_L*t_arr)
	y_ana = -1 * vel_init[0]*np.sin(Omega_L*t_arr) + vel_init[1]*np.cos(Omega_L*t_arr)
	pltPos.plot(x_ana , y_ana , 'g--' , label='Analytical path')
	pltPos.set_title('Larmor radius = ' + str(larmor_radius) + ' metres')

pltPos.set_aspect('1.0')
plt.legend(loc='best')
plt.show()



	
	


	
	





